!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! SUBMODULE: green_chebyshev
!
!> Chebyshev-specific Green-function reconstruction and moment processing.
!------------------------------------------------------------------------------

submodule (green_mod) green_chebyshev

   use chebyshev_fast_mod, only: cheb_green_fast, cheb_green_dp, cheb_green_complex
   use math_mod, only: i_unit, pi, jackson_kernel, t_polynomial, &
                       hubbard_spectral_density_element, project_spectral_density_2x2
   use basis_mod, only: norb, spin_off
   use mpi_mod, only: atoms_per_process, start_atom, end_atom, g2l_map
   use rsrec_cuda_plugin_mod, only: rsrec_cuda_backend, get_gpu_context, rsrec_cuda_plugin_compiled
   use logger_mod, only: g_logger
   use string_mod, only: fmt
   implicit none

contains

   module subroutine chebyshev_green_complex_ij(this, istart, z, g_out)
      class(green), intent(inout) :: this
      integer, intent(in) :: istart
      complex(rp), intent(in) :: z
      complex(rp), intent(out) :: g_out(:, :, :)
      real(rp) :: a, b, emin_win, emax_win
      integer :: n, n_mom

      if (size(g_out, 1) /= nb .or. size(g_out, 2) /= nb .or. size(g_out, 3) /= 4) then
         error stop 'chebyshev_green_complex_ij: output must have shape (nb,nb,4)'
      end if
      if (.not. allocated(this%recursion%mu_n) .or. istart < 1 .or. &
          istart + 3 > size(this%recursion%mu_n, 4)) then
         error stop 'chebyshev_green_complex_ij: requested pair moments are unavailable'
      end if

      call this%recursion%resolve_chebyshev_window(emin_win, emax_win)
      a = (emax_win - emin_win)/(2 - 0.3_rp)
      b = (emax_win + emin_win)/2.0_rp
      n_mom = size(this%recursion%mu_n, 3)
      g_out = cmplx(0.0_rp, 0.0_rp, rp)
      do n = 1, 4
         call cheb_green_complex(this%recursion%mu_n(:, :, :, istart + n - 1), nb, n_mom, z, a, b, g_out(:, :, n))
      end do
   end subroutine chebyshev_green_complex_ij

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Uses the moments form the Chebyshev recursions to calculate the onsite GF
   !---------------------------------------------------------------------------
   module subroutine chebyshev_green_ij(this, istart)
      class(green), intent(inout) :: this
      integer, intent(in) :: istart
      ! Local variables
      real(rp), dimension(:), allocatable :: kernel
      real(rp), dimension(:, :), allocatable :: polycheb
      real(rp), dimension(:), allocatable :: w, wscale
      real(rp) :: wstep, eps, wmin, wmax, a, b, emin_win, emax_win
      complex(rp) :: exp_factor
      integer :: ie, i, j, k, l, m, n

      this%g0 = 0.0d0

      allocate (kernel(size(this%recursion%mu_n, 3)), polycheb(this%en%channels_ldos + 10, 0:size(this%recursion%mu_n, 3)), w(this%en%channels_ldos + 10), &
                wscale(this%en%channels_ldos + 10))
      ! Defining rescaling coeficients
      call this%recursion%resolve_chebyshev_window(emin_win, emax_win)
      a = (emax_win - emin_win)/(2 - 0.3_rp)
      b = (emax_win + emin_win)/2.0_rp

      wscale(:) = (this%en%ene(:) - b)/a

      ! Calculating the Jackson Kernel
      call jackson_kernel(size(this%recursion%mu_n, 3), kernel)

      ! Calculating the Lorentz Kernel
!    call lorentz_kernel(this%control%lld, kernel, 4.0d0)

      do n = 1, 4 ! Loop on the number of on-site GFs to calculate the inter-site GFs
         ! Multiply the moments with the kernel
         do l = 1, nb
            do m = 1, nb
               this%recursion%mu_ng(l, m, :, n + istart - 1) = this%recursion%mu_n(l, m, :, n + istart - 1)*kernel(:)
            end do
         end do
         this%recursion%mu_ng(:, :, 2:size(kernel), n + istart - 1) = this%recursion%mu_ng(:, :, 2:size(kernel), n + istart - 1)*2.0_rp

         ! Calculate the Chebyshev polynomials
         call t_polynomial(size(w), size(kernel), wscale(:), polycheb)

         ! Calculate the density of states
         !$omp parallel do default(shared) private(ie, i, exp_factor, l,m)
         do ie = 1, this%en%channels_ldos + 10
            do i = 1, size(kernel)
               exp_factor = -i_unit*exp(-i_unit*(i - 1)*acos(wscale(ie)))
               do l = 1, nb
                  do m = 1, nb
                     this%g0(l, m, ie, n) = this%g0(l, m, ie, n) + this%recursion%mu_ng(l, m, i, n + istart - 1)*exp_factor
                  end do
               end do
            end do
            do l = 1, nb
               do m = 1, nb
                  this%g0(l, m, ie, n) = this%g0(l, m, ie, n)/((sqrt((a**2) - ((this%en%ene(ie) - b)**2))))
               end do
            end do
         end do
         !$omp end parallel do
      end do  ! End loop on n

      deallocate (kernel, polycheb, w, wscale)
   end subroutine chebyshev_green_ij

   !---------------------------------------------------------------------------
   !> @brief GPU drop-in for chebyshev_green_ij (intersite Chebyshev GF).
   !---------------------------------------------------------------------------
   module subroutine chebyshev_green_ij_gpu(this, istart)
      class(green), intent(inout) :: this
      integer, intent(in) :: istart
      complex(rp), allocatable :: mu_local(:, :, :, :), g0_local(:, :, :, :)
      real(rp) :: a, b, emin_win, emax_win
      integer :: nv, n_mom
      type(rsrec_cuda_backend), pointer :: gpu_backend

      if (.not. rsrec_cuda_plugin_compiled()) then
         call this%chebyshev_green_ij(istart)
         return
      end if

      this%g0 = 0.0d0
      call this%recursion%resolve_chebyshev_window(emin_win, emax_win)
      a = (emax_win - emin_win)/(2 - 0.3_rp)
      b = (emax_win + emin_win)/2.0_rp
      nv = this%en%channels_ldos + 10
      n_mom = size(this%recursion%mu_n, 3)

      allocate (mu_local(nb, nb, n_mom, 4), g0_local(nb, nb, nv, 4))
      mu_local(:, :, :, 1:4) = this%recursion%mu_n(:, :, :, istart:istart + 3)

      gpu_backend => get_gpu_context(1, nb, 1, 1, 1)
      call gpu_backend%set_precision(this%recursion%gpu_precision_code())

      call gpu_backend%chebyshev_dos(mu_local, this%en%ene(1:nv), a, b, g0_local)
      this%g0(:, :, 1:nv, 1:4) = g0_local(:, :, :, 1:4)
      deallocate (mu_local, g0_local)
   end subroutine chebyshev_green_ij_gpu

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Uses the moments form the Chebyshev recursions to calculate the onsite GF
   !> using a small complex eta into the energy channel
   !---------------------------------------------------------------------------
   module subroutine chebyshev_green_ij_eta(this, istart, eta, fermi_point, g_ef)
      class(green), intent(inout) :: this
      integer, intent(in) :: istart
      complex(rp), dimension(nb, nb, 4), intent(inout) :: g_ef
      complex(rp), intent(in) :: eta
      integer, intent(in) :: fermi_point
      real(rp) :: a, b, emin_win, emax_win
      integer :: n

      g_ef = 0.0d0

      call this%recursion%resolve_chebyshev_window(emin_win, emax_win)
      a = (emax_win - emin_win)/(2 - 0.3_rp)
      b = (emax_win + emin_win)/2.0_rp

      do n = 1, 4 ! Loop on the number of on-site GFs to calculate the inter-site GFs
         call cheb_green_complex(this%recursion%mu_n(:, :, :, n + istart - 1), &
            nb, size(this%recursion%mu_n, 3), cmplx(this%en%ene(fermi_point), 0.0_rp, rp) + eta, a, b, g_ef(:, :, n))
      end do  ! End loop on n
   end subroutine chebyshev_green_ij_eta

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Uses the moments form the Chebyshev recursions to calculate the onsite GF
   !---------------------------------------------------------------------------
   module subroutine chebyshev_green(this)
      class(green), intent(inout) :: this

      ! Deprecated wrapper retained for type-bound callers; use chebyshev_green_core.
      call chebyshev_green_core(this, use_gpu=this%control%gpu_plugin)
   end subroutine chebyshev_green

   !---------------------------------------------------------------------------
   !> @brief Deprecated wrapper; use chebyshev_green_core(this, use_gpu=.true.).
   !---------------------------------------------------------------------------
   module subroutine chebyshev_green_gpu(this)
      class(green), intent(inout) :: this

      ! Deprecated wrapper retained for type-bound callers; use chebyshev_green_core.
      call chebyshev_green_core(this, use_gpu=.true.)
   end subroutine chebyshev_green_gpu

   !---------------------------------------------------------------------------
   !> @brief Dispatcher for the on-site Chebyshev Green’s-function/DOS
   !---------------------------------------------------------------------------
   module subroutine chebyshev_dos_dispatch(this)
      class(green), intent(inout) :: this

      if (this%control%gpu_plugin) then
         if (this%control%nsp > 4) then
            call g_logger%warning('gpu_plugin requested for Chebyshev DOS, '// &
               'but nsp > 4. Using legacy path.', __FILE__, __LINE__)
            call chebyshev_green_core(this, use_gpu=.false.)
            return
         end if
         call chebyshev_green_core(this, use_gpu=.true.)
         return
      end if

      call chebyshev_green_core(this, use_gpu=.false.)

   end subroutine chebyshev_dos_dispatch

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Uses the moments from the Chebyshev recursions to calculate the onsite GF
   !> using a small complex eta into the energy channel
   !---------------------------------------------------------------------------
   module subroutine chebyshev_green_eta(this, eta, fermi_point, g_ef)
      class(green), intent(inout) :: this
      complex(rp), dimension(nb, nb, atoms_per_process), intent(inout) :: g_ef
      complex(rp), intent(in) :: eta
      integer, intent(in) :: fermi_point

      ! Deprecated wrapper retained for type-bound callers; use chebyshev_green_core.
      call chebyshev_green_core(this, eta, fermi_point, g_ef)
   end subroutine chebyshev_green_eta

   !> @brief On-site Chebyshev Green’s function, on-mesh or at a single
   !>        eta-shifted energy, CPU or GPU.
   module subroutine chebyshev_green_core(this, eta, fermi_point, g_ef, use_gpu)
      class(green), intent(inout) :: this
      complex(rp), intent(in), optional :: eta
      integer, intent(in), optional :: fermi_point
      complex(rp), dimension(nb, nb, atoms_per_process), intent(inout), optional :: g_ef
      logical, intent(in), optional :: use_gpu
      ! Local variables
      real(rp), dimension(size(this%recursion%mu_n, 3)) :: kernel
      real(rp), dimension(this%en%channels_ldos + 10, 0:size(this%recursion%mu_n, 3)) :: polycheb
      real(rp), dimension(this%en%channels_ldos + 10) :: w, wscale
      real(rp) :: wstep, eps, wmin, wmax, a, b, emin_win, emax_win
      complex(rp) :: exp_factor
      integer :: ie, i, j, k, l, m, n, nv, n_glob, n_local
      complex(rp), allocatable :: mu_local(:, :, :, :)
      complex(rp), allocatable :: g0_local(:, :, :, :)
      type(rsrec_cuda_backend), pointer :: gpu_backend
      logical :: do_gpu

      do_gpu = this%control%gpu_plugin
      if (present(use_gpu)) do_gpu = use_gpu
      if (do_gpu .and. .not. rsrec_cuda_plugin_compiled()) then
         call g_logger%warning('GPU Chebyshev Green reconstruction requested, but this executable '// &
            'was built without ENABLE_CUDA_PLUGIN. Falling back to the CPU reconstruction.', &
            __FILE__, __LINE__)
         do_gpu = .false.
      end if
      ! G1.3 keeps transport moments resident in the shared CUDA context.
      ! The DOS helper uses a different context shape and would otherwise
      ! destroy that residency before conductivity reconstruction. Keep this
      ! independent DOS stage on its canonical CPU path in that case.
      if (do_gpu .and. associated(this%recursion%gpu_backend)) then
         if (this%recursion%gpu_backend%resident_moments_available()) then
            call g_logger%info('CUDA transport moments resident; keeping Chebyshev DOS on CPU.', __FILE__, __LINE__)
            do_gpu = .false.
         end if
      end if
      if (present(eta)) then
         g_ef = 0.0d0
      else
         this%g0 = 0.0d0
      end if

      ! Defining rescaling coeficients
      call this%recursion%resolve_chebyshev_window(emin_win, emax_win)
      a = (emax_win - emin_win)/(2 - 0.3_rp)
      b = (emax_win + emin_win)/2.0_rp

      ! Number of DOS points
      nv = this%en%channels_ldos + 10

      if (.not. present(eta) .and. do_gpu) then
         call jackson_kernel(size(this%recursion%mu_n, 3), kernel)
         do n_glob = start_atom, end_atom
            n = g2l_map(n_glob)
            do l = 1, nb
               do m = 1, nb
                  this%recursion%mu_ng(l, m, :, n) = this%recursion%mu_n(l, m, :, n)*kernel(:)
               end do
            end do
            this%recursion%mu_ng(:, :, 2:size(kernel), n) = &
               this%recursion%mu_ng(:, :, 2:size(kernel), n)*2.0_rp
         end do

         gpu_backend => get_gpu_context(1, nb, 1, 1, 1)
         call gpu_backend%set_precision(this%recursion%gpu_precision_code())

         call gpu_backend%chebyshev_dos(this%recursion%mu_n(:, :, :, start_atom:end_atom), &
                                         this%en%ene(1:nv), a, b,                         &
                                         this%g0(:, :, 1:nv, start_atom:end_atom))
         call enforce_chebyshev_spectral_positivity(this, nv)
         return
      end if

      if (.not. present(eta)) then
         n_local = end_atom - start_atom + 1
         if (n_local <= 0) return

         allocate (mu_local(nb, nb, size(this%recursion%mu_n, 3), n_local))
         allocate (g0_local(nb, nb, nv, n_local))

         do n_glob = start_atom, end_atom
            n = g2l_map(n_glob)
            mu_local(:, :, :, n) = this%recursion%mu_n(:, :, :, n)
         end do

         ! fast_dp (default) and legacy reconstruct in double precision; the other backends keep fp32.
         if (trim(this%control%cheb_backend) == 'fast_dp' .or. trim(this%control%cheb_backend) == 'legacy') then
            call cheb_green_dp(mu_local, nb, size(mu_local, 3), n_local, this%en%ene(1:nv), nv, a, b, g0_local)
         else
            call cheb_green_fast(mu_local, nb, size(mu_local, 3), n_local, this%en%ene(1:nv), nv, a, b, g0_local)
         end if

         do n_glob = start_atom, end_atom
            n = g2l_map(n_glob)
            this%g0(:, :, 1:nv, n) = g0_local(:, :, :, n)
         end do

         deallocate (g0_local, mu_local)
         call enforce_chebyshev_spectral_positivity(this, nv)
         return
      end if

      wscale(:) = (this%en%ene(:) - b)/a

      ! Calculating the Jackson Kernel
      call jackson_kernel(size(this%recursion%mu_n, 3), kernel)

      ! Calculating the Lorentz Kernel
      ! call lorentz_kernel(this%control%lld, kernel, 4.0d0)

      do n_glob = start_atom, end_atom ! Loop on self-consistent atoms
         n = g2l_map(n_glob)
         ! Multiply the moments with the kernel
         do i = 1, nb
            this%recursion%mu_ng(i, i, :, n) = this%recursion%mu_n(i, i, :, n)*kernel(:)
         end do

         this%recursion%mu_ng(:, :, 2:size(kernel), n) = this%recursion%mu_ng(:, :, 2:size(kernel), n)*2.0_rp

         !do i=1, size(kernel)
         !  write(400+n, *) i, sum(this%recursion%mu_n(n, i, 1:18, 1:18))
         !end do

         ! Calculate the Chebyshev polynomials
         call t_polynomial(size(w), size(kernel), wscale(:), polycheb)

         ! Calculate the density of states
         !$omp parallel do default(shared) private(ie, i, exp_factor, l,m)
         do ie = fermi_point, fermi_point
            do i = 1, size(kernel)
               exp_factor = -i_unit*exp(-i_unit*(i - 1)*acos(((this%en%ene(ie) + eta) - b)/a))
               do l = 1, nb
                  do m = 1, nb
                     g_ef(l, m, n) = g_ef(l, m, n) + this%recursion%mu_ng(l, m, i, n)*exp_factor
                  end do
               end do
            end do
            do l = 1, nb
               do m = 1, nb
                  g_ef(l, m, n) = g_ef(l, m, n)/((sqrt((a**2) - (((this%en%ene(ie) + eta) - b)**2))))
               end do
            end do
         end do
         !$omp end parallel do
      end do  ! End loop on self-consistent atoms
   end subroutine chebyshev_green_core

   !> Restore the local spectral positivity required by the real-axis DOS
   !> route.  The finite Jackson expansion is non-negative in exact arithmetic,
   !> but a low-precision moment or GEMM path can make a tiny negative DOS.
   !> Work in each orbital's spin block so the total DOS, spin-resolved DOS,
   !> and the charge/moment integrations all consume the same physical density.
   !> Eta-shifted and intersite Green functions are intentionally excluded.
   !> With cheb_lp the clipped weight is taken out again: each atom, orbital and
   !> spin channel is rescaled so that its integrated weight equals the weight
   !> before the projection (the coherence by the geometric mean of the two
   !> factors, which keeps the 2x2 block positive semidefinite). The extended
   !> moments fix the zeroth moment, so the clipping would otherwise add charge
   !> and move the Fermi level.
   subroutine enforce_chebyshev_spectral_positivity(this, nv)
      class(green), intent(inout) :: this
      integer, intent(in) :: nv
      integer :: n_glob, n, ie, orbital, corrections
      real(rp) :: a_uu, a_dd, a_uu_before, a_dd_before, max_adjustment
      complex(rp) :: a_ud, a_ud_before, guu, gdd, gud, gdu, h_ud
      logical :: corrected
      real(rp), dimension(norb) :: dq_uu, dq_dd, q_uu, q_dd, s_uu, s_dd
      real(rp) :: s_ud, max_rescale

      corrections = 0
      max_adjustment = 0.0_rp
      max_rescale = 0.0_rp
      do n_glob = start_atom, end_atom
         n = g2l_map(n_glob)
         dq_uu = 0.0_rp
         dq_dd = 0.0_rp
         q_uu = 0.0_rp
         q_dd = 0.0_rp
         do ie = 1, nv
            do orbital = 1, norb
               guu = this%g0(orbital, orbital, ie, n)
               gdd = this%g0(orbital + spin_off, orbital + spin_off, ie, n)
               gud = this%g0(orbital, orbital + spin_off, ie, n)
               gdu = this%g0(orbital + spin_off, orbital, ie, n)

               a_uu = -aimag(guu)/pi
               a_dd = -aimag(gdd)/pi
               a_ud = hubbard_spectral_density_element(gud, gdu)
               a_uu_before = a_uu
               a_dd_before = a_dd
               a_ud_before = a_ud
               call project_spectral_density_2x2(a_uu, a_ud, a_dd, corrected)
               q_uu(orbital) = q_uu(orbital) + a_uu
               q_dd(orbital) = q_dd(orbital) + a_dd
               if (.not. corrected) cycle

               dq_uu(orbital) = dq_uu(orbital) + a_uu - a_uu_before
               dq_dd(orbital) = dq_dd(orbital) + a_dd - a_dd_before
               corrections = corrections + 1
               max_adjustment = max(max_adjustment, abs(a_uu - a_uu_before), &
                                    abs(a_dd - a_dd_before), abs(a_ud - a_ud_before))

               ! G = Hermitian(G) - i*pi*A.  Preserve the dispersive
               ! Hermitian part and replace only the non-causal spectral part.
               h_ud = 0.5_rp*(gud + conjg(gdu))
               this%g0(orbital, orbital, ie, n) = cmplx(real(guu, rp), -pi*a_uu, rp)
               this%g0(orbital + spin_off, orbital + spin_off, ie, n) = &
                  cmplx(real(gdd, rp), -pi*a_dd, rp)
               this%g0(orbital, orbital + spin_off, ie, n) = h_ud - i_unit*pi*a_ud
               this%g0(orbital + spin_off, orbital, ie, n) = conjg(h_ud) - i_unit*pi*conjg(a_ud)
            end do
         end do

         if (this%recursion%lp%enabled) then
            s_uu = 1.0_rp
            s_dd = 1.0_rp
            do orbital = 1, norb
               if (q_uu(orbital) > 0.0_rp) s_uu(orbital) = 1.0_rp - dq_uu(orbital)/q_uu(orbital)
               if (q_dd(orbital) > 0.0_rp) s_dd(orbital) = 1.0_rp - dq_dd(orbital)/q_dd(orbital)
            end do
            max_rescale = max(max_rescale, maxval(abs(1.0_rp - s_uu)), maxval(abs(1.0_rp - s_dd)))
            do ie = 1, nv
               do orbital = 1, norb
                  s_ud = sqrt(s_uu(orbital)*s_dd(orbital))
                  guu = this%g0(orbital, orbital, ie, n)
                  gdd = this%g0(orbital + spin_off, orbital + spin_off, ie, n)
                  gud = this%g0(orbital, orbital + spin_off, ie, n)
                  gdu = this%g0(orbital + spin_off, orbital, ie, n)
                  h_ud = 0.5_rp*(gud + conjg(gdu))
                  this%g0(orbital, orbital, ie, n) = cmplx(real(guu, rp), s_uu(orbital)*aimag(guu), rp)
                  this%g0(orbital + spin_off, orbital + spin_off, ie, n) = cmplx(real(gdd, rp), s_dd(orbital)*aimag(gdd), rp)
                  this%g0(orbital, orbital + spin_off, ie, n) = h_ud + s_ud*(gud - h_ud)
                  this%g0(orbital + spin_off, orbital, ie, n) = conjg(h_ud) + s_ud*(gdu - conjg(h_ud))
               end do
            end do
         end if
      end do

      if (corrections > 0) then
         call g_logger%info('Chebyshev spectral positivity projection: local_blocks='// &
            fmt('i0', corrections)//' max_adjustment='//fmt('e12.4', max_adjustment), __FILE__, __LINE__)
         if (this%recursion%lp%enabled) call g_logger%info('Chebyshev LP charge-conserving rescale: max|1-s|='// &
            fmt('e12.4', max_rescale), __FILE__, __LINE__)
      end if
   end subroutine enforce_chebyshev_spectral_positivity

end submodule green_chebyshev
