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
   use math_mod, only: i_unit, pi, jackson_kernel, lorentz_kernel, &
                       hubbard_spectral_density_element, project_spectral_density_2x2, &
                       cheb_bilinear_weights, cheb_bilinear_transform, kB_simpson_f, &
                       cheb_energy_nodes, cheb_green_point_weights, cheb_linear_weights
   use basis_mod, only: norb, spin_off
   use mpi_mod, only: atoms_per_process, start_atom, end_atom, g2l_map
   use rsrec_cuda_plugin_mod, only: rsrec_cuda_backend, get_gpu_context, rsrec_cuda_plugin_compiled
   use logger_mod, only: g_logger
   use string_mod, only: fmt
   implicit none

contains

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Uses the moments form the Chebyshev recursions to calculate the onsite GF
   !---------------------------------------------------------------------------
   module subroutine chebyshev_green_ij(this, istart)
      class(green), intent(inout) :: this
      integer, intent(in) :: istart
      integer :: n

      this%g0 = 0.0d0
      do n = 1, 4 ! Loop on the number of on-site GFs to calculate the inter-site GFs
         call chebyshev_reconstruct_block(this, this%recursion%mu_n(:, :, :, n + istart - 1), this%g0(:, :, :, n))
      end do
   end subroutine chebyshev_green_ij

   !> Real-mesh Green function g(:,:,ie) on en%ene of one block of moments mu(:,:,:) (kernel-damped Chebyshev series).
   subroutine chebyshev_reconstruct_block(this, mu, g)
      class(green), intent(inout) :: this
      complex(rp), intent(in) :: mu(:, :, :)
      complex(rp), intent(inout) :: g(:, :, :)
      ! Local variables
      real(rp), dimension(:), allocatable :: kernel
      real(rp), dimension(:), allocatable :: wscale
      complex(rp), dimension(:, :, :), allocatable :: mug
      real(rp) :: a, b, emin_win, emax_win
      complex(rp), dimension(:, :), allocatable :: ecoef, gs
      integer :: ie, i, l, m, n, nv

      allocate (kernel(size(mu, 3)), wscale(this%en%channels_ldos + 10))
      allocate (mug(nb, nb, size(mu, 3)))
      ! Defining rescaling coeficients
      call this%recursion%resolve_chebyshev_window(emin_win, emax_win)
      a = (emax_win - emin_win)/(2 - 0.3_rp)
      b = (emax_win + emin_win)/2.0_rp

      wscale(:) = (this%en%ene(:) - b)/a

      ! Calculating the damping kernel selected by control%cheb_kernel (Jackson by default)
      call cheb_kernel_vector(this, size(mu, 3), kernel)

      ! Multiply the moments with the kernel
      do l = 1, nb
         do m = 1, nb
            mug(l, m, :) = mu(l, m, :)*kernel(:)
         end do
      end do
      mug(:, :, 2:size(kernel)) = mug(:, :, 2:size(kernel))*2.0_rp

      ! g(:,ie) = (g(:,ie) + sum_m mug(:,m) E(m,ie)) / sqrt(a^2-(E_ie-b)^2), E(m,ie) = -i exp(-i (m-1) acos(x_ie))
      nv = this%en%channels_ldos + 10
      n = size(kernel)
      allocate (ecoef(n, nv), gs(nb*nb, nv))
      !$omp parallel do default(shared) private(ie, i)
      do ie = 1, nv
         do i = 1, n
            ecoef(i, ie) = -i_unit*exp(-i_unit*(i - 1)*acos(wscale(ie)))
         end do
      end do
      !$omp end parallel do
      call zgemm('N', 'N', nb*nb, nv, n, (1.0_rp, 0.0_rp), mug, nb*nb, ecoef, n, (0.0_rp, 0.0_rp), gs, nb*nb)
      do ie = 1, nv
         g(:, :, ie) = (g(:, :, ie) + reshape(gs(:, ie), [nb, nb]))/(sqrt((a**2) - ((this%en%ene(ie) - b)**2)))
      end do

      deallocate (kernel, wscale, mug, ecoef, gs)
   end subroutine chebyshev_reconstruct_block

   !> M_ij(m) = M_ji(m)^H for every moment order m.
   function moments_dagger(mu) result(md)
      complex(rp), intent(in) :: mu(:, :, :)
      complex(rp) :: md(size(mu, 2), size(mu, 1), size(mu, 3))
      integer :: m

      do m = 1, size(mu, 3)
         md(:, :, m) = conjg(transpose(mu(:, :, m)))
      end do
   end function moments_dagger

   !> Real-mesh G_ij and G_ji (nb,nb,nv) of local pair ij_loc from the one-shot moments M_ji = os_mu(:,:,:,ient),
   !> M_ij = M_ji^H. Self pair: G_ij = G_ji = G_ii.
   module subroutine chebyshev_green_ij_oneshot(this, ij_loc, gij, gji)
      class(green), intent(inout) :: this
      integer, intent(in) :: ij_loc
      complex(rp), intent(out) :: gij(:, :, :), gji(:, :, :)
      complex(rp), allocatable :: mu_local(:, :, :, :), g0_local(:, :, :, :)
      real(rp), allocatable :: damp(:)
      real(rp) :: a, b, emin_win, emax_win
      integer :: ient, nv, n_mom
      logical :: self_pair
      type(rsrec_cuda_backend), pointer :: gpu_backend

      ient = this%recursion%req%pair_entry(ij_loc)
      self_pair = this%lattice%ijpair(start_atom + ij_loc - 1, 1) == this%lattice%ijpair(start_atom + ij_loc - 1, 2)
      nv = this%en%channels_ldos + 10
      n_mom = size(this%recursion%os_mu, 3)

      if (this%control%gpu_plugin .and. rsrec_cuda_plugin_compiled()) then
         call this%recursion%resolve_chebyshev_window(emin_win, emax_win)
         a = (emax_win - emin_win)/(2 - 0.3_rp)
         b = (emax_win + emin_win)/2.0_rp
         allocate (mu_local(nb, nb, n_mom, 2), g0_local(nb, nb, nv, 2))
         mu_local(:, :, :, 1) = moments_dagger(this%recursion%os_mu(:, :, :, ient))
         mu_local(:, :, :, 2) = this%recursion%os_mu(:, :, :, ient)
         gpu_backend => get_gpu_context(1, nb, 1, 1, 1)
         call gpu_backend%set_precision(this%recursion%gpu_precision_code())
         call cheb_damping(this, n_mom, damp)
         call gpu_backend%chebyshev_dos(mu_local, this%en%ene(1:nv), a, b, g0_local, damp)
         gij = g0_local(:, :, :, 1)
         gji = g0_local(:, :, :, 2)
         deallocate (mu_local, g0_local)
      else
         gji = 0.0d0
         call chebyshev_reconstruct_block(this, this%recursion%os_mu(:, :, :, ient), gji)
         gij = 0.0d0
         if (.not. self_pair) call chebyshev_reconstruct_block(this, moments_dagger(this%recursion%os_mu(:, :, :, ient)), gij)
      end if
      if (self_pair) gij = gji
   end subroutine chebyshev_green_ij_oneshot

   !> Number of stored moments of the inter-site recursion data (cheb_lp extension included).
   function intersite_moment_count(this) result(nmom)
      class(green), intent(in) :: this
      integer :: nmom

      if (this%control%intersite_mode == 'one_shot') then
         nmom = size(this%recursion%os_mu, 3)
      else
         nmom = size(this%recursion%mu_n, 3)
      end if
   end function intersite_moment_count

   !> Weights h(0:2*nmom-2, 0:2) of the bilinear moment-space energy integral over [energy_min, fermi] (Fermi function
   !> for temperature > 0), with the window map of chebyshev_reconstruct_block, and the damping kernel(nmom).
   !> Optional lin(nmom), lin1(nmom) (lin1 needs lin): the linear weights of cheb_linear_weights for E^0 and E^1 on the same
   !> interval, window and kernel.
   module subroutine chebyshev_exchange_weights(this, h, kernel, lin, lin1)
      class(green), intent(inout) :: this
      complex(rp), allocatable, intent(out) :: h(:, :)
      real(rp), allocatable, intent(out) :: kernel(:)
      complex(rp), allocatable, intent(out), optional :: lin(:), lin1(:)
      real(rp) :: a, b, emin_win, emax_win
      integer :: nmom

      nmom = intersite_moment_count(this)
      call this%recursion%resolve_chebyshev_window(emin_win, emax_win)
      a = (emax_win - emin_win)/(2 - 0.3_rp)
      b = (emax_win + emin_win)/2.0_rp
      allocate (h(0:2*nmom - 2, 0:2), kernel(nmom))
      call cheb_bilinear_weights(nmom, a, b, this%en%energy_min, this%en%fermi, kB_simpson_f*this%en%temperature, h)
      call cheb_kernel_vector(this, nmom, kernel)
      if (present(lin)) then
         allocate (lin(nmom))
         if (present(lin1)) allocate (lin1(nmom))
         call this%chebyshev_linear_weights(nmom, lin, lin1)
      end if
   end subroutine chebyshev_exchange_weights

   !> Linear weights c(nmom) (and c1, c2 for E^1, E^2) of cheb_linear_weights at en%fermi, with the window map and the damping kernel
   !> of chebyshev_exchange_weights.
   module subroutine chebyshev_linear_weights(this, nmom, c, c1, c2)
      class(green), intent(inout) :: this
      integer, intent(in) :: nmom
      complex(rp), intent(out) :: c(nmom)
      complex(rp), intent(out), optional :: c1(nmom), c2(nmom)
      real(rp) :: a, b, emin_win, emax_win, kernel(nmom)

      call this%recursion%resolve_chebyshev_window(emin_win, emax_win)
      a = (emax_win - emin_win)/(2 - 0.3_rp)
      b = (emax_win + emin_win)/2.0_rp
      call cheb_kernel_vector(this, nmom, kernel)
      call cheb_linear_weights(kernel, a, b, this%en%energy_min, this%en%fermi, kB_simpson_f*this%en%temperature, c, c1, c2)
   end subroutine chebyshev_linear_weights

   !> Raw energy moments mp(nb, nb, 0:2, n_local) = int f(E) E^p n_ii(E) dE (about E = 0, electrons*Ry^p) of the on-site density matrix of
   !> every local recursion site (mu_n(:, :, :, 1:n_local)). With G = sum_m g_m mu_m, mu_m = <i|T_m|i> Hermitian, n = -(1/2 pi i)(G - G^H)
   !> = -(1/pi) sum_m Im(g_m) mu_m, so mp(:, :, p, i) = -(1/pi) sum_m Im(c^(p)_m) mu_m with c^(p) of chebyshev_linear_weights.
   !> herm is the largest |mp - mp^H| (zero if mu_m is Hermitian at every order); mp is returned as computed, not Hermitized.
   module subroutine chebyshev_onsite_energy_moments(this, mp, herm)
      class(green), intent(inout) :: this
      complex(rp), allocatable, intent(out) :: mp(:, :, :, :)
      real(rp), intent(out) :: herm
      complex(rp), allocatable :: c(:, :)
      integer :: nmom, n_local, n, p, m

      nmom = size(this%recursion%mu_n, 3)
      n_local = end_atom - start_atom + 1
      allocate (c(nmom, 0:2), mp(size(this%recursion%mu_n, 1), size(this%recursion%mu_n, 2), 0:2, max(n_local, 0)))
      call this%chebyshev_linear_weights(nmom, c(:, 0), c(:, 1), c(:, 2))
      mp = (0.0_rp, 0.0_rp)
      herm = 0.0_rp
      do n = 1, n_local
         do p = 0, 2
            do m = 1, nmom
               mp(:, :, p, n) = mp(:, :, p, n) - aimag(c(m, p))/pi*this%recursion%mu_n(:, :, m, n)
            end do
            herm = max(herm, maxval(abs(mp(:, :, p, n) - conjg(transpose(mp(:, :, p, n))))))
         end do
      end do
   end subroutine chebyshev_onsite_energy_moments

   !> Raw moments mij, mji(nb, nb, nmom) of G_ij and G_ji of global pair ia_glob (one_shot, pairs or self-pair source).
   subroutine pair_moment_sources(this, ia_glob, mij, mji)
      class(green), intent(in) :: this
      integer, intent(in) :: ia_glob
      complex(rp), allocatable, intent(out) :: mij(:, :, :), mji(:, :, :)
      complex(rp), allocatable :: m4(:, :, :, :)
      integer :: ia, ja_temp
      logical :: self_pair

      ia = g2l_map(ia_glob)
      ja_temp = (ia - 1)*4 + 1
      self_pair = this%lattice%ijpair(ia_glob, 1) == this%lattice%ijpair(ia_glob, 2)
      if (this%control%intersite_mode == 'one_shot') then
         mji = this%recursion%os_mu(:, :, :, this%recursion%req%pair_entry(ia))
         mij = mji
         if (.not. self_pair) mij = moments_dagger(mji)
      else if (self_pair) then
         mji = this%recursion%mu_n(:, :, :, ja_temp)
         mij = mji
      else
         m4 = this%recursion%mu_n(:, :, :, ja_temp:ja_temp + 3)
         mij = 0.5_rp*((m4(:, :, :, 1) - m4(:, :, :, 2)) + (1.0_rp/i_unit)*(m4(:, :, :, 3) - m4(:, :, :, 4)))
         mji = 0.5_rp*((m4(:, :, :, 1) - m4(:, :, :, 2)) - (1.0_rp/i_unit)*(m4(:, :, :, 3) - m4(:, :, :, 4)))
      end if
   end subroutine pair_moment_sources

   !> Raw moments of G_ij and G_ji of global pair ia_glob in gij/gji(:, :, 1:nmom, slot), from the same moments as
   !> fill_intersite_gf_pair, then the Pauli and two-index families of the slot (linear in gij, gji, so also moments).
   module subroutine fill_intersite_moments_pair(this, ia_glob, slot)
      class(green), intent(inout) :: this
      integer, intent(in) :: ia_glob, slot
      complex(rp), allocatable :: mij(:, :, :), mji(:, :, :)
      integer :: nmom

      nmom = intersite_moment_count(this)
      call pair_moment_sources(this, ia_glob, mij, mji)
      this%gij(:, :, :, slot) = 0.0_rp
      this%gji(:, :, :, slot) = 0.0_rp
      this%gij(:, :, 1:nmom, slot) = mij
      this%gji(:, :, 1:nmom, slot) = mji
      call fill_slot_families(this, slot)
   end subroutine fill_intersite_moments_pair

   !> Hamiltonian blocks H_ij = a mu_1 + b mu_0 and H_ji of global pair ia_glob from its raw moments mu_m = <i|T_m(Ht)|j>, and
   !> the window map H = a Ht + b of chebyshev_exchange_weights.
   module subroutine pair_hamiltonian_block(this, ia_glob, hij, hji, a, b)
      class(green), intent(inout) :: this
      integer, intent(in) :: ia_glob
      complex(rp), intent(out) :: hij(:, :), hji(:, :)
      real(rp), intent(out) :: a, b
      complex(rp), allocatable :: mij(:, :, :), mji(:, :, :)
      real(rp) :: emin_win, emax_win

      call this%recursion%resolve_chebyshev_window(emin_win, emax_win)
      a = (emax_win - emin_win)/(2 - 0.3_rp)
      b = (emax_win + emin_win)/2.0_rp
      call pair_moment_sources(this, ia_glob, mij, mji)
      hij = a*mij(:, :, 2) + b*mij(:, :, 1)
      hji = a*mji(:, :, 2) + b*mji(:, :, 1)
   end subroutine pair_hamiltonian_block

   !> G_ij, G_ji and their second energy derivatives (nb x nb) at E = en%fermi of global pair ia_glob, summed
   !> from the moments with the window and kernel of chebyshev_exchange_weights.
   module subroutine chebyshev_point_pair(this, ia_glob, gij0, gji0, gij2, gji2)
      class(green), intent(inout) :: this
      integer, intent(in) :: ia_glob
      complex(rp), intent(out), dimension(:, :) :: gij0, gji0, gij2, gji2
      complex(rp), allocatable :: mij(:, :, :), mji(:, :, :), g0(:), g2(:)
      real(rp), allocatable :: kernel(:)
      real(rp) :: a, b, emin_win, emax_win
      integer :: nmom, m

      nmom = intersite_moment_count(this)
      call this%recursion%resolve_chebyshev_window(emin_win, emax_win)
      a = (emax_win - emin_win)/(2 - 0.3_rp)
      b = (emax_win + emin_win)/2.0_rp
      allocate (kernel(nmom), g0(nmom), g2(nmom))
      call cheb_kernel_vector(this, nmom, kernel)
      call cheb_green_point_weights(kernel, a, b, this%en%fermi, g0, g2)
      call pair_moment_sources(this, ia_glob, mij, mji)
      gij0 = czero; gji0 = czero; gij2 = czero; gji2 = czero
      do m = 1, nmom
         gij0 = gij0 + g0(m)*mij(:, :, m)
         gji0 = gji0 + g0(m)*mji(:, :, m)
         gij2 = gij2 + g2(m)*mij(:, :, m)
         gji2 = gji2 + g2(m)*mji(:, :, m)
      end do
   end subroutine chebyshev_point_pair

   !> gij, gji of slot_t = bilinear transform (weights column hp) of the raw moments in slot_raw, and their families.
   module subroutine transform_intersite_moments(this, slot_raw, slot_t, hp, kernel)
      class(green), intent(inout) :: this
      integer, intent(in) :: slot_raw, slot_t
      complex(rp), intent(in) :: hp(0:)
      real(rp), intent(in) :: kernel(:)
      integer :: nmom

      nmom = size(kernel)
      this%gij(:, :, :, slot_t) = 0.0_rp
      this%gji(:, :, :, slot_t) = 0.0_rp
      call cheb_bilinear_transform(hp, kernel, this%gij(:, :, 1:nmom, slot_raw), this%gij(:, :, 1:nmom, slot_t))
      call cheb_bilinear_transform(hp, kernel, this%gji(:, :, 1:nmom, slot_raw), this%gji(:, :, 1:nmom, slot_t))
      call fill_slot_families(this, slot_t)
   end subroutine transform_intersite_moments

   !> Pauli and two-index families of the slot from its gij, gji.
   subroutine fill_slot_families(this, slot)
      class(green), intent(inout) :: this
      integer, intent(in) :: slot

      this%ginmag(:, :, :, slot) = 0.0_rp; this%giz(:, :, :, slot) = 0.0_rp; this%giy(:, :, :, slot) = 0.0_rp; this%gix(:, :, :, slot) = 0.0_rp
      this%gjnmag(:, :, :, slot) = 0.0_rp; this%gjx(:, :, :, slot) = 0.0_rp; this%gjy(:, :, :, slot) = 0.0_rp; this%gjz(:, :, :, slot) = 0.0_rp
      call pauli_pair(this, slot)
      call this%fill_intersite_gf_twoindex_pair(slot)
   end subroutine fill_slot_families

   !---------------------------------------------------------------------------
   !> @brief GPU drop-in for chebyshev_green_ij (intersite Chebyshev GF).
   !---------------------------------------------------------------------------
   module subroutine chebyshev_green_ij_gpu(this, istart)
      class(green), intent(inout) :: this
      integer, intent(in) :: istart
      complex(rp), allocatable :: mu_local(:, :, :, :), g0_local(:, :, :, :)
      real(rp), allocatable :: damp(:)
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

      call cheb_damping(this, n_mom, damp)
      call gpu_backend%chebyshev_dos(mu_local, this%en%ene(1:nv), a, b, g0_local, damp)
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
      real(rp), allocatable :: damp(:)
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

      call cheb_damping(this, size(this%recursion%mu_n, 3), damp)
      do n = 1, 4 ! Loop on the number of on-site GFs to calculate the inter-site GFs
         call cheb_green_complex(this%recursion%mu_n(:, :, :, n + istart - 1), &
            nb, size(this%recursion%mu_n, 3), cmplx(this%en%ene(fermi_point), 0.0_rp, rp) + eta, a, b, g_ef(:, :, n), damp)
      end do  ! End loop on n
   end subroutine chebyshev_green_ij_eta

   !> Contour G_ij and G_ji (nb,nb,size(z)) of local pair ij_loc at the complex energies z from the one-shot moments.
   !> Self pair: G_ij = G_ji = G_ii.
   module subroutine chebyshev_green_ij_eta_oneshot(this, ij_loc, z, gij, gji)
      class(green), intent(inout) :: this
      integer, intent(in) :: ij_loc
      complex(rp), intent(in) :: z(:)
      complex(rp), intent(out) :: gij(:, :, :), gji(:, :, :)
      complex(rp), allocatable :: mu_local(:, :, :, :), g0_local(:, :, :, :), mij(:, :, :)
      real(rp), allocatable :: damp(:)
      real(rp) :: a, b, emin_win, emax_win
      integer :: ient, iz, n_mom
      logical :: self_pair
      type(rsrec_cuda_backend), pointer :: gpu_backend

      ient = this%recursion%req%pair_entry(ij_loc)
      self_pair = this%lattice%ijpair(start_atom + ij_loc - 1, 1) == this%lattice%ijpair(start_atom + ij_loc - 1, 2)
      n_mom = size(this%recursion%os_mu, 3)
      call this%recursion%resolve_chebyshev_window(emin_win, emax_win)
      a = (emax_win - emin_win)/(2 - 0.3_rp)
      b = (emax_win + emin_win)/2.0_rp
      call cheb_damping(this, n_mom, damp)

      if (this%control%gpu_plugin .and. rsrec_cuda_plugin_compiled()) then
         allocate (mu_local(nb, nb, n_mom, 2), g0_local(nb, nb, size(z), 2))
         mu_local(:, :, :, 1) = moments_dagger(this%recursion%os_mu(:, :, :, ient))
         mu_local(:, :, :, 2) = this%recursion%os_mu(:, :, :, ient)
         gpu_backend => get_gpu_context(1, nb, 1, 1, 1)
         call gpu_backend%set_precision(this%recursion%gpu_precision_code())
         call gpu_backend%chebyshev_gf_eta(mu_local, real(z(1), rp), z - real(z(1), rp), a, b, g0_local, damp)
         gij = g0_local(:, :, :, 1)
         gji = g0_local(:, :, :, 2)
         deallocate (mu_local, g0_local)
      else
         if (.not. self_pair) mij = moments_dagger(this%recursion%os_mu(:, :, :, ient))
         do iz = 1, size(z)
            call cheb_green_complex(this%recursion%os_mu(:, :, :, ient), nb, n_mom, z(iz), a, b, gji(:, :, iz), damp)
            if (.not. self_pair) call cheb_green_complex(mij, nb, n_mom, z(iz), a, b, gij(:, :, iz), damp)
         end do
      end if
      if (self_pair) gij = gji
   end subroutine chebyshev_green_ij_eta_oneshot

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

   !> Moment-space SCF energy points. With t_m the diagonal on-site moments summed over atoms and spins,
   !> N(E) = -Im sum_m t_m c_m(E)/pi, c_m(E) = sum_q w_q g0_m(e_q), (e_q, w_q) = cheb_energy_nodes(e_f = E).
   !> If solve_fermi, en%fermi is the root of N(E) = valence (bisection on [energy_min, upper window edge], 1e-12 Ry).
   !> The nodes at en%fermi replace en%ene (weights in en%node_w) and the on-site G is rebuilt on them.
   module subroutine chebyshev_scf_nodes(this, t, valence, solve_fermi)
      class(green), intent(inout) :: this
      complex(rp), intent(in) :: t(:)
      real(rp), intent(in) :: valence
      logical, intent(in) :: solve_fermi
      real(rp), allocatable :: kernel(:), e(:), w(:)
      real(rp) :: a, b, emin_win, emax_win, lo, hi, mid, kbt
      integer :: n, it

      n = size(t)
      call this%recursion%resolve_chebyshev_window(emin_win, emax_win)
      a = (emax_win - emin_win)/(2 - 0.3_rp)
      b = (emax_win + emin_win)/2.0_rp
      kbt = kB_simpson_f*this%en%temperature
      allocate (kernel(n))
      call cheb_kernel_vector(this, n, kernel)

      if (solve_fermi) then
         lo = this%en%energy_min
         hi = emax_win
         if (electrons(hi) < valence) call g_logger%fatal('chebyshev_scf_nodes: valence not reached below the window top', &
                                                          __FILE__, __LINE__)
         do it = 1, 100
            if (hi - lo <= 1.0e-12_rp) exit
            mid = 0.5_rp*(lo + hi)
            if (electrons(mid) < valence) then
               lo = mid
            else
               hi = mid
            end if
         end do
         this%en%fermi = 0.5_rp*(lo + hi)
      end if

      call cheb_energy_nodes(n, a, b, this%en%energy_min, this%en%fermi, kbt, e, w)
      call move_alloc(e, this%en%ene)
      call move_alloc(w, this%en%node_w)
      call this%chebyshev_dos_dispatch()

   contains

      function electrons(ef) result(ne)
         real(rp), intent(in) :: ef
         real(rp) :: ne
         complex(rp) :: c(n)

         call cheb_linear_weights(kernel, a, b, this%en%energy_min, ef, kbt, c)
         ne = -aimag(sum(t*c))/pi
      end function electrons
   end subroutine chebyshev_scf_nodes

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
      real(rp), allocatable :: damp(:)
      complex(rp), intent(in), optional :: eta
      integer, intent(in), optional :: fermi_point
      complex(rp), dimension(nb, nb, atoms_per_process), intent(inout), optional :: g_ef
      logical, intent(in), optional :: use_gpu
      ! Local variables
      real(rp), dimension(size(this%recursion%mu_n, 3)) :: kernel
      real(rp), dimension(size(this%en%ene)) :: wscale
      real(rp) :: a, b, emin_win, emax_win
      complex(rp) :: exp_factor
      integer :: ie, i, l, m, n, nv, n_glob, n_local
      complex(rp), allocatable :: mug(:, :, :)
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
      nv = size(this%en%ene)

      if (.not. present(eta) .and. do_gpu) then
         gpu_backend => get_gpu_context(1, nb, 1, 1, 1)
         call gpu_backend%set_precision(this%recursion%gpu_precision_code())

         call cheb_damping(this, size(this%recursion%mu_n, 3), damp)
         call gpu_backend%chebyshev_dos(this%recursion%mu_n(:, :, :, start_atom:end_atom), &
                                         this%en%ene(1:nv), a, b,                         &
                                         this%g0(:, :, 1:nv, start_atom:end_atom), damp)
         call enforce_chebyshev_spectral_positivity(this, nv)
         return
      end if

      if (.not. present(eta)) then
         n_local = end_atom - start_atom + 1
         if (n_local <= 0) return

         allocate (g0_local(nb, nb, nv, n_local))

         ! The local atoms are mu_n(:,:,:,1:n_local) (g2l_map is compact): no copy of the moments.
         ! fast_dp (default) and legacy reconstruct in double precision; the other backends keep fp32.
         ! damp stays unallocated for the Jackson kernel, which the reconstruction builds inline.
         call cheb_damping(this, size(this%recursion%mu_n, 3), damp)
         if (trim(this%control%cheb_backend) == 'fast_dp' .or. trim(this%control%cheb_backend) == 'legacy') then
            call cheb_green_dp(this%recursion%mu_n(:, :, :, 1:n_local), nb, size(this%recursion%mu_n, 3), n_local, &
                               this%en%ene(1:nv), nv, a, b, g0_local, damp)
         else
            call cheb_green_fast(this%recursion%mu_n(:, :, :, 1:n_local), nb, size(this%recursion%mu_n, 3), n_local, &
                                 this%en%ene(1:nv), nv, a, b, g0_local, damp)
         end if

         do n_glob = start_atom, end_atom
            n = g2l_map(n_glob)
            this%g0(:, :, 1:nv, n) = g0_local(:, :, :, n)
         end do

         deallocate (g0_local)
         call enforce_chebyshev_spectral_positivity(this, nv)
         return
      end if

      wscale(:) = (this%en%ene(:) - b)/a

      ! Calculating the damping kernel selected by control%cheb_kernel (Jackson by default)
      call cheb_kernel_vector(this, size(this%recursion%mu_n, 3), kernel)

      allocate (mug(nb, nb, size(this%recursion%mu_n, 3)))
      do n_glob = start_atom, end_atom ! Loop on self-consistent atoms
         n = g2l_map(n_glob)
         ! Multiply the moments with the kernel (diagonal elements only)
         mug = (0.0_rp, 0.0_rp)
         do i = 1, nb
            mug(i, i, :) = this%recursion%mu_n(i, i, :, n)*kernel(:)
         end do

         mug(:, :, 2:size(kernel)) = mug(:, :, 2:size(kernel))*2.0_rp

         !do i=1, size(kernel)
         !  write(400+n, *) i, sum(this%recursion%mu_n(n, i, 1:18, 1:18))
         !end do

         ! Calculate the density of states
         !$omp parallel do default(shared) private(ie, i, exp_factor, l,m)
         do ie = fermi_point, fermi_point
            do i = 1, size(kernel)
               exp_factor = -i_unit*exp(-i_unit*(i - 1)*acos(((this%en%ene(ie) + eta) - b)/a))
               do l = 1, nb
                  do m = 1, nb
                     g_ef(l, m, n) = g_ef(l, m, n) + mug(l, m, i)*exp_factor
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
      deallocate (mug)
   end subroutine chebyshev_green_core

   !> Damping kernel g_n (n = 0..n-1) selected by control%cheb_kernel; always filled.
   subroutine cheb_kernel_vector(this, n, kernel)
      class(green), intent(in) :: this
      integer, intent(in) :: n
      real(rp), intent(out) :: kernel(n)

      select case (trim(this%control%cheb_kernel))
      case ('lorentz')
         call lorentz_kernel(n, kernel, this%control%cheb_kernel_lambda)
      case ('none')
         kernel = 1.0_rp
      case default
         call jackson_kernel(n, kernel)
      end select
   end subroutine cheb_kernel_vector

   !> Optional damping argument of the reconstruction kernels: left unallocated for the default
   !> Jackson kernel (the kernels then build it inline, bit-identical to before).
   module subroutine cheb_damping(this, n, damping)
      class(green), intent(in) :: this
      integer, intent(in) :: n
      real(rp), allocatable, intent(out) :: damping(:)

      if (trim(this%control%cheb_kernel) == 'jackson') return
      allocate (damping(n))
      call cheb_kernel_vector(this, n, damping)
   end subroutine cheb_damping

   !> Restore the local spectral positivity required by the real-axis DOS
   !> route.  The finite Jackson expansion is non-negative in exact arithmetic (the
   !> Lorentz kernel and cheb_kernel = 'none' do not guarantee it),
   !> but a low-precision moment or GEMM path can make a tiny negative DOS.
   !> Work in each orbital's spin block so the total DOS, spin-resolved DOS,
   !> and the charge/moment integrations all consume the same physical density.
   !> Eta-shifted and intersite Green functions are intentionally excluded.
   !> With cheb_lp the clipped weight is taken out again: each atom, orbital and
   !> spin channel is rescaled so that its integrated weight (node-weighted on moment-space SCF nodes) equals the weight
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
      real(rp) :: s_ud, max_rescale, wt

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
            ! Integration weight of the point: the moment-space SCF nodes carry quadrature weights, the mesh is uniform
            wt = 1.0_rp
            if (allocated(this%en%node_w)) wt = this%en%node_w(ie)
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
               q_uu(orbital) = q_uu(orbital) + wt*a_uu
               q_dd(orbital) = q_dd(orbital) + wt*a_dd
               if (.not. corrected) cycle

               dq_uu(orbital) = dq_uu(orbital) + wt*a_uu - wt*a_uu_before
               dq_dd(orbital) = dq_dd(orbital) + wt*a_dd - wt*a_dd_before
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
