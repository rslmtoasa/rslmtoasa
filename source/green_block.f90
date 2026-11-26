!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! SUBMODULE: green_block
!
!> Recursive Green-function implementations moved from green_mod.
!------------------------------------------------------------------------------

submodule (green_mod) green_block

   use mpi_mod
   use rsrec_cuda_plugin_mod, only: rsrec_cuda_backend, get_gpu_context, rsrec_cuda_plugin_compiled
   use logger_mod, only: g_logger
   implicit none

contains

   module subroutine block_green_ij_eta(this, istart, eta, fermi_point, g_ef)
      implicit none
      class(green), intent(inout) :: this
      integer, intent(in) :: istart
      complex(rp), intent(in) :: eta
      complex(rp), dimension(nb, nb, 4), intent(inout) :: g_ef
      integer, intent(in) :: fermi_point
      !
      integer :: ll, n, nw, ldim, na
      real(rp), dimension(4) :: a_inf0, b_inf0
      real(rp), dimension(nb, nb, 4) :: a_inf, b_inf
      complex(rp), dimension(nb, nb, this%en%channels_ldos + 10, 4) :: dum_g_ef
      !
      !
      ! Definitions so it is not necessary to change the code
      ll = this%control%lld
      ldim = nb
      na = 4
      nw = 10*ll

      call this%recursion%get_terminf(this%recursion%a_b(:, :, :, istart), this%recursion%b2_b(:, :, :, istart), na, ll, ldim, nw, &
                                      a_inf, b_inf, a_inf0, b_inf0)

      do n = 1, NA
         call this%bgreen(dum_g_ef(:, :, :, n), n + istart - 1, fermi_point, 1, a_inf(:, :, n), b_inf(:, :, n), eta)
         g_ef(:, :, n) = dum_g_ef(:, :, fermi_point, n)
      end do

   end subroutine block_green_ij_eta

   !> @brief One-shot block Haydock G_ji(z), G_ij(z) of local pair ij_loc for the energies z(:).
   !> @details G_ji(z) = sum_m S_m(j) X_m(z) with X the first block column of (z - T)^-1 (T = blocks 1..lld-1,
   !>          terminal self-energy Q(z) of the legacy contour route subtracted from the last diagonal block),
   !>          G_ij(z) = [G_ji(conj z)]^H. Self pair: G_ij = G_ji = G_ii.
   !>          Terminator inputs and b are built once; the z loop is threaded.
   !>          Needs zsqr to have turned the stored B^2 into B (b(:,:,m) = B_{m+1} = os_b2(:,:,m+1)).
   module subroutine block_green_ij_oneshot(this, ij_loc, z, gij, gji)
      use block_tridiag_mod, only: block_tridiag_first_column
      class(green), intent(inout) :: this
      integer, intent(in) :: ij_loc
      complex(rp), intent(in) :: z(:)
      complex(rp), intent(out) :: gij(:, :, :), gji(:, :, :)
      complex(rp), allocatable :: b(:, :, :), x(:, :, :)
      complex(rp), dimension(nb, nb) :: sigma
      real(rp), dimension(nb, nb, 1) :: a_inf, b_inf
      real(rp), dimension(1) :: a_inf0, b_inf0
      integer :: m, nblk, iseed, ient, iz, nw
      logical :: self_pair

      nblk = this%control%lld - 1
      iseed = this%recursion%req%pair_seed(ij_loc)
      ient = this%recursion%req%pair_entry(ij_loc)
      allocate (b(nb, nb, nblk))
      b(:, :, 1:nblk - 1) = this%recursion%os_b2(:, :, 2:nblk, iseed)
      b(:, :, nblk) = (0.0_rp, 0.0_rp)
      self_pair = this%lattice%ijpair(start_atom + ij_loc - 1, 1) == this%lattice%ijpair(start_atom + ij_loc - 1, 2)

      ! Terminator from the seed's coefficients exactly as in block_green_ij_eta; Q is the self-energy on block lld-1
      nw = 10*this%control%lld
      call this%recursion%get_terminf(this%recursion%os_a(:, :, :, iseed), this%recursion%os_b2(:, :, :, iseed), 1, &
                                      this%control%lld, nb, nw, a_inf, b_inf, a_inf0, b_inf0)

      !$omp parallel do default(shared) private(iz, m, x, sigma)
      do iz = 1, size(z)
         allocate (x(nb, nb, nblk))
         call bgreen_terminal_q(this%control%sym_term, .true., z(iz), a_inf(:, :, 1), b_inf(:, :, 1), sigma)
         call block_tridiag_first_column(this%recursion%os_a(:, :, 1:nblk, iseed), b, z(iz), x, sigma)
         gji(:, :, iz) = (0.0_rp, 0.0_rp)
         do m = 1, nblk
            gji(:, :, iz) = gji(:, :, iz) + matmul(this%recursion%os_s(:, :, m, ient), x(:, :, m))
         end do
         if (self_pair) then
            gij(:, :, iz) = gji(:, :, iz)
         else
            ! Advanced terminator: the legacy-mode Q keeps Re(sqrt) fixed at Re z, so Q(conj z) /= conj(Q(z)); use conj(Q(z))
            call block_tridiag_first_column(this%recursion%os_a(:, :, 1:nblk, iseed), b, conjg(z(iz)), x, conjg(sigma))
            gij(:, :, iz) = (0.0_rp, 0.0_rp)
            do m = 1, nblk
               gij(:, :, iz) = gij(:, :, iz) + matmul(this%recursion%os_s(:, :, m, ient), x(:, :, m))
            end do
            gij(:, :, iz) = conjg(transpose(gij(:, :, iz)))
         end if
         deallocate (x)
      end do
      !$omp end parallel do
   end subroutine block_green_ij_oneshot

   module subroutine block_green_ij(this, istart)
      implicit none
      class(green), intent(inout) :: this
      integer, intent(in) :: istart
      !
      complex(rp) :: eta
      integer :: ll, n, nw, ldim, na
      real(rp), dimension(4) :: a_inf0, b_inf0
      real(rp), dimension(nb, nb, 4) :: a_inf, b_inf
      !
      !
      ! Definitions so it is not necessary to change the code
      ll = this%control%lld
      ldim = nb
      na = 4
      eta = (0.0d0, 0.0d0)
      nw = 10*ll

      call this%recursion%get_terminf(this%recursion%a_b(:, :, :, istart), this%recursion%b2_b(:, :, :, istart), na, ll, ldim, nw, &
                                      a_inf, b_inf, a_inf0, b_inf0)

      do n = 1, NA
         ! Initializing Shift and Scaling for rigid band shift
         !Dfac_mat=(1.0d0,0.0d0)
         !Cshi_mat=(0.0d0,0.0d0)
         call this%bgreen(this%g0(:, :, :, n), n + istart - 1, 1, this%en%channels_ldos + 10, a_inf(:, :, n), b_inf(:, :, n), eta)

      end do
   end subroutine block_green_ij

   module subroutine block_green_ij_gpu(this, istart)
      implicit none
      class(green), intent(inout) :: this
      integer, intent(in) :: istart
      complex(rp) :: eta
      integer :: ll, n, nw, ldim, na, nv, i
      real(rp), dimension(4) :: a_inf0, b_inf0
      real(rp), dimension(nb, nb, 4) :: a_inf, b_inf
      real(rp), dimension(nb, 4) :: a_inf_d, b_inf_d
      type(rsrec_cuda_backend), pointer :: gpu_backend

      if (.not. rsrec_cuda_plugin_compiled()) then
         call this%block_green_ij(istart)
         return
      end if

      ll = this%control%lld
      ldim = nb
      na = 4
      eta = (0.0d0, 0.0d0)
      nw = 10*ll
      nv = this%en%channels_ldos + 10

      call this%recursion%get_terminf(this%recursion%a_b(:, :, :, istart), this%recursion%b2_b(:, :, :, istart), &
                                      na, ll, ldim, nw, a_inf, b_inf, a_inf0, b_inf0)
      do n = 1, na
         do i = 1, nb
            a_inf_d(i, n) = a_inf(i, i, n)
            b_inf_d(i, n) = b_inf(i, i, n)
         end do
      end do

      gpu_backend => get_gpu_context(1, nb, 1, 1, 1)
      call gpu_backend%set_precision(this%recursion%gpu_precision_code())

      call gpu_backend%block_dos(this%recursion%a_b(:, :, :, istart:istart + na - 1), &
                                 this%recursion%b2_b(:, :, :, istart:istart + na - 1), &
                                 a_inf_d, b_inf_d, this%en%ene(1:nv), eta, &
                                 this%control%sym_term, this%g0(:, :, 1:nv, 1:na))
   end subroutine block_green_ij_gpu

   module subroutine calculate_intersite_gf_twoindex(this)
      use mpi_mod
      class(green), intent(inout) :: this
      integer :: ia_glob

      do ia_glob = start_atom, end_atom
         call this%fill_intersite_gf_twoindex_pair(g2l_map(ia_glob))
      end do
   end subroutine calculate_intersite_gf_twoindex

   !> Two-index Green-function families of the pair held in slot `slot` of Ginmag, Gjnmag, Gi{x,y,z}, Gj{x,y,z}.
   module subroutine fill_intersite_gf_twoindex_pair(this, slot)
      use basis_mod, only: norb
      class(green), intent(inout) :: this
      integer, intent(in) :: slot
      integer :: j, k, l1, ie, s, pm(norb)

      this%G00ij(:, :, :, slot) = 0.0d0; this%G01ij(:, :, :, slot) = 0.0d0; this%G00ji(:, :, :, slot) = 0.0d0; this%G01ji(:, :, :, slot) = 0.0d0;
      this%Gx1ij(:, :, :, slot) = 0.0d0; this%Gy1ij(:, :, :, slot) = 0.0d0; this%Gz1ij(:, :, :, slot) = 0.0d0; this%Gx0ij(:, :, :, slot) = 0.0d0; this%Gy0ij(:, :, :, slot) = 0.0d0; this%Gz0ij(:, :, :, slot) = 0.0d0;
      this%Gx1ji(:, :, :, slot) = 0.0d0; this%Gy1ji(:, :, :, slot) = 0.0d0; this%Gz1ji(:, :, :, slot) = 0.0d0; this%Gx0ji(:, :, :, slot) = 0.0d0; this%Gy0ji(:, :, :, slot) = 0.0d0; this%Gz0ji(:, :, :, slot) = 0.0d0

      ! m -> -m partner of each orbital (l, m): pm(k) = 2 k0 - k, k0 the m = 0 index of shell l
      do k = 1, norb
         l1 = int((k - 0.9)**0.5)
         pm(k) = 2*(l1*(l1 + 1) + 1) - k
      end do
      do ie = 1, size(this%G00ij, 3)
         do j = 1, norb
            do k = 1, norb
               s = (-1)**(k + j)
               this%G00ij(k, j, ie, slot) = this%G00ij(k, j, ie, slot) + 0.5d0*(this%Ginmag(k, j, ie, slot) + s*this%Gjnmag(pm(j), pm(k), ie, slot))
               this%G01ij(k, j, ie, slot) = this%G01ij(k, j, ie, slot) + 0.5d0*(this%Ginmag(k, j, ie, slot) - s*this%Gjnmag(pm(j), pm(k), ie, slot))
               this%G00ji(k, j, ie, slot) = this%G00ji(k, j, ie, slot) + 0.5d0*(this%Gjnmag(k, j, ie, slot) + s*this%Ginmag(pm(j), pm(k), ie, slot))
               this%G01ji(k, j, ie, slot) = this%G01ji(k, j, ie, slot) + 0.5d0*(this%Gjnmag(k, j, ie, slot) - s*this%Ginmag(pm(j), pm(k), ie, slot))
               this%Gx1ij(k, j, ie, slot) = this%Gx1ij(k, j, ie, slot) + 0.5d0*(this%Gix(k, j, ie, slot) - s*this%Gjx(pm(j), pm(k), ie, slot))
               this%Gy1ij(k, j, ie, slot) = this%Gy1ij(k, j, ie, slot) + 0.5d0*(this%Giy(k, j, ie, slot) - s*this%Gjy(pm(j), pm(k), ie, slot))
               this%Gz1ij(k, j, ie, slot) = this%Gz1ij(k, j, ie, slot) + 0.5d0*(this%Giz(k, j, ie, slot) - s*this%Gjz(pm(j), pm(k), ie, slot))
               this%Gx0ij(k, j, ie, slot) = this%Gx0ij(k, j, ie, slot) + 0.5d0*(this%Gix(k, j, ie, slot) + s*this%Gjx(pm(j), pm(k), ie, slot))
               this%Gy0ij(k, j, ie, slot) = this%Gy0ij(k, j, ie, slot) + 0.5d0*(this%Giy(k, j, ie, slot) + s*this%Gjy(pm(j), pm(k), ie, slot))
               this%Gz0ij(k, j, ie, slot) = this%Gz0ij(k, j, ie, slot) + 0.5d0*(this%Giz(k, j, ie, slot) + s*this%Gjz(pm(j), pm(k), ie, slot))
               this%Gx1ji(k, j, ie, slot) = this%Gx1ji(k, j, ie, slot) + 0.5d0*(this%Gjx(k, j, ie, slot) - s*this%Gix(pm(j), pm(k), ie, slot))
               this%Gy1ji(k, j, ie, slot) = this%Gy1ji(k, j, ie, slot) + 0.5d0*(this%Gjy(k, j, ie, slot) - s*this%Giy(pm(j), pm(k), ie, slot))
               this%Gz1ji(k, j, ie, slot) = this%Gz1ji(k, j, ie, slot) + 0.5d0*(this%Gjz(k, j, ie, slot) - s*this%Giz(pm(j), pm(k), ie, slot))
               this%Gx0ji(k, j, ie, slot) = this%Gx0ji(k, j, ie, slot) + 0.5d0*(this%Gjx(k, j, ie, slot) + s*this%Gix(pm(j), pm(k), ie, slot))
               this%Gy0ji(k, j, ie, slot) = this%Gy0ji(k, j, ie, slot) + 0.5d0*(this%Gjy(k, j, ie, slot) + s*this%Giy(pm(j), pm(k), ie, slot))
               this%Gz0ji(k, j, ie, slot) = this%Gz0ji(k, j, ie, slot) + 0.5d0*(this%Gjz(k, j, ie, slot) + s*this%Giz(pm(j), pm(k), ie, slot))
            end do
         end do
      end do
   end subroutine fill_intersite_gf_twoindex_pair

   module subroutine calculate_intersite_gf(this)
      class(green), intent(inout) :: this

      ! Deprecated wrapper retained for type-bound callers; use calculate_intersite_gf_core.
      call calculate_intersite_gf_core(this, .false.)
   end subroutine calculate_intersite_gf

   module subroutine calculate_intersite_gf_eta(this)
      class(green), intent(inout) :: this

      ! Deprecated wrapper retained for type-bound callers; use calculate_intersite_gf_core.
      call calculate_intersite_gf_core(this, .true.)
   end subroutine calculate_intersite_gf_eta

   module subroutine calculate_intersite_gf_core(this, eta_mode)
      use mpi_mod
      implicit none
      class(green), intent(inout) :: this
      logical, intent(in) :: eta_mode
      integer :: ia, ja_temp, i, ia_glob, fermi_point
      complex(rp), dimension(nb, nb, 4) :: g0_ef
      complex(rp), dimension(nb, nb, 64, 2) :: g_one
      complex(rp), dimension(64) :: z
      complex(rp) :: eta
      complex(rp), dimension(64, nb, nb, 4) :: y
      real(rp) :: res, t
      real(rp), dimension(64) :: x, w

      if (eta_mode) then
         ! One-shot data is not in the four-slot form: it is reconstructed per pair below (block on the CPU, Chebyshev on the GPU if requested).
         if (this%control%gpu_plugin .and. rsrec_cuda_plugin_compiled() .and. this%control%intersite_mode /= 'one_shot' .and. &
             (this%control%recur == 'block' .or. this%control%recur == 'chebyshev')) then
            call this%calculate_intersite_gf_eta_gpu()
            return
         end if

         call gauss_legendre(64, 0.0_rp, 1.0_rp, x, w)
         call this%en%e_mesh()

         do i = 1, this%en%channels_ldos + 10
            if ((this%en%ene(i) - this%en%fermi) .le. 0.000001d0) fermi_point = i
         end do

         write (*, *) this%en%fermi, fermi_point
         this%gij_eta = 0.0d0; this%gji_eta = 0.0d0; this%ginmag_eta = 0.0d0; this%giz_eta = 0.0d0; this%giy_eta = 0.0d0; this%gix_eta = 0.0d0
         this%gjnmag_eta = 0.0d0; this%gjx_eta = 0.0d0; this%gjy_eta = 0.0d0; this%gjz_eta = 0.0d0

         if (this%control%recur == 'block') call this%recursion%zsqr()

         do ia_glob = start_atom, end_atom
            ia = g2l_map(ia_glob)
            ja_temp = (ia - 1)*4 + 1
            if (this%control%intersite_mode == 'one_shot') then
               do i = 1, 64
                  z(i) = cmplx(this%en%ene(fermi_point), 0.0_rp, rp) + cmplx(0.0_rp, (1 - x(i))/x(i))
               end do
               if (this%control%recur == 'block') then
                  call this%block_green_ij_oneshot(ia, z, g_one(:, :, :, 1), g_one(:, :, :, 2))
               else
                  call this%chebyshev_green_ij_eta_oneshot(ia, z, g_one(:, :, :, 1), g_one(:, :, :, 2))
               end if
               ! The contour combination of the pairs route returns 0.5 G_ii for i = j (the real-mesh route special-cases it to G_ii)
               if (this%lattice%ijpair(ia_glob, 1) == this%lattice%ijpair(ia_glob, 2)) g_one = 0.5_rp*g_one
               do i = 1, 64
                  this%gij_eta(i, :, :, ia) = g_one(:, :, i, 1)
                  this%gji_eta(i, :, :, ia) = g_one(:, :, i, 2)
               end do
            else
               y = (0.0_rp, 0.0_rp)
               do i = 1, 64
                  eta = (0.0_rp, 0.0_rp)
                  g0_ef = (0.0_rp, 0.0_rp)
                  res = (1 - x(i))/x(i)
                  eta = cmplx(0.0_rp, res)
                  select case (this%control%recur)
                  case ('block')
                     call this%block_green_ij_eta(ja_temp, eta, fermi_point, g0_ef)
                  case ('chebyshev')
                     call this%chebyshev_green_ij_eta(ja_temp, eta, fermi_point, g0_ef)
                  end select
                  y(i, :, :, :) = g0_ef(:, :, :)
               end do
               this%gij_eta(:, :, :, ia) = y(:, :, :, 1) - y(:, :, :, 2) + (1.0d0/i_unit*y(:, :, :, 3) - 1.0d0/i_unit*y(:, :, :, 4))
               this%gji_eta(:, :, :, ia) = y(:, :, :, 1) - y(:, :, :, 2) - (1.0d0/i_unit*y(:, :, :, 3) - 1.0d0/i_unit*y(:, :, :, 4))

               this%gij_eta(:, :, :, ia) = this%gij_eta(:, :, :, ia)*0.5d0
               this%gji_eta(:, :, :, ia) = this%gji_eta(:, :, :, ia)*0.5d0
            end if
            call pauli_eta_pair(this, ia)
         end do
         return
      end if

      call this%prepare_intersite_gf()
      do ia_glob = start_atom, end_atom
         call this%fill_intersite_gf_pair(ia_glob, g2l_map(ia_glob))
      end do
   end subroutine calculate_intersite_gf_core

   module subroutine prepare_intersite_gf(this)
      class(green), intent(inout) :: this

      if (this%control%recur == 'block') call this%recursion%zsqr()
   end subroutine prepare_intersite_gf

   !> Real-mesh G_ij and G_ji of global pair ia_glob and their Pauli families, written to slot `slot` of the
   !> green arrays. The recursion data stays indexed by the local pair. Needs prepare_intersite_gf once before.
   module subroutine fill_intersite_gf_pair(this, ia_glob, slot)
      use mpi_mod
      class(green), intent(inout) :: this
      integer, intent(in) :: ia_glob, slot
      integer :: ia, ja_temp, j, i

      ia = g2l_map(ia_glob)
      ja_temp = (ia - 1)*4 + 1
      this%gij(:, :, :, slot) = 0.0d0; this%gji(:, :, :, slot) = 0.0d0
      this%ginmag(:, :, :, slot) = 0.0d0; this%giz(:, :, :, slot) = 0.0d0; this%giy(:, :, :, slot) = 0.0d0; this%gix(:, :, :, slot) = 0.0d0
      this%gjnmag(:, :, :, slot) = 0.0d0; this%gjx(:, :, :, slot) = 0.0d0; this%gjy(:, :, :, slot) = 0.0d0; this%gjz(:, :, :, slot) = 0.0d0
      if (this%control%intersite_mode == 'one_shot') then
         if (this%control%recur == 'block') then
            call this%block_green_ij_oneshot(ia, cmplx(this%en%ene(1:this%en%channels_ldos + 10), 0.0_rp, rp), &
                                             this%gij(:, :, :, slot), this%gji(:, :, :, slot))
         else
            call this%chebyshev_green_ij_oneshot(ia, this%gij(:, :, :, slot), this%gji(:, :, :, slot))
         end if
         call pauli_pair(this, slot)
         return
      end if
      select case (this%control%recur)
      case ('block')
         if (this%control%gpu_plugin .and. rsrec_cuda_plugin_compiled()) then
            call this%block_green_ij_gpu(ja_temp)
         else
            call this%block_green_ij(ja_temp)
         end if
      case ('chebyshev')
         if (this%control%gpu_plugin .and. rsrec_cuda_plugin_compiled()) then
            call this%chebyshev_green_ij_gpu(ja_temp)
         else
            call this%chebyshev_green_ij(ja_temp)
         end if
      end select
      if (this%lattice%ijpair(ia_glob, 1) .eq. this%lattice%ijpair(ia_glob, 2)) then
         this%gij(:, :, :, slot) = this%g0(:, :, :, 1)
         this%gji(:, :, :, slot) = this%g0(:, :, :, 1)
      else
         this%gij(:, :, :, slot) = this%g0(:, :, :, 1) - this%g0(:, :, :, 2) + (1.0d0/i_unit*this%g0(:, :, :, 3) - 1.0d0/i_unit*this%g0(:, :, :, 4))
         this%gji(:, :, :, slot) = this%g0(:, :, :, 1) - this%g0(:, :, :, 2) - (1.0d0/i_unit*this%g0(:, :, :, 3) - 1.0d0/i_unit*this%g0(:, :, :, 4))
         this%gij(:, :, :, slot) = this%gij(:, :, :, slot)*0.5d0
         this%gji(:, :, :, slot) = this%gji(:, :, :, slot)*0.5d0
      end if
      call pauli_pair(this, slot)
   end subroutine fill_intersite_gf_pair

   !> Pauli families Ginmag, Giz, Giy, Gix (i side) and Gjnmag, Gjz, Gjy, Gjx (j side) of the pair in slot `slot`
   !> from gij, gji; the target slices must be zero on entry.
   module subroutine pauli_pair(this, slot)
      class(green), intent(inout) :: this
      integer, intent(in) :: slot
      integer :: i, j

      do i = 1, norb
         do j = 1, norb
            this%Ginmag(j, i, :, slot) = this%Ginmag(j, i, :, slot) + (this%gij(j, i, :, slot) + this%gij(j +spin_off, i +spin_off, :, slot))*0.5d0
            this%Giz(j, i, :, slot) = this%Giz(j, i, :, slot) + 0.5d0*(this%gij(j, i, :, slot) - this%gij(j +spin_off, i +spin_off, :, slot))
            this%Giy(j, i, :, slot) = this%Giy(j, i, :, slot) + 0.5d0*(i_unit*this%gij(j, i +spin_off, :, slot) - i_unit*this%gij(j +spin_off, i, :, slot))
            this%Gix(j, i, :, slot) = this%Gix(j, i, :, slot) + 0.5d0*(this%gij(j, i +spin_off, :, slot) + this%gij(j +spin_off, i, :, slot))
            this%Gjnmag(j, i, :, slot) = this%Gjnmag(j, i, :, slot) + (this%gji(j, i, :, slot) + this%gji(j +spin_off, i +spin_off, :, slot))*0.5d0
            this%Gjz(j, i, :, slot) = this%Gjz(j, i, :, slot) + 0.5d0*(this%gji(j, i, :, slot) - this%gji(j +spin_off, i +spin_off, :, slot))
            this%Gjy(j, i, :, slot) = this%Gjy(j, i, :, slot) + 0.5d0*(i_unit*this%gji(j, i +spin_off, :, slot) - i_unit*this%gji(j +spin_off, i, :, slot))
            this%Gjx(j, i, :, slot) = this%Gjx(j, i, :, slot) + 0.5d0*(this%gji(j, i +spin_off, :, slot) + this%gji(j +spin_off, i, :, slot))
         end do
      end do
   end subroutine pauli_pair

   !> Contour Pauli families (accumulated into Ginmag_eta, ..., Gjx_eta) of local pair ia from gij_eta, gji_eta.
   subroutine pauli_eta_pair(this, ia)
      class(green), intent(inout) :: this
      integer, intent(in) :: ia
      integer :: i, j

      do i = 1, norb
         do j = 1, norb
            this%Ginmag_eta(:, j, i, ia) = this%Ginmag_eta(:, j, i, ia) + (this%gij_eta(:, j, i, ia) + this%gij_eta(:, j +spin_off, i +spin_off, ia))*0.5d0
            this%Giz_eta(:, j, i, ia) = this%Giz_eta(:, j, i, ia) + 0.5d0*(this%gij_eta(:, j, i, ia) - this%gij_eta(:, j +spin_off, i +spin_off, ia))
            this%Giy_eta(:, j, i, ia) = this%Giy_eta(:, j, i, ia) + 0.5d0*(i_unit*this%gij_eta(:, j, i +spin_off, ia) - i_unit*this%gij_eta(:, j +spin_off, i, ia))
            this%Gix_eta(:, j, i, ia) = this%Gix_eta(:, j, i, ia) + 0.5d0*(this%gij_eta(:, j, i +spin_off, ia) + this%gij_eta(:, j +spin_off, i, ia))
            !
            this%Gjnmag_eta(:, j, i, ia) = this%Gjnmag_eta(:, j, i, ia) + (this%gji_eta(:, j, i, ia) + this%gji_eta(:, j +spin_off, i +spin_off, ia))*0.5d0
            this%Gjz_eta(:, j, i, ia) = this%Gjz_eta(:, j, i, ia) + 0.5d0*(this%gji_eta(:, j, i, ia) - this%gji_eta(:, j +spin_off, i +spin_off, ia))
            this%Gjy_eta(:, j, i, ia) = this%Gjy_eta(:, j, i, ia) + 0.5d0*(i_unit*this%gji_eta(:, j, i +spin_off, ia) - i_unit*this%gji_eta(:, j +spin_off, i, ia))
            this%Gjx_eta(:, j, i, ia) = this%Gjx_eta(:, j, i, ia) + 0.5d0*(this%gji_eta(:, j, i +spin_off, ia) + this%gji_eta(:, j +spin_off, i, ia))
         end do
      end do
   end subroutine pauli_eta_pair

   module subroutine calculate_intersite_gf_eta_gpu(this)
      use mpi_mod
      implicit none
      class(green), intent(inout) :: this
      integer :: ia, ja_temp, j, i, k, c, fermi_point, ia_glob, n_pairs, natoms, ll, ldim, nw
      complex(rp), dimension(64, nb, nb, 4) :: y
      real(rp), dimension(64) :: x, w
      complex(rp), dimension(64) :: eta_list
      real(rp), allocatable :: a_inf(:, :, :), b_inf(:, :, :), a_inf0(:), b_inf0(:)
      real(rp), allocatable :: a_inf_d(:, :), b_inf_d(:, :)
      complex(rp), allocatable :: g0_ef_all(:, :, :, :)
      type(rsrec_cuda_backend), pointer :: gpu_backend

      ll = this%control%lld
      ldim = nb
      nw = 10*ll

      call gauss_legendre(64, 0.0_rp, 1.0_rp, x, w)
      call this%en%e_mesh()
      fermi_point = 1
      do i = 1, this%en%channels_ldos + 10
         if ((this%en%ene(i) - this%en%fermi) .le. 0.000001d0) fermi_point = i
      end do

      this%gij_eta = 0.0d0; this%gji_eta = 0.0d0; this%ginmag_eta = 0.0d0; this%giz_eta = 0.0d0; this%giy_eta = 0.0d0; this%gix_eta = 0.0d0
      this%gjnmag_eta = 0.0d0; this%gjx_eta = 0.0d0; this%gjy_eta = 0.0d0; this%gjz_eta = 0.0d0

      if (this%control%recur == 'block') call this%recursion%zsqr()

      ! Gauss-Legendre eta contour (same map as the CPU path: eta = i*(1-x)/x)
      do k = 1, 64
         eta_list(k) = cmplx(0.0_rp, (1.0_rp - x(k))/x(k), kind=rp)
      end do

      ! All local pairs share the contiguous coefficient layout (4 combos/pair).
      n_pairs = end_atom - start_atom + 1
      natoms = n_pairs*4
      allocate (g0_ef_all(nb, nb, 64, natoms))

      gpu_backend => get_gpu_context(1, nb, 1, 1, 1)
      call gpu_backend%set_precision(this%recursion%gpu_precision_code())

      ! One device call: Gij at the Fermi energy for all (pair x combo) atoms
      ! and all 64 contour etas -> g0_ef_all(nb,nb,64,natoms).
      select case (this%control%recur)
      case ('block')
         allocate (a_inf(nb, nb, natoms), b_inf(nb, nb, natoms), a_inf0(natoms), b_inf0(natoms))
         allocate (a_inf_d(nb, natoms), b_inf_d(nb, natoms))
         call this%recursion%get_terminf(this%recursion%a_b, this%recursion%b2_b, &
                                         natoms, ll, ldim, nw, a_inf, b_inf, a_inf0, b_inf0)
         do k = 1, natoms
            do i = 1, nb
               a_inf_d(i, k) = a_inf(i, i, k)
               b_inf_d(i, k) = b_inf(i, i, k)
            end do
         end do
         call gpu_backend%block_gf_eta(this%recursion%a_b(:, :, :, 1:natoms), &
                                       this%recursion%b2_b(:, :, :, 1:natoms), &
                                       a_inf_d, b_inf_d, this%en%ene(fermi_point), &
                                       eta_list, this%control%sym_term, g0_ef_all)
         deallocate (a_inf, b_inf, a_inf0, b_inf0, a_inf_d, b_inf_d)
      case ('chebyshev')
         block
            real(rp) :: a, b, emin_win, emax_win
            real(rp), allocatable :: damp(:)
            call this%recursion%resolve_chebyshev_window(emin_win, emax_win)
            a = (emax_win - emin_win)/(2 - 0.3_rp)
            b = (emax_win + emin_win)/2.0_rp
            call cheb_damping(this, size(this%recursion%mu_n, 3), damp)
            call gpu_backend%chebyshev_gf_eta(this%recursion%mu_n(:, :, :, 1:natoms), &
                                              this%en%ene(fermi_point), eta_list, a, b, g0_ef_all, damp)
         end block
      end select

      ! Recombine per pair (identical algebra to the CPU loop).
      do ia_glob = start_atom, end_atom
         ia = g2l_map(ia_glob)
         ja_temp = (ia - 1)*4 + 1
         do c = 1, 4
            do k = 1, 64
               y(k, :, :, c) = g0_ef_all(:, :, k, ja_temp + c - 1)
            end do
         end do
         this%gij_eta(:, :, :, ia) = (y(:, :, :, 1) - y(:, :, :, 2) + (1.0d0/i_unit*y(:, :, :, 3) - 1.0d0/i_unit*y(:, :, :, 4)))*0.5d0
         this%gji_eta(:, :, :, ia) = (y(:, :, :, 1) - y(:, :, :, 2) - (1.0d0/i_unit*y(:, :, :, 3) - 1.0d0/i_unit*y(:, :, :, 4)))*0.5d0
         call pauli_eta_pair(this, ia)
      end do

      deallocate (g0_ef_all)
   end subroutine calculate_intersite_gf_eta_gpu

   module subroutine block_green_eta(this, eta, fermi_point, g_ef)
      use mpi_mod
      class(green), intent(inout) :: this
      complex(rp), intent(in) :: eta
      integer, intent(in) :: fermi_point
      complex(rp), dimension(nb, nb, atoms_per_process), intent(inout) :: g_ef

      ! Deprecated wrapper retained for type-bound callers; use block_green_core.
      call block_green_core(this, eta, fermi_point, g_ef)
   end subroutine block_green_eta

   module subroutine block_green(this)
      class(green), intent(inout) :: this

      ! Deprecated wrapper retained for type-bound callers; use block_green_core.
      call block_green_core(this)
   end subroutine block_green

   module subroutine block_green_core(this, eta, fermi_point, g_ef)
      use mpi_mod
      implicit none
      class(green), intent(inout) :: this
      complex(rp), intent(in), optional :: eta
      integer, intent(in), optional :: fermi_point
      complex(rp), dimension(nb, nb, atoms_per_process), intent(inout), optional :: g_ef
      integer :: nw, ll, ldim, nv, n_loc, n, n_glob, k, i
      complex(rp) :: eta0
      real(rp), dimension(this%lattice%nrec) :: a_inf0, b_inf0
      real(rp), dimension(nb, nb, this%lattice%nrec) :: a_inf, b_inf
      complex(rp), allocatable :: ab_loc(:, :, :, :), b2_loc(:, :, :, :), g0_loc(:, :, :, :)
      complex(rp), dimension(nb, nb, this%en%channels_ldos + 10, atoms_per_process) :: dum_g_ef
      real(rp), allocatable :: ainf_d(:, :), binf_d(:, :)
      type(rsrec_cuda_backend), pointer :: gpu_backend

      ll = this%control%lld
      ldim = nb
      eta0 = (0.0d0, 0.0d0)
      if (present(eta)) eta0 = eta
      nw = 10*ll
      nv = this%en%channels_ldos + 10

      call this%recursion%get_terminf(this%recursion%a_b, this%recursion%b2_b, &
                                      atoms_per_process, ll, ldim, nw, a_inf, b_inf, a_inf0, b_inf0)

      if (present(eta)) then
         do n_glob = start_atom, end_atom
            n = g2l_map(n_glob)
            call this%bgreen(dum_g_ef(:, :, :, n), n, fermi_point, 1, a_inf(:, :, n), b_inf(:, :, n), eta0)
            g_ef(:, :, n) = dum_g_ef(:, :, fermi_point, n)
         end do
         return
      end if

      if (.not. this%control%gpu_plugin .or. .not. rsrec_cuda_plugin_compiled()) then
         if (this%control%gpu_plugin .and. .not. rsrec_cuda_plugin_compiled()) then
            call g_logger%warning('GPU block Green reconstruction requested, but this executable '// &
               'was built without ENABLE_CUDA_PLUGIN. Falling back to the CPU reconstruction.', &
               __FILE__, __LINE__)
         end if
         do n_glob = start_atom, end_atom
            n = g2l_map(n_glob)
            call this%bgreen(this%g0(:, :, :, n), n, 1, nv, a_inf(:, :, n), b_inf(:, :, n), eta0)
         end do
         return
      end if

      ! Gather this rank’s local atoms into contiguous batch buffers (mirrors the
      ! n = g2l_map(n_glob) indexing of block_green’s bgreen loop).
      n_loc = end_atom - start_atom + 1
      allocate (ab_loc(nb, nb, ll, n_loc), b2_loc(nb, nb, ll, n_loc), g0_loc(nb, nb, nv, n_loc))
      allocate (ainf_d(nb, n_loc), binf_d(nb, n_loc))
      do n_glob = start_atom, end_atom
         n = g2l_map(n_glob)
         k = n_glob - start_atom + 1
         ab_loc(:, :, :, k) = this%recursion%a_b(:, :, :, n)
         b2_loc(:, :, :, k) = this%recursion%b2_b(:, :, :, n)
         do i = 1, nb
            ainf_d(i, k) = a_inf(i, i, n)
            binf_d(i, k) = b_inf(i, i, n)
         end do
      end do

      gpu_backend => get_gpu_context(1, nb, 1, 1, 1)
      call gpu_backend%set_precision(this%recursion%gpu_precision_code())

      call gpu_backend%block_dos(ab_loc, b2_loc, ainf_d, binf_d, &
                                 this%en%ene(1:nv), eta0, this%control%sym_term, g0_loc)

      ! Scatter back to the local g0 slots.
      do n_glob = start_atom, end_atom
         n = g2l_map(n_glob)
         k = n_glob - start_atom + 1
         this%g0(:, :, 1:nv, n) = g0_loc(:, :, :, k)
      end do

      deallocate (ab_loc, b2_loc, g0_loc, ainf_d, binf_d)
   end subroutine block_green_core

end submodule green_block
