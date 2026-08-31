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
      call gpu_backend%set_precision(0)  ! 0=fp32 (fast), 1=fp64 (validation)

      call gpu_backend%block_dos(this%recursion%a_b(:, :, :, istart:istart + na - 1), &
                                 this%recursion%b2_b(:, :, :, istart:istart + na - 1), &
                                 a_inf_d, b_inf_d, this%en%ene(1:nv), eta, &
                                 this%control%sym_term, this%g0(:, :, 1:nv, 1:na))
   end subroutine block_green_ij_gpu

   module subroutine calculate_intersite_gf_twoindex(this)
      use mpi_mod
   use basis_mod, only: nb, norb, spin_off
      implicit none
      class(green), intent(inout) :: this
      integer :: i, j, k, l1, l2, k0, j0, ia_glob, ia

      this%G00ij = 0.0d0; this%G01ij = 0.0d0; this%G00ji = 0.0d0; this%G01ji = 0.0d0;
      this%Gx1ij = 0.0d0; this%Gy1ij = 0.0d0; this%Gz1ij = 0.0d0; this%Gx0ij = 0.0d0; this%Gy0ij = 0.0d0; this%Gz0ij = 0.0d0;
      this%Gx1ji = 0.0d0; this%Gy1ji = 0.0d0; this%Gz1ji = 0.0d0; this%Gx0ji = 0.0d0; this%Gy0ji = 0.0d0; this%Gz0ji = 0.0d0

      do ia_glob = start_atom, end_atom
         ia = g2l_map(ia_glob)
         do j = 1, norb
            do k = 1, norb
               l1 = int((k - 0.9)**0.5)
               l2 = int((j - 0.9)**0.5)
               k0 = l1*(l1 + 1) + 1
               j0 = l2*(l2 + 1) + 1
               this%G00ij(k, j, :, IA) = this%G00ij(k, j, :, IA) + 0.5d0*(this%Ginmag(k, j, :, IA) + (-1)**(k + j)*this%Gjnmag(2*j0 - j, 2*k0 - k, :, IA))
               this%G01ij(k, j, :, IA) = this%G01ij(k, j, :, IA) + 0.5d0*(this%Ginmag(k, j, :, IA) - (-1)**(k + j)*this%Gjnmag(2*j0 - j, 2*k0 - k, :, IA))
               this%G00ji(k, j, :, IA) = this%G00ji(k, j, :, IA) + 0.5d0*(this%Gjnmag(k, j, :, IA) + (-1)**(k + j)*this%Ginmag(2*j0 - j, 2*k0 - k, :, IA))
               this%G01ji(k, j, :, IA) = this%G01ji(k, j, :, IA) + 0.5d0*(this%Gjnmag(k, j, :, IA) - (-1)**(k + j)*this%Ginmag(2*j0 - j, 2*k0 - k, :, IA))
               this%Gx1ij(k, j, :, IA) = this%Gx1ij(k, j, :, IA) + 0.5d0*(this%Gix(k, j, :, IA) - (-1)**(k + j)*this%Gjx(2*j0 - j, 2*k0 - k, :, IA))
               this%Gy1ij(k, j, :, IA) = this%Gy1ij(k, j, :, IA) + 0.5d0*(this%Giy(k, j, :, IA) - (-1)**(k + j)*this%Gjy(2*j0 - j, 2*k0 - k, :, IA))
               this%Gz1ij(k, j, :, IA) = this%Gz1ij(k, j, :, IA) + 0.5d0*(this%Giz(k, j, :, IA) - (-1)**(k + j)*this%Gjz(2*j0 - j, 2*k0 - k, :, IA))
               this%Gx0ij(k, j, :, IA) = this%Gx0ij(k, j, :, IA) + 0.5d0*(this%Gix(k, j, :, IA) + (-1)**(k + j)*this%Gjx(2*j0 - j, 2*k0 - k, :, IA))
               this%Gy0ij(k, j, :, IA) = this%Gy0ij(k, j, :, IA) + 0.5d0*(this%Giy(k, j, :, IA) + (-1)**(k + j)*this%Gjy(2*j0 - j, 2*k0 - k, :, IA))
               this%Gz0ij(k, j, :, IA) = this%Gz0ij(k, j, :, IA) + 0.5d0*(this%Giz(k, j, :, IA) + (-1)**(k + j)*this%Gjz(2*j0 - j, 2*k0 - k, :, IA))
               this%Gx1ji(k, j, :, IA) = this%Gx1ji(k, j, :, IA) + 0.5d0*(this%Gjx(k, j, :, IA) - (-1)**(k + j)*this%Gix(2*j0 - j, 2*k0 - k, :, IA))
               this%Gy1ji(k, j, :, IA) = this%Gy1ji(k, j, :, IA) + 0.5d0*(this%Gjy(k, j, :, IA) - (-1)**(k + j)*this%Giy(2*j0 - j, 2*k0 - k, :, IA))
               this%Gz1ji(k, j, :, IA) = this%Gz1ji(k, j, :, IA) + 0.5d0*(this%Gjz(k, j, :, IA) - (-1)**(k + j)*this%Giz(2*j0 - j, 2*k0 - k, :, IA))
               this%Gx0ji(k, j, :, IA) = this%Gx0ji(k, j, :, IA) + 0.5d0*(this%Gjx(k, j, :, IA) + (-1)**(k + j)*this%Gix(2*j0 - j, 2*k0 - k, :, IA))
               this%Gy0ji(k, j, :, IA) = this%Gy0ji(k, j, :, IA) + 0.5d0*(this%Gjy(k, j, :, IA) + (-1)**(k + j)*this%Giy(2*j0 - j, 2*k0 - k, :, IA))
               this%Gz0ji(k, j, :, IA) = this%Gz0ji(k, j, :, IA) + 0.5d0*(this%Gjz(k, j, :, IA) + (-1)**(k + j)*this%Giz(2*j0 - j, 2*k0 - k, :, IA))
            end do
         end do
      end do
   end subroutine calculate_intersite_gf_twoindex

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
      integer :: ia, ja_temp, j, i, ia_glob, fermi_point
      complex(rp), dimension(nb, nb, 4) :: g0_ef
      complex(rp) :: eta
      complex(rp), dimension(64, nb, nb, 4) :: y
      real(rp) :: res, t
      real(rp), dimension(64) :: x, w

      if (eta_mode) then
         if (this%control%gpu_plugin .and. rsrec_cuda_plugin_compiled() .and. &
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
            do i = 1, norb
               do j = 1, norb
                  this%Ginmag_eta(:, j, i, iA) = this%Ginmag_eta(:, j, i, iA) + (this%gij_eta(:, j, i, ia) + this%gij_eta(:, j +spin_off, i +spin_off, ia))*0.5d0
                  this%Giz_eta(:, j, i, iA) = this%Giz_eta(:, j, i, iA) + 0.5d0*(this%gij_eta(:, j, i, ia) - this%gij_eta(:, j +spin_off, i +spin_off, ia))
                  this%Giy_eta(:, j, i, iA) = this%Giy_eta(:, j, i, iA) + 0.5d0*(i_unit*this%gij_eta(:, j, i +spin_off, ia) - i_unit*this%gij_eta(:, j +spin_off, i, ia))
                  this%Gix_eta(:, j, i, iA) = this%Gix_eta(:, j, i, iA) + 0.5d0*(this%gij_eta(:, j, i +spin_off, ia) + this%gij_eta(:, j +spin_off, i, ia))
                  this%Gjnmag_eta(:, j, i, iA) = this%Gjnmag_eta(:, j, i, iA) + (this%gji_eta(:, j, i, ia) + this%gji_eta(:, j +spin_off, i +spin_off, ia))*0.5d0
                  this%Gjz_eta(:, j, i, iA) = this%Gjz_eta(:, j, i, iA) + 0.5d0*(this%gji_eta(:, j, i, ia) - this%gji_eta(:, j +spin_off, i +spin_off, ia))
                  this%Gjy_eta(:, j, i, iA) = this%Gjy_eta(:, j, i, iA) + 0.5d0*(i_unit*this%gji_eta(:, j, i +spin_off, ia) - i_unit*this%gji_eta(:, j +spin_off, i, ia))
                  this%Gjx_eta(:, j, i, iA) = this%Gjx_eta(:, j, i, iA) + 0.5d0*(this%gji_eta(:, j, i +spin_off, ia) + this%gji_eta(:, j +spin_off, i, ia))
               end do
            end do
         end do
         return
      end if

      this%gij = 0.0d0; this%gji = 0.0d0; this%ginmag = 0.0d0; this%giz = 0.0d0; this%giy = 0.0d0; this%gix = 0.0d0
      this%gjnmag = 0.0d0; this%gjx = 0.0d0; this%gjy = 0.0d0; this%gjz = 0.0d0

      if (this%control%recur == 'block') call this%recursion%zsqr()

      do ia_glob = start_atom, end_atom
         ia = g2l_map(ia_glob)
         ja_temp = (ia - 1)*4 + 1
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
         if (this%lattice%ijpair(ia, 1) .eq. this%lattice%ijpair(ia, 2)) then
            this%gij(:, :, :, ia) = this%g0(:, :, :, 1)
            this%gji(:, :, :, ia) = this%g0(:, :, :, 1)
         else
            this%gij(:, :, :, ia) = this%g0(:, :, :, 1) - this%g0(:, :, :, 2) + (1.0d0/i_unit*this%g0(:, :, :, 3) - 1.0d0/i_unit*this%g0(:, :, :, 4))
            this%gji(:, :, :, ia) = this%g0(:, :, :, 1) - this%g0(:, :, :, 2) - (1.0d0/i_unit*this%g0(:, :, :, 3) - 1.0d0/i_unit*this%g0(:, :, :, 4))
            this%gij(:, :, :, ia) = this%gij(:, :, :, ia)*0.5d0
            this%gji(:, :, :, ia) = this%gji(:, :, :, ia)*0.5d0
         end if
         do i = 1, norb
            do j = 1, norb
               this%Ginmag(j, i, :, iA) = this%Ginmag(j, i, :, iA) + (this%gij(j, i, :, ia) + this%gij(j +spin_off, i +spin_off, :, ia))*0.5d0
               this%Giz(j, i, :, iA) = this%Giz(j, i, :, iA) + 0.5d0*(this%gij(j, i, :, ia) - this%gij(j +spin_off, i +spin_off, :, ia))
               this%Giy(j, i, :, iA) = this%Giy(j, i, :, iA) + 0.5d0*(i_unit*this%gij(j, i +spin_off, :, ia) - i_unit*this%gij(j +spin_off, i, :, ia))
               this%Gix(j, i, :, iA) = this%Gix(j, i, :, iA) + 0.5d0*(this%gij(j, i +spin_off, :, ia) + this%gij(j +spin_off, i, :, ia))
               this%Gjnmag(j, i, :, iA) = this%Gjnmag(j, i, :, iA) + (this%gji(j, i, :, ia) + this%gji(j +spin_off, i +spin_off, :, ia))*0.5d0
               this%Gjz(j, i, :, iA) = this%Gjz(j, i, :, iA) + 0.5d0*(this%gji(j, i, :, ia) - this%gji(j +spin_off, i +spin_off, :, ia))
               this%Gjy(j, i, :, iA) = this%Gjy(j, i, :, iA) + 0.5d0*(i_unit*this%gji(j, i +spin_off, :, ia) - i_unit*this%gji(j +spin_off, i, :, ia))
               this%Gjx(j, i, :, iA) = this%Gjx(j, i, :, iA) + 0.5d0*(this%gji(j, i +spin_off, :, ia) + this%gji(j +spin_off, i, :, ia))
            end do
         end do
      end do
   end subroutine calculate_intersite_gf_core

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
      call gpu_backend%set_precision(0)  ! 0=fp32 (fast), 1=fp64 (validation)

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
            call this%recursion%resolve_chebyshev_window(emin_win, emax_win)
            a = (emax_win - emin_win)/(2 - 0.3_rp)
            b = (emax_win + emin_win)/2.0_rp
            call gpu_backend%chebyshev_gf_eta(this%recursion%mu_n(:, :, :, 1:natoms), &
                                              this%en%ene(fermi_point), eta_list, a, b, g0_ef_all)
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
         do i = 1, norb
            do j = 1, norb
               this%Ginmag_eta(:, j, i, iA) = this%Ginmag_eta(:, j, i, iA) + (this%gij_eta(:, j, i, ia) + this%gij_eta(:, j +spin_off, i +spin_off, ia))*0.5d0
               this%Giz_eta(:, j, i, iA) = this%Giz_eta(:, j, i, iA) + 0.5d0*(this%gij_eta(:, j, i, ia) - this%gij_eta(:, j +spin_off, i +spin_off, ia))
               this%Giy_eta(:, j, i, iA) = this%Giy_eta(:, j, i, iA) + 0.5d0*(i_unit*this%gij_eta(:, j, i +spin_off, ia) - i_unit*this%gij_eta(:, j +spin_off, i, ia))
               this%Gix_eta(:, j, i, iA) = this%Gix_eta(:, j, i, iA) + 0.5d0*(this%gij_eta(:, j, i +spin_off, ia) + this%gij_eta(:, j +spin_off, i, ia))
               !
               this%Gjnmag_eta(:, j, i, iA) = this%Gjnmag_eta(:, j, i, iA) + (this%gji_eta(:, j, i, ia) + this%gji_eta(:, j +spin_off, i +spin_off, ia))*0.5d0
               this%Gjz_eta(:, j, i, iA) = this%Gjz_eta(:, j, i, iA) + 0.5d0*(this%gji_eta(:, j, i, ia) - this%gji_eta(:, j +spin_off, i +spin_off, ia))
               this%Gjy_eta(:, j, i, iA) = this%Gjy_eta(:, j, i, iA) + 0.5d0*(i_unit*this%gji_eta(:, j, i +spin_off, ia) - i_unit*this%gji_eta(:, j +spin_off, i, ia))
               this%Gjx_eta(:, j, i, iA) = this%Gjx_eta(:, j, i, iA) + 0.5d0*(this%gji_eta(:, j, i +spin_off, ia) + this%gji_eta(:, j +spin_off, i, ia))
            end do
         end do
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

   module subroutine block_green_gpu(this)
      class(green), intent(inout) :: this

      ! Deprecated wrapper retained for type-bound callers; use block_green_core.
      call block_green_core(this)
   end subroutine block_green_gpu

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
      call gpu_backend%set_precision(0)  ! 0=fp32 (fast), 1=fp64 (validation)

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
