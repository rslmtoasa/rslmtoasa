!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Symbolic Atom
!
!> @author
!> Angela Klautau
!> Ramon Cardias
!> Lucas P. Campagna
!> S. Frota-Pessôa
!> Pascoal R. Peduto
!> Anders Bergman
!> S. B. Legoas
!> H. M. Petrilli
!> Ivan P. Miranda
!
! DESCRIPTION:
!> Module to hold symbolic_atom's parameters
!------------------------------------------------------------------------------

module symbolic_atom_mod

   use element_mod, only: element
   use potential_mod, only: potential
   use string_mod, only: sl, replace, endswith, fmt
   use precision_mod, only: rp
   use potential_mod
   use math_mod
   use logger_mod, only: g_logger
   implicit none

   private

   ! public functions
   public :: array_of_symbolic_atoms
   public :: print_state, print_state_full, print_state_formatted
   public :: save_state
   public :: load_state

   type, public :: symbolic_atom
      !> Mixture occupation in self-consistent calculation. Default: 0.01.
      !>
      !> Mixture occupation in self-consistent calculation.
      !>
      !> Default: 0.01.
      real(rp) :: mix

      !> Spin-(up/down) occupation in self-consistent calculation. Default: 0.05.
      !>
      !> Spin-(up/down) occupation in self-consistent calculation.
      !>
      !> Default: 0.05.
      real(rp) :: mixmag

      !> Specify the number of rigid band calculation used.
      !>
      !> Specify the number of rigid band calculation used
      integer :: rb
      ! TODO
      ! From common_cnstr
      real(rp), dimension(3) :: mag_cfield, mag_cfield_diff
      real(rp) :: chg_cfield, chg_cfield_diff
      real(rp) :: chg_con_val, mag_con_val
      ! TODO
      real(rp) :: a
      ! TODO (Charge?)
      real(rp) :: dq, qc, qv

      type(element) :: element
      type(potential) :: potential
   contains
      procedure :: build_from_file
      procedure :: restore_to_default
      procedure :: print_state => symbolic_atom_print_state
      procedure :: print_state_full => symbolic_atom_print_state_full
      procedure :: print_state_formatted => symbolic_atom_print_state_formatted
      procedure :: mesh_grid_size
      procedure :: B
      procedure :: rho0
      procedure :: predls
      procedure :: build_pot
      procedure :: d_matrix
      procedure :: disp_matrix
      procedure :: p_matrix
      procedure :: transform_pmatrix
      procedure :: udisp_matrix
   end type symbolic_atom

   interface symbolic_atom
      procedure :: constructor
   end interface symbolic_atom

   interface array_of_symbolic_atoms
      procedure :: array_of_symbolic_atoms_from_memory
      procedure :: array_of_symbolic_atoms_from_file
   end interface array_of_symbolic_atoms

   integer, parameter :: min_mesh_grid_size = 25

contains

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Constructor
   !
   !> @param[in] symbolic_atom Namelist file in database
   !> @param[in] database Directory to database files with 'symbolic_atom' namelist
   !> @return type(control)
   !---------------------------------------------------------------------------
   function constructor(label, database, reload) result(obj)
      type(symbolic_atom) :: obj
      character(len=*), intent(in) :: label
      character(len=*), intent(in), optional :: database
      logical, optional, intent(in) :: reload
      logical :: reload_

      !reload_ = merge(reload, .True., present(reload))
      if (present(reload)) then
         reload_ = reload
      else
         reload_ = .true.
      end if
      call obj%restore_to_default()
      obj%element = element(label, database, reload_)
      obj%potential = potential(label, database, reload_)
   end function constructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Read parameters from input file
   !
   !> @param[in] fname Namelist file
   !---------------------------------------------------------------------------
   subroutine build_from_file(this, fname)
      class(symbolic_atom), intent(out) :: this
      character(len=*), intent(in) :: fname
      call this%element%build_from_file(fname)
      call this%potential%build_from_file(fname)
   end subroutine build_from_file

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Reset all members to default
   !---------------------------------------------------------------------------
   subroutine restore_to_default(this)
      implicit none
      class(symbolic_atom), intent(out) :: this
      this%mix = -1
      this%mixmag = -1
      this%rb = -1
      this%a = 0.02D0

      call this%element%restore_to_default()
      call this%potential%restore_to_default()
   end subroutine restore_to_default

   subroutine build_pot(this)
      class(symbolic_atom), intent(inout) :: this
      integer :: i

      ! Setting the potential parameters
      ! Imaginary part is set to 0
      ! kind is rp (same as real part)
      this%potential%cx(1, :) = cmplx(this%potential%center_band(1, :), 0, rp)
      this%potential%wx(1, :) = cmplx(this%potential%width_band(1, :), 0, rp)
      this%potential%cex(1, :) = cmplx(this%potential%shifted_band(1, :), 0, rp)
      this%potential%obx(1, :) = cmplx(this%potential%obar(1, :), 0, rp)
      do i = 2, 4
         this%potential%cx(i, :) = cmplx(this%potential%center_band(2, :), 0, rp)
         this%potential%wx(i, :) = cmplx(this%potential%width_band(2, :), 0, rp)
         this%potential%cex(i, :) = cmplx(this%potential%shifted_band(2, :), 0, rp)
         this%potential%obx(i, :) = cmplx(this%potential%obar(2, :), 0, rp)
      end do
      do i = 5, 9
         this%potential%cx(i, :) = cmplx(this%potential%center_band(3, :), 0, rp)
         this%potential%wx(i, :) = cmplx(this%potential%width_band(3, :), 0, rp)
         this%potential%cex(i, :) = cmplx(this%potential%shifted_band(3, :), 0, rp)
         this%potential%obx(i, :) = cmplx(this%potential%obar(3, :), 0, rp)
      end do

      this%potential%cx0(:) = 0.5D0 * (this%potential%cx(:, 1) + this%potential%cx(:, 2))
      this%potential%cx1(:) = 0.5D0 * (this%potential%cx(:, 1) - this%potential%cx(:, 2))
      this%potential%wx0(:) = 0.5D0 * (this%potential%wx(:, 1) + this%potential%wx(:, 2))
      this%potential%wx1(:) = 0.5D0 * (this%potential%wx(:, 1) - this%potential%wx(:, 2))
      this%potential%cex0(:) = 0.5D0 * (this%potential%cex(:, 1) + this%potential%cex(:, 2))
      this%potential%cex1(:) = 0.5D0 * (this%potential%cex(:, 1) - this%potential%cex(:, 2))
      this%potential%obx0(:) = 0.5D0 * (this%potential%obx(:, 1) + this%potential%obx(:, 2))
      this%potential%obx1(:) = 0.5D0 * (this%potential%obx(:, 1) - this%potential%obx(:, 2))
   end subroutine build_pot

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Transform the potential parameters from the orthogonal to the tight-binding representation.
   !>
   !> Transform the potential parameters from the orthogonal to the tight-binding representation.
   !> @param[in] class potential variable
   !---------------------------------------------------------------------------
   subroutine predls(this, wsm)
      class(symbolic_atom), intent(inout) :: this
      real(rp) :: wsm
      real(rp) :: x, y, wow, vmad
      real(rp), dimension(this%potential%lmax + 1, 2) :: dele, qi, del, q, c, enu
      real(rp), dimension(3) :: qm
      integer :: i, j, lmax

      wow = wsm / this%potential%ws_r
      lmax = this%potential%lmax
      enu(:, :) = this%potential%enu(0:lmax, :)
      c(:, :) = this%potential%c(0:lmax, :)
      del(:, :) = this%potential%srdel(0:lmax, :)
      q(:, :) = this%potential%qpar(0:lmax, :)
      vmad = this%potential%vmad
      lmax = this%potential%lmax
      qm(:) = qm_canonical(:)

      do J = 1, 2
         do I = 1, lmax + 1
            DELE(I, J) = DEL(I, J) * WOW**(0.5D0 - I)
            QI(I, J) = Q(I, J) * WOW**(1 - 2 * I)
            X = 1 - ((QI(I, J) - QM(I)) * (C(I, J) - ENU(I, J)) / DELE(I, J) / DELE(I, J))
            Y = (QI(I, J) - QM(I)) / &
                (((C(I, J) - ENU(I, J)) * (QI(I, J) - QM(I))) - (DELE(I, J) * DELE(I, J)))
            ! Saving parameters into the tight-binding representation
            this%potential%center_band(i, j) = (C(I, J) - ENU(I, J)) * X + ENU(I, J) + VMAD
            this%potential%shifted_band(i, j) = (C(I, J) - ENU(I, J)) * X
            this%potential%width_band(i, j) = DELE(I, J) * X
            this%potential%obar(i, j) = Y
            this%potential%qi(i - 1, j) = qi(i, j)
            this%potential%dele(i - 1, j) = dele(i, j)
         end do
      end do
   end subroutine predls

   subroutine d_matrix(this, mat, e)
      class(symbolic_atom) :: this
      !
      real(rp), intent(in) :: e
      complex(rp), dimension(9, 9), intent(inout) :: mat
      integer :: k, l, m, ml
      complex(rp):: cu, cd, wu, wd, de, wuwd

      mat = 0.0D0
      do l = 0, 2
         ml = l * l + 1
         cu = cmplx(this%potential%c(l, 1) + this%potential%vmad, 0.0D0)
         cd = cmplx(this%potential%c(l, 2) + this%potential%vmad, 0.0D0)
         wu = cmplx(this%potential%dele(l, 1), 0.0D0)
         wd = cmplx(this%potential%dele(l, 2), 0.0D0)
         wuwd = wu * wd
         wu = wu * wu
         wd = wd * wd
         de = (cd * wu - cu * wd + 1.0D0 * (wd - wu) * e) / (wuwd)
         do m = 1, 2 * l + 1
            ml = l * l + m
            mat(ml, ml) = de
         end do
      end do
   end subroutine d_matrix

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculate the displacement matrix (approx. displacement of the
   !> regular solution of the Laplace equation)
   !> Implemented by Ivan Miranda on 24.09.2023
   !---------------------------------------------------------------------------
   subroutine disp_matrix(this, mat, disp_vec, l_max, ws_radius)
      !
      class(symbolic_atom) :: this
      !
      integer, intent(in) :: l_max ! maximum quantum number l
      real(rp), intent(in) :: ws_radius ! Wigner-Seitz radius of the atom in question
      real(rp), dimension(3), intent(in) :: disp_vec ! displacement vector of the atom (can be non-unitary)
      complex(rp), allocatable, dimension(:, :), intent(inout) :: mat
      !
      ! Local Variables
      integer :: i, j, k, n, m
      complex(rp) :: temp1, temp2
      integer, allocatable, dimension(:, :) :: order ! order of orbitals (l,m)
      real(rp), dimension(3) :: unit_disp
      complex(rp), allocatable, dimension(:, :) :: mat_b ! smaller block

      allocate (order(l_max + 1, (2 * l_max + 1)))
      allocate (mat_b((l_max + 1)**2, (l_max + 1)**2))

      ! Hard check for zero displacement
      if (abs(norm2(disp_vec)) .eq. 0) then
         unit_disp(:) = 0.0_RP
      else
         unit_disp(1) = disp_vec(1) / norm2(disp_vec)
         unit_disp(2) = disp_vec(2) / norm2(disp_vec)
         unit_disp(3) = disp_vec(3) / norm2(disp_vec)
      end if

      order(:, :) = 0

      do i = 0, l_max
         if (i .le. 2) then
            ! Here performs the order of the matrices in RS-LMTO-ASA
            ! Implemented as s, px, py, pz, dxy, dyz, dzx, x2-y2, 3z2-r2
            ! for the spd basis
            order(1, 1) = 1
            order(2, 1) = 3; order(2, 2) = 4; order(2, 3) = 2
            order(3, 1) = 5; order(3, 2) = 6; order(3, 3) = 9
            order(3, 4) = 7; order(3, 5) = 8
         else
            ! For the other basis, use the regular (l, m) ordering
            ! i.e., for l >= 3 (f or higer) the basis will be
            ! (3, -3), (3, -2), (3, -1), ...
            do j = -i, i, 1
               order((i + 1), (i + j + 1)) = (i**2) + (i + j + 1)
            end do
         end if
      end do

      ! In the collinear state, it's not dependent on the spin
      ! so the off-diagonal blocks are zero (up-dw, dw-up)
      mat(1:(l_max + 1)**2, ((l_max + 1)**2 + 1):(2 * (l_max + 1)**2)) = czero
      mat(((l_max + 1)**2 + 1):(2 * (l_max + 1)**2), 1:(l_max + 1)**2) = czero
      mat_b(:, :) = czero ! define the initial element as imaginary zero

      do i = 0, l_max ! quantum angular momentum index (l')
         do j = 0, l_max ! quantum angular momentum index (l)
            do k = -i, i, 1 ! magnetic quantum number index (m')
               do n = -j, j, 1 ! magnetic quantum number index (m)
                  do m = -1, 1, 1 ! magnetic quantum number index (m'')
                     if (i .gt. j) then
                        mat_b(order(i + 1, k + i + 1), order(j + 1, n + j + 1)) = czero
                     else
                        temp1 = (factorial2(2.0_RP * j - 1.0_RP) / factorial2(2.0_RP * i - 1.0_RP)) * cone
                        temp2 = realgaunt(j, i, 1, n, k, m) * real_spharm(unit_disp, 1, m)
                        mat_b(order(i + 1, k + i + 1), order(j + 1, n + j + 1)) = &
                           mat_b(order(i + 1, k + i + 1), order(j + 1, n + j + 1)) + (temp1 * temp2)
                     end if
                  end do
               end do
            end do
         end do
      end do

      mat_b(:, :) = mat_b(:, :) * ((-1.0_RP * four_pi) / (3.0_RP * ws_radius)) * cone
      mat(1:(l_max + 1)**2, 1:(l_max + 1)**2) = mat_b
      mat(((l_max + 1)**2 + 1):(2 * (l_max + 1)**2), ((l_max + 1)**2 + 1):(2 * (l_max + 1)**2)) = mat_b

      deallocate (order)
      deallocate (mat_b)

   end subroutine disp_matrix

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculate the U matrix (first order - linear - approximation to the
   !> change in the potential parameter P after an atomic displacement)
   !> Implemented by Ivan Miranda on 25.09.2023
   !---------------------------------------------------------------------------
   subroutine udisp_matrix(this, mat, pmat, disp_vec, l_max, ws_radius, energy_steps)
      !
      class(symbolic_atom) :: this
      !
      integer, intent(in) :: l_max ! maximum quantum number l
      integer, intent(in) :: energy_steps ! Energy integration points
      real(rp), intent(in) :: ws_radius ! Wigner-Seitz radius of the atom in question
      real(rp), dimension(3), intent(in) :: disp_vec ! displacement vector of the atom (can be non-unitary)
      complex(rp), allocatable, dimension(:, :, :), intent(in) :: pmat ! P potential paramter matrix in LMTO for the site
      complex(rp), allocatable, dimension(:, :, :), intent(inout) :: mat ! Final U matrix
      !
      ! Local Variables
      integer :: i
      complex(rp), allocatable, dimension(:, :) :: dmat ! Displacement matrix

      allocate (dmat(2 * (l_max + 1)**2, 2 * (l_max + 1)**2))

      dmat(:, :) = czero
      call disp_matrix(this, dmat, disp_vec, l_max, ws_radius)
      !write(*, '(9F14.9)') dmat(1:9, 1:9)
      !write(*, *) '-------'

      do i = 1, energy_steps
         mat(:, :, i) = matmul(dmat(:, :), pmat(:, :, i)) + matmul(pmat(:, :, i), transpose(dmat(:, :)))
      end do

      deallocate (dmat)

   end subroutine udisp_matrix

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Potential parameter P matrix calculated considering the
   !> orthogonal representation.
   !> Implemented by Ivan Miranda on 27.09.2023
   !---------------------------------------------------------------------------
   subroutine p_matrix(this, pmat, l_max, energy_array)
      !
      class(symbolic_atom) :: this
      !
      integer, intent(in) :: l_max ! maximum quantum number l
      real(rp), intent(in) :: energy_array(:) ! Energy array
      complex(rp), intent(inout) :: pmat(:, :, :) ! P matrix
      !
      ! Local Variables
      integer :: s, l, m, mls, ie
      complex(rp) :: temp1, temp2, temp3
      !

      pmat(:, :, :) = czero

      do s = 1, 2 ! spin index
         do l = 0, l_max ! l index
            do m = 1, 2 * l + 1
               mls = l * l + m + ((l_max + 1)**2) * (s - 1) ! Composed diagonal index
               temp1 = cmplx(this%potential%c(l, s) + this%potential%vmad, 0.0_RP)
               temp2 = cmplx(this%potential%dele(l, s), 0.0_RP)
               do ie = 1, size(energy_array)
                  temp3 = cmplx(energy_array(ie), 0.0_RP)
                  pmat(mls, mls, ie) = (temp3 - temp1) / (temp2 * temp2)
               end do
            end do
         end do
      end do

   end subroutine p_matrix

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Transformation between the screened (or not) representations
   !> of the potential function P.
   !> This assumes that the P matrices are diagonal - or spherically
   !> symmetrical (normal situation for ASA and equilibrium positions).
   !> See Turek's book Eq. (3.56)
   !> Implemented by Ivan Miranda on 24.10.2023
   !---------------------------------------------------------------------------
   subroutine transform_pmatrix(this, pmat_in, pmat_out, screening_in, screening_out)
      !
      class(symbolic_atom) :: this
      !
      real(rp), intent(in) :: screening_in(0:, :) ! Screening constants relative to pmat_in
      real(rp), intent(in) :: screening_out(0:, :) ! Screening constant relative to pmat_out
      complex(rp), intent(in) :: pmat_in(:, :, :) ! Potential P matrix (in)
      complex(rp), intent(inout) :: pmat_out(:, :, :) ! Potential P matrix (out)
      !
      ! Local Variables
      integer :: s, l, m, mls, ie
      complex(rp) :: temp1, temp2
      !

      pmat_out(:, :, :) = czero

      ! Screening_in and out should have the same size, of course
      do s = 1, size(screening_in, 2) ! spin index
         do l = 0, size(screening_in, 1) - 1 ! lmax
            do m = 1, 2 * l + 1
               mls = l * l + m + ((size(screening_in, 1))**2) * (s - 1) ! Composed diagonal index
               temp1 = cmplx(screening_in(l, s), 0.0_RP) ! Screening constants (in)
               temp2 = cmplx(screening_out(l, s), 0.0_RP) ! Screening constants (out)
               do ie = 1, size(pmat_in, 3) ! Energy channel
                  pmat_out(mls, mls, ie) = pmat_in(mls, mls, ie) / (cone + ((temp1 - temp2) * pmat_in(mls, mls, ie)))
               end do
            end do
         end do
      end do

   end subroutine transform_pmatrix

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Print class members values in namelist format
   !>
   !> Print class members values in namelist format
   !---------------------------------------------------------------------------
   subroutine symbolic_atom_print_state(this, unit, file)
      class(symbolic_atom) :: this
      integer, optional :: unit
      character(len=*), optional :: file
      integer :: newunit
      if (present(unit) .and. present(file)) then
         call g_logger%fatal('Argument error: both unit and file are present', __FILE__, __LINE__)
      else if (present(unit)) then
         call this%element%print_state(unit=unit)
         call this%potential%print_state(unit=unit)
      else if (present(file)) then
         open (newunit=newunit, file=file, action='write')
         call this%element%print_state(unit=newunit)
         call this%potential%print_state(unit=newunit)
         close (newunit)
      else
         call this%element%print_state()
         call this%potential%print_state()
      end if
   end subroutine symbolic_atom_print_state

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Print class members values in namelist format
   !>
   !> Print class members values in namelist format
   !---------------------------------------------------------------------------
   subroutine symbolic_atom_print_state_full(this, unit, file)
      class(symbolic_atom) :: this
      integer, optional :: unit
      character(len=*), optional :: file
      integer :: newunit
      if (present(unit) .and. present(file)) then
         call g_logger%fatal('Argument error: both unit and file are present', __FILE__, __LINE__)
      else if (present(unit)) then
         call this%element%print_state_full(unit=unit)
         call this%potential%print_state_full(unit=unit)
      else if (present(file)) then
         open (newunit=newunit, file=file, action='write')
         call this%element%print_state_full(unit=newunit)
         call this%potential%print_state_full(unit=newunit)
         close (newunit)
      else
         call this%element%print_state_full()
         call this%potential%print_state_full()
      end if
   end subroutine symbolic_atom_print_state_full

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Print class members values in namelist format
   !>
   !> Print class members values in namelist format
   !---------------------------------------------------------------------------
   subroutine symbolic_atom_print_state_formatted(this, unit, file)
      class(symbolic_atom) :: this
      integer, optional :: unit
      character(len=*), optional :: file
      integer :: newunit
      if (present(unit) .and. present(file)) then
         call g_logger%fatal('Argument error: both unit and file are present', __FILE__, __LINE__)
      else if (present(unit)) then
         call this%element%print_state_formatted(unit=unit)
         call this%potential%print_state_formatted(unit=unit)
      else if (present(file)) then
         open (newunit=newunit, file=file, action='write')
         call this%element%print_state_formatted(unit=newunit)
         call this%potential%print_state_formatted(unit=newunit)
         close (newunit)
      else
         call this%element%print_state_formatted()
         call this%potential%print_state_formatted()
      end if
   end subroutine symbolic_atom_print_state_formatted

   !---------------------------------------------------------------------------
   ! DESCRIPTION:

   !> TODO
   !>
   !> TODO
   !---------------------------------------------------------------------------
   function mesh_grid_size(this)
      class(symbolic_atom), intent(in) :: this
      integer :: mesh_grid_size
      real(rp) :: B, Z
      intrinsic EXP, LOG, MAX

      Z = this%element%atomic_number
      B = 1.D0 / (Z + Z + 1.D0)
      mesh_grid_size = max(min_mesh_grid_size, int(((.5D0 + log(1.D0 + this%potential%ws_r / B) / this%A) * 2.D0 - 1) / 2) * 2 + 1)
   end function mesh_grid_size

   !---------------------------------------------------------------------------
   !> DESCRIPTION:
   !> @breif
   !> TODO Calculates the quantity B whatever it is
   !---------------------------------------------------------------------------
   function B(this)
      class(symbolic_atom), intent(in) :: this
      real(rp) :: B
      B = this%potential%ws_r / (exp(this%A * this%mesh_grid_size() - this%A) - 1.D0)
   end function B

   !---------------------------------------------------------------------------
   !> DESCRIPTION:
   !> @breif
   !> TODO Calculates the initial guess for charge density
   !---------------------------------------------------------------------------
   function rho0(this, nsp)
      class(symbolic_atom), intent(in) :: this
      integer, intent(in), optional :: nsp
      real(rp), dimension(:, :), allocatable :: rho0
      real(rp) :: b, ea, sum, rpb, r, ro, fac
      integer :: ir, nr
      integer :: nsp_
      !nsp_ = merge(nsp, 1, present(nsp))
      if (present(nsp)) then
         nsp_ = nsp
      else
         nsp_ = 1
      end if
      nr = this%mesh_grid_size()
      allocate (rho0(nr, nsp_))
      b = this%B()
      ea = exp(this%a)
      sum = 0.D0
      rpb = b
      do ir = 1, nr
         r = rpb - b
         ro = exp(-5.*r) * r * r
         rho0(ir, 1) = ro
         sum = sum + this%a * rpb * ro
         rpb = rpb * ea
      end do
      fac = this%element%atomic_number / (sum * 2)
      do ir = 1, nr
         rho0(ir, 1) = rho0(ir, 1) * fac
         rho0(ir, 2) = rho0(ir, 1)
         !write(101, *) ir, rho0(ir, 1)
         !write(102, *) ir, rho0(ir, 2)
      end do
   end function rho0

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Build an array of symbolic_atoms
   !
   !> @param[in] symbolic_atom List of labels in database to build a symbolic_atom
   !> @param[in] database Directory with files used as database
   !> @return type(symbolic_atom), dimension(:), allocatable
   !---------------------------------------------------------------------------
   function array_of_symbolic_atoms_from_memory(symbolic_atoms, database)
      type(symbolic_atom), dimension(:), allocatable :: array_of_symbolic_atoms_from_memory
      character(len=*), dimension(:), intent(in) :: symbolic_atoms
      character(len=*), intent(in), optional :: database
      integer :: i, j
      allocate (array_of_symbolic_atoms_from_memory(size(symbolic_atoms)))
      if (present(database)) then
         do i = 1, size(symbolic_atoms)
            array_of_symbolic_atoms_from_memory(i) = symbolic_atom(symbolic_atoms(i), database)
         end do
      else
         do i = 1, size(symbolic_atoms)
            array_of_symbolic_atoms_from_memory(i) = symbolic_atom(symbolic_atoms(i))
         end do
      end if
   end function array_of_symbolic_atoms_from_memory

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Build an array of symbolic_atoms from namelist 'atoms' in fname
   !
   !> @param[in] fname File containg the namelist 'atoms'
   !> @param[in] size Size (integer) or array
   !> @return type(symbolic_atom), dimension(:), allocatable
   !---------------------------------------------------------------------------
   function array_of_symbolic_atoms_from_file(fname, size)
      use string_mod, only: sl
      type(symbolic_atom), dimension(:), allocatable :: array_of_symbolic_atoms_from_file
      integer, intent(in) :: size
      character(len=*), intent(in) :: fname

      ! variables associated with the reading processes
      integer :: iostatus, funit

      include 'include_codes/namelists/symbolic_atom.f90'

      allocate (label(size))

      open (newunit=funit, file=fname, action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) then
         call g_logger%fatal("file "//trim(fname)//" not found", __FILE__, __LINE__)
      end if

      read (funit, nml=atoms, iostat=iostatus)
      if (iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
         call g_logger%error('Error while reading namelist', __FILE__, __LINE__)
         call g_logger%error('iostatus = '//fmt('I0', iostatus), __FILE__, __LINE__)
      end if
      close (funit)

      array_of_symbolic_atoms_from_file = array_of_symbolic_atoms(label, database)

   end function array_of_symbolic_atoms_from_file

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Print class members values in namelist format
   !>
   !> Print class members values in namelist format. Either unit or file should be provided. If none of them are provided, then the program will write to standart output.
   !> @param[in] unit File unit used to write namelist
   !> @param[in] file File name used to write namelist
   !> @param[in] suffix If provided, prints the state to file with the same symbol added to this suffix
   !---------------------------------------------------------------------------
   subroutine print_state(array, unit, file, suffix)
      type(symbolic_atom), dimension(:), intent(in) :: array
      integer, intent(in), optional :: unit
      character(len=*), intent(in), optional :: file, suffix
      character(len=sl) :: new_file
      integer :: i, newunit
      if ((present(unit) .and. present(file)) .or. \
      (present(unit) .and. present(suffix)) .or. \
      (present(file) .and. present(suffix))) then
      call g_logger%fatal("Argument error: both unit, file and suffix are present", __FILE__, __LINE__)
      else if (present(file)) then
      open (newunit=newunit, file=file)
      do i = 1, size(array)
         call array(i)%print_state(unit=newunit)
      end do
      close (newunit)
      else if (present(suffix)) then
      do i = 1, size(array)
         new_file = trim(array(i)%element%symbol)
         if (endswith(new_file, ".nml")) then
            call replace(new_file, ".nml", "")
         end if
         new_file = trim(new_file)//trim(suffix)//".nml"
         call array(i)%print_state(file=new_file)
      end do
      else
      do i = 1, size(array)
         call array(i)%print_state(unit=unit)
      end do
      end if
      end subroutine print_state

      !---------------------------------------------------------------------------
      ! DESCRIPTION:
      !> @brief
      !> Print class members values in namelist format
      !>
      !> Print class members values in namelist format. Either unit or file should be provided. If none of them are provided, then the program will write to standart output.
      !> @param[in] unit File unit used to write namelist
      !> @param[in] file File name used to write namelist
      !> @param[in] suffix If provided, prints the state to file with the same symbol added to this suffix
      !---------------------------------------------------------------------------
      subroutine print_state_full(array, unit, file, suffix)
         type(symbolic_atom), dimension(:), intent(in) :: array
         integer, intent(in), optional :: unit
         character(len=*), intent(in), optional :: file, suffix
         character(len=sl) :: new_file
         integer :: i, newunit
         if ((present(unit) .and. present(file)) .or. \
         (present(unit) .and. present(suffix)) .or. \
         (present(file) .and. present(suffix))) then
         call g_logger%fatal("Argument error: both unit, file and suffix are present", __FILE__, __LINE__)
         else if (present(file)) then
         open (newunit=newunit, file=file, action='write')
         do i = 1, size(array)
            call array(i)%print_state_full(unit=newunit)
         end do
         close (newunit)
         else if (present(suffix)) then
         do i = 1, size(array)
            new_file = array(i)%element%symbol
            if (endswith(new_file, ".nml")) then
               call replace(new_file, ".nml", "")
            end if
            new_file = trim(new_file)//trim(suffix)//".nml"
            call array(i)%print_state_full(file=new_file)
         end do
         else
         do i = 1, size(array)
            call array(i)%print_state_full(unit=unit)
         end do
         end if
      end subroutine print_state_full

      !---------------------------------------------------------------------------
      ! DESCRIPTION:
      !> @brief
      !> Print class members values in namelist format
      !>
      !> Print class members values in namelist format. Either unit or file should be provided. If none of them are provided, then the program will write to standart output.
      !> @param[in] unit File unit used to write namelist
      !> @param[in] file File name used to write namelist
      !> @param[in] suffix If provided, prints the state to file with the same symbol added to this suffix
      !---------------------------------------------------------------------------
      subroutine print_state_formatted(array, unit, file, suffix)
         type(symbolic_atom), dimension(:), intent(in) :: array
         integer, intent(in), optional :: unit
         character(len=*), intent(in), optional :: file, suffix
         character(len=sl) :: new_file
         integer :: i, newunit
         if ((present(unit) .and. present(file)) .or. \
         (present(unit) .and. present(suffix)) .or. \
         (present(file) .and. present(suffix))) then
         call g_logger%fatal("Argument error: both unit, file and suffix are present", __FILE__, __LINE__)
         else if (present(file)) then
         open (newunit=newunit, file=file)
         do i = 1, size(array)
            call array(i)%print_state_formatted(unit=newunit)
         end do
         close (newunit)
         else if (present(suffix)) then
         do i = 1, size(array)
            new_file = trim(array(i)%element%symbol)
            if (endswith(new_file, ".nml")) then
               call replace(new_file, ".nml", "")
            end if
            new_file = trim(new_file)//trim(suffix)//".nml"
            call array(i)%print_state_formatted(file=new_file)
         end do
         else
         do i = 1, size(array)
            call array(i)%print_state_formatted(unit=unit)
         end do
         end if
      end subroutine print_state_formatted

      !---------------------------------------------------------------------------
      ! DESCRIPTION:
      !> @brief
      !> Shortcut for print_state_full with suffix="_out"
      !>
      !> Shortcut for print_state_full with suffix="_out"
      !---------------------------------------------------------------------------
      subroutine save_state(array)
         type(symbolic_atom), dimension(:), intent(in) :: array
         call print_state_formatted(array, suffix="_out")
      end subroutine save_state

      !---------------------------------------------------------------------------
      ! DESCRIPTION:
      !> @brief
      !> Checks the existence of checkpoint files (those with fname = symbol//"_out.nml" )
      !>
      !> Checks the existence of checkpoint files (those with fname = symbol//"_out.nml" )
      !---------------------------------------------------------------------------
      subroutine load_state(array)
         type(symbolic_atom), dimension(:), intent(inout) :: array
         character(len=sl) :: fname
         integer :: i, array_size
         logical :: file_exists

         array_size = size(array)
         do i = 1, array_size
            fname = trim(array(i)%element%symbol)//"_out.nml"
            inquire (FILE=fname, EXIST=file_exists)
            if (file_exists) then
               call array(i)%build_from_file(fname)
            end if
         end do
      end subroutine load_state
   end module symbolic_atom_mod
