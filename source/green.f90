!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Green
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
!> Module to handle Greens functions calculation and related routines
!------------------------------------------------------------------------------

module green_mod

   use energy_mod
   use control_mod
   use lattice_mod
   use symbolic_atom_mod
   use recursion_mod
   use density_of_states_mod
   use basis_mod, only: nb
   use precision_mod, only: rp
   use math_mod
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc
#endif
   implicit none

   private

   !> Module´s main structure
   type, public :: green
      ! General variables
      !> Recursion
      class(recursion), pointer :: recursion
      !> Lattice
      class(lattice), pointer :: lattice
      !> Symbolic atom
      class(symbolic_atom), dimension(:), pointer :: symbolic_atom
      !> Density of states
      class(dos), pointer :: dos
      !> Control
      class(control), pointer :: control
      !> Energy
      class(energy), pointer :: en
      !> Onsite Green Function
      complex(rp), dimension(:, :, :, :), allocatable :: g0
      !> Intersite Green Function
      complex(rp), dimension(:, :, :, :), allocatable :: gij, gji, ginmag, gjnmag, gix, giy, giz, gjx, gjy, gjz
      !> Intersite Green Function as a function of imaginary eta
      complex(rp), dimension(:, :, :, :), allocatable :: gij_eta, gji_eta, ginmag_eta, gjnmag_eta, gix_eta, giy_eta, giz_eta, gjx_eta, gjy_eta, gjz_eta
      !> Intersite Green Function with two indexes
      complex(rp), dimension(:, :, :, :), allocatable :: g00ij, g01ij, g00ji, g01ji, &
                                                       & gx1ij, gy1ij, gz1ij, gx0ij, gy0ij, gz0ij, &
                                                       & gx1ji, gy1ji, gz1ji, gx0ji, gy0ji, gz0ji

   contains
      procedure :: sgreen
      procedure :: bgreen
      procedure :: bgreen_complex
      procedure :: block_green
      procedure :: block_green_gpu
      procedure :: block_green_eta
      procedure :: block_green_ij
      procedure :: block_green_ij_complex
      procedure :: block_green_ij_gpu
      procedure :: block_green_ij_eta
      procedure :: calculate_intersite_gf
      procedure :: calculate_intersite_gf_complex
      procedure :: calculate_intersite_gf_twoindex
      procedure :: calculate_intersite_gf_eta
      procedure :: calculate_intersite_gf_eta_gpu
      procedure :: chebyshev_green
      procedure :: chebyshev_green_gpu
      procedure :: chebyshev_green_eta
      procedure :: chebyshev_green_ij
      procedure :: chebyshev_green_ij_complex
      procedure :: chebyshev_green_ij_gpu
      procedure :: chebyshev_green_ij_eta
      procedure :: chebyshev_dos_dispatch
      procedure :: restore_to_default
      procedure :: auxiliary_gij
      procedure :: transform_auxiliary_gij
      procedure :: gij_eta_to_gij
      final :: destructor
   end type

   interface green
      procedure :: constructor
   end interface green

   interface

   ! DESCRIPTION:
   !> @brief
   !> Constructor
   !
   !> @param[in] fname Namelist file
   !> @return type(green)
   !---------------------------------------------------------------------------
   module function constructor(dos_obj) result(obj)
      type(green) :: obj
      type(dos), target, intent(in) :: dos_obj
   end function constructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Destructor
   !---------------------------------------------------------------------------
   module subroutine destructor(this)
      type(green) :: this
   end subroutine destructor

   ! Member functions
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Reset all members to default
   !---------------------------------------------------------------------------
   module subroutine restore_to_default(this)
      class(green) :: this
   end subroutine restore_to_default

   module subroutine chebyshev_green_ij(this, istart)
      class(green), intent(inout) :: this
      integer, intent(in) :: istart
   end subroutine chebyshev_green_ij

   module subroutine chebyshev_green_ij_complex(this, istart, z_grid, g_ef)
      class(green), intent(inout) :: this
      integer, intent(in) :: istart
      complex(rp), intent(in) :: z_grid(:)
      complex(rp), intent(out) :: g_ef(:, :, :, :)
   end subroutine chebyshev_green_ij_complex

   module subroutine chebyshev_green_ij_gpu(this, istart)
      class(green), intent(inout) :: this
      integer, intent(in) :: istart
   end subroutine chebyshev_green_ij_gpu

   module subroutine chebyshev_green_ij_eta(this, istart, eta, fermi_point, g_ef)
      class(green), intent(inout) :: this
      integer, intent(in) :: istart
      complex(rp), dimension(nb, nb, 4), intent(inout) :: g_ef
      complex(rp), intent(in) :: eta
      integer, intent(in) :: fermi_point
   end subroutine chebyshev_green_ij_eta

   module subroutine chebyshev_green(this)
      class(green), intent(inout) :: this
   end subroutine chebyshev_green

   module subroutine chebyshev_green_gpu(this)
      class(green), intent(inout) :: this
   end subroutine chebyshev_green_gpu

   module subroutine chebyshev_dos_dispatch(this)
      class(green), intent(inout) :: this
   end subroutine chebyshev_dos_dispatch

   module subroutine chebyshev_green_eta(this, eta, fermi_point, g_ef)
      use mpi_mod, only: atoms_per_process
      class(green), intent(inout) :: this
      complex(rp), dimension(nb, nb, atoms_per_process), intent(inout) :: g_ef
      complex(rp), intent(in) :: eta
      integer, intent(in) :: fermi_point
   end subroutine chebyshev_green_eta

   module subroutine chebyshev_green_core(this, eta, fermi_point, g_ef, use_gpu)
      use mpi_mod, only: atoms_per_process
      class(green), intent(inout) :: this
      complex(rp), intent(in), optional :: eta
      integer, intent(in), optional :: fermi_point
      complex(rp), dimension(nb, nb, atoms_per_process), intent(inout), optional :: g_ef
      logical, intent(in), optional :: use_gpu
   end subroutine chebyshev_green_core

   ! Recursive continued-fraction Green-function implementations.
   module subroutine block_green_ij_eta(this, istart, eta, fermi_point, g_ef)
      class(green), intent(inout) :: this
      integer, intent(in) :: istart
      complex(rp), intent(in) :: eta
      complex(rp), dimension(nb, nb, 4), intent(inout) :: g_ef
      integer, intent(in) :: fermi_point
   end subroutine block_green_ij_eta

   module subroutine block_green_ij(this, istart)
      class(green), intent(inout) :: this
      integer, intent(in) :: istart
   end subroutine block_green_ij

   module subroutine block_green_ij_complex(this, istart, z_grid, g_ef)
      class(green), intent(inout) :: this
      integer, intent(in) :: istart
      complex(rp), intent(in) :: z_grid(:)
      complex(rp), intent(out) :: g_ef(:, :, :, :)
   end subroutine block_green_ij_complex

   module subroutine block_green_ij_gpu(this, istart)
      class(green), intent(inout) :: this
      integer, intent(in) :: istart
   end subroutine block_green_ij_gpu

   module subroutine calculate_intersite_gf_twoindex(this)
      class(green), intent(inout) :: this
   end subroutine calculate_intersite_gf_twoindex

   module subroutine calculate_intersite_gf(this)
      class(green), intent(inout) :: this
   end subroutine calculate_intersite_gf

   !> Build native intersite Green functions for a batch of arbitrary complex
   !> energies.  The returned arrays retain the local pair ownership and
   !> `(orbital,orbital,energy,pair)` layout of `green%gij/gji`.
   module subroutine calculate_intersite_gf_complex(this, z, g_ab, g_ba)
      use mpi_mod, only: atoms_per_process
      class(green), intent(inout) :: this
      complex(rp), intent(in) :: z(:)
      complex(rp), allocatable, intent(out) :: g_ab(:, :, :, :), g_ba(:, :, :, :)
   end subroutine calculate_intersite_gf_complex

   module subroutine calculate_intersite_gf_eta(this)
      class(green), intent(inout) :: this
   end subroutine calculate_intersite_gf_eta

   module subroutine calculate_intersite_gf_core(this, eta_mode)
      class(green), intent(inout) :: this
      logical, intent(in) :: eta_mode
   end subroutine calculate_intersite_gf_core

   module subroutine calculate_intersite_gf_eta_gpu(this)
      class(green), intent(inout) :: this
   end subroutine calculate_intersite_gf_eta_gpu

   module subroutine block_green_eta(this, eta, fermi_point, g_ef)
      use mpi_mod, only: atoms_per_process
      class(green), intent(inout) :: this
      complex(rp), intent(in) :: eta
      integer, intent(in) :: fermi_point
      complex(rp), dimension(nb, nb, atoms_per_process), intent(inout) :: g_ef
   end subroutine block_green_eta

   module subroutine block_green(this)
      class(green), intent(inout) :: this
   end subroutine block_green

   module subroutine block_green_gpu(this)
      class(green), intent(inout) :: this
   end subroutine block_green_gpu

   module subroutine block_green_core(this, eta, fermi_point, g_ef)
      use mpi_mod, only: atoms_per_process
      class(green), intent(inout) :: this
      complex(rp), intent(in), optional :: eta
      integer, intent(in), optional :: fermi_point
      complex(rp), dimension(nb, nb, atoms_per_process), intent(inout), optional :: g_ef
   end subroutine block_green_core

   module subroutine sgreen(this)
      class(green) :: this
   end subroutine sgreen

   module subroutine bgreen(this, g_out, i_site, ie_start, ie_len, a_inf, b_inf, eta)
      class(green), intent(inout) :: this
      integer, intent(in) :: i_site
      integer, intent(in) :: ie_start
      integer, intent(in) :: ie_len
      complex(rp), dimension(:, :, :), intent(inout) :: g_out
      real(rp), dimension(nb, nb), intent(in) :: a_inf
      real(rp), dimension(nb, nb), intent(in) :: b_inf
      complex(rp), intent(in) :: eta
   end subroutine bgreen

   !> Evaluate one native block-recursion Green function for a complex-energy
   !> batch.  `z_grid` is passed through the same core as the historical
   !> real-axis routine; it is not an interpolation of the DOS energy mesh.
   module subroutine bgreen_complex(this, g_out, i_site, z_grid, a_inf, b_inf, legacy_real_axis)
      class(green), intent(inout) :: this
      integer, intent(in) :: i_site
      complex(rp), intent(in) :: z_grid(:)
      complex(rp), dimension(:, :, :), intent(inout) :: g_out
      real(rp), dimension(nb, nb), intent(in) :: a_inf
      real(rp), dimension(nb, nb), intent(in) :: b_inf
      logical, intent(in), optional :: legacy_real_axis
   end subroutine bgreen_complex

   end interface

contains
   !> @brief Reshape the 64-point eta-contour intersite GF into the same
   !>        (norb,norb,echan,atom) layout as the on-mesh arrays.
   !> @details Copies Ginmag_eta/Gjnmag_eta/Gi{x,y,z}_eta/Gj{x,y,z}_eta
   !>          (shape (echan,norb,norb,atoms_per_process), filled by
   !>          calculate_intersite_gf_core’s eta_mode) into freshly
   !>          (re)allocated Ginmag/Gjnmag/Gi{x,y,z}/Gj{x,y,z} arrays of shape
   !>          (norb,norb,64,atoms_per_process), transposing the orbital and
   !>          energy-channel index order. This lets
   !>          exchange%calculate_exchange_gauss_legendre reuse the same
   !>          orbital-indexed exchange kernels (dGdG and friends) that the
   !>          on-mesh gf_route=’recursion’ path uses, without a separate
   !>          eta-indexed code path.
   subroutine gij_eta_to_gij(this)
      use mpi_mod
      !
      class(green) :: this
      integer :: echan, ijij

      if (allocated(this%Ginmag)) deallocate (this%Ginmag)
      if (allocated(this%Gjnmag)) deallocate (this%Gjnmag)
      if (allocated(this%Gix)) deallocate (this%Gix)
      if (allocated(this%Gjx)) deallocate (this%Gjx)
      if (allocated(this%Giy)) deallocate (this%Giy)
      if (allocated(this%Gjy)) deallocate (this%Gjy)
      if (allocated(this%Giz)) deallocate (this%Giz)
      if (allocated(this%Gjz)) deallocate (this%Gjz)

      allocate (this%Ginmag(norb, norb, 64, atoms_per_process))
      allocate (this%Gjnmag(norb, norb, 64, atoms_per_process))
      allocate (this%Gix(norb, norb, 64, atoms_per_process))
      allocate (this%Gjx(norb, norb, 64, atoms_per_process))
      allocate (this%Giy(norb, norb, 64, atoms_per_process))
      allocate (this%Gjy(norb, norb, 64, atoms_per_process))
      allocate (this%Giz(norb, norb, 64, atoms_per_process))
      allocate (this%Gjz(norb, norb, 64, atoms_per_process))

      do ijij = 1, atoms_per_process
         do echan = 1, 64
            this%Ginmag(:, :, echan, ijij) = this%Ginmag_eta(echan, :, :, ijij)
            this%Gjnmag(:, :, echan, ijij) = this%Gjnmag_eta(echan, :, :, ijij)
            this%Gix(:, :, echan, ijij) = this%Gix_eta(echan, :, :, ijij)
            this%Gjx(:, :, echan, ijij) = this%Gjx_eta(echan, :, :, ijij)
            this%Giy(:, :, echan, ijij) = this%Giy_eta(echan, :, :, ijij)
            this%Gjy(:, :, echan, ijij) = this%Gjy_eta(echan, :, :, ijij)
            this%Giz(:, :, echan, ijij) = this%Giz_eta(echan, :, :, ijij)
            this%Gjz(:, :, echan, ijij) = this%Gjz_eta(echan, :, :, ijij)
         end do
      end do

      return

   end subroutine gij_eta_to_gij
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Given an input Gij (physical site-resolved GF in the LMTO), the function
   !> returns the correspondent auxiliary GFs, considering here the
   !> orthogonal representation. Takes into account all njij pairs directly.
   !> It uses the fact that the phyisical GF´s are invariant with respect to
   !> screening constants (Turek´s book, page 72-73).
   !> Implemented by Ivan Miranda on 29.09.2023
   !---------------------------------------------------------------------------
   subroutine auxiliary_gij(this, green_ij, aux_gij, atom_i, atom_j)
      !
      class(green) :: this
      !
      integer, intent(in) :: atom_i, atom_j ! Input atoms for the Gij
      complex(rp), intent(in) :: green_ij(:, :, :) ! Input Gij (physical site-resolved GFs)
      complex(rp), intent(inout) :: aux_gij(:, :, :) ! Output auxiliary gij matrix (same size as Gij)
      !
      ! Local Variables
      integer :: nv, l, m, mls, s, lmaxi, lmaxj
      complex(rp), allocatable, dimension(:, :) :: cdelta_i, cdelta_j

      allocate (cdelta_i(size(green_ij, 1), size(green_ij, 2)))
      allocate (cdelta_j(size(green_ij, 1), size(green_ij, 2)))

      ! Set all matrices initially to zero
      aux_gij(:, :, :) = czero
      ! Calculate the lmax values of the input atoms
      lmaxi = this%symbolic_atom(this%lattice%iz(atom_i))%potential%lmax
      lmaxj = this%symbolic_atom(this%lattice%iz(atom_j))%potential%lmax

      do nv = 1, size(this%en%ene) ! Energy channel

         cdelta_i(:, :) = czero ! As zero in the beginning
         cdelta_j(:, :) = czero ! As zero in the beginning

         do s = 1, 2
            do l = 0, lmaxi ! Transform delta_i in complex
               do m = 1, 2*l + 1
                  mls = l*l + m + ((lmaxi + 1)**2)*(s - 1) ! composed index
                  ! Sqrt(Delta) matrix of atom i
                  cdelta_i(mls, mls) = cmplx(this%symbolic_atom(this%lattice%iz(atom_i))%potential%dele(l, s), 0.0_rp)
               end do
            end do
         end do
         do s = 1, 2
            do l = 0, lmaxj ! Transform delta_j in complex
               do m = 1, 2*l + 1
                  mls = l*l + m + ((lmaxj + 1)**2)*(s - 1) ! composed index
                  ! Sqrt(Delta) matrix of atom j
                  cdelta_j(mls, mls) = cmplx(this%symbolic_atom(this%lattice%iz(atom_j))%potential%dele(l, s), 0.0_rp)
               end do
            end do
         end do
         aux_gij(:, :, nv) = matmul(cdelta_i, matmul(green_ij(:, :, nv), cdelta_j))
      end do

      deallocate (cdelta_i)
      deallocate (cdelta_j)

   end subroutine auxiliary_gij

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Given an auxiliary GF´s (in the LMTO language), the function
   !> returns the correspondent auxiliary GF´s transformed from
   !> one representation (alpha - in) to another (beta - out).
   !> The different representations are here denoted by their
   !> set of screening constants.
   !> See Turek´s book, Eq. (3.57)
   !> Implemented by Ivan Miranda on 25.10.2023
   !---------------------------------------------------------------------------
   subroutine transform_auxiliary_gij(this, pmat_in_atom_i, pmat_out_atom_i, pmat_in_atom_j, pmat_out_atom_j, &
                                      aux_gij_in, aux_gij_out, screening_in, screening_out, atom_i, atom_j)
      !
      class(green) :: this
      !
      integer, intent(in) :: atom_i, atom_j ! Input atoms i and j
      complex(rp), intent(in) :: aux_gij_in(:, :, :) ! Auxiliary GF´s (in)
      complex(rp), intent(inout) :: aux_gij_out(:, :, :) ! Auxiliary GF´s (out), in another representation
      complex(rp), intent(in) :: pmat_in_atom_i(:, :, :) ! Potential P matrix (in) - in the same representation as GF¨s (in) for the atom i
      complex(rp), intent(in) :: pmat_out_atom_i(:, :, :) ! Potential P matrix (out) - in the same representation as GF´s (out) for the atom i
      complex(rp), intent(in) :: pmat_in_atom_j(:, :, :) ! Potential P matrix (in) - in the same representation as GF´s (in) for the atom j
      complex(rp), intent(in) :: pmat_out_atom_j(:, :, :) ! Potential P matrix (out) - in the same representation as GF´s (out) for the atom j
      real(rp), intent(in) :: screening_in(0:, :) ! Screening constants relative to the GF´s (in)
      real(rp), intent(in) :: screening_out(0:, :) ! Screening constants relative to the GF´s (out)
      !
      ! Local Variables
      integer :: l, s, m, mls, ie
      complex(rp) :: temp1, temp2
      complex(rp), allocatable, dimension(:, :, :) :: pmat_resc1, pmat_resc2, pmat_resc3 ! Re-scale P matrices

      allocate (pmat_resc1(size(pmat_in_atom_i, 1), size(pmat_in_atom_i, 2), size(pmat_in_atom_i, 3)))
      allocate (pmat_resc2(size(pmat_in_atom_i, 1), size(pmat_in_atom_i, 2), size(pmat_in_atom_i, 3)))
      allocate (pmat_resc3(size(pmat_in_atom_i, 1), size(pmat_in_atom_i, 2), size(pmat_in_atom_i, 3))) ! Only necessary if i = j

      ! First, calculate the pmat_resc1 = P^{alpha}/P^{beta} matrix (for atom i), the
      ! pmat_resc2 =  P^{alpha}/P^{beta} matrix (for atom j), and the
      ! pmat_resc3 = (beta - alpha)*P^{alpha}/P^{beta} matrix (added when i = j)

      pmat_resc1 = czero
      pmat_resc2 = czero
      pmat_resc3 = czero

      do s = 1, size(screening_in, 2) ! spin index
         do l = 0, size(screening_in, 1) - 1 ! lmax
            do m = 1, 2*l + 1
               mls = l*l + m + ((size(screening_in, 1))**2)*(s - 1) ! Composed diagonal index
               ! The screening constants can be of any atoms (i or j), because they will act only
               ! when i = j.
               temp1 = cmplx(screening_in(l, s), 0.0_rp) ! Screening constants (in)
               temp2 = cmplx(screening_out(l, s), 0.0_rp) ! Screening constants (out)
               do ie = 1, size(pmat_in_atom_i, 3) ! Energy channel
                  pmat_resc1(mls, mls, ie) = pmat_in_atom_i(mls, mls, ie)/pmat_out_atom_i(mls, mls, ie)
                  pmat_resc2(mls, mls, ie) = pmat_in_atom_j(mls, mls, ie)/pmat_out_atom_j(mls, mls, ie)
                  ! For the pmat_resc3 matrix, does not matter to consider the P potential function
                  ! for atom i or j, because it´s only acting when i = j
                  if (atom_i .eq. atom_j) then
                     pmat_resc3(mls, mls, ie) = (temp2 - temp1)*(pmat_in_atom_i(mls, mls, ie)/pmat_out_atom_i(mls, mls, ie))
                  end if
               end do
            end do
         end do
      end do

      ! Now do the rescaling of the auxiliary GF´s (in)

      do ie = 1, size(aux_gij_in, 3) ! Energy channel
         aux_gij_out(:, :, ie) = matmul(pmat_resc1(:, :, ie), matmul(aux_gij_in(:, :, ie), pmat_resc2(:, :, ie))) + &
                                 pmat_resc3(:, :, ie) ! If i != j, then the sum 0.
      end do

      deallocate (pmat_resc1)
      deallocate (pmat_resc2)
      deallocate (pmat_resc3)

   end subroutine transform_auxiliary_gij

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the Greens function
   !---------------------------------------------------------------------------

end module
