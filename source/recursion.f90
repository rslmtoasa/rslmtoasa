!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Recursion
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
!> Module to handle recursion calculations
!------------------------------------------------------------------------------

module recursion_mod

   use hamiltonian_mod
   use lattice_mod
   use control_mod
   use energy_mod
   use precision_mod, only: rp
   use math_mod
   use string_mod
   use logger_mod, only: g_logger
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc
#endif
   use timer_mod, only: g_timer
   implicit none

   private

   !> Module´s main structure
   type, public :: recursion
      !> Hamiltonian
      class(hamiltonian), pointer :: hamiltonian
      !> Lattice
      class(lattice), pointer :: lattice
      !> Energy
      class(energy), pointer :: en
      !> Control
      class(control), pointer :: control
      ! General variables

      !> Scalar recursion coefficients
      real(rp), dimension(:, :, :, :), allocatable :: a, b2
      real(rp), dimension(:), allocatable :: atemp, b2temp
      !> Block recursion coefficients
      complex(rp), dimension(:, :, :, :), allocatable :: a_b, b2_b
      complex(rp), dimension(:, :, :), allocatable :: atemp_b, b2temp_b
      !> Atom list in hopping region
      integer, dimension(:), allocatable :: izero, idum, irlist
      integer :: irnum !< Atoms within recursion shells
      !> Atom list in hopping region as a function of recursion step
      integer, dimension(:, :), allocatable :: izeroll
      !> Wave functions for recursion hopping (Lanczos)
      complex(rp), dimension(:, :), allocatable :: psi, pmn
      !> Wave functions for block recursion hopping (Lanczos)
      complex(rp), dimension(:, :, :), allocatable :: psi_b, pmn_b, hpsi, hohpsi, enupsi, socpsi
      !> Wave functions for recursion hopping (Chebyshev)
      complex(rp), dimension(:, :, :), allocatable :: psi0, psi1, psi2
      !> Chebyshev moments
      complex(rp), dimension(:, :, :, :), allocatable :: mu_n, mu_ng
      complex(rp), dimension(18, 18) :: cheb_mom_temp
      !> Variable to save H|psi>
      complex(rp), dimension(:, :), allocatable :: v
   contains
      procedure :: hop
      procedure :: crecal
      procedure :: recur
      procedure :: hop_b
      procedure :: hop_b_hoh
      procedure :: crecal_b
      procedure :: recur_b
    !!! procedure :: recur_b_par
      procedure :: recur_b_ij
      procedure :: zsqr
      procedure :: bpopt
      procedure :: get_terminf
      procedure :: get_cinf
      procedure :: emami
      procedure :: create_ll_map
      procedure :: cheb_0th_mom
      procedure :: cheb_1st_mom
      procedure :: cheb_1st_mom_hoh
      procedure :: chebyshev_recur_ij
      procedure :: chebyshev_recur_ll
      procedure :: chebyshev_recur_ll_hoh
      procedure :: chebyshev_recur_s_ll
      procedure :: chebyshev_recur
      procedure :: chebyshev_recur_full
      procedure :: restore_to_default
      final :: destructor
   end type recursion

   interface recursion
      procedure :: constructor
   end interface recursion

contains

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Constructor
   !
   !> @param[in] fname Namelist file
   !> @return type(recursion)
   !---------------------------------------------------------------------------
   function constructor(hamiltonian_obj, energy_obj) result(obj)
      type(recursion) :: obj
      type(hamiltonian), target, intent(in) :: hamiltonian_obj
      type(energy), target, intent(in) :: energy_obj

      obj%hamiltonian => hamiltonian_obj
      obj%lattice => hamiltonian_obj%charge%lattice
      obj%en => energy_obj
      obj%control => hamiltonian_obj%charge%lattice%control

      call obj%restore_to_default()
   end function constructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Destructor
   !---------------------------------------------------------------------------
   subroutine destructor(this)
      type(recursion) :: this
#ifdef USE_SAFE_ALLOC
      if (allocated(this%a)) call g_safe_alloc%deallocate('recursion.a', this%a)
      if (allocated(this%b2)) call g_safe_alloc%deallocate('recursion.b2', this%b2)
      if (allocated(this%izero)) call g_safe_alloc%deallocate('recursion.izero', this%izero)
      if (allocated(this%izeroll)) call g_safe_alloc%deallocate('recursion.izeroll', this%izeroll)
      if (allocated(this%irlist)) call g_safe_alloc%deallocate('recursion.irlist', this%irlist)
      if (allocated(this%psi)) call g_safe_alloc%deallocate('recursion.psi', this%psi)
      if (allocated(this%psi1)) call g_safe_alloc%deallocate('recursion.psi1', this%psi1)
      if (allocated(this%psi2)) call g_safe_alloc%deallocate('recursion.psi2', this%psi2)
      if (allocated(this%pmn)) call g_safe_alloc%deallocate('recursion.pmn', this%pmn)
      if (allocated(this%v)) call g_safe_alloc%deallocate('recursion.v', this%v)
      if (allocated(this%mu_n)) call g_safe_alloc%deallocate('recursion.mu_n', this%mu_n)
      if (allocated(this%mu_ng)) call g_safe_alloc%deallocate('recursion.mu_ng', this%mu_ng)
      if (allocated(this%atemp)) call g_safe_alloc%deallocate('recursion.atemp', this%atemp)
      if (allocated(this%b2temp)) call g_safe_alloc%deallocate('recursion.b2temp', this%b2temp)
      if (allocated(this%a_b)) call g_safe_alloc%deallocate('recursion.a_b', this%a_b)
      if (allocated(this%b2_b)) call g_safe_alloc%deallocate('recursion.b2_b', this%b2_b)
      if (allocated(this%psi_b)) call g_safe_alloc%deallocate('recursion.psi_b', this%psi_b)
      if (allocated(this%hpsi)) call g_safe_alloc%deallocate('recursion.hpsi', this%hpsi)
      if (allocated(this%hohpsi)) call g_safe_alloc%deallocate('recursion.hohpsi', this%hohpsi)
      if (allocated(this%enupsi)) call g_safe_alloc%deallocate('recursion.enupsi', this%enupsi)
      if (allocated(this%socpsi)) call g_safe_alloc%deallocate('recursion.socpsi', this%socpsi)
      if (allocated(this%atemp_b)) call g_safe_alloc%deallocate('recursion.atemp_b', this%atemp_b)
      if (allocated(this%b2temp_b)) call g_safe_alloc%deallocate('recursion.b2temp_b', this%b2temp_b)
      if (allocated(this%pmn_b)) call g_safe_alloc%deallocate('recursion.pmn_b', this%pmn_b)
#else
      if (allocated(this%a)) deallocate (this%a)
      if (allocated(this%b2)) deallocate (this%b2)
      if (allocated(this%izero)) deallocate (this%izero)
      if (allocated(this%izeroll)) deallocate (this%izeroll)
      if (allocated(this%irlist)) deallocate (this%irlist)
      if (allocated(this%psi)) deallocate (this%psi)
      if (allocated(this%psi1)) deallocate (this%psi1)
      if (allocated(this%psi2)) deallocate (this%psi2)
      if (allocated(this%pmn)) deallocate (this%pmn)
      if (allocated(this%v)) deallocate (this%v)
      if (allocated(this%mu_n)) deallocate (this%mu_n)
      if (allocated(this%mu_ng)) deallocate (this%mu_ng)
      if (allocated(this%atemp)) deallocate (this%atemp)
      if (allocated(this%b2temp)) deallocate (this%b2temp)
      if (allocated(this%a_b)) deallocate (this%a_b)
      if (allocated(this%b2_b)) deallocate (this%b2_b)
      if (allocated(this%psi_b)) deallocate (this%psi_b)
      if (allocated(this%hpsi)) deallocate (this%hpsi)
      if (allocated(this%hohpsi)) deallocate (this%hohpsi)
      if (allocated(this%enupsi)) deallocate (this%enupsi)
      if (allocated(this%socpsi)) deallocate (this%socpsi)
      if (allocated(this%atemp_b)) deallocate (this%atemp_b)
      if (allocated(this%b2temp_b)) deallocate (this%b2temp_b)
      if (allocated(this%pmn_b)) deallocate (this%pmn_b)
#endif
   end subroutine destructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculate the recursion coefficients a_b (block recursion). Considers the
   ! hoh part of the Hamiltonian
   !---------------------------------------------------------------------------
   subroutine hop_b_hoh(this, ll)
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m, n, nr, nnmap, nlimplus1
      integer :: ll ! Recursion step
      integer :: ino ! Atom type
      integer :: ih ! Atom number in the clust
      integer, dimension(0:this%lattice%kk) :: idum
      complex(rp), dimension(18, 18) :: summ
      complex(rp), dimension(18, 18) :: locham

      idum(:) = 0
      this%hpsi(:, :, :) = (0.0d0, 0.0d0)
      this%hohpsi(:, :, :) = (0.0d0, 0.0d0)
      this%enupsi(:, :, :) = (0.0d0, 0.0d0)
      this%socpsi(:, :, :) = (0.0d0, 0.0d0)

      nlimplus1 = this%lattice%nmax + 1
      if (this%lattice%nmax /= 0) then ! In case of impurities using the local hamiltonian
         !$omp parallel do default(shared) private(i, ino, nr, j, nnmap, locham) schedule(dynamic, 100)
         do i = 1, this%lattice%nmax ! Loop in the neighbouring
            idum(i) = this%izero(i)
            ino = this%lattice%iz(i)
            nr = this%lattice%nn(i, 1) ! Number of neighbours of atom i
            if (this%izero(i) /= 0) then
               call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%hall(:, :, 1, i), 18, this%psi_b(:, :, i), 18, cone, this%hpsi(:, :, i), 18)
               call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%enim(:, :, ino), 18, this%psi_b(:, :, i), 18, cone, this%enupsi(:, :, i), 18)
               call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%lsham(:, :, ino), 18, this%psi_b(:, :, i), 18, cone, this%socpsi(:, :, i), 18)
            end if
            if (nr >= 2) then
               do j = 2, nr ! Loop on the neighbouring
                  nnmap = this%lattice%nn(i, j)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%hall(:, :, j, i), 18, this%psi_b(:, :, nnmap), 18, cone, this%hpsi(:, :, i), 18)
                        idum(i) = 1
                     end if
                  end if
               end do ! End of loop in the neighbouring
            end if
         end do ! End of loop in the neighbouring
         !$omp end parallel do
      end if ! End of local Hamiltonian loop

      !$omp parallel do default(shared) private(i, ih, nr, j, nnmap, locham) schedule(dynamic, 100)
      do i = nlimplus1, this%lattice%kk ! Loop to find the bulk atoms using the bulk Hamiltonian
         idum(i) = this%izero(i)
         ih = this%lattice%iz(i) ! Atom type
         nr = this%lattice%nn(i, 1) ! Number of neighbours
         if (this%izero(i) /= 0) then
            call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%ee(:, :, 1, ih), 18, this%psi_b(:, :, i), 18, cone, this%hpsi(:, :, i), 18)
            call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%enim(:, :, ih), 18, this%psi_b(:, :, i), 18, cone, this%enupsi(:, :, i), 18)
            call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%lsham(:, :, ih), 18, this%psi_b(:, :, i), 18, cone, this%socpsi(:, :, i), 18)
         end if
         if (nr >= 2) then
            do j = 2, nr ! Loop on the neighbouring
               nnmap = this%lattice%nn(i, j)
               if (nnmap /= 0) then
                  if (this%izero(nnmap) /= 0) then
                     call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%ee(:, :, j, ih), 18, this%psi_b(:, :, nnmap), 18, cone, this%hpsi(:, :, i), 18)
                     idum(i) = 1
                  end if
               end if
            end do ! End of loop in the neighbouring
         end if
      end do
      !$omp end parallel do

      ! Mapping update for hoh calculation
      this%izero(:) = idum(:)

      if (this%lattice%nmax /= 0) then ! In case of impurities using the local hamiltonian
         !$omp parallel do default(shared) private(i, ino, nr, j, nnmap, locham) schedule(dynamic, 100)
         do i = 1, this%lattice%nmax ! Loop in the neighbouring
            idum(i) = this%izero(i)
            ino = this%lattice%iz(i)
            nr = this%lattice%nn(i, 1) ! Number of neighbours of atom i
            if (this%izero(i) /= 0) then
               call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%hallo(:, :, 1, i), 18, this%hpsi(:, :, i), 18, cone, this%hohpsi(:, :, i), 18)
            end if
            if (nr >= 2) then
               do j = 2, nr ! Loop on the neighbouring
                  nnmap = this%lattice%nn(i, j)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%hallo(:, :, j, i), 18, this%hpsi(:, :, nnmap), 18, cone, this%hohpsi(:, :, i), 18)
                        idum(i) = 1
                     end if
                  end if
               end do ! End of loop in the neighbouring
            end if
         end do ! End of loop in the neighbouring
         !$omp end parallel do
      end if ! End of local Hamiltonian loop

      !$omp parallel do default(shared) private(i, ih, nr, j, nnmap, locham) schedule(dynamic, 100)
      do i = nlimplus1, this%lattice%kk ! Loop to find the bulk atoms using the bulk Hamiltonian
         idum(i) = this%izero(i)
         ih = this%lattice%iz(i) ! Atom type
         nr = this%lattice%nn(i, 1) ! Number of neighbours
         if (this%izero(i) /= 0) then
            call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%eeo(:, :, 1, ih), 18, this%hpsi(:, :, i), 18, cone, this%hohpsi(:, :, i), 18)
         end if
         if (nr >= 2) then
            do j = 2, nr ! Loop on the neighbouring
               nnmap = this%lattice%nn(i, j)
               if (nnmap /= 0) then
                  if (this%izero(nnmap) /= 0) then
                     call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%eeo(:, :, j, ih), 18, this%hpsi(:, :, nnmap), 18, cone, this%hohpsi(:, :, i), 18)
                     idum(i) = 1
                  end if
               end if
            end do ! End of loop in the neighbouring
         end if
      end do
      !$omp end parallel do

      summ = 0.0d0

      ! Redefines idum for all clust atoms
      this%irnum = 0
      do i = 1, this%lattice%kk
         this%izero(i) = idum(i)
         if (this%izero(i) /= 0) then
            this%irnum = this%irnum + 1
            this%irlist(this%irnum) = i
         end if
      end do

      !$omp parallel do default(shared) private(i) reduction(+:SUMM) schedule(dynamic, 100)
      do i = 1, this%irnum
         ! H = h - hoh + e_nu + l.s
         this%hpsi(:, :, this%irlist(i)) = this%hpsi(:, :, this%irlist(i)) - this%hohpsi(:, :, this%irlist(i)) + this%enupsi(:, :, this%irlist(i)) + this%socpsi(:, :, this%irlist(i))
         !
         this%pmn_b(:, :, this%irlist(i)) = this%hpsi(:, :, this%irlist(i)) - this%pmn_b(:, :, this%irlist(i))
         !
         call zgemm('c', 'n', 18, 18, 18, cone, this%psi_b(:, :, this%irlist(i)), 18, this%hpsi(:, :, this%irlist(i)), 18, cone, summ, 18)
      end do
      !$omp end parallel do

      this%atemp_b(:, :, ll) = summ
   end subroutine
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculate the recursion coefficients a_b (block recursion)
   !---------------------------------------------------------------------------
   subroutine hop_b(this, ll)
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m, n, nr, nnmap, nlimplus1
      integer :: ll ! Recursion step
      integer :: ino ! Atom type
      integer :: ih ! Atom number in the clust
      integer, dimension(0:this%lattice%kk) :: idum
      complex(rp), dimension(18, 18) :: summ
      complex(rp), dimension(18, 18) :: locham

      idum(:) = 0
      this%hpsi(:, :, :) = (0.0d0, 0.0d0)

      nlimplus1 = this%lattice%nmax + 1
      if (this%lattice%nmax /= 0) then ! In case of impurities using the local hamiltonian
         !$omp parallel do default(shared) private(i, ino, nr, j, nnmap, locham) schedule(dynamic, 100)
         do i = 1, this%lattice%nmax ! Loop in the neighbouring
            idum(i) = this%izero(i)
            ino = this%lattice%iz(i)
            nr = this%lattice%nn(i, 1) ! Number of neighbours of atom i
            if (this%izero(i) /= 0) then
               locham = this%hamiltonian%hall(:, :, 1, i) + this%hamiltonian%lsham(:, :, ino)
               call zgemm('n', 'n', 18, 18, 18, cone, locham, 18, this%psi_b(:, :, i), 18, cone, this%hpsi(:, :, i), 18)
               !call zgemm(´n´,´n´,18,18,18,cone,this%hamiltonian%hall(:,:,1,i),18,this%psi_b(:,:,i),18,cone,this%hpsi(:,:,i),18)
               !call zgemm(´n´,´n´,18,18,18,cone,this%hamiltonian%lsham(:,:,ino),18,this%psi_b(:,:,i),18,cone,this%hpsi(:,:,I),18)
            end if
            if (nr >= 2) then
               do j = 2, nr ! Loop on the neighbouring
                  nnmap = this%lattice%nn(i, j)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%hall(:, :, j, i), 18, this%psi_b(:, :, nnmap), 18, cone, this%hpsi(:, :, i), 18)
                        idum(i) = 1
                     end if
                  end if
               end do ! End of loop in the neighbouring
            end if
         end do ! End of loop in the neighbouring
         !$omp end parallel do
      end if ! End of local Hamiltonian loop

      !$omp parallel do default(shared) private(i, ih, nr, j, nnmap, locham) schedule(dynamic, 100)
      do i = nlimplus1, this%lattice%kk ! Loop to find the bulk atoms using the bulk Hamiltonian
         idum(i) = this%izero(i)
         ih = this%lattice%iz(i) ! Atom type
         nr = this%lattice%nn(i, 1) ! Number of neighbours
         if (this%izero(i) /= 0) then
            locham = this%hamiltonian%ee(:, :, 1, ih) + this%hamiltonian%lsham(:, :, ih)
            call zgemm('n', 'n', 18, 18, 18, cone, locham, 18, this%psi_b(:, :, i), 18, cone, this%hpsi(:, :, i), 18)
            !call zgemm(´n´,´n´,18,18,18,cone,this%hamiltonian%ee(:,:,1,ih),18,this%psi_b(:,:,i),18,cone,this%hpsi(:,:,i),18)
            !call zgemm(´n´,´n´,18,18,18,cone,this%hamiltonian%lsham(:,:,ih),18,this%psi_b(:,:,i),18,cone,this%hpsi(:,:,I),18)
         end if
         if (nr >= 2) then
            do j = 2, nr ! Loop on the neighbouring
               nnmap = this%lattice%nn(i, j)
               if (nnmap /= 0) then
                  if (this%izero(nnmap) /= 0) then
                     call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%ee(:, :, j, ih), 18, this%psi_b(:, :, nnmap), 18, cone, this%hpsi(:, :, i), 18)
                     idum(i) = 1
                  end if
               end if
            end do ! End of loop in the neighbouring
         end if
      end do
      !$omp end parallel do

      summ = 0.0d0
      ! Redefines idum for all clust atoms
      this%irnum = 0
      do i = 1, this%lattice%kk
         this%izero(i) = idum(i)
         if (this%izero(i) /= 0) then
            this%irnum = this%irnum + 1
            this%irlist(this%irnum) = i
         end if
      end do

      !$omp parallel do default(shared) private(i) reduction(+:SUMM) schedule(dynamic, 100)
      do i = 1, this%irnum
         !do i=1, this%lattice%kk
         this%pmn_b(:, :, this%irlist(i)) = this%hpsi(:, :, this%irlist(i)) - this%pmn_b(:, :, this%irlist(i))
         call zgemm('c', 'n', 18, 18, 18, cone, this%psi_b(:, :, this%irlist(i)), 18, this%hpsi(:, :, this%irlist(i)), 18, cone, summ, 18)
         !call zgemm(´c´,´n´,18,18,18,cone,PSI_B(1,1,irlist(I)),18,Hpsi(1,1,irlist(I)),18,cone,SUMM,18)
      end do
      !$omp end parallel do

      this%atemp_b(:, :, ll) = summ
   end subroutine hop_b

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Perform the block recursion between sites i and j
   !---------------------------------------------------------------------------
   subroutine recur_b_ij(this)
      use mpi_mod
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, ij, j, l, ll, kk, m, reci
      integer :: llmax ! Recursion steps
      complex(rp) :: asign, bsign

      integer :: ij_loc

      llmax = this%lattice%control%lld
      !do ij=1, this%lattice%njij ! Loop on the number of pair of atoms
      do ij = start_atom, end_atom
         ij_loc = g2l_map(ij)
         i = this%lattice%ijpair(ij, 1) ! Atom number in the clust file, atom i
         j = this%lattice%ijpair(ij, 2) ! Atom number in the clust file, atom j
         call g_logger%info('Block recursion on progress between atoms '//int2str(i)//' and '//int2str(j), __FILE__, __LINE__)
         do reci = 1, 4
            ! Clear list of atoms in hopping region
            this%izero(:) = 0
            ! Initializing wave functions
            this%psi_b(:, :, :) = (0.0d0, 0.0d0)
            this%pmn_b(:, :, :) = (0.0d0, 0.0d0)

            select case (reci)
            case (1)
               asign = (1.0d0, 0.0d0)*one_over_sqrt_two
               bsign = (1.0d0, 0.0d0)*one_over_sqrt_two
               this%izero(i) = 1
               this%izero(j) = 1
            case (2)
               asign = (1.0d0, 0.0d0)*one_over_sqrt_two
               bsign = (-1.0d0, 0.0d0)*one_over_sqrt_two
               this%izero(i) = 1
               this%izero(j) = 1
            case (3)
               asign = (1.0d0, 0.0d0)*one_over_sqrt_two
               bsign = (0.0d0, 1.0d0)*one_over_sqrt_two
               this%izero(i) = 1
               this%izero(j) = 1
            case (4)
               asign = (1.0d0, 0.0d0)*one_over_sqrt_two
               bsign = (0.0d0, -1.0d0)*one_over_sqrt_two
               this%izero(i) = 1
               this%izero(j) = 1
            end select

            if ((reci .eq. 1) .and. (i .eq. j)) then
               asign = (1.0d0, 0.0d0)
               bsign = (1.0d0, 0.0d0)
            else if (i .eq. j) then
               cycle
            end if

            do l = 1, 18
               this%psi_b(l, l, i) = asign
               this%psi_b(l, l, j) = bsign
               this%atemp_b(l, l, llmax) = (0.0d0, 0.0d0)
               this%b2temp_b(l, l, 1) = (1.0d0, 0.0d0)
            end do

            call this%crecal_b()

            do ll = 1, llmax
               do l = 1, 18
                  do m = 1, 18
                     this%a_b(l, m, ll, ij_loc*4 - 4 + reci) = (this%atemp_b(l, m, ll))
                     this%b2_b(l, m, ll, ij_loc*4 - 4 + reci) = (this%b2temp_b(l, m, ll))
                  end do
               end do
            end do
         end do ! End of the loop on the basis
      end do ! End of the loop on the atom pairs
      ! For debug purposes
      !do i=1, this%lattice%nrec ! Loop on the number of atoms to be treat self-consistently
      !  do l=1, 18  ! Loop on the orbital l
      !    write(122, *) ´orbital´, l
      !    do ll=1, llmax ! Loop on the recursion steps
      !      write(122, *) this%a(ll, l, i,1), this%b2(ll,l,i,1)
      !    end do
      !  end do
      !end do
   end subroutine recur_b_ij

!!!   !----------------------------------------------------------------------------
!!!   ! DESCRIPTION:
!!!   !> @brief
!!!   !> Perform the block recursion in parallel over atoms
!!!   !----------------------------------------------------------------------------
!!!   subroutine recur_b_par(this)
!!!     use mpi_mod
!!!     class(recursion), intent(inout) :: this
!!!     ! Local variables
!!!     integer :: i, j, l, ll, m
!!!     integer :: llmax          ! Recursion steps
!!!
!!!     llmax = this%lattice%control%lld
!!!
!!!     ! Determine how many atoms each process should handle
!!!     call get_mpi_variables(rank,this%lattice%nrec)
!!!
!!!     ! Loop on the number of atoms to be treated self-consistently by this process
!!!     do i=start_atom, end_atom
!!!         j = this%lattice%irec(i) ! Atom number in the clust file
!!!         ! Clear list of atoms in hopping region
!!!         call g_logger%info(´Block recursion in progress for atom ´//int2str(j),__FILE__,__LINE__)
!!!         this%izero(:) = 0
!!!         ! Initializing wave functions
!!!         this%psi_b(:,:,:) = (0.0d0, 0.0d0)
!!!         this%pmn_b(:,:,:) = (0.0d0, 0.0d0)
!!!
!!!         do l=1,18
!!!             this%psi_b(l,l,j) = (1.0d0, 0.0d0)
!!!             this%atemp_b(l,l,llmax) = (0.0d0,0.0d0)
!!!             this%b2temp_b(l,l,1) = (1.0d0,0.0d0)
!!!         end do
!!!
!!!         this%izero(j) = 1
!!!
!!!         call this%crecal_b()
!!!
!!!         do ll= 1,llmax
!!!            do l = 1,18
!!!               do m = 1,18
!!!                  this%a_b(l,m,ll,i) = (this%atemp_b(l,m,ll))
!!!                  this%b2_b(l,m,ll,i) = (this%b2temp_b(l,m,ll))
!!!               end do
!!!               this%a(ll,l,i,1) = real(this%atemp_b(l,l,ll))
!!!               this%b2(ll,l,i,1) = real(this%b2temp_b(l,l,ll))
!!!            end do
!!!         end do
!!!     end do ! End of the loop on the nrec
!!!     ! Gather results to ensure all processes have the complete arrays
!!! #ifdef USE_MPI
!!!     call g_timer%start(´MPI recursion communication´)
!!!     call MPI_Allgather(this%a_b(:,:,:,start_atom:end_atom), (end_atom-start_atom+1)*llmax*18*18, MPI_DOUBLE_COMPLEX, &
!!!                        this%a_b, (end_atom-start_atom+1)*llmax*18*18, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
!!!     call MPI_Allgather(this%b2_b(:,:,:,start_atom:end_atom), (end_atom-start_atom+1)*llmax*18*18, MPI_DOUBLE_COMPLEX, &
!!!                        this%b2_b, (end_atom-start_atom+1)*llmax*18*18, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
!!!     !!! call MPI_Allgather(this%a_b(:,:,:,start_atom:end_atom), (end_atom-start_atom+1)*llmax*18*18, MPI_DOUBLE_COMPLEX, &
!!!     !!!                    this%a_b, (end_atom-start_atom+1)*llmax*18*18, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
!!!     !!! call MPI_Allgather(this%b2_b(:,:,:,start_atom:end_atom), (end_atom-start_atom+1)*llmax*18*18, MPI_DOUBLE_COMPLEX, &
!!!     !!!                    this%b2_b, (end_atom-start_atom+1)*llmax*18*18, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
!!!     call g_timer%stop(´MPI recursion communication´)
!!! #endif
!!!   end subroutine recur_b_par

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Perform the block recursion
   !---------------------------------------------------------------------------
   subroutine recur_b(this)
      use mpi_mod
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, l, ll, kk, m
      integer :: llmax ! Recursion steps

      ! Determine how many atoms each process should handle
      call get_mpi_variables(rank, this%lattice%nrec)

      llmax = this%lattice%control%lld
      !do i=1, this%lattice%nrec ! Loop on the number of atoms to be treat self-consistently
      do i = start_atom, end_atom ! Loop on the number of atoms to be treat self-consistently
         j = this%lattice%irec(i) ! Atom number in the clust file
         ! Clear list of atoms in hopping region
         call g_logger%info('Block recursion on progress for atom '//int2str(j), __FILE__, __LINE__)
         this%izero(:) = 0
         ! Initializing wave functions
         this%psi_b(:, :, :) = (0.0d0, 0.0d0)
         this%pmn_b(:, :, :) = (0.0d0, 0.0d0)

         ! Rotate Hamiltonian if needed
         if (this%hamiltonian%local_axis) then
            call this%hamiltonian%rotate_to_local_axis(this%lattice%symbolic_atoms(i)%potential%mom)
         end if

         do l = 1, 18
            this%psi_b(l, l, j) = (1.0d0, 0.0d0)
            this%atemp_b(l, l, llmax) = (0.0d0, 0.0d0)
            this%b2temp_b(l, l, 1) = (1.0d0, 0.0d0)
         end do

         this%izero(j) = 1

         call this%crecal_b()

         do ll = 1, llmax
            do l = 1, 18
               do m = 1, 18
                  this%a_b(l, m, ll, i - start_atom + 1) = (this%atemp_b(l, m, ll))
                  this%b2_b(l, m, ll, i - start_atom + 1) = (this%b2temp_b(l, m, ll))
               end do
               this%a(ll, l, i - start_atom + 1, 1) = real(this%atemp_b(l, l, ll))
               this%b2(ll, l, i - start_atom + 1, 1) = real(this%b2temp_b(l, l, ll))
            end do
         end do
      end do ! End of the loop on the nrec
      ! For debug purposes
      do i = start_atom, end_atom
         do l = 1, 18  ! Loop on the orbital l
            write (1000*(1 + rank) + 122, *) 'orbital', l, 'atom', i
            do ll = 1, llmax ! Loop on the recursion steps
               write (1000*(1 + rank) + 122, '(2f12.8,4x,2f12.8)') this%a(ll, l, i - start_atom + 1, 1), this%b2(ll, l, i - start_atom + 1, 1)
               !write(1000*(1+rank)+122, ´(2f12.8,4x,2f12.8)´) this%a_b(l,l,ll,i-start_atom+1), this%b2_b(l,l,ll,i-start_atom+1)
          !!!write(1000*(1+rank)+122, *) this%a(ll, l, i,1-start_atom+1), this%b2(ll,l,i,1)
            end do
         end do
      end do
   end subroutine recur_b

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the recursion coefficients b2_b (block recursion)
   !---------------------------------------------------------------------------
   subroutine crecal_b(this)
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m, info, lwork
      integer :: ll, llmax ! Recursion step
      integer :: nm1 ! LL-1
      integer :: nat ! Clust size
      integer :: hblocksize ! Hamiltonian size (18)
      complex(rp), dimension(18, 18) :: sum_b, sum_a
      complex(rp), dimension(18, 18) :: dum, b, b_i, u, lam, lam_i, u_ls, u_rs
!    complex(rp), dimension(18,18,this%lattice%kk) :: psi_t
      complex(rp), dimension(:, :, :), allocatable :: psi_t
      real(rp), dimension(18) :: ev
      complex(rp), dimension(18) :: zev
      real(rp), dimension(3*18 - 2) ::rwork
      complex(rp), dimension(3*18 - 2) ::zwork

      allocate (psi_t(18, 18, this%lattice%kk))

      sum_b = this%b2temp_b(:, :, 1)

      hblocksize = 18
      nat = this%lattice%kk
      nm1 = this%lattice%control%lld - 1
      llmax = this%lattice%control%lld
      do ll = 1, nm1
         this%atemp_b(:, :, ll) = (0.0d0, 0.0d0)
         !call HOP to get H*w = A; psi=psi ; pmn=Hpsi - B_n-1|Psi_n-1>
         !PMN = H Wn - W_(n-1)B_(n-1) in Inoue syntax  (is zero for first recursion step)
         call g_timer%start('H|PSI_n>')
         if (this%hamiltonian%hoh) then
            call this%hop_b_hoh(ll) ! a_n=<PSI_n|H|PSI_n>
         else
            call this%hop_b(ll)     ! a_n=<PSI_n|H|PSI_n>
         end if
         call g_timer%stop('H|PSI_n>')

         !psi_t(:,:,:) = this%psi_b(:,:,:)
         !call zcopy(18*18*this%lattice%kk,this%psi_b,1,psi_t,1)
         call g_timer%start('H|Psi_n-A_n|Psi_n-B_n|Psi_n-1')
         do i = 1, this%irnum
            psi_t(:, :, this%irlist(i)) = this%psi_b(:, :, this%irlist(i))
         end do
         !psi_t(:,:,:) = this%psi_b(:,:,:)

         this%b2temp_b(:, :, ll) = sum_b(:, :)

         sum_b(:, :) = 0.0d0

         !$omp parallel do default(shared) private(i) reduction(+:SUM_b)
         do i = 1, this%irnum
            !do i=1, nat
            !  ( H|Psi_n>-A_n|Psi_n>-B_n-1|Psi_n-1> = PMN-A_n|Psi_n> )
         !! H Wn - Wn An  - Wn-1 Bn-1 in Inoue syntax
            call zgemm('n', 'n', 18, 18, 18, cmone, this%psi_b(:, :, this%irlist(i)), 18, this%atemp_b(:, :, ll), 18, cone, this%pmn_b(:, :, this%irlist(i)), 18)
            !call zgemm(´n´,´n´,18,18,18,cmone,this%psi_b(:,:,i),18,this%atemp_b(:,:,ll),18,cone,this%pmn_b(:,:,i), 18)
            ! B^2_n+1 = ( H|Psi_n-A_n|Psi_n-B_n|Psi_n-1 )´*( H|Psi_n-A_n|Psi_n-B_n|Psi_n-1 )
         !! (H Wn - An Wn  - Wn-1 Bn-1)*(H Wn - An Wn - Wn-1 Bn-1)´ -in Inoue syntax
            call zgemm('c', 'n', 18, 18, 18, cone, this%pmn_b(:, :, this%irlist(i)), 18, this%pmn_b(:, :, this%irlist(i)), 18, cone, sum_b, 18)
            !call zgemm(´c´,´n´,18,18,18,cone,this%pmn_b(:,:,i),18,this%pmn_b(:,:,i),18,cone,sum_b,18)
         end do
         !$omp end parallel do
         call g_timer%stop('H|Psi_n-A_n|Psi_n-B_n|Psi_n-1')
         call g_timer%start('B_n+1')
         ! Calculate B_n+1
         u(:, :) = sum_b(:, :)
         ! Replace sqrt with eigen solver to get lamda^2 and U
         ! get lamda^2 and U in lamda=U*BB´*U*
         call zheev('v', 'u', 18, u, 18, ev, dum, 18*18, rwork, info)
         if (info /= 0) call g_logger%fatal('Diagonalization error', __FILE__, __LINE__)
         !
         !
         !

         lam = (0.0d0, 0.0d0)
         lam_i = (0.0d0, 0.0d0)
         do i = 1, 18
            lam(i, i) = cmplx(sqrt(ev(i)), kind=kind(0.0d0))
            lam_i(i, i) = 1.0d0/lam(i, i) ! 1/B_n
         end do
         !
         !  Calc. U*lamda*U´= B
         call zgemm('n', 'n', 18, 18, 18, cone, u, 18, lam, 18, czero, dum, 18)
         call zgemm('n', 'c', 18, 18, 18, cone, dum, 18, u, 18, czero, b, 18)
         !  B^-1=U*lamda^-1*U´
         call zgemm('n', 'n', 18, 18, 18, cone, u, 18, lam_i, 18, czero, dum, 18)
         call zgemm('n', 'c', 18, 18, 18, cone, dum, 18, u, 18, czero, b_i, 18)
         call g_timer%stop('B_n+1')

         call g_timer%start('<PSI|B_n+1|PSI>')
         !$omp parallel do default(shared) private(i)
         do i = 1, this%irnum
            !do i=1, nat
            call zgemm('n', 'n', 18, 18, 18, cone, this%pmn_b(:, :, this%irlist(i)), 18, b_i, 18, czero, this%psi_b(:, :, this%irlist(i)), 18)
            call zgemm('n', 'n', 18, 18, 18, cone, psi_t(:, :, this%irlist(i)), 18, b, 18, czero, this%pmn_b(:, :, this%irlist(i)), 18)
         end do
         !$omp end parallel do
         call g_timer%stop('<PSI|B_n+1|PSI>')
      end do
      this%b2temp_b(:, :, llmax) = sum_b(:, :)
   end subroutine crecal_b

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculate the square root of B2
   !---------------------------------------------------------------------------
   subroutine zsqr(this)
      use mpi_mod
      implicit none
      class(recursion), intent(inout) :: this
      !
      integer :: na, LDIM, I, J, K, L, N, M, info, n_glob
      real(rp), dimension(18) :: ev
      real(rp), dimension(3*18 - 2) ::rwork
      complex(rp), dimension(18*18) ::zwork
      complex(rp), dimension(18, 18) :: B2t, U, DUM, lam

      !
      if (this%lattice%njij == 0) then
         na = atoms_per_process !this%lattice%nrec
      else
         na = atoms_per_process*4 !4*this%lattice%njij
      end if
      !call get_mpi_variables(rank,na)
      lam = 0.0d0
      LDIM = 18
      ! do n_glob=start_atom, end_atom
      !  N = g2l_map(n_glob)
      do N = 1, na
         do L = 1, this%lattice%control%lld
            do I = 1, LDIM
               do J = 1, LDIM
                  U(J, I) = this%B2_b(J, I, L, N)
               end do
            end do
            ! replace sqrt with eigen solver to get lamda^2 and U
            ! get lamda^2 and U in lamda=U*BB´*U*
            call zheev('V', 'U', LDIM, U, LDIM, ev, zwork, LDIM*LDIM, rwork, info)
            !
            if (info /= 0) print *, 'diag', info
            ! lam=sqrt(lamda^2) ; lam_i=1/sqrt(lamda^2)
            do I = 1, LDIM
               lam(I, I) = sqrt(ev(i))
            end do
            ! calc. U*lamda*U´= B
            call zgemm('n', 'n', LDIM, LDIM, LDIM, (1.0d0, 0.0d0), U, LDIM, lam, LDIM, (0.0d0, 0.0d0), DUM, LDIM)
            call zgemm('n', 'c', LDIM, LDIM, LDIM, (1.0d0, 0.0d0), DUM, LDIM, U, LDIM, (0.0d0, 0.0d0), this%B2_b(1, 1, L, N), LDIM)
         end do
      end do
   end subroutine zsqr

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> TO DO
   !---------------------------------------------------------------------------
   subroutine get_cinf(this, alpha, beta, ll, np, nw, alpha_inf, beta_inf)
      implicit none
      class(recursion), intent(inout) :: this
      !.. Formal Arguments ..
      integer, intent(in) :: np
      integer, intent(inout) :: ll, nw
      real(rp), dimension(np, ll), intent(inout) :: alpha
      real(rp), dimension(np, ll), intent(in) :: beta
      real(rp), dimension(np), intent(out) :: alpha_inf
      real(rp), dimension(np), intent(out) :: beta_inf
      !
      !.. Local Scalars ..
      integer :: icode, iii, k, l, ll1, nb, nbp1, nl, nt, eidx, ne, nq, ifail
      real(rp) :: a1, a2, emax, emin, eps, err, e_shift, pi
      complex(rp) :: g_e
      !
      !.. Local Arrays ..
      integer, dimension(np) :: nb2
      integer, dimension(nw) :: jc
      integer, dimension(nw) :: iwk
      real(rp), dimension(ll) :: aa, am, bb, bm, sqbb
      real(rp), dimension(10) :: edge, weight, width
      real(rp), dimension(200, 10) :: bwk
      ! real(rp), dimension(np,ll) :: am2,bm2
      real(rp), dimension(np, 10) :: edge2, width2, weight2
      real(rp), dimension(nw, 2, 5) :: work
      real(rp), dimension(2) :: bandedges

      !
      !.. External Calls ..
      !
      !.. External Functions ..
      !
      !.. Intrinsic Functions ..
      intrinsic MAX, MIN
      !
      ! ... Executable Statements ...
      !
      !**************************************************************************
      icode = 0
      !*******************************************************
      eps = 1.0d-14
      err = 0.00001d0
      nbp1 = 2   !( Nband+1)
      !*************************************************************************
      do nl = 1, np
         do l = 1, ll
            aa(l) = alpha(nl, l)
            bb(l) = beta(nl, l)!**2
         end do
         ll1 = ll
         am = 0.0d0; bm = 0.0d0; edge = 0.0d0; width = 0.0d0; weight = 0.0d0
         call this%bpOPT(ll, AA, BB, LL - 1, alpha_inf(nl), beta_inf(nl), ifail)
         edge(1) = alpha_inf(nl) - 2.0d0*beta_inf(nl)
         width(1) = 4.0d0*beta_inf(nl)
      end do
   end subroutine get_cinf
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculate band-dependent terminator coefficients a_inf and b_inf
   !---------------------------------------------------------------------------
   subroutine get_terminf(this, Acoef_b, B2coef_b, na, ll, ldim, nw, a_inf, b_inf, a_inf0, b_inf0)
      implicit none
      class(recursion), intent(inout) :: this
      integer, intent(in) :: na
      integer, intent(in) :: ll
      integer, intent(inout) :: nw
      integer, intent(in) :: ldim
      complex(rp), dimension(ldim, ldim, ll, na), intent(in) :: Acoef_b, B2coef_b
      real(rp), dimension(ldim, ldim, na), intent(out) :: a_inf, b_inf
      real(rp), dimension(na), intent(out) ::  a_inf0, b_inf0
      !
      integer :: n, i, j, ll_t
      complex(rp), dimension(ldim, ldim) :: MatIn, MatOut
      real(rp), dimension(ldim, ldim, ll) :: Acoef_r, B2coef_r
      !
      a_inf = 0.0d0; b_inf = 0.0d0
      do n = 1, na
         ll_t = ll
         do i = 1, ll
            MatIn = real(Acoef_b(:, :, i, n))
            Acoef_r(:, :, i) = MatIn
            MatIn = real(B2coef_b(:, :, i, n))
            B2coef_r(:, :, i) = MatIn
         end do
         call this%get_cinf(Acoef_r, B2coef_r, ll_t, ldim*ldim, nw, a_inf(:, :, n), b_inf(:, :, n))
         do j = 1, ldim
            do i = 1, ldim
               if (IsNaN(a_inf(i, j, n))) a_inf(i, j, n) = 0.0d0
               if (IsNaN(b_inf(i, j, n))) b_inf(i, j, n) = 0.0d0
            end do
            if (a_inf(j, j, n) == 0.0d0) a_inf(j, j, n) = 0.5d0
            if (b_inf(j, j, n) == 0.0d0) b_inf(j, j, n) = 0.5d0
         end do
         a_inf0(n) = 0.0d0
         do i = 1, ldim
            a_inf0(n) = a_inf0(n) + a_inf(i, i, n)
         end do
         a_inf0(n) = a_inf0(n)/ldim
         b_inf(1, 1, n) = b_inf(1, 1, n)*1.01d0
         b_inf(10, 10, n) = b_inf(10, 10, n)*1.01d0
         b_inf0(n) = 0.0d0
         do i = 1, ldim
            b_inf0(n) = b_inf0(n) + b_inf(i, i, n)
         end do
         b_inf0(n) = b_inf0(n)/ldim
      end do
   end subroutine get_terminf

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the 0th Chebyshev moment
   !---------------------------------------------------------------------------
   subroutine cheb_0th_mom(this, psiref)
      class(recursion), intent(inout) :: this
      complex(rp), dimension(18, 18, this%lattice%kk), intent(in) :: psiref
      complex(rp), dimension(18, 18) :: dum
      integer :: nat, hblocksize, nlimplus1, k

      hblocksize = 18
      nat = this%lattice%kk
      nlimplus1 = this%lattice%nmax + 1

      this%cheb_mom_temp(:, :) = 0.0d0
      ! Write the 0th moment
      do k = 1, nat
         if (this%izero(k) /= 0) then
            call zgemm('c', 'n', 18, 18, 18, cone, psiref(:, :, k), 18, this%psi0(:, :, k), 18, cone, this%cheb_mom_temp, 18)
         end if
      end do
   end subroutine cheb_0th_mom

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the 1st Chebyshev moment
   !---------------------------------------------------------------------------
   subroutine cheb_1st_mom(this, psiref, a, b)
      class(recursion), intent(inout) :: this
      real(rp), intent(in) :: a, b ! scale and shift
      complex(rp), dimension(18, 18, this%lattice%kk), intent(in) :: psiref
      integer :: nat, hblocksize, nlimplus1, k, ih, nr, n, nb, nnmap, i

      hblocksize = 18
      nat = this%lattice%kk
      nlimplus1 = this%lattice%nmax + 1

      this%cheb_mom_temp(:, :) = 0.0d0

      ! Write |phi_1>=H|phi_0>
      nlimplus1 = this%lattice%nmax + 1
      if (this%lattice%nmax /= 0) then ! In case of impurities using the local hamiltonian
         do i = 1, this%lattice%nmax ! Loop in the neighbouring
            this%idum(i) = this%izero(i)
            ih = this%lattice%iz(i)
            nr = this%lattice%nn(i, 1) ! Number of neighbours of atom i
            if (this%izero(i) /= 0) then
               call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%hall(:, :, 1, i), 18, this%psi0(:, :, i), 18, cone, this%psi1(:, :, i), 18)
               call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%lsham(:, :, ih), 18, this%psi0(:, :, i), 18, cone, this%psi1(:, :, i), 18)
            end if
            if (nr >= 2) then
               do nb = 2, nr ! Loop on the neighbouring
                  nnmap = this%lattice%nn(i, nb)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%hall(:, :, nb, i), 18, this%psi0(:, :, nnmap), 18, cone, this%psi1(:, :, i), 18)
                        this%idum(i) = 1
                     end if
                  end if
               end do ! End of loop in the neighbouring
            end if
            ! Do the scaling and shifting
            this%psi1(:, :, i) = this%psi1(:, :, i) - b*this%psi0(:, :, i)
            this%psi1(:, :, i) = this%psi1(:, :, i)/a
         end do ! End of loop in the neighbouring
      end if ! End of local Hamiltonian loop

      ! Write |phi_1>=H|phi_0>
      do k = nlimplus1, this%lattice%kk ! Loop in the clust
         this%idum(k) = this%izero(k)
         ih = this%lattice%iz(k)
         nr = this%lattice%nn(k, 1)
         if (this%izero(k) /= 0) then
            call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%ee(1, 1, 1, ih), 18, this%psi0(:, :, k), 18, cone, this%psi1(:, :, k), 18)
            call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%lsham(:, :, ih), 18, this%psi0(:, :, k), 18, cone, this%psi1(:, :, k), 18)
         end if
         if (nr >= 2) then
            do nb = 2, nr ! Loop in the neighbouring
               nnmap = this%lattice%nn(k, nb)
               if (nnmap /= 0 .and. this%izero(nnmap) /= 0) then
                  call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%ee(1, 1, nb, ih), 18, this%psi0(:, :, nnmap), 18, cone, this%psi1(:, :, k), 18)
                  this%idum(k) = 1
               end if
            end do ! End of the loop in the neighbouring
         end if
         ! Do the scaling and shifting
         this%psi1(:, :, k) = this%psi1(:, :, k) - b*this%psi0(:, :, k)
         this%psi1(:, :, k) = this%psi1(:, :, k)/a
      end do ! End loop in the clust

      ! Write the 1st moment
      do n = 1, this%lattice%kk
         if (this%izero(n) /= 0) then
            call zgemm('c', 'n', 18, 18, 18, cone, psiref(:, :, n), 18, this%psi1(:, :, n), 18, cone, this%cheb_mom_temp, 18)
         end if
      end do
   end subroutine cheb_1st_mom

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the 1st Chebyshev moment including hoh term
   !---------------------------------------------------------------------------
   subroutine cheb_1st_mom_hoh(this, psiref, a, b)
      class(recursion), intent(inout) :: this
      real(rp), intent(in) :: a, b ! scale and shift
      complex(rp), dimension(18, 18, this%lattice%kk), intent(in) :: psiref
      integer :: nat, hblocksize, nlimplus1, k, ih, nr, n, nb, nnmap, i

      hblocksize = 18
      nat = this%lattice%kk
      nlimplus1 = this%lattice%nmax + 1

      this%cheb_mom_temp(:, :) = 0.0d0

      this%hohpsi(:, :, :) = (0.0d0, 0.0d0)
      this%enupsi(:, :, :) = (0.0d0, 0.0d0)
      this%socpsi(:, :, :) = (0.0d0, 0.0d0)

      ! Write |phi_1>=H|phi_0>
      nlimplus1 = this%lattice%nmax + 1
      if (this%lattice%nmax /= 0) then ! In case of impurities using the local hamiltonian
         do i = 1, this%lattice%nmax ! Loop in the neighbouring
            this%idum(i) = this%izero(i)
            ih = this%lattice%iz(i)
            nr = this%lattice%nn(i, 1) ! Number of neighbours of atom i
            if (this%izero(i) /= 0) then
               call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%hall(:, :, 1, i), 18, this%psi0(:, :, i), 18, cone, this%psi1(:, :, i), 18)
               call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%lsham(:, :, ih), 18, this%psi0(:, :, i), 18, cone, this%socpsi(:, :, i), 18)
               call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%enim(:, :, ih), 18, this%psi0(:, :, i), 18, cone, this%enupsi(:, :, i), 18)
            end if
            if (nr >= 2) then
               do nb = 2, nr ! Loop on the neighbouring
                  nnmap = this%lattice%nn(i, nb)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%hall(:, :, nb, i), 18, this%psi0(:, :, nnmap), 18, cone, this%psi1(:, :, i), 18)
                        this%idum(i) = 1
                     end if
                  end if
               end do ! End of loop in the neighbouring
            end if
         end do ! End of loop in the neighbouring
      end if ! End of local Hamiltonian loop

      ! Write |phi_1>=H|phi_0>
      do k = nlimplus1, this%lattice%kk ! Loop in the clust
         this%idum(k) = this%izero(k)
         ih = this%lattice%iz(k)
         nr = this%lattice%nn(k, 1)
         if (this%izero(k) /= 0) then
            call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%ee(:, :, 1, ih), 18, this%psi0(:, :, k), 18, cone, this%psi1(:, :, k), 18)
            call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%enim(:, :, ih), 18, this%psi0(:, :, k), 18, cone, this%enupsi(:, :, k), 18)
            call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%lsham(:, :, ih), 18, this%psi0(:, :, k), 18, cone, this%socpsi(:, :, k), 18)
         end if
         if (nr >= 2) then
            do nb = 2, nr ! Loop in the neighbouring
               nnmap = this%lattice%nn(k, nb)
               if (nnmap /= 0 .and. this%izero(nnmap) /= 0) then
                  call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%ee(1, 1, nb, ih), 18, this%psi0(:, :, nnmap), 18, cone, this%psi1(:, :, k), 18)
                  this%idum(k) = 1
               end if
            end do ! End of the loop in the neighbouring
         end if
         ! Do the scaling and shifting
         !this%psi1(:, :, k) = this%psi1(:, :, k) - b*this%psi0(:, :, k)
         !this%psi1(:, :, k) = this%psi1(:, :, k)/a
      end do ! End loop in the clust

      ! Mapping update for hoh calculation
      this%izero(:) = this%idum(:)

      if (this%lattice%nmax /= 0) then ! In case of impurities using the local hamiltonian
         do i = 1, this%lattice%nmax ! Loop in the neighbouring
            this%idum(i) = this%izero(i)
            ih = this%lattice%iz(i)
            nr = this%lattice%nn(i, 1) ! Number of neighbours of atom i
            if (this%izero(i) /= 0) then
               call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%hallo(:, :, 1, i), 18, this%psi1(:, :, i), 18, cone, this%hohpsi(:, :, i), 18)
            end if
            if (nr >= 2) then
               do nb = 2, nr ! Loop on the neighbouring
                  nnmap = this%lattice%nn(i, nb)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%hallo(:, :, nb, i), 18, this%psi1(:, :, nnmap), 18, cone, this%hohpsi(:, :, i), 18)
                        this%idum(i) = 1
                     end if
                  end if
               end do ! End of loop in the neighbouring
            end if
         end do ! End of loop in the neighbouring
      end if ! End of local Hamiltonian loop

      do i = nlimplus1, this%lattice%kk ! Loop to find the bulk atoms using the bulk Hamiltonian
         this%idum(i) = this%izero(i)
         ih = this%lattice%iz(i) ! Atom type
         nr = this%lattice%nn(i, 1) ! Number of neighbours
         if (this%izero(i) /= 0) then
            call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%eeo(:, :, 1, ih), 18, this%psi1(:, :, i), 18, cone, this%hohpsi(:, :, i), 18)
         end if
         if (nr >= 2) then
            do nb = 2, nr ! Loop on the neighbouring
               nnmap = this%lattice%nn(i, nb)
               if (nnmap /= 0) then
                  if (this%izero(nnmap) /= 0) then
                     call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%eeo(:, :, nb, ih), 18, this%psi1(:, :, nnmap), 18, cone, this%hohpsi(:, :, i), 18)
                     this%idum(i) = 1
                  end if
               end if
            end do ! End of loop in the neighbouring
         end if
      end do

      ! H = h - hoh + e_nu + l.s
      this%psi1(:, :, :) = this%psi1(:, :, :) - this%hohpsi(:, :, :) + this%enupsi(:, :, :) + this%socpsi(:, :, :)

      ! Do the scaling and shifting
      this%psi1(:, :, :) = this%psi1(:, :, :) - b*this%psi0(:, :, :)
      this%psi1(:, :, :) = this%psi1(:, :, :)/a

      ! Write the 1st moment
      do n = 1, this%lattice%kk
         if (this%izero(n) /= 0) then
            call zgemm('c', 'n', 18, 18, 18, cone, psiref(:, :, n), 18, this%psi1(:, :, n), 18, cone, this%cheb_mom_temp, 18)
         end if
      end do
   end subroutine cheb_1st_mom_hoh

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Recursion method to calculate Chebyshev moments for the inter-site GFs
   !---------------------------------------------------------------------------
   subroutine chebyshev_recur_ij(this)
      use mpi_mod
      class(recursion), intent(inout) :: this
      ! Local variables
      real(rp) :: a, b
      integer :: i, ij, j, l, ll, kk, m, reci
      integer :: llmax ! Recursion steps
      complex(rp) :: asign, bsign
      complex(rp), dimension(18, 18, this%lattice%kk) :: psiref

      integer :: ij_loc

      a = (this%en%energy_max - this%en%energy_min)/(2 - 0.3)
      b = (this%en%energy_max + this%en%energy_min)/2

      llmax = this%lattice%control%lld
      !do ij=1, this%lattice%njij ! Loop on the number of pair of atoms
      do ij = start_atom, end_atom
         ij_loc = g2l_map(ij)
         i = this%lattice%ijpair(ij, 1) ! Atom number in the clust file, atom i
         j = this%lattice%ijpair(ij, 2) ! Atom number in the clust file, atom j
         call g_logger%info(int2str(rank)//': Chebyshev recursion on progress between atoms '//int2str(i)//' and '//int2str(j), __FILE__, __LINE__)
         do reci = 1, 4
            ! Clear list of atoms in hopping region
            this%izero(:) = 0
            ! Initializing wave functions
            this%psi0(:, :, :) = (0.0d0, 0.0d0)
            this%psi1(:, :, :) = (0.0d0, 0.0d0)
            this%psi2(:, :, :) = (0.0d0, 0.0d0)

            psiref(:, :, :) = (0.0d0, 0.0d0)

            select case (reci)
            case (1)
               asign = (1.0d0, 0.0d0)*one_over_sqrt_two
               bsign = (1.0d0, 0.0d0)*one_over_sqrt_two
               this%izero(i) = 1
               this%izero(j) = 1
               !do l=1,18
               !this%psi0(l,l,j) = asign
               !this%psi0(l,l,j) = bsign
               !psiref(l,l,i) = asign
               !psiref(l,l,i) = bsign
               !end do
            case (2)
               asign = (1.0d0, 0.0d0)*one_over_sqrt_two
               bsign = (-1.0d0, 0.0d0)*one_over_sqrt_two
               this%izero(i) = 1
               this%izero(j) = 1
               !do l=1,18
               !this%psi0(l,l,i) = asign
               !this%psi0(l,l,j) = bsign
               !psiref(l,l,i) = asign
               !psiref(l,l,j) = bsign
               !end do
            case (3)
               asign = (1.0d0, 0.0d0)*one_over_sqrt_two
               bsign = (0.0d0, 1.0d0)*one_over_sqrt_two
               this%izero(i) = 1
               this%izero(j) = 1
            case (4)
               asign = (1.0d0, 0.0d0)*one_over_sqrt_two
               bsign = (0.0d0, -1.0d0)*one_over_sqrt_two
               this%izero(i) = 1
               this%izero(j) = 1
            end select

            do l = 1, 18
               this%psi0(l, l, i) = asign
               this%psi0(l, l, j) = bsign
               psiref(l, l, i) = asign
               psiref(l, l, j) = bsign
            end do

            ! Write the 0th moment
            call g_timer%start('<PSI_0|PSI_0>')
            call this%cheb_0th_mom(psiref)
            this%mu_n(:, :, 1, ij_loc*4 - 4 + reci) = this%cheb_mom_temp(:, :)
            call g_timer%stop('<PSI_0|PSI_0>')

            call g_timer%start('<PSI_0|PSI_1>')
            ! Write the 1st moment
            if (this%hamiltonian%hoh) then
               call this%cheb_1st_mom_hoh(psiref, a, b)
            else
               call this%cheb_1st_mom(psiref, a, b)
            end if
            this%mu_n(:, :, 2, ij_loc*4 - 4 + reci) = this%cheb_mom_temp(:, :)
            call g_timer%stop('<PSI_0|PSI_1>')

            this%izero(:) = this%idum(:)
            ! Start the recursion
            do ll = 1, this%lattice%control%lld ! Loop in the recursion steps
               call g_timer%start('<PSI_0|PSI_n>')
               if (this%hamiltonian%hoh) then
                  call this%chebyshev_recur_ll_hoh(ij_loc*4 - 4 + reci, ll, a, b)
               else
                  call this%chebyshev_recur_ll(ij_loc*4 - 4 + reci, ll, a, b)
               end if
               call g_timer%stop('<PSI_0|PSI_n>')
               if (real(sum(this%mu_n(:, :, ll + 2, ij_loc*4 - 4 + reci))) > 1000.d0) then
                  call g_logger%fatal('Chebyshev moments did not converge. Check energy limits energy_min and energy_max', __FILE__, __LINE__)
               end if
            end do ! End loop in the recursion steps
         end do
      end do
   end subroutine chebyshev_recur_ij

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Recursion method to calculate Chebyshev moments for a given number of
   !  recursion steps.
   !---------------------------------------------------------------------------
   subroutine chebyshev_recur_ll(this, i, ll, a, b)
      class(recursion), intent(inout) :: this
      integer, intent(in) :: i, ll
      real(rp), intent(in) :: a, b
      ! Local variables
      integer :: nb, ih, j, k, nr, m, n, l, hblocksize, nat, nnmap, nlimplus1
      complex(rp), dimension(18, 18) :: dum1, dum2, locham

      hblocksize = 18
      nat = this%lattice%kk

      ! Write H*|phi_1> for the local Hamiltonian
      nlimplus1 = this%lattice%nmax + 1
      if (this%lattice%nmax /= 0) then ! In case of impurities using the local hamiltonian
         !$omp parallel do default(shared) private(k, ih, nr, nb, nnmap,locham)
         do k = 1, this%lattice%nmax ! Loop in the neighbouring
            this%idum(k) = this%izero(k)
            ih = this%lattice%iz(k)
            nr = this%lattice%nn(k, 1) ! Number of neighbours of atom i
            if (this%izero(k) /= 0) then
               locham = this%hamiltonian%hall(1:18, 1:18, 1, k) + this%hamiltonian%lsham(1:18, 1:18, ih)
               call zgemm('n', 'n', 18, 18, 18, cone, locham, 18, this%psi1(:, :, k), 18, cone, this%psi2(:, :, k), 18)
            end if
            if (nr >= 2) then
               do nb = 2, nr ! Loop on the neighbouring
                  nnmap = this%lattice%nn(k, nb)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%hall(:, :, nb, k), 18, this%psi1(:, :, nnmap), 18, cone, this%psi2(:, :, k), 18)
                        this%idum(k) = 1
                     end if
                  end if
               end do ! End of loop in the neighbouring
            end if
            ! Do the scaling and shifting
            this%psi2(:, :, k) = this%psi2(:, :, k) - b*this%psi1(:, :, k)
            this%psi2(:, :, k) = this%psi2(:, :, k)/a
            ! Write 2*H*|phi_1>
            this%psi2(:, :, k) = 2*this%psi2(:, :, k)
         end do ! End of loop in the neighbouring
         !$omp end parallel do
      end if ! End of local Hamiltonian loop

      ! Write H*|phi_1> for the bulk Hamiltonian
      !$omp parallel do default(shared) private(k, ih, nr, nb, nnmap,locham)
      do k = nlimplus1, this%lattice%kk ! Loop in the clust
         this%idum(k) = this%izero(k)
         ih = this%lattice%iz(k)
         nr = this%lattice%nn(k, 1)
         if (this%izero(k) /= 0) then
            locham = this%hamiltonian%ee(1:18, 1:18, 1, ih) + this%hamiltonian%lsham(1:18, 1:18, ih)
            call zgemm('n', 'n', 18, 18, 18, cone, locham, 18, this%psi1(:, :, k), 18, cone, this%psi2(:, :, k), 18)
         end if
         if (nr >= 2) then
            do nb = 2, nr ! Loop in the neighbouring
               nnmap = this%lattice%nn(k, nb)
               if (nnmap /= 0 .and. this%izero(nnmap) /= 0) then
                  call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%ee(1, 1, nb, ih), 18, this%psi1(:, :, nnmap), 18, cone, this%psi2(:, :, k), 18)
                  this%idum(k) = 1
               end if
            end do ! End of the loop in the neighbouring
         end if
         ! Do the scaling and shifting
         this%psi2(:, :, k) = this%psi2(:, :, k) - b*this%psi1(:, :, k)
         this%psi2(:, :, k) = this%psi2(:, :, k)/a
         ! Write 2*H*|phi_1>
         this%psi2(:, :, k) = 2*this%psi2(:, :, k)
      end do ! End loop in the clust
      !$omp end parallel do

      ! Write Moments
      dum1(:, :) = (0.0d0, 0.0d0)
      dum2(:, :) = (0.0d0, 0.0d0)

      ! Redefines idum for all clust atoms
      this%irnum = 0
      do n = 1, this%lattice%kk
         this%izero(n) = this%idum(n)
         if (this%izero(n) /= 0) then
            this%irnum = this%irnum + 1
            this%irlist(this%irnum) = n
         end if
      end do

      !$omp parallel do default(shared) private(n) reduction(+:dum1) reduction(+:dum2)
      do n = 1, this%irnum
         ! Write |phi_2>=2*H*|phi_1> - |phi_0>
         this%psi2(:, :, this%irlist(n)) = this%psi2(:, :, this%irlist(n)) - this%psi0(:, :, this%irlist(n))
         call zgemm('c', 'n', 18, 18, 18, cone, this%psi1(:, :, this%irlist(n)), 18, this%psi1(:, :, this%irlist(n)), 18, cone, dum1, 18)
         call zgemm('c', 'n', 18, 18, 18, cone, this%psi2(:, :, this%irlist(n)), 18, this%psi1(:, :, this%irlist(n)), 18, cone, dum2, 18)
         this%psi0(:, :, this%irlist(n)) = this%psi1(:, :, this%irlist(n))
         this%psi1(:, :, this%irlist(n)) = this%psi2(:, :, this%irlist(n))
         this%psi2(:, :, this%irlist(n)) = (0.0d0, 0.0d0)
      end do
      !$omp end parallel  do

      this%mu_n(:, :, 2*ll + 1, i) = 2.0_rp*dum1(:, :) - this%mu_n(:, :, 1, i)
      this%mu_n(:, :, 2*ll + 2, i) = 2.0_rp*dum2(:, :) - this%mu_n(:, :, 2, i)

      if (sum(real(this%mu_n(:, :, 2*ll + 2, i))) > 1000.d0) then
         call g_logger%fatal('Chebyshev moments did not converge. Check energy limits energy_min and energy_max', __FILE__, __LINE__)
      end if
   end subroutine chebyshev_recur_ll

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Recursion method to calculate Chebyshev moments for a given number of
   !  recursion steps. Includes the hoh term
   !---------------------------------------------------------------------------
   subroutine chebyshev_recur_ll_hoh(this, i, ll, a, b)
      class(recursion), intent(inout) :: this
      integer, intent(in) :: i, ll
      real(rp), intent(in) :: a, b
      ! Local variables
      integer :: nb, ih, j, k, nr, m, n, l, hblocksize, nat, nnmap, nlimplus1
      complex(rp), dimension(18, 18) :: dum1, dum2, locham

      hblocksize = 18
      nat = this%lattice%kk

      this%hohpsi(:, :, :) = (0.0d0, 0.0d0)
      this%enupsi(:, :, :) = (0.0d0, 0.0d0)
      this%socpsi(:, :, :) = (0.0d0, 0.0d0)

      ! Write H*|phi_1> for the local Hamiltonian
      nlimplus1 = this%lattice%nmax + 1
      if (this%lattice%nmax /= 0) then ! In case of impurities using the local hamiltonian
         !$omp parallel do default(shared) private(k, ih, nr, nb, nnmap,locham)
         do k = 1, this%lattice%nmax ! Loop in the neighbouring
            this%idum(k) = this%izero(k)
            ih = this%lattice%iz(k)
            nr = this%lattice%nn(k, 1) ! Number of neighbours of atom i
            if (this%izero(k) /= 0) then
               locham = this%hamiltonian%hall(1:18, 1:18, 1, k)
               call zgemm('n', 'n', 18, 18, 18, cone, locham, 18, this%psi1(:, :, k), 18, cone, this%psi2(:, :, k), 18)
               call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%lsham(:, :, ih), 18, this%psi1(:, :, k), 18, cone, this%socpsi(:, :, k), 18)
               call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%enim(:, :, ih), 18, this%psi1(:, :, k), 18, cone, this%enupsi(:, :, k), 18)
            end if
            if (nr >= 2) then
               do nb = 2, nr ! Loop on the neighbouring
                  nnmap = this%lattice%nn(k, nb)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%hall(:, :, nb, k), 18, this%psi1(:, :, nnmap), 18, cone, this%psi2(:, :, k), 18)
                        this%idum(k) = 1
                     end if
                  end if
               end do ! End of loop in the neighbouring
            end if
         end do ! End of loop in the neighbouring
         !$omp end parallel do
      end if ! End of local Hamiltonian loop

      ! Write H*|phi_1> for the bulk Hamiltonian
      !$omp parallel do default(shared) private(k, ih, nr, nb, nnmap,locham)
      do k = nlimplus1, this%lattice%kk ! Loop in the clust
         this%idum(k) = this%izero(k)
         ih = this%lattice%iz(k)
         nr = this%lattice%nn(k, 1)
         if (this%izero(k) /= 0) then
            locham = this%hamiltonian%ee(1:18, 1:18, 1, ih)
            call zgemm('n', 'n', 18, 18, 18, cone, locham, 18, this%psi1(:, :, k), 18, cone, this%psi2(:, :, k), 18)
            call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%enim(:, :, ih), 18, this%psi1(:, :, k), 18, cone, this%enupsi(:, :, k), 18)
            call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%lsham(:, :, ih), 18, this%psi1(:, :, k), 18, cone, this%socpsi(:, :, k), 18)
         end if
         if (nr >= 2) then
            do nb = 2, nr ! Loop in the neighbouring
               nnmap = this%lattice%nn(k, nb)
               if (nnmap /= 0 .and. this%izero(nnmap) /= 0) then
                  call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%ee(1, 1, nb, ih), 18, this%psi1(:, :, nnmap), 18, cone, this%psi2(:, :, k), 18)
                  this%idum(k) = 1
               end if
            end do ! End of the loop in the neighbouring
         end if
      end do ! End loop in the clust
      !$omp end parallel do

      ! Mapping update for hoh calculation
      this%izero(:) = this%idum(:)

      if (this%lattice%nmax /= 0) then ! In case of impurities using the local hamiltonian
         !$omp parallel do default(shared) private(k, ih, nr, nb, nnmap,locham)
         do k = 1, this%lattice%nmax ! Loop in the neighbouring
            this%idum(k) = this%izero(k)
            ih = this%lattice%iz(k)
            nr = this%lattice%nn(k, 1) ! Number of neighbours of atom i
            if (this%izero(k) /= 0) then
               call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%hallo(:, :, 1, k), 18, this%psi2(:, :, k), 18, cone, this%hohpsi(:, :, k), 18)
            end if
            if (nr >= 2) then
               do nb = 2, nr ! Loop on the neighbouring
                  nnmap = this%lattice%nn(k, nb)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%hallo(:, :, nb, k), 18, this%psi2(:, :, nnmap), 18, cone, this%hohpsi(:, :, k), 18)
                        this%idum(k) = 1
                     end if
                  end if
               end do ! End of loop in the neighbouring
            end if
         end do ! End of loop in the neighbouring
         !$omp end parallel do
      end if ! End of local Hamiltonian loop

      !$omp parallel do default(shared) private(k, ih, nr, nb, nnmap,locham)
      do k = nlimplus1, this%lattice%kk ! Loop to find the bulk atoms using the bulk Hamiltonian
         this%idum(k) = this%izero(k)
         ih = this%lattice%iz(k) ! Atom type
         nr = this%lattice%nn(k, 1) ! Number of neighbours
         if (this%izero(k) /= 0) then
            call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%eeo(:, :, 1, ih), 18, this%psi2(:, :, k), 18, cone, this%hohpsi(:, :, k), 18)
         end if
         if (nr >= 2) then
            do nb = 2, nr ! Loop on the neighbouring
               nnmap = this%lattice%nn(k, nb)
               if (nnmap /= 0) then
                  if (this%izero(nnmap) /= 0) then
                     call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%eeo(:, :, nb, ih), 18, this%psi2(:, :, nnmap), 18, cone, this%hohpsi(:, :, k), 18)
                     this%idum(k) = 1
                  end if
               end if
            end do ! End of loop in the neighbouring
         end if
      end do
      !$omp end parallel do

      ! H = h - hoh + e_nu + l.s
      this%psi2(:, :, :) = this%psi2(:, :, :) - this%hohpsi(:, :, :) + this%enupsi(:, :, :) + this%socpsi(:, :, :)

      ! Do the scaling and shifting
      this%psi2(:, :, :) = this%psi2(:, :, :) - b*this%psi1(:, :, :)
      this%psi2(:, :, :) = this%psi2(:, :, :)/a
      ! Write 2*H*|phi_1>
      this%psi2(:, :, :) = 2*this%psi2(:, :, :)

      ! Write Moments
      dum1(:, :) = (0.0d0, 0.0d0)
      dum2(:, :) = (0.0d0, 0.0d0)

      ! Redefines idum for all clust atoms
      this%irnum = 0
      do n = 1, this%lattice%kk
         this%izero(n) = this%idum(n)
         if (this%izero(n) /= 0) then
            this%irnum = this%irnum + 1
            this%irlist(this%irnum) = n
         end if
      end do

      !$omp parallel do default(shared) private(n) reduction(+:dum1) reduction(+:dum2)
      do n = 1, this%irnum
         ! Write |phi_2>=2*H*|phi_1> - |phi_0>
         this%psi2(:, :, this%irlist(n)) = this%psi2(:, :, this%irlist(n)) - this%psi0(:, :, this%irlist(n))
         call zgemm('c', 'n', 18, 18, 18, cone, this%psi1(:, :, this%irlist(n)), 18, this%psi1(:, :, this%irlist(n)), 18, cone, dum1, 18)
         call zgemm('c', 'n', 18, 18, 18, cone, this%psi2(:, :, this%irlist(n)), 18, this%psi1(:, :, this%irlist(n)), 18, cone, dum2, 18)
         this%psi0(:, :, this%irlist(n)) = this%psi1(:, :, this%irlist(n))
         this%psi1(:, :, this%irlist(n)) = this%psi2(:, :, this%irlist(n))
         this%psi2(:, :, this%irlist(n)) = (0.0d0, 0.0d0)
      end do
      !$omp end parallel  do

      this%mu_n(:, :, 2*ll + 1, i) = 2.0_rp*dum1(:, :) - this%mu_n(:, :, 1, i)
      this%mu_n(:, :, 2*ll + 2, i) = 2.0_rp*dum2(:, :) - this%mu_n(:, :, 2, i)

      if (sum(real(this%mu_n(:, :, 2*ll + 2, i))) > 1000.d0) then
         call g_logger%fatal('Chebyshev moments did not converge. Check energy limits energy_min and energy_max', __FILE__, __LINE__)
      end if
   end subroutine chebyshev_recur_ll_hoh

  !!!!!!!!!!!!!!!!!!!!!!!
   !OBSOLETE TO NOT USE!!!
  !!!!!!!!!!!!!!!!!!!!!!!
   subroutine chebyshev_recur_s_ll(this, psiref, ll, a, b)
      class(recursion), intent(inout) :: this
      integer, intent(in) :: ll
      real(rp), intent(in) :: a, b
      complex(rp), dimension(18, 18, this%lattice%kk), intent(in) :: psiref
      ! Local variables
      integer :: nb, ih, j, k, nr, m, n, l, hblocksize, nat, nnmap, nlimplus1
      complex(rp), dimension(18, 18) :: dum1, dum2
      real(rp) :: inv_a

      hblocksize = 18
      nat = this%lattice%kk
      nlimplus1 = this%lattice%nmax + 1
      inv_a = 1.0_rp/a

      this%cheb_mom_temp(:, :) = 0.0d0
      ! Write H*|phi_1>
      !$omp parallel do default(shared) private(k, ih, nr, nb, nnmap)
      do k = nlimplus1, this%lattice%kk ! Loop in the clust
         this%idum(k) = this%izero(k)
         ih = this%lattice%iz(k)
         nr = this%lattice%nn(k, 1)
         if (this%izero(k) /= 0) then
            call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%ee(:, :, 1, ih), 18, this%psi1(:, :, k), 18, cone, this%psi2(:, :, k), 18)
         end if
         if (nr >= 2) then
            do nb = 2, nr ! Loop in the neighbouring
               nnmap = this%lattice%nn(k, nb)
               if (nnmap /= 0 .and. this%izero(nnmap) /= 0) then
                  call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%ee(:, :, nb, ih), 18, this%psi1(:, :, nnmap), 18, cone, this%psi2(:, :, k), 18)
                  this%idum(k) = 1
               end if
            end do ! End of the loop in the neighbouring
         end if
         ! Do the scaling and shifting
         this%psi2(:, :, k) = this%psi2(:, :, k) - b*this%psi1(:, :, k)
         this%psi2(:, :, k) = inv_a*this%psi2(:, :, k)
         ! Write 2*H*|phi_1>
         this%psi2(:, :, k) = 2.0_rp*this%psi2(:, :, k)
      end do ! End loop in the clust
      !$omp end parallel do

      ! Write Moments
      dum1(:, :) = (0.0d0, 0.0d0)
      !$omp parallel do default(shared) private(n) reduction(+:dum1)
      do n = 1, this%lattice%kk
         ! Write |phi_2>=2*H*|phi_1> - |phi_0>
         this%psi2(:, :, n) = this%psi2(:, :, n) - this%psi0(:, :, n)
         if (this%izero(n) /= 0) then
            call zgemm('c', 'n', 18, 18, 18, cone, psiref(:, :, n), 18, this%psi2(:, :, n), 18, cone, dum1, 18)
         end if
      end do
      !$omp end parallel do
      this%cheb_mom_temp(:, :) = dum1(:, :)

      this%psi0(:, :, :) = this%psi1(:, :, :)
      this%psi1(:, :, :) = this%psi2(:, :, :)
      this%psi2(:, :, :) = (0.0d0, 0.0d0)
   end subroutine chebyshev_recur_s_ll

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Recursion method using Chebyshev moments
   !---------------------------------------------------------------------------
   subroutine chebyshev_recur(this)
      use mpi_mod
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: nb, ih, i, j, k, nr, ll, m, n, l, hblocksize, nat, nnmap, nlimplus1, llcheb
      integer, dimension(0:this%lattice%kk) :: idumll
      !complex(rp) :: cone = (1.0d0, 0.0d0)
      complex(rp), dimension(18, 18) :: dum, dum1, dum2
      complex(rp), dimension(18, this%lattice%kk) :: v
      complex(rp), dimension(:, :, :, :), allocatable :: hcheb
      real(rp) :: a, b, start, finish
      ! External functions
      complex(rp), external :: zdotc

      integer :: i_loc

      hblocksize = 18
      nat = this%lattice%kk
      nlimplus1 = this%lattice%nmax + 1
      llcheb = (2*this%control%lld) + 2

      a = (this%en%energy_max - this%en%energy_min)/(2 - 0.3)
      b = (this%en%energy_max + this%en%energy_min)/2

      do i = start_atom, end_atom ! Loop on the number of atoms to be treat self-consistently by this process
         i_loc = g2l_map(i)
         j = this%lattice%irec(i) ! Atom number in the clust file
         call g_logger%info('Chebyshev recursion on progress for atom '//int2str(j), __FILE__, __LINE__)
         ! Initialize neighbouring map
         this%izeroll(:, :) = 0
         this%izeroll(j, 1) = 1

         this%izero(:) = 0
         this%izero(j) = 1
         call this%create_ll_map()

         ! Initializing wave functions
         this%psi0(:, :, :) = (0.0d0, 0.0d0)
         this%psi1(:, :, :) = (0.0d0, 0.0d0)
         this%psi2(:, :, :) = (0.0d0, 0.0d0)

         do l = 1, 18
            !> Starting state for |phi_0>
            this%psi0(l, l, j) = (1.0d0, 0.0d0)
         end do

         ! Write the 0th moment
         call g_timer%start('<PSI_0|PSI_0>')
         call this%cheb_0th_mom(this%psi0)
         this%mu_n(:, :, 1, i_loc) = (this%cheb_mom_temp(:, :))
         call g_timer%stop('<PSI_0|PSI_0>')

         call g_timer%start('<PSI_0|PSI_1>')
         if (this%hamiltonian%hoh) then
            call this%cheb_1st_mom_hoh(this%psi0, a, b)
         else
            call this%cheb_1st_mom(this%psi0, a, b)
         end if
         this%mu_n(:, :, 2, i_loc) = (this%cheb_mom_temp(:, :))
         call g_timer%stop('<PSI_0|PSI_1>')

         this%izero(:) = this%idum(:)
         ! Start the recursion
         do ll = 1, this%lattice%control%lld
            call g_timer%start('<PSI_0|PSI_n>')
            if (this%hamiltonian%hoh) then
               call this%chebyshev_recur_ll_hoh(i_loc, ll, a, b)
            else
               call this%chebyshev_recur_ll(i_loc, ll, a, b)
            end if
            call g_timer%stop('<PSI_0|PSI_n>')
         end do ! End loop in the recursion steps
      end do ! End loop on the number of atoms to be treat self-consistently
   end subroutine chebyshev_recur

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Recursion method using to find the Chebyshev moments for diagonal and
   !off-diagonal terms. Here, one cannot use the ´double-trick´. See
   ! Rev. Mod. Phys. 78, 275.
   ! OBSOLETE! DO NOT USE
   !---------------------------------------------------------------------------
   subroutine chebyshev_recur_full(this)
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: nb, ih, i, j, k, nr, ll, m, n, l, hblocksize, nat, nnmap, nlimplus1
      integer, dimension(0:this%lattice%kk) :: idumll
      complex(rp) :: cone = (1.0d0, 0.0d0)
      complex(rp), dimension(18, 18) :: dum, dum1, dum2
      complex(rp), dimension(18, this%lattice%kk) :: v
      complex(rp), dimension(:, :, :, :), allocatable :: hcheb
      complex(rp), dimension(:, :, :), allocatable :: psiref
      real(rp) :: a, b, start, finish
      ! External functions
      complex(rp), external :: zdotc

      hblocksize = 18
      nat = this%lattice%kk
      nlimplus1 = this%lattice%nmax + 1
      allocate (hcheb(18, 18, (this%lattice%nn(1, 1) + 1), this%lattice%kk), psiref(18, 18, this%lattice%kk))

      a = (this%en%energy_max - this%en%energy_min)/(2 - 0.3)
      b = (this%en%energy_max + this%en%energy_min)/2

      hcheb(:, :, :, :) = this%hamiltonian%ee(:, :, :, :) !hcheb(:, :, :, :)!/cmplx(a, 0.d0)

      do i = 1, this%lattice%nrec ! Loop on the number of atoms to be treat self-consistently
         j = this%lattice%irec(i) ! Atom number in the clust file
         ! Initialize neighbouring map
         this%izeroll(:, :) = 0
         this%izeroll(j, 1) = 1

         call this%create_ll_map()
         print '("Mapping neighrbours time = ", f10.3, " seconds.")', (finish - start)!/32

         ! Initializing wave functions
         this%psi0(:, :, :) = (0.0d0, 0.0d0)
         this%psi1(:, :, :) = (0.0d0, 0.0d0)
         this%psi2(:, :, :) = (0.0d0, 0.0d0)

         do l = 1, 18
            !> Starting state for |phi_0>
            this%psi0(l, l, j) = (1.0d0, 0.0d0)
            psiref(l, l, j) = (1.0d0, 0.0d0)
         end do

         dum(:, :) = (0.0d0, 0.0d0)
         ! Write the 0th moment
         call zgemm('c', 'n', 18, 18, 18, cone, this%psi0(:, :, j), 18, this%psi0(:, :, j), 18, cone, dum, 18)
         this%mu_n(:, :, 1, i) = (dum(:, :))

         ! Write |phi_1>=H|phi_0>
         !$omp parallel do default(shared) private(k, ih, nr, nb, nnmap)
         do k = nlimplus1, this%lattice%kk ! Loop in the clust
            idumll(k) = this%izeroll(k, 1)
            ih = this%lattice%iz(k)
            nr = this%lattice%nn(k, 1)
            if (this%izeroll(k, 1) /= 0) then
               call zgemm('n', 'n', 18, 18, 18, cone, hcheb(:, :, 1, ih), 18, this%psi0(:, :, k), 18, cone, this%psi1(:, :, k), 18)
            end if
            if (nr >= 2) then
               do nb = 2, nr ! Loop in the neighbouring
                  nnmap = this%lattice%nn(k, nb)
                  if (nnmap /= 0 .and. this%izeroll(nnmap, 1) /= 0) then
                     call zgemm('n', 'n', 18, 18, 18, cone, hcheb(:, :, nb, ih), 18, this%psi0(:, :, nnmap), 18, cone, this%psi1(:, :, k), 18)
                  end if
               end do ! End of the loop in the neighbouring
            end if
            ! Do the scaling and shifting
            this%psi1(:, :, k) = this%psi1(:, :, k) - b*this%psi0(:, :, k)
            this%psi1(:, :, k) = this%psi1(:, :, k)/a
         end do ! End loop in the clust
         !$omp end parallel do

         ! Write the 1st moment
         dum(:, :) = (0.0d0, 0.0d0)
         do n = 1, this%lattice%kk
            if (this%izeroll(n, 1) /= 0) then
               call zgemm('c', 'n', 18, 18, 18, cone, psiref(:, :, n), 18, this%psi1(:, :, n), 18, cone, dum, 18)
            end if
         end do
         this%mu_n(:, :, 2, i) = (dum(:, :))

         ! Start the recursion
         do ll = 1, this%lattice%control%lld
            ! Write H*|phi_1>
            !$omp parallel do default(shared) private(k, ih, nr, nb, nnmap)
            do k = nlimplus1, this%lattice%kk ! Loop in the clust
               idumll(k) = this%izeroll(k, ll + 1)
               ih = this%lattice%iz(k)
               nr = this%lattice%nn(k, 1)
               if (this%izeroll(k, ll + 1) /= 0) then
                  call zgemm('n', 'n', 18, 18, 18, cone, hcheb(:, :, 1, ih), 18, this%psi1(:, :, k), 18, cone, this%psi2(:, :, k), 18)
               end if
               if (nr >= 2) then
                  do nb = 2, nr ! Loop in the neighbouring
                     nnmap = this%lattice%nn(k, nb)
                     if (nnmap /= 0 .and. this%izeroll(nnmap, ll + 1) /= 0) then
                        call zgemm('n', 'n', 18, 18, 18, cone, hcheb(:, :, nb, ih), 18, this%psi1(:, :, nnmap), 18, cone, this%psi2(:, :, k), 18)
                     end if
                  end do ! End of the loop in the neighbouring
               end if
               ! Do the scaling and shifting
               this%psi2(:, :, k) = this%psi2(:, :, k) - b*this%psi1(:, :, k)
               this%psi2(:, :, k) = this%psi2(:, :, k)/a
               ! Write 2*H*|phi_1>
               this%psi2(:, :, k) = 2*this%psi2(:, :, k)
            end do ! End loop in the clust
            !$omp end parallel do

            ! Write |phi_2>=2*H*|phi_1> - |phi_0>
            this%psi2(:, :, :) = this%psi2(:, :, :) - this%psi0(:, :, :)

            ! Write Moments
            dum1(:, :) = (0.0d0, 0.0d0)
            !$omp parallel do default(shared) private(n) reduction(+:dum1)
            do n = 1, this%lattice%kk
               if (this%izeroll(n, ll + 1) /= 0) then
                  call zgemm('c', 'n', 18, 18, 18, cone, psiref(:, :, n), 18, this%psi2(:, :, n), 18, cone, dum1, 18)
               end if
            end do
            !$omp end parallel do
            this%mu_n(:, :, ll + 2, i) = (dum1(:, :))

            this%psi0(:, :, :) = this%psi1(:, :, :)
            this%psi1(:, :, :) = this%psi2(:, :, :)
            this%psi2(:, :, :) = (0.0d0, 0.0d0)
         end do ! End loop in the recursion steps
      end do ! End loop on the number of atoms to be treat self-consistently

      deallocate (hcheb)
   end subroutine chebyshev_recur_full

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the map of neighbouring map as a function of the recursion
   !steps
   !---------------------------------------------------------------------------
   subroutine create_ll_map(this)
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, nr, nnmap, ll
      integer, dimension(0:this%lattice%kk) :: idumll

      idumll(:) = 0

      do ll = 1, this%lattice%control%lld
         idumll(:) = 0
         do i = 1, this%lattice%kk ! Loop to find the bulk atoms using the bulk Hamiltonian
            idumll(i) = this%izeroll(i, ll)
            nr = this%lattice%nn(i, 1) ! Number of neighbours
            if (nr >= 2) then
               do j = 2, nr ! Loop on the neighbouring
                  nnmap = this%lattice%nn(i, j)
                  if (nnmap /= 0) then
                     if (this%izeroll(nnmap, ll) /= 0) then
                        idumll(i) = 1
                     end if
                  end if
               end do ! End of loop in the neighbouring
            end if
         end do
         this%izeroll(:, ll + 1) = idumll(:)
      end do
   end subroutine create_ll_map

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the recursion coefficients A.
   !---------------------------------------------------------------------------
   subroutine hop(this, ll)
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m, n, nr, nnmap, nlimplus1
      integer :: ll ! Recursion step
      integer :: ino ! Atom type
      integer :: ih ! Atom number in the clust
      integer, dimension(0:this%lattice%kk) :: idum
      complex(rp), dimension(18) :: dum
      real(rp) :: summ, start, finish

      summ = 0.0d0
      idum(:) = 0

      nlimplus1 = this%lattice%nmax + 1

      select case (this%lattice%control%nsp)

      case (1) ! Scalar relativistic case

         if (this%lattice%nmax /= 0) then ! In case of impurities using the local hamiltonian
            !$omp parallel do default(shared) private(i, nr, m, l, dum, j, nnmap) schedule(dynamic, 100)
            do i = 1, this%lattice%nmax ! Loop in the neighbouring
               idum(i) = this%izero(i)
               dum(:) = (0.0d0, 0.0d0)
               nr = this%lattice%nn(i, 1) ! Number of neighbours of atom i
               if (this%izero(i) /= 0) then
                  do m = 1, 9 ! Loop on the orbital m
                     do l = 1, 9 ! Loop on the orbital l
                        dum(l) = dum(l) + this%hamiltonian%hall(l, m, 1, i)*this%psi(m, i)
                        dum(l + 9) = dum(l + 9) + this%hamiltonian%hall(l + 9, m + 9, 1, i)*this%psi(m + 9, i)
                     end do ! End of loop on orbital m
                  end do ! End of loop on orbital l
               end if
               if (nr >= 2) then
                  do j = 2, nr ! Loop on the neighbouring
                     nnmap = this%lattice%nn(i, j)
                     if (nnmap /= 0) then
                        if (this%izero(nnmap) /= 0) then
                           do m = 1, 9 ! Loop in the orbital m
                              do l = 1, 9 ! Loop in the orbital l
                                 dum(l) = dum(l) + this%hamiltonian%hall(l, m, j, i)*this%psi(m, nnmap)
                                 dum(l + 9) = dum(l + 9) + this%hamiltonian%hall(l + 9, m + 9, j, i)*this%psi(m + 9, nnmap)
                              end do ! End of loop in orbital m
                           end do ! End of loop in orbital l
                           idum(i) = 1
                        end if
                     end if
                  end do ! End of loop in the neighbouring
               end if
               do l = 1, 18
                  this%v(l, i) = dum(l)
               end do
            end do ! End of loop in the neighbouring
            !$omp end parallel do
         end if ! End of local Hamiltonian loop

         !$omp parallel do default(shared) private(i, ih, nr, m, l, dum, j, nnmap) schedule(dynamic, 100)
         do i = nlimplus1, this%lattice%kk ! Loop to find the bulk atoms using the bulk Hamiltonian
            idum(i) = this%izero(i)
            ih = this%lattice%iz(i) ! Atom type
            dum(:) = (0.0d0, 0.0d0)
            nr = this%lattice%nn(i, 1) ! Number of neighbours
            if (this%izero(i) /= 0) then
               !write(125, *) ´i is ´, i
               do m = 1, 9 ! Loop on the orbital m
                  do l = 1, 9 ! Loop on the orbital l
                     dum(l) = dum(l) + this%hamiltonian%ee(l, m, 1, ih)*this%psi(m, i)
                     dum(l + 9) = dum(l + 9) + this%hamiltonian%ee(l + 9, m + 9, 1, ih)*this%psi(m + 9, i)
                  end do ! End of the loop on the orbital l
               end do ! End of loop on the orbital m
            end if
            if (nr >= 2) then
               do j = 2, nr ! Loop on the neighbouring
                  nnmap = this%lattice%nn(i, j)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        !write(125, *) i, j, this%lattice%nn(i, j)
                        do m = 1, 9 ! Loop on the orbital m
                           do l = 1, 9 ! Loop on orbital l
                              dum(l) = dum(l) + this%hamiltonian%ee(l, m, j, ih)*this%psi(m, nnmap)
                              dum(l + 9) = dum(l + 9) + this%hamiltonian%ee(l + 9, m + 9, j, ih)*this%psi(m + 9, nnmap)
                           end do ! End of loop on orbital l
                        end do ! End of loop on the orbital m
                        idum(i) = 1
                     end if
                  end if
               end do ! End of loop in the neighbouring
            end if
            do l = 1, 18
               this%v(l, i) = dum(l)
            end do
         end do
         !$omp end parallel do

         ! Redefines idum for all clust atoms
         do i = 1, this%lattice%kk
            this%izero(i) = idum(i)
            do l = 1, 18
               dum(l) = this%v(l, i)
               summ = summ + real(dum(l)*conjg(this%psi(l, i)))
               this%pmn(l, i) = dum(l) + this%pmn(l, i)
            end do
         end do
         this%atemp(ll) = summ
      end select
   end subroutine hop

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the recursion coefficients B2.
   !---------------------------------------------------------------------------
   subroutine crecal(this)
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m
      integer :: ll, llmax ! Recursion step
      integer :: nm1 ! LL-1
      integer :: nat ! Clust size
      integer :: hblocksize ! Hamiltonian size (18)
      real(rp) :: s, summ, start, finish
      complex(rp) :: dum, ajc
      complex(rp), dimension(18, this%lattice%kk) :: thpsi
      character(len=1) :: transa
      character, dimension(6) :: matdescra
      ! External functions
      complex(rp), external :: zdotc

      summ = this%b2temp(1)
      thpsi(:, :) = (0.0d0, 0.0d0)

      hblocksize = 18
      nat = this%lattice%kk
      nm1 = this%lattice%control%lld - 1
      llmax = this%lattice%control%lld
      do ll = 1, nm1
!     call g_timer%start(´hop´)
         call this%hop(ll)
!     call g_timer%stop(´hop´)
         this%b2temp(ll) = summ
         ajc = -this%atemp(ll)

         call zaxpy(nat*hblocksize, ajc, this%psi, 1, this%pmn, 1)

         summ = 0.0d0

         ! In case zdotc function lapack is not working
         do i = 1, nat
            do k = 1, 18
               summ = summ + real(conjg(this%pmn(k, i))*this%pmn(k, i))
            end do
         end do

         ! summ = real(zdotc(nat*hblocksize, this%pmn, 1, this%pmn, 1))

         s = 1.0d0/sqrt(summ)

         thpsi(:, :) = this%pmn(:, :)*s
         this%pmn(:, :) = this%psi(:, :)
         this%psi(:, :) = thpsi(:, :)

         s = sqrt(summ)

         this%pmn(:, :) = -this%pmn(:, :)*s
      end do

      this%b2temp(llmax) = summ
   end subroutine crecal

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the recursion coefficients A and B2.
   !---------------------------------------------------------------------------
   subroutine recur(this)
      use mpi_mod
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, l, ll
      integer :: llmax ! Recursion steps
      integer :: i_loc

      llmax = this%lattice%control%lld

      !do i=1, this%lattice%nrec ! Loop on the number of atoms to be treat self-consistently
      do i = start_atom, end_atom
         i_loc = g2l_map(i)
         j = this%lattice%irec(i) ! Atom number in the clust file
         do l = 1, 18 ! Loop on the orbitals
            ! Clear list of atoms in hopping region

            this%izero(:) = 0
            ! Initializing wave functions
            this%psi(:, :) = (0.0d0, 0.0d0)
            this%pmn(:, :) = (0.0d0, 0.0d0)

            this%psi(l, j) = (1.0d0, 0.0d0)
            this%izero(j) = 1
            this%atemp(:) = 0.0d0; this%b2temp(:) = 0.0d0

            this%b2temp(1) = 1.0d0
            this%atemp(llmax) = 0.0d0
            !write(125, *) ´orbital ´, l
            call this%crecal()

            do ll = 1, llmax ! Loop on the recursion steps
               this%a(ll, l, i_loc, 1) = this%atemp(ll)
               this%b2(ll, l, i_loc, 1) = this%b2temp(ll)
            end do ! End of the loop on the recursion steps
         end do ! End of the loop on the orbitals
      end do ! End of the loop on the nrec

    !!!! For debug purposes
    !!!do i=1, this%lattice%nrec ! Loop on the number of atoms to be treat self-consistently
    !!!  do l=1, 18  ! Loop on the orbital l
    !!!    write(123, *) ´orbital´, l
    !!!    do ll=1, llmax ! Loop on the recursion steps
    !!!      write(123, *) this%a(ll, l, i_loc, 1), this%b2(ll, l, i_loc, 1)
    !!!    end do
    !!!  end do
    !!!end do
   end subroutine

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Obtain the optmal values of the terminators ainf and binf by Pettifor´s
   !termination
   !---------------------------------------------------------------------------
   subroutine bpopt(this, ll, a, rb, n, ainf, rbinf, ifail)
      class(recursion), intent(inout) :: this
      ! Input
      integer, intent(in) :: n, ll
      real(rp), dimension(ll), intent(in) :: A, RB
      ! Output
      integer, intent(out) :: ifail
      real(rp), intent(out) :: ainf, rbinf
      ! Local variables
      integer :: I, JITER, NDIME
      real(rp) :: BM, BMAX, BMIN, EPS
      real(rp), dimension(ll) :: AZ, RBZ

      NDIME = ll
      IFAIL = 0
      EPS = 1.0d-05
      JITER = 0
      AINF = A(N)
      do
         JITER = JITER + 1
         AZ(1) = 0.5d0*(A(1) - AINF)
         do I = 2, N - 1
            AZ(I) = 0.5d0*(A(I) - AINF)
            RBZ(I) = 0.5d0*RB(I)
         end do
         AZ(N) = A(N) - AINF
         RBZ(N) = 1.0d0/SQRT(2.0d0)*RB(N)
         call this%EMAMI(NDIME, AZ, RBZ, N, BMAX, BMIN)
         BM = BMAX + BMIN
         BM = ABS(BM)
         AINF = AINF + (BMAX + BMIN)
         if (BM <= EPS) then
            exit
            !    elseif (JITER > 30) then
         elseif (JITER > 300) then
            IFAIL = 1
            exit
         end if
      end do
      RBINF = (BMAX - BMIN)/2.0d0
      ! write(700, *) AINF, RBINF
   end subroutine bpopt

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Obtain the max. and the min. eigenvalues of symmetric tridiagonal matrix
   !by bisection method
   !---------------------------------------------------------------------------
   subroutine emami(this, nl, as, bs, n, emax, emin)
      class(recursion), intent(inout) :: this
      ! Input
      integer, intent(in) :: N, NL
      real(rp), dimension(NL), intent(in) :: AS, BS
      ! Output
      real(rp), intent(out) :: EMAX, EMIN
      ! Local Variables
      integer :: I, ISTOP, NUM
      real(rp) :: DELE, E, E1, E2, EMAX0, EMIN0, EPS, P, RELFEH, X1, X2
      real(rp), dimension(NL) :: A, B

      EMAX0 = -1.0d6
      EMIN0 = +1.0d6
      do I = 1, N
         A(I) = AS(I)
         B(I) = BS(I)
      end do
      B(1) = 0.0d0
      B(N + 1) = 0.0d0
      do I = 1, N
         X1 = A(I) + ABS(B(I)) + ABS(B(I + 1))
         X2 = A(I) - ABS(B(I)) - ABS(B(I + 1))
         if (EMAX0 <= X1) then
            EMAX0 = X1
         end if
         if (EMIN0 > X2) then
            EMIN0 = X2
         end if
      end do
      RELFEH = 2.d0**(-39)
      EPS = 1.0d-6
      !C....CALCULATION OF EMAX
      ISTOP = 0
      EMAX = EMAX0
      EMIN = EMIN0
      do
         E = (EMAX + EMIN)/2.0d0
         ISTOP = ISTOP + 1
         if (ISTOP > 50) goto 1000
         NUM = 0
         P = A(1) - E
         if (P < 0.d0) then
            NUM = NUM + 1
         end if
         do I = 2, N
            if (P == 0.0) then
               P = (A(I) - E) - ABS(B(I))/RELFEH
               if (P < 0.d0) then
                  NUM = NUM + 1
               end if
            else
               P = (A(I) - E) - B(I)**2/P
               if (P < 0.d0) then
                  NUM = NUM + 1
               end if
            end if
         end do
         if (NUM == N) then
            EMAX = E
         end if
         if (NUM < N) then
            EMIN = E
         end if
         DELE = (EMAX - EMIN)/((EMAX + EMIN)/2.d0)
         DELE = ABS(DELE)
         if (DELE <= EPS) exit
      end do
      E1 = E
      !.....CALCULATION ON EMIN
      ISTOP = 0
      EMAX = E1
      EMIN = EMIN0
      do
         E = (EMAX + EMIN)/2.d0
         ISTOP = ISTOP + 1
         if (ISTOP > 50) goto 1000
         NUM = 0
         P = A(1) - E
         if (P < 0.d0) then
            NUM = NUM + 1
         end if
         do I = 2, N
            if (P == 0.d0) then
               P = (A(I) - E) - ABS(B(I))/RELFEH
               if (P < 0.d0) then
                  NUM = NUM + 1
               end if
            else
               P = (A(I) - E) - B(I)**2/P
               if (P < 0.d0) then
                  NUM = NUM + 1
               end if
            end if
         end do
         if (NUM == 0) then
            EMIN = E
         end if
         if (NUM > 0) then
            EMAX = E
         end if
         DELE = (EMAX - EMIN)/((EMAX + EMIN)/2.d0)
         DELE = ABS(DELE)
         if (DELE <= EPS) exit
      end do
      E2 = E
      !....
      EMAX = E1
      EMIN = E2
      return
1000  continue
      !1000 write (6, 10000)
      return
      !
      ! ... Format Declarations ...
      !
10000 format(" ", "NON-CONVERGE IN EMAMI")
   end subroutine emami

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Reset all members to default
   !---------------------------------------------------------------------------
   subroutine restore_to_default(this, full)
      class(recursion) :: this
      logical, intent(in), optional :: full
      integer :: lmax

      ! TODO: change this hard check for lmax (assuming all atoms have the same lmax)
      !lmax = this%symbolic_atom(this%lattice%iz(1))%potential%lmax
      lmax = 2

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('recursion.a', this%a, (/max(this%lattice%control%llsp, this%lattice%control%lld), 18, this%lattice%nrec, 3/))
      call g_safe_alloc%allocate('recursion.atemp', this%atemp, max(this%lattice%control%llsp, this%lattice%control%lld))
      call g_safe_alloc%allocate('recursion.b2temp', this%b2temp, max(this%lattice%control%llsp, this%lattice%control%lld))
      call g_safe_alloc%allocate('recursion.b2', this%b2, (/max(this%lattice%control%llsp, this%lattice%control%lld), 18, this%lattice%nrec, 3/))
      allocate (this%izero(0:this%lattice%kk), this%izeroll(0:this%lattice%kk, this%lattice%control%lld + 1), this%idum(0:this%lattice%kk))
      call g_safe_alloc%report_allocate('recursion.izero', this%izero)
      call g_safe_alloc%report_allocate('recursion.izeroll', this%izeroll)
      call g_safe_alloc%report_allocate('recursion.idum', this%idum)
      allocate (this%irlist(0:this%lattice%kk))
      call g_safe_alloc%report_allocate('recursion.irlist', this%irlist)

      call g_safe_alloc%allocate('recursion.psi', this%psi, (/18, this%lattice%kk/))
      call g_safe_alloc%allocate('recursion.pmn', this%pmn, (/18, this%lattice%kk/))
      call g_safe_alloc%allocate('recursion.psi1', this%psi1, (/18, 18, this%lattice%kk/))
      call g_safe_alloc%allocate('recursion.psi2', this%psi2, (/18, 18, this%lattice%kk/))
      call g_safe_alloc%allocate('recursion.psi0', this%psi0, (/18, 18, this%lattice%kk/))
      call g_safe_alloc%allocate('recursion.v', this%v, (/18, this%lattice%kk/))
      if ((this%lattice%njij == 0) .and. (this%lattice%njijk == 0)) then
         call g_safe_alloc%allocate('recursion.a_b', this%a_b, (/18, 18, this%control%lld, this%lattice%nrec/))
         call g_safe_alloc%allocate('recursion.b2_b', this%b2_b, (/18, 18, this%control%lld,, this%lattice%nrec/))
         call g_safe_alloc%allocate('recursion.mu_n', this%mu_n, (/18, 18, (2*this%lattice%control%lld) + 2, this%lattice%nrec/))
         call g_safe_alloc%allocate('recursion.mu_ng', this%mu_ng, (/18, 18, (2*this%lattice%control%lld) + 2, this%lattice%nrec/))
      else
         call g_safe_alloc%allocate('recursion.a_b', this%a_b, (/2*(lmax + 1)**2, 2*(lmax + 1)**2, this%control%lld, this%lattice%njij*4/))
         call g_safe_alloc%allocate('recursion.b2_b', this%b2_b, (/2*(lmax + 1)**2, 2*(lmax + 1)**2, this%control%lld,, this%lattice%njij*4/))
         call g_safe_alloc%allocate('recursion.mu_n', this%mu_n, (/2*(lmax + 1)**2, 2*(lmax + 1)**2, (2*this%lattice%control%lld) + 2, this%lattice%njij*4/))
         call g_safe_alloc%allocate('recursion.mu_ng', this%mu_ng, (/2*(lmax + 1)**2, 2*(lmax + 1)**2, (2*this%lattice%control%lld) + 2, this%lattice%njij*4/))
      end if
      call g_safe_alloc%allocate('recursion.psi_b', this%psi_b, (/18, 18, this%lattice%kk/))
      call g_safe_alloc%allocate('recursion.hpsi', this%hpsi, (/18, 18, this%lattice%kk/))
      call g_safe_alloc%allocate('recursion.hohpsi', this%hohpsi, (/18, 18, this%lattice%kk/))
      call g_safe_alloc%allocate('recursion.enupsi', this%enupsi, (/18, 18, this%lattice%kk/))
      call g_safe_alloc%allocate('recursion.socpsi', this%socpsi, (/18, 18, this%lattice%kk/))
      call g_safe_alloc%allocate('recursion.atemp_b', this%atemp_b, (/18, 18, this%control%lld/))
      call g_safe_alloc%allocate('recursion.b2temp_b', this%b2temp_b, (/18, 18, this%control%lld/))
      call g_safe_alloc%allocate('recursion.pmn_b', this%pmn_b, (/18, 18, this%lattice%kk/))

#else
      allocate (this%a(max(this%lattice%control%llsp, this%lattice%control%lld),&
                          &18, this%lattice%nrec, 3))
      allocate (this%atemp(max(this%lattice%control%llsp, this%lattice%control%lld)))
      allocate (this%b2temp(max(this%lattice%control%llsp, this%lattice%control%lld)))
      allocate (this%b2(max(this%lattice%control%llsp, this%lattice%control%lld),&
                          &18, this%lattice%nrec, 3))
      allocate (this%izero(0:this%lattice%kk), this%idum(0:this%lattice%kk), this%izeroll(0:this%lattice%kk, this%lattice%control%lld + 1))
      allocate (this%irlist(0:this%lattice%kk))

      allocate (this%psi(18, this%lattice%kk), this%pmn(18, this%lattice%kk))
      allocate (this%psi1(18, 18, this%lattice%kk), this%psi2(18, 18, this%lattice%kk), this%psi0(18, 18, this%lattice%kk))
      allocate (this%v(18, this%lattice%kk))
      if (this%lattice%njij == 0) then
         allocate (this%a_b(18, 18, this%control%lld, this%lattice%nrec))
         allocate (this%b2_b(18, 18, this%control%lld, this%lattice%nrec))
         allocate (this%mu_n(18, 18, (2*this%lattice%control%lld) + 2, this%lattice%nrec), &
                   this%mu_ng(18, 18, (2*this%lattice%control%lld) + 2, this%lattice%nrec))
      else
         allocate (this%a_b(18, 18, this%control%lld, this%lattice%njij*4))
         allocate (this%b2_b(18, 18, this%control%lld, this%lattice%njij*4))
         allocate (this%mu_n(18, 18, (2*this%lattice%control%lld) + 2, this%lattice%njij*4), &
                   this%mu_ng(18, 18, (2*this%lattice%control%lld) + 2, this%lattice%njij*4))
      end if
      allocate (this%psi_b(18, 18, this%lattice%kk))
      allocate (this%hpsi(18, 18, this%lattice%kk))
      allocate (this%hohpsi(18, 18, this%lattice%kk))
      allocate (this%enupsi(18, 18, this%lattice%kk))
      allocate (this%socpsi(18, 18, this%lattice%kk))
      allocate (this%atemp_b(18, 18, this%control%lld))
      allocate (this%b2temp_b(18, 18, this%control%lld))
      allocate (this%pmn_b(18, 18, this%lattice%kk))
#endif
      this%v(:, :) = 0.0d0
      this%psi(:, :) = 0.0d0
      this%psi1(:, :, :) = 0.0d0
      this%psi2(:, :, :) = 0.0d0
      this%pmn(:, :) = 0.0d0
      this%mu_n(:, :, :, :) = 0.0d0
      this%mu_ng(:, :, :, :) = 0.0d0
      this%izero(:) = 0
      this%idum(:) = 0
      this%izeroll(:, :) = 0
      this%irlist(:) = 0
      this%a(:, :, :, :) = 0.0d0
      this%b2(:, :, :, :) = 0.0d0
      this%atemp(:) = 0.0d0
      this%b2temp(:) = 0.0d0
      this%a_b(:, :, :, :) = 0.0d0
      this%b2_b(:, :, :, :) = 0.0d0
      this%psi_b(:, :, :) = 0.0d0
      this%hpsi(:, :, :) = 0.0d0
      this%hohpsi(:, :, :) = 0.0d0
      this%enupsi(:, :, :) = 0.0d0
      this%socpsi(:, :, :) = 0.0d0
      this%atemp_b(:, :, :) = 0.0d0
      this%b2temp_b(:, :, :) = 0.0d0
      this%pmn_b(:, :, :) = 0.0d0
      this%cheb_mom_temp(:, :) = 0.0d0

      if (present(full)) then
         if (full) then
            if (associated(this%hamiltonian)) call this%hamiltonian%restore_to_default()
            if (associated(this%lattice)) call this%lattice%restore_to_default()
            if (associated(this%en)) call this%en%restore_to_default()
         end if
      end if

   end subroutine restore_to_default

end module recursion_mod

