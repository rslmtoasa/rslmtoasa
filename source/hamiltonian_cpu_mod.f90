!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Hamiltonian CPU
!
!> @brief
!> CPU backend implementation for polymorphic Hamiltonian interface.
!> Handles sparse matrix-vector multiplication H|Psi> with optional
!> orthogonalization correction (HOH), spin-orbit coupling, and linearization
!> energy contributions. Manages active-site mask updates (growing shell).
!
!> @details
!> This module provides a concrete implementation of the abstract
!> hamiltonian_type interface. The apply subroutine encapsulates all
!> sparse-matrix vector multiplication logic, including:
!> - Local region (1..nmax) and bulk region (nlimplus1..kk) separation
!> - OpenMP parallelization over atoms with dynamic scheduling
!> - Optional HOH double-step with temporary buffers (hohpsi, enupsi, socpsi)
!> - Chebyshev scaling and shift: psi_out = (psi_out - b*psi_in) / a
!> - Active-site mask (izero/idum) update for recursion growth
!
!------------------------------------------------------------------------------

module hamiltonian_cpu_mod

   use hamiltonian_base_mod
   use hamiltonian_mod
   use lattice_mod
   use precision_mod, only: rp
   use basis_mod, only: nb
   implicit none

   private
   public :: cpu_hamiltonian, init_cpu_hamiltonian

   !> CPU-based Hamiltonian type that wraps the real TB-LMTO arrays
   type, extends(hamiltonian_type), public :: cpu_hamiltonian
      !> Pointer to the actual TB-LMTO Hamiltonian data
      type(hamiltonian), pointer :: h      => null()
      !> Pointer to lattice data (neighbor list, atom types, etc.)
      type(lattice),     pointer :: lat    => null()
      !> Pointer to active-site mask (modified in-place during recursion)
      integer, pointer           :: izero(:) => null()
      !> Pointer to temporary active-site mask (updated during apply)
      integer, pointer           :: idum(:)  => null()
   contains
      procedure, pass :: apply => cpu_apply
   end type cpu_hamiltonian

contains

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Initialize a CPU Hamiltonian wrapper
   !
   !> @param[out] cpu_h        Target cpu_hamiltonian object to initialize
   !> @param[in]  h            Pointer to TB-LMTO Hamiltonian
   !> @param[in]  lat          Pointer to lattice
   !> @param[in]  izero        Pointer to active-site mask
   !> @param[in]  idum         Pointer to temporary active-site mask
   !> @param[in]  hoh_enabled  Flag for HOH double-step
   !---------------------------------------------------------------------------
   subroutine init_cpu_hamiltonian(cpu_h, h, lat, izero, idum, hoh_enabled)
      type(cpu_hamiltonian), intent(out)  :: cpu_h
      type(hamiltonian), target           :: h
      type(lattice), target               :: lat
      integer, target                     :: izero(:), idum(:)
      logical, intent(in)                 :: hoh_enabled

      cpu_h%h => h
      cpu_h%lat => lat
      cpu_h%izero => izero
      cpu_h%idum => idum
      cpu_h%hoh_enabled = hoh_enabled
   end subroutine init_cpu_hamiltonian

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Apply the Hamiltonian with optional HOH correction and Chebyshev scaling
   !
   !> @param[in]     this     CPU Hamiltonian object
   !> @param[in]     psi_in   Input wave function (nb, nb, kk)
   !> @param[out]    psi_out  Output wave function after H|psi> (nb, nb, kk)
   !> @param[in]     a        Chebyshev scaling factor (typically from polynomial)
   !> @param[in]     b        Chebyshev shift parameter (typically from polynomial)
   !
   !> @details
   !> Computes psi_out = (H|psi_in> - b*psi_in) / a
   !>
   !> The H application proceeds in two phases:
   !> 1. Local region (1..nmax) + Bulk region (nlimplus1..kk) via zgemm
   !> 2. Optional HOH second-pass for orthogonalization correction
   !>
   !> Accumulates enupsi (energy parameter) and socpsi (spin-orbit) separately
   !> to allow flexible assembly: H = h - hoh + e_nu + l.s
   !>
   !> Updates izero(:) = idum(:) at the end to reflect grown active shell.
   !---------------------------------------------------------------------------
   subroutine cpu_apply(this, psi_in, psi_out, a, b)
      class(cpu_hamiltonian), intent(inout) :: this
      complex(rp), intent(in)               :: psi_in(:,:,:)
      complex(rp), intent(out)              :: psi_out(:,:,:)
      real(rp), intent(in)                  :: a, b

      ! Local variables
      integer :: i, j, k, ino, ih, nr, nnmap, nlimplus1, kk
      complex(rp), dimension(nb, nb) :: locham

      ! Temporary buffers for HOH double-step (allocated locally, not module-level)
      complex(rp), allocatable :: hohpsi(:,:,:), enupsi(:,:,:), socpsi(:,:,:)

      complex(rp), parameter :: cone = (1.0_rp, 0.0_rp)

      kk = this%lat%kk
      nlimplus1 = this%lat%nmax + 1

      ! Allocate temporary buffers
      allocate(hohpsi(nb, nb, kk))
      allocate(enupsi(nb, nb, kk))
      allocate(socpsi(nb, nb, kk))

      ! Initialize all buffers and output
      hohpsi(:,:,:) = (0.0_rp, 0.0_rp)
      enupsi(:,:,:) = (0.0_rp, 0.0_rp)
      socpsi(:,:,:) = (0.0_rp, 0.0_rp)
      psi_out(:,:,:) = (0.0_rp, 0.0_rp)

      ! Initialize idum to zero for active-site growth tracking
      this%idum(:) = 0

      ! ==========================================
      ! PHASE 1: Apply H core (h + e_nu + l.s if no HOH, or just h if HOH)
      ! ==========================================

      if (this%lat%nmax /= 0) then
         ! Local region: use site-indexed Hamiltonian (1..nmax)
         !$omp parallel do default(shared) private(i, ino, nr, j, nnmap, locham) schedule(dynamic, 100)
         do i = 1, this%lat%nmax
            this%idum(i) = this%izero(i)
            ino = this%lat%iz(i)
            nr = this%lat%nn(i, 1)

            if (this%izero(i) /= 0) then
               if (this%hoh_enabled) then
                  ! HOH case: just apply h (hopping), accumulate e_nu and l.s separately
                  call zgemm('n', 'n', nb, nb, nb, cone, this%h%hall(:, :, 1, i), nb, psi_in(:, :, i), nb, cone, psi_out(:, :, i), nb)
                  call zgemm('n', 'n', nb, nb, nb, cone, this%h%enim(:, :, ino), nb, psi_in(:, :, i), nb, cone, enupsi(:, :, i), nb)
                  call zgemm('n', 'n', nb, nb, nb, cone, this%h%lsham(:, :, ino), nb, psi_in(:, :, i), nb, cone, socpsi(:, :, i), nb)
               else
                  ! No HOH: apply full h + e_nu + l.s
                  locham = this%h%hall(:, :, 1, i) + this%h%enim(:, :, ino) + this%h%lsham(:, :, ino)
                  call zgemm('n', 'n', nb, nb, nb, cone, locham, nb, psi_in(:, :, i), nb, cone, psi_out(:, :, i), nb)
               end if
            end if

            if (nr >= 2) then
               ! Hopping terms
               do j = 2, nr
                  nnmap = this%lat%nn(i, j)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        call zgemm('n', 'n', nb, nb, nb, cone, this%h%hall(:, :, j, i), nb, psi_in(:, :, nnmap), nb, cone, psi_out(:, :, i), nb)
                        this%idum(i) = 1
                     end if
                  end if
               end do
            end if
         end do
         !$omp end parallel do
      end if

      ! Bulk region: use type-indexed Hamiltonian (nlimplus1..kk)
      !$omp parallel do default(shared) private(i, ih, nr, j, nnmap, locham) schedule(dynamic, 100)
      do i = nlimplus1, kk
         this%idum(i) = this%izero(i)
         ih = this%lat%iz(i)
         nr = this%lat%nn(i, 1)

         if (this%izero(i) /= 0) then
            if (this%hoh_enabled) then
               ! HOH case: just apply h, accumulate e_nu and l.s separately
               call zgemm('n', 'n', nb, nb, nb, cone, this%h%ee(:, :, 1, ih), nb, psi_in(:, :, i), nb, cone, psi_out(:, :, i), nb)
               call zgemm('n', 'n', nb, nb, nb, cone, this%h%enim(:, :, ih), nb, psi_in(:, :, i), nb, cone, enupsi(:, :, i), nb)
               call zgemm('n', 'n', nb, nb, nb, cone, this%h%lsham(:, :, ih), nb, psi_in(:, :, i), nb, cone, socpsi(:, :, i), nb)
            else
               ! No HOH: apply full h + e_nu + l.s
               locham = this%h%ee(:, :, 1, ih) + this%h%enim(:, :, ih) + this%h%lsham(:, :, ih)
               call zgemm('n', 'n', nb, nb, nb, cone, locham, nb, psi_in(:, :, i), nb, cone, psi_out(:, :, i), nb)
            end if
         end if

         if (nr >= 2) then
            ! Hopping terms
            do j = 2, nr
               nnmap = this%lat%nn(i, j)
               if (nnmap /= 0) then
                  if (this%izero(nnmap) /= 0) then
                     call zgemm('n', 'n', nb, nb, nb, cone, this%h%ee(:, :, j, ih), nb, psi_in(:, :, nnmap), nb, cone, psi_out(:, :, i), nb)
                     this%idum(i) = 1
                  end if
               end if
            end do
         end if
      end do
      !$omp end parallel do

      ! ==========================================
      ! PHASE 2: HOH double-step (if enabled)
      ! ==========================================
      if (this%hoh_enabled) then

         ! Update izero to reflect newly reached sites before HOH pass
         this%izero(:) = this%idum(:)
         this%idum(:) = 0

         ! Local region: apply hallo to psi_out (which contains h|psi>)
         if (this%lat%nmax /= 0) then
            !$omp parallel do default(shared) private(i, ino, nr, j, nnmap) schedule(dynamic, 100)
            do i = 1, this%lat%nmax
               this%idum(i) = this%izero(i)
               ino = this%lat%iz(i)
               nr = this%lat%nn(i, 1)

               if (this%izero(i) /= 0) then
                  call zgemm('n', 'n', nb, nb, nb, cone, this%h%hallo(:, :, 1, i), nb, psi_out(:, :, i), nb, cone, hohpsi(:, :, i), nb)
               end if

               if (nr >= 2) then
                  do j = 2, nr
                     nnmap = this%lat%nn(i, j)
                     if (nnmap /= 0) then
                        if (this%izero(nnmap) /= 0) then
                           call zgemm('n', 'n', nb, nb, nb, cone, this%h%hallo(:, :, j, i), nb, psi_out(:, :, nnmap), nb, cone, hohpsi(:, :, i), nb)
                           this%idum(i) = 1
                        end if
                     end if
                  end do
               end if
            end do
            !$omp end parallel do
         end if

         ! Bulk region: apply eeo to psi_out (which contains h|psi>)
         !$omp parallel do default(shared) private(i, ih, nr, j, nnmap) schedule(dynamic, 100)
         do i = nlimplus1, kk
            this%idum(i) = this%izero(i)
            ih = this%lat%iz(i)
            nr = this%lat%nn(i, 1)

            if (this%izero(i) /= 0) then
               call zgemm('n', 'n', nb, nb, nb, cone, this%h%eeo(:, :, 1, ih), nb, psi_out(:, :, i), nb, cone, hohpsi(:, :, i), nb)
            end if

            if (nr >= 2) then
               do j = 2, nr
                  nnmap = this%lat%nn(i, j)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        call zgemm('n', 'n', nb, nb, nb, cone, this%h%eeo(:, :, j, ih), nb, psi_out(:, :, nnmap), nb, cone, hohpsi(:, :, i), nb)
                        this%idum(i) = 1
                     end if
                  end if
               end do
            end if
         end do
         !$omp end parallel do

         ! Final assembly with HOH correction: H = h - hoh + e_nu + l.s
         psi_out(:,:,:) = psi_out(:,:,:) - hohpsi(:,:,:) + enupsi(:,:,:) + socpsi(:,:,:)

         ! Update izero one final time
         this%izero(:) = this%idum(:)

      else
         ! No HOH: just update izero from the first pass
         this%izero(:) = this%idum(:)
      end if

      ! ==========================================
      ! PHASE 3: Chebyshev scale and shift
      ! ==========================================
      psi_out(:,:,:) = (psi_out(:,:,:) - b * psi_in(:,:,:)) / a

      ! Clean up temporary buffers
      deallocate(hohpsi, enupsi, socpsi)

   end subroutine cpu_apply

end module hamiltonian_cpu_mod
