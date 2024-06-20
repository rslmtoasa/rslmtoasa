!> \brief The Depondt solver for the stochastic LLG-equation
!! \details In principle the solver is of Heun type but uses rotations to
!! keep the magnitudes of the moments.
!! \details Ref: Ph. Depondt and F.G. Mertens, J. Phys.: Condens. Matter 21, 336005 (2009)
!! \todo Replace unit length moment vectors emom with full lenght vector emomM
module Depondt
   use Parameters
!  use Profiling
   !
   implicit none
   !
   real(dblprec), dimension(:, :, :), allocatable :: mrod !< Rotated magnetic moments
   real(dblprec), dimension(:, :, :), allocatable :: btherm !< Thermal stochastic field
   real(dblprec), dimension(:, :, :), allocatable :: bloc  !< Local effective field
   real(dblprec), dimension(:, :, :), allocatable :: bdup !< Resulting effective field

   private

   public :: depondt_evolve_first, depondt_evolve_second, allocate_depondtfields

contains

   !> First step of Depond solver, calculates the stochastic field and rotates the
   !! magnetic moments according to the effective field
   subroutine depondt_evolve_first(Natom, Mensemble, lambda1_array, beff, b2eff, &
                                   btorque, emom, emom2, emomM, mmom, delta_t, Temp_array, temprescale, stt, thermal_field, do_she, she_btorque)
      !
      use Constants, only: k_bolt, gama, mub
      use RandomNumbers, only: rng_gaussian, rng_gaussianP
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3, Natom, Mensemble), intent(in) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3, Natom, Mensemble), intent(out) :: b2eff !< Temporary storage of magnetic field
      real(dblprec), dimension(3, Natom, Mensemble), intent(inout) :: thermal_field
      real(dblprec), dimension(3, Natom, Mensemble), intent(in) :: btorque     !< Spin transfer torque
      real(dblprec), dimension(3, Natom, Mensemble), intent(in) :: she_btorque !< Spin Hall effect spin transfer torque
      real(dblprec), dimension(3, Natom, Mensemble), intent(inout) :: emom     !< Current unit moment vector
      real(dblprec), dimension(3, Natom, Mensemble), intent(out) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), dimension(3, Natom, Mensemble), intent(out) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(Natom, Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      real(dblprec), intent(in) :: delta_t !< Time step
      real(dblprec), dimension(Natom), intent(in) :: Temp_array !< Temperature (array)
      real(dblprec), intent(inout) :: temprescale  !< Temperature rescaling from QHB
      character(len=1), intent(in) :: STT    !< Treat spin transfer torque?
      character(len=1), intent(in) :: do_she !< Treat the spin hall effect transfer torque

      integer :: i, k, ik
      real(dblprec) :: v, Bnorm, hx, hy, hz
      real(dblprec) :: u, cosv, sinv, lldamp
      real(dblprec) :: sigma, Dp, she_fac, stt_fac

      if (stt /= 'N') then
         stt_fac = 1.0d0
         !$omp parallel do default(shared) private(i,k)  schedule(static) collapse(2)
         do k = 1, Mensemble
            do i = 1, Natom
               ! Adding STT and SHE torques if present (prefactor instead of if-statement)
               bdup(:, i, k) = stt_fac*btorque(:, i, k)
            end do
         end do
         !$omp end parallel do
      else
         stt_fac = 0.0d0
         !$omp parallel do default(shared) private(i,k)  schedule(static) collapse(2)
         do k = 1, Mensemble
            do i = 1, Natom
               bdup(:, i, k) = 0.0d0
            end do
         end do
         !$omp end parallel do
      end if

      if (do_she /= 'N') then
         she_fac = 1.0d0
         !$omp parallel do default(shared) private(i,k)  schedule(static) collapse(2)
         do k = 1, Mensemble
            do i = 1, Natom
               ! Adding STT and SHE torques if present (prefactor instead of if-statement)
               bdup(:, i, k) = bdup(:, i, k) + she_fac*she_btorque(:, i, k)
            end do
         end do
         !$omp end parallel do
      else
         she_fac = 0.0d0
      end if

      call rng_gaussianP(btherm, 3*Natom*Mensemble, 1.d0)

      ! Dupont recipe J. Phys.: Condens. Matter 21 (2009) 336005
      !!$omp parallel do default(shared) private(ik,i,k) schedule(static) collapse(2)
      !$omp parallel do default(shared) private(ik,i,k,Dp,sigma,Bnorm,hx,hy,hz,v,lldamp,cosv,sinv,u)  schedule(static) collapse(2)
      do k = 1, Mensemble
         do i = 1, Natom

            ! Thermal field
            !   LL equations ONE universal damping
            ! Midpoint
            !  Dk(:,i)=(lambda1_array(i)/(1+lambda1_array(i)**2)*k_bolt/gama/(mub))*(gama/bn)
            !  D=Dk(j,i)*mmomi(i,j)*Temp_array(i)*temprescale
            !  sigma=sqrt(2.0d0*D)
            !  sigma_mid=sqrt(2*alfa/(1+alfa**2)*k_bolt/gama/(mub))*(gama/bn)*T/mmom
            !  sigma_dep=sqrt(2*alfa*k_bolt)/(delta_t*gama*mub)*T/mmom)

            Dp = (2.d0*lambda1_array(i)*k_bolt)/(delta_t*gama*mub)   !LLG
            !print *,Dp
            !    Dp=(2.d0*lambda1_array*k_bolt)/(delta_t*gama*mub*(1+lambda1_array**2))  !LL or is it ? (code not consistent)
            !call rng_gaussian(btherm(1:3,i,k),3,1.d0)
            sigma = sqrt(Dp*temprescale*Temp_array(i)/mmom(i, k))
            btherm(:, i, k) = btherm(:, i, k)*sigma

            !print '(a,3f12.6)','btherm', btherm(:,i,k)
            ! Construct local field
            bloc(:, i, k) = beff(:, i, k) + btherm(:, i, k)
            !print '(2x,a,i4,3g18.6)','beff  ' ,i,bloc(1:3,i,k)
            thermal_field(:, i, k) = btherm(:, i, k)

            !print '(a,3f12.6)','b_loc ', bloc(:,i,k)

            ! Construct effective field (including damping term)
            bdup(1, i, k) = bdup(1, i, k) + bloc(1, i, k) + lambda1_array(i)*emom(2, i, k)*bloc(3, i, k) - lambda1_array(i)*emom(3, i, k)*bloc(2, i, k)
            bdup(2, i, k) = bdup(2, i, k) + bloc(2, i, k) + lambda1_array(i)*emom(3, i, k)*bloc(1, i, k) - lambda1_array(i)*emom(1, i, k)*bloc(3, i, k)
            bdup(3, i, k) = bdup(3, i, k) + bloc(3, i, k) + lambda1_array(i)*emom(1, i, k)*bloc(2, i, k) - lambda1_array(i)*emom(2, i, k)*bloc(1, i, k)

            !print '(a,3f12.6)','b_bdup', bdup(:,i,k)

            ! Set up rotation matrices and perform rotations
            lldamp = 1.0D0/(1.0D0 + lambda1_array(i)**2)
            Bnorm = bdup(1, i, k)**2 + bdup(2, i, k)**2 + bdup(3, i, k)**2
            Bnorm = sqrt(Bnorm) + 1.0d-15
            hx = bdup(1, i, k)/Bnorm
            hy = bdup(2, i, k)/Bnorm
            hz = bdup(3, i, k)/Bnorm
            v = Bnorm*delta_t*gama*lldamp
            cosv = cos(v)
            sinv = sin(v)
            u = 1.0d0 - cosv
            mrod(1, i, k) = hx*hx*u*emom(1, i, k) + cosv*emom(1, i, k) + hx*hy*u*emom(2, i, k) - hz*sinv*emom(2, i, k) + &
                            hx*hz*u*emom(3, i, k) + hy*sinv*emom(3, i, k)
            mrod(2, i, k) = hy*hx*u*emom(1, i, k) + hz*sinv*emom(1, i, k) + hy*hy*u*emom(2, i, k) + cosv*emom(2, i, k) + &
                            hy*hz*u*emom(3, i, k) - hx*sinv*emom(3, i, k)
            mrod(3, i, k) = hx*hz*u*emom(1, i, k) - hy*sinv*emom(1, i, k) + hz*hy*u*emom(2, i, k) + hx*sinv*emom(2, i, k) + &
                            hz*hz*u*emom(3, i, k) + cosv*emom(3, i, k)

            !print '(a,3f12.6)','m_rod ', mrod(:,i,k)

            ! copy m(t) to emom2 and m(t+dt) to emom for heisge, save b(t)
            emom2(:, i, k) = emom(:, i, k)
            !      do k=1,Mensemble
            !      do i=1,Natom
            emomM(:, i, k) = mrod(:, i, k)*mmom(i, k)
            !      end do
            !      end do

            emom(:, i, k) = mrod(:, i, k)

            b2eff(:, i, k) = bdup(:, i, k)

         end do
      end do
      !$omp end parallel do

   end subroutine depondt_evolve_first

   !> Second step of Depond solver, calculates the corrected effective field from
   !! the predicted effective fields. Rotates the moments in the corrected field
   subroutine depondt_evolve_second(Natom, Mensemble, lambda1_array, beff, b2eff, &
                                    btorque, emom, emom2, delta_t, stt, do_she, she_btorque)
      use Constants, only: gama
      !
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3, Natom, Mensemble), intent(in) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3, Natom, Mensemble), intent(in) :: b2eff !< Temporary storage of magnetic field
      real(dblprec), dimension(3, Natom, Mensemble), intent(in) :: btorque !< Spin transfer torque
      real(dblprec), dimension(3, Natom, Mensemble), intent(in) :: she_btorque !< SHE spin transfer torque
      real(dblprec), dimension(3, Natom, Mensemble), intent(out) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3, Natom, Mensemble), intent(out) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), intent(in) :: delta_t !< Time step
      character(len=1), intent(in) :: STT    !< Treat spin transfer torque?
      character(len=1), intent(in) :: do_she !< Treat the SHE spin transfer torque
      !
      integer :: i, k, ik
      real(dblprec) :: v, Bnorm, hx, hy, hz
      real(dblprec) :: u, cosv, sinv, lldamp, she_fac, stt_fac

      if (stt /= 'N') then
         stt_fac = 1.0d0
      else
         stt_fac = 0.0d0
      end if

      if (do_she /= 'N') then
         she_fac = 1.0d0
      else
         she_fac = 0.0d0
      end if

      if (stt == 'Y' .or. do_she == 'Y') then
         !$omp parallel do default(shared) private(i,k)  schedule(static) collapse(2)
         do k = 1, Mensemble
            do i = 1, Natom
               ! Adding STT and SHE torques if present (prefactor instead of if-statement)
               bdup(:, i, k) = stt_fac*btorque(:, i, k) + she_fac*she_btorque(:, i, k)
            end do
         end do
         !$omp end parallel do
      else
         !$omp parallel do default(shared) private(i,k)  schedule(static) collapse(2)
         do k = 1, Mensemble
            do i = 1, Natom
               bdup(:, i, k) = 0.0d0
            end do
         end do
         !$omp end parallel do
      end if

      !!$omp parallel do default(shared) private(ik,i,k) schedule(static) collapse(2)
      !$omp parallel do default(shared) private(ik,i,k,Bnorm,hx,hy,hz,v,lldamp,cosv,sinv,u) schedule(static) collapse(2)
      do k = 1, Mensemble
         do i = 1, Natom

            ! Construct local field
            bloc(:, i, k) = beff(:, i, k) + btherm(:, i, k)

            ! Construct effective field (including damping term)
            bdup(1, i, k) = bdup(1, i, k) + bloc(1, i, k) + lambda1_array(i)*emom(2, i, k)*bloc(3, i, k) - lambda1_array(i)*emom(3, i, k)*bloc(2, i, k)
            bdup(2, i, k) = bdup(2, i, k) + bloc(2, i, k) + lambda1_array(i)*emom(3, i, k)*bloc(1, i, k) - lambda1_array(i)*emom(1, i, k)*bloc(3, i, k)
            bdup(3, i, k) = bdup(3, i, k) + bloc(3, i, k) + lambda1_array(i)*emom(1, i, k)*bloc(2, i, k) - lambda1_array(i)*emom(2, i, k)*bloc(1, i, k)

            ! Corrected field
            bdup(:, i, k) = 0.5d0*bdup(:, i, k) + 0.5d0*b2eff(:, i, k)
            emom(:, i, k) = emom2(:, i, k)

            ! Set up rotation matrices and perform rotations
            lldamp = 1.0D0/(1.0D0 + lambda1_array(i)**2)
            Bnorm = bdup(1, i, k)**2 + bdup(2, i, k)**2 + bdup(3, i, k)**2
            Bnorm = sqrt(Bnorm) + 1.0d-15
            hx = bdup(1, i, k)/Bnorm
            hy = bdup(2, i, k)/Bnorm
            hz = bdup(3, i, k)/Bnorm
            v = Bnorm*delta_t*gama*lldamp
            cosv = cos(v)
            sinv = sin(v)
            u = 1.0d0 - cosv
            mrod(1, i, k) = hx*hx*u*emom(1, i, k) + cosv*emom(1, i, k) + hx*hy*u*emom(2, i, k) - hz*sinv*emom(2, i, k) + &
                            hx*hz*u*emom(3, i, k) + hy*sinv*emom(3, i, k)
            mrod(2, i, k) = hy*hx*u*emom(1, i, k) + hz*sinv*emom(1, i, k) + hy*hy*u*emom(2, i, k) + cosv*emom(2, i, k) + &
                            hy*hz*u*emom(3, i, k) - hx*sinv*emom(3, i, k)
            mrod(3, i, k) = hx*hz*u*emom(1, i, k) - hy*sinv*emom(1, i, k) + hz*hy*u*emom(2, i, k) + hx*sinv*emom(2, i, k) + &
                            hz*hz*u*emom(3, i, k) + cosv*emom(3, i, k)

            ! Final update
            emom2(:, i, k) = mrod(:, i, k)
            !print '(a,3f16.4)' ,'emom2 ', emom2(:,i,k)
         end do
      end do
      !$omp end parallel do

   end subroutine depondt_evolve_second

   !> Performs a Rodrigues rotation of the magnetic moments
   !! in the effective field.
   subroutine rodrigues(Natom, Mensemble, emom, delta_t, lambda1_array)
      !
      use Constants, only: gama
      !
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3, Natom, Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), intent(in) :: delta_t !< Time step
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec) :: v, Bnorm, hx, hy, hz
      real(dblprec) :: u, cosv, sinv, lldamp

      integer :: ik, i, k

      !$omp parallel do default(shared) private(ik,i,k,Bnorm,hx,hy,hz,v,lldamp,cosv,sinv,u) schedule(static) collapse(2)
      do k = 1, Mensemble
         do i = 1, Natom
            lldamp = 1.0D0/(1.0D0 + lambda1_array(i)**2)
            Bnorm = bdup(1, i, k)**2 + bdup(2, i, k)**2 + bdup(3, i, k)**2
            Bnorm = sqrt(Bnorm)
            hx = bdup(1, i, k)/Bnorm
            hy = bdup(2, i, k)/Bnorm
            hz = bdup(3, i, k)/Bnorm
            v = Bnorm*delta_t*gama*lldamp
            cosv = cos(v)
            sinv = sin(v)
            u = 1.0d0 - cosv
            mrod(1, i, k) = hx*hx*u*emom(1, i, k) + cosv*emom(1, i, k) + hx*hy*u*emom(2, i, k) - hz*sinv*emom(2, i, k) + &
                            hx*hz*u*emom(3, i, k) + hy*sinv*emom(3, i, k)
            mrod(2, i, k) = hy*hx*u*emom(1, i, k) + hz*sinv*emom(1, i, k) + hy*hy*u*emom(2, i, k) + cosv*emom(2, i, k) + &
                            hy*hz*u*emom(3, i, k) - hx*sinv*emom(3, i, k)
            mrod(3, i, k) = hx*hz*u*emom(1, i, k) - hy*sinv*emom(1, i, k) + hz*hy*u*emom(2, i, k) + hx*sinv*emom(2, i, k) + &
                            hz*hz*u*emom(3, i, k) + cosv*emom(3, i, k)
         end do
      end do
      !$omp end parallel do

   end subroutine rodrigues

   !> Calculates stochastic field
   subroutine thermfield(Natom, Mensemble, lambda1_array, mmom, deltat, Temp_array, temprescale)
      !
      use Constants, only: k_bolt, gama, mub
      use RandomNumbers, only: rng_gaussian

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(Natom, Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      real(dblprec), intent(in) :: deltat !< Time step
      real(dblprec), dimension(Natom), intent(in) :: Temp_array  !< Temperature (array)
      real(dblprec), intent(in) :: temprescale  !< Temperature rescaling from QHB

      real(dblprec), dimension(Natom) :: Dp
      real(dblprec) :: mu, sigma

      integer :: i, k

      !   LL equations ONE universal damping
      Dp = (2.d0*lambda1_array*k_bolt)/(deltat*gama*mub)   !LLG
      !    Dp=(2.d0*lambda1_array*k_bolt)/(deltat*gama*mub*(1+lambda1_array**2))  !LL or is it ? (code not consistent)

      !   LLG equations ONE universal damping
      call rng_gaussian(btherm, 3*Natom*Mensemble, 1.d0)
      mu = 0d0

      !$omp parallel do default(shared) private(k,i,sigma) collapse(2) schedule(static)
      do k = 1, Mensemble
         do i = 1, Natom
            sigma = sqrt(Dp(i)*temprescale*Temp_array(i)/mmom(i, k))
            btherm(:, i, k) = btherm(:, i, k)*sigma
         end do
      end do
      !$omp end parallel do

   end subroutine thermfield

   !> Constructs the effective field (including damping term)
   subroutine buildbeff(Natom, Mensemble, lambda1_array, emom, btorque, stt, do_she, she_btorque)

      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3, Natom, Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3, Natom, Mensemble), intent(in) :: btorque !< Spin transfer torque
      real(dblprec), dimension(3, Natom, Mensemble), intent(in) :: she_btorque !< SHE spin transfer torque
      character(len=1), intent(in) :: STT !< Treat spin transfer torque?
      character(len=1), intent(in) :: do_she !< Treat SHE spin transfer torque
      !
      integer :: i, k, ik

      !$omp parallel do default(shared) private(ik,i,k) schedule(static) collapse(2)
      do k = 1, Mensemble
         do i = 1, Natom
            bdup(1, i, k) = bloc(1, i, k) + lambda1_array(i)*emom(2, i, k)*bloc(3, i, k) - lambda1_array(i)*emom(3, i, k)*bloc(2, i, k)
            bdup(2, i, k) = bloc(2, i, k) + lambda1_array(i)*emom(3, i, k)*bloc(1, i, k) - lambda1_array(i)*emom(1, i, k)*bloc(3, i, k)
            bdup(3, i, k) = bloc(3, i, k) + lambda1_array(i)*emom(1, i, k)*bloc(2, i, k) - lambda1_array(i)*emom(2, i, k)*bloc(1, i, k)
         end do
      end do
      !$omp end parallel do

      if (stt /= 'N') then
         !$omp parallel do default(shared) private(ik,i,k) schedule(static) collapse(2)
         do k = 1, Mensemble
            do i = 1, Natom
               bdup(:, i, k) = bdup(:, i, k) + btorque(:, i, k)
            end do
         end do
         !$omp end parallel do
      end if

      if (do_she /= 'N') then
         !$omp parallel do default(shared) private(ik,i,k) schedule(static) collapse(2)
         do k = 1, Mensemble
            do i = 1, Natom
               bdup(:, i, k) = bdup(:, i, k) + she_btorque(:, i, k)
            end do
         end do
         !$omp end parallel do

      end if

   end subroutine buildbeff

   !> Allocates work arrays for the Depondt solver
   subroutine allocate_depondtfields(Natom, Mensemble, flag)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if (flag > 0) then
         allocate (bloc(3, Natom, Mensemble), stat=i_stat)
         allocate (btherm(3, Natom, Mensemble), stat=i_stat)
         allocate (bdup(3, Natom, Mensemble), stat=i_stat)
         allocate (mrod(3, Natom, Mensemble), stat=i_stat)
      else
         deallocate (bloc, stat=i_stat)
         deallocate (btherm, stat=i_stat)
         deallocate (bdup, stat=i_stat)
         deallocate (mrod, stat=i_stat)
      end if
   end subroutine allocate_depondtfields

!!!  !> First step of Depond solver, calculates the stochastic field and rotates the
!!!  !! magnetic moments according to the effective field
!!!  subroutine depondt_evolve_first_batched(Natom, Mensemble, lambda1_array, beff, b2eff, &
!!!  btorque, emom, emom2, emomM, mmom, delta_t, Temp_array, temprescale, stt,thermal_field,do_she,she_btorque)
!!!    !
!!!    implicit none
!!!    !
!!!    integer, intent(in) :: Natom !< Number of atoms in system
!!!    integer, intent(in) :: Mensemble !< Number of ensembles
!!!    real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
!!!    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field from application of Hamiltonian
!!!    real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: b2eff !< Temporary storage of magnetic field
!!!    real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: thermal_field
!!!    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: btorque     !< Spin transfer torque
!!!    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: she_btorque !< Spin Hall effect spin transfer torque
!!!    real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom     !< Current unit moment vector
!!!    real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom2  !< Final (or temporary) unit moment vector
!!!    real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emomM  !< Current magnetic moment vector
!!!    real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
!!!    real(dblprec), intent(in) :: delta_t !< Time step
!!!    real(dblprec), dimension(Natom), intent(in) :: Temp_array !< Temperature (array)
!!!    real(dblprec), intent(inout) :: temprescale  !< Temperature rescaling from QHB
!!!    character(len=1), intent(in) :: STT    !< Treat spin transfer torque?
!!!    character(len=1), intent(in) :: do_she !< Treat the spin hall effect transfer torque
!!!
!!!    integer :: i,k,ik
!!!
!!!    ! Dupont recipe J. Phys.: Condens. Matter 21 (2009) 336005
!!!    ! Calculate stochastic field
!!!    call timing(0,'Evolution     ','OF')
!!!    call timing(0,'RNG           ','ON')
!!!
!!!    call thermfield(Natom, Mensemble, lambda1_array, mmom, delta_t, Temp_array,temprescale)
!!!
!!!    call timing(0,'RNG           ','OF')
!!!    call timing(0,'Evolution     ','ON')
!!!    ! Construct local field
!!!    !!!$omp parallel do default(shared) private(ik,i,k) schedule(static)
!!!    !!do ik=1,Natom*Mensemble
!!!    !!   i=mod(ik-1,Natom)+1
!!!    !!   k=int((ik-1)/Natom)+1
!!!    !$omp parallel do default(shared) private(ik,i,k) schedule(static) collapse(2)
!!!    do k=1,Mensemble
!!!    do i=1,Natom
!!!       !write(2630,'(1x,3f18.10,5x,4f18.10)') beff(1:3,i,k),btherm(1:3,i,k), sum(btherm(1:3,i,k)*btherm(1:3,i,k))**0.5d0
!!!       bloc(:,i,k)=beff(:,i,k)+btherm(:,i,k)
!!!       thermal_field(:,i,k)=btherm(:,i,k)
!!!    end do
!!!    end do
!!!    !$omp end parallel do
!!!
!!!    ! Construct effective field (including damping term)
!!!    call buildbeff(Natom, Mensemble,lambda1_array,emom,btorque,stt,do_she,she_btorque)
!!!
!!!    ! Set up rotation matrices and perform rotations
!!!    call rodrigues(Natom, Mensemble,emom, delta_t,lambda1_array)
!!!
!!!    ! copy m(t) to emom2 and m(t+dt) to emom for heisge, save b(t)
!!!    !!$omp parallel do default(shared) private(ik,i,k) schedule(static)
!!!    !do ik=1,Natom*Mensemble
!!!    !   i=mod(ik-1,Natom)+1
!!!    !   k=int((ik-1)/Natom)+1
!!!    !$omp parallel do default(shared) private(ik,i,k) schedule(static), collapse(2)
!!!    do k=1,Mensemble
!!!       do i=1,Natom
!!!       emom2(:,i,k)=emom(:,i,k)
!!!       !      do k=1,Mensemble
!!!       !      do i=1,Natom
!!!       emomM(:,i,k)=mrod(:,i,k)*mmom(i,k)
!!!       !      end do
!!!       !      end do
!!!
!!!       emom(:,i,k)=mrod(:,i,k)
!!!
!!!       b2eff(:,i,k)=bdup(:,i,k)
!!!
!!!       end do
!!!    end do
!!!    !$omp end parallel do
!!!
!!!  end subroutine depondt_evolve_first_batched
!!!
!!!
!!!  !> Second step of Depond solver, calculates the corrected effective field from
!!!  !! the predicted effective fields. Rotates the moments in the corrected field
!!!  subroutine depondt_evolve_second_batched(Natom, Mensemble, lambda1_array, beff, b2eff, &
!!!       btorque, emom, emom2, delta_t, stt,do_she,she_btorque)
!!!    !
!!!    implicit none
!!!    !
!!!    integer, intent(in) :: Natom !< Number of atoms in system
!!!    integer, intent(in) :: Mensemble !< Number of ensembles
!!!    real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
!!!    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field from application of Hamiltonian
!!!    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: b2eff !< Temporary storage of magnetic field
!!!    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: btorque !< Spin transfer torque
!!!    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: she_btorque !< SHE spin transfer torque
!!!    real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom   !< Current unit moment vector
!!!    real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom2  !< Final (or temporary) unit moment vector
!!!    real(dblprec), intent(in) :: delta_t !< Time step
!!!    character(len=1), intent(in) :: STT    !< Treat spin transfer torque?
!!!    character(len=1), intent(in) :: do_she !< Treat the SHE spin transfer torque
!!!    !
!!!    integer :: i,k,ik
!!!
!!!    ! Construct local field
!!!    !!!$omp parallel do default(shared) private(ik,i,k) schedule(static)
!!!    !!do ik=1,Natom*Mensemble
!!!    !!   i=mod(ik-1,Natom)+1
!!!    !!   k=int((ik-1)/Natom)+1
!!!    !$omp parallel do default(shared) private(ik,i,k) schedule(static) collapse(2)
!!!    do k=1,Mensemble
!!!       do i=1,Natom
!!!       bloc(:,i,k)=beff(:,i,k)+btherm(:,i,k)
!!!      end do
!!!    end do
!!!    !$omp end parallel do
!!!
!!!    ! Construct effective field (including damping term)
!!!    call buildbeff(Natom, Mensemble,lambda1_array,emom, btorque,stt,do_she,she_btorque)
!!!
!!!    ! Corrected field
!!!    !!!$omp parallel do default(shared) private(ik,i,k) schedule(static)
!!!    !!do ik=1,Natom*Mensemble
!!!    !!   i=mod(ik-1,Natom)+1
!!!    !!   k=int((ik-1)/Natom)+1
!!!    !$omp parallel do default(shared) private(ik,i,k) schedule(static) collapse(2)
!!!    do k=1,Mensemble
!!!    do i=1,Natom
!!!       bdup(:,i,k)=0.5d0*bdup(:,i,k)+0.5d0*b2eff(:,i,k)
!!!       emom(:,i,k)=emom2(:,i,k)
!!!    end do
!!!    end do
!!!    !$omp end parallel do
!!!
!!!    ! Final rotation
!!!    call rodrigues(Natom, Mensemble,emom, delta_t,lambda1_array)
!!!
!!!    !!!$omp parallel do default(shared) private(ik,i,k) schedule(static)
!!!    !!do ik=1,Natom*Mensemble
!!!    !!   i=mod(ik-1,Natom)+1
!!!    !!   k=int((ik-1)/Natom)+1
!!!    !$omp parallel do default(shared) private(ik,i,k) schedule(static) collapse(2)
!!!    do k=1,Mensemble
!!!      do i=1,Natom
!!!       emom2(:,i,k)=mrod(:,i,k)
!!!      end do
!!!    end do
!!!    !$omp end parallel do
!!!
!!!  end subroutine depondt_evolve_second_batched

end module depondt
