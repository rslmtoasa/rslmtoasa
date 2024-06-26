!> Library for random number generators
#ifdef VSL
include 'mkl_vsl.f90'
#endif
module RandomNumbers
   use Parameters
   use mtprng
   use stdtypes
   use omp_lib

#ifdef VSL
   USE MKL_VSL_TYPE
   USE MKL_VSL
#endif

   implicit none

   real(dblprec), dimension(:, :, :), allocatable :: ranv !< Work array for RNG
   real(dblprec), dimension(:, :, :), allocatable :: lattranv !< Work array for RNG

   logical*1 :: use_vsl       !< Use Intel Vector Statistics Library (VSL) (T/F) (also need preprocessing flag VSL)

   !States for RNGs
   type(mtprng_state), SAVE :: state_a, state_b, state_c
   type(mtprng_state), SAVE :: state_z
   !$omp threadprivate(state_a,state_b,state_c,state_z)

#ifdef VSL
   !VSL RNG parameters
   TYPE(VSL_STREAM_STATE) :: stream(0:272)
!   !$omp threadprivate(stream)
#endif

   !Work arrays for ziggurat
   integer(INT32) :: kn(128)
   real(IEEE64) :: fn(128)
   real(IEEE64) :: wn(128)

   !Parameter for ziggurat
   real(IEEE64), parameter :: r = 3.442620D+00

   !Choice of transform (0=Box-Muller, 1=Ziggurat)
   integer :: rng_trans
   integer :: rng_sphere

   private

   public :: rng_uniform, rng_gaussian, rng_norm, rannum, allocate_randomwork
   public :: rng_defaults, setup_rng_hb, rng_init, fill_rngarray, rng_uniformP
   public :: ranv, use_vsl, rng_gaussianP
   public :: lattranv, lattrannum, mt_ran_init_c

contains

   ! Sets defaults for RNGs
   !> Set default parameters for RNG
   subroutine rng_defaults()
#ifdef VSL
      use_vsl = .true.
#else
      use_vsl = .false.
#endif
   end subroutine rng_defaults

   !> First Park/Miller RNG
   FUNCTION ran2_a(idum)
      IMPLICIT NONE
      INTEGER, PARAMETER :: K4B = selected_int_kind(9)
      INTEGER(K4B), INTENT(INOUT) :: idum
      REAL(dblprec) :: ran2_a
      !   "Minimal" random number generator of Park and Miller combined with a Marsaglia shift
      !   sequence. Returns a uniform random deviate in the interval [0.0,1.0)
      !   This fully portable, scalar generator has the "traditional" (not Fortran 90) calling
      INTEGER(K4B), PARAMETER :: IA = 16807, IM = 2147483647, IQ = 127773, IR = 2836
      REAL(dblprec), SAVE :: am
      INTEGER(K4B), SAVE :: ix = -1, iy = -1, k
      !Initialize.
      if (idum <= 0 .or. iy < 0) then
         am = nearest(1.0, -1.0)/IM
         iy = ior(ieor(888889999, abs(idum)), 1)
         ix = ieor(777755555, abs(idum))
         !Set idum positive.
         idum = abs(idum) + 1
      end if
      !Marsaglia shift sequence with period 232 - 1.
      ix = ieor(ix, ishft(ix, 13))
      ix = ieor(ix, ishft(ix, -17))
      ix = ieor(ix, ishft(ix, 5))
      !Park-Miller sequence by Schrage´s method,
      k = iy/IQ
      !period 231 - 2.
      iy = IA*(iy - k*IQ) - IR*k
      if (iy < 0) iy = iy + IM
      !Combine the two generators with masking to
      ran2_a = am*ior(iand(IM, ieor(ix, iy)), 1)
      !ensure nonzero value.
   END FUNCTION ran2_a

   !.. Initilizes the function ran
   !.. seed can be a positive or negative number
   !> Initialize first Park/Miller RNG
   SUBROUTINE ran2_init_a(seed)
      IMPLICIT NONE
      INTEGER(selected_int_kind(9)), INTENT(in) :: seed
      INTEGER(selected_int_kind(9)) :: idum
      REAL(dblprec) :: rdum
      ! REAL(dblprec) :: ran_a

      !Executive statements
      idum = -abs(seed)
      rdum = ran2_a(idum)
   END SUBROUTINE ran2_init_a

   !> Second Park/Miller RNG
   FUNCTION ran2_b(idum)
      IMPLICIT NONE
      INTEGER, PARAMETER :: K4B = selected_int_kind(9)
      INTEGER(K4B), INTENT(INOUT) :: idum
      REAL(dblprec) :: ran2_b
      !   "Minimal" random number generator of Park and Miller combined with a Marsaglia shift
      !   sequence. Returns a uniform random deviate in the interval [0.0,1.0)
      !   This fully portable, scalar generator has the "traditional" (not Fortran 90) calling
      INTEGER(K4B), PARAMETER :: IA = 16807, IM = 2147483647, IQ = 127773, IR = 2836
      REAL(dblprec), SAVE :: am
      INTEGER(K4B), SAVE :: ix = -1, iy = -1, k
      !Initialize.
      if (idum <= 0 .or. iy < 0) then
         am = nearest(1.0, -1.0)/IM
         iy = ior(ieor(888889999, abs(idum)), 1)
         ix = ieor(777755555, abs(idum))
         !Set idum positive.
         idum = abs(idum) + 1
      end if
      !Marsaglia shift sequence with period 232 - 1.
      ix = ieor(ix, ishft(ix, 13))
      ix = ieor(ix, ishft(ix, -17))
      ix = ieor(ix, ishft(ix, 5))
      !Park-Miller sequence by Schrage´s method,
      k = iy/IQ
      !period 231 - 2.
      iy = IA*(iy - k*IQ) - IR*k
      if (iy < 0) iy = iy + IM
      !Combine the two generators with masking to
      ran2_b = am*ior(iand(IM, ieor(ix, iy)), 1)
      !ensure nonzero value.
   END FUNCTION ran2_b

   !.. Initilizes the function ran
   !.. seed can be a positive or negative number
   !> Initialize second Park/Miller RNG
   SUBROUTINE ran2_init_b(seed)
      IMPLICIT NONE
      INTEGER(selected_int_kind(9)), INTENT(in) :: seed
      INTEGER(selected_int_kind(9)) :: idum
      REAL(dblprec) :: rdum
      ! REAL(dblprec) :: ran_b

      !Executive statements
      idum = -abs(seed)
      rdum = ran2_b(idum)
   END SUBROUTINE ran2_init_b

   !> Third Park/Miller RNG
   FUNCTION ran2_c(idum)
      IMPLICIT NONE
      INTEGER, PARAMETER :: K4B = selected_int_kind(9)
      INTEGER(K4B), INTENT(INOUT) :: idum
      REAL(dblprec) :: ran2_c
      !   "Minimal" random number generator of Park and Miller combined with a Marsaglia shift
      !   sequence. Returns a uniform random deviate in the interval [0.0,1.0)
      !   This fully portable, scalar generator has the "traditional" (not Fortran 90) calling
      INTEGER(K4B), PARAMETER :: IA = 16807, IM = 2147483647, IQ = 127773, IR = 2836
      REAL(dblprec), SAVE :: am
      INTEGER(K4B), SAVE :: ix = -1, iy = -1, k
      !Initialize.
      if (idum <= 0 .or. iy < 0) then
         am = nearest(1.0, -1.0)/IM
         iy = ior(ieor(888889999, abs(idum)), 1)
         ix = ieor(777755555, abs(idum))
         !Set idum positive.
         idum = abs(idum) + 1
      end if
      !Marsaglia shift sequence with period 232 - 1.
      ix = ieor(ix, ishft(ix, 13))
      ix = ieor(ix, ishft(ix, -17))
      ix = ieor(ix, ishft(ix, 5))
      !Park-Miller sequence by Schrage´s method,
      k = iy/IQ
      !period 231 - 2.
      iy = IA*(iy - k*IQ) - IR*k
      if (iy < 0) iy = iy + IM
      !Combine the two generators with masking to
      ran2_c = am*ior(iand(IM, ieor(ix, iy)), 1)
      !ensure nonzero value.
   END FUNCTION ran2_c

   !.. Initilizes the function ran
   !.. seed can be a positive or negative number
   !> Initialize third Park/Miller RNG
   SUBROUTINE ran2_init_c(seed)
      IMPLICIT NONE
      INTEGER(selected_int_kind(9)), INTENT(in) :: seed
      INTEGER(selected_int_kind(9)) :: idum
      REAL(dblprec) :: rdum
      ! REAL(dblprec) :: ran_c

      !Executive statements
      idum = -abs(seed)
      rdum = ran2_c(idum)
   END SUBROUTINE ran2_init_c

   !!!!!++ Does this one even work correct? Using give strange results in basic BccFe
   !> Wrapper function, using ran_a, to return a normally distributed deviate with zero mean and unit variance
   FUNCTION gasdev(mu, sigma)
      IMPLICIT NONE
      ! INTEGER(Selected_int_kind(9)), PARAMETER :: idum = 1
      ! INTEGER, PARAMETER :: P15 = KIND(dblprec)
      REAL(dblprec), INTENT(IN) :: mu, sigma
      REAL(dblprec) :: gasdev
      REAL(dblprec), parameter :: one = 1.0_dblprec
      REAL(dblprec), parameter :: two = 2.0_dblprec
      !   Returns in harvest a normally distributed deviate with zero mean and unit variance, using
      !   ran0 as the source of uniform deviates.
      REAL(dblprec) :: rsq, v1, v2
      ! REAL(dblprec) :: ran_a
      REAL(dblprec), SAVE :: g
      LOGICAL, SAVE :: gaus_stored = .false.

      if (gaus_stored) then
         !     !We have an extra deviate handy,
         !     !so return it,
         !and unset the flag.
         gaus_stored = .false.
         gasdev = g
      else
         !We don´t have an extra deviate handy, so
         !pick two uniform numbers in the square ex-
         !tending from -1 to +1 in each direction,
         do
            v1 = mt_ran_a()
            v2 = mt_ran_a()
            !       v1=2.0_P15*v1-1.0_P15
            !       v2=2.0_P15*v2-1.0_P15
            v1 = two*v1 - one
            v2 = two*v2 - one
            !see if they are in the unit circle,
            !otherwise try again.
            rsq = v1**2 + v2**2
            if (rsq > 0.0 .and. rsq < 1.0) exit
         end do
         !Now make the Box-Muller transformation to
         !get two normal deviates. Return one and
         !save the other for next time.
         !    rsq=sqrt(-2.0_P15*log(rsq)/rsq)
         rsq = sqrt(-two*log(rsq)/rsq)
         gasdev = v1*rsq*sigma + mu
         g = v2*rsq*sigma + mu
         !Set flag.
         gaus_stored = .true.
      end if
   END FUNCTION gasdev

   !> Allocate work arrays for RNG
   subroutine allocate_randomwork(Natom, Mensemble, flag, do_ld)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)
      character(len=1) :: do_ld !< Do lattice dynamics ('Y'/'N')

      integer :: i_stat, i_all

      if (flag > 0) then
         allocate (ranv(3, Natom, Mensemble), stat=i_stat)
         if (do_ld == 'Y') then
            allocate (lattranv(3, Natom, Mensemble), stat=i_stat)
         end if
      else
         deallocate (ranv, stat=i_stat)
         if (do_ld == 'Y') then
            deallocate (lattranv, stat=i_stat)
         end if
      end if
   end subroutine allocate_randomwork

   ! Wrappers for Mersenne twister generators using VSL library
   !> Initialize the Mersenne Twister generator
   subroutine rng_init(seed)
      integer, intent(in) :: seed  !< Seed number for PRNG
#ifdef VSL
      integer :: errcode, id, nid
#endif

#ifdef VSL
      if (use_vsl) then
         id = omp_get_thread_num()
!         nid=omp_get_num_threads()
         errcode = vslnewstream(stream(id), VSL_BRNG_MT2203 + id, seed)
!             errcode=vslnewstream( stream(id), VSL_BRNG_SFMT19937, seed )
         !    errcode= vslleapfrogstream( stream(id), id, nid)
         !    errcode= vslskipaheadstream( stream(id), 16)
      else
         call mtprng_init(seed, state_c)
      end if
#else
      call mtprng_init(seed, state_c)
#endif
   end subroutine rng_init

   !> Wrapper for VSL Mersenne Twister generator - uniform distribution
   subroutine rng_uniform(out, len)
      integer, intent(in) :: len
      REAL(dblprec), dimension(len), intent(out) :: out  !< Random number in the interval [0,1)
      integer :: i
#ifdef VSL
      integer :: errcode
#endif

#ifdef VSL
      if (use_vsl) then
         errcode = vdrnguniform(VSL_RNG_METHOD_UNIFORM_STD, stream(omp_get_thread_num()), len, out, 0.d0, 1.d0)
      else
         do i = 1, len
            out(i) = mtprng_rand_real2(state_c)
         end do
      end if
#else
      do i = 1, len
         out(i) = mtprng_rand_real2(state_c)
      end do
#endif
   end subroutine rng_uniform

!!!    !> Wrapper for VSL Mersenne Twister generator - uniform distribution
   subroutine rng_uniformP(out, len)
      integer, intent(in) :: len
      REAL(dblprec), dimension(len), intent(out) :: out  !< Random number in the interval [0,1)
      integer :: i
!!! #ifdef VSL
!!!       integer :: errcode
!!! #endif
!!!       integer :: iproc,id,local_len,start_idx,stop_idx
!!!
!!! #ifdef VSL
!!!        if(use_vsl) then
!!!           !$omp parallel do default(shared) schedule(static,1) &
!!!           !$omp private(errcode,iproc,id,local_len,start_idx,stop_idx)
!!!           do iproc=1,omp_get_num_threads()
!!!              id=omp_get_thread_num()
!!!              local_len=len/omp_get_num_threads()
!!!              start_idx=id*local_len+1
!!!              stop_idx=min(len,(id+1)*local_len)
!!!              if (id==(omp_get_num_threads()-1)) stop_idx=len
!!!              local_len=stop_idx-start_idx+1
!!!              errcode=vdrnguniform(VSL_RNG_METHOD_UNIFORM_STD, stream(id), local_len, out(start_idx:stop_idx), 0.d0, 1.d0)
!!!           end do
!!!           !$omp end parallel do
!!!        else
!!!           !$omp parallel do default(shared) schedule(static,1) &
!!!           !$omp private(i,iproc,id,local_len,start_idx,stop_idx)
!!!           do iproc=1,omp_get_num_threads()
!!!              id=omp_get_thread_num()
!!!              local_len=len/omp_get_num_threads()
!!!              start_idx=id*local_len+1
!!!              stop_idx=min(len,(id+1)*local_len)
!!!              if (id==(omp_get_num_threads()-1)) stop_idx=len
!!!              local_len=stop_idx-start_idx+1
!!!              do i=start_idx,stop_idx
!!!                 out(i)= mtprng_rand_real2(state_c)
!!!              enddo
!!!           end do
!!!           !$omp end parallel do
!!!        endif
!!! #else
!!!        !$omp parallel do default(shared) schedule(static,1) &
!!!        !$omp private(i,iproc,id,local_len,start_idx,stop_idx)
!!!        do iproc=1,omp_get_num_threads()
!!!           id=omp_get_thread_num()
!!!           local_len=len/omp_get_num_threads()
!!!           start_idx=id*local_len+1
!!!           stop_idx=min(len,(id+1)*local_len)
!!!           if (id==(omp_get_num_threads()-1)) stop_idx=len
!!!           local_len=stop_idx-start_idx+1
!!!           do i=start_idx,stop_idx
!!!              out(i)= mtprng_rand_real2(state_c)
!!!           enddo
!!!        end do
!!!        !$omp end parallel do
!!! #endif
      do i = 1, len
         out(i) = mtprng_rand_real2(state_c)
      end do
   end subroutine rng_uniformP

   !> Wrapper for VSL Mersenne Twister generator - Gaussian distribution
   subroutine rng_gaussian(out, len, sigma)
      integer, intent(in) :: len
      real(dblprec), dimension(len), intent(out) :: out  !< Random number in the interval N(0,sigma)
      !REAL(dblprec),dimension(len),intent(out) :: out  !< Random number in the interval [0,1)
      real(dblprec), intent(in) :: sigma
#ifdef VSL
      integer :: errcode
#endif

#ifdef VSL
      if (use_vsl) then
         errcode = vdrnggaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream(omp_get_thread_num()), len, out, 0.d0, sigma)
      else
         call fill_rngarray(out, len)
      end if
#else
      call fill_rngarray(out, len)
#endif
   end subroutine rng_gaussian

!!!   !> Wrapper for VSL Mersenne Twister generator - uniform distribution
   subroutine rng_gaussianP(out, len, sigma)
      integer, intent(in) :: len
      REAL(dblprec), dimension(len), intent(out) :: out  !< Random number in the interval [0,1)
      real(dblprec), intent(in) :: sigma
      integer :: i
!!!#ifdef VSL
!!!      integer :: errcode
!!!#endif
!!!      integer :: iproc,id,local_len,start_idx,stop_idx
!!!
!!!#ifdef VSL
!!!       if(use_vsl) then
!!!          !$omp parallel do default(shared) schedule(static,1) &
!!!          !$omp private(errcode,iproc,id,local_len,start_idx,stop_idx)
!!!          do iproc=1,omp_get_num_threads()
!!!             id=omp_get_thread_num()
!!!             local_len=len/omp_get_num_threads()
!!!             start_idx=id*local_len+1
!!!             stop_idx=min(len,(id+1)*local_len)
!!!             if (id==(omp_get_num_threads()-1)) stop_idx=len
!!!             local_len=stop_idx-start_idx+1
!!!             errcode=vdrnggaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF,stream(id),local_len, out(start_idx:stop_idx), 0.d0, sigma)
!!!          end do
!!!          !$omp end parallel do
!!!       else
!!!         call fill_rngarray(out,len)
!!!       endif
!!!#else
      call fill_rngarray(out, len)
!!!#endif
   end subroutine rng_gaussianP

   ! Wrappers for Mersenne twister generators using mtprng library

   !> Initialize the Mersenne Twister generator #1
   subroutine mt_ran_init_a(seed)
      integer, intent(in) :: seed  !< Seed number for PRNG

      call mtprng_init(seed, state_a)
   end subroutine mt_ran_init_a

   !> Initialize the Mersenne Twister generator #2
   subroutine mt_ran_init_b(seed)
      integer, intent(in) :: seed  !< Seed number for PRNG

      call mtprng_init(seed, state_b)
   end subroutine mt_ran_init_b

   !> Initialize the Mersenne Twister generator #3
   subroutine mt_ran_init_c(seed)
      integer, intent(in) :: seed  !< Seed number for PRNG
      call mtprng_init(seed, state_c)
   end subroutine mt_ran_init_c

   !> Mersenne Twister PRNG #1. Used for heat bath
   function mt_ran_a()
      REAL(dblprec) :: mt_ran_a !< Random number in the interval [0,1)

      mt_ran_a = mtprng_rand_real2(state_a)
   end function mt_ran_a

   !> Mersenne Twister PRNG #2. Used for magnetization
   function mt_ran_b()
      REAL(dblprec) :: mt_ran_b !< Random number in the interval [0,1)

      mt_ran_b = mtprng_rand_real2(state_b)
   end function mt_ran_b

   !> Mersenne Twister PRNG #3. Used for Monte Carlo , fill array
   subroutine mt_ran_carray(out, len)
      use mtprng
      implicit none
      integer, intent(in) :: len
      REAL(dblprec), dimension(len), intent(out) :: out  !< Random number in the interval [0,1)
      integer :: i
      do i = 1, len
         out(i) = mtprng_rand_real2(state_c)
      end do
   end subroutine mt_ran_carray

   !> Mersenne Twister PRNG #3. Used for Monte Carlo
   function mt_ran_c()
      REAL(dblprec) :: mt_ran_c  !< Random number in the interval [0,1)

      mt_ran_c = mtprng_rand_real2(state_c)
   end function mt_ran_c

   !> Initialize the Mersenne Twister generator #1 (used for ziggurat)
   subroutine mt_ran_init_zg(seed)
      use mtprng
      integer, intent(in) :: seed  !< Seed number for PRNG

      call mtprng_init(seed, state_z)
   end subroutine mt_ran_init_zg

   !> Ziggurat method for generating normally distributed random number
   function r4_nor()
      !*****************************************************************************80
      !
      !! R4_NOR returns a normally distributed single precision real value.
      !
      !  Discussion:
      !
      !    The value returned is generated from a distribution with mean 0 and
      !    variance 1.
      !
      !    The underlying algorithm is the ziggurat method.
      !
      !    Before the first call to this function, the user must call R4_NOR_SETUP
      !    to determine the values of KN, FN and WN.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    04 May 2008
      !
      !  Author:
      !
      !    Original C version by George Marsaglia, Wai Wan Tsang
      !    FORTRAN90 version by John Burkardt
      !
      !  Reference:
      !
      !    George Marsaglia, Wai Wan Tsang,
      !    The Ziggurat Method for Generating Random Variables,
      !    Journal of Statistical Software,
      !    Volume 5, Number 8, October 2000, seven pages.
      !
      !  Parameters:
      !
      !    Input/output, integer(INT32) JSR, the seed.
      !
      !    Input, integer(INT32) KN(128), data computed by R4_NOR_SETUP.
      !
      !    Input, real ( kind = 8 ) FN(128), WN(128), data computed by R4_NOR_SETUP.
      !
      !    Output, real ( IEEE64 ) R4_NOR, a normally distributed random value.
      !
      implicit none

      integer(INT32) :: hz
      integer(INT32) :: iz
      real(IEEE64) :: r4_nor
      real(IEEE64) :: value
      real(IEEE64) :: x
      real(IEEE64) :: y

      hz = mtprng_rand(state_z)
      iz = iand(hz, 127)

      if (abs(hz) < kn(iz + 1)) then

         value = real(hz, IEEE64)*wn(iz + 1)

      else

         do

            if (iz == 0) then

               do
                  !         x = - 0.2904764D+00 * log (mtprng_rand_real1(state))
                  !         y = - log(mtprng_rand_real1(state))
                  x = -0.2904764D+00*log(real(mtprng_rand_real1(state_z)))
                  y = -log(real(mtprng_rand_real1(state_z)))
                  if (x*x <= y + y) then
                     exit
                  end if
               end do

               if (hz <= 0) then
                  value = -r - x
               else
                  value = +r + x
               end if

               exit

            end if

            x = real(hz, IEEE64)*wn(iz + 1)

            if (fn(iz + 1) + mtprng_rand_real1(state_z)*(fn(iz) - fn(iz + 1)) < exp(-0.5E+00*x*x)) then
               value = x
               exit
            end if

            !      hz = mtprng_rand_real1( state_z )
            hz = mtprng_rand(state_z)
            iz = iand(hz, 127)

            if (abs(hz) < kn(iz + 1)) then
               value = real(hz, IEEE64)*wn(iz + 1)
               exit
            end if

         end do

      end if

      r4_nor = value

      return
   end function r4_nor

   !< Routine for setting up data needed for r4_nor() (ziggurat routine)
   subroutine r4_nor_setup(inseed)
      !*****************************************************************************80
      !
      !! R4_NOR_SETUP sets data needed by R4_NOR.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    04 May 2008
      !
      !  Author:
      !
      !    Original C version by George Marsaglia, Wai Wan Tsang
      !    FORTRAN90 version by John Burkardt
      !
      !  Reference:
      !
      !    George Marsaglia, Wai Wan Tsang,
      !    The Ziggurat Method for Generating Random Variables,
      !    Journal of Statistical Software,
      !    Volume 5, Number 8, October 2000, seven pages.
      !
      !  Parameters:
      !
      !    Output, integer(INT32) KN(128), data needed by R4_NOR.
      !
      !    Output, real ( IEEE64 ) FN(128), WN(128), data needed by R4_NOR.
      !
      implicit none

      integer, intent(in) :: inseed
      real(IEEE64) :: dn
      integer(INT32) i
      ! real(IEEE64) :: fn(128)
      ! integer(INT32) kn(128)
      ! real(IEEE64) :: wn(128)
      real(IEEE64), parameter :: m1 = 2147483648.0D+00
      real(IEEE64) :: q
      real(IEEE64) :: tn
      real(IEEE64), parameter :: vn = 9.91256303526217D-03

      ! state=instate
      call mt_ran_init_zg(inseed)

      dn = 3.442619855899D+00
      tn = 3.442619855899D+00

      q = vn/exp(-0.5D+00*dn*dn)

      kn(1) = int((dn/q)*m1)
      kn(2) = 0

      wn(1) = real(q/m1, IEEE64)
      wn(128) = real(dn/m1, IEEE64)

      fn(1) = 1.0E+00
      fn(128) = real(exp(-0.5D+00*dn*dn), IEEE64)

      do i = 127, 2, -1
         dn = sqrt(-2.0D+00*log(vn/dn + exp(-0.5D+00*dn*dn)))
         kn(i + 1) = int((dn/tn)*m1)
         tn = dn
         fn(i) = real(exp(-0.5D+00*dn*dn), IEEE64)
         wn(i) = real(dn/m1, IEEE64)
         !print *,kn(i+1),fn(i),wn(i)
      end do

      return
   end subroutine r4_nor_setup

   !> Set up transform type and initial seed for random numbers for heat bath
   subroutine setup_rng_hb(tseed, ziggurat, rngpol)
      !
      implicit none
      !
      integer, intent(in) :: tseed !< Temperature seed
      character(len=1), intent(in) :: ziggurat !< Use ziggurat transform for RNG instead of Box-Muller transform (Y/N)
      character(len=1), intent(in) :: rngpol !< Use spherical coordinates for RNG generation (Y/N)
      !
      if (ziggurat == "Y") then
         rng_trans = 1
         !$omp parallel
         call r4_nor_setup(tseed)  ! temperature
         !$omp end parallel
      else
         rng_trans = 0
         call mt_ran_init_a(tseed) ! temperature
      end if
      if (rngpol == "Y") then
         rng_sphere = 1
      else
         rng_sphere = 0
      end if
   end subroutine setup_rng_hb

   function rng_norm(mu, sigma)
      !
      implicit none
      !
      real(dblprec), intent(in) :: mu, sigma
      real(dblprec) :: rng_norm

      if (rng_trans == 1) then
         rng_norm = r4_nor()*sigma
      else
         rng_norm = gasdev(mu, sigma)
      end if
   end function rng_norm

   subroutine fill_rngarray(ranv, arr_len)
      !
      implicit none
      !
      integer, intent(in) :: arr_len
      real(dblprec), dimension(arr_len), intent(out) :: ranv
      !
      !
      integer :: i
      real(dblprec), parameter  :: one = 1.0_dblprec
      real(dblprec), parameter  :: zero = 0.0_dblprec

      if (rng_trans == 1) then
         do i = 1, arr_len
            ranv(i) = r4_nor()
         end do
      else
         do i = 1, arr_len
            ranv(i) = gasdev(zero, one)
         end do
      end if
   end subroutine fill_rngarray

   subroutine fill_rngarray_para(ranv, arr_len)
      !
      implicit none
      !
      integer, intent(in) :: arr_len
      real(dblprec), dimension(arr_len), intent(out) :: ranv
      !
      !
      integer :: i
      real(dblprec), parameter  :: one = 1.0_dblprec
      real(dblprec), parameter  :: zero = 0.0_dblprec

      if (rng_trans == 1) then
         ranv = 0.0d0
         !!!$omp parallel do default(shared) private(i) reduction(+:ranv)
         !$omp parallel do default(shared) private(i)
         do i = 1, arr_len
            ranv(i) = r4_nor()
         end do
         !$omp end parallel do
      else
         do i = 1, arr_len
            ranv(i) = gasdev(zero, one)
         end do
      end if

   end subroutine fill_rngarray_para

   subroutine draw_three_rngs(x, y, z, mu, sigma)
      use Constants, only: pi
      !
      implicit none
      !
      real(dblprec), intent(in) :: mu, sigma
      real(dblprec), intent(out) :: x, y, z
      !
      !
      real(dblprec)  :: r, phi, theta
      real(dblprec)  :: cost, sint, cosp, sinp

      if (rng_sphere /= 1) then
         x = rng_norm(mu, sigma)
         y = rng_norm(mu, sigma)
         z = rng_norm(mu, sigma)
      else
         phi = 2.0d0*pi*(mtprng_rand_real2(state_z))
         r = rng_norm(mu, sigma)**2 + rng_norm(mu, sigma)**2 + rng_norm(mu, sigma)**2
         theta = pi*(mtprng_rand_real1(state_z))
         r = sqrt(r)
         cosp = cos(phi)
         sinp = sin(phi)
         cost = cos(theta)
         sint = sin(theta)
         x = r*cosp*sint
         y = r*sinp*sint
         z = r*cost
      end if

   end subroutine draw_three_rngs

   !> Sets up an array of random numbers
   subroutine rannum(Natom, Mensemble, NA, llg, lambda1_array, lambda2_array, &
                     compensate_drift, bn, field1, field2, mmomi, Temp_array, temprescale)

      use Constants, only: k_bolt, gama, mub
!     use InputData, only : para_rng

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: llg !< Type of equation of motion (1=LLG)
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(Natom), intent(in) :: lambda2_array !< Additional damping parameter (not used for llg=1)
      integer, intent(in) :: compensate_drift !< Correct for drift in RNG
      real(dblprec), intent(in) :: bn !< Scaling factor for LLG equation (parameter=1.0d0)
      real(dblprec), dimension(3, Mensemble), intent(in) :: field1 !< Average internal effective field
      real(dblprec), dimension(3, Mensemble), intent(in) :: field2 !< Average external effective field
      real(dblprec), dimension(Natom, Mensemble), intent(in) :: mmomi !< Inverse of magnitude of magnetic moments
      real(dblprec), dimension(Natom), intent(in) :: Temp_array  !< Temperature (array)
      real(dblprec), intent(in) :: temprescale  !< Temperature rescaling from QHB

      real(dblprec) :: D
      real(dblprec) :: mu, sigma
      real(dblprec), dimension(Mensemble, Natom) :: Dk
      real(dblprec), dimension(Mensemble) :: avf1, avf2
      real(dblprec) :: rx(NA), ry(NA), rz(NA)
      logical*1 :: para_rng = .false.
      integer :: ity
      integer :: i, j

      do j = 1, Mensemble !loop over simulations, avf1,2 different for each sim.
         avf1(j) = sqrt(field1(1, j)**2 + field1(2, j)**2 + field1(3, j)**2)
         avf2(j) = sqrt(field2(1, j)**2 + field2(2, j)**2 + field2(3, j)**2)
      end do

      !   LL equations ONE universal damping
      Dk(:, :) = 0.0d0
      if (llg == 0) then
         do i = 1, Natom
            Dk(:, i) = (lambda1_array(i)*k_bolt/gama/(mub))*(gama/bn)  !last factor for dim. less.
         end do
         !   LLG equations ONE universal damping
      else if (llg == 1) then
         do i = 1, Natom
            Dk(:, i) = (lambda1_array(i)/(1 + lambda1_array(i)**2)*k_bolt/gama/(mub))*(gama/bn)  !last factor for dim. less
         end do
         !   LL equations TWO damping parameters
      else if (llg == 2) then
         do i = 1, Natom
            do j = 1, Mensemble !loop over simulations, Dk different in each sim.
               Dk(j, i) = (lambda1_array(i)*avf1(j) + lambda2_array(i)*avf2(j))/(avf1(j) + avf2(j)) &
                          *(k_bolt/gama/(mub))*(gama/bn) !last factor for dim. less.
            end do
         end do
         !   LLG equations TWO damping parameters, fluctuations included in damping1
      else if (llg == 3) then
         do i = 1, Natom
            do j = 1, Mensemble !loop over simulations, Dk different in each sim.
               Dk(j, i) = (lambda1_array(i)*avf1(j) + lambda2_array(i)*avf2(j))/(avf1(j) + avf2(j)) &
                          /(1 + lambda1_array(i)**2)*(k_bolt/gama/(mub))*(gama/bn)!last factor for dim. les
            end do
         end do
         !   LL equations TWO damping parameters, but use thermal fluctuations corresponding to one universal damping
      else if (llg == 4) then
         do i = 1, Natom
            Dk(:, i) = (lambda1_array(i)*k_bolt/gama/(mub))*(gama/bn)  !last factor for dim. less.
         end do
         !   LLG equations TWO damping parameters, but use thermal fluctuations corresponding to one universl damping
      else if (llg == 5) then
         do i = 1, Natom
            Dk(:, i) = (lambda1_array(i)/(1 + lambda1_array(i)**2)*k_bolt/gama/(mub))*(gama/bn)  !last factor for dim. less.
         end do
         !   LL equations ONE universal damping write energies
      else if (llg == 6) then
         do i = 1, Natom
            Dk(:, i) = (lambda1_array(i)*k_bolt/gama/(mub))*(gama/bn)  !last factor for dim. less.
         end do
         !   LLG equations ONE universal damping write energies
      else if (llg == 7) then
         do i = 1, Natom
            Dk(:, i) = (lambda1_array(i)/(1 + lambda1_array(i)**2)*k_bolt/gama/(mub))*(gama/bn) !last factor for dim. less.
         end do
      end if

      ! Current RNG is not parallel!
      if (.not. para_rng) then
#ifdef VSL
         if (use_vsl) then
!            call rng_gaussian(ranv,3*Natom*Mensemble,1.d0)
            call rng_gaussianP(ranv, 3*Natom*Mensemble, 1.d0)
         else
            call fill_rngarray(ranv, 3*Natom*Mensemble)
         end if
#else
         call fill_rngarray(ranv, 3*Natom*Mensemble)
#endif
      else
#ifdef VSL
         if (use_vsl) then
!            call rng_gaussian(ranv,3*Natom*Mensemble,1.d0)
            call rng_gaussianP(ranv, 3*Natom*Mensemble, 1.d0)
         else
            call fill_rngarray_para(ranv, 3*Natom*Mensemble)
         end if
#else
         call fill_rngarray_para(ranv, 3*Natom*Mensemble)
#endif
      end if

      !!!!!++ DK mod in merge, should not be a problem. Thomas(2014/07/17)
      !$omp parallel do default(shared) private(i,j,D,sigma,mu) collapse(2)
      do j = 1, Mensemble
         do i = 1, Natom
            D = Dk(j, i)*mmomi(i, j)*Temp_array(i)*temprescale
            sigma = sqrt(2.0d0*D)
            mu = 0d0
            ranv(1, i, j) = ranv(1, i, j)*sigma
            ranv(2, i, j) = ranv(2, i, j)*sigma
            ranv(3, i, j) = ranv(3, i, j)*sigma
         end do
      end do
      !$omp end parallel do

      ! Possibility to correct the random number distribution so that
      ! the true mean is zero. Not used by default
      if (compensate_drift == 1) then
         do j = 1, Mensemble
            rx = 0.0d0; ry = 0.0d0; rz = 0.0d0
            do i = 1, Natom
               ity = mod(i - 1, NA) + 1
               rx(ity) = rx(ity) + ranv(1, i, j)
               ry(ity) = ry(ity) + ranv(2, i, j)
               rz(ity) = rz(ity) + ranv(3, i, j)
            end do
            rx = rx/Natom*NA
            ry = ry/Natom*NA
            rz = rz/Natom*NA
            do i = 1, Natom
               ity = mod(i - 1, NA) + 1
               ranv(1, i, j) = ranv(1, i, j) - rx(ity)
               ranv(2, i, j) = ranv(2, i, j) - ry(ity)
               ranv(3, i, j) = ranv(3, i, j) - rz(ity)
            end do
         end do
      end if

   end subroutine rannum

   !> Sets up an array of random numbers for the lattice dynamics
   subroutine lattrannum(Natom, Mensemble, NA, lattdampvec, compensate_drift, mioninv, Temp_array, temprescale)
      !subroutine lattrannum(Natom, Mensemble, NA,  llg, lambda1_array, lambda2_array, &
      !      compensate_drift, bn, field1, field2, mmomi, Temp_array,temprescale)

      use Constants, only: k_bolt, amu
      !use Constants, only : k_bolt, gama, mub
!     use InputData, only : para_rng

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: NA  !< Number of atoms in one cell
      !integer, intent(in) :: llg !< Type of equation of motion (1=LLG)
      !real(dblprec), dimension(Natom), intent(in) :: lattdamp !< Ionic damping parameter
      real(dblprec), dimension(Natom), intent(in) :: lattdampvec !< Ionic damping parameter
      !real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      !real(dblprec), dimension(Natom), intent(in) :: lambda2_array !< Additional damping parameter (not used for llg=1)
      integer, intent(in) :: compensate_drift !< Correct for drift in RNG
      !real(dblprec), intent(in) :: bn !< Scaling factor for LLG equation (parameter=1.0d0)
      !real(dblprec), dimension(3,Mensemble), intent(in) :: field1 !< Average internal effective field
      !real(dblprec), dimension(3,Mensemble), intent(in) :: field2 !< Average external effective field
      real(dblprec), dimension(Natom, Mensemble), intent(in) :: mioninv !< Inverse of ionic mass
      real(dblprec), dimension(Natom), intent(in) :: Temp_array  !< Temperature (array)
      real(dblprec), intent(in) :: temprescale  !< Temperature rescaling from QHB

      logical*1 :: para_rng = .false.
      real(dblprec) :: D
      real(dblprec) :: mu, sigma
      real(dblprec), dimension(Mensemble, Natom) :: Dk
      !real(dblprec), dimension(Mensemble) :: avf1, avf2
      real(dblprec) :: rx(NA), ry(NA), rz(NA)
      integer :: ity
      integer :: i, j

      !do j=1,Mensemble !loop over simulations, avf1,2 different for each sim.
      !   avf1(j)=sqrt(field1(1,j)**2+field1(2,j)**2+field1(3,j)**2)
      !   avf2(j)=sqrt(field2(1,j)**2+field2(2,j)**2+field2(3,j)**2)
      !end do

      Dk(:, :) = 0.0d0
      do i = 1, Natom
         Dk(:, i) = lattdampvec(i)*k_bolt !Check units!
         !Dk(:,i)= lattdampvec(i) * k_bolt / amu  !Check units!
         !Dk(:,i)= ( lattdamp(i) * k_bolt/gama/mub ) * ( gama/bn )  !last factor for dim. less
         !Dk(:,i)= ( lambda1_array(i) / (1+lambda1_array(i)**2) * k_bolt/gama/mub ) * ( gama/bn )  !last factor for dim. less
      end do

      !   LLG equations ONE universal damping
      !if (llg==1) then
      !   do i=1, Natom
      !      Dk(:,i)=(lambda1_array(i)/(1+lambda1_array(i)**2)*k_bolt/gama/(mub))*(gama/bn)  !last factor for dim. less
      !   enddo
      !endif

      ! Current RNG is not parallel!
      if (.not. para_rng) then
#ifdef VSL
         if (use_vsl) then
!            call rng_gaussian(lattranv,3*Natom*Mensemble,1.d0)
            call rng_gaussianP(lattranv, 3*Natom*Mensemble, 1.d0)
         else
            call fill_rngarray(lattranv, 3*Natom*Mensemble)
         end if
#else
         call fill_rngarray(lattranv, 3*Natom*Mensemble)
#endif
      else
#ifdef VSL
         if (use_vsl) then
!            call rng_gaussian(lattranv,3*Natom*Mensemble,1.d0)
            call rng_gaussianP(lattranv, 3*Natom*Mensemble, 1.d0)
         else
            call fill_rngarray_para(lattranv, 3*Natom*Mensemble)
         end if
#else
         call fill_rngarray_para(lattranv, 3*Natom*Mensemble)
#endif
      end if

      !$omp parallel do default(shared) private(i,j,D,sigma,mu) collapse(2)
      do j = 1, Mensemble
         do i = 1, Natom
            D = Dk(j, i)*mioninv(i, j)*Temp_array(i)*temprescale
            !D=Dk(j,i)*mmomi(i,j)*Temp_array(i)*temprescale
            sigma = sqrt(2.0d0*D)
            mu = 0d0
            lattranv(1, i, j) = lattranv(1, i, j)*sigma
            lattranv(2, i, j) = lattranv(2, i, j)*sigma
            lattranv(3, i, j) = lattranv(3, i, j)*sigma
         end do
      end do
      !$omp end parallel do

      ! Possibility to correct the random number distribution so that
      ! the true mean is zero. Not used by default
      if (compensate_drift == 1) then
         do j = 1, Mensemble
            rx = 0.0d0; ry = 0.0d0; rz = 0.0d0
            do i = 1, Natom
               ity = mod(i - 1, NA) + 1
               rx(ity) = rx(ity) + lattranv(1, i, j)
               ry(ity) = ry(ity) + lattranv(2, i, j)
               rz(ity) = rz(ity) + lattranv(3, i, j)
            end do
            rx = rx/Natom*NA
            ry = ry/Natom*NA
            rz = rz/Natom*NA
            do i = 1, Natom
               ity = mod(i - 1, NA) + 1
               lattranv(1, i, j) = lattranv(1, i, j) - rx(ity)
               lattranv(2, i, j) = lattranv(2, i, j) - ry(ity)
               lattranv(3, i, j) = lattranv(3, i, j) - rz(ity)
            end do
         end do
      end if

   end subroutine lattrannum

end module RandomNumbers
