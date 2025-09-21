!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Xc
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
!> Module to handle exchange-correlation processes
!------------------------------------------------------------------------------

module xc_mod
   use control_mod, only: control
   use logger_mod, only: g_logger
   use string_mod, only: int2str, real2str
   use precision_mod, only: rp
   use math_mod
   use iso_c_binding, only: c_size_t
#ifdef HAVE_LIBXC
   use xc_f03_lib_m
#endif
   implicit none
   private

   type, public :: xc
      real(rp) :: AA, ACA, ALPM, AW, BB, BCA, BLPM, BW, CCA, CW, DCA, &
                  FCA, FOURPI, FTH, OCA, OTH, PCA, QCA, RCA, SCA, TCA, &
                  UCA, XALPHA, XCCF, XCCP, XCRF, XCRP, EXCHF
      character(LEN=3) :: TXCH
      integer :: NSS, txc
      integer :: LPOT
      
      ! libXC support fields
      logical :: use_libxc = .false.
      integer :: libxc_func_id = -1
      integer :: libxc_family = -1
      integer :: libxc_nspin = -1  ! Store initialization nspin for consistency
#ifdef HAVE_LIBXC
      type(xc_f03_func_info_t) :: libxc_info
#endif
   contains
      procedure :: PBEGGA
      procedure :: CORPBE
      procedure :: EXCHPBE
      procedure :: LAGGGA
      procedure :: XCPOT
      procedure :: XCPOT_hybrid
      procedure :: xcpot_libxc_wrapper
      procedure :: exchlag
      procedure :: GCOR2
      procedure :: DIFFN
      procedure :: init_libxc
      procedure :: cleanup_libxc
      procedure :: get_libxc_functional_mapping
      procedure :: validate_libxc_compatibility
      procedure :: is_libxc_functional
      final :: destructor
   end type xc

   interface xc
      procedure :: constructor
   end interface xc

contains

   !>--------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Initialize constants and text for XCPOT.
   !>
   !> Initialize constants and text for XCPOT.
   !>
   !> TODO
   !> *On entry:
   !>
   !> EXCHF : Slater exchange factor equal to 1 for full exchange
   !> IXC   :=1        Barth-Hedin
   !>         2        Slater X-Alpha
   !>         3        Barth-Hedin-Janak
   !>         4        Vosko-Wilk-Nusair
   !>         5        Perdew-Burke-Enzerhof 96 LDA
   !>         6        Wigner exchange
   !>         7        Perdew-Zunger
   !>         8        Perdew-Burke-Enzerhof 96 GGA
   !>         9        Local Airy gas
   !> NS    : Number of spins
   !>
   !> HLS:  31-Oct-96
   !>--------------------------------------------------------------------------
   function constructor(ctrl) result(obj)
      type(control), intent(in) :: ctrl
      type(xc) :: obj

      obj%txc = ctrl%txc

      ! Value taken from SETXCP call from old source atorb.f90
      obj%EXCHF = 1.0
      obj%OTH = 1.d0/3.d0
      obj%LPOT = 1
      select case (obj%txc)
      case (1)
         !
         !        Barth-Hedin J. Phys. C5, 1629(1972)
         !
         obj%TXCH = 'B-H'
         obj%XCCP = 0.0504d0
         obj%XCCF = 0.0254d0
         obj%XCRP = 30.d0
         obj%XCRF = 75.d0
         obj%FTH = 4.d0/3.d0
         obj%AA = 0.5d0**obj%OTH
         obj%BB = 1.d0 - obj%AA
      case (11)
         !
         !        Barth-Hedin J. Phys. C5, 1629(1972)
         !
         obj%TXCH = 'ASW'
         obj%XCCP = 0.0450d0
         obj%XCCF = 0.0225d0
         obj%XCRP = 21.d0
         obj%XCRF = 52.9167d0
         obj%FTH = 4.d0/3.d0
         obj%AA = 0.5d0**obj%OTH
         obj%BB = 1.d0 - obj%AA
      case (2)
         !
         !        Slater X-Alpha
         !
         obj%TXCH = 'X-A'
         obj%XALPHA = 6.d0*obj%EXCHF*(3.d0/(4.d0*PI))**obj%OTH
         call g_logger%debug(' SETXCP: '//obj%TXCH//' Slater exchange alpha = '//real2str(obj%XALPHA), __FILE__, __LINE__)
      case (3)
         !
         !        Barth-Hedin-Janak Phys. Rev. B12, 1257(1975)
         !
         obj%TXCH = 'BHJ'
         obj%XCCP = 0.045d0
         obj%XCCF = 0.0225d0
         obj%XCRP = 21.d0
         obj%XCRF = 53.d0
         obj%FTH = 4.d0/3.d0
         obj%AA = 0.5d0**obj%OTH
         obj%BB = 1.d0 - obj%AA
      case (4)
         !
         !        Vosko-Wilk-Nusair Can. J. Phys. 58, 1200(1980)
         !
         obj%TXCH = 'VWN'
         obj%FTH = 4.d0/3.d0
         obj%AA = 2.d0**obj%FTH - 2.d0
      case (5)
         !
         !        Perdew-Burke-Enzerhof 96 LDA
         !
         obj%TXCH = 'PBE'
      case (6)
         !
         !        Wigner exchange
         !
         obj%TXCH = 'WXC'
         if (ctrl%nsp == 2) call g_logger%fatal(' SETXCP:** Spin polarization not implemented for IXC = '//int2str(obj%txc)//' Xcpot = '//obj%TXCH, __FILE__, __LINE__)
         obj%AW = 0.916d0*4.d0/3.d0
         obj%BW = 0.88d0*4.d0/3.d0
         obj%CW = 0.88d0*7.8d0/3.d0
      case (7)
         !
         !        Ceperley-Alder. Parametrization by Perdew and Zunger
         !
         obj%TXCH = 'P-Z'
         if (ctrl%nsp == 2) call g_logger%fatal(' SETXCP:** Spin polarization not implemented for IXC = '//int2str(obj%txc)//' Xcpot = '//obj%TXCH, __FILE__, __LINE__)
         obj%ACA = 1.0529d0
         obj%BCA = 0.3334d0
         obj%CCA = 7.d0*obj%ACA/6.d0
         obj%DCA = 4.d0*obj%BCA/3.d0
         obj%FCA = 4.d0/3.d0
         obj%OCA = 0.096d0
         obj%PCA = 0.0622d0
         obj%QCA = 0.0232d0
         obj%RCA = 0.004d0
         obj%SCA = obj%OCA + obj%PCA/3.d0
         obj%TCA = (2.d0*obj%QCA + obj%RCA)/3.d0
         obj%UCA = 2.d0*obj%RCA/3.d0
      case (8)
         obj%TXCH = 'GGA'
      case (9)
         !
         !        Linear Airy gas
         !
         obj%TXCH = 'LAG'
      case default
         ! Check if this is a libXC functional (txc >= 1000)
         if (obj%txc >= 1000) then
            ! libXC functional - will be handled in init_libxc
            obj%TXCH = 'LXC'  ! Mark as libXC functional
         else
            ! Unknown legacy functional
            if (ctrl%nsp == 2) call g_logger%fatal(' SETXCP:** IXC = '//int2str(obj%txc)//' not implemented', __FILE__, __LINE__)
         endif
      end select
      
      ! Initialize libXC support if needed
      call obj%init_libxc(ctrl)
   end function constructor

   !==============================================================================
   !----- FROM HERE WE HAVE SUBROUTINES AND FUCTIONS RAW COPIED FROM OLD CODE ----
   !==============================================================================

   subroutine XCPOT(this, RHO1, RHO2, RHO, RHOP, RHOPP, RR, V1, V2, EXC)
      !   ******************************************************************
      !   *                                                                *
      !   *    Calculates exchange-correlation potential according to the  *
      !   *    value of IXC. The constants used in the various expressions *
      !   *    are set in SETXCP and in the data statement below.          *
      !   *                                                                *
      !   *    HLS: 31-Oct-96                                              *
      !   ******************************************************************
      ! use xcdata
      ! use mpi_global_data
      !
      !.. Implicit Declarations ..
      ! implicit none
      !
      !.. Parameters ..
      class(xc), intent(in) :: this
      real(rp), parameter :: TOLD = 1.d-20
      real(rp), parameter :: TOLDD = 1.d-20
      !
      !.. Formal Arguments ..
      real(rp), intent(in) :: RHO, RHO1, RHO2, RR
      real(rp), intent(inout) :: EXC, V1, V2
      real(rp), dimension(2), intent(in) :: RHOP, RHOPP
      !
      !.. Local Scalars ..
      integer :: IXC
      integer :: LGGA
      real(rp) :: AP = 0.0621814
      real(rp) :: AMYX1, AMYX2, ARS, ATNF, ATNP, BETA, BRS, CNY, &
                  DBETA, DENOM1, DFS, DUC, DUC1, DUC2, EC, ECF, ECP, &
                  EPSCF, EPSCP, EPSXP, EPX, EX, FCF, FCP, FS, FX, RS78, &
                  RSF, RSF2, RSF3, RSLN, RSLOG, RSP, RSP2, RSP3, S, S4, &
                  SM, SP, SQRTRS
      real(rp) :: TF1, TP1, TRX1, TRX2, UC0, UC1, UC10, UC2, UC20, UCF, &
                  UCP, X, XFX, XPX
      real(rp) :: AF = 0.0310907, BF = 7.060428, BP = 3.72744, &
                  CF = 18.0578, CF1 = 2.9847935, &
                  CF2 = 2.7100059, CF3 = -0.1446006, &
                  CP = 12.9352, CP1 = 1.2117833, &
                  CP2 = 1.1435257, CP3 = -0.031167608, &
                  QF = 4.7309269, QP = 6.1519908, &
                  XF0 = -0.32500, XP0 = -0.10498
      real(rp) :: RS, RS1
      !
      !.. Local Arrays ..
      real(rp), dimension(2) :: RHOS
      !
      !.. External Calls ..
      ! external LAGGGA, PBEGGA
      !
      !.. Intrinsic Functions ..
      intrinsic ATAN, LOG, SQRT
      !
      ! ... Executable Statements ...
      !
      !
      if (RHO1 < TOLD .or. RHO2 < TOLDD) then
         V1 = 0.d0
         V2 = 0.d0
         EXC = 0.d0
         return
      end if
      RS1 = ((4*pi)*RHO/3.)**this%OTH
      RS = 1./RS1
      IXC = this%TXC
      select case (IXC)
      case (2)
         !
         !     Slater X-Alpha
         !
         EXC = -0.75d0*this%XALPHA*(0.5d0*RHO)**this%OTH
         V1 = -this%XALPHA*(RHO1)**this%OTH
         V2 = -this%XALPHA*(RHO2)**this%OTH
      case (4)
         !
         !     Vosko-Wilk-Nusair  Can. J. Phys. 58, 1200(1980)
         !
         X = SQRT(RS)
         XPX = X*X + BP*X + CP
         XFX = X*X + BF*X + CF
         S = (RHO2 - RHO1)/RHO
         SP = 1.d0 + S
         SM = 1.d0 - S
         S4 = S**4 - 1.d0
         FS = (SP**this%FTH + SM**this%FTH - 2.d0)/this%AA
         BETA = 1.d0/(2.74208d0 + 3.182d0*X + 0.09873d0*X*X + 0.18268d0*X**3)
         DFS = this%FTH*(SP**this%OTH - SM**this%OTH)/this%AA
         DBETA = -(0.27402d0*X + 0.09873d0 + 1.591d0/X)*BETA**2
         ATNP = ATAN(QP/(2.d0*X + BP))
         ATNF = ATAN(QF/(2.d0*X + BF))
         ECP = AP*(LOG(X*X/XPX) + CP1*ATNP - CP3*(LOG((X - XP0)**2/XPX) + CP2*ATNP))
         ECF = AF*(LOG(X*X/XFX) + CF1*ATNF - CF3*(LOG((X - XF0)**2/XFX) + CF2*ATNF))
         EC = ECP + FS*(ECF - ECP)*(1.d0 + S4*BETA)
         TP1 = (X*X + BP*X)/XPX
         TF1 = (X*X + BF*X)/XFX
         UCP = ECP - AP/3.d0*(1.d0 - TP1 - CP3*(X/(X - XP0) - TP1 - XP0*X/XPX))
         UCF = ECF - AF/3.d0*(1.d0 - TF1 - CF3*(X/(X - XF0) - TF1 - XF0*X/XFX))
         UC0 = UCP + (UCF - UCP)*FS
         UC20 = UC0 + (ECF - ECP)*SM*DFS
         UC10 = UC0 - (ECF - ECP)*SP*DFS
         DUC = (UCF - UCP)*BETA*S4*FS + (ECF - ECP)*(-RS/3.d0)*DBETA*S4*FS
         DUC2 = DUC + (ECF - ECP)*BETA*SM*(4.d0*S**3*FS + S4*DFS)
         DUC1 = DUC - (ECF - ECP)*BETA*SP*(4.d0*S**3*FS + S4*DFS)
         UC1 = UC10 + DUC1
         UC2 = UC20 + DUC2
         EPX = -0.91633059d0/RS*(1.d0 + this%FTH*FS/5.1297628d0)
         AMYX2 = -1.22177412d0/RS*SP**this%OTH
         AMYX1 = -1.22177412d0/RS*SM**this%OTH
         EXC = EC + EPX
         V1 = UC1 + AMYX1
         V2 = UC2 + AMYX2
      case (5)
         !
         !     Perdew, Burke, and Ernzerhof, submiited to PRL, May96
         !
         !     LDA version
         !
         if (RR < 1.d-10) then
!  #ifdef MPI
!         if (myid==0) then
!  #endif
            write (6, 10000) IXC
!  #ifdef MPI
!         end if
!  #endif
            stop
         else
            LGGA = 0
            RHOS(1) = RHO1
            RHOS(2) = RHO2
            call this%PBEGGA(this%LPOT, LGGA, RHOS, RHOP, RHOPP, RR, EXC, V1, V2)
         end if
      case (6)
         !
         !     Wigner expression
         !
         RS78 = 1.d0/(RS + 7.8d0)
         EXC = -0.916d0*RS1 - 0.88d0*RS78
         !
         V1 = this%CW*RS78*RS78 - this%AW*RS1 - this%BW*RS78
         V2 = V1
      case (7)
         !
         !     Ceperley-Alder. Parametrization by Perdew and Zunger.
         !
         if (RS >= 1.d0) then
            SQRTRS = SQRT(RS)
            DENOM1 = 1.d0/(1.d0 + this%ACA*SQRTRS + this%BCA*RS)
            EX = -0.9164d0*RS1
            EC = -0.2846d0*DENOM1
            EXC = EX + EC
            V1 = this%FCA*EX + EC*(1.d0 + this%CCA*SQRTRS + this%DCA*RS)*DENOM1
            V2 = V1
         else
            RSLOG = LOG(RS)
            RSLN = RS*RSLOG
            EX = -0.9164d0*RS1
            EC = -this%OCA + this%PCA*RSLOG - this%QCA*RS + this%RCA*RSLN
            EXC = EX + EC
            V1 = this%FCA*EX - this%SCA + this%PCA*RSLOG - this%TCA*RS + this%UCA*RSLN
            V2 = V1
         end if
      case (8)
         !
         !     Perdew, Burke, and Ernzerhof, submiited to PRL, May96
         !
         !     GGA version
         !
         if (RR < 1.d-10) then
!  #ifdef MPI
!         if (myid==0) then
!  #endif
            write (6, 10000) IXC
!  #ifdef MPI
!         end if
!  #endif
            stop
         else
            LGGA = 1
            RHOS(1) = RHO1
            RHOS(2) = RHO2
            call this%PBEGGA(this%LPOT, LGGA, RHOS, RHOP, RHOPP, RR, EXC, V1, V2)
         end if
      case (9)
         !
         !     Local Airy Gas plus PBE correlation
         !
         if (RR < 1.d-10) then
! #ifdef MPI
!         if (myid==0) then
! #endif
            write (6, 10000) IXC
! #ifdef MPI
!         end if
! #endif
            stop
         else
            LGGA = 0
            RHOS(1) = RHO1
            RHOS(2) = RHO2
            call this%LAGGGA(this%LPOT, LGGA, RHOS, RHOP, RHOPP, RR, EXC, V1, V2)
         end if
      case default
         !
         !     Barth-Hedin  J. PHYS. C5, 1629(1972)
         !
         RSF = RS/this%XCRF
         RSF2 = RSF*RSF
         RSF3 = RSF2*RSF
         RSP = RS/this%XCRP
         RSP2 = RSP*RSP
         RSP3 = RSP2*RSP
         FCF = (1.d0 + RSF3)*LOG(1.d0 + 1.d0/RSF) + .5d0*RSF - RSF2 - this%OTH
         FCP = (1.d0 + RSP3)*LOG(1.d0 + 1.d0/RSP) + .5d0*RSP - RSP2 - this%OTH
         EPSCP = -this%XCCP*FCP
         EPSCF = -this%XCCF*FCF
         EPSXP = -.91633059d0/RS
         CNY = 5.1297628d0*(EPSCF - EPSCP)
         X = RHO1/RHO
         FX = (X**this%FTH + (1.d0 - X)**this%FTH - this%AA)/this%BB
         EXC = EPSXP + EPSCP + FX*(CNY + this%FTH*EPSXP)/5.1297628d0
         ARS = -1.22177412d0/RS + CNY
         BRS = -this%XCCP*LOG(1.d0 + this%XCRP/RS) - CNY
         TRX1 = (2.d0*X)**this%OTH
         V1 = ARS*TRX1 + BRS
         TRX2 = (2.d0*RHO2/RHO)**this%OTH
         V2 = ARS*TRX2 + BRS
      end select
      return
      !
      ! ... Format Declarations ...
      !
10000 format(/, " XCPOT:**  RR less than 1. D-10 for IXC =", i3)
   end subroutine XCPOT

   subroutine PBEGGA(this, LPOTT, LGGA, N, ND, NDD, R, EXC, MUXC1, MUXC2)
      !   ******************************************************************
      !   *                                                                *
      !   *    Calculate the exchange-correlation energy density by        *
      !   *    Perdew-Burke-Ernzerhof generalized gradient approximation.  *
      !   *                                                                *
      !   *   *On entry:                                                   *
      !   *       n(r)   = the charge density per spin (N) ;               *
      !   *       n´(r)  = dn/dr (ND) ;                                    *´
      !   *       n´´(r) = d^2n/dr^2 (NDD) .                               *
      !   *       LGGA = 0 LDA (local density approximation) ;             *
      !   *            = 1 GGA (generalized gradient approximation) .      *
      !   *   *On exit:                                                    *
      !   *       exc    = the exchange-correlation energy density,        *
      !   *       muxc   = the exchange-correlation potential.             *
      !   *                                                                *
      !   ******************************************************************
      ! use xcdata
      ! use precision_mod, only: rp
      !
      !.. Implicit Declarations ..
      implicit none
      !
      !.. Formal Arguments ..
      class(xc), intent(in) :: this
      integer, intent(in) :: LGGA, LPOTT
      real(rp), intent(in) :: R
      real(rp), intent(out) :: EXC, MUXC1, MUXC2
      real(rp), dimension(2), intent(in) :: N, ND, NDD
      !
      !.. Local Scalars ..
      integer :: I
      real(rp) :: DVCDN, DVCUP, EC, EX, EXI, FK, G, H, KF, MUXI, NABLA, &
                  NABLA2, NDDI, NDI, NI, RS, S, SK, T, U, UU, VCDN, &
                  VCUP, VV, WW, ZET
      !
      !.. Local Arrays ..
      real(rp), dimension(2) :: VX
      !
      !.. External Calls ..
      ! external CORPBE, EXCHPBE
      !
      !.. Intrinsic Functions ..
      ! intrinsic ABS, SQRT
      !
      ! ... Executable Statements ...
      !
      !
      !PI = 0.25d0 * (4*pi)
      !
      !     Exchange energy and potential
      !
      EX = 0.d0
      do I = 1, 2
         NI = N(I) + N(I)
         EXI = 0.d0
         MUXI = 0.d0
         ! Fixed: Don't double the gradient inputs
         NDI = ND(I)
         NDDI = NDD(I)
         KF = (3.d0*PI*PI*NI)**this%OTH
         NABLA = ABS(NDI)
         S = 0.5d0*NABLA/KF/NI
         ! Fixed: Correct Laplacian in spherical coordinates
         NABLA2 = NDDI + 2.d0/R*NDI
         T = NABLA2/4.d0/KF/KF/NI
         U = NABLA*NDDI/8.d0/KF/KF/KF/NI/NI
         call this%EXCHPBE(NI, S, U, T, LGGA, LPOTT, EXI, MUXI)
         VX(I) = MUXI
         EX = EX + N(I)*EXI
      end do
      !
      !     Correlation energy and potential
      !
      NI = N(1) + N(2)
      EC = 0.d0
      VCUP = 0.d0
      DVCUP = 0.d0
      H = 0.d0
      VCDN = 0.d0
      DVCDN = 0.d0
      NDI = ND(1) + ND(2)
      NDDI = NDD(1) + NDD(2)
      ZET = (N(1) - N(2))/NI
      G = ((1.d0 + ZET)**(2.d0/3.d0) + (1.d0 - ZET)**(2.d0/3.d0))/2.d0
      NABLA = ABS(NDI)
      NABLA2 = NDDI + 2.d0/R*NDI  ! Fixed: Correct Laplacian
      FK = (3.d0*PI*PI*NI)**this%OTH
      SK = SQRT(4.d0*FK/PI)
      T = NABLA/2.d0/SK/NI/G
      UU = NABLA*NDDI/((2.d0*SK*G)**3.d0)/NI/NI
      VV = NABLA2/((2.d0*SK*G)**2.d0)/NI
      WW = (NDI*ND(1) - NDI*ND(2) - ZET*NDI*NDI)/((2.d0*SK*G)**2.d0)/NI/NI
      RS = (3.d0/(4*pi)/NI)**this%OTH
      !
      call this%CORPBE(RS, ZET, T, UU, VV, WW, LGGA, LPOTT, EC, VCUP, VCDN, H, DVCUP, DVCDN)
      !
      !     Convert to Rydberg
      !
      MUXC1 = 2.d0*(VX(1) + VCUP + DVCUP)
      MUXC2 = 2.d0*(VX(2) + VCDN + DVCDN)
      EX = 2.d0*EX/(N(1) + N(2))
      EC = 2.d0*(EC + H)
      !
      !
      EXC = EX + EC
   end subroutine PBEGGA
   subroutine EXCHPBE(this, rho, S, U, V, lgga, lpot, EX, VX)
      !----------------------------------------------------------------------
      !  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
      !  K Burkes modification of PW91 codes, May 14, 1996
      !  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
      !----------------------------------------------------------------------
      !----------------------------------------------------------------------
      !  INPUT rho : DENSITY
      !  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
      !  INPUT U:  (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KF)**3)
      !  INPUT V: (LAPLACIAN rho)/(rho*(2*KF)**2)
      !   (for U, V, see PW86(24))
      !  input lgga:  (=0=>dont put in gradient corrections, just LDA)
      !  input lpot:  (=0=>don´t get potential and don´t need U and V)
      !  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
      !----------------------------------------------------------------------
      !----------------------------------------------------------------------
      ! References:
      ! [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submiited to PRL, May96
      ! [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
      !     {\bf 40},  3399  (1989) (E).
      !----------------------------------------------------------------------
      !----------------------------------------------------------------------
      ! Formulas:
      !           e_x[unif]=ax*rho^(4/3)  [LDA]
      ! ax = -0.75*(3/pi)^(1/3)
      !        e_x[PBE]=e_x[unif]*FxPBE(s)
      !        FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
      ! uk, ul defined after [a](13)
      !----------------------------------------------------------------------
      !----------------------------------------------------------------------
      ! use precision_mod, only: rp
      !
      !.. Implicit Declarations ..
      implicit none
      !
      !.. Parameters ..
      class(xc), intent(in) :: this
      real(rp), parameter :: thrd = 1.d0/3.d0
      real(rp), parameter :: thrd4 = 4.d0/3.d0
      real(rp), parameter :: &
         ax = -0.738558766382022405884230032680836d0
      real(rp), parameter :: um = 0.2195149727645171d0
      real(rp), parameter :: uk = 0.8040d0
      real(rp), parameter :: ul = um/uk
      !
      !.. Formal Arguments ..
      integer, intent(in) :: lgga, lpot
      real(rp), intent(in) :: S, U, V, rho
      real(rp), intent(out) :: EX, VX
      !
      !.. Local Scalars ..
      real(rp) :: Fs, Fss, FxPBE, P0, S2, exunif
      !
      ! ... Executable Statements ...
      !
      !----------------------------------------------------------------------
      !----------------------------------------------------------------------
      ! construct LDA exchange energy density
      exunif = AX*rho**THRD
      if (lgga == 0) then
         ex = exunif
         vx = ex*thrd4
      else
         !----------------------------------------------------------------------
         !----------------------------------------------------------------------
         ! construct PBE enhancement factor
         S2 = S*S
         P0 = 1.d0 + ul*S2
         FxPBE = 1d0 + uk - uk/P0
         EX = exunif*FxPBE
         vx = exunif*thrd4
         if (lpot /= 0) then
            !----------------------------------------------------------------------
            !----------------------------------------------------------------------
            !  ENERGY DONE. NOW THE POTENTIAL:
            !  find first and second derivatives of Fx w.r.t s.
            !  Fs=(1/s)*d FxPBE/ ds
            !  Fss=d Fs/ds
            Fs = 2.d0*uk*ul/(P0*P0)
            Fss = -4.d0*ul*S*Fs/P0
            !----------------------------------------------------------------------
            !----------------------------------------------------------------------
            ! calculate potential from [b](24)
            VX = exunif*(THRD4*FxPBE - (U - THRD4*S2*s)*FSS - V*FS)
         end if
      end if
   end subroutine EXCHPBE

   subroutine CORPBE(this, RS, ZET, T, UU, VV, WW, lgga, lpot, ec, vcup, vcdn, H, DVCUP, DVCDN)
      !----------------------------------------------------------------------
      !  Official PBE correlation code. K. Burke, May 14, 1996.
      !  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
      !       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
      !       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
      !       : UU=(GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KS*G)**3)
      !       : VV=(LAPLACIAN rho)/(rho * (2*KS*G)**2)
      !       : WW=(GRAD rho)*(GRAD ZET)/(rho * (2*KS*G)**2
      !       :  UU, VV, WW, only needed for PBE potential
      !       : lgga=flag to do gga (0=>LSD only)
      !       : lpot=flag to do potential (0=>energy only)
      !  output: ec=lsd correlation energy from [a]
      !        : vcup=lsd up correlation potential
      !        : vcdn=lsd dn correlation potential
      !        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
      !        : dvcup=nonlocal correction to vcup
      !        : dvcdn=nonlocal correction to vcdn
      !----------------------------------------------------------------------
      !----------------------------------------------------------------------
      ! References:
      ! [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof,
      !     {\sl Generalized gradient approximation made simple}, sub.
      !     to Phys. Rev.Lett. May 1996.
      ! [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl Real-space cutoff
      !     construction of a generalized gradient approximation:  The PW91
      !     density functional}, submitted to Phys. Rev. B, Feb. 1996.
      ! [c] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
      !----------------------------------------------------------------------
      !----------------------------------------------------------------------
      ! use precision_mod, only: rp
      !
      !.. Implicit Declarations ..
      implicit none
      !
      !.. Parameters ..
      class(xc), intent(in) :: this
      real(rp), parameter :: thrd = 1.d0/3.d0
      real(rp), parameter :: thrdm = -thrd
      real(rp), parameter :: thrd2 = 2.d0*thrd
      real(rp), parameter :: sixthm = thrdm/2.d0
      real(rp), parameter :: thrd4 = 4.d0*thrd
      real(rp), parameter :: &
         GAM = 0.5198420997897463295344212145565d0
      real(rp), parameter :: fzz = 8.d0/(9.d0*GAM)
      real(rp), parameter :: &
         gamma = 0.03109069086965489503494086371273d0
      real(rp), parameter :: bet = 0.06672455060314922d0
      real(rp), parameter :: delt = bet/gamma
      real(rp), parameter :: eta = 1.d-12
      !
      !.. Formal Arguments ..
      integer, intent(in) :: lgga, lpot
      real(rp), intent(in) :: RS, T, UU, VV, WW, ZET
      real(rp), intent(inout) :: DVCDN, DVCUP, H
      real(rp), intent(out) :: ec, vcdn, vcup
      !
      !.. Local Scalars ..
      real(rp) :: ALFM, ALFRSM, B, B2, BEC, BG, COMM, ECRS, ECZET, EP, &
                  EPRS, EU, EURS, F, FAC, FACT0, FACT1, FACT2, FACT3, &
                  FACT5, FZ, G, G3, G4, GZ, PON, PREF, Q4, Q5, Q8, Q9, &
                  RSTHRD, T2, T4, T6, Z4, hB, hBT, hRS, hRST
      real(rp) :: hT, hTT, hZ, hZT, rtrs
      !
      !.. External Calls ..
      ! external gcor2
      !
      !.. Intrinsic Functions ..
      intrinsic EXP, LOG, SQRT
      !
      ! ... Executable Statements ...
      !
      !----------------------------------------------------------------------
      !----------------------------------------------------------------------
      ! find LSD energy contributions, using [c](10) and Table I[c].
      ! EU=unpolarized LSD correlation energy
      ! EURS=dEU/drs
      ! EP=fully polarized LSD correlation energy
      ! EPRS=dEP/drs
      ! ALFM=-spin stiffness, [c](3).
      ! ALFRSM=-dalpha/drs
      ! F=spin-scaling factor from [c](9).
      ! construct ec, using [c](8)
      rtrs = SQRT(rs)
      call this%gcor2(0.0310907d0, 0.21370d0, 7.5957d0, 3.5876d0, 1.6382d0, 0.49294d0, rtrs, &
                      EU, EURS)
      call this%gcor2(0.01554535d0, 0.20548d0, 14.1189d0, 6.1977d0, 3.3662d0, 0.62517d0, &
                      rtRS, EP, EPRS)
      call this%gcor2(0.0168869d0, 0.11125d0, 10.357d0, 3.6231d0, 0.88026d0, 0.49671d0, rtRS, &
                      ALFM, ALFRSM)
      Z4 = ZET**4
      F = ((1.d0 + ZET)**THRD4 + (1.d0 - ZET)**THRD4 - 2.d0)/GAM
      EC = EU*(1.d0 - F*Z4) + EP*F*Z4 - ALFM*F*(1.d0 - Z4)/FZZ
      !----------------------------------------------------------------------
      !----------------------------------------------------------------------
      ! LSD potential from [c](A1)
      ! ECRS = dEc/drs [c](A2)
      ! ECZET=dEc/dzeta [c](A3)
      ! FZ = dF/dzeta [c](A4)
      ECRS = EURS*(1.d0 - F*Z4) + EPRS*F*Z4 - ALFRSM*F*(1.d0 - Z4)/FZZ
      FZ = THRD4*((1.d0 + ZET)**THRD - (1.d0 - ZET)**THRD)/GAM
      ECZET = 4.d0*(ZET**3)*F*(EP - EU + ALFM/FZZ) + &
              FZ*(Z4*EP - Z4*EU - (1.d0 - Z4)*ALFM/FZZ)
      COMM = EC - RS*ECRS/3.d0 - ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET
      if (lgga /= 0) then
         !----------------------------------------------------------------------
         !----------------------------------------------------------------------
         ! PBE correlation energy
         ! G=phi(zeta), given after [a](3)
         ! DELT=bet/gamma
         ! B=A of [a](8)
         G = ((1.d0 + ZET)**thrd2 + (1.d0 - ZET)**thrd2)/2.d0
         G3 = G**3
         PON = -EC/(G3*gamma)
         B = DELT/(EXP(PON) - 1.d0)
         B2 = B*B
         T2 = T*T
         T4 = T2*T2
         Q4 = 1.d0 + B*T2
         Q5 = 1.d0 + B*T2 + B2*T4
         H = G3*(BET/DELT)*LOG(1.d0 + DELT*Q4*T2/Q5)
         if (lpot /= 0) then
            !----------------------------------------------------------------------
            !----------------------------------------------------------------------
            ! ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
            G4 = G3*G
            T6 = T4*T2
            RSTHRD = RS/3.d0
            GZ = (((1.d0 + zet)**2 + eta)**sixthm - ((1.d0 - zet)**2 + eta)**sixthm)/3.d0
            FAC = DELT/B + 1.d0
            BG = -3.d0*B2*EC*FAC/(BET*G4)
            BEC = B2*FAC/(BET*G3)
            Q8 = Q5*Q5 + DELT*Q4*Q5*T2
            Q9 = 1.d0 + 2.d0*B*T2
            hB = -BET*G3*B*T6*(2.d0 + B*T2)/Q8
            hRS = -RSTHRD*hB*BEC*ECRS
            FACT0 = 2.d0*DELT - 6.d0*B
            FACT1 = Q5*Q9 + Q4*Q9*Q9
            hBT = 2.d0*BET*G3*T4*((Q4*Q5*FACT0 - DELT*FACT1)/Q8)/Q8
            hRST = RSTHRD*T2*hBT*BEC*ECRS
            hZ = 3.d0*GZ*h/G + hB*(BG*GZ + BEC*ECZET)
            hT = 2.d0*BET*G3*Q9/Q8
            hZT = 3.d0*GZ*hT/G + hBT*(BG*GZ + BEC*ECZET)
            FACT2 = Q4*Q5 + B*T2*(Q4*Q9 + Q5)
            FACT3 = 2.d0*B*Q5*Q9 + DELT*FACT2
            hTT = 4.d0*BET*G3*T*(2.d0*B/Q8 - (Q9*FACT3/Q8)/Q8)
            COMM = H + HRS + HRST + T2*HT/6.d0 + 7.d0*T2*T*HTT/6.d0
            PREF = HZ - GZ*T2*HT/G
            FACT5 = GZ*(2.d0*HT + T*HTT)/G
            COMM = COMM - PREF*ZET - UU*HTT - VV*HT - WW*(HZT - FACT5)
            DVCUP = COMM + PREF
            DVCDN = COMM - PREF
         end if
      end if
   end subroutine CORPBE

   subroutine LAGGGA(this, LPOTT, LGGA, N, ND, NDD, R, EXC, MUXC1, MUXC2)
      !   ******************************************************************
      !   *                                                                *
      !   *    Calculate the exchange-correlation energy density by        *
      !   *    Local Airy Gas approximation plus Perdew-Burke-Ernzerhof    *
      !   *    correlation functional.                                     *
      !   *                                                                *
      !   *   *On entry:                                                   *
      !   *       n(r)   = the charge density per spin (N) ;               *
      !   *       n´(r)  = dn/dr (ND) ;                                    *´
      !   *       n´´(r) = d^2n/dr^2 (NDD) .                               *
      !   *       LGGA = 0 LDA (local density approximation) ;             *
      !   *            = 1 GGA (generalized gradient approximation) .      *
      !   *   *On exit:                                                    *
      !   *       exc    = the exchange-correlation energy density,        *
      !   *       muxc   = the exchange-correlation potential.             *
      !   *                                                                *
      !   ******************************************************************
      ! use xcdata
      ! use precision_mod, only: rp
      !
      !.. Implicit Declarations ..
      implicit none
      !
      !.. Formal Arguments ..
      class(xc), intent(in) :: this
      integer, intent(in) :: LGGA, LPOTT
      real(rp), intent(in) :: R
      real(rp), intent(out) :: EXC, MUXC1, MUXC2
      real(rp), dimension(2), intent(in) :: N, ND, NDD
      !
      !.. Local Scalars ..
      integer :: I
      real(rp) :: DVCDN, DVCUP, EC, EX, EXI, FK, G, H, KF, MUXI, NABLA, &
                  NABLA2, NDDI, NDI, NI, PI, RS, S, SK, T, U, UU, VCDN, &
                  VCUP, VV, WW, ZET
      !
      !.. Local Arrays ..
      real(rp), dimension(2) :: VX
      !
      !.. External Calls ..
      external CORPBE, exchlag
      !
      !.. Intrinsic Functions ..
      intrinsic ABS, SQRT
      !
      ! ... Executable Statements ...
      !
      !
      PI = 0.25d0*(4*pi)
      !
      !     Exchange energy and potential
      !
      EX = 0.d0
      do I = 1, 2
         NI = N(I) + N(I)
         EXI = 0.d0
         MUXI = 0.d0
         NDI = ND(I) + ND(I)
         if (ABS(NDI) < 1d-15) then
            NDI = 1d-15
         end if
         NDDI = NDD(I) + NDD(I)
         KF = (3.d0*PI*PI*NI)**this%OTH
         NABLA = ABS(NDI)
         S = 0.5d0*NABLA/KF/NI
         NABLA2 = NDDI + 2.d0/R*NDI  ! Fixed: Correct Laplacian
         T = NABLA2/4.d0/KF/KF/NI
         U = NABLA*NDDI/8.d0/KF/KF/KF/NI/NI
         call this%exchlag(NI, S, U, T, LGGA, LPOTT, EXI, MUXI)
         VX(I) = MUXI
         EX = EX + N(I)*EXI
      end do
      !
      !     Correlation energy and potential
      !
      NI = N(1) + N(2)
      EC = 0.d0
      VCUP = 0.d0
      DVCUP = 0.d0
      H = 0.d0
      VCDN = 0.d0
      DVCDN = 0.d0
      NDI = ND(1) + ND(2)
      NDDI = NDD(1) + NDD(2)
      ZET = (N(1) - N(2))/NI
      G = ((1.d0 + ZET)**(2.d0/3.d0) + (1.d0 - ZET)**(2.d0/3.d0))/2.d0
      NABLA = ABS(NDI)
      NABLA2 = NDDI + 2.d0/R*NDI  ! Fixed: Correct Laplacian
      FK = (3.d0*PI*PI*NI)**this%OTH
      SK = SQRT(4.d0*FK/PI)
      T = NABLA/2.d0/SK/NI/G
      UU = NABLA*NDDI/((2.d0*SK*G)**3.d0)/NI/NI
      VV = NABLA2/((2.d0*SK*G)**2.d0)/NI
      WW = (NDI*ND(1) - NDI*ND(2) - ZET*NDI*NDI)/((2.d0*SK*G)**2.d0)/NI/NI
      RS = (3.d0/(4*pi)/NI)**this%OTH
      !
      call this%CORPBE(RS, ZET, T, UU, VV, WW, LGGA, LPOTT, EC, VCUP, VCDN, H, DVCUP, DVCDN)
      !
      !     Convert to Rydberg
      !
      MUXC1 = 2.d0*(VX(1) + VCUP + DVCUP)
      MUXC2 = 2.d0*(VX(2) + VCDN + DVCDN)
      EX = 2.d0*EX/(N(1) + N(2))
      EC = 2.d0*(EC + H)
      !
      !
      EXC = EX + EC
   end subroutine LAGGGA

   subroutine exchlag(this, rho, s, u, v, lgga, lpot, ex, vx)
      !----------------------------------------------------------------------
      !  Local Airy Gas approximation
      !----------------------------------------------------------------------
      !  INPUT rho : DENSITY
      !  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
      !  INPUT U:  (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KF)**3)
      !  INPUT V: (LAPLACIAN rho)/(rho*(2*KF)**2)
      !   (for U, V, see PW86(24))
      !  input lgga:  (=0=>dont put in gradient corrections, just LDA)
      !  input lpot:  (=0=>dont get potential and dont need U and V)
      !  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
      !----------------------------------------------------------------------
      ! use precision_mod, only: rp
      !
      !.. Implicit Declarations ..
      implicit none
      !
      !.. Parameters ..
      class(xc), intent(in) :: this
      real(rp), parameter :: thrd = 1.d0/3.d0
      real(rp), parameter :: thrd4 = 4.d0/3.d0
      real(rp), parameter :: ax = -0.738558766382d0
      real(rp), parameter :: a1 = 0.041106d0
      real(rp), parameter :: a2 = 0.092070d0
      real(rp), parameter :: a3 = 0.657946d0
      real(rp), parameter :: a4 = 2.626712d0
      !
      !.. Formal Arguments ..
      integer, intent(in) :: lpot
      integer :: lgga
      real(rp), intent(in) :: rho, s, u, v
      real(rp), intent(out) :: ex, vx
      !
      !.. Local Scalars ..
      real(rp) :: exunif, fs, fss, fxlag, s4, xs, xsd, xsdd, ys, ysd, &
                  ysdd, zs, zsd, zsdd
      !
      ! ... Executable Statements ...
      !
      !----------------------------------------------------------------------
      ! construct LDA exchange energy density
      exunif = ax*(rho**thrd)
      !      IF(lgga.EQ.0) THEN
      !         ex=exunif
      !         vx=ex*thrd4
      !         RETURN
      !      ENDIF
      !----------------------------------------------------------------------
      ! construct LAG enhancement factor
      s4 = s**a4
      xs = a1*s4
      zs = 1.d0 + a2*s4
      ys = zs**a3
      fxlag = 1.d0 + xs/ys
      ex = exunif*fxlag
      vx = exunif*thrd4
      if (lpot /= 0) then
         !----------------------------------------------------------------------
         !  potential
         xsd = a4*xs/s
         xsdd = (a4 - 1.d0)*xsd/s
         zsd = a2*xsd/a1
         zsdd = a2*xsdd/a1
         ysd = a3*ys*zsd/zs
         ysdd = (a3 - 1.d0)*ysd*zsd/zs + ysd*zsdd/zsd
         fs = xsd/ys - xs*ysd/ys/ys
         fs = fs/s
         fss = xsdd/ys - 2.d0*xsd*ysd/ys/ys + 2.d0*xs*ysd*ysd/ys/ys/ys - &
               xs*ysdd/ys/ys
         fss = (fss - fs)/s
         !----------------------------------------------------------------------
         ! calculate potential
         vx = exunif*(thrd4*fxlag - (u - thrd4*s*s*s)*fss - v*fs)
      end if
   end subroutine exchlag
!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
   subroutine GCOR2(this, A, A1, B1, B2, B3, B4, rtrs, GG, GGRS)
      ! slimmed down version of GCOR used in PW91 routines, to interpolate
      ! LSD correlation energy, as given by (10) of
      ! J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
      ! K. Burke, May 11, 1996.
      ! use precision_mod, only: rp
      !
      !.. Implicit Declarations ..
      implicit none
      !
      !.. Formal Arguments ..
      class(xc), intent(in) :: this
      real(rp), intent(in) :: A, A1, B1, B2, B3, B4, rtrs
      real(rp), intent(out) :: GG, GGRS
      !
      !.. Local Scalars ..
      real(rp) :: Q0, Q1, Q2, Q3
      !
      !.. Intrinsic Functions ..
      intrinsic LOG
      !
      ! ... Executable Statements ...
      !
      Q0 = -2.d0*A*(1.d0 + A1*rtrs*rtrs)
      Q1 = 2.d0*A*rtrs*(B1 + rtrs*(B2 + rtrs*(B3 + B4*rtrs)))
      Q2 = LOG(1.d0 + 1.d0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/rtrs + 2.d0*B2 + rtrs*(3.d0*B3 + 4.d0*B4*rtrs))
      GGRS = -2.d0*A*A1*Q2 - Q0*Q3/(Q1*(1.d0 + Q1))
   end subroutine GCOR2

   subroutine DIFFN(this, RHO, RHOP, RHOPP, N, H)
      !   ******************************************************************
      !   *                                                                *
      !   *   Differentiate charge density for use in <XCPOT>.             *
      !   *                                                                *
      !   ******************************************************************
      ! use precision_mod, only: rp
      !
      !.. Implicit Declarations ..
      implicit none
      !
      !.. Formal Arguments ..
      class(xc), intent(in) :: this
      integer, intent(in) :: N
      real(rp), dimension(N), intent(in) :: RHO
      real(rp), dimension(N), intent(inout) :: RHOP, RHOPP
      real(rp), intent(in) :: H
      !
      !.. Local Scalars ..
      integer :: I, NM2
      real(rp) :: H2, SIXH, TWLWH, TWLWH2
      !
      ! ... Executable Statements ...
      !
      NM2 = N - 2
      SIXH = 6.d0*H
      TWLWH = 12.d0*H
      TWLWH2 = TWLWH*H
      H2 = H*H
      !
      !     Forward difference at the beginning of the table
      !
      RHOP(1) = ((2.d0*RHO(4) + 18.d0*RHO(2)) - (9.d0*RHO(3) + 11.d0*RHO(1)))/SIXH
      RHOP(2) = ((2.d0*RHO(5) + 18.d0*RHO(3)) - (9.d0*RHO(4) + 11.d0*RHO(2)))/SIXH
      RHOPP(1) = ((2.d0*RHO(1) + 4.d0*RHO(3)) - (5.d0*RHO(2) + RHO(4)))/H2
      RHOPP(2) = ((2.d0*RHO(2) + 4.d0*RHO(4)) - (5.d0*RHO(3) + RHO(5)))/H2
      !
      !     Central difference at the interior of the table
      !
      do I = 3, NM2
         RHOP(I) = ((RHO(I - 2) + 8.d0*RHO(I + 1)) - (8.d0*RHO(I - 1) + RHO(I + 2)))/TWLWH
         RHOPP(I) = &
            ((16.d0*RHO(I + 1) + 16.d0*RHO(I - 1)) - (RHO(I + 2) + RHO(I - 2) + 30.d0*RHO(I)))/ &
            TWLWH2
      end do
      !
      !     Backward difference at the end of the table
      !
      RHOP(N) = &
         ((11.d0*RHO(N) + 9.d0*RHO(N - 2)) - (18.d0*RHO(N - 1) + 2.d0*RHO(N - 3)))/SIXH
      RHOP(N - 1) = &
         ((11.d0*RHO(N - 1) + 9.d0*RHO(N - 3)) - (18.d0*RHO(N - 2) + 2.d0*RHO(N - 4)))/SIXH
      RHOPP(N) = ((2.d0*RHO(N) + 4.d0*RHO(N - 2)) - (5.d0*RHO(N - 1) + RHO(N - 3)))/H2
      RHOPP(N - 1) = ((2.d0*RHO(N - 1) + 4.d0*RHO(N - 3)) - (5.d0*RHO(N - 2) + RHO(N - 4)))/H2
      !do i=1, n
      !  write(808, ´(i6, 3g20.8)´) i, RHO(i), RHOP(i), RHOPP(i)
      !end do
   end subroutine DIFFN

   !==============================================================================
   !----- libXC INTERFACE PROCEDURES ---------------------------------------------
   !==============================================================================

   !>--------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Initialize libXC functional if needed based on txc value
   !> Validates compatibility with ASA before initialization
   !>--------------------------------------------------------------------------
   subroutine init_libxc(this, ctrl)
      class(xc), intent(inout) :: this
      type(control), intent(in) :: ctrl
      
      integer :: libxc_id, nspin

#ifdef HAVE_LIBXC
      type(xc_f03_func_t) :: temp_func
      character(len=256) :: func_name
#endif

      ! Check if we should use libXC for this functional
      libxc_id = this%get_libxc_functional_mapping()
      
      if (libxc_id > 0) then
         call g_logger%info('Initializing libXC functional ID: '//int2str(libxc_id), __FILE__, __LINE__)
         ! Validate functional compatibility with ASA
         if (.not. this%validate_libxc_compatibility(libxc_id)) then
            call g_logger%error('Stopping due to incompatible libXC functional (txc='//int2str(this%txc)//')', __FILE__, __LINE__)
            call g_logger%info('Valid options for ASA:', __FILE__, __LINE__)
            call g_logger%info('  Legacy functionals: txc=1-11 (LDA-PZ, PW92, PBE, BLYP, etc.)', __FILE__, __LINE__)  
            call g_logger%info('  libXC LDA functionals: txc=1001-1030', __FILE__, __LINE__)
            call g_logger%info('  libXC GGA functionals: txc=1101-1167', __FILE__, __LINE__)
            call g_logger%info('  See XC_FUNCTIONAL_GUIDE.md for complete list', __FILE__, __LINE__)
            stop 'Incompatible XC functional for ASA'
         endif
         
         this%use_libxc = .true.
         this%libxc_func_id = libxc_id
         
         if (this%use_libxc) then
            call g_logger%info('libXC initialized successfully: use_libxc = true', __FILE__, __LINE__)
         else
            call g_logger%info('libXC initialized successfully: use_libxc = false', __FILE__, __LINE__)
         endif
         
         ! Determine spin treatment - RS-LMTO always uses spin-polarized calculations
         ! (nsp=1 means non-relativistic spin-polarized in RS-LMTO convention)
         nspin = XC_POLARIZED
         
#ifdef HAVE_LIBXC
         ! Create a temporary functional to get info
         call xc_f03_func_init(temp_func, libxc_id, nspin)
         this%libxc_info = xc_f03_func_get_info(temp_func)
         this%libxc_family = xc_f03_func_info_get_family(this%libxc_info)
         this%libxc_nspin = nspin  ! Store for consistent usage
         
         ! Get the name before destroying the functional
         func_name = trim(xc_f03_func_info_get_name(this%libxc_info))
         
         call xc_f03_func_end(temp_func)
#endif
         
         ! Log information about the libXC functional
         call g_logger%info('Using libXC functional: '// &
                          trim(func_name)// &
                          ' (ID: '//int2str(libxc_id)//', Family: '//int2str(this%libxc_family)// &
                          ', nspin: '//int2str(nspin)//')', __FILE__, __LINE__)
      else
         this%use_libxc = .false.
      endif
   end subroutine init_libxc

   !>--------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Clean up libXC functional
   !>--------------------------------------------------------------------------
   subroutine cleanup_libxc(this)
      class(xc), intent(inout) :: this
      ! No cleanup needed since we use local functionals
      this%use_libxc = .false.
   end subroutine cleanup_libxc

   !>--------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Map RS-LMTO txc values to libXC functional IDs
   !> Returns -1 if no libXC equivalent exists
   !>--------------------------------------------------------------------------
   function get_libxc_functional_mapping(this) result(libxc_id)
      class(xc), intent(in) :: this
      integer :: libxc_id

#ifdef HAVE_LIBXC
      select case (this%txc)
      case (1001)              ! VWN via libXC (XC_LDA_XC_VWN = 17)
         libxc_id = 17
      case (1007)              ! VWN correlation via libXC (XC_LDA_C_VWN = 7)
         libxc_id = 7
      case (1009)              ! Perdew-Zunger via libXC (XC_LDA_C_PZ = 9)
         libxc_id = 9
      case (1101)              ! PBE exchange via libXC (XC_GGA_X_PBE = 101)
         libxc_id = 101
      case (1130)              ! PBE correlation via libXC (XC_GGA_C_PBE = 130)
         libxc_id = 130
      case (1167)              ! PBE exchange-correlation via libXC (XC_GGA_XC_PBE = 167)
         libxc_id = 167
      case default
         libxc_id = -1  ! No libXC mapping for this functional
      end select
#else
      libxc_id = -1  ! libXC not available
#endif
   end function get_libxc_functional_mapping

   !>--------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Validate that libXC functional is compatible with ASA (LDA/GGA only)
   !> Returns .true. if compatible, .false. otherwise
   !>--------------------------------------------------------------------------
   function validate_libxc_compatibility(this, libxc_id) result(is_compatible)
      class(xc), intent(in) :: this
      integer, intent(in) :: libxc_id
      logical :: is_compatible
      
#ifdef HAVE_LIBXC
      integer :: family
      character(len=256) :: family_name, func_name
      type(xc_f03_func_t) :: temp_func
      type(xc_f03_func_info_t) :: temp_info
      
      is_compatible = .false.
      
      if (libxc_id <= 0) then
         call g_logger%warning('Invalid libXC functional ID: '//int2str(libxc_id), __FILE__, __LINE__)
         return
      endif
      
      ! Initialize a temporary functional to check its family
      call xc_f03_func_init(temp_func, libxc_id, 1)  ! XC_UNPOLARIZED = 1
      temp_info = xc_f03_func_get_info(temp_func)
      
      family = xc_f03_func_info_get_family(temp_info)
      func_name = trim(xc_f03_func_info_get_name(temp_info))  ! Copy the string immediately
      
      select case(family)
      case(1)  ! XC_FAMILY_LDA = 1
         family_name = "LDA"
         is_compatible = .true.
      case(2)  ! XC_FAMILY_GGA = 2  
         family_name = "GGA"
         is_compatible = .true.
      case(3)  ! XC_FAMILY_MGGA = 3
         family_name = "meta-GGA"
         is_compatible = .false.
         call g_logger%warning('meta-GGA functional "'//trim(func_name)//'" not compatible with ASA spherical symmetry', __FILE__, __LINE__)
         call g_logger%warning('ASA lacks kinetic energy density (τ) needed for meta-GGA functionals', __FILE__, __LINE__)
      case(4)  ! XC_FAMILY_HYB_GGA = 4
         family_name = "hybrid GGA"
         is_compatible = .false.
         call g_logger%warning('Hybrid functional "'//trim(func_name)//'" not compatible with ASA implementation', __FILE__, __LINE__)
         call g_logger%warning('ASA lacks exact exchange implementation needed for hybrid functionals', __FILE__, __LINE__)
      case(5)  ! XC_FAMILY_HYB_MGGA = 5
         family_name = "hybrid meta-GGA"
         is_compatible = .false.
         call g_logger%warning('Hybrid meta-GGA functional "'//trim(func_name)//'" not compatible with ASA', __FILE__, __LINE__)
         call g_logger%warning('ASA lacks both kinetic energy density and exact exchange', __FILE__, __LINE__)
      case default
         family_name = "unknown"
         is_compatible = .false.
         call g_logger%warning('Unknown functional family ('//int2str(family)//') for "'//trim(func_name)//'"', __FILE__, __LINE__)
      end select
      
      if (is_compatible) then
         call g_logger%info('libXC functional "'//trim(func_name)//'" ('//trim(family_name)//') is compatible with ASA', __FILE__, __LINE__)
      else
         call g_logger%error('libXC functional "'//trim(func_name)//'" ('//trim(family_name)//') is NOT compatible with ASA', __FILE__, __LINE__)
         call g_logger%info('Please use LDA or GGA functionals only with ASA spherical symmetry', __FILE__, __LINE__)
      endif
      
      ! Clean up the temporary functional
      call xc_f03_func_end(temp_func)
#else
      is_compatible = .false.
      call g_logger%error('libXC not available - cannot validate functional compatibility', __FILE__, __LINE__)
#endif
   end function validate_libxc_compatibility

   !>--------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Check if txc value corresponds to a libXC functional  
   !>--------------------------------------------------------------------------
   function is_libxc_functional(this) result(is_libxc)
      class(xc), intent(in) :: this
      logical :: is_libxc
      
      ! libXC functionals use txc values >= 1000
      is_libxc = (this%txc >= 1000)
   end function is_libxc_functional

   !>--------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Hybrid XCPOT routine that dispatches to legacy or libXC backend
   !>--------------------------------------------------------------------------
   subroutine XCPOT_hybrid(this, RHO1, RHO2, RHO, RHOP, RHOPP, RR, V1, V2, EXC)
      class(xc), intent(in) :: this
      real(rp), intent(in) :: RHO, RHO1, RHO2, RR
      real(rp), intent(inout) :: EXC, V1, V2
      real(rp), dimension(2), intent(in) :: RHOP, RHOPP

      if (this%use_libxc) then
         call this%xcpot_libxc_wrapper(RHO1, RHO2, RHO, RHOP, RHOPP, RR, V1, V2, EXC)
      else
         call this%XCPOT(RHO1, RHO2, RHO, RHOP, RHOPP, RR, V1, V2, EXC)
      endif
   end subroutine XCPOT_hybrid

   !>--------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Wrapper to call libXC for exchange-correlation potential calculation
   !>--------------------------------------------------------------------------
   subroutine xcpot_libxc_wrapper(this, RHO1, RHO2, RHO, RHOP, RHOPP, RR, V1, V2, EXC)
      class(xc), intent(in) :: this
      real(rp), intent(in) :: RHO, RHO1, RHO2, RR
      real(rp), intent(inout) :: EXC, V1, V2
      real(rp), dimension(2), intent(in) :: RHOP, RHOPP

#ifdef HAVE_LIBXC
      real(rp), parameter :: TOLD = 1.d-20
      real(rp), parameter :: TOLDD = 1.d-20
      
      ! libXC arrays - use proper dimensions for spin treatment
      real(rp), dimension(2) :: rho_libxc  ! For spin-polarized: [rho_up, rho_down]
      real(rp), dimension(3) :: sigma_libxc ! For GGA: [grad_up^2, grad_up*grad_down, grad_down^2]
      real(rp), dimension(1) :: exc_libxc
      real(rp), dimension(2) :: vrho_libxc  ! Potentials for each spin
      real(rp), dimension(3) :: vsigma_libxc ! GGA gradient potentials
      type(xc_f03_func_t) :: temp_func  ! Temporary functional for testing
      integer :: nspin, family
      
      ! Initialize outputs
      V1 = 0.0d0
      V2 = 0.0d0  
      EXC = 0.0d0
      
      ! Check for negligible densities
      if (RHO1 < TOLD .or. RHO2 < TOLDD) then
         call this%XCPOT(RHO1, RHO2, RHO, RHOP, RHOPP, RR, V1, V2, EXC)
         return
      endif
      
      ! call g_logger%info('libXC called with RHO1='//real2str(RHO1)//', RHO2='//real2str(RHO2), __FILE__, __LINE__)
      
      ! Determine spin treatment - use the same as initialization
      nspin = this%libxc_nspin
      family = this%libxc_family
      
      ! Set up density arrays
      if (nspin == 2) then
         ! Spin-polarized calculation
         ! In RS-LMTO: RHO1 = spin-down, RHO2 = spin-up
         ! In libXC: rho_libxc(1) = spin-up, rho_libxc(2) = spin-down
         rho_libxc(1) = RHO2  ! spin-up density
         rho_libxc(2) = RHO1  ! spin-down density
      else
         ! Unpolarized calculation  
         rho_libxc(1) = RHO   ! total density
      endif
      
      ! Additional safety check for libXC call
      if (rho_libxc(1) <= 0.0d0) then
         call g_logger%info('Zero density detected, falling back to legacy', __FILE__, __LINE__)
         call this%XCPOT(RHO1, RHO2, RHO, RHOP, RHOPP, RR, V1, V2, EXC)
         return
      endif
      
      select case(family)
      case(1)  ! XC_FAMILY_LDA
         if (nspin == 2) then
            ! Create a local functional object for this calculation
            call xc_f03_func_init(temp_func, this%libxc_func_id, nspin)
            call xc_f03_lda_exc_vxc(temp_func, 1_c_size_t, rho_libxc, exc_libxc, vrho_libxc)
            call xc_f03_func_end(temp_func)
            
            V1 = vrho_libxc(2)  ! spin-down potential (for RHO1)
            V2 = vrho_libxc(1)  ! spin-up potential (for RHO2)
         else
            ! Create a local functional object for this calculation
            call xc_f03_func_init(temp_func, this%libxc_func_id, nspin)
            call xc_f03_lda_exc_vxc(temp_func, 1_c_size_t, rho_libxc, exc_libxc, vrho_libxc)
            call xc_f03_func_end(temp_func)
            
            V1 = vrho_libxc(1)  ! same potential for both spins
            V2 = vrho_libxc(1)
         endif
         EXC = exc_libxc(1)
         
      case(2)  ! XC_FAMILY_GGA  
         ! Set up gradient arrays
         if (nspin == 2) then
            ! RHOP(1) = d(rho_down)/dr, RHOP(2) = d(rho_up)/dr
            sigma_libxc(1) = RHOP(2)**2        ! |grad rho_up|^2
            sigma_libxc(2) = RHOP(2)*RHOP(1)   ! grad rho_up . grad rho_down  
            sigma_libxc(3) = RHOP(1)**2        ! |grad rho_down|^2
            
            ! Create local functional for GGA calculation
            call xc_f03_func_init(temp_func, this%libxc_func_id, nspin)
            call xc_f03_gga_exc_vxc(temp_func, 1_c_size_t, rho_libxc, sigma_libxc, &
                                    exc_libxc, vrho_libxc, vsigma_libxc)
            call xc_f03_func_end(temp_func)
            
            ! V1 for rho_down (RHO1), V2 for rho_up (RHO2)
            V1 = vrho_libxc(2) - 2.0d0 * vsigma_libxc(3) * RHOPP(1)
            V2 = vrho_libxc(1) - 2.0d0 * vsigma_libxc(1) * RHOPP(2)
         else
            sigma_libxc(1) = (RHOP(1) + RHOP(2))**2  ! |grad rho_total|^2
            
            ! Create local functional for GGA calculation
            call xc_f03_func_init(temp_func, this%libxc_func_id, nspin)
            call xc_f03_gga_exc_vxc(temp_func, 1_c_size_t, rho_libxc(1:1), sigma_libxc(1:1), &
                                    exc_libxc, vrho_libxc(1:1), vsigma_libxc(1:1))
            call xc_f03_func_end(temp_func)
            
            V1 = vrho_libxc(1) - 2.0d0 * vsigma_libxc(1) * (RHOPP(1) + RHOPP(2))
            V2 = V1
         endif
         EXC = exc_libxc(1)
         
      case default
         call g_logger%error('libXC functional family '//int2str(family)//' not supported in xcpot_libxc_wrapper', __FILE__, __LINE__)
         call g_logger%info('Supported families: 1 (LDA), 2 (GGA)', __FILE__, __LINE__)
         call g_logger%info('Falling back to legacy implementation', __FILE__, __LINE__)
         
         ! Fall back to legacy implementation
         call this%XCPOT(RHO1, RHO2, RHO, RHOP, RHOPP, RR, V1, V2, EXC)
         return
      end select
      
#else
      call g_logger%error('libXC not available - cannot use libXC functionals', __FILE__, __LINE__)
      stop 'libXC not available'
#endif
   end subroutine xcpot_libxc_wrapper

   ! DESCRIPTION:
   subroutine destructor(this)
      type(xc), intent(inout) :: this
      call this%cleanup_libxc()
   end subroutine destructor

end module xc_mod
