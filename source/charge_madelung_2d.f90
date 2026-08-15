submodule (charge_mod) charge_madelung_2d
   implicit none

contains

! Subroutines from surfmat library
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the Madelung potential matrix for a finite number
   !> of 2D layers of point charges and dipoles. This is used to
   !> set up the Madelung potential in the surface program INF.
   !>
   !> Matrix    (l´, l)
   !>
   !> DSS       (s, s)
   !> DSZ       (s, z)
   !> DS3Z2     (s, 3z2-1)
   !> DSX2Y2    (s, x2-y2)
   !> DSXY      (s, xy)
   !> DZZ       (z, z)
   !> DZ3Z2     (z, 3zz-1)
   !>
   !> PM    : Plate-condenser potential matrix
   !---------------------------------------------------------------------------
   module subroutine madl2d(this)

      ! Input
      class(charge), intent(inout) :: this
      ! Local variables
      integer :: I, IQ, JQ, JQP, NQM, nbas
      real(rp) :: ALPHA, AQPPZ, BETA, BETAM, BETAP, BMDL, BMG, DGI, &
                  DMDL, DMG, DRI, DRI2, DRI3, DRI4, DZ, ERFCA, ERFCM, &
                  ERFCP, EXPA, EXPAA, EXPEK, EXPP, EXPZ, FACBET, &
                  FACDIF, FACDK, FACDR, FACERF, FACG1, FACG2, &
                  FACGAU, FACGUA, FACP, FACQ1, FACQ2, FACQ3, FACQ4, &
                  FACQUA, FACQUR, FACZZ
      real(rp) :: G0MDL, GAM52, GLH, GLK, GLM, GM0, PHASE, Q0MDL, &
                  QIMDL, QM0G, QMIG, QMRG, QPPX, QPPY, QPPZ, QPX, QPY, &
                  QPZ, QRMDL, SUM0G, SUM0R, SUM1G, SUM1R, SUM20G, &
                  SUM20R, SUM2IG, SUM2IR, SUM2RG, SUM2RR, SUM30G, &
                  SUM30R, SUMM, SUMQ0, SUMQI, SUMQR, TWOLAM, TWOS, X, &
                  XI, XR
      real(rp) :: Y, Z, ZU, exf, exh, expm

      ! Process described in H. L. Skriver and N. M. Rosengaard Phys. Rev. B 43, 9538

      nbas = this%lattice%nbas
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('charge.ds3z2', this%ds3z2, (/nbas, nbas/))
      call g_safe_alloc%allocate('charge.dsx2y2', this%dsx2y2, (/nbas, nbas/))
      call g_safe_alloc%allocate('charge.dsxy', this%dsxy, (/nbas, nbas/))
      call g_safe_alloc%allocate('charge.dss', this%dss, (/nbas, nbas/))
      call g_safe_alloc%allocate('charge.dsz', this%dsz, (/nbas, nbas/))
      call g_safe_alloc%allocate('charge.dzz', this%dzz, (/nbas, nbas/))
      call g_safe_alloc%allocate('charge.dz3z2', this%dz3z2, (/nbas, nbas/))
      call g_safe_alloc%allocate('charge.am', this%am, (/nbas, nbas/))
      call g_safe_alloc%allocate('charge.bm', this%bm, (/nbas, nbas/))
      call g_safe_alloc%allocate('charge.pm', this%pm, (/nbas, nbas/))
#else
      allocate (this%ds3z2(nbas, nbas), this%dsx2y2(nbas, nbas), this%dsxy(nbas, nbas), this%dss(nbas, nbas), this%dsz(nbas, nbas) &
                , this%dzz(nbas, nbas), this%dz3z2(nbas, nbas), this%am(nbas, nbas), this%bm(nbas, nbas), this%pm(nbas, nbas))
#endif

      TWOS = 2.*this%SWS
      TWOLAM = 2.*this%ALAMDA
      FACBET = PI/this%AR2D/TWOLAM
      FACGAU = -2.*SQRT_PI/this%AR2D/this%ALAMDA
      FACQUA = 2.*this%SWS*this%SWS*SQRT(5.d0)/3.
      FACGUA = 6.*this%SWS*this%SWS*this%SWS/SQRT(15.d0)
      FACQUR = FACQUA*SQRT(3.d0/2.d0)
      FACERF = 2.*PI/this%AR2D
      FACDIF = 8.*PI*this%SWS/this%AR2D
      FACP = this%SWS*SQRT(3.d0)
      FACDK = FACP*PI/this%AR2D
      FACDR = FACP*2./SQRT_PI
      FACQ1 = 4.*SQRT_PI*this%ALAMDA/this%AR2D
      FACQ2 = 1.5*PI/this%AR2D
      FACQ3 = 2.d0/SQRT_PI
      FACQ4 = 0.5*PI/this%AR2D
      FACG1 = 8.*SQRT_PI*this%ALAMDA*this%ALAMDA/this%AR2D
      FACG2 = 2.5*PI/this%AR2D
      FACZZ = -6.d0/SQRT(5.d0)

!    write(100, ´(18f10.4)´) twos, twolam, facbet, facgau, facqua, facqur, facerf, facdif, &
!               facp, facdk, facdr, facq1, facq2, facq3, facq4, facg1, facg2, faczz

      !
      !     Layer-diagonal terms (R=R´)
      !
      SUMM = 0.d0
      SUMQ0 = -FACQ1
      SUMQR = 0.d0
      SUMQI = 0.d0
      ! Summation in the reciprocal space
      do I = 2, this%NUMVG
         DGI = this%DG(I)
         BETA = DGI/TWOLAM
         ERFCA = ERFC(BETA)
         SUMM = SUMM + 2.*ERFCA/BETA
         SUMQ0 = SUMQ0 - FACQ1*EXP(-BETA*BETA) + 2.*FACQ2*DGI*ERFCA
         SUMQR = SUMQR - 2.*FACQ4*(this%AKX(I)*this%AKX(I) - this%AKY(I)*this%AKY(I))*ERFCA/DGI
         SUMQI = SUMQI - 4.*FACQ4*this%AKX(I)*this%AKY(I)*ERFCA/DGI
      end do
      BMDL = FACBET*SUMM
      Q0MDL = FACQUA*SUMQ0
      QRMDL = FACQUR*SUMQR
      QIMDL = FACQUR*SUMQI
      !
      SUMM = 0.d0
      SUMQ0 = 0.d0
      SUMQR = 0.d0
      SUMQI = 0.d0
      ! Summation in the translational vector of the 2D lattice
      do I = 2, this%NR0
         DRI = this%DR(I)
         DRI2 = DRI*DRI
         DRI3 = DRI2*DRI
         ALPHA = this%ALAMDA*DRI
         ERFCA = ERFC(ALPHA)
         EXPAA = EXP(-ALPHA*ALPHA)
         GAM52 = 1.5*(ERFCA/FACQ3 + ALPHA*EXPAA) + EXPAA*ALPHA**3
         GAM52 = GAM52/DRI3
         SUMM = SUMM + ERFCA/DRI
         SUMQ0 = SUMQ0 - GAM52
         SUMQR = SUMQR + GAM52*(this%ASX(I)*this%ASX(I) - this%ASY(I)*this%ASY(I))/DRI2
         SUMQI = SUMQI + GAM52*2.*this%ASX(I)*this%ASY(I)/DRI2
      end do
      BMDL = BMDL + SUMM - TWOLAM/SQRT_PI
      Q0MDL = Q0MDL + FACQUA*FACQ3*SUMQ0
      QRMDL = QRMDL + FACQUR*FACQ3*SUMQR
      QIMDL = QIMDL + FACQUR*FACQ3*SUMQI
      !
      do IQ = 1, this%lattice%nbas
         this%AM(IQ, IQ) = FACGAU
         this%BM(IQ, IQ) = BMDL
         this%DSZ(IQ, IQ) = 0.d0
         this%DS3Z2(IQ, IQ) = Q0MDL
         this%DSX2Y2(IQ, IQ) = QRMDL
         this%DSXY(IQ, IQ) = QIMDL
         this%DZ3Z2(IQ, IQ) = 0.d0
      end do
      !
      !     Off-diagonal terms
      !
      NQM = this%lattice%nbas - 1
      do JQ = 1, NQM
         JQP = JQ + 1
         QPX = this%QX(JQ)
         QPY = this%QY(JQ)
         QPZ = this%QZ(JQ)
         do IQ = JQP, this%lattice%nbas
            QPPX = this%QX(IQ) - QPX
            QPPY = this%QY(IQ) - QPY
            QPPZ = this%QZ(IQ) - QPZ
            !
            !     Contribution to potential from k-parallel equal to zero
            !
            DZ = this%ALAMDA*QPPZ
            ERFCP = ERFC(DZ)
            ERFCM = 2.-ERFCP
            if (DZ > 12.) then
               EXPZ = 0.d0
            else
               EXPZ = EXP(-DZ*DZ)
            end if
            this%AM(IQ, JQ) = FACGAU*EXPZ - QPPZ*FACERF*ERFCM
            this%AM(JQ, IQ) = FACGAU*EXPZ + QPPZ*FACERF*ERFCP
            SUM0G = 0.d0
            SUM1G = ERFCM - ERFCP
            SUM20G = -FACQ1*EXPZ
            SUM2RG = 0.d0
            SUM2IG = 0.d0
            SUM30G = FACG1*DZ*EXPZ
            !
            !     Contributions from k-parallel greater than zero
            !
            do I = 2, this%NUMVG
               DGI = this%DG(I)
               PHASE = COS(this%AKX(I)*QPPX + this%AKY(I)*QPPY)
               BETA = DGI/TWOLAM
               AQPPZ = this%ALAMDA*QPPZ
               BETAP = BETA + AQPPZ
               BETAM = BETA - AQPPZ
               ERFCP = ERFC(BETAP)
               ERFCM = ERFC(BETAM)
               ! alterado por Peduto para evitar alguns casos de "overflow" que podem
               ! ocorrer quando, principalmente, o sistema se tratar de um VAX.
               ! a variavel EXPP estoura, geralmente, quando ercp e erfcm sao nulos.
               ! Mas
               ! nestes casos nao eh necessario calcular EXPP. Assim evita-se que o
               ! progra
               ! ma seja abortado.
               if (erfcp == 0) then
                  ! ERFCP underflowed to exactly zero, so EXPP*ERFCP = 0 and the
                  ! diverging EXPP need never be formed. The general branch below
                  ! then reduces to EXF = EXPM*ERFCM and EXH = -EXPM*ERFCM = -EXF.
                  expm = exp(-dgi*qppz)
                  exf = expm*erfcm
                  exh = -exf
               else
                  EXPP = EXP(DGI*QPPZ)
                  EXPM = 1./EXPP
                  EXF = EXPP*ERFCP + EXPM*ERFCM
                  EXH = EXPP*ERFCP - EXPM*ERFCM
               end if
               EXPEK = EXP(-AQPPZ*AQPPZ)*EXP(-BETA*BETA)
               BMG = PHASE*EXF/BETA
               DMG = PHASE*EXH
               QM0G = PHASE*(FACQ2*DGI*EXF - FACQ1*EXPEK)
               QMRG = PHASE*(this%AKX(I)*this%AKX(I) - this%AKY(I)*this%AKY(I))*EXF/DGI
               QMIG = PHASE*this%AKX(I)*this%AKY(I)*EXF/DGI
               GM0 = PHASE*(FACG2*DGI*DGI*EXH + FACG1*DZ*EXPEK)
               SUM0G = SUM0G + BMG
               SUM1G = SUM1G - DMG
               SUM20G = SUM20G + QM0G
               SUM2RG = SUM2RG + QMRG
               SUM2IG = SUM2IG + QMIG
               SUM30G = SUM30G + GM0
            end do
            BMDL = FACBET*SUM0G
            DMDL = FACDK*SUM1G
            Q0MDL = FACQUA*SUM20G
            QRMDL = FACQUR*FACQ4*SUM2RG
            QIMDL = FACQUR*FACQ4*2.*SUM2IG
            G0MDL = FACGUA*SUM30G
            !
            SUM0R = 0.d0
            SUM1R = 0.d0
            SUM20R = 0.d0
            SUM2RR = 0.d0
            SUM2IR = 0.d0
            SUM30R = 0.d0
            do I = 1, this%NUMVR
               X = this%ASX(I) + QPPX
               Y = this%ASY(I) + QPPY
               Z = QPPZ
               DRI2 = X*X + Y*Y + Z*Z
               DRI = SQRT(DRI2)
               DRI3 = DRI*DRI2
               DRI4 = DRI*DRI3
               ZU = -Z/DRI
               XR = (X*X - Y*Y)/DRI2
               XI = 2.*X*Y/DRI2
               if (DRI < this%RMAX) then
                  ALPHA = this%ALAMDA*DRI
                  ERFCA = ERFC(ALPHA)
                  SUM0R = SUM0R + ERFCA/DRI
                  EXPA = EXP(-ALPHA*ALPHA)
                  GLH = 0.5*SQRT_PI*ERFCA + ALPHA*EXPA
                  GLK = 1.5*GLH + EXPA*ALPHA**3
                  GLM = 2.5*GLK + EXPA*ALPHA**5
                  SUM1R = SUM1R - GLH*ZU/DRI2
                  SUM20R = SUM20R + GLK*(3.*ZU*ZU - 1.)/DRI3
                  SUM2RR = SUM2RR + GLK*XR/DRI3
                  SUM2IR = SUM2IR + GLK*XI/DRI3
                  SUM30R = SUM30R + GLM*ZU*(5.*ZU*ZU - 3.)/DRI4
               end if
            end do
            BMDL = BMDL + SUM0R
            DMDL = DMDL + FACDR*SUM1R
            Q0MDL = Q0MDL + FACQUA*FACQ3*SUM20R
            QRMDL = QRMDL + FACQUR*FACQ3*SUM2RR
            QIMDL = QIMDL + FACQUR*FACQ3*SUM2IR
            G0MDL = G0MDL + FACGUA*2.*FACQ3*SUM30R
            this%BM(IQ, JQ) = BMDL
            this%BM(JQ, IQ) = BMDL
            this%DSZ(IQ, JQ) = DMDL
            this%DSZ(JQ, IQ) = -DMDL
            this%DS3Z2(IQ, JQ) = Q0MDL
            this%DS3Z2(JQ, IQ) = Q0MDL
            this%DSX2Y2(IQ, JQ) = QRMDL
            this%DSX2Y2(JQ, IQ) = QRMDL
            this%DSXY(IQ, JQ) = QIMDL
            this%DSXY(JQ, IQ) = QIMDL
            ! DZZ is deliberately NOT set here: it is built below from the
            ! already-TWOS-scaled DS3Z2. Setting it from the unscaled Q0MDL at
            ! this point was a dead assignment (unconditionally overwritten),
            ! and a trap for whoever turns l = 2 on. See CONVENTIONS_MADELUNG.md.
            this%DZ3Z2(IQ, JQ) = G0MDL
            this%DZ3Z2(JQ, IQ) = -G0MDL
            !
         end do
      end do
      !
      !     Potential in units of e**2/(2S)
      !
      do JQ = 1, this%lattice%nbas
         do IQ = 1, this%lattice%nbas
            this%AM(IQ, JQ) = TWOS*this%AM(IQ, JQ)
            this%BM(IQ, JQ) = TWOS*this%BM(IQ, JQ)
            this%DSS(IQ, JQ) = this%AM(IQ, JQ) + this%BM(IQ, JQ)
            this%DSZ(IQ, JQ) = TWOS*this%DSZ(IQ, JQ)
            this%DS3Z2(IQ, JQ) = TWOS*this%DS3Z2(IQ, JQ)
            this%DSX2Y2(IQ, JQ) = TWOS*this%DSX2Y2(IQ, JQ)
            this%DSXY(IQ, JQ) = TWOS*this%DSXY(IQ, JQ)
            this%DZZ(IQ, JQ) = FACZZ*this%DS3Z2(IQ, JQ)
            this%DZ3Z2(IQ, JQ) = TWOS*this%DZ3Z2(IQ, JQ)
         end do
      end do
      !
      !     Plate-condenser matrix
      !
      !     Built here but consumed nowhere. This is the applied-bias / electric
      !     field hook (B7 §6): with alignment entering as a target potential
      !     step, a bias is target_step = contact_potential + eV, and PM is the
      !     matching kernel. Two caveats for whoever activates it:
      !       * PM carries the SAME lower-triangular gauge choice as AM, so it
      !         needs the same treatment (see CONVENTIONS_MADELUNG.md, C3);
      !       * a biased or field-exposed slab is a charged-slab problem, so the
      !         Q = 0 precondition no longer rescues the gauge and the setup
      !         needs its own compensation story.
      !
      do JQ = 1, this%lattice%nbas
         do IQ = 1, this%lattice%nbas
            if (IQ > JQ) then
               this%PM(IQ, JQ) = FACDIF*(this%QZ(JQ) - this%QZ(IQ))
            else
               this%PM(IQ, JQ) = 0.d0
            end if
         end do
      end do
   end subroutine madl2d

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> TODO
   !---------------------------------------------------------------------------
   module subroutine madl2r(this)
      class(charge), intent(inout) :: this
      ! Local variables
      real(rp) :: alpha, arg, beta, dgi, dri, dz, epsz, erfca, erfcm, &
                  erfcp, erfcz, erfr, expp, facbet, facgau, phase, r, &
                  rx, ry, sumg, sumr, twolam, twos, x, y, z
      real(rp), dimension(51) :: vg, vr, vz, vz0, zz
      integer :: i, iz, nz

      TWOS = 2.*this%SWS
      TWOLAM = 2.*this%ALAMDA
      FACBET = TWOS*PI/this%AR2D
      FACGAU = 1./SQRT_PI/this%ALAMDA
      NZ = 21
      DZ = TWOS/(NZ - 1)
      X = 0.
      Y = 0.
      do IZ = 1, NZ
         Z = -this%SWS + (IZ - 1)*DZ
         ZZ(IZ) = Z
         !
         !     Diagonal terms
         !
         EPSZ = this%ALAMDA*Z
         ARG = -EPSZ*EPSZ
         ERFCZ = ERFC(-EPSZ)
         SUMG = -2.*(Z*ERFCZ + FACGAU*EXP(ARG))
         VZ0(IZ) = FACBET*SUMG
         do I = 2, this%NUMVG
            DGI = this%DG(I)
            ARG = this%AKX(I)*X + this%AKY(I)*Y
            PHASE = COS(ARG)
            BETA = DGI/TWOLAM
            ERFCP = ERFC(BETA + EPSZ)
            ERFCM = ERFC(BETA - EPSZ)
            ! alteracao feita por Peduto para evitar alguns casos de "overflow" que
            ! podem
            ! ocorrer principalmente quando o sistema se tratar de um VAX
            EXPP = EXP(DGI*Z)
            SUMG = SUMG + PHASE*(EXPP*ERFCP + ERFCM/EXPP)/DGI
         end do
         VG(IZ) = FACBET*SUMG
         R = SQRT(X*X + Y*Y + Z*Z)
         ARG = this%ALAMDA*R
         if (ABS(Z) < 1.d-6) then
            ERFR = TWOLAM/SQRT_PI
         else
            ERFR = (1.-ERFC(ARG))/R
         end if
         SUMR = -ERFR
         do I = 2, this%NR0
            RX = this%ASX(I) - X
            RY = this%ASY(I) - Y
            DRI = SQRT(RX*RX + RY*RY + Z*Z)
            ALPHA = this%ALAMDA*DRI
            ERFCA = ERFC(ALPHA)
            SUMR = SUMR + ERFCA/DRI
         end do
         VR(IZ) = TWOS*SUMR
         VZ(IZ) = FACBET*SUMG + TWOS*SUMR
      end do
   end subroutine madl2r

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> TODO
   !---------------------------------------------------------------------------
   module subroutine latt2d(this)
      class(charge), intent(inout) :: this
      ! Local variables
      integer, parameter :: mr = 50000
      real(rp) :: a, amin, b, da, db, ddm, dkm, dq, dx, g1, ga, gx, gy, &
                  pqx, pqy, pqz, r1, ra, sx, sy, x, y, z
      real(rp), dimension(3) :: dd, dk
      real(rp), dimension(mr) :: csx, csy, d
      integer :: i, iq, jq, k, l, m, n, n1, ng, nr, nsh, nshl, numg, numgh, numr, numrh

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('charge.asx', this%asx, mr)
      call g_safe_alloc%allocate('charge.asy', this%asy, mr)
      call g_safe_alloc%allocate('charge.asz', this%asz, mr)
      call g_safe_alloc%allocate('charge.akx', this%akx, mr)
      call g_safe_alloc%allocate('charge.aky', this%aky, mr)
      call g_safe_alloc%allocate('charge.akz', this%akz, mr)
      call g_safe_alloc%allocate('charge.dr', this%dr, mr)
      call g_safe_alloc%allocate('charge.dg', this%dg, mr)
#else
      allocate (this%asx(mr), this%asy(mr), this%asz(mr), this%akx(mr), this%aky(mr), this%akz(mr), this%dr(mr), this%dg(mr))
#endif

      ! Calculate the radius RA of the sphere holding all the vectors  used in the
      ! lattice summations. R1 is the longest basis vector
      r1 = 1.e-06
      do IQ = 1, this%lattice%nbas
         PQX = this%QX(IQ)
         PQY = this%QY(IQ)
         PQZ = this%QZ(IQ)
         do JQ = IQ, this%lattice%nbas
            X = PQX - this%QX(JQ)
            Y = PQY - this%QY(JQ)
            Z = PQZ - this%QZ(JQ)
            DQ = SQRT(X*X + Y*Y + Z*Z)
            if (DQ >= R1) then
               R1 = DQ
            end if
         end do
      end do
      R1 = R1*1.001
      RA = this%RMAX + R1
      G1 = 0.d0
      GA = this%GMAX + G1

      do I = 1, 3
         DD(I) = SQRT(this%BSX(I)**2 + this%BSY(I)**2 + this%BSZ(I)**2)
         DK(I) = SQRT(this%BKX(I)**2 + this%BKY(I)**2 + this%BKZ(I)**2)
      end do
      DDM = MAX(DD(1), DD(2), DD(3))
      DKM = MAX(DK(1), DK(2), DK(3))
      DDM = 2*PI/DDM
      DKM = 2*PI/DKM
      NUMR = 2*(INT(RA/DKM) + 1) + 1
      NUMG = 2*(INT(GA/DDM) + 1) + 1
      NUMRH = NUMR/2 + 1
      NUMGH = NUMG/2 + 1
      ! Real space
      NR = 0
      this%NR0 = 0
      do L = 1, NUMR
         A = L - NUMRH
         do M = 1, NUMR
            B = M - NUMRH
            SX = A*this%BSX(1) + B*this%BSX(2)
            SY = A*this%BSY(1) + B*this%BSY(2)
            DX = SQRT(SX*SX + SY*SY)
            if (DX <= RA) then
               if (DX <= this%RMAX) then
                  this%NR0 = this%NR0 + 1
               end if
               NR = NR + 1
               if (NR > MR) then
                  write (6, 10005) NR, DX, RA, this%RMAX
                  stop
               else
                  D(NR) = DX
                  CSX(NR) = SX
                  CSY(NR) = SY
               end if
            end if
         end do
      end do
      !
      !     Sort vectors in order of increasing length
      !
      DA = 1.d-06
      NSH = 0
      NSHL = -1
      do K = 1, NR
         AMIN = 1000.
         do N = 1, NR
            if (D(N) < AMIN) then
               AMIN = D(N)
               N1 = N
            end if
         end do
         NSHL = NSHL + 1
         this%ASX(K) = CSX(N1)
         this%ASY(K) = CSY(N1)
         this%ASZ(K) = 0.d0
         DB = D(N1)
         this%DR(K) = DB
         if (DB > DA + 1.d-06) then
            NSH = NSH + 1
            NSHL = 0
            DA = DB
         end if
         D(N1) = 1000.
      end do
      NSH = NSH + 1
      NSHL = NSHL + 1
      this%NUMVR = NR
      !
      !     Reciprocal space
      !
      NG = 0
      do L = 1, NUMG
         A = L - NUMGH
         do M = 1, NUMG
            B = M - NUMGH
            GX = A*this%BKX(1) + B*this%BKX(2)
            GY = A*this%BKY(1) + B*this%BKY(2)
            DX = SQRT(GX*GX + GY*GY)
            if (DX <= GA) then
               NG = NG + 1
               if (NG > MR) then
                  write (6, 10006) NG, DX, GA, this%GMAX
                  stop
               else
                  D(NG) = DX
                  CSX(NG) = GX
                  CSY(NG) = GY
               end if
            end if
         end do
      end do
      !
      !     Sort vectors in order of increasing length
      !
      DA = 1.e-06
      NSH = 0
      NSHL = -1
      do K = 1, NG
         AMIN = 1000.
         do N = 1, NG
            if (D(N) < AMIN) then
               AMIN = D(N)
               N1 = N
            end if
         end do
         NSHL = NSHL + 1
         this%AKX(K) = CSX(N1)
         this%AKY(K) = CSY(N1)
         this%AKZ(K) = 0.d0
         DB = D(N1)
         this%DG(K) = DB
         if (DB > DA*1.000001) then
            NSH = NSH + 1
            NSHL = 0
            DA = DB
         end if
         D(N1) = 1000.
      end do
      NSH = NSH + 1
      NSHL = NSHL + 1
      this%NUMVG = NG

10005 format( &
         /, " LATT2D:** NR =", i5, " exceeds MR. Decrease RMAX by", &
         " changing ALAMDA or AMAX.", /, 11x, "Last vector: (ASX, ASY, ", "ASZ) =" &
         , 3f10.4)
10006 format( &
         /, " LATT2D:** NG =", i5, " exceeds MK. Decrease GMAX by", &
         " changing ALAMDA or BMAX.", /, 11x, "Last vector: (AKX, AKY, ", "AKZ) =" &
         , 3f10.4)
   end subroutine latt2d

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> TODO
   !---------------------------------------------------------------------------
   module subroutine set2d(this)
      class(charge) :: this
      ! Local variables
      real(rp) :: bx, by, bz, dperp
      real(rp), dimension(this%lattice%nbas*this%nq3) :: pz, uqx, uqy, uqz
      integer, dimension(this%lattice%nbas*this%nq3) :: index
      integer :: i, ib, inx, iq, isrf, n, nlam, nlama, nlamb, nq

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('charge.qx', this%qx, this%lattice%nbas)
      call g_safe_alloc%allocate('charge.qy', this%qy, this%lattice%nbas)
      call g_safe_alloc%allocate('charge.qz', this%qz, this%lattice%nbas)
#else
      allocate (this%qx(this%lattice%nbas), this%qy(this%lattice%nbas), this%qz(this%lattice%nbas))
#endif

      NLAM = this%lattice%nbas
      NLAMB = NLAM/2
      if (2*NLAMB == NLAM) then
         NLAMA = NLAMB - 1
      else
         NLAMA = NLAMB
      end if
      this%AR2D = this%BSX(1)*this%BSY(2) - this%BSY(1)*this%BSX(2)
      this%AR2D = sqrt(this%ar2d**2)
      DPERP = this%VOL/this%AR2D

      !
      !     Establish the interface layers
      !
      N = 0
      do IB = -NLAMA, NLAMB
         BX = IB*this%BSX(3)
         BY = IB*this%BSY(3)
         BZ = IB*this%BSZ(3)
         do IQ = 1, this%NQ3
            N = N + 1
            UQX(N) = BX + this%QX3(IQ)
            UQY(N) = BY + this%QY3(IQ)
            UQZ(N) = BZ + this%QZ3(IQ)
            PZ(N) = UQZ(N)
         end do
      end do
      call QSORT(PZ, INDEX, N)
      ISRF = 0
      do I = 1, N
         INX = INDEX(I)
         if (ABS(PZ(I)) < 1.d-6 .and. ISRF == 0) then
            ISRF = I
         end if
      end do
      IQ = 0
      do I = ISRF - NLAMA, ISRF + NLAMB
         INX = INDEX(I)
         IQ = IQ + 1
         this%QX(IQ) = UQX(INX)
         this%QY(IQ) = UQY(INX)
         this%QZ(IQ) = UQZ(INX)
      end do
   end subroutine set2d

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> TODO
   !
   !> @param[inout] m
   !> @param[inout] inx
   !> @param[inout] nl
   !---------------------------------------------------------------------------
   module subroutine qsort(m, inx, nl)
      ! sorts the vector m in increasing order
      implicit none
      ! inputs and outputs
      integer, intent(in) :: NL
      real(rp), dimension(NL), intent(inout) :: M
      integer, dimension(NL), intent(inout) :: INX
      ! Local variables
      integer :: I, IND, INIC, J, K
      real(rp) :: FIM, Z

      IND = 1
      INIC = 2
      FIM = NL
      do I = 1, NL
         INX(I) = I
      end do
      do while ((IND == 1) .and. (INIC <= FIM))
         IND = 0
         do J = INT(FIM), INIC, -1
            if (M(J) < M(J - 1)) then
               Z = M(J)
               K = INX(J)
               M(J) = M(J - 1)
               INX(J) = INX(J - 1)
               M(J - 1) = Z
               INX(J - 1) = K
            end if
         end do
         FIM = FIM - 1
         do J = INIC, INT(FIM)
            if (M(J + 1) < M(J)) then
               Z = M(J + 1)
               K = INX(J + 1)
               M(J + 1) = M(J)
               INX(J + 1) = INX(J)
               M(J) = Z
               INX(J) = K
               IND = 1
            end if
         end do
         INIC = INIC + 1
      end do
   end subroutine qsort

end submodule charge_madelung_2d
