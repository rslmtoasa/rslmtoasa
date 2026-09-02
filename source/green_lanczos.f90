!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! SUBMODULE: green_lanczos
!
!> Recursive Green-function implementations moved from green_mod.
!------------------------------------------------------------------------------

submodule (green_mod) green_lanczos

   use mpi_mod
   use logger_mod, only: g_logger
   use string_mod, only: fmt, int2str, real2str, log2str
   implicit none

contains

   module subroutine sgreen(this)
      use mpi_mod
      !this, bdos, g0, e, nv, ldim, ll, na, nmdir, mom, ef)
      class(green) :: this
      ! Input
      !integer, intent(in) :: nv, ldim, ll, na, nmdir
      !real(rp), intent(in) :: ef
      ! Output
      !real(rp), dimension(nv), intent(inout) :: e
      !real(rp), dimension(nv, ldim, ldim, na), intent(out) :: bdos
      !complex(rp), dimension(nv, ldim, ldim, na), intent(out) :: g0
      !real(rp), dimension(na, 3), intent(inout) :: mom
      ! Local variables
      integer :: ia, ja, mdir, nw, ll_t, ie, j, i
      real(rp), dimension(this%en%channels_ldos + 10, this%lattice%ntype) :: dx, dy, dz
      real(rp), dimension(nb, this%en%channels_ldos + 10) :: doso
      real(rp), dimension(nb, this%en%channels_ldos + 10) :: dmag, dnmag
      complex(rp) :: dfac, sfac, impi
      complex(rp), dimension(4) :: gspinor
      complex(rp), dimension(3) :: gmask, lmask
      complex(rp), dimension(2, 3) :: gfac
      integer, dimension(4, 3) :: goff
      real(rp), dimension(this%control%lld, this%control%lld)  :: Sm
      real(rp), dimension(nb)  :: q_int
      real(rp) :: e_start, e_stop

      integer :: ia_glob

      impi = (pi, 0.0d0)
      gfac(1, 1) = 1.0d0; gfac(1, 2) = -i_unit; gfac(1, 3) = 1.0d0
      gfac(2, 1) = 1.0d0; gfac(2, 2) = i_unit; gfac(2, 3) = -1.0d0
      gmask(1) = 0.0d0; gmask(2) = 0.0d0; gmask(3) = 1.0d0
      !   if(lrot) then
      !      lmask(1)=0.0d0;lmask(2)=0.0d0;lmask(3)=1.0d0
      !   else
      lmask(1) = 1.0d0/3.0d0; lmask(2) = 1.0d0/3.0d0; lmask(3) = 1.0d0/3.0d0
      !   end if

      goff(1, 1) = 0; goff(2, 1) = 9; goff(1, 2) = 0; goff(2, 2) = 9; goff(1, 3) = 0; goff(2, 3) = 0
      goff(3, 1) = 9; goff(4, 1) = 0; goff(3, 2) = 9; goff(4, 2) = 0; goff(3, 3) = 9; goff(4, 3) = 9
      dfac = i_unit*impi/(2.0d0, 0.0d0)
      dx = 0; dy = 0; dz = 0
      nw = 10*this%lattice%ntype*this%control%lld
      ll_t = this%control%lld

      doso = 0.0d0
      this%g0 = 0.0d0
      !do ia = 1, this%lattice%nrec
      do ia_glob = start_atom, end_atom
         ia = g2l_map(ia_glob)
         do mdir = 1, this%control%nmdir
            doso = 0.0d0
            call this%dos%density(doso, ia, mdir)
            if (this%control%nmdir == 1) then
               do ie = 1, this%en%channels_ldos + 10
                  do j = 1, nb
                     this%g0(j, j, ie, ia) = -i_unit*doso(j, ie)*impi
                  end do
                  write (300 + ia, *) this%en%ene(ie), sum(doso(1:norb, ie)), sum(doso(norb+1:nb, ie))
               end do
            else
               do ie = 1, this%en%channels_ldos + 10
                  do j = 1, norb
                     ! Charge, from main direction.. (not z-component)
                     this%g0(j, j, ie, ia) = this%g0(j, j, ie, ia) - (doso(j, ie) + doso(j +spin_off, ie))*dfac*lmask(mdir)!*1.0d0 /3.0d0 !mom(ja, mdir)**2
                     this%g0(j +spin_off, j +spin_off, ie, ia) = this%g0(j +spin_off, j +spin_off, ie, ia) - (doso(j, ie) + doso(j +spin_off, ie))*dfac*lmask(mdir)!*1.0d0 /3.0d0 !mom(ja, mdir)**2
                     ! Spin dependent part
                     this%g0(j + goff(1, mdir), j + goff(2, mdir), ie, ia) = &
                        this%g0(j + goff(1, mdir), j + goff(2, mdir), ie, ia) - (doso(j, ie) - doso(j +spin_off, ie))*gfac(1, mdir)*dfac

                     this%g0(j + goff(3, mdir), j + goff(4, mdir), ie, ia) = &
                        this%g0(j + goff(3, mdir), j + goff(4, mdir), ie, ia) - (doso(j, ie) - doso(j +spin_off, ie))*gfac(2, mdir)*dfac
                  end do
               end do
            end if
         end do
      end do
   end subroutine sgreen

   module subroutine bgreen(this, g_out, i_site, ie_start, ie_len, a_inf, b_inf, eta)
      class(green), intent(inout) :: this
      integer, intent(in) :: i_site, ie_start, ie_len
      complex(rp), dimension(:, :, :), intent(inout) :: g_out
      real(rp), dimension(nb, nb), intent(in) :: a_inf, b_inf
      complex(rp), intent(in) :: eta
      complex(rp), allocatable :: z_grid(:)
      integer :: ie

      if (ie_start < 1 .or. ie_len < 1 .or. ie_start+ie_len-1 > size(g_out, 3) .or. &
          ie_start+ie_len-1 > size(this%en%ene)) then
         error stop 'bgreen: requested energy range is outside the Green-function storage'
      end if
      allocate(z_grid(ie_len))
      do ie = 1, ie_len
         z_grid(ie) = cmplx(this%en%ene(ie_start+ie-1), 0.0_rp, rp) + eta
      end do
      g_out = cmplx(0.0_rp, 0.0_rp, rp)
      call this%bgreen_complex(g_out(:, :, ie_start:ie_start+ie_len-1), i_site, z_grid, a_inf, b_inf, .true.)
   end subroutine bgreen

   !> Shared block-recursion core for the legacy energy mesh and arbitrary
   !> complex energies.  The optional legacy flag preserves the historical
   !> real-axis terminator and output spelling for existing consumers; the
   !> public complex path evaluates the terminator analytically at every z.
   module subroutine bgreen_complex(this, g_out, i_site, z_grid, a_inf, b_inf, legacy_real_axis)
      implicit none
      class(green), intent(inout) :: this
      integer, intent(in) :: i_site
      complex(rp), dimension(:, :, :), intent(inout) :: g_out
      real(rp), dimension(nb, nb), intent(in) :: a_inf
      real(rp), dimension(nb, nb), intent(in) :: b_inf
      complex(rp), intent(in) :: z_grid(:)
      logical, intent(in), optional :: legacy_real_axis
      !
      integer :: nv, ldim, ll, ll_t, llinf
      real(rp) :: factor_z
      integer :: i, j, l, ei, info, ln, lwork
      integer :: ii, jj
      real(rp) :: recMaxA, recMaxB, valrec
      real(rp) :: maxQ, maxB2r
      real(rp) :: qMax, wMax, b2Max
      real(rp) :: qval, wval, b2val
      logical :: qHasNaN
      logical :: recNaN
      logical :: found_nan
      complex(rp) :: vv
      integer, dimension(nb) :: ipiv
      complex(rp) :: etop, ebot, ea, eb
      real(rp), dimension(this%lattice%nrec) :: a_inf0, b_inf0
      complex(rp), dimension(nb, nb) :: Q, Qp, Q2p, Z, one, W, B2z, Qt, P
      complex(rp) :: zoff, cone, czero, ze, zterm, det, im
      complex(rp) :: coff
      complex(rp), dimension(nb*nb) :: work
      real(rp), dimension(this%en%channels_ldos + 10) :: ene
      complex(rp), dimension(nb, nb) :: Dfac_mat, Cshi_mat
      real(rp) :: a_diag, b_diag
      complex(rp) :: z_eval, eta_local
      logical :: legacy_mode
      !
      integer, dimension(nb) :: m_tab
      !
      integer :: n_glob
      logical, save :: dump_done = .false.
      character(len=200) :: dump_fname
      integer :: dump_unit, dump_iostat
      !
      g_out = (0.0d0, 0.0d0)
      !
      ! Definitions so it is not necessary to change the code
      ll = this%control%lld
      nv = this%en%channels_ldos
      ldim = nb
      factor_z = 1.0d0
      legacy_mode = .false.
      if (present(legacy_real_axis)) legacy_mode = legacy_real_axis
      if (size(g_out, 1) /= nb .or. size(g_out, 2) /= nb .or. size(g_out, 3) < size(z_grid)) then
         error stop 'bgreen_complex: output and complex-energy batch dimensions are incompatible'
      end if

      ! Lightweight runtime shapes/logging (once per green instance call)
      ! call g_logger%info('DEBUG:bgreen shapes nb='//int2str(nb)//' ldim='//int2str(ldim), __FILE__, __LINE__)
      ! call g_logger%info('DEBUG:bgreen a_b dims=('//int2str(size(this%recursion%a_b,1))//','//int2str(size(this%recursion%a_b,2))//','//int2str(size(this%recursion%a_b,3))//','//int2str(size(this%recursion%a_b,4))//')', __FILE__, __LINE__)

      ! Unchanged code
      ll_t = ll
      llinf = ll!+ll/2!0!3*ll!/2!+30 !300
      cone = (1.0d0, 0.0d0)
      czero = (0.0d0, 0.0d0)
      im = (0.0d0, 1.0d0)
      zoff = (0.0d0, 0.01d0)*0.0d0
      coff = 0.000d0*(0.0d0, 1.0d0)
      one = (0.0d0, 0.0d0)

      ! Initializing Shift and Scaling for rigid band shift
      Dfac_mat = (1.0d0, 0.0d0)
      Cshi_mat = (0.0d0, 0.0d0)

      do i = 1, ldim
         m_tab(i) = i
      end do

      do i = 1, ldim
         one(i, i) = (1.0d0, 0.0d0)
      end do
      if (ldim >= 10) then
         a_diag = (a_inf(1, 1) + a_inf(10, 10))*0.5d0
         b_diag = (b_inf(1, 1) + b_inf(10, 10))*0.5d0
      else
         a_diag = a_inf(1, 1)
         b_diag = b_inf(1, 1)
      end if
      ! Standard semi-circle construction of Greens functions
      !write(12334,´(18f8.4)´) a_inf
      !write(12335,´(18f8.4)´) b_inf
      !write(12336,´(18f8.4)´) a_inf
      !write(12337,´(18f8.4)´) this%recursion%a_b(:,:,2,i_site)
      !write(12337,´(18f8.4)´) this%recursion%b2_b(:,:,2,i_site)
    !!!$omp parallel do default(shared) private(ei, Z, Ze, Q, B2z, i, ea, eb, det, zoff, l, ln, j, P , ipiv, info, work, lwork, W)
      !$omp parallel do default(shared) &
      !$omp          private(ei, Z, Ze, Q, B2z, i, etop, ebot, ea, eb, det, zoff, l, ln, j, P, ipiv, info, work, lwork, W, &
      !$omp                  found_nan, qHasNaN, qMax, wMax, b2Max, qval, wval, b2val, ii, jj, &
      !$omp                  recMaxA, recMaxB, recNaN, valrec, vv, maxQ, maxB2r, z_eval, eta_local)
      do ei = 1, size(z_grid)
         z_eval = z_grid(ei)
         if (legacy_mode) then
            Z = real(z_eval, rp)*one
            ze = cmplx(real(z_eval, rp), 0.0_rp, rp)
            eta_local = z_eval-ze
         else
            Z = z_eval*one
            ze = z_eval
            eta_local = cmplx(0.0_rp, 0.0_rp, rp)
         end if
         Q = (0.0d0, 0.0d0)
         B2z = 0.0d0
         found_nan = .false.
         if (this%control%sym_term) then
            do i = 1, ldim !  Orbital-independent
               etop = a_diag + 2.0d0*b_diag
               ebot = a_diag - 2.0d0*b_diag
               if (legacy_mode) then
                  ea = real(z_eval, rp) - etop
                  eb = real(z_eval, rp) - ebot
               else
                  ea = z_eval - etop
                  eb = z_eval - ebot
               end if
               det = ea*eb
               zoff = sqrt(det)
               Q(i, i) = (ze + eta_local - a_diag - zoff)*0.5d0
            end do
         else
            do i = 1, ldim  ! Orbital-dependent
               if (i == 1 .or. (ldim >= 10 .and. i == 10)) then  !(now the s-bands are broadened earlier)
                  etop = a_inf(i, i) - Cshi_mat(i, i) + 2*b_inf(i, i)*1.025d0/Dfac_mat(i, i)
                  ebot = a_inf(i, i) - Cshi_mat(i, i) - 2*b_inf(i, i)*1.025d0/Dfac_mat(i, i)
               else
                  etop = a_inf(i, i) - Cshi_mat(i, i) + 2*b_inf(i, i)!*1.01d0/Dfac_mat(i,i)
                  ebot = a_inf(i, i) - Cshi_mat(i, i) - 2*b_inf(i, i)!*1.01d0/Dfac_mat(i,i)
               end if
               if (legacy_mode) then
                  ea = real(z_eval, rp)/Dfac_mat(i, i) - etop - 1.0d0*Cshi_mat(i, i)
                  eb = real(z_eval, rp)/Dfac_mat(i, i) - ebot - 1.0d0*Cshi_mat(i, i)
               else
                  ea = z_eval/Dfac_mat(i, i) - etop - 1.0d0*Cshi_mat(i, i)
                  eb = z_eval/Dfac_mat(i, i) - ebot - 1.0d0*Cshi_mat(i, i)
               end if
               det = ea*eb
               zoff = sqrt(det)!*0.5d0
               zoff = zoff*(1.0)**ei
           !!! if(nv<=10) zoff=zoff*1.0d-3   ! Hack for Efermi evaluation
               ! Below for orbital dependent terminator
               if (legacy_mode) then
                  Q(i, i) = ((real(z_eval, rp) + eta_local)/Dfac_mat(i, i) - 1.0d0*Cshi_mat(i, i) - &
                     a_inf(i, i) - real(zoff) - aimag(zoff)*factor_z*im)*0.5d0
               else
                  Q(i, i) = (z_eval/Dfac_mat(i, i) - 1.0d0*Cshi_mat(i, i) - a_inf(i, i) - zoff)*0.5d0
               end if
            end do
         end if
         do l = llinf - 1, 1, -1
            ln = min(l, ll)
            do j = 1, ldim
               do i = 1, ldim
                  if (real(Z(i, j)) .ne. 0) then
                     P(i, j) = (Z(i, j) + eta_local)
                  else
                     P(i, j) = Z(i, j)
                  end if
                  if ((abs(real(Q(i, j))) .lt. 1.0d-12) .and. (abs(aimag(Q(i, j))) .lt. 1.0d-12)) then
                     Q(i, j) = (0.0d0, 0.0d0)
                  end if
                  Q(i, j) = P(i, j)/Dfac_mat(i, j) - Cshi_mat(i, j) - cone*this%recursion%a_b(i, j, ln, i_site) - Q(i, j)
                  B2z(i, j) = cone*this%recursion%b2_b(i, j, ln, i_site)
               end do
            end do
            ! Check Q for NaNs before LU factorization
            qMax = 0.0_rp
            qHasNaN = .false.
            do jj = 1, ldim
               do ii = 1, ldim
                  qval = abs(Q(ii, jj))
                  if (qval == qval) then
                     if (qval > qMax) qMax = qval
                  else
                     qHasNaN = .true.
                  end if
               end do
            end do
            if (qHasNaN) then
               call g_logger%error('Pre-LU Q contains NaN (site='//int2str(i_site)//' ie='//int2str(ei)//' Q_max='//real2str(qMax)//')', __FILE__, __LINE__)
               !$omp critical
               if (.not. dump_done) then
                  dump_done = .true.
                  dump_fname = 'debug_bgreen_preLU_site'//int2str(i_site)//'_ie'//int2str(ei)//'.txt'
                  dump_unit = 901
                  open(dump_unit, file=trim(dump_fname), status='replace', action='write', iostat=dump_iostat)
                  if (dump_iostat == 0) then
                     write(dump_unit, '(A)') 'PRE-LU DEBUG DUMP: site='//int2str(i_site)//' ie='//int2str(ei)
                        write(dump_unit, '(A)') 'SHAPES: nb='//int2str(nb)//' ldim='//int2str(ldim)//' ll='//int2str(ll)
                        write(dump_unit, '(A)') 'Small P (real imag) first 4x4:'
                        do jj = 1, min(4,ldim)
                           do ii = 1, min(4,ldim)
                              write(dump_unit, '(I4,1X,I4,1X,ES24.16,1X,ES24.16)') ii, jj, real(P(ii,jj)), aimag(P(ii,jj))
                           end do
                        end do
                        write(dump_unit, '(A)') 'Small Dfac_mat (real imag) first 4x4:'
                        do jj = 1, min(4,ldim)
                           do ii = 1, min(4,ldim)
                              write(dump_unit, '(I4,1X,I4,1X,ES24.16,1X,ES24.16)') ii, jj, real(Dfac_mat(ii,jj)), aimag(Dfac_mat(ii,jj))
                           end do
                        end do
                        write(dump_unit, '(A)') 'Small Cshi_mat (real imag) first 4x4:'
                        do jj = 1, min(4,ldim)
                           do ii = 1, min(4,ldim)
                              write(dump_unit, '(I4,1X,I4,1X,ES24.16,1X,ES24.16)') ii, jj, real(Cshi_mat(ii,jj)), aimag(Cshi_mat(ii,jj))
                           end do
                        end do
                        write(dump_unit, '(A)') 'Sample recursion contributions (ln=1..min(4,ll)):'
                        do ln = 1, min(4,ll)
                           write(dump_unit, '(A,I0)') 'ln=', ln
                           do jj = 1, min(4,ldim)
                              do ii = 1, min(4,ldim)
                                 write(dump_unit, '(I4,1X,I4,1X,ES24.16,1X,ES24.16,1X,ES24.16,1X,ES24.16)') ii, jj, real(this%recursion%a_b(ii,jj,ln,i_site)), aimag(this%recursion%a_b(ii,jj,ln,i_site)), real(this%recursion%b2_b(ii,jj,ln,i_site)), aimag(this%recursion%b2_b(ii,jj,ln,i_site))
                              end do
                           end do
                        end do
                        write(dump_unit, '(A)') '--- continuing with full Q/B2z dump below ---'
                     write(dump_unit, '(A)') 'Q matrix (real imag):'
                     do jj = 1, ldim
                        do ii = 1, ldim
                           write(dump_unit, '(I4,1X,I4,1X,ES24.16,1X,ES24.16)') ii, jj, real(Q(ii,jj)), aimag(Q(ii,jj))
                        end do
                     end do
                     write(dump_unit, '(A)') 'B2z matrix (real imag):'
                     do jj = 1, ldim
                        do ii = 1, ldim
                           write(dump_unit, '(I4,1X,I4,1X,ES24.16,1X,ES24.16)') ii, jj, real(B2z(ii,jj)), aimag(B2z(ii,jj))
                        end do
                     end do
                     close(dump_unit)
                  end if
               end if
               !$omp end critical
               found_nan = .true.
               exit
            end if

            ! LU factorization
            call zgetrf(ldim, ldim, Q, ldim, ipiv, info)
            if (info /= 0) then
               ! Collect diagnostics and attempt tiny regularization then retry once
               maxQ = 0.0_rp
               maxB2r = 0.0_rp
               do jj = 1, ldim
                  do ii = 1, ldim
                     if (abs(Q(ii, jj)) > maxQ) maxQ = abs(Q(ii, jj))
                     if (abs(B2z(ii, jj)) > maxB2r) maxB2r = abs(B2z(ii, jj))
                  end do
               end do
               recMaxA = 0.0_rp
               recMaxB = 0.0_rp
               recNaN = .false.
               do ln = 1, ll
                  do ii = 1, ldim
                     do jj = 1, ldim
                        valrec = abs(this%recursion%a_b(ii, jj, ln, i_site))
                        if (valrec > recMaxA) recMaxA = valrec
                        if (IsNaN(real(this%recursion%a_b(ii, jj, ln, i_site))) .or. IsNaN(aimag(this%recursion%a_b(ii, jj, ln, i_site)))) recNaN = .true.
                        valrec = abs(this%recursion%b2_b(ii, jj, ln, i_site))
                        if (valrec > recMaxB) recMaxB = valrec
                        if (IsNaN(real(this%recursion%b2_b(ii, jj, ln, i_site))) .or. IsNaN(aimag(this%recursion%b2_b(ii, jj, ln, i_site)))) recNaN = .true.
                     end do
                  end do
               end do
               call g_logger%warning('LU factorization failed (zgetrf info='//fmt('I0', info)//') site='//int2str(i_site)//' ie='//int2str(ei)//' maxQ='//real2str(maxQ)//' maxB2='//real2str(maxB2r)//' recMaxA='//real2str(recMaxA)//' recMaxB='//real2str(recMaxB), __FILE__, __LINE__)
               ! Try tiny regularization
               Q = Q + (1.0d-12, 0.0d0)*one
               call zgetrf(ldim, ldim, Q, ldim, ipiv, info)
               if (info /= 0) then
                  call g_logger%error('zgetrf retry failed (info='//fmt('I0', info)//') skipping energy slice site='//int2str(i_site)//' ie='//int2str(ei), __FILE__, __LINE__)
                  found_nan = .true.
                  exit
               end if
            end if
            ! Inverse of Q from LU factorization
            lwork = ldim*ldim
            call zgetri(ldim, Q, ldim, ipiv, work, lwork, info)
            ! Check Q for NaNs after inversion
            qMax = 0.0_rp
            qHasNaN = .false.
            do jj = 1, ldim
               do ii = 1, ldim
                  qval = abs(Q(ii, jj))
                  if (qval == qval) then
                     if (qval > qMax) qMax = qval
                  else
                     qHasNaN = .true.
                  end if
               end do
            end do
            if (qHasNaN) then
               call g_logger%error('Post-inverse Q contains NaN (site='//int2str(i_site)//' ie='//int2str(ei)//' Q_max='//real2str(qMax)//' zgetri_info='//fmt('I0', info)//')', __FILE__, __LINE__)
               !$omp critical
               if (.not. dump_done) then
                  dump_done = .true.
                  dump_fname = 'debug_bgreen_postInv_site'//int2str(i_site)//'_ie'//int2str(ei)//'.txt'
                  dump_unit = 902
                  open(dump_unit, file=trim(dump_fname), status='replace', action='write', iostat=dump_iostat)
                  if (dump_iostat == 0) then
                     write(dump_unit, '(A)') 'POST-INV DEBUG DUMP: site='//int2str(i_site)//' ie='//int2str(ei)
                     write(dump_unit, '(A)') 'Q matrix (real imag):'
                     do jj = 1, ldim
                        do ii = 1, ldim
                           write(dump_unit, '(I4,1X,I4,1X,ES24.16,1X,ES24.16)') ii, jj, real(Q(ii,jj)), aimag(Q(ii,jj))
                        end do
                     end do
                     write(dump_unit, '(A)') 'B2z matrix (real imag):'
                     do jj = 1, ldim
                        do ii = 1, ldim
                           write(dump_unit, '(I4,1X,I4,1X,ES24.16,1X,ES24.16)') ii, jj, real(B2z(ii,jj)), aimag(B2z(ii,jj))
                        end do
                     end do
                     close(dump_unit)
                  end if
               end if
               !$omp end critical
               found_nan = .true.
               exit
            end if
            if (info /= 0) then
               call g_logger%warning('zgetri (inverse) failed (info='//fmt('I0', info)//') site='//int2str(i_site)//' ie='//int2str(ei), __FILE__, __LINE__)
               ! Attempt small regularization and retry inverse
               Q = Q + (1.0d-12, 0.0d0)*one
               call zgetrf(ldim, ldim, Q, ldim, ipiv, info)
               if (info == 0) then
                  call zgetri(ldim, Q, ldim, ipiv, work, lwork, info)
               end if
               if (info /= 0) then
                  call g_logger%error('zgetri retry failed (info='//fmt('I0', info)//') skipping energy slice site='//int2str(i_site)//' ie='//int2str(ei), __FILE__, __LINE__)
                  found_nan = .true.
                  exit
               end if
            end if
            call zgemm('n', 'n', ldim, ldim, ldim, cone, Q, ldim, B2z, ldim, czero, W, ldim)
            call zgemm('c', 'n', ldim, ldim, ldim, cone, B2z, ldim, W, ldim, czero, Q, ldim)
         end do
         !
         if (.not. found_nan) then
            do j = 1, ldim
               do i = 1, ldim
                  !bdos(ei,m_tab(i),m_tab(j),n)=bdos(ei,m_tab(i),m_tab(j),n) + abs(-aimag(Q(i,j))/3.14159265359d0)/Dfac_mat(i,j)
                  if (legacy_mode .and. aimag(eta_local) .eq. 0) then
                     g_out(m_tab(i), m_tab(j), ei) = g_out(m_tab(i), m_tab(j), ei) + &
                                                     real(Q(i, j)/Dfac_mat(i, j)) + aimag(Q(i, j)/Dfac_mat(i, j))*(1.0d0)**ei*(0.0d0, 1.0d0)
                  else
                     g_out(m_tab(i), m_tab(j), ei) = g_out(m_tab(i), m_tab(j), ei) + (Q(i, j)/Dfac_mat(i, j))
                  end if
               end do
            end do
         end if
         ! Quick NaN check for this energy slice (log first occurrence)
         found_nan = .false.
         do j = 1, ldim
            do i = 1, ldim
               vv = g_out(m_tab(i), m_tab(j), ei)
               if (real(vv) /= real(vv) .or. aimag(vv) /= aimag(vv)) then
                  call g_logger%warning('NaN detected in g_out at site '//fmt('I0', i_site)//' ie='//fmt('I0', ei), __FILE__, __LINE__)
                  found_nan = .true.
                        ! Collect recursion coefficient summaries for this site
                        recMaxA = 0.0_rp
                        recMaxB = 0.0_rp
                        recNaN = .false.
                        do ln = 1, ll
                           do ii = 1, ldim
                              do jj = 1, ldim
                                 valrec = abs(this%recursion%a_b(ii, jj, ln, i_site))
                                 if (valrec > recMaxA) recMaxA = valrec
                                 if (IsNaN(real(this%recursion%a_b(ii, jj, ln, i_site))) .or. IsNaN(aimag(this%recursion%a_b(ii, jj, ln, i_site)))) recNaN = .true.
                                 valrec = abs(this%recursion%b2_b(ii, jj, ln, i_site))
                                 if (valrec > recMaxB) recMaxB = valrec
                                 if (IsNaN(real(this%recursion%b2_b(ii, jj, ln, i_site))) .or. IsNaN(aimag(this%recursion%b2_b(ii, jj, ln, i_site)))) recNaN = .true.
                              end do
                           end do
                        end do
                        ! call g_logger%info('DEBUG:bgreen recursion site='//int2str(i_site)//' ie='//int2str(ei)//' maxA='//real2str(recMaxA)//' maxB='//real2str(recMaxB)//' recNaN='//log2str(recNaN), __FILE__, __LINE__)
                        ! ! Log diagonal of a_inf/b_inf and Dfac_mat for quick inspection
                        ! do ii = 1, ldim
                        !    call g_logger%info('DEBUG:bgreen a_inf('//int2str(ii)//')='//real2str(a_inf(ii,ii))//' b_inf='//real2str(b_inf(ii,ii))//' Dfac='//real2str(real(Dfac_mat(ii,ii))), __FILE__, __LINE__)
                        ! end do
                        ! Additional diagnostics: check Q, W and B2z for NaNs/large values
                        qMax = 0.0_rp
                        qHasNaN = .false.
                        wMax = 0.0_rp
                        b2Max = 0.0_rp
                        do jj = 1, ldim
                           do ii = 1, ldim
                              qval = abs(Q(ii,jj))
                              if (qval == qval) then
                                 if (qval > qMax) qMax = qval
                              else
                                 qHasNaN = .true.
                              end if
                              wval = abs(W(ii,jj))
                              if (wval == wval) then
                                 if (wval > wMax) wMax = wval
                              end if
                              b2val = abs(B2z(ii,jj))
                              if (b2val == b2val) then
                                 if (b2val > b2Max) b2Max = b2val
                              end if
                           end do
                        end do
                        ! call g_logger%info('DEBUG:bgreen Q_max='//real2str(qMax)//' Q_hasNaN='//log2str(qHasNaN)//' W_max='//real2str(wMax)//' B2z_max='//real2str(b2Max), __FILE__, __LINE__)
                       !$omp critical
                       if (.not. dump_done) then
                          dump_done = .true.
                          dump_fname = 'debug_bgreen_site'//int2str(i_site)//'_ie'//int2str(ei)//'.txt'
                          dump_unit = 900
                          open(dump_unit, file=trim(dump_fname), status='replace', action='write', iostat=dump_iostat)
                          if (dump_iostat == 0) then
                             write(dump_unit, '(A)') 'DEBUG DUMP: site='//int2str(i_site)//' ie='//int2str(ei)
                             write(dump_unit, '(A)') 'Q matrix (real imag):'
                             do jj = 1, ldim
                                do ii = 1, ldim
                                   write(dump_unit, '(I4,1X,I4,1X,ES24.16,1X,ES24.16)') ii, jj, real(Q(ii,jj)), aimag(Q(ii,jj))
                                end do
                             end do
                             write(dump_unit, '(A)') 'W matrix (real imag):'
                             do jj = 1, ldim
                                do ii = 1, ldim
                                   write(dump_unit, '(I4,1X,I4,1X,ES24.16,1X,ES24.16)') ii, jj, real(W(ii,jj)), aimag(W(ii,jj))
                                end do
                             end do
                             write(dump_unit, '(A)') 'B2z matrix (real imag):'
                             do jj = 1, ldim
                                do ii = 1, ldim
                                   write(dump_unit, '(I4,1X,I4,1X,ES24.16,1X,ES24.16)') ii, jj, real(B2z(ii,jj)), aimag(B2z(ii,jj))
                                end do
                             end do
                             write(dump_unit, '(A)') 'recursion a_b (ln, i, j, real imag):'
                             do ln = 1, ll
                                do jj = 1, ldim
                                   do ii = 1, ldim
                                      write(dump_unit, '(I4,1X,I4,1X,I4,1X,ES24.16,1X,ES24.16)') ln, ii, jj, real(this%recursion%a_b(ii,jj,ln,i_site)), aimag(this%recursion%a_b(ii,jj,ln,i_site))
                                   end do
                                end do
                             end do
                             write(dump_unit, '(A)') 'recursion b2_b (ln, i, j, real imag):'
                             do ln = 1, ll
                                do jj = 1, ldim
                                   do ii = 1, ldim
                                      write(dump_unit, '(I4,1X,I4,1X,I4,1X,ES24.16,1X,ES24.16)') ln, ii, jj, real(this%recursion%b2_b(ii,jj,ln,i_site)), aimag(this%recursion%b2_b(ii,jj,ln,i_site))
                                   end do
                                end do
                             end do
                             close(dump_unit)
                          end if
                       end if
                       !$omp end critical
                        exit
               end if
            end do
            if (found_nan) exit
         end do
         !write(12333,´(19f8.4)´) e(ei), (real(g_out(i,i,ei)),i=1,18)
      end do
      !$omp end parallel do

   end subroutine bgreen_complex

end submodule green_lanczos
