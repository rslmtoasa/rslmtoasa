module cfd
   !
   use Parameters
   !
   implicit none
   !
   integer :: i_cons
   integer :: stdout = 6
   integer, parameter :: ncomp = 3
   real(dblprec), parameter :: b2t = 235298.924212429_dblprec
   !
   real(dblprec) :: constrained_mom_err
   real(dblprec) :: lambda = 10
   real(dblprec) :: lambda_t
   real(dblprec) :: induced_mom_thresh = 0.5_dblprec
   real(dblprec) :: cfd_prefac = 1.0_dblprec
   real(dblprec) :: bfield_beta = 1.0_dblprec
   integer :: b_constr_iter
   !
   ! arrays for PID
   real(dblprec), dimension(:, :), allocatable :: d_delta
   real(dblprec), dimension(:, :), allocatable :: s_delta
   real(dblprec), dimension(:, :), allocatable :: dd_delta
   !
   private
   !
   public :: initialize_cfd, constrain
   !
contains
   !
   subroutine initialize_cfd(ndim, flag, i_cons_in, code_prefac)
      !
      implicit none
      !
      integer, intent(in) :: ndim
      integer, intent(in) :: flag
      integer, intent(in), optional :: i_cons_in
      integer, intent(in), optional :: code_prefac
      !
      if (flag > 0) then
         allocate (d_delta(ncomp, ndim))
         allocate (s_delta(ncomp, ndim))
         allocate (dd_delta(ncomp, ndim))
         b_constr_iter = 0
         lambda_t = lambda
         if (present(i_cons_in)) i_cons = i_cons_in
         if (present(code_prefac)) cfd_prefac = code_prefac
      else
         deallocate (d_delta)
         deallocate (s_delta)
         deallocate (dd_delta)
      end if
      !
   end subroutine initialize_cfd
   !
   subroutine constrain(mom_in, mom_ref, bfield, ndim)
      !
      use abspinlib
      !
      implicit none
      !
      ! arguments
      integer, intent(in) :: ndim
      real(dblprec), intent(in) :: mom_in(ncomp, ndim)
      real(dblprec), intent(in) :: mom_ref(ncomp, ndim)
      real(dblprec), intent(inout) :: bfield(ncomp, ndim)
      !
      !
      integer :: iidim, na
      real(dblprec), dimension(:, :), allocatable :: mom_tmp
      real(dblprec), dimension(3) :: e_i, e_out, c_in
      real(dblprec), dimension(3) :: m_delta, bfield_new
      real(dblprec) :: etcon, ma, mnorm

      !!!   cfd_prefac=b2t*omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3)

      allocate (mom_tmp(ncomp, ndim))

      do na = 1, ndim
         if (i_cons == 2) then
            ! Lagrange multiplier without orthogonalization
            !
            mom_tmp(:, na) = mom_in(:, na)/norm2(mom_in(:, na)) - mom_ref(:, na)
            ! Gramm-Schmidt step
            etcon = etcon + lambda_t*sum(mom_tmp(:, na)**2)
         else if (i_cons == 3) then
            ! Lagrange multiplier with orthogonalization (b _|_ m)
            !
            mom_tmp(:, na) = mom_in(:, na)/norm2(mom_in(:, na)) - mom_ref(:, na)
            ! Gramm-Schmidt step
            mom_tmp(:, na) = mom_tmp(:, na) - sum(mom_tmp(:, na)*mom_ref(:, na))*mom_ref(:, na)
            etcon = etcon + lambda_t*sum(mom_tmp(:, na)**2)
         else if (i_cons == 4) then
            ! i_cons = 4 means that we try to use a PID regulator
            !
            !
            !
            write (stdout, '(4x,a)') ' | AMN-PID noncolinear constraints '
            !
            ! Check moment magnitude
            mnorm = sqrt(mom_in(1, na)**2 + mom_in(2, na)**2 + mom_in(3, na)**2)
            write (stdout, '(4x,a,i4)') " | - atom: ", na
            if (mnorm .lt. induced_mom_thresh) then
               write (stdout, '(2x,a,i4,a,f10.4)') ' | Local magnetization for atom ', na, ' is less than threshold', ma
               m_delta = 0.0_dblprec
            else
               c_in = bfield(:, na)
               ! Direction only
               e_out = mom_in(:, na)/ma
               e_i = mom_ref(:, na)/norm2(mom_ref(:, na))
               !! Direction and magnitude
               !e_out=mom_in(:,na)
               !e_i =mom_ref(:,na)
               ! P I D
               ! Full direction (untested)
               !m_delta=(e_i-e_out)
               ! Perp direction (works for bcc fe)
               m_delta = -(mom_in(:, na) - sum(mom_in(:, na)*e_i)*e_i)
            end if
            !
            ! Reducing the effect for first iteration (ie when d_delta=0)
            if (norm2(d_delta(:, na)) < 1e-15) m_delta = 0.1_dblprec*m_delta
            ! e) m_delta=-lambda_t*(mom_in(:,na)-sum(mom_in(:,na)*e_i)*e_i)*10.0_dblprec
            ! others:lambda_t=0.1
            !gs
            !m_delta=-(e_out-norm2(e_out*e_i)*e_i)
            !
            ! Check to don't mix first iteration
            if (norm2(d_delta(:, na)) > 1e-15) dd_delta(:, na) = m_delta - d_delta(:, na)
            !
            write (stdout, '(4x,a,i4,3f15.8)') " | Output moments     for atom ", na, mom_in(1:3, na)
            write (stdout, '(4x,a,i4,3f15.8)') " | Input direction    for atom ", na, mom_ref(:, na)
            write (stdout, '(4x,a,i4,3f15.8)') " | Outut direction    for atom ", na, e_out
            write (stdout, '(4x,a,i4,3f15.8)') " | Input field        for atom ", na, bfield(1:3, na)
            !
            ! Check to don't mix first iteration
            if (norm2(d_delta(:, na)) > 1e-15) s_delta(:, na) = s_delta(:, na) + m_delta

            !bfield(:,na)=lambda_t*(1.20_dblprec*m_delta+0.35_dblprec*s_delta(:,na)+0.10_dblprec*dd_delta(:,na))
            bfield(:, na) = lambda*(1.30_dblprec*m_delta + 0.35_dblprec*s_delta(:, na) - 0.10_dblprec*dd_delta(:, na))   !<-- use this for atoms

            !bfield_pts(:,ir)=lambda_t*(1.00_dblprec*m_delta+0.12_dblprec*s_delta_pts(:,ir)+0.10_dblprec*dd_delta(:,na))   !ok for grids

            ! Calculate Zeeman-like constraining energy cost
            etcon = etcon + sum(bfield(1:3, na)*mom_in(:, na))
            d_delta(:, na) = m_delta
            !
            write (stdout, '(4x,a,i4,3f15.8)') " | P  contribution    for atom ", na, d_delta(:, na)
            write (stdout, '(4x,a,i4,3f15.8)') " | I  contribution    for atom ", na, s_delta(:, na)
            write (stdout, '(4x,a,i4,3f15.8)') " | D  contribution    for atom ", na, dd_delta(:, na)
            write (stdout, '(4x,a,i4,4f15.8)') " | Constraining field for atom ", na, bfield(1:3, na)
            write (stdout, '(4x,a,i4,3f15.4)') " | Constraining field for atom (t)", na, cfd_prefac*bfield(1:3, na)

         else if (i_cons == 5) then
            ! i_cons = 5 means that we try the Ma-Dudarev approach
            ! which is very analogous to the normal Lagrange approach
            !
            !
            write (stdout, '(2x,a)') ' Ma-Dudarev constraints '
            constrained_mom_err = 0.0_dblprec
            !
            ! Check moment magnitude
            ma = dsqrt(mom_in(1, na)**2 + mom_in(2, na)**2 + mom_in(3, na)**2)
            write (stdout, '(4x,a,i4)') " | - Atom: ", na
            !
            if (ma .lt. induced_mom_thresh) then
               write (stdout, '(2x,a,i4,a,f10.4)') ' | Local magnetization for atom ', na, ' is less than threshold', ma
               bfield(:, na) = 0.0_dblprec
            else
               write (stdout, '(4x,a,i4,3f15.8)') " | Output moments     for atom ", na, mom_in(1:3, na)
               write (stdout, '(4x,a,i4,3f15.8)') " | Input direction    for atom ", na, mom_ref(:, na)
               write (stdout, '(4x,a,i4,3f15.8)') " | Input field        for atom ", na, bfield(1:3, na)
               !
               e_out = mom_in(:, na)
               e_i = mom_ref(:, na)/norm2(mom_ref(:, na))*ma
               !e_dot=e_out(1)*e_i (1)+e_out(2)*e_i (2)+e_out(3)*e_i (3)
               ! new constraining field
               !bfield_new=-10.0_dblprec*lambda_t*(e_out-e_i)
               bfield_new = lambda_t*(e_out - e_i)
               ! gram-schmidt orthogonalization ( b _|_ m)
               !bfield_new=bfield_new-sum(e_out*bfield_new)*e_out/ma**2
               bfield_new = bfield_new - sum(e_i*bfield_new)*e_i/ma**2
               ! mixing field
               bfield(:, na) = (1.0_dblprec - bfield_beta)*bfield(:, na) - bfield_beta*bfield_new
               !
               ! calculate zeeman-like constraining energy cost
               etcon = etcon + lambda_t*(sqrt(sum(mom_in(1:3, na)*mom_in(:, na))) - sum(mom_in(1:3, na)*e_i/ma))
               constrained_mom_err = constrained_mom_err + sum(e_out - e_i)**2
               ! write (stdout,'(4x,a,i4,3f15.8)' ) " | new field          for atom ",na,bfield_new
               write (stdout, '(4x,a,i4,4f15.8)') " | Output field       for atom ", na, bfield(1:3, na), &
                  sum(bfield(:, na)*e_i)
               write (stdout, '(4x,a,i4,4f15.8)') " | Output field (t)   for atom ", na, &
                  bfield(1:3, na)*cfd_prefac, &
                  sum(bfield(:, na)*e_i)
               write (stdout, '(4x,a,i4,3f15.8)') " | Output direction    for atom ", na, e_out/ma
            end if
         end if
      end do ! na
      !
      constrained_mom_err = sqrt(constrained_mom_err)/ndim

      b_constr_iter = b_constr_iter + 1
      !if(i_cons==5) lambda_t=min(lambda_t+4_dblprec,100.0_dblprec)
      if (i_cons == 5) lambda_t = min(lambda_t + 1.0_dblprec - lambda_t/lambda, lambda)
      !if(i_cons==3) lambda_t=min(lambda_t+1.0_dblprec-lambda_t/lambda,lambda)
      ! works for moderate lambdas
      !if(i_cons==3) then
      !   if (b_constr_iter<=30.0_dblprec) then
      !      lambda_t=min(lambda_t*(2.0_dblprec-b_constr_iter/30.0_dblprec),lambda)
      !   else
      !      lambda_t=lambda
      !   end if
      !end if
      !
      !
      ! scale up lambda in case of Lagrangian formulation
      if (etcon < 1.0d-2) lambda_t = min(1.2_dblprec*lambda_t, 1.0e4_dblprec)
      if (i_cons == 3) lambda_t = min(lambda_t*(2.0_dblprec - lambda_t/lambda), lambda)
      !if(i_cons==5) lambda_t=min(lambda_t*(2.0_dblprec-lambda_t/lambda),lambda)
      if (i_cons == 5) lambda_t = min(lambda_t + 2.0_dblprec, lambda)
      ! if(i_cons==4) lambda_t=min(lambda_t+1.0_dblprec,25.0_dblprec)
      !if(i_cons==4) lambda_t=lambda_t+1.0_dblprec
      !if(i_cons==5) lambda_t=lambda_t+lambda_t*min(0.1_dblprec,etcon**2)
      !if(i_cons==5) lambda_t=lambda_t*(1.0_dblprec+0.5_dblprec*(etcon))
      !if(i_cons==5) lambda_t=lambda_t*(1.0_dblprec+2.0_dblprec*min(constrained_mom_err,0.1_dblprec))
      if (i_cons == 5) write (stdout, '(4x,a,f12.4a,g10.2)') " | New lambda_t: ", lambda_t, "     error: ", constrained_mom_err
      write (stdout, '(4x,a)') " | -  "
      deallocate (mom_tmp)
      return
   end subroutine constrain

end module cfd
