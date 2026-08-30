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
   logical :: diagnostics_enabled = .false.

   !> Machine-readable diagnostics for one controller update.  All vectors are
   !> in the same spin frame as the arguments to constrain(): the ordinary
   !> Cartesian frame for real-space/explicit-texture calculations and the
   !> primitive rotating frame for gbt_single_q.
   type, public :: constraint_diagnostics
      integer :: iteration = 0
      integer :: nsite = 0
      real(dblprec) :: convergence_metric = 0.0_dblprec
      real(dblprec) :: max_angular_error = 0.0_dblprec
      !> Instantaneous residual merit assembled by the selected controller.
      !> This is a numerical penalty/merit quantity, not a DFT total energy.
      real(dblprec) :: controller_penalty_energy = 0.0_dblprec
      !> Expectation value of the constraining Pauli term in the measured state,
      !> sum_i B_con(i) dot m(i).  This is a band/eigenvalue bookkeeping
      !> diagnostic only; it is not added to the physical DFT energy.
      real(dblprec) :: field_coupling_energy = 0.0_dblprec
      !> The controller does not implement the variational cDFT Lagrange
      !> functional.  Keep this explicit so callers cannot infer one from
      !> controller_penalty_energy.
      logical :: lagrange_functional_represented = .false.
      real(dblprec), allocatable :: target(:, :)
      real(dblprec), allocatable :: actual(:, :)
      real(dblprec), allocatable :: moment_magnitude(:)
      real(dblprec), allocatable :: angular_error(:)
      real(dblprec), allocatable :: transverse_residual(:, :)
      real(dblprec), allocatable :: bfield(:, :)
      real(dblprec), allocatable :: bfield_magnitude(:)
      real(dblprec), allocatable :: bfield_longitudinal(:)
      real(dblprec), allocatable :: torque(:, :)
   contains
      procedure :: clear => clear_constraint_diagnostics
   end type constraint_diagnostics
   !
   ! arrays for PID
   real(dblprec), dimension(:, :), allocatable :: d_delta
   real(dblprec), dimension(:, :), allocatable :: s_delta
   real(dblprec), dimension(:, :), allocatable :: dd_delta
   !
   private
   !
   public :: initialize_cfd, constrain, set_constraint_diagnostics, set_constraint_controller_parameters
   !
contains
   !
   subroutine initialize_cfd(ndim, flag, i_cons_in, code_prefac, diagnostics_in)
      !
      implicit none
      !
      integer, intent(in) :: ndim
      integer, intent(in) :: flag
      integer, intent(in), optional :: i_cons_in
      integer, intent(in), optional :: code_prefac
      logical, intent(in), optional :: diagnostics_in
      !
      if (flag > 0) then
         if (present(diagnostics_in)) diagnostics_enabled = diagnostics_in
         if (allocated(d_delta)) then
            if (size(d_delta, 2) == ndim) then
               if (present(i_cons_in)) i_cons = i_cons_in
               if (present(code_prefac)) cfd_prefac = code_prefac
               return
            end if
         end if
         if (allocated(d_delta)) deallocate(d_delta)
         if (allocated(s_delta)) deallocate(s_delta)
         if (allocated(dd_delta)) deallocate(dd_delta)
         allocate (d_delta(ncomp, ndim))
         allocate (s_delta(ncomp, ndim))
         allocate (dd_delta(ncomp, ndim))
         b_constr_iter = 0
         lambda_t = lambda
         ! Zero-initialize PID arrays to avoid uninitialized values
         d_delta = 0.0_dblprec
         s_delta = 0.0_dblprec
         dd_delta = 0.0_dblprec
         if (present(i_cons_in)) i_cons = i_cons_in
         if (present(code_prefac)) cfd_prefac = code_prefac
      else
         diagnostics_enabled = .false.
         if (allocated(d_delta)) deallocate (d_delta)
         if (allocated(s_delta)) deallocate (s_delta)
         if (allocated(dd_delta)) deallocate (dd_delta)
      end if
      !
   end subroutine initialize_cfd

   !> Enable or disable the verbose legacy controller trace.  Structured
   !> diagnostics returned by constrain() are independent of this switch.
   subroutine set_constraint_diagnostics(enabled)
      logical, intent(in) :: enabled

      diagnostics_enabled = enabled
   end subroutine set_constraint_diagnostics

   !> Set the controller nominal penalty/gain strength.  This is primarily
   !> an audit and experiment hook; it never changes the physical-energy
   !> bookkeeping.  A subsequent initialize_cfd call starts lambda_t from it.
   subroutine set_constraint_controller_parameters(lambda_in)
      real(dblprec), intent(in) :: lambda_in

      if (lambda_in <= 0.0_dblprec) error stop 'constraint controller strength must be positive'
      lambda = lambda_in
      if (allocated(d_delta)) lambda_t = lambda
   end subroutine set_constraint_controller_parameters
   !
   subroutine constrain(mom_in, mom_ref, bfield, ndim, etcon_out, update_lambda, diagnostics_out)
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
      real(dblprec), intent(out), optional :: etcon_out
      logical, intent(in), optional :: update_lambda
      type(constraint_diagnostics), intent(out), optional :: diagnostics_out
      !
      !
      integer :: na
      real(dblprec), dimension(:, :), allocatable :: mom_tmp
      real(dblprec), dimension(3) :: e_i, e_out, c_in
      real(dblprec), dimension(3) :: m_delta, bfield_new, transverse, torque
      real(dblprec) :: etcon, ma, mnorm, angle, cosine, metric_sum, field_coupling
      logical :: update_lambda_

      !!!   cfd_prefac=b2t*omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3)

      allocate (mom_tmp(ncomp, ndim))
      etcon = 0.0_dblprec
      field_coupling = 0.0_dblprec
      constrained_mom_err = 0.0_dblprec
      mom_tmp = 0.0_dblprec
      update_lambda_ = .true.
      if (present(update_lambda)) update_lambda_ = update_lambda

      if (present(diagnostics_out)) then
         call diagnostics_out%clear()
         diagnostics_out%iteration = b_constr_iter + 1
         diagnostics_out%nsite = ndim
         allocate(diagnostics_out%target(3, ndim), diagnostics_out%actual(3, ndim), &
                  diagnostics_out%moment_magnitude(ndim), diagnostics_out%angular_error(ndim), &
                  diagnostics_out%transverse_residual(3, ndim), diagnostics_out%bfield(3, ndim), &
                  diagnostics_out%bfield_magnitude(ndim), diagnostics_out%bfield_longitudinal(ndim), &
                  diagnostics_out%torque(3, ndim))
         diagnostics_out%target = 0.0_dblprec
         diagnostics_out%actual = 0.0_dblprec
         diagnostics_out%moment_magnitude = 0.0_dblprec
         diagnostics_out%angular_error = 0.0_dblprec
         diagnostics_out%transverse_residual = 0.0_dblprec
         diagnostics_out%bfield = 0.0_dblprec
         diagnostics_out%bfield_magnitude = 0.0_dblprec
         diagnostics_out%bfield_longitudinal = 0.0_dblprec
         diagnostics_out%torque = 0.0_dblprec
      end if
      metric_sum = 0.0_dblprec

      do na = 1, ndim
         if (i_cons == 2) then
            ! Historical multiplier-style penalty controller without
            ! orthogonalization.
            !
            ma = norm2(mom_in(:, na))
            if (ma <= tiny(1.0_dblprec) .or. norm2(mom_ref(:, na)) <= tiny(1.0_dblprec)) then
               e_out = 0.0_dblprec
               if (norm2(mom_ref(:, na)) > tiny(1.0_dblprec)) then
                  e_i = mom_ref(:, na)/norm2(mom_ref(:, na))
               else
                  e_i = [0.0_dblprec, 0.0_dblprec, 1.0_dblprec]
               end if
               mom_tmp(:, na) = 0.0_dblprec
            else
               e_out = mom_in(:, na)/ma
               e_i = mom_ref(:, na)/norm2(mom_ref(:, na))
               mom_tmp(:, na) = e_out - e_i
               ! Gramm-Schmidt step
               etcon = etcon + 0.5_dblprec*lambda_t*sum(mom_tmp(:, na)**2)
               bfield(:, na) = bfield(:, na) - lambda_t*mom_tmp(:, na)
            end if
         else if (i_cons == 3) then
            ! Historical multiplier-style penalty controller with
            ! orthogonalization (b _|_ m).
            !
            ma = norm2(mom_in(:, na))
            if (ma <= tiny(1.0_dblprec) .or. norm2(mom_ref(:, na)) <= tiny(1.0_dblprec)) then
               e_out = 0.0_dblprec
               if (norm2(mom_ref(:, na)) > tiny(1.0_dblprec)) then
                  e_i = mom_ref(:, na)/norm2(mom_ref(:, na))
               else
                  e_i = [0.0_dblprec, 0.0_dblprec, 1.0_dblprec]
               end if
               mom_tmp(:, na) = 0.0_dblprec
            else
               e_out = mom_in(:, na)/ma
               e_i = mom_ref(:, na)/norm2(mom_ref(:, na))
               mom_tmp(:, na) = e_out - e_i
               ! Gramm-Schmidt step
               mom_tmp(:, na) = mom_tmp(:, na) - sum(mom_tmp(:, na)*e_i)*e_i
               etcon = etcon + 0.5_dblprec*lambda_t*sum(mom_tmp(:, na)**2)
               bfield(:, na) = bfield(:, na) - lambda_t*mom_tmp(:, na)
               ! Method 3 is a transverse Lagrange multiplier.  Remove a
               ! longitudinal restart component as well as round-off drift.
               bfield(:, na) = bfield(:, na) - sum(bfield(:, na)*e_i)*e_i
            end if
         else if (i_cons == 4) then
            ! i_cons = 4 means that we try to use a PID regulator
            !
            !
            !
            if (diagnostics_enabled) write (stdout, '(4x,a)') ' | AMN-PID noncolinear constraints '
            !
            ! Check moment magnitude
            mnorm = sqrt(mom_in(1, na)**2 + mom_in(2, na)**2 + mom_in(3, na)**2)
            c_in = bfield(:, na)
            ma = mnorm
            e_out = 0.0_dblprec
            if (norm2(mom_ref(:, na)) > tiny(1.0_dblprec)) then
               e_i = mom_ref(:, na)/norm2(mom_ref(:, na))
            else
               e_i = [0.0_dblprec, 0.0_dblprec, 1.0_dblprec]
            end if
            if (diagnostics_enabled) write (stdout, '(4x,a,i4)') " | - atom: ", na
            if (mnorm .lt. induced_mom_thresh) then
               if (diagnostics_enabled) write (stdout, '(2x,a,i4,a,f10.4)') ' | Local magnetization for atom ', na, ' is less than threshold', mnorm
               m_delta = 0.0_dblprec
            else
               ! Direction only
               e_out = mom_in(:, na)/mnorm
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
            ! Check to don´t mix first iteration
            dd_delta(:, na) = 0.0_dblprec
            if (norm2(d_delta(:, na)) > 1e-15) dd_delta(:, na) = m_delta - d_delta(:, na)
            !
            if (diagnostics_enabled) then
               write (stdout, '(4x,a,i4,3f15.8)') " | Output moments     for atom ", na, mom_in(1:3, na)
               write (stdout, '(4x,a,i4,3f15.8)') " | Input direction    for atom ", na, mom_ref(:, na)
               write (stdout, '(4x,a,i4,3f15.8)') " | Outut direction    for atom ", na, e_out
               write (stdout, '(4x,a,i4,3f15.8)') " | Input field        for atom ", na, bfield(1:3, na)
            end if
            !
            ! Check to don´t mix first iteration
            if (norm2(d_delta(:, na)) > 1e-15) s_delta(:, na) = s_delta(:, na) + m_delta

            !bfield(:,na)=lambda_t*(1.20_dblprec*m_delta+0.35_dblprec*s_delta(:,na)+0.10_dblprec*dd_delta(:,na))
            if (b_constr_iter == 0) then
               bfield(:, na) = c_in + lambda*(1.30_dblprec*m_delta + 0.35_dblprec*s_delta(:, na) - 0.10_dblprec*dd_delta(:, na))
            else
               bfield(:, na) = lambda*(1.30_dblprec*m_delta + 0.35_dblprec*s_delta(:, na) - 0.10_dblprec*dd_delta(:, na))
            end if

            !bfield_pts(:,ir)=lambda_t*(1.00_dblprec*m_delta+0.12_dblprec*s_delta_pts(:,ir)+0.10_dblprec*dd_delta(:,na))   !ok for grids

            ! Calculate Zeeman-like constraining energy cost
            etcon = etcon + 0.5_dblprec*lambda*sum(m_delta**2)
            d_delta(:, na) = m_delta
            !
            if (diagnostics_enabled) then
               write (stdout, '(4x,a,i4,3f15.8)') " | P  contribution    for atom ", na, d_delta(:, na)
               write (stdout, '(4x,a,i4,3f15.8)') " | I  contribution    for atom ", na, s_delta(:, na)
               write (stdout, '(4x,a,i4,3f15.8)') " | D  contribution    for atom ", na, dd_delta(:, na)
               write (stdout, '(4x,a,i4,4f15.8)') " | Constraining field for atom ", na, bfield(1:3, na)
               write (stdout, '(4x,a,i4,3f15.4)') " | Constraining field for atom (t)", na, cfd_prefac*bfield(1:3, na)
            end if

            ! PID is also a transverse controller: a longitudinal restart
            ! component is not a torque and must not leak into the next
            ! Hamiltonian build.
            bfield(:, na) = bfield(:, na) - sum(bfield(:, na)*e_i)*e_i

         else if (i_cons == 5) then
            ! i_cons = 5 means that we try the Ma-Dudarev approach
            ! which is very analogous to the normal Lagrange approach
            !
            !
            if (diagnostics_enabled) write (stdout, '(2x,a)') ' Ma-Dudarev constraints '
            constrained_mom_err = 0.0_dblprec
            !
            ! Check moment magnitude
            ma = dsqrt(mom_in(1, na)**2 + mom_in(2, na)**2 + mom_in(3, na)**2)
            if (diagnostics_enabled) write (stdout, '(4x,a,i4)') " | - Atom: ", na
            !
            if (ma .lt. induced_mom_thresh) then
               if (diagnostics_enabled) write (stdout, '(2x,a,i4,a,f10.4)') ' | Local magnetization for atom ', na, ' is less than threshold', ma
               bfield(:, na) = 0.0_dblprec
            else
               if (diagnostics_enabled) then
                  write (stdout, '(4x,a,i4,3f15.8)') " | Output moments     for atom ", na, mom_in(1:3, na)
                  write (stdout, '(4x,a,i4,3f15.8)') " | Input direction    for atom ", na, mom_ref(:, na)
                  write (stdout, '(4x,a,i4,3f15.8)') " | Input field        for atom ", na, bfield(1:3, na)
               end if
               !
               e_out = mom_in(:, na)
               if (norm2(mom_ref(:, na)) <= tiny(1.0_dblprec)) then
                  e_i = [0.0_dblprec, 0.0_dblprec, ma]
               else
                  e_i = mom_ref(:, na)/norm2(mom_ref(:, na))*ma
               end if
               !e_dot=e_out(1)*e_i (1)+e_out(2)*e_i (2)+e_out(3)*e_i (3)
               ! new constraining field
               !bfield_new=-10.0_dblprec*lambda_t*(e_out-e_i)
               bfield_new = lambda_t*(e_out - e_i)
               ! gram-schmidt orthogonalization ( b _|_ m)
               !bfield_new=bfield_new-sum(e_out*bfield_new)*e_out/ma**2
               bfield_new = bfield_new - sum(e_i*bfield_new)*e_i/ma**2
               ! mixing field
               bfield(:, na) = (1.0_dblprec - bfield_beta)*bfield(:, na) - bfield_beta*bfield_new
               bfield(:, na) = bfield(:, na) - sum(bfield(:, na)*e_i)*e_i/ma**2
               !
               ! calculate zeeman-like constraining energy cost
               etcon = etcon + lambda_t*(ma - sum(mom_in(1:3, na)*e_i/ma))
               constrained_mom_err = constrained_mom_err + sum(e_out - e_i)**2
               ! write (stdout,´(4x,a,i4,3f15.8)´ ) " | new field          for atom ",na,bfield_new
               if (diagnostics_enabled) then
                  write (stdout, '(4x,a,i4,4f15.8)') " | Output field       for atom ", na, bfield(1:3, na), &
                     sum(bfield(:, na)*e_i)
                  write (stdout, '(4x,a,i4,4f15.8)') " | Output field (t)   for atom ", na, &
                     bfield(:, na)*cfd_prefac, sum(bfield(:, na)*e_i)
                  write (stdout, '(4x,a,i4,3f15.8)') " | Output direction    for atom ", na, e_out/ma
               end if
            end if
         end if

         ! Keep all transverse controller variants free of a longitudinal
         ! restart component, including the low-moment PID branch.
         if ((i_cons == 3 .or. i_cons == 4) .and. norm2(mom_ref(:, na)) > tiny(1.0_dblprec)) then
            e_i = mom_ref(:, na)/norm2(mom_ref(:, na))
            bfield(:, na) = bfield(:, na) - sum(bfield(:, na)*e_i)*e_i
         end if

         ! Diagnostics are derived from the measured moment and the field that
         ! will be consumed by the next Hamiltonian build.  Keep this outside
         ! the controller branches so every supported algorithm reports the
         ! same quantities.
         ma = norm2(mom_in(:, na))
         if (norm2(mom_ref(:, na)) > tiny(1.0_dblprec)) then
            e_i = mom_ref(:, na)/norm2(mom_ref(:, na))
         else
            e_i = [0.0_dblprec, 0.0_dblprec, 1.0_dblprec]
         end if
         if (ma > tiny(1.0_dblprec)) then
            e_out = mom_in(:, na)/ma
            cosine = max(-1.0_dblprec, min(1.0_dblprec, sum(e_out*e_i)))
            transverse = e_out - sum(e_out*e_i)*e_i
            ! atan2 retains angular resolution when the residual is below
            ! 1e-6 rad; acos(dot) loses those digits through cancellation.
            angle = atan2(norm2(transverse), cosine)
            torque = cross3(e_i, e_out)
         else
            e_out = 0.0_dblprec
            angle = 0.5_dblprec*acos(-1.0_dblprec)
            transverse = 0.0_dblprec
            torque = 0.0_dblprec
         end if
         metric_sum = metric_sum + angle*angle
         if (present(diagnostics_out)) then
            diagnostics_out%target(:, na) = e_i
            diagnostics_out%actual(:, na) = e_out
            diagnostics_out%moment_magnitude(na) = ma
            diagnostics_out%angular_error(na) = angle
            diagnostics_out%transverse_residual(:, na) = transverse
            diagnostics_out%bfield(:, na) = bfield(:, na)
            diagnostics_out%bfield_magnitude(na) = norm2(bfield(:, na))
            diagnostics_out%bfield_longitudinal(na) = sum(bfield(:, na)*e_i)
            diagnostics_out%torque(:, na) = torque
         end if
      end do ! na
      !
      constrained_mom_err = sqrt(constrained_mom_err/max(1, ndim))
      ! `bfield` is the field that will be consumed by the next Hamiltonian
      ! build.  Since mom_in is the measured (un-normalized) state supplied by
      ! the SCF caller, this is the expectation value of the added Pauli term.
      ! It is deliberately kept separate from etcon: the latter is a
      ! controller residual and neither quantity is a physical DFT energy.
      field_coupling = sum(bfield*mom_in)
      if (present(diagnostics_out)) then
         diagnostics_out%convergence_metric = sqrt(metric_sum/real(max(1, ndim), dblprec))
         diagnostics_out%max_angular_error = maxval(diagnostics_out%angular_error)
         diagnostics_out%controller_penalty_energy = etcon
         diagnostics_out%field_coupling_energy = field_coupling
         diagnostics_out%lagrange_functional_represented = .false.
      end if

      if (.not. update_lambda_) then
         if (present(etcon_out)) etcon_out = etcon
         deallocate (mom_tmp)
         return
      end if

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
      ! Scale up the controller strength when the merit value is small.
      if (etcon < 1.0d-2) lambda_t = min(1.2_dblprec*lambda_t, 1.0e4_dblprec)
      if (i_cons == 3) lambda_t = min(lambda_t*(2.0_dblprec - lambda_t/lambda), lambda)
      !if(i_cons==5) lambda_t=min(lambda_t*(2.0_dblprec-lambda_t/lambda),lambda)
      if (i_cons == 5) lambda_t = min(lambda_t + 2.0_dblprec, lambda)
      ! if(i_cons==4) lambda_t=min(lambda_t+1.0_dblprec,25.0_dblprec)
      !if(i_cons==4) lambda_t=lambda_t+1.0_dblprec
      !if(i_cons==5) lambda_t=lambda_t+lambda_t*min(0.1_dblprec,etcon**2)
      !if(i_cons==5) lambda_t=lambda_t*(1.0_dblprec+0.5_dblprec*(etcon))
      !if(i_cons==5) lambda_t=lambda_t*(1.0_dblprec+2.0_dblprec*min(constrained_mom_err,0.1_dblprec))
      if (diagnostics_enabled .and. i_cons == 5) write (stdout, '(4x,a,f12.4a,g10.2)') " | New lambda_t: ", lambda_t, "     error: ", constrained_mom_err
      if (diagnostics_enabled) write (stdout, '(4x,a)') " | -  "
      if (present(etcon_out)) etcon_out = etcon
      deallocate (mom_tmp)
      return
   end subroutine constrain

   subroutine clear_constraint_diagnostics(this)
      class(constraint_diagnostics), intent(inout) :: this

      if (allocated(this%target)) deallocate(this%target)
      if (allocated(this%actual)) deallocate(this%actual)
      if (allocated(this%moment_magnitude)) deallocate(this%moment_magnitude)
      if (allocated(this%angular_error)) deallocate(this%angular_error)
      if (allocated(this%transverse_residual)) deallocate(this%transverse_residual)
      if (allocated(this%bfield)) deallocate(this%bfield)
      if (allocated(this%bfield_magnitude)) deallocate(this%bfield_magnitude)
      if (allocated(this%bfield_longitudinal)) deallocate(this%bfield_longitudinal)
      if (allocated(this%torque)) deallocate(this%torque)
      this%iteration = 0
      this%nsite = 0
      this%convergence_metric = 0.0_dblprec
      this%max_angular_error = 0.0_dblprec
      this%controller_penalty_energy = 0.0_dblprec
      this%field_coupling_energy = 0.0_dblprec
      this%lagrange_functional_represented = .false.
   end subroutine clear_constraint_diagnostics

   pure function cross3(a, b) result(c)
      real(dblprec), intent(in) :: a(3), b(3)
      real(dblprec) :: c(3)

      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)
   end function cross3

end module cfd
