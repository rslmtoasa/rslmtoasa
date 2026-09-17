!------------------------------------------------------------------------------
! DRESP-03G -- direct finite-H complex-energy resolvent bridge.
!
! This module deliberately knows only about ordinary finite matrices.  It does
! not import reciprocal eigenpairs, Lehmann kernels, LMTO P(E), or the native
! exchange implementation.  The production caller supplies the live H(k),
! H(k+q), T(q), T(-q), and C matrices.
!------------------------------------------------------------------------------
module lr_kl_contour_mod
   use precision_mod, only: rp
   use math_mod, only: i_unit, pi
   implicit none
   private

   type, public :: finite_h_contour_options
      integer :: contour_points = 32
      character(len=16) :: contour_shape = 'ellipse'
      real(rp) :: contour_margin = 0.25_rp
      real(rp) :: contour_height_fraction = 0.35_rp
      logical :: account_fermi_poles = .true.
   end type finite_h_contour_options

   type, public :: finite_h_contour_report
      integer :: contour_points = 0
      integer :: fermi_poles = 0
      real(rp) :: solve_seconds = 0.0_rp
      real(rp) :: contour_seconds = 0.0_rp
      real(rp) :: pole_seconds = 0.0_rp
   end type finite_h_contour_report

   public :: build_finite_temperature_contour
   public :: build_zero_temperature_occupied_contour
   public :: finite_temperature_complex_fermi
   public :: finite_temperature_regularized_fermi
   public :: force_theorem_finite_q_hessian_from_resolvent
   public :: force_theorem_finite_q_hessian_from_resolvent_batch
   public :: force_theorem_finite_q_hessian_from_zero_temperature_contour

   interface
      subroutine zgetrf(n, m, a, lda, ipiv, info)
         import rp
         integer, intent(in) :: n, m, lda
         complex(rp), intent(inout) :: a(lda, *)
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
      end subroutine zgetrf
      subroutine zgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
         import rp
         character(len=1), intent(in) :: trans
         integer, intent(in) :: n, nrhs, lda, ldb
         complex(rp), intent(in) :: a(lda, *)
         integer, intent(in) :: ipiv(*)
         complex(rp), intent(inout) :: b(ldb, *)
         integer, intent(out) :: info
      end subroutine zgetrs
   end interface

contains

   !> Construct the counter-clockwise ellipse used by the finite-T bridge.
   !>
   !> The ellipse encloses the real-axis spectra of both supplied matrices.
   !> Its semiminor axis is intentionally allowed to cross Matsubara poles.
   !> Those poles are returned explicitly so that their residues can be added
   !> back.  This permits a well-conditioned, broad contour at low temperature
   !> without silently changing the Fermi operator.
   subroutine build_finite_temperature_contour(h_source, h_endpoint, fermi, kT, options, nodes, weights, fermi_poles)
      complex(rp), intent(in) :: h_source(:, :), h_endpoint(:, :)
      real(rp), intent(in) :: fermi, kT
      type(finite_h_contour_options), intent(in) :: options
      complex(rp), allocatable, intent(out) :: nodes(:), weights(:), fermi_poles(:)
      real(rp) :: lower, upper, center, semimajor, semiminor, theta, dtheta
      real(rp) :: desired_height, pole_y, ellipse_value
      integer :: i, l, max_l, npoint, npole

      if (kT <= 0.0_rp) error stop 'build_finite_temperature_contour: kT must be positive'
      call validate_options(options)
      call hermitian_bounds(h_source, h_endpoint, lower, upper)
      center = 0.5_rp*(lower + upper)
      semimajor = max(0.5_rp*(upper-lower) + options%contour_margin, options%contour_margin)
      desired_height = max(options%contour_height_fraction*semimajor, 1.0e-8_rp)
      if (.not. options%account_fermi_poles) desired_height = min(desired_height, 0.45_rp*pi*kT)
      semiminor = desired_height

      npoint = options%contour_points
      allocate(nodes(npoint), weights(npoint))
      dtheta = 2.0_rp*pi/real(npoint,rp)
      do i = 1, npoint
         theta = dtheta*real(i-1,rp)
         nodes(i) = cmplx(center + semimajor*cos(theta), semiminor*sin(theta), rp)
         ! z(theta) is counter-clockwise for theta increasing from zero.
         weights(i) = ((-semimajor*sin(theta) + i_unit*semiminor*cos(theta))*dtheta)/(2.0_rp*pi*i_unit)
      end do

      if (options%account_fermi_poles) then
         max_l = max(2, int(ceiling(semiminor/(pi*kT))) + 2)
         npole = 0
         do l = -max_l, max_l
            pole_y = pi*kT*real(2*l+1,rp)
            ellipse_value = ((fermi-center)/semimajor)**2 + (pole_y/semiminor)**2
            if (ellipse_value < 1.0_rp-1.0e-12_rp) npole = npole + 1
         end do
         allocate(fermi_poles(npole))
         npole = 0
         do l = -max_l, max_l
            pole_y = pi*kT*real(2*l+1,rp)
            ellipse_value = ((fermi-center)/semimajor)**2 + (pole_y/semiminor)**2
            if (ellipse_value < 1.0_rp-1.0e-12_rp) then
               npole = npole + 1
               fermi_poles(npole) = cmplx(fermi,pole_y,rp)
            end if
         end do
      else
         allocate(fermi_poles(0))
      end if
   end subroutine build_finite_temperature_contour

   !> Construct an occupied-state contour for the T=0 limit.
   !>
   !> `occupied_bounds=[emin,emax]` must lie strictly inside the insulating
   !> gap, with emax below the first unoccupied eigenvalue.  The weight is
   !> f(z)=1, and no finite-T Fermi poles are present.  The routine does not
   !> diagonalize either matrix; the bounds are an explicit fixture/driver
   !> contract for this classical limit.
   subroutine build_zero_temperature_occupied_contour(occupied_bounds, options, nodes, weights)
      real(rp), intent(in) :: occupied_bounds(2)
      type(finite_h_contour_options), intent(in) :: options
      complex(rp), allocatable, intent(out) :: nodes(:), weights(:)
      real(rp) :: center, semimajor, semiminor, theta, dtheta
      integer :: i, npoint

      call validate_options(options)
      if (occupied_bounds(2) <= occupied_bounds(1)) then
         error stop 'build_zero_temperature_occupied_contour: invalid occupied bounds'
      end if
      center = 0.5_rp*sum(occupied_bounds)
      semimajor = 0.5_rp*(occupied_bounds(2)-occupied_bounds(1)) + options%contour_margin
      semiminor = max(options%contour_height_fraction*semimajor, 1.0e-8_rp)
      npoint = options%contour_points
      allocate(nodes(npoint), weights(npoint))
      dtheta = 2.0_rp*pi/real(npoint,rp)
      do i = 1, npoint
         theta = dtheta*real(i-1,rp)
         nodes(i) = cmplx(center + semimajor*cos(theta), semiminor*sin(theta), rp)
         weights(i) = ((-semimajor*sin(theta) + i_unit*semiminor*cos(theta))*dtheta)/(2.0_rp*pi*i_unit)
      end do
   end subroutine build_zero_temperature_occupied_contour

   !> Stable complex Fermi function on a contour node.
   pure complex(rp) function finite_temperature_complex_fermi(z, fermi, kT) result(value)
      complex(rp), intent(in) :: z
      real(rp), intent(in) :: fermi, kT
      complex(rp) :: argument, reduced

      if (kT <= 0.0_rp) error stop 'finite_temperature_complex_fermi: kT must be positive'
      argument = (z-fermi)/kT
      if (real(argument,rp) > 40.0_rp) then
         reduced = exp(-argument)
         value = reduced/(1.0_rp+reduced)
      else if (real(argument,rp) < -40.0_rp) then
         reduced = exp(argument)
         value = 1.0_rp/(1.0_rp+reduced)
      else
         value = 1.0_rp/(exp(argument)+1.0_rp)
      end if
   end function finite_temperature_complex_fermi

   !> Fermi function with all poles enclosed by the selected contour removed.
   !>
   !> The residue of f at z_l=mu+i*pi*(2l+1)kT is -kT.  Adding kT/(z-z_l)
   !> cancels that pole.  Cauchy theorem then gives
   !>
   !>   integral_C f_reg(z) R(z) dz/(2*pi*i)
   !>     = integral_C f(z) R(z) dz/(2*pi*i) + kT sum_l R(z_l),
   !>
   !> which is precisely the physical-spectrum contour after the Fermi-pole
   !> residues have been removed.  This regularization is numerically much
   !> better than integrating a meromorphic weight close to its last enclosed
   !> pole, while remaining an exact finite-dimensional identity.
   pure complex(rp) function finite_temperature_regularized_fermi(z, fermi, kT, fermi_poles) result(value)
      complex(rp), intent(in) :: z, fermi_poles(:)
      real(rp), intent(in) :: fermi, kT
      integer :: i
      value = finite_temperature_complex_fermi(z,fermi,kT)
      do i = 1, size(fermi_poles)
         value = value + kT/(z-fermi_poles(i))
      end do
   end function finite_temperature_regularized_fermi

   !> Finite-T direct-resolvent Hessian for one k -> k+q pair.
   !>
   !> H_source is H(k), H_endpoint is H(k+q).  T_q has rows at k+q and
   !> columns at k; T_minus_q has rows at k and columns at k+q.  The returned
   !> TT term is the explicit ordered-pair average, while contact is
   !> Tr[f(H_k) C(q,-q;k)].
   subroutine force_theorem_finite_q_hessian_from_resolvent(h_source, h_endpoint, fermi, kT, torques_q, torques_minus_q, mixed, &
      options, hessian, torque_torque, mixed_contact, complete, report)
      complex(rp), intent(in) :: h_source(:, :), h_endpoint(:, :)
      real(rp), intent(in) :: fermi, kT
      complex(rp), intent(in) :: torques_q(:, :, :), torques_minus_q(:, :, :), mixed(:, :, :, :)
      type(finite_h_contour_options), intent(in) :: options
      complex(rp), intent(out) :: hessian(:, :), torque_torque(:, :), mixed_contact(:, :), complete(:, :)
      type(finite_h_contour_report), intent(out), optional :: report
      complex(rp), allocatable :: nodes(:), weights(:), poles(:)

      call build_finite_temperature_contour(h_source, h_endpoint, fermi, kT, options, nodes, weights, poles)
      call resolvent_hessian_core(h_source, h_endpoint, fermi, kT, .true., nodes, weights, poles, torques_q, torques_minus_q, mixed, &
         hessian, torque_torque, mixed_contact, complete, report)
      deallocate(nodes, weights, poles)
   end subroutine force_theorem_finite_q_hessian_from_resolvent

   !> Brillouin-zone average of the direct-resolvent Hessian.
   subroutine force_theorem_finite_q_hessian_from_resolvent_batch(h_source, h_endpoint, fermi, kT, k_weights, torques_q, &
      torques_minus_q, mixed, options, hessian, torque_torque, mixed_contact, complete, report)
      complex(rp), intent(in) :: h_source(:, :, :), h_endpoint(:, :, :)
      real(rp), intent(in) :: fermi, kT, k_weights(:)
      complex(rp), intent(in) :: torques_q(:, :, :, :), torques_minus_q(:, :, :, :), mixed(:, :, :, :, :)
      type(finite_h_contour_options), intent(in) :: options
      complex(rp), intent(out) :: hessian(:, :), torque_torque(:, :), mixed_contact(:, :), complete(:, :)
      type(finite_h_contour_report), intent(out), optional :: report
      complex(rp), allocatable :: h(:, :), tt(:, :), cc(:, :), allh(:, :)
      type(finite_h_contour_report) :: one_report
      real(rp) :: weight_sum
      integer :: ik, nk, nsite

      nk = size(h_source,3); nsite = size(torques_q,3)
      if (size(h_endpoint,3) /= nk .or. size(k_weights) /= nk .or. size(torques_q,4) /= nk .or. &
          size(torques_minus_q,4) /= nk .or. size(mixed,5) /= nk) error stop 'resolvent batch: k dimension mismatch'
      weight_sum = sum(k_weights)
      if (weight_sum <= tiny(1.0_rp)) error stop 'resolvent batch: zero k-weight sum'
      allocate(h(nsite,nsite), tt(nsite,nsite), cc(nsite,nsite), allh(nsite,nsite))
      hessian = 0.0_rp; torque_torque = 0.0_rp; mixed_contact = 0.0_rp; complete = 0.0_rp
      if (present(report)) report = finite_h_contour_report()
      do ik = 1, nk
         call force_theorem_finite_q_hessian_from_resolvent(h_source(:,:,ik), h_endpoint(:,:,ik), fermi, kT, torques_q(:,:,:,ik), &
            torques_minus_q(:,:,:,ik), mixed(:,:,:,:,ik), options, h, tt, cc, allh, one_report)
         hessian = hessian + k_weights(ik)*h/weight_sum
         torque_torque = torque_torque + k_weights(ik)*tt/weight_sum
         mixed_contact = mixed_contact + k_weights(ik)*cc/weight_sum
         complete = complete + k_weights(ik)*allh/weight_sum
         if (present(report)) then
            report%contour_points = report%contour_points + one_report%contour_points
            report%fermi_poles = report%fermi_poles + one_report%fermi_poles
            report%solve_seconds = report%solve_seconds + one_report%solve_seconds
            report%contour_seconds = report%contour_seconds + one_report%contour_seconds
            report%pole_seconds = report%pole_seconds + one_report%pole_seconds
         end if
      end do
      deallocate(h, tt, cc, allh)
   end subroutine force_theorem_finite_q_hessian_from_resolvent_batch

   !> Zero-temperature occupied-contour version for a gapped fixture.
   subroutine force_theorem_finite_q_hessian_from_zero_temperature_contour(h_source, h_endpoint, occupied_bounds, torques_q, &
      torques_minus_q, mixed, options, hessian, torque_torque, mixed_contact, complete, report)
      complex(rp), intent(in) :: h_source(:, :), h_endpoint(:, :)
      real(rp), intent(in) :: occupied_bounds(2)
      complex(rp), intent(in) :: torques_q(:, :, :), torques_minus_q(:, :, :), mixed(:, :, :, :)
      type(finite_h_contour_options), intent(in) :: options
      complex(rp), intent(out) :: hessian(:, :), torque_torque(:, :), mixed_contact(:, :), complete(:, :)
      type(finite_h_contour_report), intent(out), optional :: report
      complex(rp), allocatable :: nodes(:), weights(:), poles(:)

      call build_zero_temperature_occupied_contour(occupied_bounds, options, nodes, weights)
      allocate(poles(0))
      call resolvent_hessian_core(h_source, h_endpoint, 0.0_rp, 0.0_rp, .false., nodes, weights, poles, torques_q, torques_minus_q, mixed, &
         hessian, torque_torque, mixed_contact, complete, report)
      deallocate(nodes, weights, poles)
   end subroutine force_theorem_finite_q_hessian_from_zero_temperature_contour

   ! The common contraction is intentionally expressed in orbital matrix space
   ! rather than in an eigenbasis.  Each energy node performs direct LU
   ! factorization of z I-H and solves for the vertex/contact right-hand sides.
   subroutine resolvent_hessian_core(h_source, h_endpoint, fermi, kT, finite_temperature, nodes, weights, poles, torques_q, &
      torques_minus_q, mixed, hessian, torque_torque, mixed_contact, complete, report)
      complex(rp), intent(in) :: h_source(:, :), h_endpoint(:, :), nodes(:), weights(:), poles(:)
      real(rp), intent(in) :: fermi, kT
      logical, intent(in) :: finite_temperature
      complex(rp), intent(in) :: torques_q(:, :, :), torques_minus_q(:, :, :), mixed(:, :, :, :)
      complex(rp), intent(out) :: hessian(:, :), torque_torque(:, :), mixed_contact(:, :), complete(:, :)
      type(finite_h_contour_report), intent(out), optional :: report
      complex(rp), allocatable :: a0(:, :), a1(:, :), rhs0(:, :), rhs1(:, :)
      complex(rp) :: coefficient, trace_first, trace_second, trace_contact
      complex(rp), allocatable :: contour_tt(:,:), contour_contact(:,:), pole_tt(:,:), pole_contact(:,:)
      integer, allocatable :: piv0(:), piv1(:)
      integer :: nmat, nsite, nnode, npole, a, b, info, clock_start, clock_stop, clock_rate
      integer :: contour_start, contour_stop, pole_start, pole_stop
      real(rp) :: solve_seconds, contour_seconds, pole_seconds

      nmat = size(h_source,1); nsite = size(torques_q,3); nnode = size(nodes); npole = size(poles)
      if (size(h_source,2) /= nmat .or. any(shape(h_endpoint) /= [nmat,nmat]) .or. size(weights) /= nnode .or. &
          size(torques_q,1) /= nmat .or. size(torques_q,2) /= nmat .or. size(torques_minus_q,1) /= nmat .or. &
          size(torques_minus_q,2) /= nmat .or. size(torques_q,3) /= nsite .or. size(torques_minus_q,3) /= nsite .or. &
          size(mixed,1) /= nmat .or. size(mixed,2) /= nmat .or. size(mixed,3) /= nsite .or. size(mixed,4) /= nsite .or. &
          size(hessian,1) /= nsite .or. size(hessian,2) /= nsite .or. size(torque_torque,1) /= nsite .or. &
          size(mixed_contact,1) /= nsite .or. size(mixed_contact,2) /= nsite .or. size(complete,1) /= nsite .or. &
          size(complete,2) /= nsite) error stop 'resolvent Hessian: shape mismatch'

      allocate(a0(nmat,nmat), a1(nmat,nmat), rhs0(nmat,nmat), rhs1(nmat,nmat), piv0(nmat), piv1(nmat), &
         contour_tt(nsite,nsite), contour_contact(nsite,nsite), pole_tt(nsite,nsite), pole_contact(nsite,nsite))
      contour_tt = 0.0_rp; contour_contact = 0.0_rp; pole_tt = 0.0_rp; pole_contact = 0.0_rp
      solve_seconds = 0.0_rp; contour_seconds = 0.0_rp; pole_seconds = 0.0_rp
      call system_clock(count_rate=clock_rate)

      call system_clock(contour_start)
      do nnode = 1, size(nodes)
         a0 = -h_source; a1 = -h_endpoint
         do a = 1, nmat
            a0(a,a) = a0(a,a) + nodes(nnode)
            a1(a,a) = a1(a,a) + nodes(nnode)
         end do
         call system_clock(clock_start)
         call zgetrf(nmat, nmat, a0, nmat, piv0, info)
         if (info /= 0) error stop 'resolvent Hessian: source LU factorization failed'
         call zgetrf(nmat, nmat, a1, nmat, piv1, info)
         if (info /= 0) error stop 'resolvent Hessian: endpoint LU factorization failed'
         call system_clock(clock_stop)
         solve_seconds = solve_seconds + elapsed_seconds(clock_start,clock_stop,clock_rate)
         coefficient = weights(nnode)
         if (finite_temperature) coefficient = coefficient*finite_temperature_regularized_fermi(nodes(nnode),fermi,kT,poles)
         call accumulate_at_factorized_node(a0,a1,piv0,piv1,coefficient,contour_tt,contour_contact)
      end do
      call system_clock(contour_stop)
      contour_seconds = elapsed_seconds(contour_start,contour_stop,clock_rate)

      ! The regularized contour integrand has no enclosed Fermi poles, but the
      ! physical-spectrum contour is recovered by adding their residues.  The
      ! explicit solve here is part of the exact finite-T identity, not a
      ! quadrature correction or a spectral reconstruction.
      call system_clock(pole_start)
      do nnode = 1, size(poles)
         a0 = -h_source; a1 = -h_endpoint
         do a = 1, nmat
            a0(a,a) = a0(a,a) + poles(nnode)
            a1(a,a) = a1(a,a) + poles(nnode)
         end do
         call system_clock(clock_start)
         call zgetrf(nmat, nmat, a0, nmat, piv0, info)
         if (info /= 0) error stop 'resolvent Hessian: source pole LU factorization failed'
         call zgetrf(nmat, nmat, a1, nmat, piv1, info)
         if (info /= 0) error stop 'resolvent Hessian: endpoint pole LU factorization failed'
         call system_clock(clock_stop)
         solve_seconds = solve_seconds + elapsed_seconds(clock_start,clock_stop,clock_rate)
         call accumulate_at_factorized_node(a0,a1,piv0,piv1,cmplx(kT,0.0_rp,rp),pole_tt,pole_contact)
      end do
      call system_clock(pole_stop)
      pole_seconds = elapsed_seconds(pole_start,pole_stop,clock_rate)

      torque_torque = contour_tt + pole_tt
      mixed_contact = contour_contact + pole_contact
      hessian = torque_torque + mixed_contact
      complete = hessian
      if (present(report)) then
         report%contour_points = size(nodes)
         report%fermi_poles = size(poles)
         report%solve_seconds = solve_seconds
         report%contour_seconds = contour_seconds
         report%pole_seconds = pole_seconds
      end if
      deallocate(a0,a1,rhs0,rhs1,piv0,piv1,contour_tt,contour_contact,pole_tt,pole_contact)

   contains

      subroutine accumulate_at_factorized_node(lu0,lu1,ip0,ip1,weight,tt_sum,contact_sum)
         complex(rp), intent(in) :: lu0(:, :), lu1(:, :), weight
         integer, intent(in) :: ip0(:), ip1(:)
         complex(rp), intent(inout) :: tt_sum(:,:), contact_sum(:,:)
         integer :: ia, ib, solve_start, solve_stop, local_info

         do ia = 1, nsite
            do ib = 1, nsite
               rhs0 = torques_minus_q(:,:,ib)
               call system_clock(solve_start)
               call zgetrs('N',nmat,nmat,lu0,nmat,ip0,rhs0,nmat,local_info)
               if (local_info /= 0) error stop 'resolvent Hessian: source vertex solve failed'
               rhs1 = matmul(torques_q(:,:,ia),rhs0)
               call zgetrs('N',nmat,nmat,lu1,nmat,ip1,rhs1,nmat,local_info)
               if (local_info /= 0) error stop 'resolvent Hessian: endpoint vertex solve failed'
               trace_first = trace_matrix(rhs1)

               rhs0 = torques_minus_q(:,:,ia)
               call zgetrs('N',nmat,nmat,lu0,nmat,ip0,rhs0,nmat,local_info)
               if (local_info /= 0) error stop 'resolvent Hessian: reverse source vertex solve failed'
               rhs1 = matmul(torques_q(:,:,ib),rhs0)
               call zgetrs('N',nmat,nmat,lu1,nmat,ip1,rhs1,nmat,local_info)
               if (local_info /= 0) error stop 'resolvent Hessian: reverse endpoint vertex solve failed'
               trace_second = trace_matrix(rhs1)

               rhs0 = mixed(:,:,ia,ib)
               call zgetrs('N',nmat,nmat,lu0,nmat,ip0,rhs0,nmat,local_info)
               if (local_info /= 0) error stop 'resolvent Hessian: contact solve failed'
               trace_contact = trace_matrix(rhs0)
               tt_sum(ia,ib) = tt_sum(ia,ib) + weight*0.5_rp*(trace_first+trace_second)
               contact_sum(ia,ib) = contact_sum(ia,ib) + weight*trace_contact
               call system_clock(solve_stop)
               solve_seconds = solve_seconds + elapsed_seconds(solve_start,solve_stop,clock_rate)
            end do
         end do
      end subroutine accumulate_at_factorized_node
   end subroutine resolvent_hessian_core

   subroutine validate_options(options)
      type(finite_h_contour_options), intent(in) :: options
      if (options%contour_points < 8) error stop 'finite-H contour: contour_points must be at least 8'
      if (trim(options%contour_shape) /= 'ellipse') error stop 'finite-H contour: only contour_shape=ellipse is implemented'
      if (options%contour_margin <= 0.0_rp) error stop 'finite-H contour: contour_margin must be positive'
      if (options%contour_height_fraction <= 0.0_rp) error stop 'finite-H contour: height fraction must be positive'
   end subroutine validate_options

   subroutine hermitian_bounds(h_source, h_endpoint, lower, upper)
      complex(rp), intent(in) :: h_source(:, :), h_endpoint(:, :)
      real(rp), intent(out) :: lower, upper
      real(rp) :: lo0, hi0, lo1, hi1
      call matrix_gershgorin_bounds(h_source,lo0,hi0)
      call matrix_gershgorin_bounds(h_endpoint,lo1,hi1)
      lower = min(lo0,lo1); upper = max(hi0,hi1)
   end subroutine hermitian_bounds

   subroutine matrix_gershgorin_bounds(matrix, lower, upper)
      complex(rp), intent(in) :: matrix(:, :)
      real(rp), intent(out) :: lower, upper
      real(rp) :: radius
      integer :: i
      lower = huge(1.0_rp); upper = -huge(1.0_rp)
      do i = 1, size(matrix,1)
         radius = sum(abs(matrix(i,:))) - abs(matrix(i,i))
         lower = min(lower,real(matrix(i,i),rp)-radius)
         upper = max(upper,real(matrix(i,i),rp)+radius)
      end do
   end subroutine matrix_gershgorin_bounds

   pure function trace_matrix(matrix) result(value)
      complex(rp), intent(in) :: matrix(:, :)
      complex(rp) :: value
      integer :: i
      value = 0.0_rp
      do i = 1, min(size(matrix,1),size(matrix,2))
         value = value + matrix(i,i)
      end do
   end function trace_matrix

   pure function elapsed_seconds(start_count, end_count, rate) result(seconds)
      integer, intent(in) :: start_count, end_count, rate
      real(rp) :: seconds
      if (rate > 0) then
         seconds = real(end_count-start_count,rp)/real(rate,rp)
      else
         seconds = 0.0_rp
      end if
   end function elapsed_seconds

end module lr_kl_contour_mod
