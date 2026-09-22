! DRESP-09U independent production representation-map tangent fixture.
!
! It exercises the scalar predls derivative, the post-predls spin lift, the
! live bond algebra, a commuting structure-constant case, and the
! nonmagnetic zero-response control.  Central differences are used only as
! independent oracles; no finite difference is used by the service itself.
program test_lr_lmto_representation_tangent

   use precision_mod, only: rp
   use basis_mod, only: basis_init, nb, norb
   use math_mod, only: init_math_operators, hcpx
   use logger_mod, only: g_logger
   use lmto_magnetic_tangent_mod, only: lmto_bond_value, lmto_hhmag_to_spinor
   use lr_dresp08_native_mapping_mod, only: dresp08_rotate_collinear_operator
   use lr_lmto_representation_tangent_mod, only: lmto_predls_transform, lmto_predls_transform_tangent, &
      lmto_spin_parameter_tangent, lmto_representation_bond_tangent, lmto_representation_relative_residual
   implicit none

   integer, parameter :: nchannel = 3
   real(rp), parameter :: tol = 5.0e-11_rp
   real(rp) :: c(nchannel,2), enu(nchannel,2), srdel(nchannel,2), qpar(nchannel,2)
   real(rp) :: dc(nchannel,2), denu(nchannel,2), dsrdel(nchannel,2), dqpar(nchannel,2), qm(nchannel)
   real(rp) :: center(nchannel,2), shifted(nchannel,2), width(nchannel,2), obar(nchannel,2)
   real(rp) :: dcenter(nchannel,2), dshifted(nchannel,2), dwidth(nchannel,2), dobar(nchannel,2)
   real(rp) :: cp(nchannel,2), cm(nchannel,2), sp(nchannel,2), sm(nchannel,2), wp(nchannel,2), wm(nchannel,2)
   real(rp) :: op(nchannel,2), om(nchannel,2), h, current_error, error_fd(4), ratio(4)
   complex(rp), allocatable :: up(:), down(:), dup(:), ddown(:), matrix(:,:), delta_matrix(:,:)
   complex(rp), allocatable :: up_matrix(:,:), down_matrix(:,:), plus(:,:), minus(:,:)
   complex(rp), allocatable :: hhh(:,:), hcomm(:,:), c0(:), c1(:), wx0i(:), wx1i(:), wx0j(:), wx1j(:)
   complex(rp), allocatable :: zero(:), dbond(:,:), bond_plus(:,:), bond_minus(:,:), hhmag(:,:,:)
   real(rp) :: moment(3), dmoment(3), zero_moment(3), error_spin, error_bond, error_nonmag, error_commuting, sitewise_error
   integer :: i, j, branch
   logical :: failed

   call basis_init(2)
   call init_math_operators()
   call g_logger%init()
   allocate(up(norb),down(norb),dup(norb),ddown(norb),matrix(nb,nb),delta_matrix(nb,nb),up_matrix(norb,norb), &
      down_matrix(norb,norb),plus(nb,nb),minus(nb,nb),hhh(norb,norb),hcomm(norb,norb),c0(norb),c1(norb), &
      wx0i(norb),wx1i(norb),wx0j(norb),wx1j(norb),zero(norb),dbond(nb,nb),bond_plus(nb,nb),bond_minus(nb,nb), &
      hhmag(norb,norb,4))

   do i = 1, nchannel
      do j = 1, 2
         c(i,j) = 0.31_rp + 0.07_rp*real(i,rp) - 0.023_rp*real(j,rp)
         enu(i,j) = -0.19_rp + 0.041_rp*real(i,rp) + 0.017_rp*real(j,rp)
         srdel(i,j) = 0.77_rp + 0.05_rp*real(i,rp) + 0.021_rp*real(j,rp)
         qpar(i,j) = 0.13_rp + 0.031_rp*real(i,rp) - 0.011_rp*real(j,rp)
         dc(i,j) = 0.009_rp*real(i+2*j,rp)
         denu(i,j) = -0.006_rp*real(2*i+j,rp)
         dsrdel(i,j) = 0.004_rp*real(i+j,rp)
         dqpar(i,j) = -0.003_rp*real(2*i+j,rp)
      end do
      qm(i) = 0.04_rp + 0.013_rp*real(i,rp)
   end do

   call lmto_predls_transform_tangent(c, enu, srdel, qpar, dc, denu, dsrdel, dqpar, 0.17_rp, 1.23_rp, 1.41_rp, qm, &
      center, shifted, width, obar, dcenter, dshifted, dwidth, dobar)
   ratio = 0.0_rp
   do branch = 1, 4
      h = 2.0e-3_rp/2.0_rp**real(branch-1,rp)
      call lmto_predls_transform(c+h*dc, enu+h*denu, srdel+h*dsrdel, qpar+h*dqpar, 0.17_rp, 1.23_rp, 1.41_rp, qm, &
         cp, sp, wp, op)
      call lmto_predls_transform(c-h*dc, enu-h*denu, srdel-h*dsrdel, qpar-h*dqpar, 0.17_rp, 1.23_rp, 1.41_rp, qm, &
         cm, sm, wm, om)
      select case (branch)
      case (1)
         current_error = maxval(abs((cp-cm)/(2.0_rp*h)-dcenter))
      case (2)
         current_error = maxval(abs((sp-sm)/(2.0_rp*h)-dshifted))
      case (3)
         current_error = maxval(abs((wp-wm)/(2.0_rp*h)-dwidth))
      case (4)
         current_error = maxval(abs((op-om)/(2.0_rp*h)-dobar))
      end select
      h = h/2.0_rp
      call lmto_predls_transform(c+h*dc, enu+h*denu, srdel+h*dsrdel, qpar+h*dqpar, 0.17_rp, 1.23_rp, 1.41_rp, qm, &
         cp, sp, wp, op)
      call lmto_predls_transform(c-h*dc, enu-h*denu, srdel-h*dsrdel, qpar-h*dqpar, 0.17_rp, 1.23_rp, 1.41_rp, qm, &
         cm, sm, wm, om)
      select case (branch)
      case (1)
         error_fd(branch) = maxval(abs((cp-cm)/(2.0_rp*h)-dcenter))
      case (2)
         error_fd(branch) = maxval(abs((sp-sm)/(2.0_rp*h)-dshifted))
      case (3)
         error_fd(branch) = maxval(abs((wp-wm)/(2.0_rp*h)-dwidth))
      case (4)
         error_fd(branch) = maxval(abs((op-om)/(2.0_rp*h)-dobar))
      end select
      ratio(branch) = current_error/max(error_fd(branch),tiny(1.0_rp))
   end do

   do i = 1, norb
      up(i) = cmplx(0.11_rp+0.007_rp*real(i,rp),0.0_rp,rp)
      down(i) = cmplx(-0.08_rp+0.005_rp*real(i,rp),0.0_rp,rp)
      dup(i) = cmplx(0.0_rp,0.0_rp,rp)
      ddown(i) = cmplx(0.0_rp,0.0_rp,rp)
      up_matrix(i,:) = cmplx(0.0_rp,0.0_rp,rp)
      down_matrix(i,:) = cmplx(0.0_rp,0.0_rp,rp)
      up_matrix(i,i) = up(i)
      down_matrix(i,i) = down(i)
   end do
   moment = [0.0_rp,0.0_rp,1.0_rp]
   dmoment = [1.0_rp,0.0_rp,0.0_rp]
   zero_moment = [0.0_rp,0.0_rp,0.0_rp]
   zero = cmplx(0.0_rp,0.0_rp,rp)
   call lmto_spin_parameter_tangent(up, down, dup, ddown, moment, dmoment, matrix, delta_matrix)
   call hcpx(up_matrix, 'cart2sph')
   call hcpx(down_matrix, 'cart2sph')
   call dresp08_rotate_collinear_operator(up_matrix, down_matrix, 1.0e-5_rp, plus)
   call dresp08_rotate_collinear_operator(up_matrix, down_matrix, -1.0e-5_rp, minus)
   error_spin = lmto_representation_relative_residual((plus-minus)/(2.0e-5_rp), delta_matrix)

   hhh = cmplx(0.0_rp,0.0_rp,rp)
   hcomm = cmplx(0.0_rp,0.0_rp,rp)
   do i = 1, norb
      hhh(i,i) = cmplx(0.17_rp+0.011_rp*real(i,rp),0.0_rp,rp)
      hcomm(i,i) = hhh(i,i)
      c0(i) = cmplx(0.03_rp+0.002_rp*real(i,rp),0.0_rp,rp)
      c1(i) = cmplx(0.06_rp+0.003_rp*real(i,rp),0.0_rp,rp)
      wx0i(i) = cmplx(0.83_rp+0.01_rp*real(i,rp),0.0_rp,rp)
      wx1i(i) = cmplx(0.12_rp+0.004_rp*real(i,rp),0.0_rp,rp)
      wx0j(i) = cmplx(0.91_rp-0.006_rp*real(i,rp),0.0_rp,rp)
      wx1j(i) = cmplx(0.08_rp+0.002_rp*real(i,rp),0.0_rp,rp)
      do j = i+1, norb
         hhh(i,j) = cmplx(0.004_rp*real(i+j,rp),0.001_rp*real(j-i,rp),rp)
         hhh(j,i) = conjg(hhh(i,j))
      end do
   end do
   call lmto_representation_bond_tangent(hhh, wx0i, wx1i, wx0j, wx1j, c0, c1, zero, zero, zero, zero, zero, zero, &
      moment, moment, dmoment, dmoment, .false., dbond)
   call build_bond(hhh, wx0i, wx1i, wx0j, wx1j, c0, c1, moment, bond_plus)
   call build_bond(hhh, wx0i, wx1i, wx0j, wx1j, c0, c1, rotate_moment(moment, 1.0e-5_rp), bond_plus)
   call build_bond(hhh, wx0i, wx1i, wx0j, wx1j, c0, c1, rotate_moment(moment, -1.0e-5_rp), bond_minus)
   error_bond = lmto_representation_relative_residual((bond_plus-bond_minus)/(2.0e-5_rp), dbond)

   ! Spin scalar/nonmagnetic negative control: all vector channels vanish.
   call lmto_representation_bond_tangent(hhh, wx0i, zero, wx0j, zero, c0, zero, zero, zero, zero, zero, zero, zero, &
      moment, moment, dmoment, dmoment, .false., dbond)
   error_nonmag = sqrt(sum(abs(dbond)**2))

   ! The same rotation oracle with a diagonal commuting structure block is an
   ! independent closed-form check of the endpoint product algebra.
   call lmto_representation_bond_tangent(hcomm, wx0i, wx1i, wx0j, wx1j, c0, c1, zero, zero, zero, zero, zero, zero, &
      moment, moment, dmoment, dmoment, .false., dbond)
   call build_bond(hcomm, wx0i, wx1i, wx0j, wx1j, c0, c1, rotate_moment(moment, 2.0e-6_rp), plus)
   call build_bond(hcomm, wx0i, wx1i, wx0j, wx1j, c0, c1, rotate_moment(moment, -2.0e-6_rp), minus)
   error_commuting = lmto_representation_relative_residual((plus-minus)/(4.0e-6_rp), dbond)

   ! DRESP-12 sitewise fixture: the production endpoint tangent is linear in
   ! the two endpoint rotations, so two independent site rotations must add
   ! without introducing a new radial solution or a fitted connection.
   call lmto_representation_bond_tangent(hhh, wx0i, wx1i, wx0j, wx1j, c0, c1, zero, zero, zero, zero, zero, zero, &
      moment, moment, dmoment, zero_moment, .false., bond_plus)
   call lmto_representation_bond_tangent(hhh, wx0i, wx1i, wx0j, wx1j, c0, c1, zero, zero, zero, zero, zero, zero, &
      moment, moment, zero_moment, dmoment, .false., bond_minus)
   call lmto_representation_bond_tangent(hhh, wx0i, wx1i, wx0j, wx1j, c0, c1, zero, zero, zero, zero, zero, zero, &
      moment, moment, dmoment, dmoment, .false., dbond)
   sitewise_error = lmto_representation_relative_residual(dbond, bond_plus + bond_minus)

   failed = maxval(error_fd) > 2.0e-8_rp .or. minval(ratio(2:4)) < 3.5_rp .or. error_spin > tol .or. &
      error_bond > 2.0e-9_rp .or. error_nonmag > tol .or. error_commuting > 2.0e-9_rp .or. sitewise_error > tol
   if (failed) then
      write (*,'(a)') 'UnitLrLmtoRepresentationTangent: FAIL'
      write (*,'(a,4es12.4)') '  predls_fd=', error_fd
      write (*,'(a,4es12.4)') '  predls_ratio=', ratio
      write (*,'(a,es12.4)') '  spin_matrix_fd=', error_spin
      write (*,'(a,es12.4)') '  bond_fd=', error_bond
      write (*,'(a,es12.4)') '  nonmagnetic=', error_nonmag
      write (*,'(a,es12.4)') '  commuting=', error_commuting
      write (*,'(a,es12.4)') '  sitewise_superposition=', sitewise_error
      error stop 1
   end if
   write (*,'(a,4es12.4)') '  predls_fd_center_shifted_width_obar=', error_fd
   write (*,'(a,4es12.4)') '  predls_theta2_ratio=', ratio
   write (*,'(a,es12.4)') '  spin_matrix_fd=', error_spin
   write (*,'(a,es12.4)') '  bond_fd=', error_bond
   write (*,'(a,es12.4)') '  nonmagnetic_zero_tangent=', error_nonmag
   write (*,'(a,es12.4)') '  commuting_structure=', error_commuting
   write (*,'(a,es12.4)') '  sitewise_superposition=', sitewise_error
   write (*,'(a)') 'UnitLrLmtoRepresentationTangent: PASS (analytic predls, spin lift, bond, commuting, nonmagnetic, sitewise controls)'

contains

   subroutine build_bond(hhh, wx0i, wx1i, wx0j, wx1j, c0, c1, moment, result)
      complex(rp), intent(in) :: hhh(:,:), wx0i(:), wx1i(:), wx0j(:), wx1j(:), c0(:), c1(:)
      real(rp), intent(in) :: moment(3)
      complex(rp), intent(out) :: result(:,:)
      integer :: idir
      call lmto_bond_value(hhh, wx0i, wx1i, wx0j, wx1j, c0, c1, moment, moment, .false., hhmag)
      do idir = 1, 4
         call hcpx(hhmag(:,:,idir), 'cart2sph')
      end do
      call lmto_hhmag_to_spinor(hhmag, result)
   end subroutine build_bond

   pure function rotate_moment(value, angle) result(rotated)
      real(rp), intent(in) :: value(3), angle
      real(rp) :: rotated(3)
      rotated = [sin(angle), 0.0_rp, cos(angle)]
   end function rotate_moment

end program test_lr_lmto_representation_tangent
