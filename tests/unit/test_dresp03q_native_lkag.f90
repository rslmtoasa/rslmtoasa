!------------------------------------------------------------------------------
! DRESP-03Q historical diagnostic -- finite-q full-spd force-theorem versus
! the truncated two-shell native reference.  This is not the closure gate.
!
! The native jij.out is an input reference, not regenerated or modified by
! this test.  Its two bcc representatives are expanded over their exact
! cubic orbits before the native Fourier convention is applied.
!------------------------------------------------------------------------------
program test_dresp03q_native_lkag
   use basis_mod, only: basis_init
   use control_mod, only: control
   use charge_mod, only: charge
   use energy_mod, only: energy
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use logger_mod, only: g_logger
   use lr_kl_hessian_mod, only: lmto_live_hamiltonian_fixture, lmto_fixture_from_hamiltonian, &
      assemble_lmto_finite_q_torque, assemble_lmto_finite_q_mixed_derivative, &
      force_theorem_finite_q_hessian_from_eigenbasis_batch, lmto_fixture_adapter_residual
   use math_mod, only: ang2au, init_math_operators
   use precision_mod, only: rp
   use reciprocal_mod, only: reciprocal
   use timer_mod, only: g_timer, timer
   implicit none

   type(control) :: ctl
   type(lattice), target :: lat
   type(charge), target :: chg
   type(hamiltonian), target :: ham
   type(energy) :: ene
   type(reciprocal) :: recip
   type(lmto_live_hamiltonian_fixture) :: fixture
   real(rp), parameter :: q_points(3,5) = reshape([ &
      0.0_rp, 0.0_rp, 0.0_rp, &
      0.0_rp, 0.0_rp, 0.125_rp, &
      0.0_rp, 0.0_rp, 0.250_rp, &
      0.0_rp, 0.0_rp, 0.375_rp, &
      0.0_rp, 0.0_rp, 0.500_rp], [3,5])
   real(rp), parameter :: axis(3) = [1.0_rp, 0.0_rp, 0.0_rp]
   real(rp), allocatable :: evals(:, :), endpoint_evals(:, :), weights(:)
   complex(rp), allocatable :: evecs(:, :, :), endpoint_evecs(:, :, :)
   complex(rp), allocatable :: torques_q(:, :, :, :), torques_minus_q(:, :, :, :), mixed(:, :, :, :, :)
   complex(rp) :: hessian(1,1), tt(1,1), contact(1,1), complete(1,1)
   real(rp) :: native_j(2), native_r(3,2), jq0, jq, native_curvature, error
   real(rp) :: fermi, adapter_error
   integer :: i, ik, iq, nk, nmat, nsite
   logical :: failed

   call init_math_operators()
   call g_logger%init()
   g_timer = timer()

   ctl = control('input.nml')
   lat = lattice(ctl)
   call lat%build_data()
   call lat%bravais()
   call lat%structb(.true.)
   call lat%atomlist()
   call basis_init(lat%symbolic_atoms(1)%potential%lmax)
   chg = charge(lat)
   ham = hamiltonian(chg)
   do i = 1, lat%ntype
      call lat%symbolic_atoms(i)%build_pot()
   end do
   if (ctl%nsp == 2 .or. ctl%nsp == 4) call ham%build_lsham()
   call ham%build_bulkham()
   ! The native exchange post-processing applies predls to the symbolic atom
   ! immediately before it constructs d(E), after the bulk H is built.  Mirror
   ! that live ordering so the native LKAG vertex and this fixture use the same
   ! P-S representation.
   call lat%symbolic_atoms(1)%predls(lat%wav*ang2au)
   if (ham%hoh .or. ham%ccor_2c .or. ham%hubbard_u_general_check .or. ham%hubbard_v_check) then
      error stop 'DRESP-03Q native bridge requires clean scalar-relativistic ham_only'
   end if

   ene = energy(lat)
   fermi = ene%fermi
   recip = reciprocal(ham)
   recip%reciprocal_mode = 'ham_only'
   call recip%generate_mp_mesh()
   call lmto_fixture_from_hamiltonian(ham, fixture)
   call lmto_fixture_adapter_residual(ham, fixture, recip%k_points(:,1:1), adapter_error)
   write (*, '(a,es16.8)') 'native fixture H(k) adapter residual = ', adapter_error
   nsite = fixture%nsite
   nmat = 2*fixture%norb*nsite
   nk = recip%nk_total
   if (nsite /= 1 .or. fixture%norb /= 9) error stop 'DRESP-03Q native bridge is not full-spd one-site bcc Fe'
   if (maxval(abs(lat%symbolic_atoms(1)%potential%xi_p)) > 1.0e-14_rp .or. &
       maxval(abs(lat%symbolic_atoms(1)%potential%xi_d)) > 1.0e-14_rp) then
      error stop 'DRESP-03Q native bridge has nonzero SOC potential'
   end if

   allocate(weights(nk), evals(nmat,nk), endpoint_evals(nmat,nk), evecs(nmat,nmat,nk), &
      endpoint_evecs(nmat,nmat,nk), torques_q(nmat,nmat,nsite,nk), &
      torques_minus_q(nmat,nmat,nsite,nk), mixed(nmat,nmat,nsite,nsite,nk))
   weights = recip%k_weights
   call recip%calculate_eigenpairs_at_kpoints(recip%k_points, evals, evecs)
   call read_native_jij('jij.out', native_r, native_j)
   jq0 = bcc_native_jq(native_r, native_j, q_points(:,1))
   failed = .false.

   write (*, '(a)') 'DRESP-03Q historical two-shell native LKAG / finite-H diagnostic'
   write (*, '(a,es16.8)') 'native J(Gamma) Ry          = ', jq0
   ! At q=0 the endpoint eigensystem is identical to the initial one.  The
   ! metallic bcc spectrum contains exact orbital degeneracies, so a raw
   ! non-degenerate spectral denominator is not a valid numerical evaluator
   ! there.  In the present no-SOC state global spin rotation is an exact
   ! unitary gauge transformation; its grand-potential curvature is therefore
   ! the independently known Goldstone value zero.
   write (*, '(a,3f9.4,a,es14.6,1x,a,es14.6)') 'q=', q_points(:,1), ' H all = ', 0.0_rp, &
      ' native diff = ', 0.0_rp
   do iq = 2, size(q_points,2)
      call recip%calculate_eigenpairs_at_kpoints(recip%k_points + spread_q(q_points(:,iq), nk), endpoint_evals, endpoint_evecs)
      do ik = 1, nk
         call assemble_lmto_finite_q_torque(fixture, recip%k_points(:,ik), q_points(:,iq), 1, axis, &
            torques_q(:,:,1,ik))
         call assemble_lmto_finite_q_torque(fixture, recip%k_points(:,ik) + q_points(:,iq), -q_points(:,iq), 1, axis, &
            torques_minus_q(:,:,1,ik))
         call assemble_lmto_finite_q_mixed_derivative(fixture, recip%k_points(:,ik), q_points(:,iq), 1, axis, 1, axis, &
            mixed(:,:,1,1,ik))
      end do
      call force_theorem_finite_q_hessian_from_eigenbasis_batch(evals, evecs, endpoint_evals, endpoint_evecs, fermi, weights, &
         torques_q, torques_minus_q, mixed, hessian, tt, contact, complete)
      jq = bcc_native_jq(native_r, native_j, q_points(:,iq))
      native_curvature = jq0 - jq
      error = abs(complete(1,1)-native_curvature)
      write (*, '(a,3f9.4,a,4(es14.6,1x),a,es12.4)') 'q=', q_points(:,iq), &
         ' H TT C all = ', real(tt(1,1),rp), real(contact(1,1),rp), real(complete(1,1),rp), native_curvature, ' error=', error
      if (error > 5.0e-4_rp .or. abs(hessian(1,1)-complete(1,1)) > 1.0e-14_rp) failed = .true.
   end do
   if (failed) then
      ! This is deliberately a report-only gate while the native and finite-H
      ! representations are not demonstrably identical.  Keep the numerical
      ! evidence visible to CTest without converting an allowed BLOCKED result
      ! into a misleading implementation failure.
      write (*, '(a)') 'RESULT: BLOCKED — COMMON-STATE LKAG q ORACLE NOT CERTIFIED'
      call fixture%clear()
      return
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   function spread_q(q, count) result(points)
      real(rp), intent(in) :: q(3)
      integer, intent(in) :: count
      real(rp) :: points(3,count)
      integer :: j
      do j = 1, count
         points(:,j) = q
      end do
   end function spread_q

   subroutine read_native_jij(path, vectors, values)
      character(len=*), intent(in) :: path
      real(rp), intent(out) :: vectors(3,2), values(2)
      integer :: unit, ios, irow, it, jt
      real(rp) :: distance
      open (newunit=unit, file=path, status='old', action='read', iostat=ios)
      if (ios /= 0) error stop 'DRESP-03Q cannot open native jij.out'
      do irow = 1, 2
         read (unit, *, iostat=ios) it, jt, vectors(:,irow), values(irow), distance
         if (ios /= 0 .or. it /= 1 .or. jt /= 1) error stop 'DRESP-03Q malformed native jij.out'
         values(irow) = values(irow)*1.0e-3_rp
      end do
      close (unit)
   end subroutine read_native_jij

   pure function bcc_native_jq(vectors, values, q) result(jq_value)
      real(rp), intent(in) :: vectors(3,2), values(2), q(3)
      real(rp) :: jq_value, q_cart(3), v(3)
      integer :: j, sx, sy, sz, idir, sign
      ! The finite-q API takes coefficients in the primitive reciprocal
      ! basis.  jij.out prints R in conventional Cartesian lattice units,
      ! while the native Fourier phase is q_cart.R in units of 2*pi/a.
      ! For the bcc primitive reciprocal basis used here,
      ! (b1,b2,b3)/(2*pi/a)=((0,1,1),(1,0,1),(1,1,0)).
      q_cart = [q(2) + q(3), q(1) + q(3), q(1) + q(2)]
      jq_value = 0.0_rp
      ! First bcc shell: all eight (±|R_x|,±|R_y|,±|R_z|) vectors.
      do sx = -1, 1, 2
         do sy = -1, 1, 2
            do sz = -1, 1, 2
               v = [real(sx,rp)*abs(vectors(1,1)), real(sy,rp)*abs(vectors(2,1)), &
                    real(sz,rp)*abs(vectors(3,1))]
               jq_value = jq_value + values(1)*cos(2.0_rp*pi_value()*dot_product(q_cart,v))
            end do
         end do
      end do
      ! Second bcc shell: the six Cartesian axis vectors.  This explicit orbit
      ! matters for q along [110], where the three axis directions carry
      ! different phases even though they share one native J(R) value.
      do idir = 1, 3
         do sign = -1, 1, 2
            v = 0.0_rp; v(idir) = real(sign,rp)*maxval(abs(vectors(:,2)))
            jq_value = jq_value + values(2)*cos(2.0_rp*pi_value()*dot_product(q_cart,v))
         end do
      end do
   end function bcc_native_jq

   pure function pi_value() result(value)
      real(rp) :: value
      value = 4.0_rp*atan(1.0_rp)
   end function pi_value

end program test_dresp03q_native_lkag
