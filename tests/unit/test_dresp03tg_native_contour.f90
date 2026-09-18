! DRESP-03TG-CLOSE -- native Turek contour and Fourier closure gates.
program test_dresp03tg_native_contour
   use basis_mod, only: basis_init
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use charge_mod, only: charge
   use control_mod, only: control
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use logger_mod, only: g_logger
   use lr_lmto_turek_contour_mod, only: native_turek_contour_options, native_turek_contour_report, &
      native_build_contour, native_exchange_q_contour, native_exchange_jij_contour, &
      native_fourier_jq_to_jij, native_fourier_jij_to_jq
   use math_mod, only: ang2au, init_math_operators
   use precision_mod, only: rp
   use timer_mod, only: g_timer, timer
   implicit none

   type(control) :: ctl
   type(lattice), target :: lat
   type(charge), target :: chg
   type(hamiltonian), target :: ham
   type(native_turek_contour_options) :: options
   type(native_turek_contour_report) :: report
   complex(rp), allocatable :: nodes(:), weights(:), poles(:)
   real(rp), allocatable :: jq(:, :, :), q_points(:, :), q_mesh(:, :), k_points(:, :), k_weights(:)
   real(rp), allocatable :: real_space(:, :), jij(:, :, :), jij_roundtrip(:, :, :), jq_roundtrip(:, :, :)
   real(rp), allocatable :: native_jq(:, :, :), native_jij(:, :, :), native_jq_roundtrip(:, :, :)
   real(rp) :: bounds(2), contour_weight_error, fourier_error, native_fourier_error, q_value_error
   integer :: i, ix, iy, iz, iq, ik, ir

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
   call ham%build_bulkham()
   do i = 1, lat%ntype
      call lat%symbolic_atoms(i)%predls(lat%wav*ang2au)
   end do
   if (lat%nrec /= 1) error stop 'DRESP-03TG-CLOSE contour fixture requires one reciprocal site'

   options%contour_points = 32
   options%contour_shape = 'ellipse'
   options%contour_margin = 0.25_rp
   options%contour_height_fraction = 0.35_rp
   options%account_fermi_poles = .true.
   ! The native path poles are the full accepted spd spectrum; the compact
   ! atomic-input DOS window is not a safe contour bound for this gate.
   bounds = [-3.0_rp, 3.0_rp]
   call native_build_contour(bounds, -0.070393_rp, 300.0_rp*6.3336814e-6_rp, options, nodes, weights, poles)
   contour_weight_error = abs(sum(weights))
   if (contour_weight_error > 2.0e-13_rp .or. size(nodes) /= 32 .or. size(poles) < 1) then
      error stop 'DRESP-03TG-CLOSE native contour construction failed'
   end if

   allocate(k_points(3,1), k_weights(1), q_points(3,3), jq(1,1,3))
   k_points(:,1) = 0.0_rp; k_weights(1) = 1.0_rp
   q_points = 0.0_rp
   q_points(1,2) = 0.25_rp
   q_points(1,3) = -0.25_rp
   call native_exchange_q_contour(lat, k_points, k_weights, q_points, -0.070393_rp, 300.0_rp*6.3336814e-6_rp, &
      bounds, options, jq, report)
   q_value_error = 0.0_rp
   do iq = 1, 3
      if (.not. ieee_is_finite(jq(1,1,iq))) q_value_error = huge(1.0_rp)
   end do
   if (report%contour_points /= 32*3 .or. report%fermi_poles < 1 .or. q_value_error >= huge(1.0_rp)) then
      error stop 'DRESP-03TG-CLOSE native contour evaluation failed'
   end if
   deallocate(q_points)
   allocate(q_mesh(3,8), real_space(3,8), jij(1,1,8), jij_roundtrip(1,1,8), jq_roundtrip(1,1,8), &
      native_jq(1,1,8), native_jij(1,1,8), native_jq_roundtrip(1,1,8))
   iq = 0
   do iz = 0, 1
      do iy = 0, 1
         do ix = 0, 1
            iq = iq+1
            q_mesh(:,iq) = [real(ix,rp)/2.0_rp, real(iy,rp)/2.0_rp, real(iz,rp)/2.0_rp]
            real_space(:,iq) = [real(ix,rp), real(iy,rp), real(iz,rp)]
            jij(1,1,iq) = 0.1_rp*real(iq,rp)
         end do
      end do
   end do
   call native_fourier_jij_to_jq(q_mesh, real_space, jij, jq_roundtrip)
   call native_fourier_jq_to_jij(q_mesh, jq_roundtrip, real_space, jij_roundtrip)
   fourier_error = maxval(abs(jij_roundtrip-jij))
   if (fourier_error > 3.0e-13_rp) error stop 'DRESP-03TG-CLOSE Jij/J(q) Fourier identity failed'

   ! Close a complete live native q mesh as well.  This exercises the
   ! production contour result, not only the algebraic Fourier helper.
   call native_exchange_jij_contour(lat, k_points, k_weights, q_mesh, real_space, -0.070393_rp, &
      300.0_rp*6.3336814e-6_rp, bounds, options, native_jq, native_jij, report)
   call native_fourier_jij_to_jq(q_mesh, real_space, native_jij, native_jq_roundtrip)
   native_fourier_error = maxval(abs(native_jq_roundtrip-native_jq))
   if (native_fourier_error > 3.0e-13_rp) error stop 'DRESP-03TG-CLOSE native Jij/J(q) closure failed'

   write (*,'(a,es14.6)') 'DRESP-03TG-CLOSE contour weight residual = ', contour_weight_error
   write (*,'(a,i0,a,i0)') 'DRESP-03TG-CLOSE contour nodes/poles = ', report%contour_points, '/', report%fermi_poles
   write (*,'(a,3(es14.6,1x))') 'DRESP-03TG-CLOSE native J(q) = ', jq(1,1,1), jq(1,1,2), jq(1,1,3)
   write (*,'(a,es14.6)') 'DRESP-03TG-CLOSE Fourier roundtrip residual = ', fourier_error
   write (*,'(a,es14.6)') 'DRESP-03TG-CLOSE native Jij/J(q) residual = ', native_fourier_error
   write (*,'(a)') 'DRESP-03TG-CLOSE: PASS'

   deallocate(nodes,weights,poles,k_points,k_weights,q_mesh,real_space,jij,jij_roundtrip,jq,jq_roundtrip, &
      native_jq,native_jij,native_jq_roundtrip)
end program test_dresp03tg_native_contour
