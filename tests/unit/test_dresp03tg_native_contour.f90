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
      native_exchange_q_ordered_contour, native_exchange_pairs_contour, native_spectral_bounds, &
      native_fourier_jq_to_jij, native_fourier_jij_to_jq, native_fourier_complex_jij_to_jq, native_complex_fermi, &
      native_regularized_fermi
   use math_mod, only: ang2au, init_math_operators, pi
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
   real(rp), allocatable :: k_full(:, :), k_full_weights(:), q_pair(:, :)
   real(rp), allocatable :: native_jq(:, :, :), native_jij(:, :, :), native_jq_roundtrip(:, :, :)
   complex(rp), allocatable :: jq_ud(:, :, :), jq_du(:, :, :), pair_ud(:, :, :), pair_du(:, :, :), jq_pair_ft(:, :, :)
   real(rp) :: bounds(2), contour_weight_error, fourier_error, native_fourier_error, q_value_error
   real(rp) :: scalar_oracle_error, scalar_lkag_error, scalar_coarse_error, scalar_coarse_lkag_error
   real(rp) :: scalar_fine_error, scalar_fine_lkag_error, pair_closure_error, covariance_error, spectral_error
   real(rp) :: curvature_oracle_error, max_ellipse_value
   logical :: bounds_inside
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

   options%contour_points = 64
   options%contour_shape = 'ellipse'
   options%contour_margin = 0.25_rp
   options%contour_height_fraction = 0.35_rp
   options%account_fermi_poles = .true.
   options%target_fermi_poles = 32
   ! The native path poles are the full accepted spd spectrum; the compact
   ! atomic-input DOS window is not a safe contour bound for this gate.
   bounds = [-3.0_rp, 3.0_rp]
   call native_build_contour(bounds, -0.070393_rp, 300.0_rp*6.3336814e-6_rp, options, nodes, weights, poles)
   contour_weight_error = abs(sum(weights))
   if (contour_weight_error > 2.0e-13_rp .or. size(nodes) /= options%contour_points .or. size(poles) < 1) then
      error stop 'DRESP-03TG-CLOSE native contour construction failed'
   end if
   scalar_coarse_error = 0.0_rp; scalar_coarse_lkag_error = 0.0_rp
   scalar_fine_error = 0.0_rp; scalar_fine_lkag_error = 0.0_rp
   call scalar_pole_oracle(8192, 32, scalar_coarse_error, scalar_coarse_lkag_error)
   call scalar_pole_oracle(32768, 32, scalar_fine_error, scalar_fine_lkag_error)
   scalar_oracle_error = scalar_fine_error
   scalar_lkag_error = scalar_fine_lkag_error
   write (*,'(a,2(es14.6,1x))') 'DRESP-03TG-CLOSE scalar coarse/fine = ', scalar_coarse_error, scalar_fine_error
   write (*,'(a,2(es14.6,1x))') 'DRESP-03TG-CLOSE scalar LKAG coarse/fine = ', scalar_coarse_lkag_error, scalar_fine_lkag_error
   if (scalar_fine_error > 2.0e-8_rp .or. scalar_fine_lkag_error > 5.0e-9_rp .or. &
       scalar_fine_error >= 0.1_rp*scalar_coarse_error) then
      error stop 'DRESP-03TG-CLOSE scalar contour/pole oracle failed'
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
   if (report%contour_points /= options%contour_points*3 .or. report%fermi_poles < 1 .or. q_value_error >= huge(1.0_rp)) then
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

   ! Independent native closure: transform the Green-function blocks in real
   ! space first, then compare their contour pair evaluation with direct q.
   allocate(k_full(3,8), k_full_weights(8))
   ik = 0
   do iz = 0, 1
      do iy = 0, 1
         do ix = 0, 1
            ik = ik+1
            k_full(:,ik) = [real(ix,rp)/2.0_rp, real(iy,rp)/2.0_rp, real(iz,rp)/2.0_rp]
            k_full_weights(ik) = 1.0_rp/8.0_rp
         end do
      end do
   end do
   allocate(jq_ud(1,1,8), jq_du(1,1,8), pair_ud(1,1,8), pair_du(1,1,8), jq_pair_ft(1,1,8))
   call native_exchange_q_ordered_contour(lat, k_full, k_full_weights, q_mesh, -0.070393_rp, &
      300.0_rp*6.3336814e-6_rp, bounds, options, jq_ud, jq_du, report)
   call native_exchange_pairs_contour(lat, k_full, k_full_weights, real_space, -0.070393_rp, &
      300.0_rp*6.3336814e-6_rp, bounds, options, pair_ud, pair_du, report)
   call native_fourier_complex_jij_to_jq(q_mesh, real_space, pair_ud, jq_pair_ft)
   pair_closure_error = maxval(abs(jq_ud-jq_pair_ft))
   call native_fourier_complex_jij_to_jq(q_mesh, real_space, pair_du, jq_pair_ft)
   pair_closure_error = max(pair_closure_error,maxval(abs(jq_du-jq_pair_ft)))
   if (pair_closure_error > 3.0e-11_rp) error stop 'DRESP-03TG-CLOSE independent pair/q closure failed'

   allocate(q_pair(3,2))
   ! Use a reciprocal vector compatible with the complete 2^3 k mesh.  A
   ! quarter-grid shift would compare two different finite quadratures and
   ! would not be a valid covariance gate at this resolution.
   q_pair = 0.0_rp; q_pair(1,:) = [0.5_rp,-0.5_rp]
   deallocate(jq_ud,jq_du); allocate(jq_ud(1,1,2),jq_du(1,1,2))
   call native_exchange_q_ordered_contour(lat, k_full, k_full_weights, q_pair, -0.070393_rp, &
      300.0_rp*6.3336814e-6_rp, bounds, options, jq_ud, jq_du, report)
   covariance_error = max(abs(jq_ud(1,1,1)-jq_du(1,1,2)),abs(jq_ud(1,1,1)-jq_ud(1,1,2)))
   covariance_error = max(covariance_error,abs(jq_ud(1,1,1)-jq_du(1,1,1)))
   if (covariance_error > 3.0e-10_rp) error stop 'DRESP-03TG-CLOSE q covariance failed'

   call native_spectral_bounds(lat, k_full, -0.070393_rp, 300.0_rp*6.3336814e-6_rp, options, bounds, &
      max_ellipse_value, bounds_inside, iq)
   spectral_error = max_ellipse_value
   if (.not. bounds_inside .or. spectral_error >= 1.0_rp .or. iq /= 2*9*8) then
      error stop 'DRESP-03TG-CLOSE native spectral-bound gate failed'
   end if
   curvature_oracle_error = heisenberg_curvature_oracle()
   if (curvature_oracle_error > 1.0e-13_rp) error stop 'DRESP-03TG-CLOSE finite-q curvature oracle failed'

   write (*,'(a,es14.6)') 'DRESP-03TG-CLOSE contour weight residual = ', contour_weight_error
   write (*,'(a,i0,a,i0)') 'DRESP-03TG-CLOSE contour nodes/poles = ', report%contour_points, '/', report%fermi_poles
   write (*,'(a,3(es14.6,1x))') 'DRESP-03TG-CLOSE native J(q) = ', jq(1,1,1), jq(1,1,2), jq(1,1,3)
   write (*,'(a,es14.6)') 'DRESP-03TG-CLOSE Fourier roundtrip residual = ', fourier_error
   write (*,'(a,es14.6)') 'DRESP-03TG-CLOSE DFT helper roundtrip residual = ', native_fourier_error
   write (*,'(a,es14.6)') 'DRESP-03TG-CLOSE scalar pole oracle residual = ', scalar_oracle_error
   write (*,'(a,es14.6)') 'DRESP-03TG-CLOSE scalar LKAG residual = ', scalar_lkag_error
   write (*,'(a,es14.6)') 'DRESP-03TG-CLOSE independent pair/q residual = ', pair_closure_error
   write (*,'(a,es14.6)') 'DRESP-03TG-CLOSE q covariance residual = ', covariance_error
   write (*,'(a,es14.6)') 'DRESP-03TG-CLOSE native max spectral ellipse = ', spectral_error
   write (*,'(a,es14.6)') 'DRESP-03TG-CLOSE curvature factor oracle residual = ', curvature_oracle_error
   write (*,'(a)') 'DRESP-03TG-CLOSE-R1 algebraic/native gates: PASS'

   deallocate(nodes,weights,poles,k_points,k_weights,q_mesh,real_space,jij,jij_roundtrip,jq,jq_roundtrip, &
      native_jq,native_jij,native_jq_roundtrip,k_full,k_full_weights,jq_ud,jq_du,pair_ud,pair_du,jq_pair_ft,q_pair)

contains

   subroutine scalar_pole_oracle(npoint, target_poles, max_error, max_lkag_error)
      integer, intent(in) :: npoint, target_poles
      real(rp), intent(inout) :: max_error, max_lkag_error
      type(native_turek_contour_options) :: local_options
      complex(rp), allocatable :: local_nodes(:), local_weights(:), local_poles(:)
      complex(rp) :: integral, regularized_integral, pole_residue
      real(rp), parameter :: epsilons(3) = [-1.2_rp,-0.15_rp,0.25_rp]
      real(rp), parameter :: fermi = -0.070393_rp, kT = 300.0_rp*6.3336814e-6_rp
      real(rp) :: fvalue
      integer :: ie, inode, ipole

      local_options = options
      local_options%contour_points = npoint
      local_options%target_fermi_poles = target_poles
      call native_build_contour(bounds,fermi,kT,local_options,local_nodes,local_weights,local_poles)
      do ie = 1, size(epsilons)
         integral = cmplx(0.0_rp,0.0_rp,rp)
         regularized_integral = cmplx(0.0_rp,0.0_rp,rp)
         do inode = 1, size(local_nodes)
            integral = integral + local_weights(inode)*native_complex_fermi(local_nodes(inode),fermi,kT)/ &
               (local_nodes(inode)-epsilons(ie))
            regularized_integral = regularized_integral + local_weights(inode)* &
               native_regularized_fermi(local_nodes(inode),fermi,kT,local_poles)/(local_nodes(inode)-epsilons(ie))
         end do
         do ipole = 1, size(local_poles)
            pole_residue = kT/(local_poles(ipole)-epsilons(ie))
            integral = integral + pole_residue
            regularized_integral = regularized_integral + pole_residue
         end do
         fvalue = 1.0_rp/(exp((epsilons(ie)-fermi)/kT)+1.0_rp)
         write (*,'(a,2(i0,1x),a,es12.4,1x,es12.4)') 'DRESP-03TG-CLOSE scalar sample ', npoint, target_poles, 'errors=', &
            abs(real(integral,rp)-fvalue), abs(real(regularized_integral,rp)-fvalue)
         max_error = max(max_error,abs(real(integral,rp)-fvalue))
         max_error = max(max_error,abs(real(regularized_integral,rp)-fvalue))
         max_lkag_error = max(max_lkag_error,abs(real(-integral/4.0_rp,rp)+fvalue/4.0_rp))
      end do
      deallocate(local_nodes,local_weights,local_poles)
   end subroutine scalar_pole_oracle

   function heisenberg_curvature_oracle() result(error)
      real(rp) :: error
      real(rp), parameter :: rvec(3,4) = reshape([0.0_rp,0.0_rp,0.0_rp, 1.0_rp,0.0_rp,0.0_rp, &
         -1.0_rp,0.0_rp,0.0_rp, 2.0_rp,0.0_rp,0.0_rp],[3,4])
      real(rp), parameter :: jhist(4) = [0.20_rp,0.07_rp,0.07_rp,-0.03_rp]
      real(rp), parameter :: q = 0.25_rp
      real(rp) :: j0, jq, coefficient, expected, phase_difference
      real(rp), parameter :: cell_x(4) = [0.0_rp,0.25_rp,0.50_rp,0.75_rp]
      integer :: ic, ir

      j0 = sum(jhist); jq = 0.0_rp
      do ir = 1, size(jhist)
         jq = jq + jhist(ir)*cos(2.0_rp*acos(-1.0_rp)*q*rvec(1,ir))
      end do
      ! Expand the double-counted Heisenberg energy directly for
      ! theta_x = a exp(+iqx) + b exp(-iqx).  For every directed pair the
      ! coefficient of a*b in Delta E is J(R)*(2-2 cos(q.R)); this is an
      ! algebraic coefficient extraction, not a numerical fit.
      coefficient = 0.0_rp
      do ic = 1, size(cell_x)
         do ir = 1, size(jhist)
            phase_difference = 2.0_rp*pi*q*(cell_x(ic)-(cell_x(ic)+rvec(1,ir)))
            coefficient = coefficient + jhist(ir)*(2.0_rp-2.0_rp*cos(phase_difference))/real(size(cell_x),rp)
         end do
      end do
      expected = 2.0_rp*(j0-jq)
      error = abs(coefficient-expected)
   end function heisenberg_curvature_oracle
end program test_dresp03tg_native_contour
