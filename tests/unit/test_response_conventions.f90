!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
!> @brief Pins response operators, circular normalization, and the XC-provider
!> data contract without implementing chi_KS or a TDDFT Dyson solve.
program test_response_conventions
   use precision_mod, only: rp
   use math_mod, only: i_unit, rho2nm
   use response_components_mod, only: RESPONSE_CHARGE, RESPONSE_MX, RESPONSE_MY, RESPONSE_MZ, &
      RESPONSE_PLUS, RESPONSE_MINUS
   use response_basis_mod, only: response_operator, ladder_operator
   use response_vertices_mod, only: response_channel
   use tddft_chi0_mod, only: tddft_chi0_options, tddft_chi0_result, build_chi_ks_from_eigenpairs
   use xc_response_kernel_mod, only: xc_response_sample, xc_response_kernel_provider, xc_response_radial_projection
   implicit none

   real(rp), parameter :: tol = 1.0e-13_rp
   logical :: failed

   failed = .false.
   call test_pauli_and_circular_algebra()
   call test_two_level_spin_split_fixture()
   call test_xc_provider_contract()
   call test_radial_alsda_projection()
   call test_transverse_goldstone_kernel_normalization()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine test_pauli_and_circular_algebra()
      complex(rp), allocatable :: s0(:, :), sx(:, :), sy(:, :), sz(:, :), sp(:, :), sm(:, :)
      complex(rp), allocatable :: sx_two_orb(:, :)
      complex(rp) :: comm(2, 2)

      s0 = response_operator(RESPONSE_CHARGE, 1)
      sx = response_operator(RESPONSE_MX, 1)
      sy = response_operator(RESPONSE_MY, 1)
      sz = response_operator(RESPONSE_MZ, 1)
      sp = ladder_operator(RESPONSE_PLUS, 1)
      sm = ladder_operator(RESPONSE_MINUS, 1)

      call assert_matrix('sigma_x sigma_y = i sigma_z', matmul(sx, sy), i_unit*sz)
      call assert_matrix('sigma_i squared = identity', matmul(sz, sz), s0)
      comm = matmul(sp, sm) - matmul(sm, sp)
      call assert_matrix('standard ladder commutator', comm, sz)
      call assert_matrix('m-plus is twice sigma-plus', response_operator(RESPONSE_PLUS, 1), 2.0_rp*sp)
      call assert_matrix('m-minus is twice sigma-minus', response_operator(RESPONSE_MINUS, 1), 2.0_rp*sm)

      ! Exact two-orbital ordering: (1-up, 2-up, 1-down, 2-down).
      sx_two_orb = response_operator(RESPONSE_MX, 2)
      call assert_scalar('spin-major first orbital up/down', sx_two_orb(1, 3), cmplx(1.0_rp, 0.0_rp, rp))
      call assert_scalar('spin-major second orbital up/down', sx_two_orb(2, 4), cmplx(1.0_rp, 0.0_rp, rp))
      call assert_scalar('no orbital-major spin mixing', sx_two_orb(1, 2), cmplx(0.0_rp, 0.0_rp, rp))
   end subroutine test_pauli_and_circular_algebra

   subroutine test_two_level_spin_split_fixture()
      real(rp), parameter :: e0 = 1.25_rp, b = 0.37_rp
      complex(rp), allocatable :: s0(:, :), sz(:, :), mplus(:, :), sigma_plus(:, :)
      complex(rp) :: h(2, 2), rho(2, 2), expectation
      real(rp) :: n, m(3)

      s0 = response_operator(RESPONSE_CHARGE, 1)
      sz = response_operator(RESPONSE_MZ, 1)
      mplus = response_operator(RESPONSE_PLUS, 1)
      sigma_plus = ladder_operator(RESPONSE_PLUS, 1)
      h = e0*s0 + b*sz
      call assert_scalar('H = H0 I + Hz sigma has E_up', h(1, 1), cmplx(e0 + b, 0.0_rp, rp))
      call assert_scalar('H = H0 I + Hz sigma has E_down', h(2, 2), cmplx(e0 - b, 0.0_rp, rp))
      call assert_scalar('m-plus spin-flip vertex', mplus(1, 2), cmplx(2.0_rp, 0.0_rp, rp))
      call assert_scalar('sigma-plus spin-flip vertex', sigma_plus(1, 2), cmplx(1.0_rp, 0.0_rp, rp))

      ! A coherent transverse one-electron density pins m+ = mx + i my.
      rho = reshape([cmplx(0.5_rp, 0.0_rp, rp), cmplx(0.5_rp, 0.0_rp, rp), &
                     cmplx(0.5_rp, 0.0_rp, rp), cmplx(0.5_rp, 0.0_rp, rp)], [2, 2])
      call rho2nm(rho, n, m)
      expectation = sum(rho*transpose(mplus))
      call assert_scalar('m-plus density normalization', expectation, cmplx(m(1), m(2), rp))
   end subroutine test_two_level_spin_split_fixture

   subroutine test_xc_provider_contract()
      type(xc_response_sample) :: sample
      type(xc_response_kernel_provider) :: provider

      sample%vxc_up = -0.80_rp
      sample%vxc_down = -0.60_rp
      sample%vxc_scalar = 0.5_rp*(sample%vxc_up + sample%vxc_down)
      sample%bxc_energy = 0.5_rp*(sample%vxc_up - sample%vxc_down)
      call provider%initialize(1, 'test-xc')
      call provider%record_ground_state_site(1, 2.2_rp, sample)
      call provider%set_site_derivatives(1, dvxc_dn=1.0_rp, dvxc_dm=2.0_rp, &
         dbxc_dn=3.0_rp, dbxc_dm=4.0_rp, k_perp_circular=5.0_rp)
      if (abs(provider%site(1)%bxc_energy + 0.10_rp) > tol .or. &
          abs(provider%site(1)%spin_population - 2.2_rp) > tol .or. &
          .not. provider%site(1)%has_k_perp_circular .or. abs(provider%site(1)%k_perp_circular - 5.0_rp) > tol) then
         write (*, '(a)') 'FAIL XC response-provider ground-state contract'
         failed = .true.
      end if
   end subroutine test_xc_provider_contract

   subroutine test_radial_alsda_projection()
      type(xc_response_kernel_provider) :: provider
      type(xc_response_radial_projection) :: projection

      ! m(r) = [2,1], B_xc(r) = [1,2]; all values are energy/population
      ! quantities in the SCF convention.  The response projector population
      ! is M_site=2.2, so K = (2 + 2)/(2*2.2^2) and B_site = 4/2.2.
      call projection%clear()
      call projection%accumulate(1.0_rp, rho_down=1.0_rp, rho_up=3.0_rp, &
         vxc_down=1.0_rp, vxc_up=3.0_rp)
      call projection%accumulate(1.0_rp, rho_down=2.0_rp, rho_up=3.0_rp, &
         vxc_down=0.0_rp, vxc_up=4.0_rp)
      call provider%initialize(1, 'test-xc')
      call provider%record_radial_projection(1, projection)
      call provider%set_site_spin_population(1, 2.2_rp)
      if (.not. provider%site(1)%has_k_perp_circular .or. &
          abs(provider%site(1)%radial_spin_population - 3.0_rp) > tol .or. &
          abs(provider%site(1)%bxc_energy - 4.0_rp/2.2_rp) > tol .or. &
          abs(provider%site(1)%k_perp_circular - 4.0_rp/(2.0_rp*2.2_rp*2.2_rp)) > tol .or. &
          abs(provider%site(1)%spin_population - 2.2_rp) > tol) then
         write (*, '(a)') 'FAIL radial ALSDA XC projection'
         failed = .true.
      end if
   end subroutine test_radial_alsda_projection

   ! This self-consistent two-level ferromagnet pins the complete transverse
   ! convention.  m+/- vertices carry 2, while the induced Hamiltonian field
   ! couples with one half; changing either the sign or the half breaks Xi=1.
   subroutine test_transverse_goldstone_kernel_normalization()
      type(xc_response_kernel_provider) :: provider
      type(xc_response_radial_projection) :: projection
      type(tddft_chi0_options) :: options
      type(tddft_chi0_result) :: response
      type(response_channel) :: left(1), right(1)
      real(rp) :: weights(1), eigenvalues(2, 1), omega(1)
      complex(rp) :: eigenvectors(2, 2, 1), xi

      ! H_xc = B sigma_z with B=-0.05 Ry: up is occupied, M=+1,
      ! chi_KS^{+-}(0)=-4/0.10=-40 Ry^-1 and K=B/(2M)=-0.025 Ry.
      call projection%clear()
      call projection%accumulate(1.0_rp, rho_down=0.0_rp, rho_up=1.0_rp, &
         vxc_down=0.05_rp, vxc_up=-0.05_rp)
      call provider%initialize(1, 'two-level-goldstone')
      call provider%record_radial_projection(1, projection)
      call provider%set_site_spin_population(1, 1.0_rp)
      call assert_scalar('transverse kernel includes signed circular half', &
         cmplx(provider%site(1)%k_perp_circular, 0.0_rp, rp), cmplx(-0.025_rp, 0.0_rp, rp))

      left(1) = response_channel(1, RESPONSE_PLUS)
      right(1) = response_channel(1, RESPONSE_MINUS)
      weights = 1.0_rp
      eigenvalues(:, 1) = [-0.05_rp, 0.05_rp]
      eigenvectors = cmplx(0.0_rp, 0.0_rp, rp)
      eigenvectors(1, 1, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      eigenvectors(2, 2, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      omega = 0.0_rp
      options%eta = 1.0e-15_rp
      options%fermi_level = 0.0_rp
      options%electronic_temperature = 0.0_rp
      call build_chi_ks_from_eigenpairs(weights, eigenvalues, eigenvectors, eigenvalues, eigenvectors, [1], &
         left, right, omega, options, response)
      xi = response%chi(1, 1, 1)*provider%site(1)%k_perp_circular
      call assert_scalar('two-level transverse Goldstone identity', xi, cmplx(1.0_rp, 0.0_rp, rp))
   end subroutine test_transverse_goldstone_kernel_normalization

   subroutine assert_matrix(label, actual, expected)
      character(len=*), intent(in) :: label
      complex(rp), intent(in) :: actual(:, :), expected(:, :)
      if (maxval(abs(actual - expected)) > tol) then
         write (*, '(a,1x,a)') 'FAIL', label
         failed = .true.
      end if
   end subroutine assert_matrix

   subroutine assert_scalar(label, actual, expected)
      character(len=*), intent(in) :: label
      complex(rp), intent(in) :: actual, expected
      if (abs(actual - expected) > tol) then
         write (*, '(a,1x,a)') 'FAIL', label
         failed = .true.
      end if
   end subroutine assert_scalar

end program test_response_conventions
