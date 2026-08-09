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
   use xc_response_kernel_mod, only: xc_response_sample, xc_response_kernel_provider
   implicit none

   real(rp), parameter :: tol = 1.0e-13_rp
   logical :: failed

   failed = .false.
   call test_pauli_and_circular_algebra()
   call test_two_level_spin_split_fixture()
   call test_xc_provider_contract()

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
         dbxc_dn=3.0_rp, dbxc_dm=4.0_rp, k_perp=5.0_rp)
      if (abs(provider%site(1)%bxc_energy + 0.10_rp) > tol .or. &
          abs(provider%site(1)%spin_population - 2.2_rp) > tol .or. &
          .not. provider%site(1)%has_k_perp .or. abs(provider%site(1)%k_perp - 5.0_rp) > tol) then
         write (*, '(a)') 'FAIL XC response-provider ground-state contract'
         failed = .true.
      end if
   end subroutine test_xc_provider_contract

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
