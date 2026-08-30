! RS-LMTO-ASA -- WP05 rotating/lab-frame covariance oracles
program test_gbt_frame_contract
   use precision_mod, only: rp
   use gbt_structure_mod, only: gbt_frame_t, gbt_frame_from_phase, gbt_frame_for_cell, gbt_frame_link, &
                                gbt_endpoint_link, gbt_rotating_to_lab_vector, gbt_lab_to_rotating_vector, &
                                gbt_rotating_to_lab_spinor, gbt_lab_to_rotating_spinor, &
                                gbt_rotating_to_lab_density, gbt_lab_to_rotating_density, &
                                gbt_reference_moment_is_collinear
   use spin_density_mod, only: sd_matrix_from_cartesian, sd_cartesian_from_matrix
   implicit none

   real(rp), parameter :: tol = 2.0e-12_rp
   real(rp), parameter :: pi = acos(-1.0_rp)
   real(rp) :: max_error
   integer :: failures

   max_error = 0.0_rp
   failures = 0
   call test_su2_so3_and_density()
   call test_commensurate_spiral_sites()
   call test_frame_links_and_reference_guard()

   write (*, '(a,es12.4)') 'GBT frame-contract maximum error: ', max_error
   if (failures > 0 .or. max_error >= tol) then
      write (*, '(a,i0)') 'RESULT: FAIL (checks failed: ', failures
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine check(value, label)
      real(rp), intent(in) :: value
      character(len=*), intent(in) :: label

      max_error = max(max_error, value)
      if (value >= tol) then
         failures = failures + 1
         write (*, '(a,a,a,es12.4)') 'FAIL ', trim(label), ': ', value
      end if
   end subroutine check

   subroutine check_true(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label

      if (.not. condition) then
         failures = failures + 1
         write (*, '(a,a)') 'FAIL ', trim(label)
      end if
   end subroutine check_true

   subroutine test_su2_so3_and_density()
      type(gbt_frame_t) :: frame
      real(rp) :: m_rot(3), m_lab(3), m_back(3), m_from_density(3), n, n_lab
      complex(rp) :: rho_rot(2, 2), rho_lab(2, 2), rho_back(2, 2)
      complex(rp) :: psi_rot(2), psi_lab(2), psi_back(2), psi_density(2, 2)
      complex(rp) :: expected_density(2, 2)
      integer :: i, j

      call gbt_frame_from_phase(0.71_rp, 0.37_rp, -0.83_rp, frame)
      m_rot = [0.31_rp, -0.27_rp, 0.84_rp]
      call gbt_rotating_to_lab_vector(frame, m_rot, m_lab)
      call gbt_lab_to_rotating_vector(frame, m_lab, m_back)
      call check(maxval(abs(m_back - m_rot)), 'SO(3) vector inverse')
      call check(maxval(abs(matmul(transpose(frame%R), frame%R) - identity_real())), 'SO(3) orthogonality')

      call sd_matrix_from_cartesian(1.7_rp, m_rot, rho_rot)
      call gbt_rotating_to_lab_density(frame, rho_rot, rho_lab)
      call sd_cartesian_from_matrix(rho_lab, n_lab, m_from_density)
      call check(abs(n_lab - 1.7_rp), 'density charge invariance')
      call check(maxval(abs(m_from_density - m_lab)), 'SU(2)/SO(3) density covariance')
      call check(abs(sqrt(sum(m_from_density**2)) - sqrt(sum(m_rot**2))), 'density moment magnitude')
      call gbt_lab_to_rotating_density(frame, rho_lab, rho_back)
      call check(maxval(abs(rho_back - rho_rot)), 'density matrix inverse')

      psi_rot = [cmplx(0.37_rp, -0.22_rp, rp), cmplx(-0.19_rp, 0.41_rp, rp)]
      call gbt_rotating_to_lab_spinor(frame, psi_rot, psi_lab)
      call gbt_lab_to_rotating_spinor(frame, psi_lab, psi_back)
      call check(maxval(abs(psi_back - psi_rot)), 'SU(2) spinor inverse')
      do i = 1, 2
         do j = 1, 2
            psi_density(i, j) = psi_rot(i)*conjg(psi_rot(j))
         end do
      end do
      call gbt_rotating_to_lab_density(frame, psi_density, rho_lab)
      do i = 1, 2
         do j = 1, 2
            expected_density(i, j) = psi_lab(i)*conjg(psi_lab(j))
         end do
      end do
      call check(maxval(abs(rho_lab - expected_density)), 'spinor/density transform agreement')
      n = real(rho_rot(1, 1) + rho_rot(2, 2), rp)
      call check(abs(n - 1.7_rp), 'rotating density trace')
   end subroutine test_su2_so3_and_density

   subroutine test_commensurate_spiral_sites()
      type(gbt_frame_t) :: frame
      real(rp) :: q(3), r(3), m_rot(3), m_lab(3), expected(3), phase
      real(rp) :: theta(2), phi(2)
      integer :: isub, n

      q = [0.25_rp, 0.0_rp, 0.0_rp]
      theta = [0.62_rp, 1.11_rp]
      phi = [0.18_rp, -0.47_rp]
      m_rot = [0.0_rp, 0.0_rp, 1.0_rp]

      ! This is the site-by-site oracle for the explicit period-four
      ! supercell: each sublattice has its own reference phase offset, while
      ! translation advances the lab azimuth by q dot R exactly once.
      do isub = 1, 2
         do n = 0, 3
            r = [real(n, rp)*4.0_rp, 0.0_rp, 0.0_rp]
            phase = 2.0_rp*pi*q(1)*real(n, rp)
            call gbt_frame_for_cell(q, r, 4.0_rp, theta(isub), phi(isub), frame)
            call gbt_rotating_to_lab_vector(frame, m_rot, m_lab)
            expected = [sin(theta(isub))*cos(phi(isub) + phase), &
                        sin(theta(isub))*sin(phi(isub) + phase), cos(theta(isub))]
            call check(maxval(abs(m_lab - expected)), 'commensurate lab moment site')
            call check(abs(sqrt(sum(m_lab**2)) - 1.0_rp), 'commensurate moment norm')
         end do
      end do
   end subroutine test_commensurate_spiral_sites

   subroutine test_frame_links_and_reference_guard()
      type(gbt_frame_t) :: frame_a, frame_b
      complex(rp) :: link(2, 2), endpoint(2, 2), reverse(2, 2)
      real(rp) :: q(3), d(3)

      q = [0.19_rp, -0.07_rp, 0.11_rp]
      d = [1.3_rp, -0.8_rp, 0.6_rp]
      call gbt_frame_from_phase(0.41_rp, -0.26_rp, 0.0_rp, frame_a)
      call gbt_frame_from_phase(0.93_rp, 0.58_rp, 2.0_rp*pi*0.19_rp*1.3_rp + &
                                2.0_rp*pi*(-0.07_rp)*(-0.8_rp) + 2.0_rp*pi*0.11_rp*0.6_rp, frame_b)
      call gbt_frame_link(frame_a, frame_b, link)
      call gbt_endpoint_link(q, d, 1.0_rp, 0.41_rp, -0.26_rp, 0.93_rp, 0.58_rp, endpoint)
      call check(maxval(abs(link - endpoint)), 'frame link/endpoint-link agreement')
      call gbt_frame_link(frame_b, frame_a, reverse)
      call check(maxval(abs(reverse - transpose(conjg(link)))), 'frame link reversal')

      call check_true(gbt_reference_moment_is_collinear([0.0_rp, 0.0_rp, -1.0_rp], 1.0e-10_rp), &
                      'collinear rotating reference accepted')
      call check_true(.not. gbt_reference_moment_is_collinear([1.0e-3_rp, 0.0_rp, 1.0_rp], 1.0e-10_rp), &
                      'transverse rotating reference rejected')
   end subroutine test_frame_links_and_reference_guard

   pure function identity_real() result(eye)
      real(rp) :: eye(3, 3)

      eye = 0.0_rp
      eye(1, 1) = 1.0_rp
      eye(2, 2) = 1.0_rp
      eye(3, 3) = 1.0_rp
   end function identity_real

end program test_gbt_frame_contract
