program test_gbt_wp02_gauge_shifted_k
   ! WP02 is a fixed-potential operator test.  The GBT side uses the
   ! production primitive-link/lift/contract path.  The shifted-k reference
   ! is assembled independently from ordinary spin-resolved scalar blocks.
   use precision_mod, only: rp
   use gbt_structure_mod, only: gbt_bond_phase, gbt_endpoint_link, &
      gbt_lift_orbital_block, gbt_contract_collinear
   implicit none

   integer, parameter :: norb = 2
   integer, parameter :: nb = 2*norb
   integer, parameter :: nbond_one = 5
   integer, parameter :: nloop = 3
   real(rp), parameter :: tol = 1.0e-12_rp
   real(rp), parameter :: negative_tol = 1.0e-3_rp
   real(rp), parameter :: alat = 3.4_rp
   real(rp), parameter :: pi = acos(-1.0_rp)
   real(rp), parameter :: q_test(3) = [0.173_rp, -0.117_rp, 0.086_rp]
   real(rp), parameter :: k_test(3) = [0.137_rp, -0.219_rp, 0.083_rp]
   real(rp), parameter :: q_points(3,3) = reshape([ &
      0.173_rp, -0.117_rp, 0.086_rp, &
     -0.271_rp,  0.314_rp, 0.119_rp, &
      0.223_rp,  0.071_rp, -0.307_rp], [3,3])
   real(rp), parameter :: k_points(3,3) = reshape([ &
      0.137_rp, -0.219_rp, 0.083_rp, &
     -0.241_rp,  0.193_rp, 0.157_rp, &
      0.311_rp,  0.047_rp, -0.283_rp], [3,3])
   real(rp), parameter :: one_translation(3,nbond_one) = reshape([ &
      0.0_rp, 0.0_rp, 0.0_rp, &
      1.0_rp, 0.0_rp, 0.0_rp, &
     -1.0_rp, 0.0_rp, 0.0_rp, &
      0.0_rp, 1.0_rp, 0.0_rp, &
      0.0_rp,-1.0_rp, 0.0_rp], [3,nbond_one])

   complex(rp) :: raw_one(norb,norb,nbond_one)
   complex(rp) :: raw_inversion(norb,norb,nbond_one)
   complex(rp) :: w0(norb), w1(norb), c0(norb), c1(norb)
   logical :: failed
   real(rp) :: largest_abs, largest_rel, largest_eigen

   failed = .false.
   largest_abs = 0.0_rp
   largest_rel = 0.0_rp
   largest_eigen = 0.0_rp
   call initialize_fixture()
   call test_endpoint_link_algebra()
   call test_closed_loop_gauge()
   call test_realspace_bond_reversal()
   call test_reciprocal_hermiticity()
   call test_shifted_k_identity()
   call test_q_reversal_symmetry()
   call test_negative_controls()

   write (*, '(a,3(es12.4,1x))') 'WP02 overall maxima (abs, rel, eigen): ', &
      largest_abs, largest_rel, largest_eigen
   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine initialize_fixture()
      ! Complex directed blocks expose conjugation and endpoint ordering.
      raw_one(:,:,1) = reshape([ &
         cmplx(0.41_rp, 0.00_rp, rp), cmplx(0.07_rp, 0.03_rp, rp), &
         cmplx(0.07_rp,-0.03_rp, rp), cmplx(0.53_rp, 0.00_rp, rp)], [norb,norb])
      raw_one(:,:,2) = reshape([ &
         cmplx(0.70_rp, 0.20_rp, rp), cmplx(-0.30_rp, 0.40_rp, rp), &
         cmplx(0.10_rp,-0.50_rp, rp), cmplx(0.80_rp,-0.10_rp, rp)], [norb,norb])
      raw_one(:,:,3) = transpose(conjg(raw_one(:,:,2)))
      raw_one(:,:,4) = reshape([ &
         cmplx(-0.22_rp, 0.17_rp, rp), cmplx(0.15_rp,-0.26_rp, rp), &
         cmplx(-0.31_rp, 0.09_rp, rp), cmplx(0.44_rp, 0.05_rp, rp)], [norb,norb])
      raw_one(:,:,5) = transpose(conjg(raw_one(:,:,4)))

      ! Inversion-symmetric control: the two opposite translations have the
      ! same real symmetric block.  This is used only for q -> -q symmetry.
      raw_inversion(:,:,1) = raw_one(:,:,1)
      raw_inversion(:,:,2) = reshape([ &
         cmplx(0.62_rp, 0.0_rp, rp), cmplx(-0.18_rp, 0.0_rp, rp), &
         cmplx(-0.18_rp, 0.0_rp, rp), cmplx(0.37_rp, 0.0_rp, rp)], [norb,norb])
      raw_inversion(:,:,3) = raw_inversion(:,:,2)
      raw_inversion(:,:,4) = reshape([ &
         cmplx(-0.29_rp, 0.0_rp, rp), cmplx(0.11_rp, 0.0_rp, rp), &
         cmplx(0.11_rp, 0.0_rp, rp), cmplx(0.51_rp, 0.0_rp, rp)], [norb,norb])
      raw_inversion(:,:,5) = raw_inversion(:,:,4)

      w0 = [cmplx(1.08_rp, 0.0_rp, rp), cmplx(0.74_rp, 0.0_rp, rp)]
      w1 = [cmplx(0.23_rp, 0.0_rp, rp), cmplx(-0.17_rp, 0.0_rp, rp)]
      c0 = [cmplx(-0.19_rp, 0.0_rp, rp), cmplx(0.27_rp, 0.0_rp, rp)]
      c1 = [cmplx(0.39_rp, 0.0_rp, rp), cmplx(-0.13_rp, 0.0_rp, rp)]
   end subroutine initialize_fixture

   subroutine test_endpoint_link_algebra()
      call check_link_case('one-sublattice endpoint frames', q_test, [1.2_rp,-0.7_rp,0.8_rp], &
         0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp)
      call check_link_case('multi-sublattice endpoint frames', q_test, [-0.4_rp,1.1_rp,0.6_rp], &
         0.63_rp, 0.27_rp, 1.11_rp,-0.58_rp)
      call check_link_case('second multi-sublattice endpoint frames', q_points(:,2), [0.2_rp,0.9_rp,-1.4_rp], &
         1.27_rp,-0.31_rp, 0.44_rp, 0.82_rp)
   end subroutine test_endpoint_link_algebra

   subroutine check_link_case(label, q, d, theta_a, phi_a, theta_b, phi_b)
      character(len=*), intent(in) :: label
      real(rp), intent(in) :: q(3), d(3), theta_a, phi_a, theta_b, phi_b
      complex(rp) :: link(2,2), reverse(2,2), ua(2,2), ub(2,2), rz(2,2), oracle(2,2)
      real(rp) :: unitary_error, reversal_error, oracle_error, alpha

      call gbt_endpoint_link(q, d, alat, theta_a, phi_a, theta_b, phi_b, link)
      call gbt_endpoint_link(q, -d, alat, theta_b, phi_b, theta_a, phi_a, reverse)
      call independent_frame(theta_a, phi_a, ua)
      call independent_frame(theta_b, phi_b, ub)
      alpha = gbt_bond_phase(q, d, alat)
      rz = cmplx(0.0_rp, 0.0_rp, rp)
      rz(1,1) = exp(cmplx(0.0_rp,-0.5_rp*alpha,rp))
      rz(2,2) = exp(cmplx(0.0_rp, 0.5_rp*alpha,rp))
      oracle = matmul(transpose(conjg(ua)), matmul(rz, ub))
      unitary_error = maxval(abs(matmul(transpose(conjg(link)), link) - identity2()))
      reversal_error = maxval(abs(reverse - transpose(conjg(link))))
      oracle_error = maxval(abs(link - oracle))
      call report('link '//trim(label)//' unitarity', unitary_error, unitary_error, 0.0_rp)
      call report('link '//trim(label)//' reversal', reversal_error, reversal_error, 0.0_rp)
      call report('link '//trim(label)//' independent oracle', oracle_error, oracle_error, 0.0_rp)
   end subroutine check_link_case

   subroutine test_closed_loop_gauge()
      real(rp) :: d1(3), d2(3), d3(3)
      real(rp) :: loop_error
      complex(rp) :: l12(2,2), l23(2,2), l31(2,2), product(2,2)

      d1 = [1.0_rp, 0.0_rp, 0.0_rp]*alat
      d2 = [0.0_rp, 1.0_rp, 0.0_rp]*alat
      d3 = [-1.0_rp,-1.0_rp,0.0_rp]*alat
      call gbt_endpoint_link(q_test, d1, alat, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, l12)
      call gbt_endpoint_link(q_test, d2, alat, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, l23)
      call gbt_endpoint_link(q_test, d3, alat, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, l31)
      product = matmul(l12, matmul(l23, l31))
      loop_error = maxval(abs(product - identity2()))
      call report('one-sublattice translated closed loop', loop_error, loop_error, 0.0_rp)

      d1 = [1.0_rp, 0.5_rp, 0.0_rp]*alat
      d2 = [-0.25_rp, 0.5_rp, 0.75_rp]*alat
      d3 = -d1-d2
      call gbt_endpoint_link(q_points(:,2), d1, alat, 0.63_rp, 0.27_rp, 1.11_rp,-0.58_rp, l12)
      call gbt_endpoint_link(q_points(:,2), d2, alat, 1.11_rp,-0.58_rp, 0.63_rp, 0.27_rp, l23)
      call gbt_endpoint_link(q_points(:,2), d3, alat, 0.63_rp, 0.27_rp, 0.63_rp, 0.27_rp, l31)
      product = matmul(l12, matmul(l23, l31))
      loop_error = maxval(abs(product - identity2()))
      call report('multi-sublattice translated closed loop', loop_error, loop_error, 0.0_rp)
   end subroutine test_closed_loop_gauge

   subroutine test_realspace_bond_reversal()
      complex(rp) :: forward(nb,nb), reverse(nb,nb)
      real(rp) :: error

      call contracted_bond(raw_one(:,:,2), q_test, one_translation(:,2)*alat, &
         0.63_rp, 0.27_rp, 1.11_rp,-0.58_rp, w0, w1, w0, w1, forward)
      call contracted_bond(transpose(conjg(raw_one(:,:,2))), q_test, -one_translation(:,2)*alat, &
         1.11_rp,-0.58_rp, 0.63_rp, 0.27_rp, w0, w1, w0, w1, reverse)
      error = maxval(abs(reverse - transpose(conjg(forward))))
      call report('contracted real-space bond reversal', error, error, 0.0_rp)
   end subroutine test_realspace_bond_reversal

   subroutine test_reciprocal_hermiticity()
      complex(rp) :: h(nb,nb)
      real(rp) :: error
      integer :: ip

      do ip = 1, size(k_points,2)
         call build_gbt_one(raw_one, q_points(:,ip), k_points(:,ip), h)
         error = maxval(abs(h - transpose(conjg(h))))
         call report('reciprocal H(k) Hermiticity '//trim(index_string(ip)), error, error, 0.0_rp)
      end do
   end subroutine test_reciprocal_hermiticity

   subroutine test_shifted_k_identity()
      complex(rp) :: h_gbt(nb,nb), h_expected(nb,nb), h_up(norb,norb), h_down(norb,norb)
      real(rp) :: abs_error, rel_error, eigen_error
      integer :: ip

      do ip = 1, size(k_points,2)
         call build_gbt_one(raw_one, q_points(:,ip), k_points(:,ip), h_gbt)
         call build_ordinary_sector(raw_one, k_points(:,ip)-0.5_rp*q_points(:,ip), 1, h_up)
         call build_ordinary_sector(raw_one, k_points(:,ip)+0.5_rp*q_points(:,ip), 2, h_down)
         h_expected = cmplx(0.0_rp, 0.0_rp, rp)
         h_expected(1:norb,1:norb) = h_up
         h_expected(norb+1:nb,norb+1:nb) = h_down
         abs_error = maxval(abs(h_gbt-h_expected))
         rel_error = norm_error(h_gbt, h_expected)
         eigen_error = eigenvalue_error(h_gbt, h_expected)
         call report('shifted-k matrix identity '//trim(index_string(ip)), abs_error, rel_error, eigen_error)
      end do
   end subroutine test_shifted_k_identity

   subroutine test_q_reversal_symmetry()
      complex(rp) :: h_plus(nb,nb), h_minus(nb,nb)
      real(rp) :: error

      ! For a centrosymmetric SOC-free fixture the useful pointwise statement
      ! is spec H(k,q) = spec H(-k,-q), not generally H(k,q)=H(k,-q).
      call build_gbt_one(raw_inversion, q_test, k_test, h_plus)
      call build_gbt_one(raw_inversion, -q_test, -k_test, h_minus)
      error = eigenvalue_error(h_plus, h_minus)
      call report('centrosymmetric q <-> -q spectral symmetry', 0.0_rp, 0.0_rp, error)
   end subroutine test_q_reversal_symmetry

   subroutine test_negative_controls()
      complex(rp) :: correct(nb,nb), wrong_half(nb,nb), wrong_sign(nb,nb)
      complex(rp) :: correct_bond(nb,nb), wrong_endpoint(nb,nb)
      real(rp) :: half_error, sign_error, endpoint_error

      call build_gbt_one(raw_one, q_test, k_test, correct)
      call build_gbt_one(raw_one, 2.0_rp*q_test, k_test, wrong_half)
      call build_gbt_one(raw_one, -q_test, k_test, wrong_sign)
      half_error = maxval(abs(wrong_half-correct))
      sign_error = maxval(abs(wrong_sign-correct))
      call report_negative('negative control: missing half phase', half_error)
      call report_negative('negative control: q sign flip', sign_error)

      call contracted_bond(raw_one(:,:,2), q_test, one_translation(:,2)*alat, &
         0.63_rp, 0.27_rp, 1.11_rp,-0.58_rp, w0, w1, w0, w1, correct_bond)
      call contracted_bond(raw_one(:,:,2), q_test, one_translation(:,2)*alat, &
         1.11_rp,-0.58_rp, 0.63_rp, 0.27_rp, w0, w1, w0, w1, wrong_endpoint)
      endpoint_error = maxval(abs(wrong_endpoint-correct_bond))
      call report_negative('negative control: reversed endpoint order', endpoint_error)
   end subroutine test_negative_controls

   subroutine build_gbt_one(raw, q, k, h)
      complex(rp), intent(in) :: raw(norb,norb,nbond_one)
      real(rp), intent(in) :: q(3), k(3)
      complex(rp), intent(out) :: h(nb,nb)
      complex(rp) :: link(2,2), lifted(nb,nb), pair(nb,nb), phase
      integer :: bond, i

      h = cmplx(0.0_rp, 0.0_rp, rp)
      do bond = 1, nbond_one
         call gbt_endpoint_link(q, one_translation(:,bond)*alat, alat, 0.0_rp, 0.0_rp, &
            0.0_rp, 0.0_rp, link)
         call gbt_lift_orbital_block(raw(:,:,bond), link, lifted)
         call gbt_contract_collinear(lifted, w0, w1, w0, w1, 1.0_rp, 1.0_rp, pair)
         if (bond == 1) then
            do i = 1, norb
               pair(i,i) = pair(i,i) + c0(i) + c1(i)
               pair(i+norb,i+norb) = pair(i+norb,i+norb) + c0(i) - c1(i)
            end do
         end if
         phase = exp(cmplx(0.0_rp, 2.0_rp*pi*dot_product(k, one_translation(:,bond)), rp))
         h = h + pair*phase
      end do
   end subroutine build_gbt_one

   subroutine build_ordinary_sector(raw, k, sector, h)
      complex(rp), intent(in) :: raw(norb,norb,nbond_one)
      real(rp), intent(in) :: k(3)
      integer, intent(in) :: sector
      complex(rp), intent(out) :: h(norb,norb)
      complex(rp) :: phase, source_factor(norb), target_factor(norb)
      integer :: bond, i, j

      if (sector == 1) then
         source_factor = w0+w1
         target_factor = w0+w1
      else if (sector == 2) then
         source_factor = w0-w1
         target_factor = w0-w1
      else
         failed = .true.
         h = cmplx(0.0_rp,0.0_rp,rp)
         return
      end if
      h = cmplx(0.0_rp, 0.0_rp, rp)
      do bond = 1, nbond_one
         phase = exp(cmplx(0.0_rp, 2.0_rp*pi*dot_product(k, one_translation(:,bond)), rp))
         do i = 1, norb
            do j = 1, norb
               h(i,j) = h(i,j) + source_factor(i)*raw(i,j,bond)*target_factor(j)*phase
            end do
         end do
         if (bond == 1) then
            if (sector == 1) then
               do i = 1, norb
                  h(i,i) = h(i,i) + c0(i)+c1(i)
               end do
            else
               do i = 1, norb
                  h(i,i) = h(i,i) + c0(i)-c1(i)
               end do
            end if
         end if
      end do
   end subroutine build_ordinary_sector

   subroutine contracted_bond(raw, q, d, theta_a, phi_a, theta_b, phi_b, w0_a, w1_a, w0_b, w1_b, h)
      complex(rp), intent(in) :: raw(norb,norb), w0_a(norb), w1_a(norb), w0_b(norb), w1_b(norb)
      real(rp), intent(in) :: q(3), d(3), theta_a, phi_a, theta_b, phi_b
      complex(rp), intent(out) :: h(nb,nb)
      complex(rp) :: link(2,2), lifted(nb,nb)

      call gbt_endpoint_link(q, d, alat, theta_a, phi_a, theta_b, phi_b, link)
      call gbt_lift_orbital_block(raw, link, lifted)
      call gbt_contract_collinear(lifted, w0_a, w1_a, w0_b, w1_b, 1.0_rp, 1.0_rp, h)
   end subroutine contracted_bond

   subroutine independent_frame(theta, phi, frame)
      real(rp), intent(in) :: theta, phi
      complex(rp), intent(out) :: frame(2,2)
      complex(rp) :: eiphi
      real(rp) :: c, s

      c = cos(0.5_rp*theta)
      s = sin(0.5_rp*theta)
      eiphi = exp(cmplx(0.0_rp, phi, rp))
      frame(1,1) = cmplx(c,0.0_rp,rp)
      frame(1,2) = -conjg(eiphi)*s
      frame(2,1) = eiphi*s
      frame(2,2) = cmplx(c,0.0_rp,rp)
   end subroutine independent_frame

   function identity2() result(identity)
      complex(rp) :: identity(2,2)
      identity = cmplx(0.0_rp,0.0_rp,rp)
      identity(1,1) = cmplx(1.0_rp,0.0_rp,rp)
      identity(2,2) = identity(1,1)
   end function identity2

   function norm_error(actual, reference) result(error)
      complex(rp), intent(in) :: actual(:,:), reference(:,:)
      real(rp) :: error, denominator
      denominator = max(1.0_rp, sqrt(sum(abs(reference)**2)))
      error = sqrt(sum(abs(actual-reference)**2))/denominator
   end function norm_error

   function eigenvalue_error(first, second) result(error)
      complex(rp), intent(in) :: first(:,:), second(:,:)
      complex(rp), allocatable :: first_copy(:,:), second_copy(:,:), work(:), query(:)
      real(rp), allocatable :: first_eigen(:), second_eigen(:), rwork(:)
      real(rp) :: error
      integer :: n, info, lwork
      external zheev

      n = size(first,1)
      allocate(first_copy(n,n), second_copy(n,n), first_eigen(n), second_eigen(n), &
         rwork(max(1,3*n-2)), query(1))
      first_copy = first
      second_copy = second
      call zheev('N','U',n,first_copy,n,first_eigen,query,-1,rwork,info)
      if (info /= 0) then
         failed = .true.; error = huge(1.0_rp); return
      end if
      lwork = max(1,int(real(query(1),rp)))
      allocate(work(lwork))
      call zheev('N','U',n,first_copy,n,first_eigen,work,lwork,rwork,info)
      if (info /= 0) then
         failed = .true.; error = huge(1.0_rp); return
      end if
      call zheev('N','U',n,second_copy,n,second_eigen,work,lwork,rwork,info)
      if (info /= 0) then
         failed = .true.; error = huge(1.0_rp); return
      end if
      error = maxval(abs(first_eigen-second_eigen))
      deallocate(first_copy,second_copy,first_eigen,second_eigen,rwork,query,work)
   end function eigenvalue_error

   function index_string(index) result(label)
      integer, intent(in) :: index
      character(len=8) :: label
      write(label,'(i0)') index
   end function index_string

   subroutine report(label, abs_error, rel_error, eigen_error)
      character(len=*), intent(in) :: label
      real(rp), intent(in) :: abs_error, rel_error, eigen_error
      write (*, '(a,3(es12.4,1x))') trim(label)//' (abs, rel, eigen): ', abs_error, rel_error, eigen_error
      largest_abs = max(largest_abs, abs_error)
      largest_rel = max(largest_rel, rel_error)
      largest_eigen = max(largest_eigen, eigen_error)
      if (abs_error > tol .or. rel_error > tol .or. eigen_error > tol) then
         write (*, '(a,a)') 'FAIL ', trim(label)
         failed = .true.
      end if
   end subroutine report

   subroutine report_negative(label, error)
      character(len=*), intent(in) :: label
      real(rp), intent(in) :: error
      write (*, '(a,es12.4)') trim(label)//' residual: ', error
      if (error <= negative_tol) then
         write (*, '(a,a)') 'FAIL ', trim(label)//' was not detected'
         failed = .true.
      end if
   end subroutine report_negative

end program test_gbt_wp02_gauge_shifted_k
