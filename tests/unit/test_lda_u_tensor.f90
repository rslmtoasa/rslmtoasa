program test_lda_u_tensor
   use precision_mod, only: rp
   use math_mod, only: Coulomb_mat_spherical, hcpx
   implicit none

   real(rp), parameter :: tol = 2.0e-12_rp
   logical :: failed
   character(len=32) :: mode

   call get_command_argument(1, mode)
   if (trim(mode) == '--dump') then
      call dump_reference_grid()
      stop
   end if
   failed = .false.
   call test_monopole()
   call test_independent_anisotropic_anchors()
   call test_tensor_symmetries(1)
   call test_tensor_symmetries(2)
   call test_hcpx_spd_mapping()
   call test_simultaneous_u5_covariance()
   call test_d_shell_rotation()
   call test_liechtenstein_directional_derivative()
   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine dump_reference_grid()
      integer :: l, a, b, c, d, channel
      complex(rp) :: value
      real(rp) :: f0, f2, f4

      do channel = 0, 4, 2
         f0 = 0.0_rp
         f2 = 0.0_rp
         f4 = 0.0_rp
         select case (channel)
         case (0)
            f0 = 1.0_rp
         case (2)
            f2 = 0.37_rp
         case (4)
            f4 = 0.19_rp
         end select
         do l = 1, 2
            do a = -l, l
               do b = -l, l
                  do c = -l, l
                     do d = -l, l
                        value = Coulomb_mat_spherical(l, a, b, c, d, f0, f2, f4, 0.0_rp)
                        write (*, '(a,6(1x,i3),2(1x,es25.16))') 'TENSOR', channel, l, a, b, c, d, real(value, rp), aimag(value)
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine dump_reference_grid

   subroutine test_monopole()
      integer :: a, b, c, d
      real(rp) :: expected
      complex(rp) :: actual

      do a = -2, 2
         do b = -2, 2
            do c = -2, 2
               do d = -2, 2
                  expected = 0.0_rp
                  if (a == b .and. c == d) expected = 1.37_rp
                  actual = Coulomb_mat_spherical(2, a, b, c, d, 1.37_rp, 0.0_rp, 0.0_rp, 0.0_rp)
                  call require(abs(actual - cmplx(expected, 0.0_rp, rp)) < tol, 'F0 tensor is delta_ab delta_cd')
               end do
            end do
         end do
      end do
   end subroutine test_monopole

   subroutine test_independent_anisotropic_anchors()
      complex(rp) :: p_element, d_element

      ! Values generated independently with SymPy Wigner-3j reference in
      ! tests/validation/ldau_tensor_reference.py.
      p_element = Coulomb_mat_spherical(1, -1, 0, 0, -1, 1.0_rp, 0.37_rp, 0.0_rp, 0.0_rp)
      d_element = Coulomb_mat_spherical(2, -2, -1, -1, -2, 1.0_rp, 0.37_rp, 0.19_rp, 0.0_rp)
      call require(abs(p_element - cmplx(0.0444_rp, 0.0_rp, rp)) < tol, 'independent p-shell F2 tensor anchor')
      call require(abs(d_element - cmplx(0.047460317460317474_rp, 0.0_rp, rp)) < tol, &
                    'independent d-shell F2/F4 off-diagonal tensor anchor')
   end subroutine test_independent_anisotropic_anchors

   subroutine test_tensor_symmetries(l)
      integer, intent(in) :: l
      integer :: a, b, c, d
      complex(rp) :: abcd, cdab, ba_dc

      do a = -l, l
         do b = -l, l
            do c = -l, l
               do d = -l, l
                  abcd = Coulomb_mat_spherical(l, a, b, c, d, 0.81_rp, 0.37_rp, 0.19_rp, 0.11_rp)
                  cdab = Coulomb_mat_spherical(l, c, d, a, b, 0.81_rp, 0.37_rp, 0.19_rp, 0.11_rp)
                  ba_dc = Coulomb_mat_spherical(l, b, a, d, c, 0.81_rp, 0.37_rp, 0.19_rp, 0.11_rp)
                  call require(abs(abcd - cdab) < tol, 'pair-exchange symmetry')
                  call require(abs(conjg(abcd) - ba_dc) < tol, 'Hermitian pair-density symmetry')
                  if (abs(abcd) > tol) call require(a + c == b + d, 'complex spherical tensor magnetic selection rule')
               end do
            end do
         end do
      end do
   end subroutine test_tensor_symmetries

   subroutine test_hcpx_spd_mapping()
      complex(rp) :: transform(9, 9), identity(9, 9), h_cart(9, 9), h_expected(9, 9), h_actual(9, 9)
      integer :: i, j
      real(rp) :: q

      q = 1.0_rp/sqrt(2.0_rp)
      transform = (0.0_rp, 0.0_rp)
      transform(1, 1) = (1.0_rp, 0.0_rp)
      transform(2, 4) = cmplx(-q, 0.0_rp, rp)
      transform(2, 2) = cmplx(q, 0.0_rp, rp)
      transform(3, 4) = cmplx(0.0_rp, q, rp)
      transform(3, 2) = cmplx(0.0_rp, q, rp)
      transform(4, 3) = (1.0_rp, 0.0_rp)
      transform(5, 5) = cmplx(0.0_rp, q, rp)
      transform(5, 9) = cmplx(0.0_rp, -q, rp)
      transform(6, 6) = cmplx(0.0_rp, q, rp)
      transform(6, 8) = cmplx(0.0_rp, q, rp)
      transform(7, 6) = cmplx(q, 0.0_rp, rp)
      transform(7, 8) = cmplx(-q, 0.0_rp, rp)
      transform(8, 5) = cmplx(q, 0.0_rp, rp)
      transform(8, 9) = cmplx(q, 0.0_rp, rp)
      transform(9, 7) = (1.0_rp, 0.0_rp)

      identity = (0.0_rp, 0.0_rp)
      do i = 1, 9
         identity(i, i) = (1.0_rp, 0.0_rp)
      end do
      call require(maxval(abs(matmul(conjg(transpose(transform)), transform) - identity)) < tol, &
                    'hcpx spd orbital transformation is unitary')

      h_cart = (0.0_rp, 0.0_rp)
      do i = 1, 9
         h_cart(i, i) = cmplx(0.11_rp*real(i, rp), 0.0_rp, rp)
         do j = i + 1, 9
            h_cart(i, j) = cmplx(0.013_rp*real(i + 2*j, rp), 0.007_rp*real(j - i, rp), rp)
            h_cart(j, i) = conjg(h_cart(i, j))
         end do
      end do
      h_expected = matmul(conjg(transpose(transform)), matmul(h_cart, transform))
      h_actual = h_cart
      call hcpx(h_actual, 'cart2sph')
      call require(maxval(abs(h_actual - h_expected)) < tol, 'hcpx cart2sph matches explicit spd orbital map')
      call hcpx(h_actual, 'sph2cart')
      call require(maxval(abs(h_actual - h_cart)) < tol, 'hcpx spd Hamiltonian round trip')
   end subroutine test_hcpx_spd_mapping

   subroutine test_simultaneous_u5_covariance()
      complex(rp) :: t0(-2:2, -2:2, -2:2, -2:2), t1(-2:2, -2:2, -2:2, -2:2)
      complex(rp) :: n0(-2:2, -2:2), n1(-2:2, -2:2), c(-2:2, -2:2)
      complex(rp) :: total
      complex(rp) :: zero(-2:2, -2:2), v0(-2:2, -2:2), v1(-2:2, -2:2), vexpected(-2:2, -2:2)
      complex(rp) :: vdn0(-2:2, -2:2), vdn1(-2:2, -2:2)
      real(rp) :: e0, e1
      integer :: a, b, cidx, d, i, j, k, m

      do a = -2, 2
         do b = -2, 2
            do cidx = -2, 2
               do d = -2, 2
                  t0(a, b, cidx, d) = Coulomb_mat_spherical(2, a, b, cidx, d, 1.0_rp, 0.37_rp, 0.19_rp, 0.0_rp)
               end do
            end do
         end do
      end do
      do a = -2, 2
         do b = -2, 2
            n0(a, b) = cmplx(0.013_rp*real(1 + a*b, rp), 0.017_rp*real(a - b, rp), rp)
         end do
      end do
      ! Complex DFT is a deterministic, non-real unitary matrix.
      do a = -2, 2
         do b = -2, 2
            c(a, b) = exp(cmplx(0.0_rp, 2.0_rp*acos(-1.0_rp)*real((a + 2)*(b + 2), rp)/5.0_rp, rp))/sqrt(5.0_rp)
         end do
      end do
      n1 = (0.0_rp, 0.0_rp)
      do a = -2, 2
         do b = -2, 2
            do i = -2, 2
               do j = -2, 2
                  n1(a, b) = n1(a, b) + c(i, a)*n0(i, j)*conjg(c(j, b))
               end do
            end do
         end do
      end do
      t1 = (0.0_rp, 0.0_rp)
      do a = -2, 2
         do b = -2, 2
            do cidx = -2, 2
               do d = -2, 2
                  total = (0.0_rp, 0.0_rp)
                  do i = -2, 2
                     do j = -2, 2
                        do k = -2, 2
                           do m = -2, 2
                              total = total + conjg(c(i, a))*c(j, b)*conjg(c(k, cidx))*c(m, d)*t0(i, j, k, m)
                           end do
                        end do
                     end do
                  end do
                  t1(a, b, cidx, d) = total
               end do
            end do
         end do
      end do
      e0 = interaction_energy(t0, n0)
      e1 = interaction_energy(t1, n1)
      call require(abs(e0 - e1) < 5.0e-12_rp, 'U(5) simultaneous tensor and stored-density transform')
      zero = (0.0_rp, 0.0_rp)
      call liechtenstein_potential(t0, n0, zero, 0.0_rp, 0.0_rp, v0, vdn0)
      call liechtenstein_potential(t1, n1, zero, 0.0_rp, 0.0_rp, v1, vdn1)
      vexpected = (0.0_rp, 0.0_rp)
      do a = -2, 2
         do b = -2, 2
            do i = -2, 2
               do j = -2, 2
                  vexpected(a, b) = vexpected(a, b) + conjg(c(i, a))*v0(i, j)*c(j, b)
               end do
            end do
         end do
      end do
      call require(maxval(abs(v1 - vexpected)) < 5.0e-12_rp, 'U(5) Hubbard potential covariance')
   end subroutine test_simultaneous_u5_covariance

   subroutine test_d_shell_rotation()
      complex(rp) :: tensor(-2:2, -2:2, -2:2, -2:2), density(-2:2, -2:2)
      complex(rp) :: rotated(-2:2, -2:2), rotation(-2:2, -2:2), ly(-2:2, -2:2)
      complex(rp) :: term
      real(rp) :: e0, e1, theta
      integer :: a, b, i, j, m

      do a = -2, 2
         do b = -2, 2
            do i = -2, 2
               do j = -2, 2
                  tensor(a, b, i, j) = Coulomb_mat_spherical(2, a, b, i, j, 0.73_rp, 0.29_rp, 0.14_rp, 0.0_rp)
               end do
            end do
            density(a, b) = cmplx(0.021_rp*real(1 + a*b, rp), 0.012_rp*real(a - b, rp), rp)
         end do
      end do

      ! Construct Ly=(L+ - L-)/(2i) for l=2, then exp(-i theta Ly).
      ly = (0.0_rp, 0.0_rp)
      do m = -2, 1
         term = cmplx(sqrt(real(6 - m*(m + 1), rp)), 0.0_rp, rp)
         ly(m + 1, m) = term/cmplx(0.0_rp, 2.0_rp, rp)
         ly(m, m + 1) = -term/cmplx(0.0_rp, 2.0_rp, rp)
      end do
      theta = 0.413_rp
      rotation = (0.0_rp, 0.0_rp)
      do a = -2, 2
         rotation(a, a) = (1.0_rp, 0.0_rp)
      end do
      ! Taylor series converges rapidly for this bounded 5x5 generator.
      block
         complex(rp) :: power(-2:2, -2:2), next_power(-2:2, -2:2)
         integer :: order
         power = rotation
         do order = 1, 48
            call multiply_matrix(power, ly, next_power)
            power = next_power*cmplx(0.0_rp, -theta/real(order, rp), rp)
            rotation = rotation + power
         end do
      end block
      rotated = (0.0_rp, 0.0_rp)
      do a = -2, 2
         do b = -2, 2
            do i = -2, 2
               do j = -2, 2
                  rotated(a, b) = rotated(a, b) + rotation(i, a)*density(i, j)*conjg(rotation(j, b))
               end do
            end do
         end do
      end do
      e0 = interaction_energy(tensor, density)
      e1 = interaction_energy(tensor, rotated)
      call require(abs(e0 - e1) < 5.0e-12_rp, 'rotationally invariant d-shell tensor under physical y rotation')
   end subroutine test_d_shell_rotation

   subroutine test_liechtenstein_directional_derivative()
      complex(rp) :: tensor(-2:2, -2:2, -2:2, -2:2)
      complex(rp) :: nup(-2:2, -2:2), ndn(-2:2, -2:2)
      complex(rp) :: dup(-2:2, -2:2), ddn(-2:2, -2:2)
      complex(rp) :: vup(-2:2, -2:2), vdn(-2:2, -2:2)
      real(rp), parameter :: steps(3) = [1.0e-4_rp, 1.0e-5_rp, 1.0e-6_rp]
      real(rp) :: h, eplus, eminus, numeric, analytic, error
      integer :: a, b, c, d, istep

      do a = -2, 2
         do b = -2, 2
            do c = -2, 2
               do d = -2, 2
                  tensor(a, b, c, d) = Coulomb_mat_spherical(2, a, b, c, d, 1.0_rp, 0.37_rp, 0.19_rp, 0.0_rp)
               end do
            end do
            if (a == b) then
               nup(a, b) = cmplx(0.34_rp + 0.025_rp*real(a, rp), 0.0_rp, rp)
               ndn(a, b) = cmplx(0.27_rp - 0.018_rp*real(a, rp), 0.0_rp, rp)
               dup(a, b) = cmplx(0.013_rp*real(a + 3, rp), 0.0_rp, rp)
               ddn(a, b) = cmplx(0.009_rp*real(2 - a, rp), 0.0_rp, rp)
            else
               nup(a, b) = cmplx(0.004_rp*real(1 + a*b, rp), 0.003_rp*real(a - b, rp), rp)
               ndn(a, b) = cmplx(0.003_rp*real(2 + a*b, rp), -0.002_rp*real(a - b, rp), rp)
               dup(a, b) = cmplx(0.002_rp*real(a + b, rp), 0.004_rp*real(a - b, rp), rp)
               ddn(a, b) = cmplx(0.003_rp*real(a + b, rp), 0.002_rp*real(a - b, rp), rp)
            end if
         end do
      end do
      call liechtenstein_potential(tensor, nup, ndn, 1.0_rp, (0.37_rp + 0.19_rp)/14.0_rp, vup, vdn)
      analytic = 0.0_rp
      do a = -2, 2
         do b = -2, 2
            analytic = analytic + real(vup(a, b)*dup(a, b) + vdn(a, b)*ddn(a, b), rp)
         end do
      end do
      call require(maxval(abs(vup - conjg(transpose(vup)))) < 5.0e-12_rp, 'complex d-shell up potential is Hermitian')
      call require(maxval(abs(vdn - conjg(transpose(vdn)))) < 5.0e-12_rp, 'complex d-shell down potential is Hermitian')
      do istep = 1, size(steps)
         h = steps(istep)
         eplus = liechtenstein_energy(tensor, nup + h*dup, ndn + h*ddn, 1.0_rp, (0.37_rp + 0.19_rp)/14.0_rp)
         eminus = liechtenstein_energy(tensor, nup - h*dup, ndn - h*ddn, 1.0_rp, (0.37_rp + 0.19_rp)/14.0_rp)
         numeric = (eplus - eminus)/(2.0_rp*h)
         error = abs(numeric - analytic)
         call require(error < 2.0e-8_rp, 'complex off-diagonal Liechtenstein directional derivative')
      end do
   end subroutine test_liechtenstein_directional_derivative

   subroutine liechtenstein_potential(tensor, nup, ndn, u, jbar, vup, vdn)
      complex(rp), intent(in) :: tensor(-2:2, -2:2, -2:2, -2:2), nup(-2:2, -2:2), ndn(-2:2, -2:2)
      real(rp), intent(in) :: u, jbar
      complex(rp), intent(out) :: vup(-2:2, -2:2), vdn(-2:2, -2:2)
      integer :: a, b, c, d
      real(rp) :: nup_trace, ndn_trace

      nup_trace = 0.0_rp
      ndn_trace = 0.0_rp
      do a = -2, 2
         nup_trace = nup_trace + real(nup(a, a), rp)
         ndn_trace = ndn_trace + real(ndn(a, a), rp)
      end do
      vup = (0.0_rp, 0.0_rp)
      vdn = (0.0_rp, 0.0_rp)
      do a = -2, 2
         do b = -2, 2
            do c = -2, 2
               do d = -2, 2
                  vup(a, b) = vup(a, b) + tensor(a, b, c, d)*ndn(c, d) + &
                              (tensor(a, b, c, d) - tensor(a, d, c, b))*nup(c, d)
                  vdn(a, b) = vdn(a, b) + tensor(a, b, c, d)*nup(c, d) + &
                              (tensor(a, b, c, d) - tensor(a, d, c, b))*ndn(c, d)
               end do
            end do
         end do
      end do
      do a = -2, 2
         vup(a, a) = vup(a, a) - u*(nup_trace + ndn_trace - 0.5_rp) + jbar*(nup_trace - 0.5_rp)
         vdn(a, a) = vdn(a, a) - u*(nup_trace + ndn_trace - 0.5_rp) + jbar*(ndn_trace - 0.5_rp)
      end do
   end subroutine liechtenstein_potential

   function liechtenstein_energy(tensor, nup, ndn, u, jbar) result(energy)
      complex(rp), intent(in) :: tensor(-2:2, -2:2, -2:2, -2:2), nup(-2:2, -2:2), ndn(-2:2, -2:2)
      real(rp), intent(in) :: u, jbar
      real(rp) :: energy, direct, exchange, nup_trace, ndn_trace
      complex(rp) :: densities(-2:2, -2:2, 2)
      integer :: a, b, c, d, s, t

      densities(:, :, 1) = nup
      densities(:, :, 2) = ndn
      direct = 0.0_rp
      exchange = 0.0_rp
      do s = 1, 2
         do t = 1, 2
            do a = -2, 2
               do b = -2, 2
                  do c = -2, 2
                     do d = -2, 2
                        direct = direct + real(tensor(a, b, c, d)*densities(a, b, s)*densities(c, d, t), rp)
                     end do
                  end do
               end do
            end do
         end do
         do a = -2, 2
            do b = -2, 2
               do c = -2, 2
                  do d = -2, 2
                     exchange = exchange + real(tensor(a, b, c, d)*densities(a, d, s)*densities(c, b, s), rp)
                  end do
               end do
            end do
         end do
      end do
      nup_trace = sum([(real(nup(a, a), rp), a=-2,2)])
      ndn_trace = sum([(real(ndn(a, a), rp), a=-2,2)])
      energy = 0.5_rp*(direct - exchange) - 0.5_rp*u*(nup_trace + ndn_trace)*(nup_trace + ndn_trace - 1.0_rp) &
               + 0.5_rp*jbar*(nup_trace*(nup_trace - 1.0_rp) + ndn_trace*(ndn_trace - 1.0_rp))
   end function liechtenstein_energy

   function interaction_energy(tensor, density) result(energy)
      complex(rp), intent(in) :: tensor(-2:2, -2:2, -2:2, -2:2), density(-2:2, -2:2)
      real(rp) :: energy, direct, exchange
      integer :: a, b, c, d

      direct = 0.0_rp
      exchange = 0.0_rp
      do a = -2, 2
         do b = -2, 2
            do c = -2, 2
               do d = -2, 2
                  direct = direct + real(tensor(a, b, c, d)*density(a, b)*density(c, d), rp)
                  exchange = exchange + real(tensor(a, b, c, d)*density(a, d)*density(c, b), rp)
               end do
            end do
         end do
      end do
      energy = 0.5_rp*(direct - exchange)
   end function interaction_energy

   subroutine multiply_matrix(left, right, output)
      complex(rp), intent(in) :: left(-2:2, -2:2), right(-2:2, -2:2)
      complex(rp), intent(out) :: output(-2:2, -2:2)
      integer :: a, b, k

      output = (0.0_rp, 0.0_rp)
      do a = -2, 2
         do b = -2, 2
            do k = -2, 2
               output(a, b) = output(a, b) + left(a, k)*right(k, b)
            end do
         end do
      end do
   end subroutine multiply_matrix

   subroutine require(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label

      if (.not. condition) then
         write (*, '(a,a)') 'FAIL: ', trim(label)
         failed = .true.
      end if
   end subroutine require

end program test_lda_u_tensor
