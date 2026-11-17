! Unit tests of the full-spinor d-shell Liechtenstein + FLL functional (hubbard_spinor_mod).
!
! Expected values come from: the existing collinear production routines (a), exact algebraic
! identities and invariances (b, d, e, g), finite differences of the functional (c), and the
! committed output of the independent R08 oracle (f, spinor_oracle_reference.txt). Every check
! prints measured value and tolerance; a control must exceed its tolerance to pass.
program test_hubbard_spinor
   use precision_mod, only: rp
   use basis_mod, only: basis_init, nb, spin_off
   use control_mod, only: control
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use logger_mod, only: g_logger
   use hubbard_spinor_mod, only: hubbard_d_tensor, hubbard_spinor_functional
   implicit none

   integer, parameter :: nd = 5, n2 = 10
   real(rp), parameter :: pi = 3.14159265358979323846_rp
   logical :: failed = .false.
   integer(8) :: seed = 20261026_8
   character(len=512) :: reference_path

   call get_command_argument(1, reference_path)
   call g_logger%init()
   call test_collinear_reduction()
   call test_su2_invariance()
   call test_derivative()
   call test_dudarev_limit()
   call test_single_projector()
   call test_r08_oracle()
   call test_hermitian_potential()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains

   function uniform() result(x)
      real(rp) :: x
      seed = mod(seed*1103515245_8 + 12345_8, 2147483648_8)
      x = real(seed, rp)/2147483648.0_rp
   end function uniform

   !> measured <= tol passes; for a control, measured > tol passes.
   subroutine check(name, measured, tol, control)
      character(len=*), intent(in) :: name
      real(rp), intent(in) :: measured, tol
      logical, intent(in), optional :: control
      logical :: ok, is_control

      is_control = .false.
      if (present(control)) is_control = control
      if (is_control) then
         ok = measured > tol
      else
         ok = measured <= tol
      end if
      if (is_control) then
         write (*, '(a,1x,a,2(2x,es10.3))') merge('PASS', 'FAIL', ok), trim(name)//' [control: must exceed]', measured, tol
      else
         write (*, '(a,1x,a,2(2x,es10.3))') merge('PASS', 'FAIL', ok), trim(name), measured, tol
      end if
      if (.not. ok) failed = .true.
   end subroutine check

   function rand_hermitian(scale) result(h)
      real(rp), intent(in) :: scale
      complex(rp) :: h(n2, n2)
      integer :: i, j

      do i = 1, n2
         h(i, i) = cmplx(scale*(uniform() - 0.5_rp), 0.0_rp, rp)
         do j = i + 1, n2
            h(i, j) = scale*cmplx(uniform() - 0.5_rp, uniform() - 0.5_rp, rp)
            h(j, i) = conjg(h(i, j))
         end do
      end do
   end function rand_hermitian

   !> Random unitary by Gram-Schmidt of a random complex matrix.
   function rand_unitary(n) result(q)
      integer, intent(in) :: n
      complex(rp) :: q(n, n)
      integer :: i, k

      do i = 1, n
         do k = 1, n
            q(k, i) = cmplx(uniform() - 0.5_rp, uniform() - 0.5_rp, rp)
         end do
         do k = 1, i - 1
            q(:, i) = q(:, i) - dot_product(q(:, k), q(:, i))*q(:, k)
         end do
         q(:, i) = q(:, i)/sqrt(real(dot_product(q(:, i), q(:, i)), rp))
      end do
   end function rand_unitary

   !> Physical spinor density: unitary times occupations in [0,1].
   function rand_physical() result(p)
      complex(rp) :: p(n2, n2), q(n2, n2)
      integer :: i

      q = rand_unitary(n2)
      p = (0.0_rp, 0.0_rp)
      do i = 1, n2
         p = p + uniform()*spread(q(:, i), 2, n2)*spread(conjg(q(:, i)), 1, n2)
      end do
   end function rand_physical

   function max_abs(a) result(m)
      complex(rp), intent(in) :: a(:, :)
      real(rp) :: m

      m = maxval(abs(a))
   end function max_abs

   subroutine test_collinear_reduction()
      real(rp), parameter :: u = 1.4_rp, jbar = 0.22_rp
      type(hamiltonian), target :: ham
      type(lattice), target :: lat
      type(control), target :: ctl
      complex(rp) :: t(nd, nd, nd, nd), p(n2, n2), pup(nd, nd), pdn(nd, nd), v(n2, n2), vref(n2, n2)
      complex(rp) :: density(1, 4, 2, 7, 7), potential(1, 4, 2, 7, 7), pshuf(n2, n2)
      real(rp) :: e_int, e_dc, e_hub, ci(1, 4), cd(1, 4), ch(1, 4), ct(1, 4)
      real(rp) :: s_int, s_dc, s_hub

      call basis_init(2)
      call ctl%restore_to_default()
      ctl%lmax = 2
      ctl%nsp = 2
      lat%ntype = 1
      allocate (lat%symbolic_atoms(1))
      call lat%symbolic_atoms(1)%restore_to_default()
      ham%lattice => lat
      ham%control => ctl
      allocate (ham%hubbard_u_pot(nb, nb, 1))
      ham%hubbard_u_pot = 0.0_rp
      ham%hubbard_u_potential_form = 'liechtenstein'
      lat%symbolic_atoms(1)%potential%hubbard_u(3) = u
      lat%symbolic_atoms(1)%potential%hubbard_j(3) = jbar

      ! Block-diagonal P: random Hermitian spin blocks with different up and down content.
      p = (0.0_rp, 0.0_rp)
      pup = rand_physical_block(0.9_rp)
      pdn = rand_physical_block(0.5_rp)
      p(1:nd, 1:nd) = pup
      p(nd + 1:n2, nd + 1:n2) = pdn
      lat%symbolic_atoms(1)%potential%ldm(3, 1, 1:5, 1:5) = pup
      lat%symbolic_atoms(1)%potential%ldm(3, 2, 1:5, 1:5) = pdn
      call ham%calculate_hubbard_u_potential_general()
      vref = (0.0_rp, 0.0_rp)
      vref(1:nd, 1:nd) = ham%hubbard_u_pot(5:9, 5:9, 1)
      vref(nd + 1:n2, nd + 1:n2) = ham%hubbard_u_pot(5 + spin_off:9 + spin_off, 5 + spin_off:9 + spin_off, 1)
      density = (0.0_rp, 0.0_rp)
      potential = (0.0_rp, 0.0_rp)
      density(1, 3, 1, 1:5, 1:5) = pup
      density(1, 3, 2, 1:5, 1:5) = pdn
      potential(1, 3, 1, 1:5, 1:5) = vref(1:nd, 1:nd)
      potential(1, 3, 2, 1:5, 1:5) = vref(nd + 1:n2, nd + 1:n2)
      call ham%evaluate_hubbard_energy(density, potential, ci, cd, ch, ct)

      call hubbard_d_tensor(u, jbar, t)
      call hubbard_spinor_functional(t, u, jbar, p, e_int, e_dc, e_hub, v)
      call check('(a) collinear reduction: V (production potential)', max_abs(v - vref), 1.0e-13_rp)
      call check('(a) collinear reduction: E_int', abs(e_int - ci(1, 3)), 1.0e-13_rp)
      call check('(a) collinear reduction: E_DC', abs(e_dc - cd(1, 3)), 1.0e-13_rp)
      call check('(a) collinear reduction: E_Hub', abs(e_hub - ch(1, 3)), 1.0e-13_rp)

      pshuf = p
      pshuf(1:nd, 1:nd) = pdn
      pshuf(nd + 1:n2, nd + 1:n2) = pup
      call hubbard_spinor_functional(t, u, jbar, pshuf, s_int, s_dc, s_hub, v)
      call check('(a) control, spin blocks swapped: V', max_abs(v - vref), 1.0e-13_rp, .true.)
      ! swapping the two spin blocks is a symmetry of E, so the energy control permutes the
      ! orbital labels of the down block instead
      pshuf = p
      pshuf(nd + 1:n2, nd + 1:n2) = pdn([2, 3, 4, 5, 1], [2, 3, 4, 5, 1])
      call hubbard_spinor_functional(t, u, jbar, pshuf, s_int, s_dc, s_hub, v)
      call check('(a) control, down block orbitals permuted: E_Hub', abs(s_hub - ch(1, 3)), 1.0e-13_rp, .true.)
   end subroutine test_collinear_reduction

   function rand_physical_block(scale) result(b)
      real(rp), intent(in) :: scale
      complex(rp) :: b(nd, nd), q(nd, nd)
      integer :: i

      q = rand_unitary(nd)
      b = (0.0_rp, 0.0_rp)
      do i = 1, nd
         b = b + scale*uniform()*spread(q(:, i), 2, nd)*spread(conjg(q(:, i)), 1, nd)
      end do
   end function rand_physical_block

   subroutine test_su2_invariance()
      real(rp), parameter :: u = 1.4_rp, jbar = 0.22_rp
      complex(rp) :: t(nd, nd, nd, nd), p(n2, n2), pr(n2, n2), v(n2, n2), vr(n2, n2), rot(n2, n2), s(2, 2)
      real(rp) :: e_int, e_dc, e_hub, er_int, er_dc, er_hub, worst_e, worst_v, worst_c_v, theta, phi, chi
      integer :: trial, a, i, k

      call hubbard_d_tensor(u, jbar, t)
      worst_e = 0.0_rp
      worst_v = 0.0_rp
      worst_c_v = 0.0_rp
      do trial = 1, 5
         theta = pi*uniform()
         phi = 2.0_rp*pi*uniform()
         chi = 2.0_rp*pi*uniform()
         s(1, 1) = cmplx(cos(theta/2.0_rp)*cos(phi), cos(theta/2.0_rp)*sin(phi), rp)
         s(1, 2) = cmplx(sin(theta/2.0_rp)*cos(chi), sin(theta/2.0_rp)*sin(chi), rp)
         s(2, 1) = -conjg(s(1, 2))
         s(2, 2) = conjg(s(1, 1))
         rot = (0.0_rp, 0.0_rp)
         do i = 1, 2
            do k = 1, 2
               do a = 1, nd
                  rot((i - 1)*nd + a, (k - 1)*nd + a) = s(i, k)
               end do
            end do
         end do
         p = rand_physical()
         pr = matmul(rot, matmul(p, conjg(transpose(rot))))
         call hubbard_spinor_functional(t, u, jbar, p, e_int, e_dc, e_hub, v)
         call hubbard_spinor_functional(t, u, jbar, pr, er_int, er_dc, er_hub, vr)
         worst_e = max(worst_e, abs(er_hub - e_hub))
         worst_v = max(worst_v, max_abs(vr - matmul(rot, matmul(v, conjg(transpose(rot))))))
         ! control: V of the rotated density compared with the unrotated V
         worst_c_v = max(worst_c_v, max_abs(vr - v))
      end do
      call check('(b) SU(2) invariance: E_Hub unchanged', worst_e, 1.0e-12_rp)
      call check('(b) SU(2) covariance: V -> R V R^dagger', worst_v, 1.0e-12_rp)
      call check('(b) control, rotate P without rotating V', worst_c_v, 1.0e-12_rp, .true.)
      call test_su2_energy_control(t, u, jbar)
   end subroutine test_su2_invariance

   !> Control for the energy: rotating only the spin-up block (a non-SU(2) transformation)
   !> must change E_Hub.
   subroutine test_su2_energy_control(t, u, jbar)
      complex(rp), intent(in) :: t(nd, nd, nd, nd)
      real(rp), intent(in) :: u, jbar
      complex(rp) :: p(n2, n2), pr(n2, n2), v(n2, n2), rot(n2, n2)
      real(rp) :: e_int, e_dc, e_hub, er_int, er_dc, er_hub
      integer :: a

      p = rand_physical()
      rot = (0.0_rp, 0.0_rp)
      do a = 1, n2
         rot(a, a) = (1.0_rp, 0.0_rp)
      end do
      do a = 1, nd
         rot(a, a) = cmplx(0.0_rp, 1.0_rp, rp)      ! phase on the up spin only
      end do
      rot(1, 1) = (0.0_rp, 0.0_rp)
      rot(1, 2) = (1.0_rp, 0.0_rp)
      rot(2, 1) = (1.0_rp, 0.0_rp)
      rot(2, 2) = (0.0_rp, 0.0_rp)                   ! plus an orbital swap within the up block
      pr = matmul(rot, matmul(p, conjg(transpose(rot))))
      call hubbard_spinor_functional(t, u, jbar, p, e_int, e_dc, e_hub, v)
      call hubbard_spinor_functional(t, u, jbar, pr, er_int, er_dc, er_hub, v)
      call check('(b) control, up-block-only orbital swap is not a symmetry: E_Hub', abs(er_hub - e_hub), 1.0e-12_rp, .true.)
   end subroutine test_su2_energy_control

   subroutine test_derivative()
      real(rp), parameter :: u = 1.4_rp, jbar = 0.22_rp
      real(rp), parameter :: steps(3) = [3.0e-1_rp, 1.0e-1_rp, 3.0e-2_rp]
      complex(rp) :: t(nd, nd, nd, nd), p(n2, n2), dp(n2, n2), v(n2, n2), vscratch(n2, n2)
      real(rp) :: e0(3), ep(3), em(3), slope, slope_t, fd, worst_fd, quad(3), worst_quad, worst_ctrl, dummy1, dummy2
      integer :: n, i, k

      call hubbard_d_tensor(u, jbar, t)
      p = rand_physical()
      dp = rand_hermitian(1.0_rp)
      dp = dp/sqrt(sum(abs(dp)**2))
      call hubbard_spinor_functional(t, u, jbar, p, dummy1, dummy2, e0(1), v)
      slope = 0.0_rp
      slope_t = 0.0_rp
      do i = 1, n2
         do k = 1, n2
            slope = slope + real(v(i, k)*dp(k, i), rp)
            slope_t = slope_t + real(v(k, i)*dp(k, i), rp)       ! Tr(V^T dP)
         end do
      end do
      worst_fd = 0.0_rp
      worst_ctrl = huge(1.0_rp)
      do n = 1, 3
         call hubbard_spinor_functional(t, u, jbar, p + steps(n)*dp, dummy1, dummy2, ep(n), vscratch)
         call hubbard_spinor_functional(t, u, jbar, p - steps(n)*dp, dummy1, dummy2, em(n), vscratch)
         fd = (ep(n) - em(n))/(2.0_rp*steps(n))
         worst_fd = max(worst_fd, abs(fd - slope)/abs(slope))
         worst_ctrl = min(worst_ctrl, abs(fd - slope_t)/abs(slope))
         ! E is exactly quadratic in P: E(P+h dP) - E(P) - h Tr(V dP) = h^2 Q for every h
         quad(n) = (ep(n) - e0(1) - steps(n)*slope)/steps(n)**2
      end do
      worst_quad = max(abs(quad(2) - quad(1)), abs(quad(3) - quad(1)))/max(abs(quad(1)), tiny(1.0_rp))
      call check('(c) central difference vs Tr(V dP), relative, h=0.3,0.1,0.03', worst_fd, 1.0e-8_rp)
      call check('(c) remainder/h^2 constant across the three steps, relative', worst_quad, 1.0e-8_rp)
      call check('(c) control, transposed V: relative mismatch', worst_ctrl, 1.0e-8_rp, .true.)
   end subroutine test_derivative

   subroutine test_dudarev_limit()
      real(rp), parameter :: u = 1.9_rp
      complex(rp) :: t(nd, nd, nd, nd), p(n2, n2), v(n2, n2), vref(n2, n2), eye(n2, n2)
      real(rp) :: e_int, e_dc, e_hub, eref, eblocks
      integer :: i

      call hubbard_d_tensor(u, 0.0_rp, t)
      eye = (0.0_rp, 0.0_rp)
      do i = 1, n2
         eye(i, i) = (1.0_rp, 0.0_rp)
      end do
      p = rand_hermitian(1.2_rp)        ! arbitrary Hermitian spinor density, not necessarily physical
      call hubbard_spinor_functional(t, u, 0.0_rp, p, e_int, e_dc, e_hub, v)
      eref = 0.5_rp*u*real(sum([(p(i, i), i=1, n2)]) - sum(p*transpose(p)), rp)
      vref = u*(0.5_rp*eye - p)
      call check('(d) J=0 Dudarev: E_Hub = U/2 Tr(P - P^2)', abs(e_hub - eref), 1.0e-13_rp)
      call check('(d) J=0 Dudarev: V = U(1/2 - P)', max_abs(v - vref), 1.0e-13_rp)
      ! control: the same formula applied to the two 5x5 spin blocks separately (drops the coherences)
      eblocks = 0.5_rp*u*real(sum([(p(i, i), i=1, n2)]) - sum(p(1:nd, 1:nd)*transpose(p(1:nd, 1:nd))) &
         - sum(p(nd + 1:n2, nd + 1:n2)*transpose(p(nd + 1:n2, nd + 1:n2))), rp)
      call check('(d) control, block-wise Dudarev ignoring spin coherences: E_Hub', abs(e_hub - eblocks), 1.0e-13_rp, .true.)
   end subroutine test_dudarev_limit

   subroutine test_single_projector()
      real(rp), parameter :: u = 1.9_rp
      complex(rp) :: t(nd, nd, nd, nd), psi(n2), p(n2, n2), v(n2, n2), q(n2, n2)
      real(rp) :: e_int, e_dc, e_hub, f_int, f_dc, f_hub
      integer :: i, k

      call hubbard_d_tensor(u, 0.0_rp, t)
      q = rand_unitary(n2)
      psi = q(:, 1)
      do i = 1, n2
         do k = 1, n2
            p(i, k) = psi(i)*conjg(psi(k))
         end do
      end do
      call hubbard_spinor_functional(t, u, 0.0_rp, p, e_int, e_dc, e_hub, v)
      call check('(e) single occupied spinor, J=0: E_Hub', abs(e_hub), 1.0e-14_rp)
      call check('(e) single occupied spinor, J=0: E_int', abs(e_int), 1.0e-14_rp)
      call hubbard_spinor_functional(t, u, 0.0_rp, 0.5_rp*p, f_int, f_dc, f_hub, v)
      call check('(e) control, half-occupied spinor (E_Hub = U/8): E_Hub', abs(f_hub), 1.0e-14_rp, .true.)
   end subroutine test_single_projector

   subroutine test_r08_oracle()
      complex(rp) :: t(nd, nd, nd, nd), p(n2, n2), vref(n2, n2), v(n2, n2)
      real(rp) :: u, jbar, eref(3), e_int, e_dc, e_hub, ctrl_e, ctrl_v
      real(rp) :: re, im
      integer :: unit, ncase, icase, i, k, ii, kk, ios
      character(len=48) :: name

      open (newunit=unit, file=trim(reference_path), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write (*, '(a)') 'FAIL (f) cannot open R08 reference file '//trim(reference_path)
         failed = .true.
         return
      end if
      read (unit, *) ncase
      do icase = 1, ncase
         read (unit, *) u, jbar
         do i = 1, n2
            do k = 1, n2
               read (unit, *) ii, kk, re, im
               p(ii, kk) = cmplx(re, im, rp)
            end do
         end do
         read (unit, *) eref
         do i = 1, n2
            do k = 1, n2
               read (unit, *) ii, kk, re, im
               vref(ii, kk) = cmplx(re, im, rp)
            end do
         end do
         call hubbard_d_tensor(u, jbar, t)
         call hubbard_spinor_functional(t, u, jbar, p, e_int, e_dc, e_hub, v)
         write (name, '(a,i0,a,f5.2)') '(f) R08 case ', icase, ' J=', jbar
         call check(trim(name)//': E_int, E_DC, E_Hub', maxval(abs([e_int, e_dc, e_hub] - eref)), 1.0e-12_rp)
         call check(trim(name)//': V', max_abs(v - vref), 1.0e-12_rp)
         if (icase == 2) then
            ! control: evaluate this density with the J of the first case (J=0) against this reference
            call hubbard_d_tensor(u, 0.0_rp, t)
            call hubbard_spinor_functional(t, u, 0.0_rp, p, e_int, e_dc, e_hub, v)
            ctrl_e = maxval(abs([e_int, e_dc, e_hub] - eref))
            ctrl_v = max_abs(v - vref)
         end if
      end do
      close (unit)
      call check('(f) control, wrong J against the J>0 reference: energies', ctrl_e, 1.0e-12_rp, .true.)
      call check('(f) control, wrong J against the J>0 reference: V', ctrl_v, 1.0e-12_rp, .true.)
   end subroutine test_r08_oracle

   subroutine test_hermitian_potential()
      real(rp), parameter :: u = 2.3_rp, jbar = 0.41_rp
      complex(rp) :: t(nd, nd, nd, nd), p(n2, n2), v(n2, n2), skew(n2, n2)
      real(rp) :: e_int, e_dc, e_hub, worst
      integer :: trial

      call hubbard_d_tensor(u, jbar, t)
      worst = 0.0_rp
      do trial = 1, 5
         p = rand_physical()
         call hubbard_spinor_functional(t, u, jbar, p, e_int, e_dc, e_hub, v)
         worst = max(worst, max_abs(v - conjg(transpose(v))))
      end do
      call check('(g) V Hermitian for physical spinor P', worst, 1.0e-14_rp)
      skew = rand_hermitian(0.3_rp)
      p = rand_physical() + cmplx(0.0_rp, 1.0_rp, rp)*skew      ! add an anti-Hermitian part
      call hubbard_spinor_functional(t, u, jbar, p, e_int, e_dc, e_hub, v)
      call check('(g) control, non-Hermitian P gives non-Hermitian V', max_abs(v - conjg(transpose(v))), 1.0e-14_rp, .true.)
   end subroutine test_hermitian_potential

end program test_hubbard_spinor
