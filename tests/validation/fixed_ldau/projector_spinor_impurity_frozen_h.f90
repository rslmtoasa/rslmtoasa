! S6 T2(a): nsp=2 impurity local moments, one frozen U-carrying host plus one recursion atom.
! S6c T2(a): the same checks on a surface-like layout (calctype S, second half of the program):
! a frozen U host (type 1), an empty sphere without U (type 2) and two surface types with U (types 3, 4),
! three recursion atoms with irec(k)=k+1 and the legacy ib map deliberately pointing at the host.
!
! The recursion atom (type 2) has H = d crystal field + exchange + xi L.S + spin mixing + the spinor V_U
! built from its own P_in; the frozen host (type 1) carries a different P and U. Deliberately ib(1)=1
! while irec(1)=2 and iz(2)=2, so the legacy host map would pick the host's V. The production RS routine
! bands%calculate_moments is run in spectral and in local mode on the same Green function. Expected
! c'_s and q2'_s: the spinor formulas of SPINOR_DM_CONTRACT.md evaluated here on the spectral run's
! P, M1, c, q2 with V from the S1 functional (1e-12 / 1e-10); c'_s also from the exact eigenpairs of H
! with Lorentzian occupations (RS quadrature budget). The host must be untouched.
program projector_spinor_impurity_frozen_h
   use precision_mod, only: rp
   use basis_mod, only: basis_init, nb, norb, spin_off
   use math_mod, only: init_math_operators, pi
   use hubbard_spinor_mod, only: hubbard_d_tensor, hubbard_spinor_functional
   use logger_mod, only: g_logger
   use mpi_mod, only: rank, ierr, atoms_per_process, get_mpi_variables, parallel_context, assignment(=)
   use control_mod, only: control
   use lattice_mod, only: lattice
   use charge_mod, only: charge
   use hamiltonian_mod, only: hamiltonian
   use energy_mod, only: energy
   use recursion_mod, only: recursion
   use density_of_states_mod, only: dos
   use green_mod, only: green
   use bands_mod, only: bands
#ifdef USE_MPI
   use mpi
#endif
   implicit none

   integer, parameter :: nch = 6000
   real(rp), parameter :: ef = 0.0_rp, eta = 0.004_rp, emin = -3.0_rp
   real(rp), parameter :: u_d = 0.25_rp, j_d = 0.04_rp, tol_rs_c = 3.0e-4_rp

   type(parallel_context) :: parallel
   type(control), target :: ctl
   type(lattice), target :: lat
   type(charge), target :: chg
   type(hamiltonian), target :: ham
   type(energy), target :: en
   type(recursion), target :: rec
   type(dos), target :: dos_state
   type(green), target :: gr
   type(bands), allocatable :: rs
   complex(rp), allocatable :: q(:, :), h(:, :), qa(:, :, :)
   real(rp), allocatable :: eig(:), eiga(:, :)
   complex(rp) :: p_rec(10, 10), p_host(10, 10), v_rec(10, 10), v_host(10, 10), v_ctl(10, 10)
   complex(rp) :: v_prod_rec(10, 10), v_prod_host(10, 10), p_scratch(10, 10)
   complex(rp) :: p_spec(10, 10), m1_spec(10, 10), t5(5, 5, 5, 5)
   real(rp) :: c_spec(2), q2_spec(2), q0(2), c_loc(2), q2_loc(2)
   real(rp) :: c_exp(2), q2_exp(2), c_ctl(2), q2_ctl(2), c_dir(2)
   real(rp) :: gap, gap_best, w
   complex(rp) :: host_p_before(10, 10)
   real(rp) :: host_c_before(2), host_q0_before(2), host_q2_before(2)
   logical :: failed = .false.
   integer :: s
   !> Model parameters (xi, spin coupling, d level shift) of types 1..4 of the surface layout.
   real(rp), parameter :: model_xi(4) = [0.0_rp, 0.05_rp, 0.35_rp, 0.2_rp]
   real(rp), parameter :: model_cpl(4) = [0.0_rp, 0.3_rp, 1.0_rp, 0.7_rp]
   real(rp), parameter :: model_shift(4) = [0.0_rp, 0.2_rp, 0.0_rp, 0.15_rp]

   external :: zheev
#ifdef USE_MPI
   call MPI_INIT(ierr)
#endif
   parallel = parallel_context()
   call g_logger%init()
   call basis_init(2)
   call init_math_operators()
   allocate (q(nb, nb), h(nb, nb), eig(nb), qa(nb, nb, 3), eiga(nb, 3))
   call setup_objects('I', 2)
   call hubbard_d_tensor(u_d, j_d, t5)

   ! P_in of the recursion atom: the full model; of the host: the model without spin mixing.
   call build_model(0.35_rp, 1.0_rp)
   call exact_references(p_rec)
   call build_model(0.0_rp, 0.0_rp)
   call exact_references(p_host)
   lat%symbolic_atoms(1)%potential%ldm_spinor = p_host
   lat%symbolic_atoms(2)%potential%ldm_spinor = p_rec
   call ham%calculate_hubbard_u_potential_general()
   call hubbard_spinor_functional(t5, u_d, j_d, p_rec, gap, gap_best, w, v_rec)
   call hubbard_spinor_functional(t5, u_d, j_d, p_host, gap, gap_best, w, v_host)
   call spinor_block(ham%hubbard_u_pot(:, :, 2), v_prod_rec)
   call spinor_block(ham%hubbard_u_pot(:, :, 1), v_prod_host)
   call report('inserted V of the recursion atom vs S1 functional of its P_in', maxval(abs(v_prod_rec - v_rec)), 1.0e-14_rp)
   call report('inserted V of the host vs S1 functional of its P', maxval(abs(v_prod_host - v_host)), 1.0e-14_rp)
   call report('max |V_up,dn| of the recursion atom (nonzero)', maxval(abs(v_rec(1:5, 6:10))), 1.0e-3_rp, .true.)
   call report('host V differs from recursion V (site map is observable)', maxval(abs(v_host - v_rec)), 1.0e-3_rp, .true.)

   ! H of the recursion atom including its inserted V; exact eigenpairs, E_F in the widest gap.
   call build_model(0.35_rp, 1.0_rp)
   call add_v(v_prod_rec)
   call exact_references(p_scratch)
   qa(:, :, 1) = q
   eiga(:, 1) = eig
   call setup_rs_green()
   host_p_before = lat%symbolic_atoms(1)%potential%ldm_spinor
   host_c_before = lat%symbolic_atoms(1)%potential%gravity_center(3, :)
   host_q0_before = lat%symbolic_atoms(1)%potential%ql(1, 2, :)
   host_q2_before = lat%symbolic_atoms(1)%potential%ql(3, 2, :)

   ham%hubbard_radial_moments = 'spectral'
   call rs%calculate_moments()
   p_spec = lat%symbolic_atoms(2)%potential%ldm_spinor
   m1_spec = lat%symbolic_atoms(2)%potential%m1_spinor
   do s = 1, 2
      q0(s) = lat%symbolic_atoms(2)%potential%ql(1, 2, s)
      c_spec(s) = lat%symbolic_atoms(2)%potential%gravity_center(3, s) + lat%symbolic_atoms(2)%potential%vmad
      q2_spec(s) = lat%symbolic_atoms(2)%potential%ql(3, 2, s)
   end do
   ham%hubbard_radial_moments = 'local'
   call rs%calculate_moments()
   do s = 1, 2
      c_loc(s) = lat%symbolic_atoms(2)%potential%gravity_center(3, s) + lat%symbolic_atoms(2)%potential%vmad
      q2_loc(s) = lat%symbolic_atoms(2)%potential%ql(3, 2, s)
   end do

   call formulas(v_rec, p_spec, m1_spec, q0, c_spec, q2_spec, c_exp, q2_exp)
   call report('(a) local c prime vs spinor formula on the spectral quantities', maxval(abs(c_loc - c_exp)), 1.0e-12_rp)
   call report('(a) local q2 prime vs spinor formula on the spectral quantities', maxval(abs(q2_loc - q2_exp)), 1.0e-10_rp)
   call report('(a) size of the c prime shift (resolvable)', maxval(abs(c_loc - c_spec)), 1.0e-3_rp, .true.)
   call exact_c_prime(v_rec, 1, c_dir)
   call report('(a) local c prime vs exact eigenpairs with Lorentzian occupations (RS quadrature)', &
      maxval(abs(c_loc - c_dir)), tol_rs_c)

   call report('(a) host skipped: P unchanged', maxval(abs(lat%symbolic_atoms(1)%potential%ldm_spinor - host_p_before)), 0.0_rp)
   call report('(a) host skipped: centre unchanged', &
      maxval(abs(lat%symbolic_atoms(1)%potential%gravity_center(3, :) - host_c_before)), 0.0_rp)
   call report('(a) host skipped: q0, q2 unchanged', &
      max(maxval(abs(lat%symbolic_atoms(1)%potential%ql(1, 2, :) - host_q0_before)), &
          maxval(abs(lat%symbolic_atoms(1)%potential%ql(3, 2, :) - host_q2_before))), 0.0_rp)

   call formulas(v_host, p_spec, m1_spec, q0, c_spec, q2_spec, c_ctl, q2_ctl)
   call report('(a) control, host V (legacy ib map): c prime', maxval(abs(c_loc - c_ctl)), 1.0e-6_rp, .true.)
   call report('(a) control, host V (legacy ib map): q2 prime', maxval(abs(q2_loc - q2_ctl)), 1.0e-6_rp, .true.)
   v_ctl = v_rec
   v_ctl(1:5, 6:10) = (0.0_rp, 0.0_rp)
   v_ctl(6:10, 1:5) = (0.0_rp, 0.0_rp)
   call formulas(v_ctl, p_spec, m1_spec, q0, c_spec, q2_spec, c_ctl, q2_ctl)
   call report('(a) control, spin-off-diagonal V blocks dropped: c prime', maxval(abs(c_loc - c_ctl)), 1.0e-6_rp, .true.)
   call report('(a) control, spin-off-diagonal V blocks dropped: q2 prime', maxval(abs(q2_loc - q2_ctl)), 1.0e-6_rp, .true.)

   call surface_case()

#ifdef USE_MPI
   call reduce_failure()
#endif
   if (rank == 0) then
      if (failed) then
         write (*, '(a)') 'projector_spinor_impurity_frozen_h: FAIL'
      else
         write (*, '(a)') 'projector_spinor_impurity_frozen_h: PASS'
      end if
   end if
#ifdef USE_MPI
   call parallel%restore_to_default()
   call MPI_FINALIZE(ierr)
#endif
   if (failed) error stop 1

contains

   subroutine report(name, measured, tol, control_check)
      character(len=*), intent(in) :: name
      real(rp), intent(in) :: measured, tol
      logical, intent(in), optional :: control_check
      logical :: ok, is_control

      is_control = .false.
      if (present(control_check)) is_control = control_check
      ok = merge(measured > tol, measured <= tol, is_control)
      if (.not. ok) failed = .true.
      if (rank == 0) then
         if (is_control) then
            write (*, '(a,1x,a,2(2x,es10.3))') merge('PASS', 'FAIL', ok), trim(name)//' [must exceed]', measured, tol
         else
            write (*, '(a,1x,a,2(2x,es10.3))') merge('PASS', 'FAIL', ok), trim(name), measured, tol
         end if
      end if
   end subroutine report

   !> d block of an nb x nb Hamiltonian-space array, index (sigma-1)*5+(m+3).
   subroutine spinor_block(a, b)
      complex(rp), intent(in) :: a(nb, nb)
      complex(rp), intent(out) :: b(10, 10)
      integer :: s, t

      do s = 1, 2
         do t = 1, 2
            b((s - 1)*5 + 1:s*5, (t - 1)*5 + 1:t*5) = a((s - 1)*spin_off + 5:(s - 1)*spin_off + 9, &
                                                       (t - 1)*spin_off + 5:(t - 1)*spin_off + 9)
         end do
      end do
   end subroutine spinor_block

   subroutine add_v(v)
      complex(rp), intent(in) :: v(10, 10)
      integer :: s, t

      do s = 1, 2
         do t = 1, 2
            h((s - 1)*norb + 5:(s - 1)*norb + 9, (t - 1)*norb + 5:(t - 1)*norb + 9) = &
               h((s - 1)*norb + 5:(s - 1)*norb + 9, (t - 1)*norb + 5:(t - 1)*norb + 9) + &
               v((s - 1)*5 + 1:s*5, (t - 1)*5 + 1:t*5)
         end do
      end do
   end subroutine add_v

   !> c'_s = c_s - t_s/q0_s and q2'_s = q2_s - 2 Re Tr(Pi_s V M1c_s) + Tr(V Pi_s V P) - t_s^2/q0_s.
   subroutine formulas(v, p, m1, q0s, cs, q2s, cp, q2p)
      complex(rp), intent(in) :: v(10, 10), p(10, 10), m1(10, 10)
      real(rp), intent(in) :: q0s(2), cs(2), q2s(2)
      real(rp), intent(out) :: cp(2), q2p(2)
      complex(rp) :: vp(10, 10), vpv(10, 10), vm1c(10, 10)
      real(rp) :: t
      integer :: s, k

      vp = matmul(v, p)
      vpv = matmul(vp, v)
      do s = 1, 2
         vm1c = matmul(v, m1 - cs(s)*p)
         t = 0.0_rp
         q2p(s) = q2s(s)
         do k = (s - 1)*5 + 1, s*5
            t = t + real(vp(k, k), rp)
            q2p(s) = q2p(s) - 2.0_rp*real(vm1c(k, k), rp) + real(vpv(k, k), rp)
         end do
         cp(s) = cs(s) - t/q0s(s)
         q2p(s) = q2p(s) - t*t/q0s(s)
      end do
   end subroutine formulas

   !> c'_s = [sum_n occ1_n <a_n|Pi_s|a_n> eps-weighted - sum_n occ_n Re <a_n|Pi_s V|a_n>] / sum_n occ_n <a_n|Pi_s|a_n> from the eigenpairs of H
   !> (H already contains V), with the finite-window Lorentzian occupations of the RS quadrature.
   subroutine exact_c_prime(v, ka, cp)
      complex(rp), intent(in) :: v(10, 10)
      integer, intent(in) :: ka
      real(rp), intent(out) :: cp(2)
      complex(rp) :: a(10), va(10)
      real(rp) :: occ, occ_m1, num(2), den(2)
      integer :: k, s, sr, i

      num = 0.0_rp
      den = 0.0_rp
      do k = 1, nb
         occ = (atan((ef - eiga(k, ka))/eta) - atan((emin - eiga(k, ka))/eta))/pi
         occ_m1 = eiga(k, ka)*occ + eta/(2.0_rp*pi)*log(((ef - eiga(k, ka))**2 + eta**2)/((emin - eiga(k, ka))**2 + eta**2))
         do sr = 1, 2
            do i = 1, 5
               a((sr - 1)*5 + i) = qa((sr - 1)*norb + 4 + i, k, ka)
            end do
         end do
         va = matmul(v, a)
         do s = 1, 2
            den(s) = den(s) + occ*sum(abs(a((s - 1)*5 + 1:s*5))**2)
            num(s) = num(s) + occ_m1*sum(abs(a((s - 1)*5 + 1:s*5))**2) - &
                     occ*real(dot_product(a((s - 1)*5 + 1:s*5), va((s - 1)*5 + 1:s*5)), rp)
         end do
      end do
      cp = num/den
   end subroutine exact_c_prime

   subroutine build_model(xi, coupling, dshift)
      real(rp), intent(in) :: xi, coupling
      real(rp), intent(in), optional :: dshift
      complex(rp) :: hd(5, 5), lz(5, 5), lp(5, 5), x(5, 5), e5(5, 5)
      real(rp) :: lm
      real(rp), parameter :: es(4) = [-1.4_rp, -1.0_rp, -0.8_rp, -0.6_rp]
      real(rp), parameter :: ed(5) = [-0.34_rp, -0.19_rp, 0.16_rp, 0.39_rp, 0.61_rp]
      integer :: i, j, m

      h = (0.0_rp, 0.0_rp)
      e5 = (0.0_rp, 0.0_rp)
      do i = 1, 4
         h(i, i) = es(i)
         h(norb + i, norb + i) = es(i) + 0.05_rp
      end do
      hd = (0.0_rp, 0.0_rp)
      lz = (0.0_rp, 0.0_rp)
      lp = (0.0_rp, 0.0_rp)
      x = (0.0_rp, 0.0_rp)
      do i = 1, 5
         e5(i, i) = (1.0_rp, 0.0_rp)
         hd(i, i) = ed(i)
         lz(i, i) = real(i - 3, rp)
         do j = i + 1, 5
            hd(i, j) = cmplx(0.021_rp*real(i + j, rp), 0.013_rp*real(j - i, rp), rp)
            hd(j, i) = conjg(hd(i, j))
         end do
         do j = 1, 5
            x(i, j) = cmplx(0.031_rp*real(i - j, rp) + 0.017_rp, 0.023_rp*real(i + 2*j, rp) - 0.05_rp, rp)
         end do
      end do
      do m = -2, 1
         lm = sqrt(6.0_rp - real(m*(m + 1), rp))
         lp(m + 4, m + 3) = cmplx(lm, 0.0_rp, rp)
      end do
      if (present(dshift)) hd = hd + dshift*e5
      h(5:9, 5:9) = hd - 0.07_rp*e5 + 0.5_rp*xi*lz
      h(norb + 5:nb, norb + 5:nb) = hd + 0.07_rp*e5 - 0.5_rp*xi*lz
      h(5:9, norb + 5:nb) = 0.5_rp*xi*conjg(transpose(lp)) + coupling*x
      h(norb + 5:nb, 5:9) = 0.5_rp*xi*lp + coupling*conjg(transpose(x))
   end subroutine build_model

   !> Eigenpairs of the current h (E_F placed in the widest gap between levels 6 and 12) and the exact
   !> occupied 10x10 d projector.
   subroutine exact_references(p_out)
      complex(rp), intent(out) :: p_out(10, 10)
      complex(rp), allocatable :: work(:)
      real(rp), allocatable :: rwork(:)
      real(rp) :: gap_b, shift
      integer :: info, a, b, k, row, col, sa, sb

      q = h
      allocate (work(2*nb), rwork(3*nb - 2))
      call zheev('V', 'U', nb, q, nb, eig, work, 2*nb, rwork, info)
      if (info /= 0) failed = .true.
      gap_b = 0.0_rp
      shift = 0.0_rp
      do k = 6, 11
         if (eig(k + 1) - eig(k) > gap_b) then
            gap_b = eig(k + 1) - eig(k)
            shift = 0.5_rp*(eig(k) + eig(k + 1))
         end if
      end do
      eig = eig - shift
      if (gap_b < 0.1_rp) failed = .true.
      p_out = (0.0_rp, 0.0_rp)
      do sa = 1, 2
         do sb = 1, 2
            do a = 1, 5
               do b = 1, 5
                  row = (sa - 1)*norb + 4 + a
                  col = (sb - 1)*norb + 4 + b
                  do k = 1, nb
                     if (eig(k) < ef) p_out((sa - 1)*5 + a, (sb - 1)*5 + b) = &
                        p_out((sa - 1)*5 + a, (sb - 1)*5 + b) + q(row, k)*conjg(q(col, k))
                  end do
               end do
            end do
         end do
      end do
      deallocate (work, rwork)
   end subroutine exact_references

   !> Surface-like layout: type 1 frozen U host, type 2 empty sphere (no U), types 3 and 4 surface layers with U,
   !> each recursion atom k (type k+1) with its own H (crystal field, exchange, xi L.S, spin mixing, own V).
   subroutine surface_case()
      complex(rp) :: p_in(10, 10, 4), v_t(10, 10, 4), p_sp(10, 10, 3), m1_sp(10, 10, 3), vv(10, 10), vctl(10, 10)
      complex(rp) :: host_p_before(10, 10)
      real(rp) :: q0a(2, 3), ca(2, 3), q2a(2, 3), cl(2, 3), q2l(2, 3), ce(2), q2e(2), cd(2), cc(2), q2c(2)
      real(rp) :: host_c_before(2), host_q0_before(2), host_q2_before(2), e_i, e_d, e_h
      character(len=40) :: tag
      integer :: k, t, sp

      call setup_objects('S', 4)
      do t = 1, 4
         call build_model(model_xi(t), model_cpl(t), model_shift(t))
         call exact_references(p_in(:, :, t))
         lat%symbolic_atoms(t)%potential%ldm_spinor = p_in(:, :, t)
      end do
      call ham%calculate_hubbard_u_potential_general()
      do t = 1, 4
         write (tag, '(a,i0)') '(S) type ', t
         if (t == 2) then
            call report(trim(tag)//' empty sphere: inserted V is zero', maxval(abs(ham%hubbard_u_pot(:, :, 2))), 0.0_rp)
            v_t(:, :, 2) = (0.0_rp, 0.0_rp)
            cycle
         end if
         call hubbard_spinor_functional(t5, u_d, j_d, p_in(:, :, t), e_i, e_d, e_h, v_t(:, :, t))
         call spinor_block(ham%hubbard_u_pot(:, :, t), vv)
         call report(trim(tag)//' inserted V vs S1 functional of its P_in', maxval(abs(vv - v_t(:, :, t))), 1.0e-14_rp)
      end do
      call report('(S) max |V_up,dn| of type 3 (nonzero)', maxval(abs(v_t(1:5, 6:10, 3))), 1.0e-3_rp, .true.)
      call report('(S) max |V_up,dn| of type 4 (nonzero)', maxval(abs(v_t(1:5, 6:10, 4))), 1.0e-3_rp, .true.)
      call report('(S) V of type 3 differs from V of type 4 (site map is observable)', &
         maxval(abs(v_t(:, :, 3) - v_t(:, :, 4))), 1.0e-3_rp, .true.)

      ! H of recursion atom k (type k+1) with its own inserted V; exact eigenpairs, E_F in the widest gap.
      do k = 1, 3
         t = k + 1
         call build_model(model_xi(t), model_cpl(t), model_shift(t))
         if (t /= 2) call add_v(v_t(:, :, t))
         call exact_references(vv)
         qa(:, :, k) = q
         eiga(:, k) = eig
      end do
      call setup_rs_green()
      host_p_before = lat%symbolic_atoms(1)%potential%ldm_spinor
      host_c_before = lat%symbolic_atoms(1)%potential%gravity_center(3, :)
      host_q0_before = lat%symbolic_atoms(1)%potential%ql(1, 2, :)
      host_q2_before = lat%symbolic_atoms(1)%potential%ql(3, 2, :)

      ham%hubbard_radial_moments = 'spectral'
      call rs%calculate_moments()
      do k = 1, 3
         p_sp(:, :, k) = lat%symbolic_atoms(k + 1)%potential%ldm_spinor
         m1_sp(:, :, k) = lat%symbolic_atoms(k + 1)%potential%m1_spinor
         do sp = 1, 2
            q0a(sp, k) = lat%symbolic_atoms(k + 1)%potential%ql(1, 2, sp)
            ca(sp, k) = lat%symbolic_atoms(k + 1)%potential%gravity_center(3, sp) + lat%symbolic_atoms(k + 1)%potential%vmad
            q2a(sp, k) = lat%symbolic_atoms(k + 1)%potential%ql(3, 2, sp)
         end do
      end do
      ham%hubbard_radial_moments = 'local'
      call rs%calculate_moments()
      do k = 1, 3
         do sp = 1, 2
            cl(sp, k) = lat%symbolic_atoms(k + 1)%potential%gravity_center(3, sp) + lat%symbolic_atoms(k + 1)%potential%vmad
            q2l(sp, k) = lat%symbolic_atoms(k + 1)%potential%ql(3, 2, sp)
         end do
      end do

      do k = 1, 3
         t = k + 1
         write (tag, '(a,i0)') '(S) type ', t
         if (t == 2) then
            call report(trim(tag)//' empty sphere skipped: c, q2 equal the spectral run', &
               max(maxval(abs(cl(:, k) - ca(:, k))), maxval(abs(q2l(:, k) - q2a(:, k)))), 0.0_rp)
            call formulas(v_t(:, :, 1), p_sp(:, :, k), m1_sp(:, :, k), q0a(:, k), ca(:, k), q2a(:, k), cc, q2c)
            call report(trim(tag)//' control, host V applied to the empty sphere: c prime', &
               maxval(abs(cl(:, k) - cc)), 1.0e-6_rp, .true.)
            cycle
         end if
         call formulas(v_t(:, :, t), p_sp(:, :, k), m1_sp(:, :, k), q0a(:, k), ca(:, k), q2a(:, k), ce, q2e)
         call report(trim(tag)//' local c prime vs spinor formula on the spectral quantities', &
            maxval(abs(cl(:, k) - ce)), 1.0e-12_rp)
         call report(trim(tag)//' local q2 prime vs spinor formula on the spectral quantities', &
            maxval(abs(q2l(:, k) - q2e)), 1.0e-10_rp)
         call report(trim(tag)//' size of the c prime shift (resolvable)', maxval(abs(cl(:, k) - ca(:, k))), 1.0e-3_rp, .true.)
         call exact_c_prime(v_t(:, :, t), k, cd)
         call report(trim(tag)//' local c prime vs exact eigenpairs with Lorentzian occupations (RS quadrature)', &
            maxval(abs(cl(:, k) - cd)), tol_rs_c)
         ! site swap: the V of the other surface type
         call formulas(v_t(:, :, 7 - t), p_sp(:, :, k), m1_sp(:, :, k), q0a(:, k), ca(:, k), q2a(:, k), cc, q2c)
         call report(trim(tag)//' control, site swap (V of the other surface type): c prime', &
            maxval(abs(cl(:, k) - cc)), 1.0e-6_rp, .true.)
         call report(trim(tag)//' control, site swap (V of the other surface type): q2 prime', &
            maxval(abs(q2l(:, k) - q2c)), 1.0e-6_rp, .true.)
         call formulas(v_t(:, :, 1), p_sp(:, :, k), m1_sp(:, :, k), q0a(:, k), ca(:, k), q2a(:, k), cc, q2c)
         call report(trim(tag)//' control, host V (legacy ib map): c prime', maxval(abs(cl(:, k) - cc)), 1.0e-6_rp, .true.)
         vctl = v_t(:, :, t)
         vctl(1:5, 6:10) = (0.0_rp, 0.0_rp)
         vctl(6:10, 1:5) = (0.0_rp, 0.0_rp)
         call formulas(vctl, p_sp(:, :, k), m1_sp(:, :, k), q0a(:, k), ca(:, k), q2a(:, k), cc, q2c)
         call report(trim(tag)//' control, spin-off-diagonal V blocks dropped: c prime', &
            maxval(abs(cl(:, k) - cc)), 1.0e-6_rp, .true.)
         call report(trim(tag)//' control, spin-off-diagonal V blocks dropped: q2 prime', &
            maxval(abs(q2l(:, k) - q2c)), 1.0e-6_rp, .true.)
      end do

      call report('(S) host skipped: P unchanged', maxval(abs(lat%symbolic_atoms(1)%potential%ldm_spinor - host_p_before)), 0.0_rp)
      call report('(S) host skipped: centre unchanged', &
         maxval(abs(lat%symbolic_atoms(1)%potential%gravity_center(3, :) - host_c_before)), 0.0_rp)
      call report('(S) host skipped: q0, q2 unchanged', &
         max(maxval(abs(lat%symbolic_atoms(1)%potential%ql(1, 2, :) - host_q0_before)), &
             maxval(abs(lat%symbolic_atoms(1)%potential%ql(3, 2, :) - host_q2_before))), 0.0_rp)
   end subroutine surface_case

   !> Symbolic types 1..ntype; type 1 is the frozen host (not a recursion site), types 2..ntype are recursion atoms
   !> (irec(k) = k+1). The legacy ib map points every recursion atom at the host. U on every type but 2 when
   !> calctype is S (type 2 is the empty sphere); U on both types when calctype is I.
   subroutine setup_objects(calctype, ntype)
      character(len=1), intent(in) :: calctype
      integer, intent(in) :: ntype
      integer :: i, nrec

      nrec = ntype - 1
      call get_mpi_variables(rank, nrec)
      call ctl%restore_to_default()
      ctl%lmax = 2
      ctl%nsp = 2
      ctl%calctype = calctype
      ctl%density_policy = 'constrained_spiral'

      lat%nrec = nrec; lat%nbulk = 1; lat%nbulk_bulk = 1; lat%ntype = ntype
      lat%nn_max = 1; lat%kk = 1; lat%nmax = ntype; lat%alat = 1.0_rp
      lat%r2 = 1.0_rp; lat%a = 0.0_rp; lat%a(1, 1) = 1.0_rp
      lat%a(2, 2) = 1.0_rp; lat%a(3, 3) = 1.0_rp
      lat%a_cart_inv = lat%a; lat%a_cart_inv_ready = .true.
      lat%control => ctl
      if (allocated(lat%ib)) deallocate (lat%ib, lat%irec, lat%atlist, lat%iz, lat%num, lat%nn, lat%sbar, lat%symbolic_atoms)
      allocate (lat%ib(nrec), lat%irec(nrec), lat%atlist(ntype), lat%iz(ntype), lat%num(ntype), lat%nn(ntype, 1), &
                lat%sbar(norb, norb, 1, ntype), lat%symbolic_atoms(ntype))
      ! Deliberately ib(k)=1 while irec(k)=k+1 and iz(k+1)=k+1.
      lat%ib = 1; lat%num = 1; lat%nn = 1
      lat%sbar = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, ntype
         lat%atlist(i) = i
         lat%iz(i) = i
         if (i > 1) lat%irec(i - 1) = i
      end do
      do i = 1, ntype
         call lat%symbolic_atoms(i)%restore_to_default()
         if (.not. (calctype == 'S' .and. i == 2)) then
            lat%symbolic_atoms(i)%potential%hubbard_u(3) = u_d
            lat%symbolic_atoms(i)%potential%hubbard_j(3) = j_d
         end if
         lat%symbolic_atoms(i)%potential%mom = [0.0_rp, 0.0_rp, 1.0_rp]
         lat%symbolic_atoms(i)%potential%qpar = 1.0_rp
         lat%symbolic_atoms(i)%potential%srdel = 0.1_rp
         lat%symbolic_atoms(i)%potential%c = 0.0_rp
         lat%symbolic_atoms(i)%potential%gravity_center = 0.0_rp
         lat%symbolic_atoms(i)%potential%pl = 0.0_rp
         allocate (lat%symbolic_atoms(i)%potential%ldm_spinor(10, 10))
      end do
      ! A distinctive host state, so that "untouched" is a real statement.
      lat%symbolic_atoms(1)%potential%gravity_center(3, :) = [0.123_rp, 0.234_rp]
      lat%symbolic_atoms(1)%potential%ql(1, 2, :) = [4.5_rp, 4.1_rp]
      lat%symbolic_atoms(1)%potential%ql(3, 2, :) = [0.31_rp, 0.29_rp]

      chg%lattice => lat
      chg%symbolic_atom => lat%symbolic_atoms
      ham%charge => chg; ham%lattice => lat; ham%control => ctl
      ham%local_axis = .false.
      ham%operator_generation = 0
      ham%hubbard_u_general_check = .true.
      ham%hubbard_soc_experimental = .true.
      ham%hubbard_u_potential_form = 'liechtenstein'
      ham%hubbard_u_sc_check = .false.; ham%hubbard_v_check = .false.
      if (allocated(ham%hubbard_u_pot)) deallocate (ham%hubbard_u_pot)
      allocate (ham%hubbard_u_pot(nb, nb, ntype))
      ham%hubbard_u_pot = cmplx(0.0_rp, 0.0_rp, rp)

      en%lattice => lat
      call en%restore_to_default()
      en%channels_ldos = nch
      en%energy_min = emin
      en%energy_max = 3.0_rp
      en%fermi = ef
      en%temperature = 0.0_rp
      call en%e_mesh()

      rec%hamiltonian => ham; rec%lattice => lat; rec%en => en; rec%control => ctl
      dos_state%recursion => rec; dos_state%en => en; dos_state%lattice => lat
      dos_state%control => ctl; dos_state%symbolic_atom => lat%symbolic_atoms
      gr%recursion => rec; gr%lattice => lat; gr%symbolic_atom => lat%symbolic_atoms
      gr%dos => dos_state; gr%control => ctl; gr%en => en
      if (allocated(gr%g0)) deallocate (gr%g0)
      allocate (gr%g0(nb, nb, nch + 10, atoms_per_process))
      gr%g0 = cmplx(0.0_rp, 0.0_rp, rp)

      if (allocated(rs)) deallocate (rs)
      allocate (rs)
      rs%green => gr; rs%lattice => lat; rs%symbolic_atom => lat%symbolic_atoms
      rs%dos => dos_state; rs%en => en; rs%control => ctl; rs%recursion => rec
      call rs%restore_to_default()
      rs%nv1 = nint((en%fermi - en%energy_min)/en%edel)
      rs%e1 = en%energy_min + real(rs%nv1 - 1, rp)*en%edel
   end subroutine setup_objects

   subroutine setup_rs_green()
      complex(rp), allocatable :: gmat(:, :)
      integer :: ie, band, a, b
      complex(rp) :: z
      integer :: ka

      gr%g0 = cmplx(0.0_rp, 0.0_rp, rp)
      allocate (gmat(nb, nb))
      do ka = 1, size(gr%g0, 4)
         do ie = 1, nch + 10
            gmat = cmplx(0.0_rp, 0.0_rp, rp)
            do band = 1, nb
               z = cmplx(en%ene(ie) - eiga(band, ka), eta, rp)
               do a = 1, nb
                  do b = 1, nb
                     gmat(a, b) = gmat(a, b) + qa(a, band, ka)*conjg(qa(b, band, ka))/z
                  end do
               end do
            end do
            gr%g0(:, :, ie, ka) = gmat
         end do
      end do
      deallocate (gmat)
   end subroutine setup_rs_green

#ifdef USE_MPI
   subroutine reduce_failure()
      integer :: local_flag, global_flag

      local_flag = merge(1, 0, failed)
      call MPI_ALLREDUCE(local_flag, global_flag, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      failed = global_flag /= 0
   end subroutine reduce_failure
#endif

end program projector_spinor_impurity_frozen_h
