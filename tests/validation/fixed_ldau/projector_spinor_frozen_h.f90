! Source-connected spinor (nsp=2) occupation validation on one frozen, finite 18x18 operator.
!
! H = d-shell crystal field + exchange splitting + on-site xi L.S + a spin-mixing (up-down)
! coupling. Its exact eigenpairs come from LAPACK zheev, so the oracles are exact eigenpair sums
! (KS, T=0) and the analytic finite-window Lorentzian (RS), both for the full 10x10 spinor d block
! P(a,b) = sum_n f_n psi_n(a) conj(psi_n(b)) and M1 = sum_n f_n eps_n psi_n(a) conj(psi_n(b)), with
! the index (sigma-1)*5+(m+3). The production KS routine reciprocal%calculate_spinor_dm and the
! production RS extractor bands%calculate_moments (diagnostic flag on) are compared with them.
program projector_spinor_frozen_h
   use precision_mod, only: rp
   use basis_mod, only: basis_init, nb, norb
   use math_mod, only: init_math_operators, pi, simpson_m, hubbard_spectral_density_element, &
                       hubbard_density_diagnostics
   use logger_mod, only: g_logger
   use mpi_mod, only: rank, ierr, atoms_per_process, get_mpi_variables, parallel_context, assignment(=), numprocs
   use control_mod, only: control
   use lattice_mod, only: lattice
   use charge_mod, only: charge
   use hamiltonian_mod, only: hamiltonian
   use energy_mod, only: energy
   use recursion_mod, only: recursion
   use density_of_states_mod, only: dos
   use green_mod, only: green
   use bands_mod, only: bands
   use reciprocal_mod, only: reciprocal
   use kpoint_workset_mod, only: kpoint_workset, make_kpoint_workset
#ifdef USE_MPI
   use mpi
#endif
   implicit none

   integer, parameter :: nk = 2, nch = 6000
   real(rp), parameter :: ef = 0.0_rp, eta = 0.004_rp, emin = -3.0_rp
   real(rp), parameter :: tol_ks = 5.0e-12_rp, tol_rs_p = 1.0e-4_rp, tol_rs_m1 = 3.0e-4_rp
   real(rp), parameter :: tol_zero = 1.0e-12_rp, tol_herm = 1.0e-14_rp, tol_eig = 1.0e-8_rp
   integer, parameter :: meshes(4) = [750, 1500, 3000, 6000]

   type(parallel_context) :: parallel
   type(control), target :: ctl
   type(lattice), target :: lat
   type(charge), target :: chg
   type(hamiltonian), target :: ham
   type(energy), target :: en
   type(recursion), target :: rec
   type(dos), target :: dos_state
   type(green), target :: gr
   type(bands) :: rs
   type(reciprocal) :: ks
   complex(rp), allocatable :: q(:, :), h(:, :)
   real(rp), allocatable :: eig(:)
   complex(rp) :: p_exact(10, 10), m1_exact(10, 10), p_lor(10, 10), m1_lor(10, 10)
   complex(rp) :: p_ks(10, 10), m1_ks(10, 10), p_rs(10, 10), m1_rs(10, 10), p_col(2, 5, 5)
   complex(rp) :: p_wrong(10, 10), m1_wrong(10, 10), swapped(10, 10)
   logical :: failed = .false.
   integer :: local_nk, ik

   external :: zheev
#ifdef USE_MPI
   call MPI_INIT(ierr)
#endif
   parallel = parallel_context()
   call g_logger%init()
   call basis_init(2)
   call init_math_operators()
   call get_mpi_variables(rank, 1)
   allocate (q(nb, nb), h(nb, nb), eig(nb))
   call setup_objects()

   ! ---- (a), (c): full model, xi=0.35, spin-mixing coupling on
   call build_model(0.35_rp, 1.0_rp)
   call exact_references()
   call setup_ks()
   call ks%calculate_spinor_dm(lat)
   p_ks = lat%symbolic_atoms(1)%potential%ldm_spinor
   m1_ks = lat%symbolic_atoms(1)%potential%m1_spinor
   call report('(a) KS P vs exact eigenpair sum', maxval(abs(p_ks - p_exact)), tol_ks)
   call report('(a) KS M1 vs exact eigenpair sum [Ry]', maxval(abs(m1_ks - m1_exact)), tol_ks)
   call report('(a) max |P_up,dn| of the model (nonzero)', maxval(abs(p_exact(1:5, 6:10))), 1.0e-3_rp, .true.)

   ks%eigenvectors = conjg(ks%eigenvectors)
   call ks%calculate_spinor_dm(lat)
   p_wrong = lat%symbolic_atoms(1)%potential%ldm_spinor
   m1_wrong = lat%symbolic_atoms(1)%potential%m1_spinor
   ks%eigenvectors = conjg(ks%eigenvectors)
   call report('(a) control, KS conjugated eigenvectors: P', maxval(abs(p_wrong - p_exact)), 0.1_rp, .true.)
   call report('(a) control, KS conjugated eigenvectors: M1 [Ry]', maxval(abs(m1_wrong - m1_exact)), 3.0e-3_rp, .true.)
   call report('(a) control, KS transposed P', maxval(abs(transpose(p_ks) - p_exact)), 0.1_rp, .true.)
   swapped(1:5, 1:5) = p_ks(6:10, 6:10)
   swapped(6:10, 6:10) = p_ks(1:5, 1:5)
   swapped(1:5, 6:10) = p_ks(6:10, 1:5)
   swapped(6:10, 1:5) = p_ks(1:5, 6:10)
   call report('(a) control, KS spin blocks swapped', maxval(abs(swapped - p_exact)), 1.0e-3_rp, .true.)

   call setup_rs_green(.true.)
   call rs%calculate_moments()
   p_rs = lat%symbolic_atoms(1)%potential%ldm_spinor
   m1_rs = lat%symbolic_atoms(1)%potential%m1_spinor
   call report('(a) RS P vs finite-window Lorentzian', maxval(abs(p_rs - p_lor)), tol_rs_p)
   call report('(a) RS M1 vs finite-window Lorentzian [Ry]', maxval(abs(m1_rs - m1_lor)), tol_rs_m1)
   call mesh_refinement()
   call setup_rs_green(.false.)
   call rs%calculate_moments()
   p_wrong = lat%symbolic_atoms(1)%potential%ldm_spinor
   m1_wrong = lat%symbolic_atoms(1)%potential%m1_spinor
   call report('(a) control, RS Green function of H^T: P', maxval(abs(p_wrong - p_exact)), 0.1_rp, .true.)
   call report('(a) control, RS Green function of H^T: M1 [Ry]', maxval(abs(m1_wrong - m1_exact)), 3.0e-3_rp, .true.)
   swapped(1:5, 1:5) = p_rs(6:10, 6:10)
   swapped(6:10, 6:10) = p_rs(1:5, 1:5)
   swapped(1:5, 6:10) = p_rs(6:10, 1:5)
   swapped(6:10, 1:5) = p_rs(1:5, 6:10)
   call report('(a) control, RS spin blocks swapped', maxval(abs(swapped - p_lor)), 1.0e-3_rp, .true.)

   ! ---- (c) Hermiticity, eigenvalues, trace against the radial d charge
   call checks_physical('KS', p_ks)
   call setup_rs_green(.true.)
   call rs%calculate_moments()
   call checks_physical('RS', lat%symbolic_atoms(1)%potential%ldm_spinor)
   call ks%calculate_band_moments()
   call report('(c) KS Tr P_ss vs band-moment q0 (both spins)', &
      max(abs(real(trace_block(p_ks, 1), rp) - ks%band_moments(1, 3, 1, 1)), &
          abs(real(trace_block(p_ks, 2), rp) - ks%band_moments(1, 3, 2, 1))), tol_ks)
   call report('(c) KS Tr P vs summed d charge', &
      abs(real(trace_block(p_ks, 1) + trace_block(p_ks, 2), rp) - ks%band_moments(1, 3, 1, 1) - &
          ks%band_moments(1, 3, 2, 1)), tol_ks)
   call report('(c) RS Tr P_ss vs ql(1,d,s) (same Green function and quadrature, by construction)', &
      max(abs(real(trace_block(lat%symbolic_atoms(1)%potential%ldm_spinor, 1), rp) - &
              lat%symbolic_atoms(1)%potential%ql(1, 2, 1)), &
          abs(real(trace_block(lat%symbolic_atoms(1)%potential%ldm_spinor, 2), rp) - &
              lat%symbolic_atoms(1)%potential%ql(1, 2, 2))), 1.0e-10_rp)

   ! ---- (b) xi -> 0 and no spin-mixing coupling
   call build_model(0.0_rp, 0.0_rp)
   call exact_references()
   call setup_ks()
   call ks%calculate_spinor_dm(lat)
   p_ks = lat%symbolic_atoms(1)%potential%ldm_spinor
   call report('(b) xi=0: KS spin-off-diagonal blocks', max(maxval(abs(p_ks(1:5, 6:10))), maxval(abs(p_ks(6:10, 1:5)))), tol_zero)
   ctl%nsp = 1
   call ks%calculate_ldm_from_projected_dos(lat, fixed_liechtenstein=.true.)
   p_col(1, :, :) = lat%symbolic_atoms(1)%potential%ldm(3, 1, 1:5, 1:5)
   p_col(2, :, :) = lat%symbolic_atoms(1)%potential%ldm(3, 2, 1:5, 1:5)
   ctl%nsp = 2
   call report('(b) xi=0: KS diagonal blocks vs the nsp=1 extraction', &
      max(maxval(abs(p_ks(1:5, 1:5) - p_col(1, :, :))), maxval(abs(p_ks(6:10, 6:10) - p_col(2, :, :)))), tol_zero)
   call setup_rs_green(.true.)
   call rs%calculate_moments()
   p_rs = lat%symbolic_atoms(1)%potential%ldm_spinor
   call report('(b) xi=0: RS spin-off-diagonal blocks', max(maxval(abs(p_rs(1:5, 6:10))), maxval(abs(p_rs(6:10, 1:5)))), tol_zero)
   ctl%nsp = 1
   call rs%calculate_moments()
   ctl%nsp = 2
   call report('(b) xi=0: RS diagonal blocks vs the nsp=1 extraction', &
      max(maxval(abs(p_rs(1:5, 1:5) - lat%symbolic_atoms(1)%potential%ldm(3, 1, 1:5, 1:5))), &
          maxval(abs(p_rs(6:10, 6:10) - lat%symbolic_atoms(1)%potential%ldm(3, 2, 1:5, 1:5)))), tol_zero)
   call build_model(0.35_rp, 0.0_rp)
   call exact_references()
   call setup_ks()
   call ks%calculate_spinor_dm(lat)
   call report('(b) control, xi=0.35 without the coupling: KS spin-off-diagonal blocks', &
      maxval(abs(lat%symbolic_atoms(1)%potential%ldm_spinor(1:5, 6:10))), tol_zero, .true.)

#ifdef USE_MPI
   call reduce_failure()
#endif
   if (rank == 0) then
      if (failed) then
         write (*, '(a)') 'projector_spinor_frozen_h: FAIL'
      else
         write (*, '(a)') 'projector_spinor_frozen_h: PASS'
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

   function trace_block(p, spin) result(t)
      complex(rp), intent(in) :: p(10, 10)
      integer, intent(in) :: spin
      complex(rp) :: t
      integer :: i

      t = (0.0_rp, 0.0_rp)
      do i = 1, 5
         t = t + p((spin - 1)*5 + i, (spin - 1)*5 + i)
      end do
   end function trace_block

   subroutine checks_physical(label, p)
      character(len=*), intent(in) :: label
      complex(rp), intent(in) :: p(10, 10)
      real(rp) :: trace_value, min_eig, max_eig, antiherm
      integer :: info

      call hubbard_density_diagnostics(p, trace_value, min_eig, max_eig, antiherm, info)
      call report('(c) '//label//' P Hermiticity (by construction in RS)', antiherm, tol_herm)
      call report('(c) '//label//' P eigenvalues in [-1e-8, 1+1e-8]', max(-min_eig, max_eig - 1.0_rp, 0.0_rp), tol_eig)
   end subroutine checks_physical

   !> H = d crystal field (complex orbital mixing) + exchange splitting + xi L.S + coupling*X (up-down).
   subroutine build_model(xi, coupling)
      real(rp), intent(in) :: xi, coupling
      complex(rp) :: hd(5, 5), lz(5, 5), lp(5, 5), x(5, 5)
      real(rp) :: lm
      real(rp), parameter :: es(4) = [-1.4_rp, -1.0_rp, -0.8_rp, -0.6_rp]
      real(rp), parameter :: ed(5) = [-0.34_rp, -0.19_rp, 0.16_rp, 0.39_rp, 0.61_rp]
      integer :: i, j, m

      h = (0.0_rp, 0.0_rp)
      do i = 1, 4
         h(i, i) = es(i)
         h(norb + i, norb + i) = es(i) + 0.05_rp
      end do
      hd = (0.0_rp, 0.0_rp)
      lz = (0.0_rp, 0.0_rp)
      lp = (0.0_rp, 0.0_rp)
      x = (0.0_rp, 0.0_rp)
      do i = 1, 5
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
      do m = -2, 1                                       ! L+ |m> = sqrt(l(l+1) - m(m+1)) |m+1>
         lm = sqrt(6.0_rp - real(m*(m + 1), rp))
         lp(m + 4, m + 3) = cmplx(lm, 0.0_rp, rp)
      end do
      h(5:9, 5:9) = hd - 0.07_rp*identity5() + 0.5_rp*xi*lz
      h(norb + 5:nb, norb + 5:nb) = hd + 0.07_rp*identity5() - 0.5_rp*xi*lz
      h(5:9, norb + 5:nb) = 0.5_rp*xi*conjg(transpose(lp)) + coupling*x     ! (up,down): L-/2
      h(norb + 5:nb, 5:9) = 0.5_rp*xi*lp + coupling*conjg(transpose(x))     ! (down,up): L+/2
   end subroutine build_model

   function identity5() result(e)
      complex(rp) :: e(5, 5)
      integer :: i

      e = (0.0_rp, 0.0_rp)
      do i = 1, 5
         e(i, i) = (1.0_rp, 0.0_rp)
      end do
   end function identity5

   subroutine exact_references()
      complex(rp), allocatable :: work(:)
      real(rp), allocatable :: rwork(:)
      real(rp) :: occ, occ_m1, gap_best, shift
      integer :: info, a, b, k, row, col, sa, sb

      q = h
      allocate (work(2*nb), rwork(3*nb - 2))
      call zheev('V', 'U', nb, q, nb, eig, work, 2*nb, rwork, info)
      if (info /= 0) failed = .true.
      ! Put E_F in the middle of the widest gap between the 6th and 12th levels so that no level
      ! sits near the Fermi endpoint (the finite-window quadrature budget assumes this).
      gap_best = 0.0_rp
      do k = 6, 11
         if (eig(k + 1) - eig(k) > gap_best) then
            gap_best = eig(k + 1) - eig(k)
            shift = 0.5_rp*(eig(k) + eig(k + 1))
         end if
      end do
      eig = eig - shift
      if (gap_best < 0.1_rp) failed = .true.
      p_exact = (0.0_rp, 0.0_rp)
      m1_exact = (0.0_rp, 0.0_rp)
      p_lor = (0.0_rp, 0.0_rp)
      m1_lor = (0.0_rp, 0.0_rp)
      do sa = 1, 2
         do sb = 1, 2
            do a = 1, 5
               do b = 1, 5
                  row = (sa - 1)*norb + 4 + a
                  col = (sb - 1)*norb + 4 + b
                  do k = 1, nb
                     if (eig(k) < ef) then
                        p_exact((sa - 1)*5 + a, (sb - 1)*5 + b) = p_exact((sa - 1)*5 + a, (sb - 1)*5 + b) + q(row, k)*conjg(q(col, k))
                        m1_exact((sa - 1)*5 + a, (sb - 1)*5 + b) = m1_exact((sa - 1)*5 + a, (sb - 1)*5 + b) + &
                           eig(k)*q(row, k)*conjg(q(col, k))
                     end if
                     occ = (atan((ef - eig(k))/eta) - atan((emin - eig(k))/eta))/pi
                     occ_m1 = eig(k)*occ + eta/(2.0_rp*pi)*log(((ef - eig(k))**2 + eta**2)/((emin - eig(k))**2 + eta**2))
                     p_lor((sa - 1)*5 + a, (sb - 1)*5 + b) = p_lor((sa - 1)*5 + a, (sb - 1)*5 + b) + occ*q(row, k)*conjg(q(col, k))
                     m1_lor((sa - 1)*5 + a, (sb - 1)*5 + b) = m1_lor((sa - 1)*5 + a, (sb - 1)*5 + b) + &
                        occ_m1*q(row, k)*conjg(q(col, k))
                  end do
               end do
            end do
         end do
      end do
      deallocate (work, rwork)
   end subroutine exact_references

   !> Independent Simpson integration of the Lorentzian spectral function at four meshes.
   subroutine mesh_refinement()
      real(rp), allocatable :: mesh(:), yr(:), yi(:)
      complex(rp) :: spectral, gab, gba, integ(10, 10), integ_m1(10, 10)
      real(rp) :: step, endpoint, vr, vi, mr, mi, errs(4), m1errs(4), egrid
      integer :: level, n, npts, a, b, sa, sb, row, col, band_id, ie

      do level = 1, 4
         n = meshes(level)
         step = (ef - emin)/real(n/2, rp)
         npts = nint((ef - emin)/step)
         endpoint = emin + real(npts - 1, rp)*step
         allocate (mesh(n + 10), yr(n + 10), yi(n + 10))
         do ie = 1, n + 10
            mesh(ie) = emin + real(ie - 1, rp)*step
         end do
         do sa = 1, 2
            do sb = 1, 2
               do a = 1, 5
                  do b = 1, 5
                     row = (sa - 1)*norb + 4 + a
                     col = (sb - 1)*norb + 4 + b
                     do ie = 1, n + 10
                        spectral = (0.0_rp, 0.0_rp)
                        do band_id = 1, nb
                           gab = q(row, band_id)*conjg(q(col, band_id))/cmplx(mesh(ie) - eig(band_id), eta, rp)
                           gba = q(col, band_id)*conjg(q(row, band_id))/cmplx(mesh(ie) - eig(band_id), eta, rp)
                           spectral = spectral + hubbard_spectral_density_element(gab, gba)
                        end do
                        yr(ie) = real(spectral, rp)
                        yi(ie) = aimag(spectral)
                     end do
                     call simpson_m(vr, step, ef, npts, yr, endpoint, 0, mesh)
                     call simpson_m(vi, step, ef, npts, yi, endpoint, 0, mesh)
                     call simpson_m(mr, step, ef, npts, yr, endpoint, 1, mesh)
                     call simpson_m(mi, step, ef, npts, yi, endpoint, 1, mesh)
                     integ((sa - 1)*5 + a, (sb - 1)*5 + b) = cmplx(vr, vi, rp)
                     integ_m1((sa - 1)*5 + a, (sb - 1)*5 + b) = cmplx(mr, mi, rp)
                  end do
               end do
            end do
         end do
         errs(level) = maxval(abs(integ - p_lor))
         m1errs(level) = maxval(abs(integ_m1 - m1_lor))
         deallocate (mesh, yr, yi)
         egrid = 0.0_rp
      end do
      call report('(a) RS Simpson error vs Lorentzian, nchannels=6000', errs(4), tol_rs_p)
      call report('(a) RS M1 Simpson error vs Lorentzian, nchannels=6000 [Ry]', m1errs(4), tol_rs_m1)
      call report('(a) refinement: P error decreases 750->1500->3000->6000 (violations)', &
         real(count(errs(2:4) >= errs(1:3)), rp), 0.5_rp)
      call report('(a) refinement: M1 error decreases 750->1500->3000->6000 (violations)', &
         real(count(m1errs(2:4) >= m1errs(1:3)), rp), 0.5_rp)
      call report('(a) RS production vs the 6000-point Simpson reference', maxval(abs(p_rs - integ)), 5.0e-12_rp)
   end subroutine mesh_refinement

   subroutine setup_objects()
      call ctl%restore_to_default()
      ctl%lmax = 2
      ctl%nsp = 2
      ctl%calctype = 'B'
      ctl%density_policy = 'constrained_spiral'     ! radial moments resolved along the global z axis

      lat%nrec = 1; lat%nbulk = 0; lat%nbulk_bulk = 1; lat%ntype = 1
      lat%nn_max = 1; lat%kk = 1; lat%nmax = 1; lat%alat = 1.0_rp
      lat%r2 = 1.0_rp; lat%a = 0.0_rp; lat%a(1, 1) = 1.0_rp
      lat%a(2, 2) = 1.0_rp; lat%a(3, 3) = 1.0_rp
      lat%a_cart_inv = lat%a; lat%a_cart_inv_ready = .true.
      lat%control => ctl
      allocate (lat%ib(1), lat%atlist(1), lat%iz(1), lat%num(1), lat%nn(1, 1), &
                lat%sbar(norb, norb, 1, 1), lat%symbolic_atoms(1))
      lat%ib = 1; lat%atlist = 1; lat%iz = 1; lat%num = 1; lat%nn = 1
      lat%sbar = cmplx(0.0_rp, 0.0_rp, rp)
      call lat%symbolic_atoms(1)%restore_to_default()
      lat%symbolic_atoms(1)%potential%hubbard_u(3) = 0.25_rp
      lat%symbolic_atoms(1)%potential%hubbard_j(3) = 0.04_rp
      lat%symbolic_atoms(1)%potential%mom = [0.0_rp, 0.0_rp, 1.0_rp]
      lat%symbolic_atoms(1)%potential%qpar = 1.0_rp
      lat%symbolic_atoms(1)%potential%srdel = 0.1_rp
      lat%symbolic_atoms(1)%potential%c = 0.0_rp
      lat%symbolic_atoms(1)%potential%gravity_center = 0.0_rp
      lat%symbolic_atoms(1)%potential%pl = 0.0_rp

      chg%lattice => lat
      chg%symbolic_atom => lat%symbolic_atoms
      ham%charge => chg; ham%lattice => lat; ham%control => ctl
      ham%local_axis = .false.
      ham%operator_generation = 0
      ham%hubbard_u_general_check = .true.
      ham%hubbard_spinor_dm_diagnostic = .true.
      ham%hubbard_u_potential_form = 'liechtenstein'
      ham%hubbard_u_sc_check = .false.; ham%hubbard_v_check = .false.
      allocate (ham%hubbard_u_pot(nb, nb, 1))
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
      allocate (gr%g0(nb, nb, nch + 10, atoms_per_process))
      gr%g0 = cmplx(0.0_rp, 0.0_rp, rp)

      rs%green => gr; rs%lattice => lat; rs%symbolic_atom => lat%symbolic_atoms
      rs%dos => dos_state; rs%en => en; rs%control => ctl; rs%recursion => rec
      call rs%restore_to_default()
      rs%nv1 = nint((en%fermi - en%energy_min)/en%edel)
      rs%e1 = en%energy_min + real(rs%nv1 - 1, rp)*en%edel
   end subroutine setup_objects

   subroutine setup_ks()
      real(rp) :: all_points(3, nk), weights(nk)
      type(kpoint_workset) :: workset

      all_points = 0.0_rp
      weights = 1.0_rp/real(nk, rp)
      workset = make_kpoint_workset(all_points, weights, parallel, numprocs > 1, .true.)
      local_nk = workset%nk_local
      ks%hamiltonian => ham; ks%lattice => lat; ks%control => ctl
      ks%n_spin_components = 2
      ks%n_sites = 1; ks%n_orb_types = 3
      ks%nk_mesh = [nk, 1, 1]
      ks%nk_total = nk; ks%nk_local = local_nk
      ks%k_mesh_distributed_active = numprocs > 1
      ks%use_symmetry_reduction = .false.
      ks%temperature = 0.0_rp; ks%fermi_level = ef
      ks%k_workset = workset
      if (allocated(ks%eigenvalues)) deallocate (ks%eigenvalues, ks%eigenvectors, ks%k_weights, ks%k_l2g_map, ks%k_g2l_map)
      allocate (ks%eigenvalues(nb, local_nk), ks%eigenvectors(nb, nb, local_nk))
      allocate (ks%k_weights(nk), ks%k_l2g_map(local_nk), ks%k_g2l_map(nk))
      ks%k_weights = weights
      ks%k_points = all_points
      ks%k_l2g_map = workset%local_to_global
      ks%k_g2l_map = workset%global_to_local
      ks%cached_operator_generation = ham%operator_generation
      do ik = 1, local_nk
         ks%eigenvalues(:, ik) = eig
         ks%eigenvectors(:, :, ik) = q
      end do
   end subroutine setup_ks

   subroutine setup_rs_green(correct)
      logical, intent(in) :: correct
      complex(rp), allocatable :: gmat(:, :)
      integer :: ie, band, a, b
      complex(rp) :: z

      gr%g0 = cmplx(0.0_rp, 0.0_rp, rp)
      allocate (gmat(nb, nb))
      do ie = 1, nch + 10
         gmat = cmplx(0.0_rp, 0.0_rp, rp)
         do band = 1, nb
            z = cmplx(en%ene(ie) - eig(band), eta, rp)
            do a = 1, nb
               do b = 1, nb
                  if (correct) then
                     gmat(a, b) = gmat(a, b) + q(a, band)*conjg(q(b, band))/z
                  else
                     gmat(a, b) = gmat(a, b) + conjg(q(a, band))*q(b, band)/z
                  end if
               end do
            end do
         end do
         gr%g0(:, :, ie, 1) = gmat
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

end program projector_spinor_frozen_h
