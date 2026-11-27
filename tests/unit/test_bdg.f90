!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!> @brief BdG kernels (bdg_assemble, bdg_accumulate, bdg_gap, bdg_iterate) on the developer oracles O1-O4.
!> @details Argument after the case name: oracle directory. Model Hamiltonians are diagonal in the index
!>          L + (sigma-1)*norb, built from band arrays eps(nk, nmat) on the mesh k_n = (n-1/2)/nk (w = 1/nk).
!>          Tolerances come from <dir>/tolerances.nml. Real parts are compared signed, relative error |x-ref|/|ref|.
!------------------------------------------------------------------------------
module test_bdg_mod
   use precision_mod, only: rp
   use bdg_mod
   use math_mod, only: hermitian_eigensystem
   implicit none

   integer, parameter :: nmax = 4096, max_iter_budget = 3000
   real(rp) :: delta_rel, delta_near_tc_rel, chi_rel, delta_imag_abs, nel_abs, moment_abs, normal_abs, spectrum_abs
   logical :: failed = .false.

   !> Model: H(k_ik) = diag(eps(ik, :)) + h0, H(-k_ik) = diag(eps(nk+1-ik, :)) + h0 (h0 optional, k-independent).
   type, extends(bdg_hamiltonian_source) :: model_source
      real(rp), allocatable :: eps(:, :)
      complex(rp), allocatable :: h0(:, :) !< optional k-independent matrix added to H(k) and H(-k)
   contains
      procedure :: hk => model_hk
   end type model_source

contains

   subroutine model_hk(this, ik, hk, hmk)
      class(model_source), intent(in) :: this
      integer, intent(in) :: ik
      complex(rp), intent(out) :: hk(:, :), hmk(:, :)
      integer :: a

      hk = (0.0_rp, 0.0_rp)
      hmk = (0.0_rp, 0.0_rp)
      do a = 1, size(hk, 1)
         hk(a, a) = this%eps(ik, a)
         hmk(a, a) = this%eps(size(this%eps, 1) + 1 - ik, a)
      end do
      if (allocated(this%h0)) then
         hk = hk + this%h0
         hmk = hmk + this%h0
      end if
   end subroutine model_hk

   subroutine read_tolerances(dir)
      character(len=*), intent(in) :: dir
      integer :: iu
      namelist /bdg_tolerances/ delta_rel, delta_near_tc_rel, chi_rel, delta_imag_abs, nel_abs, moment_abs, normal_abs, &
         spectrum_abs
      open (newunit=iu, file=dir//'/tolerances.nml', status='old', action='read')
      read (iu, nml=bdg_tolerances)
      close (iu)
      write (*, '(a)') 'test | measured | tolerance'
   end subroutine read_tolerances

   subroutine report(name, err, tol)
      character(len=*), intent(in) :: name
      real(rp), intent(in) :: err, tol
      write (*, '(a,a,es10.2,a,es10.2,a)') name, ' | ', err, ' | ', tol, merge('      ', ' FAIL ', err < tol)
      if (.not. (err < tol)) failed = .true.
   end subroutine report

   !> Convergence row: the last iteration's max|delta_new - delta_old| must be below the iteration tolerance.
   subroutine report_conv(name, niter, dmax, tol)
      character(len=*), intent(in) :: name
      integer, intent(in) :: niter
      real(rp), intent(in) :: dmax, tol
      write (*, '(a,a,i0,a,i0)') name, ' niter = ', niter, ' of budget ', max_iter_budget
      call report(name//' dmax', dmax, tol)
   end subroutine report_conv

   real(rp) function rel(x, ref)
      real(rp), intent(in) :: x, ref
      rel = abs(x - ref)/abs(ref)
   end function rel

   subroutine finish()
      if (failed) error stop 1
   end subroutine finish

   !> One site, norb orbitals; runs bdg_iterate on bands eps(nk, 2 norb). delta holds the seed on input.
   subroutine solve(eps, norb, mu, kt, kmat, mix, tol, max_iter, delta, chi, nel, niter, dmax, h0)
      real(rp), intent(in) :: eps(:, :), mu, kt, kmat(:, :), mix, tol
      complex(rp), intent(in), optional :: h0(:, :)
      integer, intent(in) :: norb, max_iter
      complex(rp), intent(inout) :: delta(:, :)
      complex(rp), intent(out) :: chi(:, :)
      real(rp), intent(out) :: nel(2), dmax
      integer, intent(out) :: niter
      type(model_source) :: src
      integer :: a, site_of(2*norb), orb_of(2*norb), spin_of(2*norb), nk

      nk = size(eps, 1)
      src%eps = eps
      if (present(h0)) src%h0 = h0
      site_of = 1
      do a = 1, 2*norb
         orb_of(a) = mod(a - 1, norb) + 1
         spin_of(a) = merge(1, 2, a <= norb)
      end do
      call bdg_iterate(src, spread(1.0_rp/nk, 1, nk), mu, kt, reshape(kmat, [norb, norb, 1]), [1], site_of, orb_of, &
                       spin_of, mix, tol, max_iter, delta, chi, nel, niter, dmax)
   end subroutine solve

   !> Reads the O1 band e1(nk), mu and g from the o1 file.
   subroutine read_o1_band(dir, e1, mu, g)
      character(len=*), intent(in) :: dir
      real(rp), allocatable, intent(out) :: e1(:)
      real(rp), intent(out) :: mu, g
      character(len=512) :: mesh_rule, basis_rule, band_rule
      integer :: nk, norb, iu
      real(rp) :: t, k(nmax), eps(nmax)
      namelist /o1_inputs/ nk, mu, mesh_rule, basis_rule, norb, t, g, band_rule, k, eps
      open (newunit=iu, file=dir//'/o1_single_band.nml', status='old', action='read')
      read (iu, nml=o1_inputs)
      close (iu)
      e1 = eps(1:nk)
   end subroutine read_o1_band

   !> Catches: overall sign of D, f in place of (1 - f) (dmax row), missing (1 - f), wrong k weights, wrong N formula.
   !> Not caught here: the sign of D(dn,up) alone (D = Delta x 1), which spin-diagonal H hides; see test_bdg_o2x.
   subroutine test_bdg_o1(dir)
      character(len=*), intent(in) :: dir
      character(len=512) :: mesh_rule, basis_rule, band_rule
      integer :: nk, norb, iu, it, niter
      real(rp) :: mu, t, g, k(nmax), eps(nmax), nel(2), dmax, mix(3), tol(3), kt(4)
      real(rp) :: delta_ref(4), chi_ref(3), nel_ref(3), nnor_ref(3)
      real(rp) :: ktc, kt_t0, delta_t0, chi_t0, nel_t0, nel_normal_t0, kt_t05, delta_t05, chi_t05, nel_t05, nel_normal_t05
      real(rp) :: kt_t09, delta_t09, chi_t09, nel_t09, nel_normal_t09, kt_t11, delta_t11, lin_rate_t11
      complex(rp) :: delta(1, 1), chi(1, 1)
      character(len=4), parameter :: lab(3) = ['t0 ', 't05', 't09']
      namelist /o1_inputs/ nk, mu, mesh_rule, basis_rule, norb, t, g, band_rule, k, eps
      namelist /o1_reference/ ktc, kt_t0, delta_t0, chi_t0, nel_t0, nel_normal_t0, kt_t05, delta_t05, chi_t05, nel_t05, &
         nel_normal_t05, kt_t09, delta_t09, chi_t09, nel_t09, nel_normal_t09, kt_t11, delta_t11, lin_rate_t11

      call read_tolerances(dir)
      open (newunit=iu, file=dir//'/o1_single_band.nml', status='old', action='read')
      read (iu, nml=o1_inputs)
      read (iu, nml=o1_reference)
      close (iu)
      kt = [kt_t0, kt_t05, kt_t09, kt_t11]
      delta_ref = [delta_t0, delta_t05, delta_t09, delta_t11]
      chi_ref = [chi_t0, chi_t05, chi_t09]
      nel_ref = [nel_t0, nel_t05, nel_t09]
      nnor_ref = [nel_normal_t0, nel_normal_t05, nel_normal_t09]
      mix = [1.0_rp, 1.0_rp, 4.0_rp]
      tol = [1.0e-15_rp, 1.0e-15_rp, 1.0e-15_rp]
      do it = 1, 3
         delta = 0.003_rp
         call solve(reshape(spread(eps(1:nk), 2, 2), [nk, 2]), 1, mu, kt(it), reshape([g], [1, 1]), mix(it), tol(it), max_iter_budget, &
                    delta, chi, nel, niter, dmax)
         write (*, '(a,a,i0,a,es9.2,a,f4.1)') lab(it), ' niter = ', niter, ' tol = ', tol(it), ' mix = ', mix(it)
         call report_conv('O1 '//trim(lab(it)), niter, dmax, tol(it))
         call report('O1 delta '//trim(lab(it)), rel(real(delta(1, 1), rp), delta_ref(it)), &
                     merge(delta_near_tc_rel, delta_rel, it == 3))
         call report('O1 chi '//trim(lab(it)), rel(real(chi(1, 1), rp), chi_ref(it)), chi_rel)
         call report('O1 sum nel '//trim(lab(it)), abs(sum(nel) - nel_ref(it)), nel_abs)
         call report('O1 |Im delta| '//trim(lab(it)), abs(aimag(delta(1, 1))), delta_imag_abs)
         delta = 0.0_rp
         call solve(reshape(spread(eps(1:nk), 2, 2), [nk, 2]), 1, mu, kt(it), reshape([g], [1, 1]), 1.0_rp, 1.0_rp, 0, &
                    delta, chi, nel, niter, dmax)
         call report('O1 normal nel '//trim(lab(it)), abs(sum(nel) - nnor_ref(it)), nel_abs)
      end do
      delta = 0.003_rp
      call solve(reshape(spread(eps(1:nk), 2, 2), [nk, 2]), 1, mu, kt(4), reshape([g], [1, 1]), 10.0_rp, 1.0e-14_rp, max_iter_budget, &
                 delta, chi, nel, niter, dmax)
      write (*, '(a,i0,a)') 't11 niter = ', niter, ' tol = 1.0e-14 mix = 10.0'
      call report_conv('O1 t11', niter, dmax, 1.0e-14_rp)
      call report('O1 |delta| t11', abs(delta(1, 1)), normal_abs)
      call finish()
   end subroutine test_bdg_o1

   !> Catches: K applied as a diagonal g/3 * I instead of the full matrix (then Delta = 0).
   subroutine test_bdg_o1b(dir)
      character(len=*), intent(in) :: dir
      character(len=512) :: mesh_rule, basis_rule, band_rule
      integer :: nk, norb, iu, niter, l
      real(rp) :: mu, g, kmat(nmax), kt, nel(2), dmax, delta(nmax), chi(nmax)
      real(rp), allocatable :: e1(:)
      complex(rp) :: dl(3, 1), ch(3, 1)
      namelist /o1b_inputs/ nk, mu, mesh_rule, basis_rule, norb, kmat, kt, band_rule
      namelist /o1b_reference/ delta, chi

      call read_tolerances(dir)
      call read_o1_band(dir, e1, mu, g)
      open (newunit=iu, file=dir//'/o1b_uniform_k.nml', status='old', action='read')
      read (iu, nml=o1b_inputs)
      read (iu, nml=o1b_reference)
      close (iu)
      dl(:, 1) = [0.003_rp, 0.004_rp, 0.005_rp]
      call solve(spread(e1, 2, 6), 3, mu, kt, reshape(kmat(1:9), [3, 3]), 1.0_rp, 1.0e-15_rp, max_iter_budget, dl, ch, nel, niter, dmax)
      write (*, '(a,i0,a)') 'o1b niter = ', niter, ' tol = 1.0e-15 mix = 1.0'
      call report_conv('O1b', niter, dmax, 1.0e-15_rp)
      do l = 1, 3
         call report('O1b delta L='//char(48 + l), rel(real(dl(l, 1), rp), delta(l)), delta_rel)
         call report('O1b chi L='//char(48 + l), rel(real(ch(l, 1), rp), chi(l)), chi_rel)
         call report('O1b |Im delta| L='//char(48 + l), abs(aimag(dl(l, 1))), delta_imag_abs)
      end do
      call finish()
   end subroutine test_bdg_o1b

   !> Catches: Zeeman field acting as a pair breaker at kT = 0, or dropped from the pair energies, or occupations mixed up.
   subroutine test_bdg_o2(dir)
      character(len=*), intent(in) :: dir
      character(len=512) :: mesh_rule, basis_rule, band_rule
      integer :: nk, norb, iu, niter, ic
      real(rp) :: mu, g, nel(2), dmax, h(2), kt(2), dref(2), cref(2), nu(2), nd(2)
      real(rp) :: kt_a, h_a, delta_a, chi_a, n_up_a, n_dn_a, kt_b, h_b, delta_b, chi_b, n_up_b, n_dn_b
      real(rp), allocatable :: e1(:), eps2(:, :)
      complex(rp) :: delta(1, 1), chi(1, 1)
      character(len=1), parameter :: lab(2) = ['a', 'b']
      namelist /o2_inputs/ nk, mu, mesh_rule, basis_rule, norb, g, band_rule
      namelist /o2_reference/ kt_a, h_a, delta_a, chi_a, n_up_a, n_dn_a, kt_b, h_b, delta_b, chi_b, n_up_b, n_dn_b

      call read_tolerances(dir)
      call read_o1_band(dir, e1, mu, g)
      open (newunit=iu, file=dir//'/o2_zeeman.nml', status='old', action='read')
      read (iu, nml=o2_inputs)
      read (iu, nml=o2_reference)
      close (iu)
      kt = [kt_a, kt_b]
      h = [h_a, h_b]
      dref = [delta_a, delta_b]
      cref = [chi_a, chi_b]
      nu = [n_up_a, n_up_b]
      nd = [n_dn_a, n_dn_b]
      allocate (eps2(nk, 2))
      do ic = 1, 2
         eps2(:, 1) = e1 - h(ic)
         eps2(:, 2) = e1 + h(ic)
         delta = 0.003_rp
         call solve(eps2, 1, mu, kt(ic), reshape([g], [1, 1]), 1.0_rp, 1.0e-15_rp, max_iter_budget, delta, chi, nel, niter, dmax)
         write (*, '(a,a,a,i0,a)') 'o2', lab(ic), ' ', niter, ' iterations, tol = 1.0e-15 mix = 1.0'
         call report_conv('O2'//lab(ic), niter, dmax, 1.0e-15_rp)
         call report('O2'//lab(ic)//' delta', rel(real(delta(1, 1), rp), dref(ic)), delta_rel)
         call report('O2'//lab(ic)//' chi', rel(real(chi(1, 1), rp), cref(ic)), chi_rel)
         call report('O2'//lab(ic)//' n_up', abs(nel(1) - nu(ic)), nel_abs)
         call report('O2'//lab(ic)//' n_dn', abs(nel(2) - nd(ic)), nel_abs)
         call report('O2'//lab(ic)//' |Im delta|', abs(aimag(delta(1, 1))), delta_imag_abs)
         if (ic == 1) call report('O2a |n_up - n_dn|', abs(nel(1) - nel(2)), moment_abs)
      end do
      call finish()
   end subroutine test_bdg_o2

   !> O2 with the field along x: H = eps(k) I - h sigma_x. On-site same-orbital pairing is a spin singlet, invariant under a
   !> global spin rotation, so Delta, chi, N_up + N_dn equal the O2 (field along z) references and nel(1) = nel(2).
   !> Catches: the sign of D(dn,up) (D = Delta x 1 instead of Delta x i sigma_y), invisible when H is spin-diagonal.
   subroutine test_bdg_o2x(dir)
      character(len=*), intent(in) :: dir
      character(len=512) :: mesh_rule, basis_rule, band_rule
      integer :: nk, norb, iu, niter, ic
      real(rp) :: mu, g, nel(2), dmax, h(2), kt(2), dref(2), cref(2), nu(2), nd(2)
      real(rp) :: kt_a, h_a, delta_a, chi_a, n_up_a, n_dn_a, kt_b, h_b, delta_b, chi_b, n_up_b, n_dn_b
      real(rp), allocatable :: e1(:), eps2(:, :)
      complex(rp) :: delta(1, 1), chi(1, 1), h0(2, 2)
      character(len=1), parameter :: lab(2) = ['a', 'b']
      namelist /o2_inputs/ nk, mu, mesh_rule, basis_rule, norb, g, band_rule
      namelist /o2_reference/ kt_a, h_a, delta_a, chi_a, n_up_a, n_dn_a, kt_b, h_b, delta_b, chi_b, n_up_b, n_dn_b

      call read_tolerances(dir)
      call read_o1_band(dir, e1, mu, g)
      open (newunit=iu, file=dir//'/o2_zeeman.nml', status='old', action='read')
      read (iu, nml=o2_inputs)
      read (iu, nml=o2_reference)
      close (iu)
      kt = [kt_a, kt_b]
      h = [h_a, h_b]
      dref = [delta_a, delta_b]
      cref = [chi_a, chi_b]
      nu = [n_up_a, n_up_b]
      nd = [n_dn_a, n_dn_b]
      allocate (eps2(nk, 2))
      eps2(:, 1) = e1
      eps2(:, 2) = e1
      do ic = 1, 2
         h0 = (0.0_rp, 0.0_rp)
         h0(1, 2) = -h(ic)
         h0(2, 1) = -h(ic)
         delta = 0.003_rp
         call solve(eps2, 1, mu, kt(ic), reshape([g], [1, 1]), 1.0_rp, 1.0e-15_rp, max_iter_budget, delta, chi, nel, niter, dmax, h0)
         write (*, '(a,a,a,i0,a)') 'o2x', lab(ic), ' ', niter, ' iterations, tol = 1.0e-15 mix = 1.0'
         call report_conv('O2x'//lab(ic), niter, dmax, 1.0e-15_rp)
         call report('O2x'//lab(ic)//' delta', rel(real(delta(1, 1), rp), dref(ic)), delta_rel)
         call report('O2x'//lab(ic)//' chi', rel(real(chi(1, 1), rp), cref(ic)), chi_rel)
         call report('O2x'//lab(ic)//' N_up + N_dn', abs(nel(1) + nel(2) - nu(ic) - nd(ic)), nel_abs)
         call report('O2x'//lab(ic)//' nel(1)', abs(nel(1) - 0.5_rp*(nu(ic) + nd(ic))), nel_abs)
         call report('O2x'//lab(ic)//' nel(2)', abs(nel(2) - 0.5_rp*(nu(ic) + nd(ic))), nel_abs)
         call report('O2x'//lab(ic)//' |Im delta|', abs(aimag(delta(1, 1))), delta_imag_abs)
      end do
      call finish()
   end subroutine test_bdg_o2x

   !> Catches: hole block built from -H(k)* instead of -H(-k)* (spectra off by 4e-3 Ry, Delta collapses to 0).
   subroutine test_bdg_o3(dir)
      character(len=*), intent(in) :: dir
      character(len=512) :: mesh_rule, basis_rule, band_rule, spectra_dims, wrong_note
      integer :: nk, norb, iu, niter, j, nidx, ik
      real(rp) :: mu, g, kt, alpha, delta_fixed, spectrum_k_index(16), delta, chi, n_up, n_dn
      real(rp) :: spectra(64), wrong_spectra(64), wrong_delta, nel(2), dmax, ev(4), eps_up(nmax), eps_dn(nmax)
      complex(rp) :: dl(1, 1), ch(1, 1), hk(2, 2), hmk(2, 2), hbdg(4, 4)
      type(model_source) :: src
      namelist /o3_inputs/ nk, mu, mesh_rule, basis_rule, norb, g, kt, alpha, band_rule, eps_up, eps_dn, delta_fixed, &
         spectrum_k_index
      namelist /o3_reference/ delta, chi, n_up, n_dn, spectra, spectra_dims, wrong_spectra, wrong_delta, wrong_note

      call read_tolerances(dir)
      spectrum_k_index = 0.0_rp
      open (newunit=iu, file=dir//'/o3_kodd.nml', status='old', action='read')
      read (iu, nml=o3_inputs)
      read (iu, nml=o3_reference)
      close (iu)
      allocate (src%eps(nk, 2))
      src%eps(:, 1) = eps_up(1:nk)
      src%eps(:, 2) = eps_dn(1:nk)
      dl = 0.003_rp
      call solve(src%eps, 1, mu, kt, reshape([g], [1, 1]), 1.0_rp, 1.0e-15_rp, max_iter_budget, dl, ch, nel, niter, dmax)
      write (*, '(a,i0,a)') 'o3 niter = ', niter, ' tol = 1.0e-15 mix = 1.0'
      call report_conv('O3', niter, dmax, 1.0e-15_rp)
      call report('O3 delta', rel(real(dl(1, 1), rp), delta), delta_rel)
      call report('O3 chi', rel(real(ch(1, 1), rp), chi), chi_rel)
      call report('O3 n_up', abs(nel(1) - n_up), nel_abs)
      call report('O3 n_dn', abs(nel(2) - n_dn), nel_abs)
      call report('O3 |Im delta|', abs(aimag(dl(1, 1))), delta_imag_abs)
      nidx = count(spectrum_k_index > 0.5_rp)
      do j = 1, nidx
         ik = nint(spectrum_k_index(j))
         call src%hk(ik, hk, hmk)
         call bdg_assemble(hk, hmk, mu, reshape([cmplx(delta_fixed, 0.0_rp, rp)], [1, 1]), [1, 1], [1, 1], [1, 2], hbdg)
         call hermitian_eigensystem(hbdg, ev)
         call report('O3 spectrum k#'//char(48 + j), maxval(abs(ev - spectra(4*j - 3:4*j))), spectrum_abs)
      end do
      call finish()
   end subroutine test_bdg_o3

   !> Catches: cross-channel coupling dropped (O4a), K transposed (O4b).
   subroutine test_bdg_o4(dir)
      character(len=*), intent(in) :: dir
      character(len=512) :: mesh_rule, basis_rule, band_rule, kmat_note
      integer :: nk, norb, iu, niter, ic, l
      real(rp) :: mu, kt, t_b, c_b, eps_a(nmax), eps_b(nmax), kmat_a(4), kmat_b(4), nel(2), dmax
      real(rp) :: delta_a(2), chi_a(2), delta_b(2), chi_b(2), kmat(4, 2), dref(2, 2), cref(2, 2)
      real(rp), allocatable :: eps4(:, :)
      complex(rp) :: dl(2, 1), ch(2, 1)
      character(len=1), parameter :: lab(2) = ['a', 'b']
      namelist /o4_inputs/ nk, mu, mesh_rule, basis_rule, norb, kt, t_b, c_b, band_rule, eps_a, eps_b, kmat_a, kmat_b, kmat_note
      namelist /o4_reference/ delta_a, chi_a, delta_b, chi_b

      call read_tolerances(dir)
      open (newunit=iu, file=dir//'/o4_two_band.nml', status='old', action='read')
      read (iu, nml=o4_inputs)
      read (iu, nml=o4_reference)
      close (iu)
      kmat(:, 1) = kmat_a
      kmat(:, 2) = kmat_b
      dref(:, 1) = delta_a
      dref(:, 2) = delta_b
      cref(:, 1) = chi_a
      cref(:, 2) = chi_b
      allocate (eps4(nk, 4))
      eps4(:, 1) = eps_a(1:nk)
      eps4(:, 2) = eps_b(1:nk)
      eps4(:, 3:4) = eps4(:, 1:2)
      do ic = 1, 2
         dl(:, 1) = [0.005_rp, 0.001_rp]
         call solve(eps4, 2, mu, kt, reshape(kmat(:, ic), [2, 2]), 1.0_rp, 1.0e-15_rp, max_iter_budget, dl, ch, nel, niter, dmax)
         write (*, '(a,a,a,i0,a)') 'o4', lab(ic), ' ', niter, ' iterations, tol = 1.0e-15 mix = 1.0'
         call report_conv('O4'//lab(ic), niter, dmax, 1.0e-15_rp)
         do l = 1, 2
            call report('O4'//lab(ic)//' delta L='//char(48 + l), rel(real(dl(l, 1), rp), dref(l, ic)), delta_rel)
            call report('O4'//lab(ic)//' chi L='//char(48 + l), rel(real(ch(l, 1), rp), cref(l, ic)), chi_rel)
            call report('O4'//lab(ic)//' |Im delta| L='//char(48 + l), abs(aimag(dl(l, 1))), delta_imag_abs)
         end do
      end do
      call finish()
   end subroutine test_bdg_o4

   !> Internal plumbing check, not validation. Catches: wrong L->l map, transposed l/l', missing eV->Ry.
   subroutine test_bdg_kmat()
      integer, parameter :: norb = 9, ntype = 2
      integer, parameter :: l_of(norb) = [0, 1, 1, 1, 2, 2, 2, 2, 2]
      real(rp) :: lambda(3, 3, ntype), kmat(norb, norb, ntype), err, ref
      integer :: i, j, t, l, lp

      do t = 1, ntype
         do lp = 0, 2
            do l = 0, 2
               lambda(l + 1, lp + 1, t) = 10.0_rp*l + lp + 100.0_rp*t
            end do
         end do
      end do
      call bdg_kmat_from_lambda(lambda, kmat)
      err = 0.0_rp
      do t = 1, ntype
         do j = 1, norb
            do i = 1, norb
               ref = (10.0_rp*l_of(i) + l_of(j) + 100.0_rp*t)/13.605703976_rp
               err = max(err, abs(kmat(i, j, t) - ref)/ref)
            end do
         end do
      end do
      write (*, '(a)') 'test | measured | tolerance'
      write (*, '(a,es10.2,a)') 'kmat max relative error | ', err, ' | 1.00E-14'
      if (err > 1.0e-14_rp) error stop 'kmat mismatch'
   end subroutine test_bdg_kmat

end module test_bdg_mod
