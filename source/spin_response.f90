!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: spin_response
!
!> @author
!> Anders Bergman
!
! DESCRIPTION:
!> Transverse spin response: pure kernels on plain arrays (C0, d amplitudes) and
!> the spin_response type that hands the frozen ground state to them (C1, C2).
!> Conventions are those of cleanup/spin_response_restart_spec.md:
!>
!>   chi0_ij(w) = sum_k w_k sum_{n,m} (f_nk - f_m,k+q) T_i T_j^* / (w + e_nk - e_m,k+q + i eta)
!>   T_i^{nm}   = sum_mu conj(a^{nk}_{i mu up}) a^{m,k+q}_{i mu dn}
!>   (I + chi0 U) chi = chi0,   chi0(0,0) U M = -M,   L = -(chi - chi^H)/(2 pi i)
!>
!> Spin index of the amplitudes: 1 = up (majority), 2 = down. The kernels take plain
!> arrays; only the spin_response type uses the reciprocal and lattice objects.
!------------------------------------------------------------------------------
module spin_response_mod
   use precision_mod, only: rp
   use math_mod, only: pi, i_unit
   use basis_mod, only: nb, spin_off
   use lattice_mod, only: lattice
   use reciprocal_mod, only: reciprocal, fermi_dirac_occupation, kB_Ry_per_K
   use kpoint_workset_mod, only: kpoint_workset
   use logger_mod, only: g_logger
   use string_mod, only: int2str, real2str
   implicit none
   private

   public :: spin_response
   public :: d_amplitudes_coefficient
   public :: d_amplitudes_juelich
   public :: accumulate_d_occupation
   public :: accumulate_chi0
   public :: u_juelich
   public :: u_mills
   public :: solve_dyson
   public :: spectral_trace
   public :: pole_estimates

   !> Energy difference (Ry) below which the static value uses the Fermi divided difference.
   real(rp), parameter :: degenerate_tol = 1.0e-10_rp

   !> Angular momentum of the response orbitals; d orbitals are ml = l*l + 1 .. (l + 1)**2 of each spin block.
   integer, parameter :: l_d = 2
   !> k points per eigenpair chunk.
   integer, parameter :: k_chunk = 64
   !> Largest accepted |N - target| (electrons) on the response mesh; the electron_count_abs of tests/spin_response/oracles/tolerances.nml.
   real(rp), parameter :: electron_count_tol = 1.0e-6_rp

   type :: spin_response
      class(reciprocal), pointer :: reciprocal => null()
      class(lattice), pointer :: lattice => null()
      !> 'juelich' | 'mills'
      character(len=16) :: method
      character(len=256) :: output_prefix
      real(rp) :: fermi_input, fermi_used, kt, electron_count
      !> Symmetry reduction or time reversal was requested in &reciprocal and overridden by prepare.
      logical :: reduction_overridden
      !> Spin labels were exchanged because the d moment was negative.
      logical :: spin_swapped
      !> Physical spin block of each label: [1, 2], or [2, 1] when spin_swapped.
      integer :: spin_map(2)
      !> d occupation per (site, spin label) with coefficient and Juelich amplitudes.
      real(rp), allocatable :: d_occ_mills(:, :), d_occ_juelich(:, :)
      !> Juelich d moment per site (mu_B), the one that enters the Goldstone condition.
      real(rp), allocatable :: moment(:)
   contains
      procedure :: build_from_file
      procedure :: restore_to_default
      procedure :: prepare
      procedure :: eigenpairs_chunk
      procedure :: write_state
      final :: destructor
   end type spin_response

   interface spin_response
      procedure :: constructor
   end interface spin_response

contains

   !> @brief Constructor: collaborators plus defaults and the &spin_response namelist.
   function constructor(reciprocal_obj) result(obj)
      type(spin_response) :: obj
      type(reciprocal), target, intent(inout) :: reciprocal_obj

      obj%reciprocal => reciprocal_obj
      obj%lattice => reciprocal_obj%lattice
      call obj%restore_to_default()
      call obj%build_from_file()
   end function constructor

   subroutine destructor(this)
      type(spin_response) :: this

      if (allocated(this%d_occ_mills)) deallocate (this%d_occ_mills)
      if (allocated(this%d_occ_juelich)) deallocate (this%d_occ_juelich)
      if (allocated(this%moment)) deallocate (this%moment)
   end subroutine destructor

   subroutine restore_to_default(this)
      class(spin_response), intent(inout) :: this

      this%method = 'juelich'
      this%output_prefix = 'spin_response'
      this%fermi_input = 0.0_rp
      this%fermi_used = 0.0_rp
      this%kt = 0.0_rp
      this%electron_count = 0.0_rp
      this%reduction_overridden = .false.
      this%spin_swapped = .false.
      this%spin_map = [1, 2]
   end subroutine restore_to_default

   subroutine build_from_file(this)
      class(spin_response), intent(inout) :: this

      integer :: iostatus, funit

      include 'include_codes/namelists/spin_response.f90'

      method = this%method
      output_prefix = this%output_prefix

      open (newunit=funit, file=this%lattice%control%fname, action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) then
         call g_logger%fatal('file '//trim(this%lattice%control%fname)//' not found', __FILE__, __LINE__)
      end if
      read (funit, nml=spin_response, iostat=iostatus)
      if (iostatus /= 0 .and. .not. is_iostat_end(iostatus)) then
         call g_logger%fatal('error while reading namelist spin_response, iostat = '//int2str(iostatus), __FILE__, __LINE__)
      end if
      close (funit)

      if (trim(method) /= 'juelich' .and. trim(method) /= 'mills') then
         call g_logger%fatal('spin_response method must be juelich or mills, got '//trim(method), __FILE__, __LINE__)
      end if
      this%method = method
      this%output_prefix = output_prefix
   end subroutine build_from_file

   !> @brief Eigenpairs at k+q for mesh points first:last, from the reciprocal arbitrary-k service.
   !> @param[in]  q_direct  q in fractional coordinates of b1, b2, b3.
   !> @param[out] e         Eigenvalues (Ry), shape (nmat, last-first+1).
   !> @param[out] v         Eigenvectors, shape (nmat, nmat, last-first+1).
   subroutine eigenpairs_chunk(this, q_direct, first, last, e, v)
      class(spin_response), intent(inout) :: this
      real(rp), intent(in) :: q_direct(3)
      integer, intent(in) :: first, last
      real(rp), allocatable, intent(out) :: e(:, :)
      complex(rp), allocatable, intent(out) :: v(:, :, :)

      type(kpoint_workset) :: kq

      kq = this%reciprocal%k_workset%shifted(q_direct)
      call this%reciprocal%calculate_eigenpairs_at_kpoints(kq%points(:, first:last), e, v)
   end subroutine eigenpairs_chunk

   !> @brief Frozen-state handoff (C1) and d occupations and moments (C2).
   !> @details Pass 1 diagonalizes the full-BZ mesh in chunks, keeps the eigenvalues and fixes E_F
   !>          (the &energy value in reciprocal%fermi_level, or the solve on this mesh when
   !>          auto_find_fermi). The electron count is summed from the occupation array that pass 2
   !>          feeds to the d occupations. M_mills < 0 sets spin_swapped and relabels spins.
   subroutine prepare(this)
      class(spin_response), intent(inout) :: this

      class(reciprocal), pointer :: rec
      integer :: nk, nmat, nsite, first, last, ik, n, i, s
      real(rp), allocatable :: f(:, :), e(:, :), e_kn(:, :), enu(:, :), p(:, :), n_mills(:, :), n_juelich(:, :)
      complex(rp), allocatable :: v(:, :, :), a(:, :, :, :, :), aj(:, :, :, :, :)

      rec => this%reciprocal
      nsite = rec%lattice%nrec
      nmat = nb*nsite

      this%reduction_overridden = rec%use_symmetry_reduction .or. rec%use_time_reversal
      rec%use_symmetry_reduction = .false.
      rec%use_time_reversal = .false.
      call rec%generate_mp_mesh()
      nk = rec%k_workset%nk_local
      call g_logger%info('spin_response%prepare: full-BZ mesh '//int2str(rec%nk_mesh(1))//'x'//int2str(rec%nk_mesh(2))// &
                         'x'//int2str(rec%nk_mesh(3))//', symmetry reduction and time reversal off'// &
                         merge(' (overriding &reciprocal)', '                         ', this%reduction_overridden), &
                         __FILE__, __LINE__)

      if (allocated(rec%eigenvalues)) deallocate (rec%eigenvalues)
      allocate (rec%eigenvalues(nmat, nk))
      do first = 1, nk, k_chunk
         last = min(nk, first + k_chunk - 1)
         call this%eigenpairs_chunk([0.0_rp, 0.0_rp, 0.0_rp], first, last, e, v)
         rec%eigenvalues(:, first:last) = e
      end do

      this%kt = rec%temperature*kB_Ry_per_K
      this%fermi_input = rec%fermi_level
      if (rec%auto_find_fermi) rec%fermi_level = rec%find_fermi_level_from_eigenvalues(rec%total_electrons)
      this%fermi_used = rec%fermi_level

      allocate (f(nk, nmat))
      do ik = 1, nk
         do n = 1, nmat
            f(ik, n) = fermi_dirac_occupation(rec%eigenvalues(n, ik), this%fermi_used, this%kt)
         end do
      end do
      this%electron_count = sum(rec%k_workset%weights*sum(f, dim=2))
      call g_logger%info('spin_response%prepare: E_F input '//real2str(this%fermi_input, '(F12.8)')//' Ry, used '// &
                         real2str(this%fermi_used, '(F12.8)')//' Ry, N = '//real2str(this%electron_count, '(F14.8)')// &
                         ' for target '//real2str(rec%total_electrons, '(F10.4)'), __FILE__, __LINE__)
      if (abs(this%electron_count - rec%total_electrons) > electron_count_tol) then
         call g_logger%fatal('spin_response%prepare: electron count on the response mesh differs from the target', &
                             __FILE__, __LINE__)
      end if

      allocate (enu(nsite, 2), p(nsite, 2))
      do i = 1, nsite
         enu(i, :) = rec%lattice%symbolic_atoms(rec%lattice%nbulk + i)%potential%enu(l_d, :)
         p(i, :) = 1.0_rp/rec%lattice%symbolic_atoms(rec%lattice%nbulk + i)%potential%ppar(l_d, :)**2
      end do

      allocate (n_mills(nsite, 2), n_juelich(nsite, 2))
      n_mills = 0.0_rp
      n_juelich = 0.0_rp
      do first = 1, nk, k_chunk
         last = min(nk, first + k_chunk - 1)
         call this%eigenpairs_chunk([0.0_rp, 0.0_rp, 0.0_rp], first, last, e, v)
         e_kn = transpose(e)
         call d_amplitudes_coefficient(v, nsite, a)
         call accumulate_d_occupation(rec%k_workset%weights(first:last), f(first:last, :), a, n_mills)
         aj = a
         call d_amplitudes_juelich(aj, e_kn, this%fermi_used, enu, p)
         call accumulate_d_occupation(rec%k_workset%weights(first:last), f(first:last, :), aj, n_juelich)
      end do

      this%spin_swapped = sum(n_mills(:, 1) - n_mills(:, 2)) < 0.0_rp
      this%spin_map = merge([2, 1], [1, 2], this%spin_swapped)
      if (allocated(this%d_occ_mills)) deallocate (this%d_occ_mills, this%d_occ_juelich, this%moment)
      allocate (this%d_occ_mills(nsite, 2), this%d_occ_juelich(nsite, 2), this%moment(nsite))
      do s = 1, 2
         this%d_occ_mills(:, s) = n_mills(:, this%spin_map(s))
         this%d_occ_juelich(:, s) = n_juelich(:, this%spin_map(s))
      end do
      this%moment = this%d_occ_juelich(:, 1) - this%d_occ_juelich(:, 2)
   end subroutine prepare

   !> @brief Write <output_prefix>_state.dat: the frozen-state values the response used.
   subroutine write_state(this)
      class(spin_response), intent(in) :: this

      integer :: funit, i

      open (newunit=funit, file=trim(this%output_prefix)//'_state.dat', action='write', status='replace')
      write (funit, '(a)') '# quantity value ; energies in Ry, moments in mu_B, occupations in electrons'
      write (funit, '(a,a)') 'method ', trim(this%method)
      write (funit, '(a,es24.16)') 'fermi_input ', this%fermi_input
      write (funit, '(a,es24.16)') 'fermi_used ', this%fermi_used
      write (funit, '(a,es24.16)') 'electron_count ', this%electron_count
      write (funit, '(a,es24.16)') 'kT ', this%kt
      write (funit, '(a,3i6)') 'mesh_nk ', this%reciprocal%nk_mesh
      write (funit, '(a,l2)') 'mesh_full_bz T ; reduction_overridden ', this%reduction_overridden
      write (funit, '(a,l2)') 'spin_swapped ', this%spin_swapped
      write (funit, '(a)') 'u_juelich n/a'
      write (funit, '(a)') 'u_mills n/a'
      write (funit, '(a)') '# site m_mills m_juelich d_up_mills d_dn_mills d_up_juelich d_dn_juelich'
      do i = 1, size(this%moment)
         write (funit, '(i6,6es24.16)') i, this%d_occ_mills(i, 1) - this%d_occ_mills(i, 2), this%moment(i), &
            this%d_occ_mills(i, :), this%d_occ_juelich(i, :)
      end do
      close (funit)
   end subroutine write_state

   !> @brief d block of the eigenvectors as amplitudes v_{i mu sigma}, with spin 1 = the up block.
   !> @param[in]  v      Eigenvectors, shape (nmat, nband, nk), site-major with up block then down block per site.
   !> @param[in]  nsite  Number of sites in the cell.
   !> @param[out] a      Amplitudes, shape (nk, nband, 2*l_d+1, nsite, 2).
   pure subroutine d_amplitudes_coefficient(v, nsite, a)
      complex(rp), intent(in) :: v(:, :, :)
      integer, intent(in) :: nsite
      complex(rp), allocatable, intent(out) :: a(:, :, :, :, :)

      integer :: i, s, mu

      allocate (a(size(v, 3), size(v, 2), 2*l_d + 1, nsite, 2))
      do s = 1, 2
         do i = 1, nsite
            do mu = 1, 2*l_d + 1
               a(:, :, mu, i, s) = transpose(v((i - 1)*nb + (s - 1)*spin_off + l_d**2 + mu, :, :))
            end do
         end do
      end do
   end subroutine d_amplitudes_coefficient

   !> @brief Multiply coefficient amplitudes by the Juelich factor
   !>        (1 + (e - enu)(ef - enu) p)/sqrt(1 + (ef - enu)^2 p), with enu and p of the spin block.
   !> @param[inout] a    Amplitudes from d_amplitudes_coefficient, in physical spin blocks.
   !> @param[in]    e    Eigenvalues (Ry), shape (nk, nband).
   !> @param[in]    ef   Fermi level (Ry).
   !> @param[in]    enu  Linearization energies of the d channel (Ry), shape (nsite, 2).
   !> @param[in]    p    <phidot^2> of the d channel (Ry^-2), shape (nsite, 2).
   pure subroutine d_amplitudes_juelich(a, e, ef, enu, p)
      complex(rp), intent(inout) :: a(:, :, :, :, :)
      real(rp), intent(in) :: e(:, :), ef, enu(:, :), p(:, :)

      integer :: i, s
      real(rp) :: x

      do s = 1, 2
         do i = 1, size(a, 4)
            x = ef - enu(i, s)
            a(:, :, :, i, s) = a(:, :, :, i, s)* &
               spread((1.0_rp + (e - enu(i, s))*x*p(i, s))/sqrt(1.0_rp + x*x*p(i, s)), 3, size(a, 3))
         end do
      end do
   end subroutine d_amplitudes_juelich

   !> @brief Add sum_k w_k sum_n f_kn sum_mu |a_{i mu sigma}|^2 to n_d(i, sigma).
   !> @param[in]    w    k weights, shape (nk).
   !> @param[in]    f    Occupations, shape (nk, nband).
   !> @param[in]    a    Amplitudes, shape (nk, nband, norb, nsite, 2).
   !> @param[inout] n_d  d occupation, shape (nsite, 2).
   pure subroutine accumulate_d_occupation(w, f, a, n_d)
      real(rp), intent(in) :: w(:), f(:, :)
      complex(rp), intent(in) :: a(:, :, :, :, :)
      real(rp), intent(inout) :: n_d(:, :)

      integer :: i, s

      do s = 1, 2
         do i = 1, size(a, 4)
            n_d(i, s) = n_d(i, s) + sum(spread(w, 2, size(f, 2))*f*sum(abs(a(:, :, :, i, s))**2, dim=3))
         end do
      end do
   end subroutine accumulate_d_occupation

   !> @brief Add the Lehmann sum of one k set to chi0 (the caller zeroes chi0 and may chunk k).
   !> @details eta = 0 is valid only at omega = 0, where a pair with |e_n - e_m| < 1e-10 Ry
   !>          takes the Fermi divided difference -f(1-f)/kt. Every other pair is pruned only
   !>          when f_n - f_m is exactly zero.
   !> @param[in]    kt     Temperature (Ry).
   !> @param[in]    w      k weights, shape (nk).
   !> @param[in]    e_k    Energies at k, shape (nk, nb); f_k occupations, same shape.
   !> @param[in]    a_k    d amplitudes at k, shape (nk, nb, norb, nsite, 2).
   !> @param[in]    e_kq   Energies at k+q, shape (nk, nb); f_kq, a_kq likewise.
   !> @param[in]    omega  Energies, shape (nw) (Ry).
   !> @param[in]    eta    Broadening (Ry).
   !> @param[inout] chi0   Bare response, shape (nsite, nsite, nw).
   pure subroutine accumulate_chi0(kt, w, e_k, f_k, a_k, e_kq, f_kq, a_kq, omega, eta, chi0)
      real(rp), intent(in) :: kt, w(:), e_k(:, :), f_k(:, :), e_kq(:, :), f_kq(:, :), omega(:), eta
      complex(rp), intent(in) :: a_k(:, :, :, :, :), a_kq(:, :, :, :, :)
      complex(rp), intent(inout) :: chi0(:, :, :)

      integer :: ik, n, m, i, j, iw
      real(rp) :: df, de
      complex(rp) :: t(size(a_k, 4)), c

      do ik = 1, size(w)
         do n = 1, size(e_k, 2)
            do m = 1, size(e_kq, 2)
               df = f_k(ik, n) - f_kq(ik, m)
               de = e_k(ik, n) - e_kq(ik, m)
               do i = 1, size(t)
                  t(i) = sum(conjg(a_k(ik, n, :, i, 1))*a_kq(ik, m, :, i, 2))
               end do
               do iw = 1, size(omega)
                  if (eta == 0.0_rp .and. omega(iw) == 0.0_rp .and. abs(de) < degenerate_tol) then
                     c = cmplx(-f_k(ik, n)*(1.0_rp - f_k(ik, n))/kt, 0.0_rp, rp)
                  else if (df == 0.0_rp) then
                     cycle
                  else
                     c = df/(omega(iw) + de + i_unit*eta)
                  end if
                  do j = 1, size(t)
                     do i = 1, size(t)
                        chi0(i, j, iw) = chi0(i, j, iw) + w(ik)*c*t(i)*conjg(t(j))
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine accumulate_chi0

   !> @brief Juelich U from the Goldstone condition sum_j chi0_ij(0,0) M_j U_j = -M_i.
   !> @param[in]  chi0s  Static bare response chi0(0,0), real, shape (nsite, nsite).
   !> @param[in]  m      Site moments, shape (nsite).
   !> @param[out] u      Interaction (Ry), shape (nsite).
   subroutine u_juelich(chi0s, m, u)
      real(rp), intent(in) :: chi0s(:, :), m(:)
      real(rp), intent(out) :: u(:)

      real(rp) :: amat(size(m), size(m))
      integer :: j, info, ipiv(size(m))

      do j = 1, size(m)
         amat(:, j) = chi0s(:, j)*m(j)
      end do
      u = -m
      call dgesv(size(m), 1, amat, size(m), ipiv, u, size(m), info)
      if (info /= 0) error stop 'u_juelich: singular Goldstone matrix'
   end subroutine u_juelich

   !> @brief Mills U = (C_dn - C_up)/M per site, from the d band centres.
   elemental function u_mills(c_up, c_dn, m) result(u)
      real(rp), intent(in) :: c_up, c_dn, m
      real(rp) :: u

      u = (c_dn - c_up)/m
   end function u_mills

   !> @brief Solve (I + chi0 U) chi = chi0 at each omega by LU, with U diagonal in sites.
   !> @details Never call at eta = 0 for q = 0: the Goldstone condition makes I + chi0 U singular there.
   !> @param[in]  chi0  Bare response, shape (nsite, nsite, nw).
   !> @param[in]  u     Interaction (Ry), shape (nsite).
   !> @param[out] chi   Enhanced response, shape (nsite, nsite, nw).
   subroutine solve_dyson(chi0, u, chi)
      complex(rp), intent(in) :: chi0(:, :, :)
      real(rp), intent(in) :: u(:)
      complex(rp), intent(out) :: chi(:, :, :)

      integer :: ns, iw, i, j, info, ipiv(size(u))
      complex(rp) :: amat(size(u), size(u))

      ns = size(u)
      do iw = 1, size(chi0, 3)
         do j = 1, ns
            amat(:, j) = chi0(:, j, iw)*u(j)
            amat(j, j) = amat(j, j) + 1.0_rp
         end do
         chi(:, :, iw) = chi0(:, :, iw)
         call zgesv(ns, ns, amat, ns, ipiv, chi(:, :, iw), ns, info)
         if (info /= 0) error stop 'solve_dyson: singular I + chi0 U'
      end do
   end subroutine solve_dyson

   !> @brief tr L(omega) with L = -(chi - chi^H)/(2 pi i), i.e. -Im tr chi / pi.
   pure function spectral_trace(chi) result(trl)
      complex(rp), intent(in) :: chi(:, :, :)
      real(rp) :: trl(size(chi, 3))

      integer :: iw, i

      do iw = 1, size(chi, 3)
         trl(iw) = 0.0_rp
         do i = 1, size(chi, 1)
            trl(iw) = trl(iw) - aimag(chi(i, i, iw))/pi
         end do
      end do
   end function spectral_trace

   !> @brief One-site pole estimates: parabolic peak of tr L and first zero of Re(1 + chi0 U) with interpolated omega > 0.
   !> @details A maximum on either grid end gives has_peak = .false.; no sign change gives
   !>          has_crossing = .false. peak and crossing are undefined unless their flag is set.
   !> @param[in]  omega  Increasing grid, shape (nw).
   !> @param[in]  trl    tr L on the grid, shape (nw).
   !> @param[in]  chi0   Bare response on the grid, shape (nw).
   !> @param[in]  u      Interaction (Ry).
   pure subroutine pole_estimates(omega, trl, chi0, u, peak, crossing, has_peak, has_crossing)
      real(rp), intent(in) :: omega(:), trl(:), u
      complex(rp), intent(in) :: chi0(:)
      real(rp), intent(out) :: peak, crossing
      logical, intent(out) :: has_peak, has_crossing

      integer :: imax, iw
      real(rp) :: h, g0, g1, root

      peak = 0.0_rp
      crossing = 0.0_rp
      imax = maxloc(trl, 1)
      has_peak = imax > 1 .and. imax < size(trl)
      if (has_peak) then
         h = 0.5_rp*(omega(imax + 1) - omega(imax - 1))
         peak = omega(imax) + h*(trl(imax - 1) - trl(imax + 1)) &
                /(2.0_rp*(trl(imax - 1) - 2.0_rp*trl(imax) + trl(imax + 1)))
      end if

      has_crossing = .false.
      do iw = 1, size(omega) - 1
         g0 = 1.0_rp + u*real(chi0(iw), rp)
         g1 = 1.0_rp + u*real(chi0(iw + 1), rp)
         if (g0*g1 < 0.0_rp) then
            root = omega(iw) - g0*(omega(iw + 1) - omega(iw))/(g1 - g0)
            if (root > 0.0_rp) then
               crossing = root
               has_crossing = .true.
               exit
            end if
         end if
      end do
   end subroutine pole_estimates

end module spin_response_mod
