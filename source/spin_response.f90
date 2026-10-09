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
   use math_mod, only: pi, i_unit, inverse_3x3
   use basis_mod, only: nb, spin_off
   use lattice_mod, only: lattice
   use reciprocal_mod, only: reciprocal, fermi_dirac_occupation, kB_Ry_per_K
   use logger_mod, only: g_logger
   use string_mod, only: int2str, real2str
   implicit none
   private

   public :: spin_response
   public :: d_amplitudes_coefficient
   public :: d_amplitudes_juelich
   public :: accumulate_d_occupation
   public :: accumulate_chi0
   public :: mills_residual
   public :: interaction_values
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
   !> Largest q list the namelist can carry.
   integer, parameter :: max_q = 64

   type :: spin_response
      class(reciprocal), pointer :: reciprocal => null()
      class(lattice), pointer :: lattice => null()
      !> 'juelich' | 'mills'
      character(len=16) :: method
      character(len=256) :: output_prefix
      integer :: n_omega
      real(rp) :: omega_min, omega_max, eta
      !> q points (direct coordinates, shape (3, nq)) and the omega grid (Ry).
      real(rp), allocatable :: q_direct(:, :), omega(:)
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
      !> Juelich and Mills U (Ry), Mills Goldstone residual, and the d exchange splitting C_dn - C_up (Ry), per site.
      real(rp), allocatable :: uj(:), um(:), delta_mills(:), gap_d(:)
      !> Position of each site relative to site 1, in direct (lattice) coordinates, shape (3, nsite).
      real(rp), allocatable :: site_offset(:, :)
   contains
      procedure :: build_from_file
      procedure :: restore_to_default
      procedure :: prepare
      procedure :: eigenpairs_chunk
      procedure, private :: set_site_offsets
      procedure :: bare_response
      procedure :: interaction
      procedure :: run
      procedure :: write_state
      procedure, private :: d_channel
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
      if (allocated(this%uj)) deallocate (this%uj, this%um, this%delta_mills, this%gap_d)
   end subroutine destructor

   subroutine restore_to_default(this)
      class(spin_response), intent(inout) :: this

      this%method = 'juelich'
      this%output_prefix = 'spin_response'
      this%n_omega = 201
      this%omega_min = 0.0_rp
      this%omega_max = 0.05_rp
      this%eta = 0.002_rp
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

      integer :: iostatus, funit, i, nlines
      real(rp) :: qtmp(3)

      include 'include_codes/namelists/spin_response.f90'

      method = this%method
      output_prefix = this%output_prefix
      n_q = 1
      q_list = 0.0_rp
      q_file = ''
      omega_min = this%omega_min
      omega_max = this%omega_max
      n_omega = this%n_omega
      eta = this%eta

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
      if (n_omega < 1 .or. (n_omega > 1 .and. omega_max <= omega_min) .or. eta < 0.0_rp) then
         call g_logger%fatal('spin_response needs n_omega >= 1, omega_max > omega_min for n_omega > 1, and eta >= 0', &
                             __FILE__, __LINE__)
      end if
      this%method = method
      this%output_prefix = output_prefix
      this%n_omega = n_omega
      this%omega_min = omega_min
      this%omega_max = omega_max
      this%eta = eta
      this%omega = [(omega_min + (i - 1)*(omega_max - omega_min)/real(max(n_omega - 1, 1), rp), i=1, n_omega)]

      if (len_trim(q_file) > 0) then
         open (newunit=funit, file=trim(q_file), action='read', iostat=iostatus, status='old')
         if (iostatus /= 0) call g_logger%fatal('q_file '//trim(q_file)//' not found', __FILE__, __LINE__)
         nlines = 0
         do
            read (funit, *, iostat=iostatus) qtmp
            if (iostatus /= 0) exit
            nlines = nlines + 1
         end do
         rewind (funit)
         allocate (this%q_direct(3, nlines))
         do i = 1, nlines
            read (funit, *) this%q_direct(:, i)
         end do
         close (funit)
      else
         if (n_q < 1 .or. n_q > max_q) call g_logger%fatal('spin_response n_q must be in 1..'//int2str(max_q), __FILE__, __LINE__)
         this%q_direct = q_list(:, 1:n_q)
      end if
   end subroutine build_from_file

   !> @brief Eigenpairs at k+q for mesh points first:last, from the reciprocal arbitrary-k service, in the gauge of H(k+q).
   !> @details H(k) is assembled at the folded k with phases exp(i 2 pi k.(r_j - r_i)), r_j - r_i the neighbour minus centre
   !>          vector in direct coordinates (reciprocal_fourier.f90, the `angle` and `phase` lines and the block placement
   !>          h(site i, site j) += ee*phase; neighbour vectors from clusba). Folding k+q by an integer vector G therefore
   !>          gives H(k_folded) = D H(k+q) D^+ with D = diag(exp(2 pi i G.tau_i)), and the eigenvector at k+q is
   !>          exp(-2 pi i G.tau_i) times the one returned. Site 1 is the reference, so one site is left untouched.
   !> @param[in]  q_direct  q in fractional coordinates of b1, b2, b3.
   !> @param[out] e         Eigenvalues (Ry), shape (nmat, last-first+1).
   !> @param[out] v         Eigenvectors, shape (nmat, nmat, last-first+1).
   subroutine eigenpairs_chunk(this, q_direct, first, last, e, v)
      class(spin_response), intent(inout) :: this
      real(rp), intent(in) :: q_direct(3)
      integer, intent(in) :: first, last
      real(rp), allocatable, intent(out) :: e(:, :)
      complex(rp), allocatable, intent(out) :: v(:, :, :)

      real(rp), allocatable :: k_folded(:, :), kq(:, :)
      real(rp) :: g(3)
      integer :: ik, i

      ! k+q is folded for this chunk only: kpoint_workset%shifted copies the whole workset on every call (quadratic in nk).
      kq = this%reciprocal%k_workset%points(:, first:last) + spread(q_direct, dim=2, ncopies=last - first + 1)
      kq = kq - floor(kq + 0.5_rp)
      call this%reciprocal%calculate_eigenpairs_at_kpoints(kq, e, v, k_folded)
      do ik = 1, last - first + 1
         g = real(nint(this%reciprocal%k_workset%points(:, first + ik - 1) + q_direct - k_folded(:, ik)), rp)
         do i = 2, this%lattice%nrec
            v((i - 1)*nb + 1:i*nb, :, ik) = v((i - 1)*nb + 1:i*nb, :, ik)* &
                                            exp(-2.0_rp*pi*i_unit*dot_product(g, this%site_offset(:, i)))
         end do
      end do
   end subroutine eigenpairs_chunk

   !> @brief Site positions relative to site 1 in direct coordinates, from the cluster positions the Hamiltonian uses.
   subroutine set_site_offsets(this)
      class(spin_response), intent(inout) :: this

      integer :: i
      real(rp) :: cell_inverse(3, 3)

      if (allocated(this%site_offset)) deallocate (this%site_offset)
      allocate (this%site_offset(3, this%lattice%nrec))
      cell_inverse = merge(this%lattice%a_cart_inv, inverse_3x3(this%lattice%a), this%lattice%a_cart_inv_ready)
      do i = 1, this%lattice%nrec
         this%site_offset(:, i) = matmul(cell_inverse, this%lattice%cr(:, this%lattice%atlist(this%lattice%ib(i))) &
                                         - this%lattice%cr(:, this%lattice%atlist(this%lattice%ib(1))))
      end do
   end subroutine set_site_offsets

   !> @brief Frozen-state handoff (C1) and d occupations and moments (C2).
   !> @details Pass 1 diagonalizes the full-BZ mesh in chunks, keeps the eigenvalues and fixes E_F
   !>          (the &energy value in reciprocal%fermi_level, or the solve on this mesh when
   !>          auto_find_fermi). The electron count is summed from the occupation array that pass 2
   !>          feeds to the d occupations. M_mills < 0 sets spin_swapped and relabels spins.
   subroutine prepare(this)
      class(spin_response), intent(inout) :: this

      class(reciprocal), pointer :: rec
      integer :: nk, nmat, nsite, first, last, ik, n, s
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
      call this%set_site_offsets()
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

      call this%d_channel(enu, p)

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

   !> @brief Linearization energy and p = 1/ppar^2 of the d channel per (site, physical spin).
   subroutine d_channel(this, enu, p)
      class(spin_response), intent(in) :: this
      real(rp), allocatable, intent(out) :: enu(:, :), p(:, :)

      integer :: i

      allocate (enu(this%lattice%nrec, 2), p(this%lattice%nrec, 2))
      do i = 1, this%lattice%nrec
         enu(i, :) = this%lattice%symbolic_atoms(this%lattice%nbulk + i)%potential%enu(l_d, :)
         p(i, :) = 1.0_rp/this%lattice%symbolic_atoms(this%lattice%nbulk + i)%potential%ppar(l_d, :)**2
      end do
   end subroutine d_channel

   !> @brief Amplitudes (spin-relabelled, Mills or Juelich) and occupations at one endpoint of a chunk.
   !> @param[out] et  Eigenvalues, shape (nk, nband); f the occupations; a the amplitudes.
   subroutine endpoint(this, e, v, juelich, enu, p, et, f, a)
      class(spin_response), intent(in) :: this
      real(rp), intent(in) :: e(:, :), enu(:, :), p(:, :)
      complex(rp), intent(in) :: v(:, :, :)
      logical, intent(in) :: juelich
      real(rp), allocatable, intent(out) :: et(:, :), f(:, :)
      complex(rp), allocatable, intent(out) :: a(:, :, :, :, :)

      integer :: ik, n

      et = transpose(e)
      call d_amplitudes_coefficient(v, this%lattice%nrec, a)
      if (juelich) call d_amplitudes_juelich(a, et, this%fermi_used, enu, p)
      a = a(:, :, :, :, this%spin_map)
      allocate (f(size(et, 1), size(et, 2)))
      do n = 1, size(et, 2)
         do ik = 1, size(et, 1)
            f(ik, n) = fermi_dirac_occupation(et(ik, n), this%fermi_used, this%kt)
         end do
      end do
   end subroutine endpoint

   !> @brief Bare response chi0(q, omega) over the full-BZ mesh, k chunk by chunk, both endpoints from eigenpairs_chunk.
   !> @param[in]  q_direct  q in fractional coordinates of b1, b2, b3.
   !> @param[in]  omega     Energies (Ry); eta (Ry) the broadening, 0 allowed at omega = 0 (static value).
   !> @param[in]  juelich   Juelich amplitudes if .true., coefficient (Mills) amplitudes otherwise.
   !> @param[out] chi0      shape (nsite, nsite, size(omega)).
   subroutine bare_response(this, q_direct, omega, eta, juelich, chi0)
      class(spin_response), intent(inout) :: this
      real(rp), intent(in) :: q_direct(3), omega(:), eta
      logical, intent(in) :: juelich
      complex(rp), intent(out) :: chi0(:, :, :)

      integer :: first, last, nk
      real(rp), allocatable :: enu(:, :), p(:, :), e_k(:, :), e_kq(:, :), et_k(:, :), et_kq(:, :), f_k(:, :), f_kq(:, :)
      complex(rp), allocatable :: v_k(:, :, :), v_kq(:, :, :), a_k(:, :, :, :, :), a_kq(:, :, :, :, :)

      call this%d_channel(enu, p)
      nk = this%reciprocal%k_workset%nk_local
      chi0 = (0.0_rp, 0.0_rp)
      do first = 1, nk, k_chunk
         last = min(nk, first + k_chunk - 1)
         call this%eigenpairs_chunk([0.0_rp, 0.0_rp, 0.0_rp], first, last, e_k, v_k)
         call this%eigenpairs_chunk(q_direct, first, last, e_kq, v_kq)
         call endpoint(this, e_k, v_k, juelich, enu, p, et_k, f_k, a_k)
         call endpoint(this, e_kq, v_kq, juelich, enu, p, et_kq, f_kq, a_kq)
         call accumulate_chi0(this%kt, this%reciprocal%k_workset%weights(first:last), et_k, f_k, a_k, &
                              et_kq, f_kq, a_kq, omega, eta, chi0)
      end do
   end subroutine bare_response

   !> @brief Juelich and Mills U per site from the static chi0(0, 0) (omega = 0, eta = 0) of each amplitude kind.
   !> @details Mills C_up, C_dn are the d entries of center_band, the on-site energy of H(k), in label order.
   subroutine interaction(this)
      class(spin_response), intent(inout) :: this

      integer :: i, ns
      real(rp), allocatable :: c_up(:)
      complex(rp), allocatable :: chi0_mills(:, :, :), chi0_juelich(:, :, :)

      ns = this%lattice%nrec
      allocate (chi0_mills(ns, ns, 1), chi0_juelich(ns, ns, 1), c_up(ns))
      allocate (this%uj(ns), this%um(ns), this%delta_mills(ns), this%gap_d(ns))
      call this%bare_response([0.0_rp, 0.0_rp, 0.0_rp], [0.0_rp], 0.0_rp, .false., chi0_mills)
      call this%bare_response([0.0_rp, 0.0_rp, 0.0_rp], [0.0_rp], 0.0_rp, .true., chi0_juelich)
      do i = 1, ns
         c_up(i) = this%lattice%symbolic_atoms(this%lattice%nbulk + i)%potential%center_band(l_d + 1, this%spin_map(1))
         this%gap_d(i) = this%lattice%symbolic_atoms(this%lattice%nbulk + i)%potential%center_band(l_d + 1, this%spin_map(2)) &
                         - c_up(i)
      end do
      call interaction_values(real(chi0_mills(:, :, 1), rp), real(chi0_juelich(:, :, 1), rp), &
                              this%d_occ_mills(:, 1) - this%d_occ_mills(:, 2), this%moment, c_up, c_up + this%gap_d, &
                              this%uj, this%um, this%delta_mills)
   end subroutine interaction

   !> @brief Loop over the q list: chi0, Dyson, tr L and pole estimates; writes <prefix>_q<NNN>.dat and <prefix>_dispersion.dat.
   subroutine run(this)
      class(spin_response), intent(inout) :: this

      integer :: iq, iw, nq, nw, ns, funit, fdisp
      real(rp) :: peak, crossing, qcart(3)
      real(rp), allocatable :: u(:), trl(:)
      complex(rp), allocatable :: chi0(:, :, :), chi(:, :, :)
      logical :: has_peak, has_crossing
      character(len=256) :: fname
      character(len=24) :: peak_text, crossing_text

      ns = this%lattice%nrec
      u = merge(this%um, this%uj, trim(this%method) == 'mills')
      nq = size(this%q_direct, 2)
      nw = size(this%omega)
      allocate (chi0(ns, ns, nw), chi(ns, ns, nw))
      open (newunit=fdisp, file=trim(this%output_prefix)//'_dispersion.dat', action='write', status='replace')
      write (fdisp, '(a)') '# q_direct(3) q_cartesian(3, units 2pi/a) |q|(1/Angstrom) peak_omega(Ry) crossing_omega(Ry)'
      do iq = 1, nq
         call this%bare_response(this%q_direct(:, iq), this%omega, this%eta, trim(this%method) == 'juelich', chi0)
         call solve_dyson(chi0, u, chi)
         trl = spectral_trace(chi)
         has_peak = .false.
         has_crossing = .false.
         if (ns == 1) call pole_estimates(this%omega, trl, chi0(1, 1, :), u(1), peak, crossing, has_peak, has_crossing)

         write (fname, '(a,a,i3.3,a)') trim(this%output_prefix), '_q', iq, '.dat'
         open (newunit=funit, file=trim(fname), action='write', status='replace')
         write (funit, '(a)') '# omega(Ry) re_tr_chi0(1/Ry) im_tr_chi0(1/Ry) re_tr_chi(1/Ry) im_tr_chi(1/Ry) tr_L(1/Ry) min_abs_eig(I+chi0*U)'
         do iw = 1, nw
            write (funit, '(7es22.14)') this%omega(iw), sum_diagonal(chi0(:, :, iw)), sum_diagonal(chi(:, :, iw)), trl(iw), &
               min_abs_eig(chi0(:, :, iw), u)
         end do
         close (funit)

         qcart = matmul(this%reciprocal%reciprocal_vectors, this%q_direct(:, iq))/(2.0_rp*pi)
         peak_text = 'n/a'
         crossing_text = 'n/a'
         if (has_peak) write (peak_text, '(es22.14)') peak
         if (has_crossing) write (crossing_text, '(es22.14)') crossing
         write (fdisp, '(7es16.8,2a24)') this%q_direct(:, iq), qcart, 2.0_rp*pi*sqrt(sum(qcart**2))/this%lattice%alat, &
            peak_text, crossing_text
      end do
      close (fdisp)
   end subroutine run

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
      write (funit, '(a)') '# site m_mills m_juelich d_up_mills d_dn_mills d_up_juelich d_dn_juelich'
      do i = 1, size(this%moment)
         write (funit, '(i6,6es24.16)') i, this%d_occ_mills(i, 1) - this%d_occ_mills(i, 2), this%moment(i), &
            this%d_occ_mills(i, :), this%d_occ_juelich(i, :)
      end do
      if (allocated(this%uj)) then
         write (funit, '(a)') '# site u_juelich(Ry) u_mills(Ry) juelich_goldstone_residual mills_residual_delta implied_gap_delta*Delta_d(Ry)'
         do i = 1, size(this%uj)
            write (funit, '(i6,2es24.16,a,es24.16,es24.16)') i, this%uj(i), this%um(i), '  0 (by construction)', &
               this%delta_mills(i), this%delta_mills(i)*this%gap_d(i)
         end do
      else
         write (funit, '(a)') 'u n/a'
      end if
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
   subroutine accumulate_chi0(kt, w, e_k, f_k, a_k, e_kq, f_kq, a_kq, omega, eta, chi0)
      real(rp), intent(in) :: kt, w(:), e_k(:, :), f_k(:, :), e_kq(:, :), f_kq(:, :), omega(:), eta
      complex(rp), intent(in) :: a_k(:, :, :, :, :), a_kq(:, :, :, :, :)
      complex(rp), intent(inout) :: chi0(:, :, :)

      integer :: ik, n, m, i, j, iw
      real(rp) :: df, de
      complex(rp) :: c
      complex(rp) :: t(size(a_k, 4), size(e_kq, 2), size(e_k, 2), size(w))

      ! vertex products do not depend on omega; each omega point then runs the whole k, n, m sum in the serial order
      do ik = 1, size(w)
         do n = 1, size(e_k, 2)
            do m = 1, size(e_kq, 2)
               do i = 1, size(t, 1)
                  t(i, m, n, ik) = sum(conjg(a_k(ik, n, :, i, 1))*a_kq(ik, m, :, i, 2))
               end do
            end do
         end do
      end do
      !$omp parallel do default(shared) private(iw, ik, n, m, i, j, df, de, c) schedule(static) if (size(omega) > 1)
      do iw = 1, size(omega)
         do ik = 1, size(w)
            do n = 1, size(e_k, 2)
               do m = 1, size(e_kq, 2)
                  df = f_k(ik, n) - f_kq(ik, m)
                  de = e_k(ik, n) - e_kq(ik, m)
                  if (eta == 0.0_rp .and. omega(iw) == 0.0_rp .and. abs(de) < degenerate_tol) then
                     c = cmplx(-f_k(ik, n)*(1.0_rp - f_k(ik, n))/kt, 0.0_rp, rp)
                  else if (df == 0.0_rp) then
                     cycle
                  else
                     c = df/(omega(iw) + de + i_unit*eta)
                  end if
                  do j = 1, size(t, 1)
                     do i = 1, size(t, 1)
                        chi0(i, j, iw) = chi0(i, j, iw) + w(ik)*c*t(i, m, n, ik)*conjg(t(j, m, n, ik))
                     end do
                  end do
               end do
            end do
         end do
      end do
      !$omp end parallel do
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

   !> @brief Mills Goldstone residual per site, delta_i = 1 + sum_j chi0s_ij u_j m_j / m_i (1 + chi0(0,0) U for one site).
   !> @param[in] chi0s  Static bare response chi0(0,0), real, shape (nsite, nsite).
   !> @param[in] u      Interaction (Ry), shape (nsite).
   !> @param[in] m      Site moments, shape (nsite).
   pure function mills_residual(chi0s, u, m) result(delta)
      real(rp), intent(in) :: chi0s(:, :), u(:), m(:)
      real(rp) :: delta(size(m))

      delta = 1.0_rp + matmul(chi0s, u*m)/m
   end function mills_residual

   !> @brief Juelich U (Goldstone solve with the Juelich chi0 and M), Mills U (centre difference over the Mills M),
   !>        and the Mills residual (Mills chi0 and M with the Mills U).
   !> @param[in]  chi0s_mills, chi0s_juelich  Static chi0(0,0) with each amplitude kind, shape (nsite, nsite).
   !> @param[in]  m_mills, m_juelich          Site moments of each kind, shape (nsite).
   !> @param[in]  c_up, c_dn                  d band centres of the majority and minority label (Ry), shape (nsite).
   subroutine interaction_values(chi0s_mills, chi0s_juelich, m_mills, m_juelich, c_up, c_dn, uj, um, delta_mills)
      real(rp), intent(in) :: chi0s_mills(:, :), chi0s_juelich(:, :), m_mills(:), m_juelich(:), c_up(:), c_dn(:)
      real(rp), intent(out) :: uj(:), um(:), delta_mills(:)

      call u_juelich(chi0s_juelich, m_juelich, uj)
      um = u_mills(c_up, c_dn, m_mills)
      delta_mills = mills_residual(chi0s_mills, um, m_mills)
   end subroutine interaction_values

   !> @brief Trace of a square matrix.
   pure function sum_diagonal(a) result(t)
      complex(rp), intent(in) :: a(:, :)
      complex(rp) :: t
      integer :: i

      t = (0.0_rp, 0.0_rp)
      do i = 1, size(a, 1)
         t = t + a(i, i)
      end do
   end function sum_diagonal

   !> @brief Smallest |eigenvalue| of I + chi0 U at one omega (gauge invariant for several sites).
   function min_abs_eig(chi0, u) result(m)
      complex(rp), intent(in) :: chi0(:, :)
      real(rp), intent(in) :: u(:)
      real(rp) :: m

      integer :: n, j, info
      complex(rp) :: a(size(u), size(u)), w(size(u)), vl(1, 1), vr(1, 1), work(4*size(u))
      real(rp) :: rwork(2*size(u))

      n = size(u)
      do j = 1, n
         a(:, j) = chi0(:, j)*u(j)
         a(j, j) = a(j, j) + 1.0_rp
      end do
      call zgeev('N', 'N', n, a, n, w, vl, 1, vr, 1, work, size(work), rwork, info)
      if (info /= 0) error stop 'min_abs_eig: zgeev failed'
      m = minval(abs(w))
   end function min_abs_eig

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
