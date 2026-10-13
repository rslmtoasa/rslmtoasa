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
!> The response space is indexed by kappa = (site - 1)*nch + channel, channel = s, p, d, f as chosen in
!> &spin_response channels; the kernels treat kappa as the site index and never see the channel.
!------------------------------------------------------------------------------
module spin_response_mod
   use precision_mod, only: rp
   use math_mod, only: pi, i_unit, inverse_3x3
   use basis_mod, only: nb, spin_off, lmax_basis
   use lattice_mod, only: lattice
   use reciprocal_mod, only: reciprocal, fermi_dirac_occupation, kB_Ry_per_K
   use logger_mod, only: g_logger
   use string_mod, only: int2str, real2str
   implicit none
   private

   public :: spin_response
   public :: channel_amplitudes_coefficient
   public :: d_amplitudes_coefficient
   public :: d_amplitudes_juelich
   public :: accumulate_d_occupation
   public :: accumulate_chi0
   public :: mills_residual
   public :: interaction_values
   public :: u_juelich
   public :: u_mills
   public :: radial_kernel_matrix
   public :: solve_dyson
   public :: adiabatic_frequencies
   public :: peak_estimate
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
   !> Omega points per matrix-multiplication block of accumulate_chi0; fixed so results do not depend on the thread count.
   integer, parameter :: omega_block = 32
   !> Points with |m| below this fraction of max|m| take K = dV/m from their neighbours in radial_kernel_matrix.
   real(rp), parameter :: m_floor = 1.0e-12_rp

   type :: spin_response
      class(reciprocal), pointer :: reciprocal => null()
      class(lattice), pointer :: lattice => null()
      !> 'juelich' | 'mills'
      character(len=16) :: method
      !> 'goldstone' (scalar U on d) | 'alsda' | 'alsda_scaled' (ALSDA matrix times the Goldstone factor lambda).
      character(len=16) :: kernel
      character(len=256) :: output_prefix
      !> l of each response channel, ascending (0 = s .. 3 = f); always contains d.
      integer, allocatable :: channel_l(:)
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
      !> The same occupations for every (site, channel) row kappa and spin label.
      real(rp), allocatable :: occ_mills(:, :), occ_juelich(:, :)
      !> Juelich d moment per site (mu_B), the one that enters the Goldstone condition.
      real(rp), allocatable :: moment(:)
      !> Juelich and Mills U (Ry), Mills Goldstone residual, and the d exchange splitting C_dn - C_up (Ry), per site.
      real(rp), allocatable :: uj(:), um(:), delta_mills(:), gap_d(:)
      !> With SOC: Juelich U and d moment per site of the pass without the L.S term (see keep_nosoc_kernel).
      real(rp), allocatable :: u_nosoc(:), moment_nosoc(:)
      !> Position of each site relative to site 1, in direct (lattice) coordinates, shape (3, nsite).
      real(rp), allocatable :: site_offset(:, :)
      !> ALSDA kernel K over kappa (Ry), unscaled, block diagonal in sites; U = -K enters the Dyson equation.
      real(rp), allocatable :: k_alsda(:, :)
      !> Interaction matrix over kappa (Ry) that run puts into the Dyson equation.
      real(rp), allocatable :: kmat(:, :)
      !> Goldstone factor of the ALSDA kernel, and the q = 0 Rayleigh eigenvalue of 1 + chi0 U along M for the unscaled kernel.
      real(rp) :: kernel_lambda, kernel_goldstone_eig
   contains
      procedure :: build_from_file
      procedure :: restore_to_default
      procedure :: prepare
      procedure :: eigenpairs_chunk
      procedure, private :: set_site_offsets
      procedure :: bare_response
      procedure :: bare_response_list
      procedure :: set_alsda_kernel
      procedure :: interaction
      procedure :: keep_nosoc_kernel
      procedure :: run
      procedure :: write_state
      procedure, private :: channel_parameters, d_rows
      final :: destructor
   end type spin_response

   interface spin_response
      procedure :: constructor
   end interface spin_response

   !> Dyson solve with the kernel as a vector (diagonal kernel) or as a matrix over kappa.
   interface solve_dyson
      procedure :: solve_dyson_diagonal, solve_dyson_matrix
   end interface solve_dyson

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
      if (allocated(this%occ_mills)) deallocate (this%occ_mills, this%occ_juelich)
      if (allocated(this%moment)) deallocate (this%moment)
      if (allocated(this%uj)) deallocate (this%uj, this%um, this%delta_mills, this%gap_d)
      if (allocated(this%k_alsda)) deallocate (this%k_alsda)
      if (allocated(this%kmat)) deallocate (this%kmat)
      if (allocated(this%u_nosoc)) deallocate (this%u_nosoc, this%moment_nosoc)
   end subroutine destructor

   subroutine restore_to_default(this)
      class(spin_response), intent(inout) :: this

      this%method = 'juelich'
      this%kernel = 'goldstone'
      this%kernel_lambda = 0.0_rp
      this%kernel_goldstone_eig = 0.0_rp
      this%channel_l = [l_d]
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

      integer :: iostatus, funit, i, nlines, l
      real(rp) :: qtmp(3)

      include 'include_codes/namelists/spin_response.f90'

      method = this%method
      kernel = this%kernel
      channels = 'd'
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
      if (trim(kernel) /= 'goldstone' .and. trim(kernel) /= 'alsda' .and. trim(kernel) /= 'alsda_scaled') then
         call g_logger%fatal('spin_response kernel must be goldstone, alsda or alsda_scaled, got '//trim(kernel), __FILE__, __LINE__)
      end if
      if (trim(kernel) /= 'goldstone' .and. this%lattice%control%txc /= 1) then
         call g_logger%fatal('spin_response kernel '//trim(kernel)//' needs txc = 1 (LDA): for GGA and libXC B_xc/m is not the transverse kernel', &
                             __FILE__, __LINE__)
      end if
      if (this%lattice%control%has_soc() .and. (trim(method) /= 'juelich' .or. trim(kernel) /= 'goldstone')) then
         call g_logger%fatal('spin_response with SOC (nsp = 2 or 4) needs method juelich and kernel goldstone', __FILE__, __LINE__)
      end if
      if (n_omega < 1 .or. (n_omega > 1 .and. omega_max <= omega_min) .or. eta < 0.0_rp) then
         call g_logger%fatal('spin_response needs n_omega >= 1, omega_max > omega_min for n_omega > 1, and eta >= 0', &
                             __FILE__, __LINE__)
      end if
      if (verify(trim(channels), 'spdf') /= 0 .or. index(channels, 'd') == 0) then
         call g_logger%fatal('spin_response channels must be letters of spdf and contain d, got '//trim(channels), __FILE__, __LINE__)
      end if
      this%channel_l = pack([0, 1, 2, 3], [(index(channels, 'spdf'(l + 1:l + 1)) > 0, l=0, 3)])
      if (maxval(this%channel_l) > lmax_basis) then
         call g_logger%fatal('spin_response channels '//trim(channels)//' need lmax >= '//int2str(maxval(this%channel_l)), &
                             __FILE__, __LINE__)
      end if
      this%method = method
      this%kernel = kernel
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
   !>          from the Fermi solve on this mesh (auto_find_fermi must be on; the &energy value is only
   !>          reported). The electron count is summed from the occupation array that pass 2
   !>          feeds to the d occupations. M_mills < 0 sets spin_swapped and relabels spins.
   subroutine prepare(this)
      class(spin_response), intent(inout) :: this

      class(reciprocal), pointer :: rec
      integer :: nk, nmat, nsite, nsc, first, last, ik, n, s
      integer, allocatable :: id(:)
      real(rp), allocatable :: f(:, :), e(:, :), e_kn(:, :), enu(:, :), p(:, :), n_mills(:, :), n_juelich(:, :)
      complex(rp), allocatable :: v(:, :, :), a(:, :, :, :, :), aj(:, :, :, :, :)

      rec => this%reciprocal
      nsite = rec%lattice%nrec
      nsc = nsite*size(this%channel_l)
      id = this%d_rows()
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
      call g_logger%info('spin_response%prepare: channels '//channel_letters(this%channel_l), __FILE__, __LINE__)

      if (allocated(rec%eigenvalues)) deallocate (rec%eigenvalues)
      allocate (rec%eigenvalues(nmat, nk))
      do first = 1, nk, k_chunk
         last = min(nk, first + k_chunk - 1)
         call this%eigenpairs_chunk([0.0_rp, 0.0_rp, 0.0_rp], first, last, e, v)
         rec%eigenvalues(:, first:last) = e
      end do

      this%kt = rec%temperature*kB_Ry_per_K
      this%fermi_input = rec%fermi_level
      if (.not. rec%auto_find_fermi) call g_logger%fatal('spin_response%prepare: E_F is always solved on the response mesh; auto_find_fermi must be .true.', &
                                                         __FILE__, __LINE__)
      rec%fermi_level = rec%find_fermi_level_from_eigenvalues(rec%total_electrons)
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

      call this%channel_parameters(enu, p)

      allocate (n_mills(nsc, 2), n_juelich(nsc, 2))
      n_mills = 0.0_rp
      n_juelich = 0.0_rp
      do first = 1, nk, k_chunk
         last = min(nk, first + k_chunk - 1)
         call this%eigenpairs_chunk([0.0_rp, 0.0_rp, 0.0_rp], first, last, e, v)
         e_kn = transpose(e)
         call channel_amplitudes_coefficient(v, nsite, this%channel_l, a)
         call accumulate_d_occupation(rec%k_workset%weights(first:last), f(first:last, :), a, n_mills)
         aj = a
         call d_amplitudes_juelich(aj, e_kn, this%fermi_used, enu, p)
         call accumulate_d_occupation(rec%k_workset%weights(first:last), f(first:last, :), aj, n_juelich)
      end do

      this%spin_swapped = sum(n_mills(id, 1) - n_mills(id, 2)) < 0.0_rp
      if (this%spin_swapped .and. this%lattice%control%has_soc()) then
         call g_logger%fatal('spin_response%prepare: with SOC the d moment must point along +z (spin_swapped is .true.)', __FILE__, __LINE__)
      end if
      this%spin_map = merge([2, 1], [1, 2], this%spin_swapped)
      if (allocated(this%d_occ_mills)) deallocate (this%d_occ_mills, this%d_occ_juelich, this%moment, this%occ_mills, this%occ_juelich)
      allocate (this%d_occ_mills(nsite, 2), this%d_occ_juelich(nsite, 2), this%moment(nsite), &
                this%occ_mills(nsc, 2), this%occ_juelich(nsc, 2))
      do s = 1, 2
         this%occ_mills(:, s) = n_mills(:, this%spin_map(s))
         this%occ_juelich(:, s) = n_juelich(:, this%spin_map(s))
      end do
      this%d_occ_mills = this%occ_mills(id, :)
      this%d_occ_juelich = this%occ_juelich(id, :)
      this%moment = this%d_occ_juelich(:, 1) - this%d_occ_juelich(:, 2)
   end subroutine prepare

   !> @brief Linearization energy and p = 1/ppar^2 of each (site, channel) row kappa and physical spin.
   subroutine channel_parameters(this, enu, p)
      class(spin_response), intent(in) :: this
      real(rp), allocatable, intent(out) :: enu(:, :), p(:, :)

      integer :: i, ic, nch

      nch = size(this%channel_l)
      allocate (enu(this%lattice%nrec*nch, 2), p(this%lattice%nrec*nch, 2))
      do i = 1, this%lattice%nrec
         do ic = 1, nch
            enu((i - 1)*nch + ic, :) = this%lattice%symbolic_atoms(this%lattice%nbulk + i)%potential%enu(this%channel_l(ic), :)
            p((i - 1)*nch + ic, :) = 1.0_rp/this%lattice%symbolic_atoms(this%lattice%nbulk + i)%potential%ppar(this%channel_l(ic), :)**2
         end do
      end do
   end subroutine channel_parameters

   !> @brief Row kappa of the d channel of each site.
   function d_rows(this) result(id)
      class(spin_response), intent(in) :: this
      integer, allocatable :: id(:)

      integer :: i

      id = [((i - 1)*size(this%channel_l) + findloc(this%channel_l, l_d, 1), i=1, this%lattice%nrec)]
   end function d_rows

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
      call channel_amplitudes_coefficient(v, this%lattice%nrec, this%channel_l, a)
      if (juelich) call d_amplitudes_juelich(a, et, this%fermi_used, enu, p)
      a = a(:, :, :, :, this%spin_map)
      allocate (f(size(et, 1), size(et, 2)))
      do n = 1, size(et, 2)
         do ik = 1, size(et, 1)
            f(ik, n) = fermi_dirac_occupation(et(ik, n), this%fermi_used, this%kt)
         end do
      end do
   end subroutine endpoint

   !> @brief Bare response chi0(q, omega) for each q of a list over the full-BZ mesh, k chunk by chunk.
   !> @details The k endpoint of a chunk is computed once and shared by all q; each q still sums the chunks in order.
   !> @param[in]  q_direct  q in fractional coordinates of b1, b2, b3, shape (3, nq).
   !> @param[in]  omega     Energies (Ry); eta (Ry) the broadening, 0 allowed at omega = 0 (static value).
   !> @param[in]  juelich   Juelich amplitudes if .true., coefficient (Mills) amplitudes otherwise.
   !> @param[out] chi0      shape (nsc, nsc, size(omega), nq), nsc = nsite*nchannel.
   !> @param[out] chi0_static  Optional chi0(q, omega = 0, eta = 0) from the same eigenpairs, shape (nsc, nsc, nq).
   subroutine bare_response_list(this, q_direct, omega, eta, juelich, chi0, chi0_static)
      class(spin_response), intent(inout) :: this
      real(rp), intent(in) :: q_direct(:, :), omega(:), eta
      logical, intent(in) :: juelich
      complex(rp), intent(out) :: chi0(:, :, :, :)
      complex(rp), intent(out), optional :: chi0_static(:, :, :)

      integer :: first, last, nk, iq
      real(rp), allocatable :: enu(:, :), p(:, :), e_k(:, :), e_kq(:, :), et_k(:, :), et_kq(:, :), f_k(:, :), f_kq(:, :)
      complex(rp), allocatable :: v_k(:, :, :), v_kq(:, :, :), a_k(:, :, :, :, :), a_kq(:, :, :, :, :)

      call this%channel_parameters(enu, p)
      nk = this%reciprocal%k_workset%nk_local
      chi0 = (0.0_rp, 0.0_rp)
      if (present(chi0_static)) chi0_static = (0.0_rp, 0.0_rp)
      do first = 1, nk, k_chunk
         last = min(nk, first + k_chunk - 1)
         call this%eigenpairs_chunk([0.0_rp, 0.0_rp, 0.0_rp], first, last, e_k, v_k)
         call endpoint(this, e_k, v_k, juelich, enu, p, et_k, f_k, a_k)
         do iq = 1, size(q_direct, 2)
            call this%eigenpairs_chunk(q_direct(:, iq), first, last, e_kq, v_kq)
            call endpoint(this, e_kq, v_kq, juelich, enu, p, et_kq, f_kq, a_kq)
            call accumulate_chi0(this%kt, this%reciprocal%k_workset%weights(first:last), et_k, f_k, a_k, &
                                 et_kq, f_kq, a_kq, omega, eta, chi0(:, :, :, iq))
            if (present(chi0_static)) then
               call accumulate_chi0(this%kt, this%reciprocal%k_workset%weights(first:last), et_k, f_k, a_k, &
                                    et_kq, f_kq, a_kq, [0.0_rp], 0.0_rp, chi0_static(:, :, iq:iq))
            end if
         end do
      end do
   end subroutine bare_response_list

   !> @brief Bare response chi0(q, omega) for one q; see bare_response_list.
   !> @param[out] chi0  shape (nsc, nsc, size(omega)).
   subroutine bare_response(this, q_direct, omega, eta, juelich, chi0)
      class(spin_response), intent(inout) :: this
      real(rp), intent(in) :: q_direct(3), omega(:), eta
      logical, intent(in) :: juelich
      complex(rp), intent(out) :: chi0(:, :, :)

      complex(rp), allocatable :: c4(:, :, :, :)

      allocate (c4(size(chi0, 1), size(chi0, 2), size(chi0, 3), 1))
      call this%bare_response_list(reshape(q_direct, [3, 1]), omega, eta, juelich, c4)
      chi0 = c4(:, :, :, 1)
   end subroutine bare_response

   !> @brief Store the ALSDA kernel block K_ab (Ry) of one site in the (site, channel) matrix.
   subroutine set_alsda_kernel(this, isite, kblock)
      class(spin_response), intent(inout) :: this
      integer, intent(in) :: isite
      real(rp), intent(in) :: kblock(:, :)

      integer :: nch, nsc

      nch = size(this%channel_l)
      nsc = this%lattice%nrec*nch
      if (.not. allocated(this%k_alsda)) then
         allocate (this%k_alsda(nsc, nsc))
         this%k_alsda = 0.0_rp
      end if
      this%k_alsda((isite - 1)*nch + 1:isite*nch, (isite - 1)*nch + 1:isite*nch) = kblock
   end subroutine set_alsda_kernel

   !> @brief Juelich and Mills U per site from the static chi0(0, 0) (omega = 0, eta = 0) of each amplitude kind, d rows only.
   !> @details Mills C_up, C_dn are the d entries of center_band, the on-site energy of H(k), in label order.
   subroutine interaction(this)
      class(spin_response), intent(inout) :: this

      integer :: i, ns, nsc
      integer, allocatable :: id(:)
      real(rp), allocatable :: c_up(:), mom(:), w(:), chi0s(:, :)
      complex(rp), allocatable :: chi0_mills(:, :, :), chi0_juelich(:, :, :)
      logical :: mills

      ns = this%lattice%nrec
      nsc = ns*size(this%channel_l)
      id = this%d_rows()
      allocate (chi0_mills(nsc, nsc, 1), chi0_juelich(nsc, nsc, 1), c_up(ns))
      if (allocated(this%uj)) deallocate (this%uj, this%um, this%delta_mills, this%gap_d)
      allocate (this%uj(ns), this%um(ns), this%delta_mills(ns), this%gap_d(ns))
      call this%bare_response([0.0_rp, 0.0_rp, 0.0_rp], [0.0_rp], 0.0_rp, .false., chi0_mills)
      call this%bare_response([0.0_rp, 0.0_rp, 0.0_rp], [0.0_rp], 0.0_rp, .true., chi0_juelich)
      do i = 1, ns
         c_up(i) = this%lattice%symbolic_atoms(this%lattice%nbulk + i)%potential%center_band(l_d + 1, this%spin_map(1))
         this%gap_d(i) = this%lattice%symbolic_atoms(this%lattice%nbulk + i)%potential%center_band(l_d + 1, this%spin_map(2)) &
                         - c_up(i)
      end do
      call interaction_values(real(chi0_mills(id, id, 1), rp), real(chi0_juelich(id, id, 1), rp), &
                              this%d_occ_mills(:, 1) - this%d_occ_mills(:, 2), this%moment, c_up, c_up + this%gap_d, &
                              this%uj, this%um, this%delta_mills)

      mills = trim(this%method) == 'mills'
      if (allocated(this%kmat)) deallocate (this%kmat)
      allocate (this%kmat(nsc, nsc))
      this%kmat = 0.0_rp
      do i = 1, ns
         this%kmat(id(i), id(i)) = merge(this%um(i), this%uj(i), mills)
         ! SOC: K = Delta/m with Delta = U m of the pass without L.S, written as U (m0/m) so that m0 = m gives U exactly.
         if (allocated(this%u_nosoc)) this%kmat(id(i), id(i)) = this%u_nosoc(i)*(this%moment_nosoc(i)/this%moment(i))
      end do
      if (trim(this%kernel) == 'goldstone') return
      ! ALSDA: U = -K. lambda is the least-squares factor along M for chi0(0,0) lambda U M = -M (equals U_Juelich for one d channel).
      mom = merge(this%occ_mills(:, 1) - this%occ_mills(:, 2), this%occ_juelich(:, 1) - this%occ_juelich(:, 2), mills)
      chi0s = real(merge(chi0_mills(:, :, 1), chi0_juelich(:, :, 1), mills), rp)
      w = matmul(chi0s, matmul(-this%k_alsda, mom))
      this%kernel_lambda = -dot_product(w, mom)/dot_product(w, w)
      this%kernel_goldstone_eig = dot_product(mom, mom + w)/dot_product(mom, mom)
      this%kmat = -this%k_alsda*merge(this%kernel_lambda, 1.0_rp, trim(this%kernel) == 'alsda_scaled')
   end subroutine interaction

   !> @brief Keep the Juelich U and d moment of the pass just run (no L.S term); the next interaction call scales U by m0/m.
   subroutine keep_nosoc_kernel(this)
      class(spin_response), intent(inout) :: this

      this%u_nosoc = this%uj
      this%moment_nosoc = this%moment
   end subroutine keep_nosoc_kernel

   !> @brief Loop over the q list: chi0, Dyson, tr L and pole estimates; writes <prefix>_q<NNN>.dat, <prefix>_dispersion.dat and <prefix>_adiabatic.dat.
   subroutine run(this)
      class(spin_response), intent(inout) :: this

      integer :: iq, iw, nq, nw, ns, nsc, funit, fdisp, fadia
      real(rp) :: peak, crossing, qcart(3), imag_max
      real(rp), allocatable :: trl(:), moment_row(:), omega_ad(:), omega_force(:)
      complex(rp), allocatable :: chi0(:, :, :, :), chi(:, :, :), chi0_static(:, :, :)
      logical :: has_peak, has_crossing
      character(len=256) :: fname
      character(len=24) :: peak_text, crossing_text

      ns = this%lattice%nrec
      nsc = ns*size(this%channel_l)
      nq = size(this%q_direct, 2)
      nw = size(this%omega)
      allocate (chi0(nsc, nsc, nw, nq), chi(nsc, nsc, nw), chi0_static(nsc, nsc, nq), omega_ad(nsc), omega_force(nsc))
      call this%bare_response_list(this%q_direct, this%omega, this%eta, trim(this%method) == 'juelich', chi0, chi0_static)
      ! signed moment of each row kappa, as in interaction
      moment_row = merge(this%occ_mills(:, 1) - this%occ_mills(:, 2), this%occ_juelich(:, 1) - this%occ_juelich(:, 2), &
                         trim(this%method) == 'mills')
      open (newunit=fadia, file=trim(this%output_prefix)//'_adiabatic.dat', action='write', status='replace')
      write (fadia, '(a)') '# q_direct(3) q_cartesian(3, units 2pi/a) |q|(1/Angstrom) omega_adiabatic(Ry, ascending, one per row of the response space) max_abs_imag(Ry) omega_force(Ry, ascending, N U (I + chi0 U) from the same static chi0)'
      open (newunit=fdisp, file=trim(this%output_prefix)//'_dispersion.dat', action='write', status='replace')
      write (fdisp, '(a)') '# q_direct(3) q_cartesian(3, units 2pi/a) |q|(1/Angstrom) peak_omega(Ry) crossing_omega(Ry) min_re_eig(I+chi0*U)_at_omega_1'
      do iq = 1, nq
         call solve_dyson(chi0(:, :, :, iq), this%kmat, chi)
         trl = spectral_trace(chi)
         has_peak = .false.
         has_crossing = .false.
         if (nsc == 1) then
            call pole_estimates(this%omega, trl, chi0(1, 1, :, iq), this%kmat(1, 1), peak, crossing, has_peak, has_crossing)
         else if (ns == 1) then
            call peak_estimate(this%omega, trl, peak, has_peak)
         end if

         write (fname, '(a,a,i3.3,a)') trim(this%output_prefix), '_q', iq, '.dat'
         open (newunit=funit, file=trim(fname), action='write', status='replace')
         write (funit, '(a)') '# omega(Ry) re_tr_chi0(1/Ry) im_tr_chi0(1/Ry) re_tr_chi(1/Ry) im_tr_chi(1/Ry) tr_L(1/Ry) min_abs_eig(I+chi0*U)'
         do iw = 1, nw
            write (funit, '(7es22.14)') this%omega(iw), sum_diagonal(chi0(:, :, iw, iq)), sum_diagonal(chi(:, :, iw)), trl(iw), &
               min_abs_eig(chi0(:, :, iw, iq), this%kmat)
         end do
         close (funit)

         qcart = matmul(this%reciprocal%reciprocal_vectors, this%q_direct(:, iq))/(2.0_rp*pi)
         peak_text = 'n/a'
         crossing_text = 'n/a'
         if (has_peak) write (peak_text, '(es22.14)') peak
         if (has_crossing) write (crossing_text, '(es22.14)') crossing
         write (fdisp, '(7es16.8,2a24,es16.8)') this%q_direct(:, iq), qcart, 2.0_rp*pi*sqrt(sum(qcart**2))/this%lattice%alat, &
            peak_text, crossing_text, min_real_eig(chi0(:, :, 1, iq), this%kmat)
         call adiabatic_frequencies(chi0_static(:, :, iq), this%kmat, moment_row, omega_ad, imag_max, omega_force)
         write (fadia, '(7es16.8,*(es22.14))') this%q_direct(:, iq), qcart, 2.0_rp*pi*sqrt(sum(qcart**2))/this%lattice%alat, &
            omega_ad, imag_max, omega_force
      end do
      close (fdisp)
      close (fadia)
   end subroutine run

   !> @brief Write <output_prefix>_state.dat: the frozen-state values the response used.
   subroutine write_state(this)
      class(spin_response), intent(in) :: this

      integer :: funit, i

      open (newunit=funit, file=trim(this%output_prefix)//'_state.dat', action='write', status='replace')
      write (funit, '(a)') '# quantity value ; energies in Ry, moments in mu_B, occupations in electrons'
      write (funit, '(a,a)') 'method ', trim(this%method)
      if (size(this%channel_l) > 1) write (funit, '(a,a)') 'channels ', channel_letters(this%channel_l)
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
      if (allocated(this%u_nosoc)) then
         write (funit, '(a)') '# SOC: K = U m0/m with U (Ry) and m0 (mu_B) from the pass without L.S, m the d moment above'
         do i = 1, size(this%u_nosoc)
            write (funit, '(a,i6,2es24.16)') 'u_m0_nosoc ', i, this%u_nosoc(i), this%moment_nosoc(i)
         end do
      end if
      if (trim(this%kernel) /= 'goldstone') then
         write (funit, '(a,a)') 'kernel ', trim(this%kernel)
         write (funit, '(a,es24.16)') 'kernel_lambda ', this%kernel_lambda
         write (funit, '(a,es24.16)') 'kernel_goldstone_eig ', this%kernel_goldstone_eig
         write (funit, '(a)') '# unscaled ALSDA kernel K, rows and columns (site, channel) (Ry)'
         do i = 1, size(this%k_alsda, 1)
            write (funit, '(*(es24.16))') this%k_alsda(i, :)
         end do
      end if
      close (funit)
   end subroutine write_state

   !> @brief Channel letters of an l list, for example 'spd'.
   pure function channel_letters(l_list) result(text)
      integer, intent(in) :: l_list(:)
      character(len=size(l_list)) :: text

      integer :: ic

      do ic = 1, size(l_list)
         text(ic:ic) = 'spdf'(l_list(ic) + 1:l_list(ic) + 1)
      end do
   end function channel_letters

   !> @brief Channel blocks of the eigenvectors as amplitudes v_{kappa mu sigma}, with spin 1 = the up block.
   !> @details Row kappa = (site - 1)*nch + channel; channel l uses orbitals l**2 + 1 .. (l + 1)**2 of the spin block,
   !>          and orbital slots beyond 2l + 1 are exact zeros.
   !> @param[in]  v       Eigenvectors, shape (nmat, nband, nk), site-major with up block then down block per site.
   !> @param[in]  nsite   Number of sites in the cell.
   !> @param[in]  l_list  Angular momenta of the channels, shape (nch).
   !> @param[out] a       Amplitudes, shape (nk, nband, 2*maxval(l_list)+1, nsite*nch, 2).
   pure subroutine channel_amplitudes_coefficient(v, nsite, l_list, a)
      complex(rp), intent(in) :: v(:, :, :)
      integer, intent(in) :: nsite, l_list(:)
      complex(rp), allocatable, intent(out) :: a(:, :, :, :, :)

      integer :: i, ic, s, mu, nch

      nch = size(l_list)
      allocate (a(size(v, 3), size(v, 2), 2*maxval(l_list) + 1, nsite*nch, 2))
      a = (0.0_rp, 0.0_rp)
      do s = 1, 2
         do i = 1, nsite
            do ic = 1, nch
               do mu = 1, 2*l_list(ic) + 1
                  a(:, :, mu, (i - 1)*nch + ic, s) = transpose(v((i - 1)*nb + (s - 1)*spin_off + l_list(ic)**2 + mu, :, :))
               end do
            end do
         end do
      end do
   end subroutine channel_amplitudes_coefficient

   !> @brief d amplitudes: channel_amplitudes_coefficient for the d channel alone, shape (nk, nband, 2*l_d+1, nsite, 2).
   pure subroutine d_amplitudes_coefficient(v, nsite, a)
      complex(rp), intent(in) :: v(:, :, :)
      integer, intent(in) :: nsite
      complex(rp), allocatable, intent(out) :: a(:, :, :, :, :)

      call channel_amplitudes_coefficient(v, nsite, [l_d], a)
   end subroutine d_amplitudes_coefficient

   !> @brief Multiply coefficient amplitudes by the Juelich factor
   !>        (1 + (e - enu)(ef - enu) p)/sqrt(1 + (ef - enu)^2 p), with enu and p of the spin block.
   !> @param[inout] a    Amplitudes from d_amplitudes_coefficient, in physical spin blocks.
   !> @param[in]    e    Eigenvalues (Ry), shape (nk, nband).
   !> @param[in]    ef   Fermi level (Ry).
   !> @param[in]    enu  Linearization energies of each row kappa (Ry), shape (nsc, 2).
   !> @param[in]    p    <phidot^2> of each row kappa (Ry^-2), shape (nsc, 2).
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
   !> @details Pairs (n, m) are listed once per call: band n needs a nonzero spin-1 amplitude at k and band m a nonzero
   !>          spin-2 amplitude at k+q (otherwise the vertex is exactly zero), and f_n - f_m must be nonzero.
   !>          eta = 0 is valid only at omega = 0, where a pair with |e_n - e_m| < 1e-10 Ry takes the Fermi divided
   !>          difference -f(1-f)/kt and stays in the list even when f_n = f_m. Nothing is pruned by a threshold.
   !>          Omega is cut into blocks of omega_block; each block is one matrix product over the pair list, so the
   !>          result does not depend on the number of threads.
   !> @param[in]    kt     Temperature (Ry).
   !> @param[in]    w      k weights, shape (nk).
   !> @param[in]    e_k    Energies at k, shape (nk, nb); f_k occupations, same shape.
   !> @param[in]    a_k    Amplitudes at k, shape (nk, nb, norb, nsc, 2), row kappa of the response space.
   !> @param[in]    e_kq   Energies at k+q, shape (nk, nb); f_kq, a_kq likewise.
   !> @param[in]    omega  Energies, shape (nw) (Ry).
   !> @param[in]    eta    Broadening (Ry).
   !> @param[inout] chi0   Bare response, shape (nsc, nsc, nw).
   subroutine accumulate_chi0(kt, w, e_k, f_k, a_k, e_kq, f_kq, a_kq, omega, eta, chi0)
      real(rp), intent(in) :: kt, w(:), e_k(:, :), f_k(:, :), e_kq(:, :), f_kq(:, :), omega(:), eta
      complex(rp), intent(in) :: a_k(:, :, :, :, :), a_kq(:, :, :, :, :)
      complex(rp), intent(inout) :: chi0(:, :, :)

      integer :: ik, n, m, i, j, p, np, nsc, ib, iw0, nbl, iw
      real(rp) :: de, x
      logical :: up(size(w), size(e_k, 2)), down(size(w), size(e_kq, 2))
      real(rp), allocatable :: dfp(:), dep(:), fdp(:)
      complex(rp), allocatable :: pm(:, :), cm(:, :), blk(:, :)
      complex(rp) :: t(size(a_k, 4))

      nsc = size(a_k, 4)
      do n = 1, size(e_k, 2)
         up(:, n) = [(any(a_k(ik, n, :, :, 1) /= (0.0_rp, 0.0_rp)), ik=1, size(w))]
      end do
      do m = 1, size(e_kq, 2)
         down(:, m) = [(any(a_kq(ik, m, :, :, 2) /= (0.0_rp, 0.0_rp)), ik=1, size(w))]
      end do
      np = sum(count(up, 2)*count(down, 2))
      allocate (pm(np, nsc*nsc), dfp(np), dep(np), fdp(np))
      np = 0
      do ik = 1, size(w)
         do n = 1, size(e_k, 2)
            if (.not. up(ik, n)) cycle
            do m = 1, size(e_kq, 2)
               if (.not. down(ik, m)) cycle
               de = e_k(ik, n) - e_kq(ik, m)
               if (f_k(ik, n) == f_kq(ik, m) .and. .not. (eta == 0.0_rp .and. abs(de) < degenerate_tol)) cycle
               np = np + 1
               dfp(np) = f_k(ik, n) - f_kq(ik, m)
               dep(np) = de
               fdp(np) = merge(-f_k(ik, n)*(1.0_rp - f_k(ik, n))/kt, 0.0_rp, abs(de) < degenerate_tol)
               do i = 1, nsc
                  t(i) = sum(conjg(a_k(ik, n, :, i, 1))*a_kq(ik, m, :, i, 2))
               end do
               do j = 1, nsc
                  do i = 1, nsc
                     pm(np, (j - 1)*nsc + i) = w(ik)*t(i)*conjg(t(j))
                  end do
               end do
            end do
         end do
      end do
      if (np == 0) return

      !$omp parallel do default(shared) private(ib, iw0, nbl, iw, p, x, cm, blk) schedule(static) if (size(omega) > omega_block)
      do ib = 1, (size(omega) + omega_block - 1)/omega_block
         iw0 = (ib - 1)*omega_block + 1
         nbl = min(omega_block, size(omega) - iw0 + 1)
         allocate (cm(np, nbl), blk(nbl, nsc*nsc))
         do iw = 1, nbl
            do p = 1, np
               x = omega(iw0 + iw - 1) + dep(p)
               if (eta == 0.0_rp .and. omega(iw0 + iw - 1) == 0.0_rp .and. abs(dep(p)) < degenerate_tol) then
                  cm(p, iw) = cmplx(fdp(p), 0.0_rp, rp)
               else if (dfp(p) == 0.0_rp) then
                  cm(p, iw) = (0.0_rp, 0.0_rp)
               else
                  cm(p, iw) = cmplx(dfp(p)*x, -dfp(p)*eta, rp)/(x*x + eta*eta)
               end if
            end do
         end do
         call zgemm('T', 'N', nbl, nsc*nsc, np, (1.0_rp, 0.0_rp), cm, np, pm, size(pm, 1), (0.0_rp, 0.0_rp), blk, nbl)
         do iw = 1, nbl
            chi0(:, :, iw0 + iw - 1) = chi0(:, :, iw0 + iw - 1) + reshape(blk(iw, :), [nsc, nsc])
         end do
         deallocate (cm, blk)
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

   !> @brief ALSDA transverse kernel matrix of one site, K_ab = (1/4 pi) int dr P_a P_b K(r)/r^2 over channels a, b, spherical vertex (L = 0).
   !> @details K(r) = dV/m with m = n_up - n_dn (core included), P_a = sqrt(gfac_up gfac_dn) G_up G_dn + S_up S_dn. The sum
   !>          starts at the second mesh point (dV = 0 at r = 0 while m is not), with the Simpson weights of the full mesh.
   !>          Where |m| <= m_floor max|m| (m = 0 on a mesh point), K is the cubic interpolation from the four nearest
   !>          points in mesh index. The Dyson equation takes U = -K.
   !> @param[in]  r, wq              Mesh (bohr) and quadrature weights including dr/di, shape (nr).
   !> @param[in]  m, dv              n_up - n_dn (bohr^-3) and V_xc,up - V_xc,dn (Ry), shape (nr).
   !> @param[in]  g_up, g_dn         Large components times r of each channel, shape (nr, nch); s_* the small components.
   !> @param[in]  gfac_up, gfac_dn   Scalar-relativistic metric factors, shape (nr, nch).
   !> @param[out] kmat               K_ab (Ry), shape (nch, nch).
   subroutine radial_kernel_matrix(r, wq, m, dv, g_up, g_dn, s_up, s_dn, gfac_up, gfac_dn, kmat)
      real(rp), intent(in) :: r(:), wq(:), m(:), dv(:), g_up(:, :), g_dn(:, :), s_up(:, :), s_dn(:, :), gfac_up(:, :), gfac_dn(:, :)
      real(rp), intent(out) :: kmat(:, :)

      integer :: i, j, nr
      logical :: guarded(size(r))
      real(rp) :: k(size(r)), p(size(r), size(g_up, 2))

      nr = size(r)
      p = sqrt(gfac_up*gfac_dn)*g_up*g_dn + s_up*s_dn
      guarded = abs(m) <= m_floor*maxval(abs(m))
      k = 0.0_rp
      do i = 2, nr
         if (.not. guarded(i)) k(i) = dv(i)/m(i)
      end do
      do i = 2, nr
         if (.not. guarded(i)) cycle
         if (i < 4 .or. i > nr - 2) error stop 'radial_kernel_matrix: m = 0 within two mesh points of an end'
         if (count(guarded(i - 2:i + 2)) > 1) error stop 'radial_kernel_matrix: m = 0 on neighbouring mesh points'
         k(i) = (-k(i - 2) + 4.0_rp*k(i - 1) + 4.0_rp*k(i + 1) - k(i + 2))/6.0_rp
      end do
      do i = 1, size(kmat, 1)
         do j = 1, size(kmat, 2)
            kmat(i, j) = sum(wq(2:)*p(2:, i)*p(2:, j)*k(2:)/r(2:)**2)/(4.0_rp*pi)
         end do
      end do
   end subroutine radial_kernel_matrix

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

   !> @brief Eigenvalues of I + chi0 K at one omega (gauge invariant for several sites).
   function eig_one_plus(chi0, kmat) result(w)
      complex(rp), intent(in) :: chi0(:, :)
      real(rp), intent(in) :: kmat(:, :)
      complex(rp) :: w(size(kmat, 1))

      integer :: n, j, info
      complex(rp) :: a(size(kmat, 1), size(kmat, 1)), vl(1, 1), vr(1, 1), work(4*size(kmat, 1))
      real(rp) :: rwork(2*size(kmat, 1))

      n = size(kmat, 1)
      a = matmul(chi0, kmat)
      do j = 1, n
         a(j, j) = a(j, j) + 1.0_rp
      end do
      call zgeev('N', 'N', n, a, n, w, vl, 1, vr, 1, work, size(work), rwork, info)
      if (info /= 0) error stop 'eig_one_plus: zgeev failed'
   end function eig_one_plus

   !> @brief Smallest |eigenvalue| of I + chi0 K at one omega.
   function min_abs_eig(chi0, kmat) result(m)
      complex(rp), intent(in) :: chi0(:, :)
      real(rp), intent(in) :: kmat(:, :)
      real(rp) :: m

      m = minval(abs(eig_one_plus(chi0, kmat)))
   end function min_abs_eig

   !> @brief Smallest real part of the eigenvalues of I + chi0 K at one omega (signed; the static delta_q at omega = 0, eta = 0).
   function min_real_eig(chi0, kmat) result(m)
      complex(rp), intent(in) :: chi0(:, :)
      real(rp), intent(in) :: kmat(:, :)
      real(rp) :: m

      m = minval(real(eig_one_plus(chi0, kmat), rp))
   end function min_real_eig

   !> @brief Solve (I + chi0 K) chi = chi0 at each omega by LU, with K a matrix over kappa.
   !> @details Never call at eta = 0 for q = 0: the Goldstone condition makes I + chi0 K singular there.
   !> @param[in]  chi0  Bare response, shape (nsc, nsc, nw).
   !> @param[in]  kmat  Interaction kernel (Ry), shape (nsc, nsc).
   !> @param[out] chi   Enhanced response, shape (nsc, nsc, nw).
   subroutine solve_dyson_matrix(chi0, kmat, chi)
      complex(rp), intent(in) :: chi0(:, :, :)
      real(rp), intent(in) :: kmat(:, :)
      complex(rp), intent(out) :: chi(:, :, :)

      integer :: ns, iw, j, info, ipiv(size(kmat, 1))
      complex(rp) :: amat(size(kmat, 1), size(kmat, 1))

      ns = size(kmat, 1)
      do iw = 1, size(chi0, 3)
         amat = matmul(chi0(:, :, iw), kmat)
         do j = 1, ns
            amat(j, j) = amat(j, j) + 1.0_rp
         end do
         chi(:, :, iw) = chi0(:, :, iw)
         call zgesv(ns, ns, amat, ns, ipiv, chi(:, :, iw), ns, info)
         if (info /= 0) error stop 'solve_dyson: singular I + chi0 K'
      end do
   end subroutine solve_dyson_matrix

   !> @brief solve_dyson_matrix with the diagonal kernel K = diag(u).
   subroutine solve_dyson_diagonal(chi0, u, chi)
      complex(rp), intent(in) :: chi0(:, :, :)
      real(rp), intent(in) :: u(:)
      complex(rp), intent(out) :: chi(:, :, :)

      real(rp) :: kmat(size(u), size(u))
      integer :: i

      kmat = 0.0_rp
      do i = 1, size(u)
         kmat(i, i) = u(i)
      end do
      call solve_dyson_matrix(chi0, kmat, chi)
   end subroutine solve_dyson_diagonal

   !> @brief Adiabatic frequencies Omega = -N (chi0^-1 + K) from the static chi0(q, 0), N = diag(m) the signed moments.
   !> @details For chi(w) = (w - Omega)^-1 N the static value is -Omega^-1 N, and (I + chi0 K) chi = chi0 gives
   !>          chi^-1 = chi0^-1 + K. Omega is similar to a Hermitian matrix, so its eigenvalues are real: positive
   !>          where m > 0 and negative where m < 0. At q = 0 the Goldstone condition makes Omega m = 0.
   !> @param[in]  chi0s       Static chi0(q, 0), shape (nsc, nsc).
   !> @param[in]  kmat        Interaction kernel (Ry), shape (nsc, nsc).
   !> @param[in]  m           Signed moment of each row kappa, shape (nsc).
   !> @param[out] omega       Real parts of the eigenvalues of Omega (Ry), ascending, shape (nsc).
   !> @param[out] imag_max    Largest |imaginary part| of those eigenvalues (Ry).
   !> @param[out] omega_force Optional: real parts of the eigenvalues of the force-theorem matrix N K (I + chi0 K) (Ry), ascending,
   !>                         from the same chi0s; -chi0^-1 of Omega is replaced by K.
   subroutine adiabatic_frequencies(chi0s, kmat, m, omega, imag_max, omega_force)
      complex(rp), intent(in) :: chi0s(:, :)
      real(rp), intent(in) :: kmat(:, :), m(:)
      real(rp), intent(out) :: omega(:), imag_max
      real(rp), intent(out), optional :: omega_force(:)

      integer :: n, i, info, ipiv(size(m))
      complex(rp) :: a(size(m), size(m)), b(size(m), size(m)), w(size(m)), vl(1, 1), vr(1, 1), work(4*size(m))
      real(rp) :: rwork(2*size(m))

      n = size(m)
      a = chi0s
      b = (0.0_rp, 0.0_rp)
      do i = 1, n
         b(i, i) = (1.0_rp, 0.0_rp)
      end do
      call zgesv(n, n, a, n, ipiv, b, n, info)
      if (info /= 0) error stop 'adiabatic_frequencies: singular static chi0'
      do i = 1, n
         a(i, :) = -m(i)*(b(i, :) + kmat(i, :))
      end do
      call zgeev('N', 'N', n, a, n, w, vl, 1, vr, 1, work, size(work), rwork, info)
      if (info /= 0) error stop 'adiabatic_frequencies: zgeev failed'
      omega = real(w, rp)
      call sort_ascending(omega)
      imag_max = maxval(abs(aimag(w)))
      if (present(omega_force)) then
         a = cmplx(kmat, 0.0_rp, rp)
         a = a + matmul(a, matmul(chi0s, a))
         do i = 1, n
            a(i, :) = m(i)*a(i, :)
         end do
         call zgeev('N', 'N', n, a, n, w, vl, 1, vr, 1, work, size(work), rwork, info)
         if (info /= 0) error stop 'adiabatic_frequencies: zgeev failed (force theorem)'
         omega_force = real(w, rp)
         call sort_ascending(omega_force)
      end if
   end subroutine adiabatic_frequencies

   !> @brief In-place ascending selection sort of a short list.
   pure subroutine sort_ascending(x)
      real(rp), intent(inout) :: x(:)

      integer :: i, j
      real(rp) :: t

      do i = 1, size(x) - 1
         j = i - 1 + minloc(x(i:), 1)
         t = x(i)
         x(i) = x(j)
         x(j) = t
      end do
   end subroutine sort_ascending

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

      integer :: iw
      real(rp) :: g0, g1, root

      crossing = 0.0_rp
      call peak_estimate(omega, trl, peak, has_peak)

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

   !> @brief Parabolic peak of tr L; a maximum on either grid end gives has_peak = .false. and peak undefined (0).
   pure subroutine peak_estimate(omega, trl, peak, has_peak)
      real(rp), intent(in) :: omega(:), trl(:)
      real(rp), intent(out) :: peak
      logical, intent(out) :: has_peak

      integer :: imax
      real(rp) :: h

      peak = 0.0_rp
      imax = maxloc(trl, 1)
      has_peak = imax > 1 .and. imax < size(trl)
      if (has_peak) then
         h = 0.5_rp*(omega(imax + 1) - omega(imax - 1))
         peak = omega(imax) + h*(trl(imax - 1) - trl(imax + 1)) &
                /(2.0_rp*(trl(imax - 1) - 2.0_rp*trl(imax) + trl(imax + 1)))
      end if
   end subroutine peak_estimate

end module spin_response_mod
