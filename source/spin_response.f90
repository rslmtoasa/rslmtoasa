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
   public :: read_q_points
   public :: read_omega_grid

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

   !> Pure kernels on plain arrays, implemented in spin_response_kernels.f90.
   interface
      module subroutine accumulate_chi0(kt, w, e_k, f_k, a_k, e_kq, f_kq, a_kq, omega, eta, chi0)
         real(rp), intent(in) :: kt, w(:), e_k(:, :), f_k(:, :), e_kq(:, :), f_kq(:, :), omega(:), eta
         complex(rp), intent(in) :: a_k(:, :, :, :, :), a_kq(:, :, :, :, :)
         complex(rp), intent(inout) :: chi0(:, :, :)
      end subroutine accumulate_chi0
      module subroutine u_juelich(chi0s, m, u)
         real(rp), intent(in) :: chi0s(:, :), m(:)
         real(rp), intent(out) :: u(:)
      end subroutine u_juelich
      elemental module function u_mills(c_up, c_dn, m) result(u)
         real(rp), intent(in) :: c_up, c_dn, m
         real(rp) :: u
      end function u_mills
      module subroutine radial_kernel_matrix(r, wq, m, dv, g_up, g_dn, s_up, s_dn, gfac_up, gfac_dn, kmat)
         real(rp), intent(in) :: r(:), wq(:), m(:), dv(:), g_up(:, :), g_dn(:, :), s_up(:, :), s_dn(:, :), gfac_up(:, :), gfac_dn(:, :)
         real(rp), intent(out) :: kmat(:, :)
      end subroutine radial_kernel_matrix
      pure module function mills_residual(chi0s, u, m) result(delta)
         real(rp), intent(in) :: chi0s(:, :), u(:), m(:)
         real(rp) :: delta(size(m))
      end function mills_residual
      module subroutine interaction_values(chi0s_mills, chi0s_juelich, m_mills, m_juelich, c_up, c_dn, uj, um, delta_mills)
         real(rp), intent(in) :: chi0s_mills(:, :), chi0s_juelich(:, :), m_mills(:), m_juelich(:), c_up(:), c_dn(:)
         real(rp), intent(out) :: uj(:), um(:), delta_mills(:)
      end subroutine interaction_values
      pure module function sum_diagonal(a) result(t)
         complex(rp), intent(in) :: a(:, :)
         complex(rp) :: t
         integer :: i
      end function sum_diagonal
      module function min_abs_eig(chi0, kmat) result(m)
         complex(rp), intent(in) :: chi0(:, :)
         real(rp), intent(in) :: kmat(:, :)
         real(rp) :: m
      end function min_abs_eig
      module function min_real_eig(chi0, kmat) result(m)
         complex(rp), intent(in) :: chi0(:, :)
         real(rp), intent(in) :: kmat(:, :)
         real(rp) :: m
      end function min_real_eig
      module subroutine solve_dyson_matrix(chi0, kmat, chi)
         complex(rp), intent(in) :: chi0(:, :, :)
         real(rp), intent(in) :: kmat(:, :)
         complex(rp), intent(out) :: chi(:, :, :)
      end subroutine solve_dyson_matrix
      module subroutine solve_dyson_diagonal(chi0, u, chi)
         complex(rp), intent(in) :: chi0(:, :, :)
         real(rp), intent(in) :: u(:)
         complex(rp), intent(out) :: chi(:, :, :)
      end subroutine solve_dyson_diagonal
      module subroutine adiabatic_frequencies(chi0s, kmat, m, omega, imag_max, omega_force)
         complex(rp), intent(in) :: chi0s(:, :)
         real(rp), intent(in) :: kmat(:, :), m(:)
         real(rp), intent(out) :: omega(:), imag_max
         real(rp), intent(out), optional :: omega_force(:)
      end subroutine adiabatic_frequencies
      pure module function spectral_trace(chi) result(trl)
         complex(rp), intent(in) :: chi(:, :, :)
         real(rp) :: trl(size(chi, 3))
      end function spectral_trace
      pure module subroutine pole_estimates(omega, trl, chi0, u, peak, crossing, has_peak, has_crossing)
         real(rp), intent(in) :: omega(:), trl(:), u
         complex(rp), intent(in) :: chi0(:)
         real(rp), intent(out) :: peak, crossing
         logical, intent(out) :: has_peak, has_crossing
      end subroutine pole_estimates
      pure module subroutine peak_estimate(omega, trl, peak, has_peak)
         real(rp), intent(in) :: omega(:), trl(:)
         real(rp), intent(out) :: peak
         logical, intent(out) :: has_peak
      end subroutine peak_estimate
   end interface

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

      integer :: iostatus, funit, i, l

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

      call q_points(n_q, q_list, q_file, this%q_direct)
   end subroutine build_from_file

   !> @brief q points (direct coordinates, shape (3, nq)) from q_file (one q per line) when it is set, else q_list(:, 1:n_q).
   subroutine q_points(n_q, q_list, q_file, q_direct)
      integer, intent(in) :: n_q
      real(rp), intent(in) :: q_list(:, :)
      character(len=*), intent(in) :: q_file
      real(rp), allocatable, intent(out) :: q_direct(:, :)

      integer :: funit, iostatus, nlines, i
      real(rp) :: qtmp(3)

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
         allocate (q_direct(3, nlines))
         do i = 1, nlines
            read (funit, *) q_direct(:, i)
         end do
         close (funit)
      else
         if (n_q < 1 .or. n_q > max_q) call g_logger%fatal('spin_response n_q must be in 1..'//int2str(max_q), __FILE__, __LINE__)
         q_direct = q_list(:, 1:n_q)
      end if
   end subroutine q_points

   !> @brief q points of the &spin_response namelist of file fname (n_q, q_list, q_file; default q = 0), see q_points.
   subroutine read_q_points(fname, q_direct)
      character(len=*), intent(in) :: fname
      real(rp), allocatable, intent(out) :: q_direct(:, :)

      integer :: funit, iostatus

      include 'include_codes/namelists/spin_response.f90'

      n_q = 1
      q_list = 0.0_rp
      q_file = ''
      open (newunit=funit, file=fname, action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) call g_logger%fatal('file '//trim(fname)//' not found', __FILE__, __LINE__)
      read (funit, nml=spin_response, iostat=iostatus)
      if (iostatus /= 0 .and. .not. is_iostat_end(iostatus)) then
         call g_logger%fatal('error while reading namelist spin_response, iostat = '//int2str(iostatus), __FILE__, __LINE__)
      end if
      close (funit)
      call q_points(n_q, q_list, q_file, q_direct)
   end subroutine read_q_points

   !> @brief omega grid (Ry) of the &spin_response namelist of file fname: n_omega points from omega_min to omega_max (defaults 0,
   !>        0.05, 201).
   subroutine read_omega_grid(fname, omega)
      character(len=*), intent(in) :: fname
      real(rp), allocatable, intent(out) :: omega(:)

      integer :: funit, iostatus, i

      include 'include_codes/namelists/spin_response.f90'

      omega_min = 0.0_rp
      omega_max = 0.05_rp
      n_omega = 201
      open (newunit=funit, file=fname, action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) call g_logger%fatal('file '//trim(fname)//' not found', __FILE__, __LINE__)
      read (funit, nml=spin_response, iostat=iostatus)
      if (iostatus /= 0 .and. .not. is_iostat_end(iostatus)) then
         call g_logger%fatal('error while reading namelist spin_response, iostat = '//int2str(iostatus), __FILE__, __LINE__)
      end if
      close (funit)
      if (n_omega < 1 .or. (n_omega > 1 .and. omega_max <= omega_min)) then
         call g_logger%fatal('spin_response needs n_omega >= 1 and omega_max > omega_min for n_omega > 1', __FILE__, __LINE__)
      end if
      omega = [(omega_min + (i - 1)*(omega_max - omega_min)/real(max(n_omega - 1, 1), rp), i=1, n_omega)]
   end subroutine read_omega_grid

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
      if (this%kt <= 0.0_rp) call g_logger%fatal('spin_response%prepare: needs &reciprocal temperature > 0 (dos_method = ''blochl'' '// &
                                                 'forces 0 K): the static chi0 of degenerate pairs is -f(1-f)/kT', __FILE__, __LINE__)
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
         write (fdisp, '(7es16.8,2(1x,a24),es16.8)') this%q_direct(:, iq), qcart, 2.0_rp*pi*sqrt(sum(qcart**2))/this%lattice%alat, &
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

end module spin_response_mod
