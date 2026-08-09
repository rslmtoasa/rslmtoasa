!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Bare finite-temperature Kohn-Sham response from spinor eigenpairs.
!>
!> This is the reference eigenpair backend for TDDFT.  Frequencies and eta are
!> energies in Rydberg and the denominator is retarded,
!> `omega + e_n - e_m + i*eta`.  The right transition factor is evaluated as
!> `<m|B|n>` rather than obtained by assuming that B is Hermitian.  That
!> distinction is essential for the circular `chi^{+-}` channel, where the
!> response operator and perturbing operator are adjoints of one another.
module tddft_chi0_mod
   use precision_mod, only: rp
   use math_mod, only: pi
   use response_vertices_mod, only: response_channel, response_transition_vertex
   implicit none

   private

   ! Keep this exactly aligned with reciprocal_occupations.f90.  TDDFT-06 can
   ! later obtain these values directly from the reciprocal/configuration
   ! objects; the array provider intentionally has no reciprocal-state side
   ! effects.
   real(rp), parameter, public :: tddft_kB_Ry_per_K = 6.3336814e-6_rp
   real(rp), parameter, public :: tddft_occupation_kT_floor = 1.0e-10_rp

   !> Inputs controlling the exact all-band reference path.  A zero
   !> occupation_prune_tolerance disables pruning; a positive value is an
   !> explicitly approximate performance option.
   type, public :: tddft_chi0_options
      real(rp) :: eta = 0.0_rp
      real(rp) :: fermi_level = 0.0_rp
      real(rp) :: electronic_temperature = 0.0_rp
      integer :: band_first = 1
      integer :: band_last = 0
      real(rp) :: occupation_prune_tolerance = 0.0_rp
      integer :: k_mesh_shape(3) = 0
   end type tddft_chi0_options

   !> Reproducibility metadata written with every chi_KS output.
   type, public :: tddft_chi0_metadata
      character(len=32) :: backend = 'eigenpairs'
      character(len=32) :: energy_integration = 'not applicable'
      character(len=16) :: energy_unit = 'Rydberg'
      character(len=32) :: susceptibility_unit = '1/Rydberg'
      character(len=96) :: frequency_convention = 'retarded: omega is energy; denominator omega+en-em+i*eta'
      character(len=96) :: spectral_convention = 'Stoner spectral weight = -Im chi_KS^{+-}/pi (positive excitation positive)'
      real(rp) :: eta = 0.0_rp
      real(rp) :: fermi_level = 0.0_rp
      real(rp) :: electronic_temperature = 0.0_rp
      real(rp) :: electronic_kT = 0.0_rp
      real(rp) :: k_weight_sum = 0.0_rp
      integer :: k_mesh_shape(3) = 0
      integer :: nk = 0
      integer :: available_band_count = 0
      integer :: band_first = 0
      integer :: band_last = 0
      real(rp) :: occupation_prune_tolerance = 0.0_rp
      real(rp) :: green_eta = 0.0_rp
      real(rp) :: integration_energy_min = 0.0_rp
      real(rp) :: integration_energy_max = 0.0_rp
      integer :: integration_energy_points = 0
   end type tddft_chi0_metadata

   !> Response and directly consumable KS/Stoner spectral products.  The
   !> matrix dimensions are (left channel, right channel, frequency).  The
   !> site-diagonal array is (site channel, frequency), and is populated only
   !> for matching left/right site channels.
   type, public :: tddft_chi0_result
      complex(rp), allocatable :: chi(:, :, :)
      real(rp), allocatable :: re_chi(:, :, :)
      real(rp), allocatable :: im_chi(:, :, :)
      real(rp), allocatable :: site_diagonal_spectrum(:, :)
      real(rp), allocatable :: trace_spectrum(:)
      real(rp), allocatable :: stoner_spectral_map(:, :)
      type(tddft_chi0_metadata) :: metadata
   end type tddft_chi0_result

   public :: build_chi_ks_from_eigenpairs
   public :: tddft_fermi_occupation
   public :: write_chi_ks_text

contains

   !> Stable Fermi occupation using the existing reciprocal-space temperature
   !> convention (temperature in K, energies in Ry).  It intentionally shares
   !> reciprocal_occupations.f90's kT floor and exponential cutoffs.
   pure real(rp) function tddft_fermi_occupation(eigenvalue, fermi_level, temperature) result(occupation)
      real(rp), intent(in) :: eigenvalue, fermi_level, temperature
      real(rp) :: argument, kT

      kT = max(temperature*tddft_kB_Ry_per_K, tddft_occupation_kT_floor)
      argument = (eigenvalue - fermi_level)/kT
      if (argument >= 50.0_rp) then
         occupation = 0.0_rp
      else if (argument <= -50.0_rp) then
         occupation = 1.0_rp
      else
         occupation = 1.0_rp/(exp(argument) + 1.0_rp)
      end if
   end function tddft_fermi_occupation

   !> Build the generalized bare KS response for batches of energy-valued
   !> frequencies.  `eigenvalues_k` and `eigenvalues_kq` hold the n and m
   !> states, respectively; callers obtain the latter from TDDFT-01's exact
   !> arbitrary-k service.  k weights may be normalized probabilities or raw
   !> multiplicities and are normalized by their explicit sum.
   subroutine build_chi_ks_from_eigenpairs(k_weights, eigenvalues_k, eigenvectors_k, eigenvalues_kq, &
      eigenvectors_kq, site_orbital_counts, left_channels, right_channels, omega, options, result)
      real(rp), intent(in) :: k_weights(:)
      real(rp), intent(in) :: eigenvalues_k(:, :), eigenvalues_kq(:, :), omega(:)
      complex(rp), intent(in) :: eigenvectors_k(:, :, :), eigenvectors_kq(:, :, :)
      integer, intent(in) :: site_orbital_counts(:)
      type(response_channel), intent(in) :: left_channels(:), right_channels(:)
      type(tddft_chi0_options), intent(in) :: options
      type(tddft_chi0_result), intent(out) :: result

      integer :: nk, nbands, nspinor, nleft, nright, nw, ik, n, m, iw
      integer :: band_first, band_last
      real(rp) :: weight_sum, occupation_difference, transition_energy, prefactor
      complex(rp) :: denominator
      complex(rp), allocatable :: left_vertex(:), right_vertex(:)

      nk = size(k_weights)
      nbands = size(eigenvalues_k, 1)
      nleft = size(left_channels)
      nright = size(right_channels)
      nw = size(omega)
      nspinor = 2*sum(site_orbital_counts)
      call validate_chi_ks_inputs(nk, nbands, nspinor, nleft, nright, nw, k_weights, eigenvalues_k, &
         eigenvectors_k, eigenvalues_kq, eigenvectors_kq, options)

      band_first = options%band_first
      band_last = options%band_last
      if (band_last == 0) band_last = nbands
      if (band_first < 1 .or. band_last < band_first .or. band_last > nbands) then
         error stop 'build_chi_ks_from_eigenpairs: invalid selected band window'
      end if
      weight_sum = sum(k_weights)

      allocate(result%chi(nleft, nright, nw), result%re_chi(nleft, nright, nw), &
         result%im_chi(nleft, nright, nw))
      result%chi = cmplx(0.0_rp, 0.0_rp, rp)
      allocate(left_vertex(nleft), right_vertex(nright))

      ! This is intentionally the all-band f_n-f_m reference expression by
      ! default.  Frequency-independent vertices and numerators are built once
      ! per transition and immediately accumulated over the supplied batch.
      do ik = 1, nk
         prefactor = k_weights(ik)/weight_sum
         do n = band_first, band_last
            do m = band_first, band_last
               occupation_difference = tddft_fermi_occupation(eigenvalues_k(n, ik), options%fermi_level, &
                  options%electronic_temperature) - tddft_fermi_occupation(eigenvalues_kq(m, ik), &
                  options%fermi_level, options%electronic_temperature)
               if (options%occupation_prune_tolerance > 0.0_rp) then
                  if (abs(occupation_difference) <= options%occupation_prune_tolerance) cycle
               end if

               do iw = 1, nleft
                  left_vertex(iw) = response_transition_vertex(left_channels(iw), site_orbital_counts, &
                     eigenvectors_k(:, n, ik), eigenvectors_kq(:, m, ik))
               end do
               ! Do not replace this with conjg(vertex(right,n,m)): B is
               ! generally non-Hermitian in a circular channel.  The retarded
               ! Kubo factor is <n|A|m><m|B|n>.
               do iw = 1, nright
                  right_vertex(iw) = response_transition_vertex(right_channels(iw), site_orbital_counts, &
                     eigenvectors_kq(:, m, ik), eigenvectors_k(:, n, ik))
               end do
               transition_energy = eigenvalues_k(n, ik) - eigenvalues_kq(m, ik)
               do iw = 1, nw
                  denominator = cmplx(omega(iw) + transition_energy, options%eta, rp)
                  result%chi(:, :, iw) = result%chi(:, :, iw) + prefactor*occupation_difference* &
                     outer_product(left_vertex, right_vertex)/denominator
               end do
            end do
         end do
      end do
      deallocate(left_vertex, right_vertex)

      result%re_chi = real(result%chi, rp)
      result%im_chi = aimag(result%chi)
      call build_spectral_products(left_channels, right_channels, result)
      result%metadata%eta = options%eta
      result%metadata%fermi_level = options%fermi_level
      result%metadata%electronic_temperature = options%electronic_temperature
      result%metadata%electronic_kT = max(options%electronic_temperature*tddft_kB_Ry_per_K, &
         tddft_occupation_kT_floor)
      result%metadata%k_weight_sum = weight_sum
      result%metadata%k_mesh_shape = options%k_mesh_shape
      result%metadata%nk = nk
      result%metadata%available_band_count = nbands
      result%metadata%band_first = band_first
      result%metadata%band_last = band_last
      result%metadata%occupation_prune_tolerance = options%occupation_prune_tolerance
   end subroutine build_chi_ks_from_eigenpairs

   !> Write a self-describing plain-text chi_KS result.  The `matrix` records
   !> contain Re and Im chi; `site_diagonal` and `trace` are the explicitly
   !> labelled KS/Stoner products using -Im(chi)/pi.
   subroutine write_chi_ks_text(filename, omega, result)
      character(len=*), intent(in) :: filename
      real(rp), intent(in) :: omega(:)
      type(tddft_chi0_result), intent(in) :: result
      integer :: unit, ios, iw, ileft, iright, idiag

      if (.not. allocated(result%chi) .or. size(omega) /= size(result%chi, 3)) then
         error stop 'write_chi_ks_text: omega/result shape mismatch'
      end if
      open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
      if (ios /= 0) error stop 'write_chi_ks_text: cannot open output file'
      write(unit, '(a)') '# quantity = bare chi_KS; no TDDFT enhancement; eta is numerical broadening'
      write(unit, '(a,a)') '# energy_unit = ', trim(result%metadata%energy_unit)
      write(unit, '(a,a)') '# susceptibility_unit = ', trim(result%metadata%susceptibility_unit)
      write(unit, '(a,a)') '# frequency_convention = ', trim(result%metadata%frequency_convention)
      write(unit, '(a,a)') '# spectral_convention = ', trim(result%metadata%spectral_convention)
      write(unit, '(a,a)') '# chi0_backend = ', trim(result%metadata%backend)
      write(unit, '(a,a)') '# energy_integration = ', trim(result%metadata%energy_integration)
      write(unit, '(a,es24.16)') '# eta_Ry = ', result%metadata%eta
      write(unit, '(a,es24.16)') '# fermi_level_Ry = ', result%metadata%fermi_level
      write(unit, '(a,es24.16)') '# electronic_temperature_K = ', result%metadata%electronic_temperature
      write(unit, '(a,es24.16)') '# electronic_kT_Ry = ', result%metadata%electronic_kT
      write(unit, '(a,3(1x,i0))') '# k_mesh_shape =', result%metadata%k_mesh_shape
      write(unit, '(a,i0)') '# nk = ', result%metadata%nk
      write(unit, '(a,es24.16)') '# raw_k_weight_sum = ', result%metadata%k_weight_sum
      write(unit, '(a,i0)') '# available_band_count = ', result%metadata%available_band_count
      write(unit, '(a,2(1x,i0))') '# band_window_first_last = ', result%metadata%band_first, result%metadata%band_last
      write(unit, '(a,es24.16)') '# occupation_prune_tolerance = ', result%metadata%occupation_prune_tolerance
      if (result%metadata%integration_energy_points > 0) then
         write(unit, '(a,es24.16)') '# green_eta_Ry = ', result%metadata%green_eta
         write(unit, '(a,2(1x,es24.16))') '# integration_energy_window_Ry = ', &
            result%metadata%integration_energy_min, result%metadata%integration_energy_max
         write(unit, '(a,i0)') '# integration_energy_points = ', result%metadata%integration_energy_points
      end if
      write(unit, '(a)') '# record omega_Ry kind left_channel right_channel Re_chi_Ry^-1 Im_chi_Ry^-1 spectral_Ry^-1'
      do iw = 1, size(omega)
         do ileft = 1, size(result%chi, 1)
            do iright = 1, size(result%chi, 2)
               write(unit, '(es24.16,1x,a,2(1x,i0),3(1x,es24.16))') omega(iw), 'matrix', ileft, iright, &
                  result%re_chi(ileft, iright, iw), result%im_chi(ileft, iright, iw), &
                  -result%im_chi(ileft, iright, iw)/pi
            end do
         end do
         do idiag = 1, size(result%site_diagonal_spectrum, 1)
            write(unit, '(es24.16,1x,a,2(1x,i0),3(1x,es24.16))') omega(iw), 'site_diagonal', idiag, idiag, &
               0.0_rp, 0.0_rp, result%stoner_spectral_map(idiag, iw)
         end do
         write(unit, '(es24.16,1x,a,2(1x,i0),3(1x,es24.16))') omega(iw), 'trace', 0, 0, &
            0.0_rp, 0.0_rp, result%trace_spectrum(iw)
      end do
      close(unit)
   end subroutine write_chi_ks_text

   function outer_product(left, right) result(product)
      complex(rp), intent(in) :: left(:), right(:)
      complex(rp) :: product(size(left), size(right))
      integer :: ileft, iright

      do iright = 1, size(right)
         do ileft = 1, size(left)
            product(ileft, iright) = left(ileft)*right(iright)
         end do
      end do
   end function outer_product

   subroutine build_spectral_products(left_channels, right_channels, result)
      type(response_channel), intent(in) :: left_channels(:), right_channels(:)
      type(tddft_chi0_result), intent(inout) :: result
      integer :: nsite_diag, idiag, iw

      nsite_diag = min(size(left_channels), size(right_channels))
      allocate(result%site_diagonal_spectrum(nsite_diag, size(result%chi, 3)), &
         result%stoner_spectral_map(nsite_diag, size(result%chi, 3)), result%trace_spectrum(size(result%chi, 3)))
      result%site_diagonal_spectrum = 0.0_rp
      result%stoner_spectral_map = 0.0_rp
      do idiag = 1, nsite_diag
         if (left_channels(idiag)%site /= right_channels(idiag)%site) cycle
         do iw = 1, size(result%chi, 3)
            result%site_diagonal_spectrum(idiag, iw) = -aimag(result%chi(idiag, idiag, iw))/pi
            result%stoner_spectral_map(idiag, iw) = result%site_diagonal_spectrum(idiag, iw)
         end do
      end do
      do iw = 1, size(result%chi, 3)
         result%trace_spectrum(iw) = -sum(aimag([(result%chi(idiag, idiag, iw), idiag=1, nsite_diag)]))/pi
      end do
   end subroutine build_spectral_products

   subroutine validate_chi_ks_inputs(nk, nbands, nspinor, nleft, nright, nw, k_weights, eigenvalues_k, &
      eigenvectors_k, eigenvalues_kq, eigenvectors_kq, options)
      integer, intent(in) :: nk, nbands, nspinor, nleft, nright, nw
      real(rp), intent(in) :: k_weights(:), eigenvalues_k(:, :), eigenvalues_kq(:, :)
      complex(rp), intent(in) :: eigenvectors_k(:, :, :), eigenvectors_kq(:, :, :)
      type(tddft_chi0_options), intent(in) :: options

      if (nk <= 0 .or. nbands <= 0 .or. nspinor <= 0 .or. nleft <= 0 .or. nright <= 0 .or. nw <= 0) then
         error stop 'build_chi_ks_from_eigenpairs: empty input is not valid'
      end if
      if (size(eigenvalues_k, 2) /= nk .or. any(shape(eigenvalues_kq) /= shape(eigenvalues_k))) then
         error stop 'build_chi_ks_from_eigenpairs: eigenvalue shapes are incompatible'
      end if
      if (size(eigenvectors_k, 1) /= nspinor .or. size(eigenvectors_k, 2) /= nbands .or. &
          size(eigenvectors_k, 3) /= nk .or. any(shape(eigenvectors_kq) /= shape(eigenvectors_k))) then
         error stop 'build_chi_ks_from_eigenpairs: eigenvector shapes are incompatible'
      end if
      if (any(k_weights < 0.0_rp) .or. sum(k_weights) <= tiny(1.0_rp)) then
         error stop 'build_chi_ks_from_eigenpairs: k weights must have a positive sum'
      end if
      if (options%eta <= 0.0_rp) error stop 'build_chi_ks_from_eigenpairs: eta must be positive'
      if (options%electronic_temperature < 0.0_rp) then
         error stop 'build_chi_ks_from_eigenpairs: electronic temperature must be non-negative'
      end if
      if (options%occupation_prune_tolerance < 0.0_rp) then
         error stop 'build_chi_ks_from_eigenpairs: occupation pruning tolerance must be non-negative'
      end if
   end subroutine validate_chi_ks_inputs

end module tddft_chi0_mod
