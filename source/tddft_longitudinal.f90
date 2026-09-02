!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Coupled charge/longitudinal-spin response and legacy static helpers.
!>
!> The production adapters below reuse the generalized TDDFT chi_KS engine and
!> tddft_dyson_mod.  The original independently converged +/- Bz calibration
!> and low-frequency fit are retained as TDDFT-08 compatibility APIs, while
!> TDDFT-13 provides the coupled charge/longitudinal kernel assembly.
module tddft_longitudinal_mod
   use precision_mod, only: rp
   use response_components_mod, only: RESPONSE_CHARGE, RESPONSE_MZ
   use response_vertices_mod, only: response_channel
   use tddft_chi0_mod, only: tddft_chi0_options, tddft_chi0_result, build_chi_ks_from_eigenpairs
   use tddft_chi0_green_mod, only: green_function_provider, green_chi0_options, build_chi_ks_from_green_functions
   use xc_response_kernel_mod, only: xc_response_kernel_provider
   implicit none

   private

   type, public :: tddft_longitudinal_options
      real(rp) :: pair_tolerance = 1.0e-10_rp
      real(rp) :: linearity_tolerance = 5.0e-2_rp
      real(rp) :: static_agreement_tolerance = 5.0e-2_rp
      real(rp) :: fit_omega_min = 0.0_rp
      real(rp) :: fit_omega_max = huge(1.0_rp)
   end type tddft_longitudinal_options

   !> Static response reconstructed column-by-column from symmetric fields.
   !> `linearity_error` compares every available |Delta B| estimate with the
   !> smallest-field one, so a user cannot silently accept a nonlinear field.
   type, public :: tddft_longitudinal_static_result
      real(rp), allocatable :: chi(:, :)
      real(rp), allocatable :: field_steps(:, :)
      real(rp), allocatable :: linearity_error(:)
      logical :: linearity_passed = .false.
   end type tddft_longitudinal_static_result

   type, public :: tddft_longitudinal_result
      real(rp), allocatable :: m0(:)
      complex(rp), allocatable :: chi_static(:, :)
      complex(rp), allocatable :: u_parallel(:, :)
      real(rp), allocatable :: t_parallel(:)
      real(rp), allocatable :: gamma_parallel(:)
      real(rp), allocatable :: fit_residual(:)
      real(rp), allocatable :: dynamic_static_error(:)
      integer :: fit_first = 0
      integer :: fit_last = 0
      logical :: static_agreement_passed = .false.
   end type tddft_longitudinal_result

   public :: read_longitudinal_static_fields
   public :: build_longitudinal_static_response
   public :: build_longitudinal_kernel
   public :: calibrate_longitudinal_response
   public :: write_longitudinal_report
   public :: longitudinal_index
   public :: build_charge_longitudinal_channels
   public :: build_charge_longitudinal_chi_ks
   public :: build_charge_longitudinal_chi_ks_from_green_functions
   public :: build_charge_longitudinal_kernel
   public :: append_longitudinal_response_metadata

contains

   !> Site-major index for the coupled longitudinal basis:
   !> (charge(site=1), m_z(site=1), charge(site=2), m_z(site=2), ...).
   pure integer function longitudinal_index(site, component) result(index)
      integer, intent(in) :: site, component

      if (site < 1 .or. (component /= RESPONSE_CHARGE .and. component /= RESPONSE_MZ)) then
         error stop 'longitudinal_index: invalid site or component'
      end if
      index = 2*(site-1) + merge(1, 2, component == RESPONSE_CHARGE)
   end function longitudinal_index

   subroutine build_charge_longitudinal_channels(nsite, channels)
      integer, intent(in) :: nsite
      type(response_channel), allocatable, intent(out) :: channels(:)
      integer :: isite

      if (nsite < 1) error stop 'build_charge_longitudinal_channels: nsite must be positive'
      allocate(channels(2*nsite))
      do isite = 1, nsite
         channels(longitudinal_index(isite, RESPONSE_CHARGE)) = response_channel(isite, RESPONSE_CHARGE)
         channels(longitudinal_index(isite, RESPONSE_MZ)) = response_channel(isite, RESPONSE_MZ)
      end do
   end subroutine build_charge_longitudinal_channels

   !> Common eigenpair adapter for the coupled charge/longitudinal basis.
   subroutine build_charge_longitudinal_chi_ks(k_weights, eigenvalues_k, eigenvectors_k, eigenvalues_kq, eigenvectors_kq, &
      site_orbital_counts, omega, options, result)
      real(rp), intent(in) :: k_weights(:), eigenvalues_k(:, :), eigenvalues_kq(:, :), omega(:)
      complex(rp), intent(in) :: eigenvectors_k(:, :, :), eigenvectors_kq(:, :, :)
      integer, intent(in) :: site_orbital_counts(:)
      type(tddft_chi0_options), intent(in) :: options
      type(tddft_chi0_result), intent(out) :: result
      type(response_channel), allocatable :: channels(:)

      call build_charge_longitudinal_channels(size(site_orbital_counts), channels)
      call build_chi_ks_from_eigenpairs(k_weights, eigenvalues_k, eigenvectors_k, eigenvalues_kq, eigenvectors_kq, &
         site_orbital_counts, channels, channels, omega, options, result)
   end subroutine build_charge_longitudinal_chi_ks

   !> Common K-space Green-function adapter.  Native real-space GF remains a
   !> deliberate capability guard in the production driver until its
   !> multi-component Fourier source is derived.
   subroutine build_charge_longitudinal_chi_ks_from_green_functions(one_particle, k_weights, site_orbital_counts, omega, &
      options, result)
      class(green_function_provider), intent(in) :: one_particle
      real(rp), intent(in) :: k_weights(:), omega(:)
      integer, intent(in) :: site_orbital_counts(:)
      type(green_chi0_options), intent(in) :: options
      type(tddft_chi0_result), intent(out) :: result
      type(response_channel), allocatable :: channels(:)

      call build_charge_longitudinal_channels(size(site_orbital_counts), channels)
      call build_chi_ks_from_green_functions(one_particle, k_weights, site_orbital_counts, channels, channels, omega, &
         options, result)
   end subroutine build_charge_longitudinal_chi_ks_from_green_functions

   !> Assemble K=f_H+f_xc in the same site-major basis as the bare response.
   !> Hartree is charge-charge only; the local XC block is the ground-state
   !> XC functional projected (n,m_z) ALSDA Hessian.  No Goldstone projection
   !> is applied because charge/longitudinal fluctuations are not rigid-spin
   !> rotations.
   subroutine build_charge_longitudinal_kernel(provider, coulomb_site, kernel)
      type(xc_response_kernel_provider), intent(in) :: provider
      real(rp), intent(in) :: coulomb_site(:, :)
      complex(rp), intent(out) :: kernel(:, :)
      integer :: nsite, isite, jsite, icharge, ispin, jcharge
      logical :: supported
      character(len=256) :: reason

      if (.not. allocated(provider%site)) error stop 'build_charge_longitudinal_kernel: XC provider is not initialized'
      nsite = size(provider%site)
      if (nsite < 1 .or. any(shape(coulomb_site) /= [nsite, nsite]) .or. &
          any(shape(kernel) /= [2*nsite, 2*nsite])) then
         error stop 'build_charge_longitudinal_kernel: incompatible site or response dimensions'
      end if
      call provider%longitudinal_response_capability(supported, reason)
      if (.not. supported) error stop 'build_charge_longitudinal_kernel: '//trim(reason)

      kernel = cmplx(0.0_rp, 0.0_rp, rp)
      do isite = 1, nsite
         icharge = longitudinal_index(isite, RESPONSE_CHARGE)
         ispin = longitudinal_index(isite, RESPONSE_MZ)
         kernel(icharge, icharge) = cmplx(provider%site(isite)%dvxc_dn, 0.0_rp, rp)
         kernel(icharge, ispin) = cmplx(provider%site(isite)%dvxc_dm, 0.0_rp, rp)
         kernel(ispin, icharge) = cmplx(provider%site(isite)%dbxc_dn, 0.0_rp, rp)
         kernel(ispin, ispin) = cmplx(provider%site(isite)%dbxc_dm, 0.0_rp, rp)
      end do
      do isite = 1, nsite
         icharge = longitudinal_index(isite, RESPONSE_CHARGE)
         do jsite = 1, nsite
            jcharge = longitudinal_index(jsite, RESPONSE_CHARGE)
            kernel(icharge, jcharge) = kernel(icharge, jcharge) + cmplx(coulomb_site(isite, jsite), 0.0_rp, rp)
         end do
      end do
   end subroutine build_charge_longitudinal_kernel

   subroutine append_longitudinal_response_metadata(filename, kernel, functional_label)
      character(len=*), intent(in) :: filename
      complex(rp), intent(in) :: kernel(:, :)
      character(len=*), intent(in) :: functional_label
      integer :: unit, ios

      open(newunit=unit, file=filename, status='old', position='append', action='write', iostat=ios)
      if (ios /= 0) error stop 'append_longitudinal_response_metadata: cannot open output file'
      write(unit, '(a)') '# longitudinal_response_metadata_begin'
      write(unit, '(a)') '# response_basis = site-major (charge,m_z) coupled longitudinal block'
      write(unit, '(a,i0)') '# response_dimension = ', size(kernel, 1)
      write(unit, '(a)') '# response_index = 2*(site-1)+1 charge; 2*(site-1)+2 m_z'
      write(unit, '(a)') '# kernel_provenance = ground-state XC local ALSDA Hessian plus projected Hartree'
      write(unit, '(a)') '# hartree_block = charge-charge only; supplied in Rydberg response metric'
      write(unit, '(a,a)') '# xc_functional = ', trim(functional_label)
      write(unit, '(a)') '# goldstone_constraint = none; longitudinal fluctuations are not rigid rotations'
      write(unit, '(a)') '# eta_role = numerical broadening only; no damping or linewidth is inferred'
      write(unit, '(a)') '# llb_mapping = not performed; susceptibility output is the input to future dissipative analysis'
      write(unit, '(a)') '# longitudinal_response_metadata_end'
      close(unit)
   end subroutine append_longitudinal_response_metadata

   !> Read independently converged SCF results.  The first non-comment record
   !> is `nsite nrecords`; each following record is
   !> `perturbed_site signed_DeltaB_Ry m_z(site=1) ... m_z(site=nsite)`.
   !> Every source site must have matched +B and -B records, with one or more
   !> magnitudes.  Repeating this for all source sites gives a full site-site
   !> response rather than only its on-site diagonal.
   subroutine read_longitudinal_static_fields(filename, nsite, source, field, moments)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nsite
      integer, allocatable, intent(out) :: source(:)
      real(rp), allocatable, intent(out) :: field(:), moments(:, :)
      integer :: unit, ios, nfile, nrecords, ir
      character(len=1024) :: line

      open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) error stop 'read_longitudinal_static_fields: cannot open static-field file'
      do
         read(unit, '(A)', iostat=ios) line
         if (ios /= 0) error stop 'read_longitudinal_static_fields: missing header'
         line = adjustl(line)
         if (len_trim(line) > 0 .and. line(1:1) /= '#') exit
      end do
      read(line, *, iostat=ios) nfile, nrecords
      if (ios /= 0 .or. nfile /= nsite .or. nrecords < 2) then
         error stop 'read_longitudinal_static_fields: invalid nsite/nrecords header'
      end if
      allocate(source(nrecords), field(nrecords), moments(nsite, nrecords))
      ir = 0
      do while (ir < nrecords)
         read(unit, '(A)', iostat=ios) line
         if (ios /= 0) error stop 'read_longitudinal_static_fields: truncated static-field data'
         line = adjustl(line)
         if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
         ir = ir + 1
         read(line, *, iostat=ios) source(ir), field(ir), moments(:, ir)
         if (ios /= 0) error stop 'read_longitudinal_static_fields: invalid field record'
      end do
      close(unit)
   end subroutine read_longitudinal_static_fields

   subroutine build_longitudinal_static_response(source, field, moments, options, result)
      integer, intent(in) :: source(:)
      real(rp), intent(in) :: field(:), moments(:, :)
      type(tddft_longitudinal_options), intent(in) :: options
      type(tddft_longitudinal_static_result), intent(out) :: result
      integer :: nsite, nrecord, isource, ir, jr, npair, ipair, min_pair
      integer, allocatable :: plus_record(:), minus_record(:)
      real(rp), allocatable :: slopes(:, :)
      real(rp) :: scale

      nsite = size(moments, 1); nrecord = size(moments, 2)
      if (size(source) /= nrecord .or. size(field) /= nrecord .or. nsite < 1) then
         error stop 'build_longitudinal_static_response: incompatible field data'
      end if
      if (options%pair_tolerance <= 0.0_rp .or. options%linearity_tolerance < 0.0_rp) then
         error stop 'build_longitudinal_static_response: invalid tolerances'
      end if
      allocate(result%chi(nsite, nsite), result%field_steps(nsite, nrecord/2), result%linearity_error(nsite))
      result%chi = 0.0_rp; result%field_steps = 0.0_rp; result%linearity_error = 0.0_rp
      do isource = 1, nsite
         npair = 0
         do ir = 1, nrecord
            if (source(ir) /= isource .or. field(ir) <= 0.0_rp) cycle
            do jr = 1, nrecord
               if (source(jr) /= isource .or. field(jr) >= 0.0_rp) cycle
               if (abs(field(ir) + field(jr)) <= options%pair_tolerance*max(1.0_rp, abs(field(ir)))) then
                  npair = npair + 1
                  exit
               end if
            end do
         end do
         if (npair < 1) error stop 'build_longitudinal_static_response: each site needs matched +DeltaB/-DeltaB data'
         allocate(plus_record(npair), minus_record(npair), slopes(nsite, npair))
         ipair = 0
         do ir = 1, nrecord
            if (source(ir) /= isource .or. field(ir) <= 0.0_rp) cycle
            do jr = 1, nrecord
               if (source(jr) /= isource .or. field(jr) >= 0.0_rp) cycle
               if (abs(field(ir) + field(jr)) <= options%pair_tolerance*max(1.0_rp, abs(field(ir)))) then
                  ipair = ipair + 1; plus_record(ipair) = ir; minus_record(ipair) = jr
                  slopes(:, ipair) = (moments(:, ir) - moments(:, jr))/(2.0_rp*field(ir))
                  result%field_steps(isource, ipair) = field(ir)
                  exit
               end if
            end do
         end do
         min_pair = minloc(field(plus_record), dim=1)
         result%chi(:, isource) = sum(slopes, dim=2)/real(npair, rp)
         do ipair = 1, npair
            scale = max(1.0_rp, maxval(abs(slopes(:, min_pair))))
            result%linearity_error(isource) = max(result%linearity_error(isource), &
               maxval(abs(slopes(:, ipair) - slopes(:, min_pair)))/scale)
         end do
         deallocate(plus_record, minus_record, slopes)
      end do
      result%linearity_passed = all(result%linearity_error <= options%linearity_tolerance)
      if (.not. result%linearity_passed) error stop 'build_longitudinal_static_response: DeltaB->0 linearity check failed'
   end subroutine build_longitudinal_static_response

   !> Calibrate U = chi_KS(0)^-1 - chi_static^-1, then fit each on-site
   !> response to chi(0)/(1 + i omega T).  The plus sign is the translation of
   !> the requested relaxational form to this code’s retarded +i*eta Kubo
   !> convention, where a positive-frequency absorptive diagonal has Im chi<0.
   subroutine calibrate_longitudinal_response(m0, chi_ks_static, chi_static, omega, chi_dynamic, options, result)
      real(rp), intent(in) :: m0(:), omega(:)
      complex(rp), intent(in) :: chi_ks_static(:, :), chi_static(:, :), chi_dynamic(:, :, :)
      type(tddft_longitudinal_options), intent(in) :: options
      type(tddft_longitudinal_result), intent(out) :: result
      complex(rp), allocatable :: model(:)
      integer :: nsite, nw, i, iw, nfit
      real(rp) :: denom, residual_scale

      nsite = size(m0); nw = size(omega)
      if (nsite < 1 .or. any(shape(chi_ks_static) /= [nsite, nsite]) .or. &
          any(shape(chi_static) /= [nsite, nsite]) .or. any(shape(chi_dynamic) /= [nsite, nsite, nw])) then
         error stop 'calibrate_longitudinal_response: incompatible response dimensions'
      end if
      if (options%fit_omega_max < options%fit_omega_min .or. options%static_agreement_tolerance < 0.0_rp) then
         error stop 'calibrate_longitudinal_response: invalid fit/static tolerances'
      end if
      allocate(result%m0(nsite), result%chi_static(nsite, nsite), result%u_parallel(nsite, nsite), &
         result%t_parallel(nsite), result%gamma_parallel(nsite), result%fit_residual(nsite), &
         result%dynamic_static_error(nsite))
      call build_longitudinal_kernel(chi_ks_static, chi_static, result%u_parallel)
      result%m0 = m0; result%chi_static = chi_static
      result%t_parallel = 0.0_rp; result%gamma_parallel = 0.0_rp; result%fit_residual = 0.0_rp
      result%dynamic_static_error = 0.0_rp
      result%fit_first = 0; result%fit_last = 0
      do iw = 1, nw
         if (omega(iw) >= options%fit_omega_min .and. omega(iw) <= options%fit_omega_max .and. omega(iw) > 0.0_rp) then
            if (result%fit_first == 0) result%fit_first = iw
            result%fit_last = iw
         end if
      end do
      if (result%fit_first == 0) error stop 'calibrate_longitudinal_response: fit needs positive frequencies in range'
      nfit = result%fit_last - result%fit_first + 1
      do i = 1, nsite
         denom = 0.0_rp
         do iw = result%fit_first, result%fit_last
            denom = denom + omega(iw)**2
            result%t_parallel(i) = result%t_parallel(i) + omega(iw)*aimag(chi_static(i, i)/chi_dynamic(i, i, iw))
         end do
         result%t_parallel(i) = result%t_parallel(i)/denom
         if (result%t_parallel(i) <= 0.0_rp) error stop 'calibrate_longitudinal_response: nonpositive longitudinal relaxation time'
         result%gamma_parallel(i) = 1.0_rp/result%t_parallel(i)
         allocate(model(nfit))
         do iw = result%fit_first, result%fit_last
            model(iw-result%fit_first+1) = chi_static(i, i)/cmplx(1.0_rp, omega(iw)*result%t_parallel(i), rp)
         end do
         residual_scale = max(1.0_rp, abs(chi_static(i, i)))
         result%fit_residual(i) = sqrt(sum(abs(chi_dynamic(i, i, result%fit_first:result%fit_last) - model)**2)/real(nfit, rp))/residual_scale
         deallocate(model)
         iw = minloc(abs(omega), dim=1)
         result%dynamic_static_error(i) = abs(chi_dynamic(i, i, iw) - chi_static(i, i))/max(1.0_rp, abs(chi_static(i, i)))
      end do
      result%static_agreement_passed = all(result%dynamic_static_error <= options%static_agreement_tolerance)
      if (.not. result%static_agreement_passed) error stop 'calibrate_longitudinal_response: dynamic/static acceptance tolerance failed'
   end subroutine calibrate_longitudinal_response

   !> Matrix form of U_parallel = chi_KS(0)^-1 - chi_parallel(0)^-1.
   !> Keeping this separate makes the static finite-field calibration directly
   !> usable by the existing generalized Dyson engine before a fit is requested.
   subroutine build_longitudinal_kernel(chi_ks_static, chi_static, u_parallel)
      complex(rp), intent(in) :: chi_ks_static(:, :), chi_static(:, :)
      complex(rp), intent(out) :: u_parallel(:, :)
      complex(rp), allocatable :: inv_ks(:, :), inv_static(:, :)
      integer :: nsite

      nsite = size(chi_ks_static, 1)
      if (nsite < 1 .or. size(chi_ks_static, 2) /= nsite .or. any(shape(chi_static) /= [nsite, nsite]) .or. &
          any(shape(u_parallel) /= [nsite, nsite])) error stop 'build_longitudinal_kernel: incompatible matrices'
      allocate(inv_ks(nsite, nsite), inv_static(nsite, nsite))
      call invert_complex_matrix(chi_ks_static, inv_ks)
      call invert_complex_matrix(chi_static, inv_static)
      u_parallel = inv_ks - inv_static
   end subroutine build_longitudinal_kernel

   subroutine write_longitudinal_report(filename, omega, eta, static_result, result)
      character(len=*), intent(in) :: filename
      real(rp), intent(in) :: omega(:), eta
      type(tddft_longitudinal_static_result), intent(in) :: static_result
      type(tddft_longitudinal_result), intent(in) :: result
      integer :: unit, ios, i, j

      open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
      if (ios /= 0) error stop 'write_longitudinal_report: cannot open output file'
      write(unit, '(a)') '# TDDFT longitudinal response: retarded denominator omega+DeltaE+i*eta'
      write(unit, '(a)') '# Gamma_parallel=1/T_parallel is a microscopic inverse relaxation time, NOT an LLB alpha_parallel.'
      write(unit, '(a,es24.16)') '# eta_Ry = ', eta
      write(unit, '(a,2(1x,i0))') '# fit_frequency_indices = ', result%fit_first, result%fit_last
      write(unit, '(a,2(1x,es24.16))') '# fit_frequency_range_Ry = ', omega(result%fit_first), omega(result%fit_last)
      write(unit, '(a,l1)') '# static_field_linearity_passed = ', static_result%linearity_passed
      write(unit, '(a,l1)') '# dynamic_static_agreement_passed = ', result%static_agreement_passed
      write(unit, '(a)') '# site m0 chi_parallel_0 U_parallel_diagonal T_parallel_Ry^-1 Gamma_parallel_Ry fit_residual dynamic_static_relative_error linearity_relative_error'
      do i = 1, size(result%m0)
         write(unit, '(i0,8(1x,es24.16))') i, result%m0(i), real(result%chi_static(i, i), rp), &
            real(result%u_parallel(i, i), rp), result%t_parallel(i), result%gamma_parallel(i), result%fit_residual(i), &
            result%dynamic_static_error(i), static_result%linearity_error(i)
      end do
      write(unit, '(a)') '# U_parallel_matrix: i j Re Im'
      do j = 1, size(result%u_parallel, 2)
         do i = 1, size(result%u_parallel, 1)
            write(unit, '(2(i0,1x),2(es24.16,1x))') i, j, real(result%u_parallel(i, j), rp), aimag(result%u_parallel(i, j))
         end do
      end do
      close(unit)
   end subroutine write_longitudinal_report

   subroutine invert_complex_matrix(matrix, inverse)
      complex(rp), intent(in) :: matrix(:, :)
      complex(rp), intent(out) :: inverse(:, :)
      complex(rp), allocatable :: system(:, :)
      integer, allocatable :: ipiv(:)
      integer :: n, i, info

      n = size(matrix, 1)
      if (n < 1 .or. size(matrix, 2) /= n .or. any(shape(inverse) /= [n, n])) then
         error stop 'invert_complex_matrix: incompatible matrix dimensions'
      end if
      allocate(system(n, n), ipiv(n))
      system = matrix; inverse = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, n
         inverse(i, i) = cmplx(1.0_rp, 0.0_rp, rp)
      end do
      call zgesv(n, n, system, n, ipiv, inverse, n, info)
      if (info /= 0) error stop 'invert_complex_matrix: singular susceptibility matrix'
   end subroutine invert_complex_matrix

end module tddft_longitudinal_mod
