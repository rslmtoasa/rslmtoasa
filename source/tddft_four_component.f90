!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Full site-resolved charge-spin LR-TDDFT assembly.
!>
!> The response basis is site-major `(n,mx,my,mz)`.  Bare response remains
!> the existing exact eigenpair/vertex summation; this module supplies the
!> unambiguous four-component channel list, rotates each local ALSDA kernel
!> into that global basis, adds a caller-supplied charge Hartree matrix, and
!> reports (rather than forces) no-SOC rigid-spin zero modes.
module tddft_four_component_mod
   use precision_mod, only: rp
   use response_components_mod, only: RESPONSE_CHARGE, RESPONSE_MX, RESPONSE_MY, RESPONSE_MZ
   use response_vertices_mod, only: response_channel, build_site_charge_spin_channels
   use tddft_chi0_mod, only: tddft_chi0_options, tddft_chi0_result, build_chi_ks_from_eigenpairs
   use xc_response_kernel_mod, only: xc_response_kernel_provider, cartesian_transverse_kernel
   implicit none

   private

   type, public :: tddft_four_component_zero_mode_diagnostics
      logical :: applicable = .false.
      logical :: symmetry_broken = .false.
      integer :: number_of_modes = 0
      real(rp), allocatable :: residual(:)
      real(rp), allocatable :: norm(:)
   end type tddft_four_component_zero_mode_diagnostics

   public :: build_four_component_chi_ks
   public :: build_four_component_kernel
   public :: evaluate_four_component_zero_modes
   public :: four_component_index

contains

   pure integer function four_component_index(site, component) result(index)
      integer, intent(in) :: site, component

      if (site < 1 .or. component < RESPONSE_CHARGE .or. component > RESPONSE_MZ) then
         error stop 'four_component_index: invalid site or Cartesian response component'
      end if
      index = 4*(site-1) + component + 1
   end function four_component_index

   !> Build every chi_KS^{mu,nu} through the same generalized spinor vertices
   !> as the validated transverse and longitudinal routes.  This accepts SOC
   !> and non-collinear spinors without a collinearity-specific branch.
   subroutine build_four_component_chi_ks(k_weights, eigenvalues_k, eigenvectors_k, eigenvalues_kq, eigenvectors_kq, &
      site_orbital_counts, omega, options, result)
      real(rp), intent(in) :: k_weights(:), eigenvalues_k(:, :), eigenvalues_kq(:, :), omega(:)
      complex(rp), intent(in) :: eigenvectors_k(:, :, :), eigenvectors_kq(:, :, :)
      integer, intent(in) :: site_orbital_counts(:)
      type(tddft_chi0_options), intent(in) :: options
      type(tddft_chi0_result), intent(out) :: result
      type(response_channel), allocatable :: channels(:)

      call build_site_charge_spin_channels(size(site_orbital_counts), channels)
      call build_chi_ks_from_eigenpairs(k_weights, eigenvalues_k, eigenvectors_k, eigenvalues_kq, eigenvectors_kq, &
         site_orbital_counts, channels, channels, omega, options, result)
   end subroutine build_four_component_chi_ks

   !> Assemble K = f_H + f_xc in global response coordinates.  `coulomb_site`
   !> is the physically projected charge-charge Coulomb kernel for this q;
   !> no spin component is populated by Hartree.  The local ALSDA derivative
   !> convention is (n,m1,m2,m3), with m3 along the local ground-state moment.
   subroutine build_four_component_kernel(provider, coulomb_site, kernel)
      type(xc_response_kernel_provider), intent(in) :: provider
      real(rp), intent(in) :: coulomb_site(:, :)
      complex(rp), intent(out) :: kernel(:, :)
      integer :: nsite, isite, jsite, i0, j0
      logical :: full_response_supported
      character(len=256) :: capability_reason
      real(rp) :: local_kernel(4, 4), rotation(3, 3), transform(4, 4), global_kernel(4, 4)

      if (.not. allocated(provider%site)) error stop 'build_four_component_kernel: XC provider is not initialized'
      nsite = size(provider%site)
      if (nsite < 1 .or. any(shape(coulomb_site) /= [nsite, nsite]) .or. any(shape(kernel) /= [4*nsite, 4*nsite])) then
         error stop 'build_four_component_kernel: incompatible site or response dimensions'
      end if
      call provider%full_response_capability(full_response_supported, capability_reason)
      if (.not. full_response_supported) then
         error stop 'build_four_component_kernel: full response unavailable: '//trim(capability_reason)
      end if
      kernel = cmplx(0.0_rp, 0.0_rp, rp)

      do isite = 1, nsite
         call build_local_alsda_kernel(provider, isite, local_kernel)
         call local_frame(provider%site(isite)%magnetization_direction, rotation)
         transform = 0.0_rp
         transform(1, 1) = 1.0_rp
         transform(2:4, 2:4) = rotation
         global_kernel = matmul(transform, matmul(local_kernel, transpose(transform)))
         i0 = 4*(isite-1)
         kernel(i0+1:i0+4, i0+1:i0+4) = cmplx(global_kernel, 0.0_rp, rp)
      end do
      do jsite = 1, nsite
         j0 = 4*(jsite-1)
         do isite = 1, nsite
            i0 = 4*(isite-1)
            kernel(i0+1, j0+1) = kernel(i0+1, j0+1) + cmplx(coulomb_site(isite, jsite), 0.0_rp, rp)
         end do
      end do
   end subroutine build_four_component_kernel

   !> Diagnose the three global rigid rotations of a no-SOC, field-free
   !> state.  Zero-norm generators (e.g. rotation about a collinear moment)
   !> are omitted.  With SOC or an external field, the physically allowed
   !> anisotropy gap is retained and diagnostics explicitly report that no
   !> zero-mode condition applies.
   subroutine evaluate_four_component_zero_modes(chi_ks_static, kernel, magnetization, has_soc, has_external_field, diagnostics)
      complex(rp), intent(in) :: chi_ks_static(:, :), kernel(:, :)
      real(rp), intent(in) :: magnetization(:, :)
      logical, intent(in) :: has_soc, has_external_field
      type(tddft_four_component_zero_mode_diagnostics), intent(out) :: diagnostics
      complex(rp), allocatable :: xi(:, :), mode(:), residual_vector(:)
      real(rp) :: axis(3), mode_norm
      integer :: nsite, iaxis, isite, nmode

      nsite = size(magnetization, 2)
      if (size(magnetization, 1) /= 3 .or. nsite < 1 .or. any(shape(chi_ks_static) /= [4*nsite, 4*nsite]) .or. &
         any(shape(kernel) /= [4*nsite, 4*nsite])) then
         error stop 'evaluate_four_component_zero_modes: incompatible response dimensions'
      end if
      diagnostics%symmetry_broken = has_soc .or. has_external_field
      if (diagnostics%symmetry_broken) return
      diagnostics%applicable = .true.
      allocate(xi(4*nsite, 4*nsite), mode(4*nsite), residual_vector(4*nsite), diagnostics%residual(3), diagnostics%norm(3))
      xi = matmul(chi_ks_static, kernel)
      diagnostics%residual = 0.0_rp
      diagnostics%norm = 0.0_rp
      nmode = 0
      do iaxis = 1, 3
         axis = 0.0_rp; axis(iaxis) = 1.0_rp
         mode = cmplx(0.0_rp, 0.0_rp, rp)
         do isite = 1, nsite
            mode(4*(isite-1)+2:4*isite) = cmplx(cross(magnetization(:, isite), axis), 0.0_rp, rp)
         end do
         mode_norm = sqrt(sum(abs(mode)**2))
         if (mode_norm <= tiny(1.0_rp)) cycle
         nmode = nmode + 1
         residual_vector = matmul(xi, mode) - mode
         diagnostics%norm(nmode) = mode_norm
         diagnostics%residual(nmode) = sqrt(sum(abs(residual_vector)**2))/mode_norm
      end do
      diagnostics%number_of_modes = nmode
      if (nmode < 3) then
         call shrink_diagnostics(diagnostics, nmode)
      end if
   end subroutine evaluate_four_component_zero_modes

   subroutine build_local_alsda_kernel(provider, isite, local_kernel)
      type(xc_response_kernel_provider), intent(in) :: provider
      integer, intent(in) :: isite
      real(rp), intent(out) :: local_kernel(4, 4)

      local_kernel = 0.0_rp
      local_kernel(1, 1) = provider%site(isite)%dvxc_dn
      local_kernel(1, 4) = provider%site(isite)%dvxc_dm
      local_kernel(4, 1) = provider%site(isite)%dbxc_dn
      local_kernel(4, 4) = provider%site(isite)%dbxc_dm
      ! k_perp_circular is Bxc/(2M) for unhalved sigma+/- vertices.  The
      ! Cartesian x/y block couples directly to sigma_x/y, so it is Bxc/M.
      local_kernel(2, 2) = cartesian_transverse_kernel(provider, isite)
      local_kernel(3, 3) = cartesian_transverse_kernel(provider, isite)
   end subroutine build_local_alsda_kernel

   subroutine local_frame(direction, rotation)
      real(rp), intent(in) :: direction(3)
      real(rp), intent(out) :: rotation(3, 3)
      real(rp) :: e1(3), e2(3), e3(3), reference(3), norm_e3, norm_e1

      norm_e3 = sqrt(sum(direction**2))
      if (norm_e3 <= tiny(1.0_rp)) error stop 'local_frame: zero ground-state magnetization direction'
      e3 = direction/norm_e3
      reference = [1.0_rp, 0.0_rp, 0.0_rp]
      if (abs(dot_product(reference, e3)) > 0.9_rp) reference = [0.0_rp, 1.0_rp, 0.0_rp]
      e1 = reference - dot_product(reference, e3)*e3
      norm_e1 = sqrt(sum(e1**2))
      e1 = e1/norm_e1
      e2 = cross(e3, e1)
      rotation(:, 1) = e1
      rotation(:, 2) = e2
      rotation(:, 3) = e3
   end subroutine local_frame

   pure function cross(left, right) result(value)
      real(rp), intent(in) :: left(3), right(3)
      real(rp) :: value(3)

      value = [left(2)*right(3)-left(3)*right(2), left(3)*right(1)-left(1)*right(3), &
         left(1)*right(2)-left(2)*right(1)]
   end function cross

   subroutine shrink_diagnostics(diagnostics, nmode)
      type(tddft_four_component_zero_mode_diagnostics), intent(inout) :: diagnostics
      integer, intent(in) :: nmode
      real(rp), allocatable :: residual(:), norm(:)

      allocate(residual(nmode), norm(nmode))
      if (nmode > 0) then
         residual = diagnostics%residual(1:nmode)
         norm = diagnostics%norm(1:nmode)
      end if
      call move_alloc(residual, diagnostics%residual)
      call move_alloc(norm, diagnostics%norm)
   end subroutine shrink_diagnostics

end module tddft_four_component_mod
