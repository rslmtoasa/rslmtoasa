!------------------------------------------------------------------------------
! RSGF-00 native-GF endpoint augmentation oracle
!------------------------------------------------------------------------------
!
! The test compares the response-owned four-branch adapter with two
! independent constructions on a finite orthogonal LMTO-like model:
! dense A (zI-H)^(-1) A^dagger and a Lehmann sum built from H eigenvectors.
! It also keeps three deliberately incomplete adapters as negative controls.
!------------------------------------------------------------------------------
program test_lr_gf_endpoint_augmentation
   use precision_mod, only: rp
   use dyson_kernel_mod, only: dyson_kspace_inverse
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, lmto_orbital_l
   use lr_gf_endpoint_augmentation_mod, only: lr_gf_augmented_block, lr_gf_endpoint_capabilities, &
      lr_gf_dense_hamiltonian_provider, augment_lr_gf_endpoint, augment_lr_gf_endpoint_pair
   use response_angular_basis_mod, only: response_harmonic
   implicit none

   integer, parameter :: lmax = 1, norb = (lmax + 1)**2, nlocal = 2*norb
   integer, parameter :: nsite = 2, nmat = nsite*nlocal, npoint = 7, nz = 4
   real(rp), parameter :: tolerance = 2.0e-11_rp, negative_floor = 1.0e-6_rp
   complex(rp), parameter :: z_values(nz) = [ &
      cmplx(-0.47_rp, 0.19_rp, rp), cmplx(0.13_rp, 0.23_rp, rp), &
      cmplx(0.81_rp, 0.17_rp, rp), cmplx(1.37_rp, 0.31_rp, rp)]
   type(lmto_radial_basis) :: radial(nsite)
   type(lr_gf_endpoint_capabilities) :: capabilities
   type(lr_gf_augmented_block) :: endpoint_off, endpoint_on, endpoint_provider, endpoint_pair_forward, endpoint_pair_reverse
   type(lr_gf_dense_hamiltonian_provider) :: provider
   complex(rp) :: hgamma(nmat, nmat), h_eff(nmat, nmat), g_full(nmat, nmat)
   complex(rp) :: sigma_zero(nmat, nmat), g_ab(nlocal, nlocal), g_ba(nlocal, nlocal), g_aa(nlocal, nlocal)
   complex(rp) :: a_left(2, nmat), a_right(2, nmat), reference_off(2, 2), &
      reference_on(2, 2), candidate(2, 2)
   complex(rp) :: branches(2, 2, 4), wrong_contacts(2, 2), wrong_hgh(2, 2), wrong_phidot(2, 2)
   complex(rp) :: eigenvectors(nmat, nmat), spectral_reference(2, 2)
   real(rp) :: eigenvalues(nmat)
   integer :: iz, info
   logical :: failed

   failed = .false.
   call setup_radial_bases(radial)
   call setup_hamiltonian(radial, hgamma, h_eff)
   sigma_zero = cmplx(0.0_rp, 0.0_rp, rp)
   call provider%initialize(h_eff, nlocal)

   do iz = 1, nz
      call dyson_kspace_inverse(h_eff, z_values(iz), sigma_zero, g_full)
      g_ab = g_full(1:nlocal, nlocal + 1:nmat)
      g_ba = g_full(nlocal + 1:nmat, 1:nlocal)
      g_aa = g_full(1:nlocal, 1:nlocal)

      call augment_lr_gf_endpoint(radial(1), radial(2), 1, 2, z_values(iz), g_ab, endpoint_off, capabilities, &
         hgamma(1:nlocal, nlocal + 1:nmat))
      call augment_lr_gf_endpoint(radial(1), radial(1), 1, 1, z_values(iz), g_aa, endpoint_on, capabilities, &
         hgamma(1:nlocal, 1:nlocal))
      call augment_lr_gf_endpoint(radial(1), radial(2), 1, 2, z_values(iz), g_ab, endpoint_provider, capabilities, &
         hamiltonian_provider=provider)
      call augment_lr_gf_endpoint_pair(radial(1), radial(2), 1, 2, z_values(iz), g_ab, g_ba, endpoint_pair_forward, &
         endpoint_pair_reverse, capabilities, hgamma(1:nlocal, nlocal + 1:nmat), &
         hgamma(nlocal + 1:nmat, 1:nlocal))
      call report('forward/reverse endpoint pair contract', &
         maxval(abs(endpoint_pair_forward%h_g_h - endpoint_off%h_g_h)), tolerance, failed)

      call endpoint_maps(radial(1), hgamma(1:nlocal, :), 1, 4, 0.71_rp, 0.39_rp, a_left)
      call endpoint_maps(radial(2), hgamma(nlocal + 1:nmat, :), nlocal + 1, 5, 1.17_rp, 2.01_rp, a_right)
      reference_off = matmul(a_left, matmul(g_full, conjg(transpose(a_right))))
      call endpoint_off%evaluate_point(4, 0.71_rp, 0.39_rp, 5, 1.17_rp, 2.01_rp, candidate)
      call report('dense offsite augmented-GF equivalence', maxval(abs(candidate - reference_off)), tolerance, failed)
      call endpoint_provider%evaluate_point(4, 0.71_rp, 0.39_rp, 5, 1.17_rp, 2.01_rp, candidate)
      call report('effective-action provider equivalence', maxval(abs(candidate - reference_off)), tolerance, failed)

      call endpoint_maps(radial(1), hgamma(1:nlocal, :), 1, 4, 0.71_rp, 0.39_rp, a_left)
      call endpoint_maps(radial(1), hgamma(1:nlocal, :), 1, 4, 0.71_rp, 0.39_rp, a_right)
      reference_on = matmul(a_left, matmul(g_full, conjg(transpose(a_right))))
      call endpoint_on%evaluate_point(4, 0.71_rp, 0.39_rp, 4, 0.71_rp, 0.39_rp, candidate)
      call report('dense onsite augmented-GF equivalence', maxval(abs(candidate - reference_on)), tolerance, failed)

      call endpoint_off%branch_at_point(4, 0.71_rp, 0.39_rp, 5, 1.17_rp, 2.01_rp, branches)
      wrong_phidot = branches(:, :, 1)
      call report('omitting phidot augmentation fails', maxval(abs(wrong_phidot - reference_off)), negative_floor, failed, .true.)

      call endpoint_on%branch_at_point(4, 0.71_rp, 0.39_rp, 4, 0.71_rp, 0.39_rp, branches)
      wrong_contacts = onsite_without_contacts(endpoint_on, 4, 0.71_rp, 0.39_rp)
      call report('omitting onsite -I contacts fails', maxval(abs(wrong_contacts - reference_on)), negative_floor, failed, .true.)

      wrong_hgh = branch_without_hgamma(endpoint_off, 4, 0.71_rp, 0.39_rp, 5, 1.17_rp, 2.01_rp)
      call report('omitting offsite -h_ab fails', maxval(abs(wrong_hgh - reference_off)), negative_floor, failed, .true.)

      call hermitian_eig(h_eff, eigenvalues, eigenvectors, info)
      if (info /= 0) then
         write (*, '(a,i0)') 'zheev info = ', info
         failed = .true.
      else
         spectral_reference = spectral_endpoint(eigenvalues, eigenvectors, z_values(iz), a_left, a_right)
         call report('Lehmann augmented-GF equivalence', maxval(abs(candidate - spectral_reference)), tolerance, failed)
      end if
   end do

   if (failed) then
      write (*, '(a)') 'UnitLrGfEndpointAugmentation: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'UnitLrGfEndpointAugmentation: PASS (dense, spectral, provider, negative controls)'

contains

   subroutine setup_radial_bases(bases)
      type(lmto_radial_basis), intent(out) :: bases(:)
      integer :: isite, ir, iorb, ispin, l
      real(rp) :: radius

      do isite = 1, size(bases)
         call bases(isite)%initialize(npoint, lmax, 2)
         bases(isite)%rofi(1) = 0.0_rp
         do ir = 2, npoint
            bases(isite)%rofi(ir) = 0.09_rp*real(ir - 1, rp) + 0.013_rp*real(isite, rp)
         end do
         bases(isite)%mesh_a = 0.03_rp
         bases(isite)%mesh_b = 0.10_rp
         bases(isite)%channel_present = .true.
         do ispin = 1, 2
            do l = 0, lmax
               bases(isite)%enu_work(l + 1, ispin) = -0.31_rp + 0.19_rp*real(isite, rp) + &
                  0.07_rp*real(l, rp) + 0.043_rp*real(ispin, rp)
               do ir = 1, npoint
                  radius = bases(isite)%rofi(ir)
                  ! Regular U_l numerators, with nonzero first energy response.
                  bases(isite)%phi_large(ir, l + 1, ispin) = (0.82_rp + 0.03_rp*real(isite, rp) + &
                     0.02_rp*real(ispin, rp))*(radius**(l + 1))
                  bases(isite)%phidot_large(ir, l + 1, ispin) = (0.11_rp + 0.017_rp*real(l, rp))* &
                     (radius**(l + 1)) + 0.008_rp*radius**(l + 2)
               end do
            end do
         end do
      end do
   end subroutine setup_radial_bases

   subroutine setup_hamiltonian(bases, hgamma_out, h_eff_out)
      type(lmto_radial_basis), intent(in) :: bases(:)
      complex(rp), intent(out) :: hgamma_out(:, :), h_eff_out(:, :)
      integer :: site_i, i, j, gi, gj, l, ispin
      real(rp) :: enu

      hgamma_out = cmplx(0.0_rp, 0.0_rp, rp)
      do site_i = 1, nsite
         do ispin = 1, 2
            do i = 1, norb
               gi = (site_i - 1)*nlocal + (ispin - 1)*norb + i
               l = lmto_orbital_l(i)
               hgamma_out(gi, gi) = cmplx(0.061_rp + 0.009_rp*real(gi, rp) + &
                  0.003_rp*real(site_i, rp), 0.0_rp, rp)
               do j = i + 1, norb
                  gj = (site_i - 1)*nlocal + (ispin - 1)*norb + j
                  hgamma_out(gi, gj) = cmplx(0.004_rp*real(i + j, rp), &
                     0.002_rp*real(j - i, rp), rp)
                  hgamma_out(gj, gi) = conjg(hgamma_out(gi, gj))
               end do
               enu = bases(site_i)%enu_work(l + 1, ispin)
               h_eff_out(gi, gi) = cmplx(enu, 0.0_rp, rp)
            end do
         end do
      end do
      do i = 1, nlocal
         do j = 1, nlocal
            hgamma_out(i, nlocal + j) = cmplx(0.013_rp*real(i + 2*j, rp), &
               0.0015_rp*real(i - j, rp), rp)
            hgamma_out(nlocal + j, i) = conjg(hgamma_out(i, nlocal + j))
         end do
      end do
      h_eff_out = hgamma_out
      do site_i = 1, nsite
         do ispin = 1, 2
            do i = 1, norb
               gi = (site_i - 1)*nlocal + (ispin - 1)*norb + i
               l = lmto_orbital_l(i)
               h_eff_out(gi, gi) = h_eff_out(gi, gi) + &
                  cmplx(bases(site_i)%enu_work(l + 1, ispin), 0.0_rp, rp)
            end do
         end do
      end do
      ! The assignment above used h_eff_out before its explicit zeroing only
      ! through hgamma_out, so all non-diagonal blocks are retained.
   end subroutine setup_hamiltonian

   subroutine endpoint_maps(radial, hgamma_rows, column_offset, ir, theta, azimuth, map)
      type(lmto_radial_basis), intent(in) :: radial
      complex(rp), intent(in) :: hgamma_rows(:, :)
      integer, intent(in) :: column_offset, ir
      real(rp), intent(in) :: theta, azimuth
      complex(rp), intent(out) :: map(:, :)
      complex(rp) :: phi(2, nlocal), phidot(2, nlocal)
      integer :: iorb, ispin, l, index

      map = cmplx(0.0_rp, 0.0_rp, rp)
      phi = cmplx(0.0_rp, 0.0_rp, rp)
      phidot = cmplx(0.0_rp, 0.0_rp, rp)
      do ispin = 1, 2
         do iorb = 1, norb
            l = lmto_orbital_l(iorb)
            index = (ispin - 1)*norb + iorb
            phi(ispin, index) = cmplx(radial%phi_large(ir, l + 1, ispin)/radial%rofi(ir), 0.0_rp, rp)* &
               response_harmonic(l, iorb - l*l - l - 1, theta, azimuth)
            phidot(ispin, index) = cmplx(radial%phidot_large(ir, l + 1, ispin)/radial%rofi(ir), 0.0_rp, rp)* &
               response_harmonic(l, iorb - l*l - l - 1, theta, azimuth)
         end do
      end do
      map(:, :) = cmplx(0.0_rp, 0.0_rp, rp)
      map(:, 1:size(hgamma_rows, 2)) = matmul(phi, identity_matrix(nlocal, size(hgamma_rows, 2), column_offset)) + &
         matmul(phidot, hgamma_rows)
   end subroutine endpoint_maps

   function identity_matrix(nrow, ncol, offset) result(identity)
      integer, intent(in) :: nrow, ncol, offset
      complex(rp) :: identity(nrow, ncol)
      integer :: i

      identity = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, min(nrow, ncol)
         if (i + offset - 1 <= ncol) identity(i, i + offset - 1) = cmplx(1.0_rp, 0.0_rp, rp)
      end do
   end function identity_matrix

   function onsite_without_contacts(endpoint, ir, theta, azimuth) result(value)
      type(lr_gf_augmented_block), intent(in) :: endpoint
      integer, intent(in) :: ir
      real(rp), intent(in) :: theta, azimuth
      complex(rp) :: value(2, 2), local_branches(2, 2, 4)

      call endpoint%branch_at_point(ir, theta, azimuth, ir, theta, azimuth, local_branches)
      ! Recompute only the two contact-sensitive branches from the stored maps:
      ! hG -> D G and Gh -> G D.  Their difference from the production branch
      ! is precisely the missing onsite identity contribution.
      value = local_branches(:, :, 1) + local_branches(:, :, 2) + local_branches(:, :, 3) + local_branches(:, :, 4)
      value = value + contact_branch(endpoint, ir, theta, azimuth, .true.)
   end function onsite_without_contacts

   function contact_branch(endpoint, ir, theta, azimuth, onsite) result(value)
      type(lr_gf_augmented_block), intent(in) :: endpoint
      integer, intent(in) :: ir
      real(rp), intent(in) :: theta, azimuth
      logical, intent(in) :: onsite
      complex(rp) :: value(2, 2)
      complex(rp) :: left(2, nlocal), right(2, nlocal)

      value = cmplx(0.0_rp, 0.0_rp, rp)
      if (.not. onsite) return
      call point_map(endpoint, ir, theta, azimuth, endpoint%phi_left, left)
      call point_map(endpoint, ir, theta, azimuth, endpoint%phi_right, right)
      ! The negative contact in each single-sided branch is +Phi I Phi^dagger
      ! when reconstructed from the contact-free D*G/G*D branches.
      value = matmul(left, conjg(transpose(right))) + matmul(left, conjg(transpose(right)))
   end function contact_branch

   subroutine point_map(endpoint, ir, theta, azimuth, radial_map, map)
      type(lr_gf_augmented_block), intent(in) :: endpoint
      integer, intent(in) :: ir
      real(rp), intent(in) :: theta, azimuth
      real(rp), intent(in) :: radial_map(:, :, :)
      complex(rp), intent(out) :: map(:, :)
      integer :: iorb, ispin, l, index

      map = cmplx(0.0_rp, 0.0_rp, rp)
      do ispin = 1, 2
         do iorb = 1, norb
            l = lmto_orbital_l(iorb)
            index = (ispin - 1)*norb + iorb
            map(ispin, index) = cmplx(radial_map(ir, iorb, ispin)/endpoint%radius_left(ir), 0.0_rp, rp)* &
               response_harmonic(l, iorb - l*l - l - 1, theta, azimuth)
         end do
      end do
   end subroutine point_map

   function branch_without_hgamma(endpoint, ir_left, theta_left, az_left, ir_right, theta_right, az_right) result(value)
      type(lr_gf_augmented_block), intent(in) :: endpoint
      integer, intent(in) :: ir_left, ir_right
      real(rp), intent(in) :: theta_left, az_left, theta_right, az_right
      complex(rp) :: value(2, 2), branches(2, 2, 4)

      call endpoint%branch_at_point(ir_left, theta_left, az_left, ir_right, theta_right, az_right, branches)
      value = sum(branches, dim=3)
      ! The fourth branch used D G D - delta D, with the real hgamma term
      ! removed, so add the omitted Phidot h_ab Phidot^dagger contribution.
      value = value + omitted_hgamma_branch(endpoint, ir_left, theta_left, az_left, ir_right, theta_right, az_right)
   end function branch_without_hgamma

   function omitted_hgamma_branch(endpoint, ir_left, theta_left, az_left, ir_right, theta_right, az_right) result(value)
      type(lr_gf_augmented_block), intent(in) :: endpoint
      integer, intent(in) :: ir_left, ir_right
      real(rp), intent(in) :: theta_left, az_left, theta_right, az_right
      complex(rp) :: value(2, 2), left(2, nlocal), right(2, nlocal)
      call point_map(endpoint, ir_left, theta_left, az_left, endpoint%phidot_left, left)
      call point_map(endpoint, ir_right, theta_right, az_right, endpoint%phidot_right, right)
      value = matmul(left, matmul(endpoint%hgamma_ab, conjg(transpose(right))))
   end function omitted_hgamma_branch

   function spectral_endpoint(eigenvalues, eigenvectors, z, left_map, right_map) result(value)
      real(rp), intent(in) :: eigenvalues(:)
      complex(rp), intent(in) :: eigenvectors(:, :), z
      complex(rp), intent(in) :: left_map(:, :), right_map(:, :)
      complex(rp) :: value(2, 2), left_state(2), right_state(2)
      integer :: n

      value = cmplx(0.0_rp, 0.0_rp, rp)
      do n = 1, size(eigenvalues)
         left_state = matmul(left_map, eigenvectors(:, n))
         right_state = matmul(right_map, eigenvectors(:, n))
         value = value + matmul(reshape(left_state, [2, 1]), &
            reshape(conjg(right_state), [1, 2]))/(z - eigenvalues(n))
      end do
   end function spectral_endpoint

   subroutine hermitian_eig(matrix_in, eigenvalues, eigenvectors, info)
      complex(rp), intent(in) :: matrix_in(:, :)
      real(rp), intent(out) :: eigenvalues(:)
      complex(rp), intent(out) :: eigenvectors(:, :)
      integer, intent(out) :: info
      complex(rp) :: work_query(1)
      complex(rp), allocatable :: work(:)
      real(rp) :: rwork(max(1, 3*size(eigenvalues) - 2))
      integer :: lwork

      eigenvectors = matrix_in
      call zheev('V', 'U', size(eigenvalues), eigenvectors, size(eigenvalues), eigenvalues, work_query, -1, rwork, info)
      if (info /= 0) return
      lwork = max(1, nint(real(work_query(1), rp)))
      allocate(work(lwork))
      call zheev('V', 'U', size(eigenvalues), eigenvectors, size(eigenvalues), eigenvalues, work, lwork, rwork, info)
      deallocate(work)
   end subroutine hermitian_eig

   subroutine report(label, value, limit, failed, require_greater)
      character(len=*), intent(in) :: label
      real(rp), intent(in) :: value, limit
      logical, intent(inout) :: failed
      logical, intent(in), optional :: require_greater
      logical :: greater

      greater = .false.
      if (present(require_greater)) greater = require_greater
      write (*, '(a,es14.6)') trim(label)//' = ', value
      if ((.not. greater .and. value > limit) .or. (greater .and. value <= limit)) failed = .true.
   end subroutine report

end program test_lr_gf_endpoint_augmentation
