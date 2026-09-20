program test_dresp09s_scalar_relativistic

   use precision_mod, only: rp
   use self_mod, only: legacy_gintsr
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use lr_response_space_mod, only: response_space_layout
   use lr_dresp09s_scalar_relativistic_mod, only: sr_l0_lower_transverse_factor, sr_same_spin_numerator, &
      sr_mixed_transverse_point, sr_l0_source_components, sr_l0_source, sr_l0_density_from_matrix
   use response_angular_basis_mod, only: response_gaunt
   implicit none

   type(lmto_radial_basis) :: radial(1)
   type(response_space_layout) :: space
   complex(rp), allocatable :: field(:), components(:, :, :), components_minus(:, :, :), density_matrix(:, :), density(:, :), source(:, :)
   real(rp), allocatable :: radius(:), potential(:)
   real(rp) :: g(18), gp(18), gpp(18), expected, got, max_same, max_mixed, duality
   real(rp) :: a, b, z, energy, energy_down, r, tmc_up, tmc_down, field_coefficient
   real(rp) :: legacy_sum, recovered_sum, max_legacy, max_hermitian
   complex(rp) :: duality_complex
   integer :: npoint, lmax, ir, l, spin, p, q, iorb, flat, unit, response_flat
   integer :: item_l, item_m, item_site, item_ir, item_channel

   npoint = 9
   lmax = 2
   a = 0.08_rp
   b = 0.7_rp
   z = 26.0_rp
   energy = -0.30_rp
   energy_down = -0.17_rp
   allocate(radius(npoint), potential(npoint))
   do ir = 1, npoint
      radius(ir) = real(ir - 1, rp)*0.12_rp
      potential(ir) = -1.5_rp + 0.04_rp*real(ir - 1, rp)
   end do
   call radial(1)%initialize(npoint, lmax, 2)
   radial(1)%rofi = radius
   radial(1)%mesh_a = a
   radial(1)%mesh_b = b
   radial(1)%nuclear_z = z
   do spin = 1, 2
      do l = 0, lmax
         do ir = 1, npoint
            r = radius(ir)
            g = 0.0_rp
            gp = 0.0_rp
            gpp = 0.0_rp
            g(ir) = (0.82_rp + 0.07_rp*spin + 0.03_rp*l)*max(r, 0.02_rp)**(l + 1)
            gp(ir) = (0.14_rp + 0.02_rp*spin + 0.01_rp*l)*max(r, 0.02_rp)**(l + 1)
            gpp(ir) = (0.02_rp + 0.01_rp*spin)*max(r, 0.02_rp)**(l + 1)
         end do
         ! Capture expects a complete mesh-packed pair.  The fixture uses
         ! nonzero distinct lower numerators and nonzero energy derivatives.
         do ir = 1, npoint
            g(ir + npoint) = (0.11_rp + 0.02_rp*spin + 0.01_rp*l)*max(radius(ir), 0.02_rp)**(l + 1)
            gp(ir + npoint) = (0.03_rp + 0.004_rp*spin)*max(radius(ir), 0.02_rp)**(l + 1)
            gpp(ir + npoint) = (0.006_rp + 0.001_rp*spin)*max(radius(ir), 0.02_rp)**(l + 1)
         end do
         call radial(1)%capture_channel(l, spin, merge(energy, energy_down, spin == 1), radius, &
            potential + 0.02_rp*real(spin - 1, rp), a, b, z, g, gp, gpp, merge(energy, energy_down, spin == 1))
      end do
   end do
   call space%initialize(1, 0, radius, a, b, 1)

   max_same = 0.0_rp
   max_legacy = 0.0_rp
   do spin = 1, 2
      do l = 0, lmax
         do ir = 2, npoint
            do p = 0, 1
         do q = 0, 1
                  tmc_up = radial(1)%tmc(ir, l + 1, spin)
                  expected = sr_endpoint(radial(1), ir, l, spin, p, 1)* &
                     sr_endpoint(radial(1), ir, l, spin, q, 1)*(1.0_rp + real(l*(l + 1), rp)/ &
                     (tmc_up*tmc_up*radius(ir)**2)) + sr_endpoint(radial(1), ir, l, spin, p, 2)* &
                     sr_endpoint(radial(1), ir, l, spin, q, 2)
                  got = sr_same_spin_numerator(radial(1), ir, l, spin, p, q)
                  max_same = max(max_same, abs(got - expected))
               end do
            end do
            recovered_sum = recovered_gintsr_integral(radial(1), l, spin, p, q)
            call legacy_gintsr(legacy_pair(radial(1), l, spin, p), legacy_pair(radial(1), l, spin, q), a, b, npoint, z, &
               radial(1)%enu_radial(l + 1, spin), l, radial(1)%potential(:, spin), radius, legacy_sum)
            max_legacy = max(max_legacy, abs(legacy_sum - recovered_sum))
         end do
      end do
   end do
   if (max_same > 2.0e-14_rp .or. max_legacy > 2.0e-14_rp) error stop 'DRESP-09S same-spin GINTSR algebra failed'

   ! Independent mixed-spin reference: the two endpoint TMC factors must be
   ! retained separately, and the transverse minor term is -1/3.
   max_mixed = 0.0_rp
   do l = 0, lmax
      do ir = 2, npoint
         do p = 0, 1
            do q = 0, 1
               tmc_up = radial(1)%tmc(ir, l + 1, 1)
               tmc_down = radial(1)%tmc(ir, l + 1, 2)
               expected = (sr_endpoint(radial(1), ir, l, 1, p, 1)*sr_endpoint(radial(1), ir, l, 2, q, 1) + &
                  sr_l0_lower_transverse_factor*(sr_endpoint(radial(1), ir, l, 1, p, 2)* &
                  sr_endpoint(radial(1), ir, l, 2, q, 2) + real(l*(l + 1), rp)* &
                  sr_endpoint(radial(1), ir, l, 1, p, 1)*sr_endpoint(radial(1), ir, l, 2, q, 1)/ &
                  (tmc_up*tmc_down*radius(ir)**2)))/radius(ir)**2
               got = sr_mixed_transverse_point(radial(1), ir, l, 1, p, 2, q)
               max_mixed = max(max_mixed, abs(got - expected))
            end do
         end do
      end do
   end do
   if (max_mixed > 2.0e-14_rp) error stop 'DRESP-09S unlike-spin mixed bilinear failed'

   ! Symmetric limit: collapse the two channels and compare to the derived
   ! angularly averaged transverse formula, not to a heuristic GFAC average.
   radial(1)%phi_large(:, :, 2) = radial(1)%phi_large(:, :, 1)
   radial(1)%phi_small(:, :, 2) = radial(1)%phi_small(:, :, 1)
   radial(1)%phidot_large(:, :, 2) = radial(1)%phidot_large(:, :, 1)
   radial(1)%phidot_small(:, :, 2) = radial(1)%phidot_small(:, :, 1)
   radial(1)%phiddot_large(:, :, 2) = radial(1)%phiddot_large(:, :, 1)
   radial(1)%phiddot_small(:, :, 2) = radial(1)%phiddot_small(:, :, 1)
   radial(1)%enu_radial(:, 2) = radial(1)%enu_radial(:, 1)
   radial(1)%enu_work(:, 2) = radial(1)%enu_work(:, 1)
   radial(1)%tmc(:, :, 2) = radial(1)%tmc(:, :, 1)
   radial(1)%gfac(:, :, 2) = radial(1)%gfac(:, :, 1)
   max_mixed = 0.0_rp
   do l = 0, lmax
      do ir = 2, npoint
         got = sr_mixed_transverse_point(radial(1), ir, l, 1, 0, 2, 0)
         expected = (sr_endpoint(radial(1), ir, l, 1, 0, 1)**2 + sr_l0_lower_transverse_factor* &
            (sr_endpoint(radial(1), ir, l, 1, 0, 2)**2 + real(l*(l + 1), rp)* &
            sr_endpoint(radial(1), ir, l, 1, 0, 1)**2/(radial(1)%tmc(ir, l + 1, 1)**2*radius(ir)**2)))/radius(ir)**2
         max_mixed = max(max_mixed, abs(got - expected))
      end do
   end do
   if (max_mixed > 2.0e-14_rp) error stop 'DRESP-09S symmetric mixed-spin limit failed'

   allocate(field(space%ndim), components(2*(lmax + 1)**2, 2*(lmax + 1)**2, 4), &
      components_minus(2*(lmax + 1)**2, 2*(lmax + 1)**2, 4), &
      density_matrix(2*(lmax + 1)**2, 2*(lmax + 1)**2), density(1, npoint), source(2*(lmax + 1)**2, 2*(lmax + 1)**2))
   field_coefficient = 0.7_rp
   field = cmplx(0.0_rp, 0.0_rp, rp)
   do ir = 1, npoint
      call response_flatten_local(1, 0, 0, ir, response_flat)
      field(response_flat) = cmplx(field_coefficient, 0.0_rp, rp)
   end do
   call sr_l0_source_components(space, radial, field, 1, components)
   call sr_l0_source_components(space, radial, field, 2, components_minus)
   max_hermitian = 0.0_rp
   do p = 0, 1
      do q = 0, 1
         max_hermitian = max(max_hermitian, maxval(abs(components_minus(:, :, 1 + q + 2*p) - &
            transpose(conjg(components(:, :, 1 + p + 2*q))))))
      end do
   end do
   if (max_hermitian > 3.0e-13_rp) error stop 'DRESP-09S plus/minus Hermitian relation failed'
   do p = 0, 1
      do q = 0, 1
         call branch_reference(radial(1), space, field(response_flat), p, q, expected)
         got = real(components(1 + (lmax + 1)**2, 1, 1 + p + 2*q), rp)
         ! The selected entry is s(up/down); only its scale is checked here.
         if (abs(got - expected) > 3.0e-13_rp) error stop 'DRESP-09S branch endpoint contraction failed'
      end do
   end do

   density_matrix = cmplx(0.0_rp, 0.0_rp, rp)
   density_matrix(1, 1 + (lmax + 1)**2) = cmplx(0.2_rp, 0.1_rp, rp)
   call sr_l0_density_from_matrix(space, radial, density_matrix, 1, density)
   call sr_l0_source(space, radial, field, 1, cmplx(0.0_rp, 0.0_rp, rp)*density_matrix, source)
   duality_complex = cmplx(0.0_rp, 0.0_rp, rp)
   do ir = 1, npoint
      call response_flatten_local(1, 0, 0, ir, response_flat)
      duality_complex = duality_complex + field(response_flat)*density(1, ir)*space%radial_weights(ir)
   end do
   if (abs(duality_complex - sum(density_matrix*transpose(source))) > 3.0e-12_rp) then
      error stop 'DRESP-09S field/density pairing failed'
   end if

   write (*, '(a,es16.8)') 'DRESP-09S same-spin closure = ', max_same
   write (*, '(a,es16.8)') 'DRESP-09S legacy GINTSR closure = ', max_legacy
   write (*, '(a,es16.8)') 'DRESP-09S mixed-spin closure = ', max_mixed
   write (*, '(a,es16.8)') 'DRESP-09S plus/minus Hermitian closure = ', max_hermitian
   write (*, '(a,es16.8)') 'DRESP-09S field-density duality = ', abs(duality_complex - sum(density_matrix*transpose(source)))

contains

   pure real(rp) function sr_endpoint(radial, ir, l, spin, power, component) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, spin, power, component
      if (component == 1) then
         value = merge(radial%phi_large(ir, l + 1, spin) - radial%enu_work(l + 1, spin)* &
            radial%phidot_large(ir, l + 1, spin), radial%phidot_large(ir, l + 1, spin), power == 0)
      else
         value = merge(radial%phi_small(ir, l + 1, spin) - radial%enu_work(l + 1, spin)* &
            radial%phidot_small(ir, l + 1, spin), radial%phidot_small(ir, l + 1, spin), power == 0)
      end if
   end function sr_endpoint

   subroutine response_flatten_local(site, response_l, response_m, ir, flat)
      integer, intent(in) :: site, response_l, response_m, ir
      integer, intent(out) :: flat
      flat = (((site - 1)*(space%response_lmax + 1)**2 + response_l*response_l + response_l + response_m)* &
         space%npoint + ir - 1)*space%nchannel + 1
   end subroutine response_flatten_local

   subroutine branch_reference(radial, response, coefficient, p, q, value)
      type(lmto_radial_basis), intent(in) :: radial
      type(response_space_layout), intent(in) :: response
      complex(rp), intent(in) :: coefficient
      integer, intent(in) :: p, q
      real(rp), intent(out) :: value
      integer :: ir
      value = 0.0_rp
      do ir = 2, response%npoint
         value = value + response%radial_weights(ir)*real(coefficient, rp)*response_gaunt(0, 0, 0, 0, 0, 0)* &
            sr_mixed_transverse_point(radial, ir, 0, 2, q, 1, p)
      end do
   end subroutine branch_reference

   function legacy_pair(radial, l, spin, power) result(pair)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: l, spin, power
      real(rp) :: pair(radial%npoint, 2)
      integer :: ir

      do ir = 1, radial%npoint
         pair(ir, 1) = sr_endpoint(radial, ir, l, spin, power, 1)
         pair(ir, 2) = sr_endpoint(radial, ir, l, spin, power, 2)
      end do
   end function legacy_pair

   function recovered_gintsr_integral(radial, l, spin, p, q) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: l, spin, p, q
      real(rp) :: value
      integer :: ir

      value = 0.0_rp
      do ir = 2, radial%npoint - 1, 2
         value = value + (radial%rofi(ir) + radial%mesh_b)*sr_same_spin_numerator(radial, ir, l, spin, p, q)
      end do
      value = 2.0_rp*value
      do ir = 3, radial%npoint - 2, 2
         value = value + (radial%rofi(ir) + radial%mesh_b)*sr_same_spin_numerator(radial, ir, l, spin, p, q)
      end do
      value = 2.0_rp*value
      ir = radial%npoint
      value = value + (radial%rofi(ir) + radial%mesh_b)*sr_same_spin_numerator(radial, ir, l, spin, p, q)
      value = value*radial%mesh_a/3.0_rp
   end function recovered_gintsr_integral

end program test_dresp09s_scalar_relativistic
