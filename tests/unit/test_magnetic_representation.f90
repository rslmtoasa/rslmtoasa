program test_magnetic_representation
   use precision_mod, only: rp
   use magnetic_representation_mod, only: periodic_nc, gbt_single_q, explicit_texture, &
      normalize_magnetic_representation, select_endpoint_moments, texture_bulk_reuse_is_valid
   implicit none

   logical :: failed

   failed = .false.
   call test_mode_normalization()
   call test_same_species_site_identity()
   call test_skyrmion_like_texture()
   call test_bulk_reuse_guard()
   call test_gbt_is_not_nc_moment_routing()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine test_mode_normalization()
      if (trim(normalize_magnetic_representation(' GBT_SINGLE_Q')) /= trim(gbt_single_q)) then
         write (*, '(a)') 'FAIL representation normalization'
         failed = .true.
      end if
   end subroutine test_mode_normalization

   subroutine test_same_species_site_identity()
      real(rp) :: site_mom(3, 2), mi(3), mj(3), type_mom(3)
      logical :: valid

      type_mom = [0.0_rp, 0.0_rp, 1.0_rp]
      site_mom(:, 1) = [1.0_rp, 0.0_rp, 0.0_rp]
      site_mom(:, 2) = [0.0_rp, 1.0_rp, 0.0_rp]

      call select_endpoint_moments(explicit_texture, 1, 2, type_mom, type_mom, site_mom, mi, mj, valid)
      if (.not. valid .or. maxval(abs(mi - site_mom(:, 1))) > 0.0_rp .or. &
          maxval(abs(mj - site_mom(:, 2))) > 0.0_rp .or. maxval(abs(mi - mj)) < 1.0_rp) then
         write (*, '(a)') 'FAIL same-species explicit texture lost site identity'
         failed = .true.
      end if

      call select_endpoint_moments(periodic_nc, 1, 2, type_mom, type_mom, site_mom, mi, mj, valid)
      if (.not. valid .or. maxval(abs(mi - type_mom)) > 0.0_rp .or. maxval(abs(mj - type_mom)) > 0.0_rp) then
         write (*, '(a)') 'FAIL periodic_nc did not retain per-type moments'
         failed = .true.
      end if
   end subroutine test_same_species_site_identity

   subroutine test_skyrmion_like_texture()
      real(rp), parameter :: pi_local = acos(-1.0_rp)
      real(rp) :: site_mom(3, 5), mi(3), mj(3), type_mom(3), phi
      logical :: valid
      integer :: i, j

      type_mom = [0.0_rp, 0.0_rp, 1.0_rp]
      site_mom(:, 1) = [0.0_rp, 0.0_rp, -1.0_rp]
      do i = 2, 5
         phi = 0.5_rp*pi_local*real(i - 2, rp)
         site_mom(:, i) = [cos(phi), sin(phi), 0.0_rp]
      end do
      do i = 1, 5
         j = modulo(i, 5) + 1
         call select_endpoint_moments(explicit_texture, i, j, type_mom, type_mom, site_mom, mi, mj, valid)
         if (.not. valid .or. abs(sqrt(sum(mi*mi)) - 1.0_rp) > 2.0e-15_rp .or. &
             abs(sqrt(sum(mj*mj)) - 1.0_rp) > 2.0e-15_rp) then
            write (*, '(a,i0)') 'FAIL skyrmion-like texture endpoint ', i
            failed = .true.
         end if
      end do
   end subroutine test_skyrmion_like_texture

   subroutine test_bulk_reuse_guard()
      integer :: site_types(4), representatives(2)
      real(rp) :: moments(3, 4)

      site_types = [1, 2, 1, 2]
      representatives = [1, 2]
      moments(:, 1) = [0.0_rp, 0.0_rp, 1.0_rp]
      moments(:, 2) = [0.0_rp, 0.0_rp, -1.0_rp]
      moments(:, 3) = moments(:, 1)
      moments(:, 4) = moments(:, 2)
      if (.not. texture_bulk_reuse_is_valid(site_types, representatives, moments, 1.0e-12_rp)) then
         write (*, '(a)') 'FAIL true magnetic-supercell reuse was rejected'
         failed = .true.
      end if
      moments(:, 3) = [1.0_rp, 0.0_rp, 0.0_rp]
      if (texture_bulk_reuse_is_valid(site_types, representatives, moments, 1.0e-12_rp)) then
         write (*, '(a)') 'FAIL invalid per-type explicit-texture reuse was accepted'
         failed = .true.
      end if
   end subroutine test_bulk_reuse_guard

   subroutine test_gbt_is_not_nc_moment_routing()
      real(rp) :: site_mom(3, 1), mi(3), mj(3), type_mom(3)
      logical :: valid

      type_mom = [0.0_rp, 0.0_rp, 1.0_rp]
      site_mom(:, 1) = type_mom
      call select_endpoint_moments(gbt_single_q, 1, 1, type_mom, type_mom, site_mom, mi, mj, valid)
      if (valid) then
         write (*, '(a)') 'FAIL gbt_single_q entered ordinary NC moment provider'
         failed = .true.
      end if
   end subroutine test_gbt_is_not_nc_moment_routing

end program test_magnetic_representation
