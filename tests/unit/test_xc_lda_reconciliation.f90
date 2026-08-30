! XCR-01: fixed-density reconciliation of the legacy and libXC LDA paths.
!
! This executable deliberately calls only the two production XC kernels.  It
! does not enter the SCF loop, construct a potential, or determine a Fermi
! level.  The first command-line argument, when present, is an output
! directory.  Three CSV files are written there:
!
!   xc_lda_reconciliation.csv
!   xc_lda_stoner_response.csv
!   xc_lda_exchange_validation.csv
!
! The source of truth for the density convention is kept here next to the
! driver: densities are electrons/bohr^3, rs is in bohr, and the libXC
! Hartree results are converted by the production wrapper to Rydberg.
program test_xc_lda_reconciliation
   use precision_mod, only: rp
   use control_mod, only: control
   use logger_mod, only: g_logger
   use xc_mod, only: xc
   implicit none

   integer, parameter :: n_rs = 9
   integer, parameter :: n_zeta = 12
   integer, parameter :: n_delta = 5
   real(rp), parameter :: pi = 4.0_rp*atan(1.0_rp)
   real(rp), parameter :: rs_grid(n_rs) = [ &
      0.5_rp, 1.0_rp, 1.5_rp, 2.0_rp, 2.5_rp, 3.0_rp, 4.0_rp, 5.0_rp, 6.0_rp]
   real(rp), parameter :: zeta_grid(n_zeta) = [ &
      0.0_rp, 1.0e-4_rp, 1.0e-3_rp, 1.0e-2_rp, 0.05_rp, 0.1_rp, &
      0.25_rp, 0.5_rp, 0.8_rp, 0.95_rp, 0.99_rp, 1.0_rp]
   real(rp), parameter :: response_delta(n_delta) = [ &
      1.0e-2_rp, 1.0e-3_rp, 1.0e-4_rp, 1.0e-5_rp, 1.0e-6_rp]
   real(rp), parameter :: radius = 1.0_rp
   real(rp), parameter :: zero_gradient(2) = 0.0_rp
   real(rp), parameter :: zero_laplacian(2) = 0.0_rp
   character(len=*), parameter :: status_regular = 'REGULAR'
   character(len=*), parameter :: status_full_polarization = 'FULLY_POLARIZED_LEGACY_GUARD'

   type(control) :: ctl
   type(xc) :: legacy_bh, libxc_vbh
   type(xc) :: legacy_pw, libxc_pw
   type(xc) :: legacy_vwn, libxc_vwn
   type(xc) :: libxc_exchange
   character(len=1024) :: output_dir, path
   integer :: ios, unit_data, unit_response, unit_exchange, unit_gradient
   integer :: irs, izeta, idelta
   logical :: failed
   real(rp) :: max_exchange_abs_error, max_exchange_rel_error

   real(rp) :: n, rho_up, rho_down
   real(rp) :: legacy_exc, legacy_vup, legacy_vdn
   real(rp) :: libxc_exc, libxc_vup, libxc_vdn
   character(len=64) :: status, classification

   call get_command_argument(1, output_dir)
   if (len_trim(output_dir) == 0) output_dir = '.'

   call g_logger%init()
   call ctl%restore_to_default()
   ctl%nsp = 2

   ctl%txc = 1
   legacy_bh = xc(ctl)
   ctl%txc = 101
   libxc_vbh = xc(ctl)
   ctl%txc = 5
   legacy_pw = xc(ctl)
   ctl%txc = 105
   libxc_pw = xc(ctl)
   ctl%txc = 4
   legacy_vwn = xc(ctl)
   ctl%txc = 106
   libxc_vwn = xc(ctl)
   ctl%txc = 1001
   libxc_exchange = xc(ctl)

   path = trim(output_dir)//'/xc_lda_reconciliation.csv'
   open (newunit=unit_data, file=trim(path), status='replace', action='write', iostat=ios)
   if (ios /= 0) error stop 'Could not open xc_lda_reconciliation.csv'
   write (unit_data, '(a)') &
      'functional_pair,rs,zeta,rho_up,rho_down,exc_legacy,exc_libxc,'// &
      'vup_legacy,vup_libxc,vdn_legacy,vdn_libxc,v0_legacy,v0_libxc,'// &
      'bxc_legacy,bxc_libxc,delta_vxc_legacy,delta_vxc_libxc,'// &
      'exc_abs_difference,exc_relative_difference,vup_abs_difference,'// &
      'vup_relative_difference,vdn_abs_difference,vdn_relative_difference,'// &
      'bxc_abs_difference,bxc_relative_difference,relative_difference,'// &
      'evaluation_status,classification'

   do irs = 1, n_rs
      n = density_from_rs(rs_grid(irs))
      do izeta = 1, n_zeta
         rho_up = 0.5_rp*n*(1.0_rp + zeta_grid(izeta))
         rho_down = 0.5_rp*n*(1.0_rp - zeta_grid(izeta))
         if (zeta_grid(izeta) == 1.0_rp) then
            status = status_full_polarization
         else
            status = status_regular
         endif

         call evaluate_pair(legacy_bh, libxc_vbh, rho_down, rho_up, &
            legacy_exc, legacy_vup, legacy_vdn, libxc_exc, libxc_vup, libxc_vdn)
         classification = classification_for_pair('BH_VBH', status, &
            max(abs(legacy_vup - libxc_vup), abs(legacy_vdn - libxc_vdn)))
         call write_comparison_row(unit_data, 'BH_VBH', rs_grid(irs), zeta_grid(izeta), &
            rho_up, rho_down, legacy_exc, libxc_exc, legacy_vup, libxc_vup, &
            legacy_vdn, libxc_vdn, status, classification)

         call evaluate_pair(legacy_pw, libxc_pw, rho_down, rho_up, &
            legacy_exc, legacy_vup, legacy_vdn, libxc_exc, libxc_vup, libxc_vdn)
         classification = classification_for_pair('PW92_PW', status, &
            max(abs(legacy_vup - libxc_vup), abs(legacy_vdn - libxc_vdn)))
         call write_comparison_row(unit_data, 'PW92_PW', rs_grid(irs), zeta_grid(izeta), &
            rho_up, rho_down, legacy_exc, libxc_exc, legacy_vup, libxc_vup, &
            legacy_vdn, libxc_vdn, status, classification)

         call evaluate_pair(legacy_vwn, libxc_vwn, rho_down, rho_up, &
            legacy_exc, legacy_vup, legacy_vdn, libxc_exc, libxc_vup, libxc_vdn)
         classification = classification_for_pair('VWN_VWN', status, &
            max(abs(legacy_vup - libxc_vup), abs(legacy_vdn - libxc_vdn)))
         call write_comparison_row(unit_data, 'VWN_VWN', rs_grid(irs), zeta_grid(izeta), &
            rho_up, rho_down, legacy_exc, libxc_exc, legacy_vup, libxc_vup, &
            legacy_vdn, libxc_vdn, status, classification)
      enddo
   enddo
   close (unit_data)

   ! Audit the reported legacy potential against a finite difference of its
   ! own energy density.  This localizes a possible common charge-channel
   ! offset without modifying the legacy expression.
   path = trim(output_dir)//'/xc_lda_energy_gradient.csv'
   open (newunit=unit_gradient, file=trim(path), status='replace', action='write', iostat=ios)
   if (ios /= 0) error stop 'Could not open xc_lda_energy_gradient.csv'
   write (unit_gradient, '(a)') &
      'functional_pair,rs,zeta,step,rho_up,rho_down,vup_reported,vup_from_energy,'// &
      'vdn_reported,vdn_from_energy,vup_residual,vdn_residual,max_residual'
   do irs = 1, n_rs
      n = density_from_rs(rs_grid(irs))
      do izeta = 1, n_zeta - 1
         rho_up = 0.5_rp*n*(1.0_rp + zeta_grid(izeta))
         rho_down = 0.5_rp*n*(1.0_rp - zeta_grid(izeta))
         call write_energy_gradient_row(unit_gradient, 'BH_VBH', rs_grid(irs), zeta_grid(izeta), &
            rho_down, rho_up, legacy_bh)
         call write_energy_gradient_row(unit_gradient, 'PW92_PW', rs_grid(irs), zeta_grid(izeta), &
            rho_down, rho_up, legacy_pw)
         call write_energy_gradient_row(unit_gradient, 'VWN_VWN', rs_grid(irs), zeta_grid(izeta), &
            rho_down, rho_up, legacy_vwn)
      enddo
   enddo
   close (unit_gradient)

   path = trim(output_dir)//'/xc_lda_stoner_response.csv'
   open (newunit=unit_response, file=trim(path), status='replace', action='write', iostat=ios)
   if (ios /= 0) error stop 'Could not open xc_lda_stoner_response.csv'
   write (unit_response, '(a)') &
      'functional_pair,rs,delta,i_xc_legacy,i_xc_libxc,absolute_difference,'// &
      'relative_difference,classification'
   do irs = 1, n_rs
      do idelta = 1, n_delta
         call write_response_row(unit_response, 'BH_VBH', rs_grid(irs), response_delta(idelta), &
            legacy_bh, libxc_vbh)
         call write_response_row(unit_response, 'PW92_PW', rs_grid(irs), response_delta(idelta), &
            legacy_pw, libxc_pw)
         call write_response_row(unit_response, 'VWN_VWN', rs_grid(irs), response_delta(idelta), &
            legacy_vwn, libxc_vwn)
      enddo
   enddo
   close (unit_response)

   path = trim(output_dir)//'/xc_lda_exchange_validation.csv'
   open (newunit=unit_exchange, file=trim(path), status='replace', action='write', iostat=ios)
   if (ios /= 0) error stop 'Could not open xc_lda_exchange_validation.csv'
   write (unit_exchange, '(a)') &
      'rs,zeta,rho_up,rho_down,exc_libxc,exc_analytic,vup_libxc,vup_analytic,'// &
      'vdn_libxc,vdn_analytic,exc_abs_error,vup_abs_error,vdn_abs_error,'// &
      'max_abs_error,max_relative_error,status'
   max_exchange_abs_error = 0.0_rp
   max_exchange_rel_error = 0.0_rp
   do irs = 1, n_rs
      do izeta = 1, n_zeta - 1
         n = density_from_rs(rs_grid(irs))
         rho_up = 0.5_rp*n*(1.0_rp + zeta_grid(izeta))
         rho_down = 0.5_rp*n*(1.0_rp - zeta_grid(izeta))
         call evaluate_exchange(libxc_exchange, rho_down, rho_up, &
            libxc_exc, libxc_vup, libxc_vdn)
         call write_exchange_row(unit_exchange, rs_grid(irs), zeta_grid(izeta), &
            rho_up, rho_down, libxc_exc, libxc_vup, libxc_vdn, &
            max_exchange_abs_error, max_exchange_rel_error)
      enddo
   enddo
   close (unit_exchange)

   failed = max_exchange_abs_error > 2.0e-12_rp .or. max_exchange_rel_error > 2.0e-12_rp
   write (*, '(a,es16.8)') 'exchange_max_abs_error=', max_exchange_abs_error
   write (*, '(a,es16.8)') 'exchange_max_relative_error=', max_exchange_rel_error
   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   endif
   write (*, '(a)') 'RESULT: PASS'

contains

   pure function density_from_rs(rs) result(density)
      real(rp), intent(in) :: rs
      real(rp) :: density
      density = 3.0_rp/(4.0_rp*pi*rs**3)
   end function density_from_rs

   subroutine evaluate_pair(legacy, libxc, rho_dn, rho_up, legacy_e, legacy_vup, legacy_vdn, &
      libxc_e, libxc_vup, libxc_vdn)
      type(xc), intent(in) :: legacy, libxc
      real(rp), intent(in) :: rho_dn, rho_up
      real(rp), intent(out) :: legacy_e, legacy_vup, legacy_vdn
      real(rp), intent(out) :: libxc_e, libxc_vup, libxc_vdn
      real(rp) :: legacy_v1, legacy_v2, libxc_v1, libxc_v2
      real(rp) :: rho_total

      rho_total = rho_dn + rho_up
      call legacy%XCPOT(rho_dn, rho_up, rho_total, zero_gradient, zero_laplacian, radius, &
         legacy_v1, legacy_v2, legacy_e)
      call libxc%xcpot_libxc_wrapper(rho_dn, rho_up, rho_total, zero_gradient, zero_laplacian, radius, &
         libxc_v1, libxc_v2, libxc_e)
      legacy_vdn = legacy_v1
      legacy_vup = legacy_v2
      libxc_vdn = libxc_v1
      libxc_vup = libxc_v2
   end subroutine evaluate_pair

   subroutine evaluate_exchange(functional, rho_dn, rho_up, energy, vup, vdn)
      type(xc), intent(in) :: functional
      real(rp), intent(in) :: rho_dn, rho_up
      real(rp), intent(out) :: energy, vup, vdn
      real(rp) :: v1, v2, rho_total

      rho_total = rho_dn + rho_up
      call functional%xcpot_libxc_wrapper(rho_dn, rho_up, rho_total, zero_gradient, zero_laplacian, radius, &
         v1, v2, energy)
      vdn = v1
      vup = v2
   end subroutine evaluate_exchange

   subroutine write_comparison_row(unit, pair_name, rs, zeta, rho_up, rho_down, &
      legacy_e, libxc_e, legacy_vup, libxc_vup, legacy_vdn, libxc_vdn, status, classification)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: pair_name, status, classification
      real(rp), intent(in) :: rs, zeta, rho_up, rho_down
      real(rp), intent(in) :: legacy_e, libxc_e, legacy_vup, libxc_vup, legacy_vdn, libxc_vdn
      real(rp) :: legacy_v0, libxc_v0, legacy_b, libxc_b
      real(rp) :: legacy_delta, libxc_delta
      real(rp) :: e_abs, e_rel, up_abs, up_rel, dn_abs, dn_rel, b_abs, b_rel, rel

      legacy_v0 = 0.5_rp*(legacy_vup + legacy_vdn)
      libxc_v0 = 0.5_rp*(libxc_vup + libxc_vdn)
      legacy_b = 0.5_rp*(legacy_vdn - legacy_vup)
      libxc_b = 0.5_rp*(libxc_vdn - libxc_vup)
      legacy_delta = legacy_vup - legacy_vdn
      libxc_delta = libxc_vup - libxc_vdn
      e_abs = abs(legacy_e - libxc_e)
      up_abs = abs(legacy_vup - libxc_vup)
      dn_abs = abs(legacy_vdn - libxc_vdn)
      b_abs = abs(legacy_b - libxc_b)
      e_rel = relative_difference(legacy_e, libxc_e)
      up_rel = relative_difference(legacy_vup, libxc_vup)
      dn_rel = relative_difference(legacy_vdn, libxc_vdn)
      b_rel = relative_difference(legacy_b, libxc_b)
      rel = max(e_rel, up_rel, dn_rel, b_rel)

      write (unit, '(a,",",25(es24.16e3,","),a,",",a)') trim(pair_name), &
         rs, zeta, rho_up, rho_down, legacy_e, libxc_e, legacy_vup, libxc_vup, &
         legacy_vdn, libxc_vdn, legacy_v0, libxc_v0, legacy_b, libxc_b, &
         legacy_delta, libxc_delta, e_abs, e_rel, up_abs, up_rel, dn_abs, dn_rel, &
         b_abs, b_rel, rel, trim(status), trim(classification)
   end subroutine write_comparison_row

   subroutine write_response_row(unit, pair_name, rs, delta, legacy, libxc)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: pair_name
      real(rp), intent(in) :: rs, delta
      type(xc), intent(in) :: legacy, libxc
      real(rp) :: density, dn_plus, up_plus, dn_minus, up_minus
      real(rp) :: le_plus, lvup_plus, lvdn_plus, le_minus, lvup_minus, lvdn_minus
      real(rp) :: ce_plus, cvup_plus, cvdn_plus, ce_minus, cvup_minus, cvdn_minus
      real(rp) :: i_legacy, i_libxc, abs_diff, rel_diff

      density = density_from_rs(rs)
      up_plus = 0.5_rp*density*(1.0_rp + delta)
      dn_plus = 0.5_rp*density*(1.0_rp - delta)
      up_minus = 0.5_rp*density*(1.0_rp - delta)
      dn_minus = 0.5_rp*density*(1.0_rp + delta)
      call evaluate_pair(legacy, libxc, dn_plus, up_plus, le_plus, lvup_plus, lvdn_plus, &
         ce_plus, cvup_plus, cvdn_plus)
      call evaluate_pair(legacy, libxc, dn_minus, up_minus, le_minus, lvup_minus, lvdn_minus, &
         ce_minus, cvup_minus, cvdn_minus)
      i_legacy = ((lvdn_plus - lvup_plus) - (lvdn_minus - lvup_minus))/(2.0_rp*density*delta)
      i_libxc = ((cvdn_plus - cvup_plus) - (cvdn_minus - cvup_minus))/(2.0_rp*density*delta)
      abs_diff = abs(i_legacy - i_libxc)
      rel_diff = relative_difference(i_legacy, i_libxc)
      write (unit, '(a,",",4(es24.16e3,","),es24.16e3,",",es24.16e3,",",a)') &
         trim(pair_name), rs, delta, i_legacy, i_libxc, abs_diff, rel_diff, &
         trim(classification_for_pair(pair_name, status_regular, 0.0_rp))
   end subroutine write_response_row

   function classification_for_pair(pair_name, status, potential_difference) result(classification)
      character(len=*), intent(in) :: pair_name, status
      real(rp), intent(in) :: potential_difference
      character(len=64) :: classification

      if (status == status_full_polarization) then
         classification = 'LEGACY_KERNEL_DEFECT'
      else if (trim(pair_name) == 'BH_VBH' .and. potential_difference > 1.0e-7_rp) then
         ! The BH energy and magnetic derivative agree closely, but the
         ! legacy charge-average potential develops a reproducible common
         ! offset as polarization grows.  Keep this as a measured defect;
         ! no production formula is changed by this diagnostic.
         classification = 'LEGACY_KERNEL_DEFECT'
      else
         classification = 'EXPECTED_PARAMETRIZATION_DIFFERENCE'
      endif
   end function classification_for_pair

   subroutine write_energy_gradient_row(unit, pair_name, rs, zeta, rho_dn, rho_up, functional)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: pair_name
      real(rp), intent(in) :: rs, zeta, rho_dn, rho_up
      type(xc), intent(in) :: functional
      real(rp) :: step, e0, e_dn_plus, e_dn_minus, e_up_plus, e_up_minus
      real(rp) :: vdn_reported, vup_reported, vup_from_energy, vdn_from_energy
      real(rp) :: vup_residual, vdn_residual, max_residual

      step = 1.0e-6_rp*max(rho_dn, rho_up)
      call evaluate_legacy_energy(functional, rho_dn, rho_up, e0, vdn_reported, vup_reported)
      call evaluate_legacy_energy(functional, rho_dn + step, rho_up, e_dn_plus, vdn_reported, vup_reported)
      call evaluate_legacy_energy(functional, rho_dn - step, rho_up, e_dn_minus, vdn_reported, vup_reported)
      vdn_from_energy = ((rho_dn + step + rho_up)*e_dn_plus - &
         (rho_dn - step + rho_up)*e_dn_minus)/(2.0_rp*step)

      call evaluate_legacy_energy(functional, rho_dn, rho_up + step, e_up_plus, vdn_reported, vup_reported)
      call evaluate_legacy_energy(functional, rho_dn, rho_up - step, e_up_minus, vdn_reported, vup_reported)
      vup_from_energy = ((rho_dn + rho_up + step)*e_up_plus - &
         (rho_dn + rho_up - step)*e_up_minus)/(2.0_rp*step)

      call evaluate_legacy_energy(functional, rho_dn, rho_up, e0, vdn_reported, vup_reported)
      vdn_residual = vdn_reported - vdn_from_energy
      vup_residual = vup_reported - vup_from_energy
      max_residual = max(abs(vup_residual), abs(vdn_residual))
      write (unit, '(a,",",11(es24.16e3,","),es24.16e3)') trim(pair_name), rs, zeta, step, &
         rho_up, rho_dn, vup_reported, vup_from_energy, vdn_reported, vdn_from_energy, &
         vup_residual, vdn_residual, max_residual
   end subroutine write_energy_gradient_row

   subroutine evaluate_legacy_energy(functional, rho_dn, rho_up, energy, vdn, vup)
      type(xc), intent(in) :: functional
      real(rp), intent(in) :: rho_dn, rho_up
      real(rp), intent(out) :: energy, vdn, vup
      real(rp) :: rho_total

      rho_total = rho_dn + rho_up
      call functional%XCPOT(rho_dn, rho_up, rho_total, zero_gradient, zero_laplacian, radius, &
         vdn, vup, energy)
   end subroutine evaluate_legacy_energy

   subroutine write_exchange_row(unit, rs, zeta, rho_up, rho_down, libxc_e, libxc_vup, libxc_vdn, &
      max_abs_error, max_rel_error)
      integer, intent(in) :: unit
      real(rp), intent(in) :: rs, zeta, rho_up, rho_down, libxc_e, libxc_vup, libxc_vdn
      real(rp), intent(inout) :: max_abs_error, max_rel_error
      real(rp) :: analytic_e, analytic_vup, analytic_vdn
      real(rp) :: e_abs, up_abs, dn_abs, row_abs, row_rel

      analytic_e = -1.5_rp*(6.0_rp/pi)**(1.0_rp/3.0_rp)* &
         (rho_up**(4.0_rp/3.0_rp) + rho_down**(4.0_rp/3.0_rp))/(rho_up + rho_down)
      analytic_vup = -2.0_rp*(6.0_rp/pi)**(1.0_rp/3.0_rp)*rho_up**(1.0_rp/3.0_rp)
      analytic_vdn = -2.0_rp*(6.0_rp/pi)**(1.0_rp/3.0_rp)*rho_down**(1.0_rp/3.0_rp)
      e_abs = abs(libxc_e - analytic_e)
      up_abs = abs(libxc_vup - analytic_vup)
      dn_abs = abs(libxc_vdn - analytic_vdn)
      row_abs = max(e_abs, up_abs, dn_abs)
      row_rel = max(relative_difference(libxc_e, analytic_e), &
         relative_difference(libxc_vup, analytic_vup), relative_difference(libxc_vdn, analytic_vdn))
      max_abs_error = max(max_abs_error, row_abs)
      max_rel_error = max(max_rel_error, row_rel)
      write (unit, '(15(es24.16e3,","),a)') rs, zeta, rho_up, rho_down, &
         libxc_e, analytic_e, libxc_vup, analytic_vup, libxc_vdn, analytic_vdn, &
         e_abs, up_abs, dn_abs, row_abs, row_rel, trim(status_regular)
   end subroutine write_exchange_row

   pure function relative_difference(first, second) result(value)
      real(rp), intent(in) :: first, second
      real(rp) :: value
      value = abs(first - second)/max(abs(first), abs(second), 1.0e-30_rp)
   end function relative_difference

end program test_xc_lda_reconciliation
