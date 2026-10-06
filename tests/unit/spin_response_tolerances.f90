!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test helper
!> @brief Every key of &spin_response_tolerances (tests/spin_response/oracles/tolerances.nml), declared once.
!> @details A new key needs one declaration and one namelist entry here and a line in the file. A value of
!>          -1 means the key is absent from the file; `require` stops a test that needs it.
!------------------------------------------------------------------------------
module spin_response_tolerances_mod
   use precision_mod, only: rp
   implicit none

   real(rp) :: a1_rel = -1.0_rp, a2_chi_rel = -1.0_rp, a2_identity_rel = -1.0_rp, static_eta0_rel = -1.0_rp
   real(rp) :: a3_rel = -1.0_rp, window_moment_rel = -1.0_rp, dyson_u0_rel = -1.0_rp
   real(rp) :: pole_peak_rel = -1.0_rp, pole_crossing_rel = -1.0_rp
   real(rp) :: electron_count_abs = -1.0_rp, eigen_contract_abs = -1.0_rp, moment_rel = -1.0_rp
   real(rp) :: s1_rel = -1.0_rp, s2_rel = -1.0_rp, e1_fe_rel = -1.0_rp, q0_pole_abs = -1.0_rp
   real(rp) :: region1_rel = -1.0_rp, region2_rel = -1.0_rp, baseline_repro_rel = -1.0_rp
   namelist /spin_response_tolerances/ a1_rel, a2_chi_rel, a2_identity_rel, static_eta0_rel, a3_rel, &
      window_moment_rel, dyson_u0_rel, pole_peak_rel, pole_crossing_rel, &
      electron_count_abs, eigen_contract_abs, moment_rel, s1_rel, s2_rel, e1_fe_rel, q0_pole_abs, &
      region1_rel, region2_rel, baseline_repro_rel

contains

   subroutine read_tolerance_file(oracle_dir)
      character(*), intent(in) :: oracle_dir

      integer :: u

      open (newunit=u, file=trim(oracle_dir)//'/tolerances.nml', status='old', action='read')
      read (u, nml=spin_response_tolerances)
      close (u)
   end subroutine read_tolerance_file

   subroutine require(key, tol)
      character(*), intent(in) :: key
      real(rp), intent(in) :: tol

      if (.not. tol > 0.0_rp) then
         write (*, '(a)') 'tolerances.nml: key '//key//' missing or <= 0'
         error stop 1
      end if
   end subroutine require

end module spin_response_tolerances_mod
