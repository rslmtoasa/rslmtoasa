!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!> @brief spin_response C3-C6 on bcc Fe: sum rule, q -> -q, q dependence, U values and the driver's q files.
!> @details Usage, from a scratch copy of tests/spin_response/fe_bcc in which `rslmto.x input_driver.nml`
!>          has run (see fe_scratch.py): test_spin_response_c3 <oracle directory>. The in-process state is
!>          built from input_driver.nml, so the q list, omega grid and eta are the driver's.
!------------------------------------------------------------------------------
program test_spin_response_c3
   use precision_mod, only: rp
   use math_mod, only: pi
   use control_mod, only: control
   use lattice_mod, only: lattice
   use charge_mod, only: charge
   use energy_mod, only: energy
   use hamiltonian_mod, only: hamiltonian
   use reciprocal_mod, only: reciprocal, fermi_dirac_occupation, kB_Ry_per_K
   use spin_response_mod, only: spin_response, d_amplitudes_coefficient, accumulate_chi0
   use spin_response_fe_state
   implicit none

   type(control), target :: control_obj
   type(lattice), target :: lattice_obj
   type(charge), target :: charge_obj
   type(energy), target :: energy_obj
   type(hamiltonian), target :: ham
   type(reciprocal), target :: rec, rec_ref
   type(spin_response) :: sr
   character(len=512) :: dir
   real(rp), allocatable :: bm(:, :, :, :), omega(:)
   complex(rp), allocatable :: c_a(:, :, :), c_b(:, :, :)
   real(rp) :: mref, q(3), eta, w, dev, qs(3, 2)
   integer :: nw, iw, i

   call get_command_argument(1, dir)
   if (len_trim(dir) == 0) error stop 'usage: test_spin_response_c3 <oracle directory>'
   call read_tolerances(trim(dir))
   call require('s1_rel', s1_rel)
   call require('e1_fe_rel', e1_fe_rel)
   call require('q0_pole_abs', q0_pole_abs)
   call build_fe_state(.false., control_obj, lattice_obj, charge_obj, energy_obj, ham, 'input_driver.nml')
   rec = reciprocal(ham)
   rec%fermi_level = energy_obj%fermi
   sr = spin_response(rec)
   call sr%prepare()
   call sr%interaction()
   call reference_band_moments(ham, rec_ref, bm)
   mref = bm(1, 3, 1, 1) - bm(1, 3, 2, 1)

   ! C3.3 (E1): integral of tr L over |omega| <= 5 Ry, eta = 0.02 Ry, d omega = eta/4, Mills amplitudes, vs the
   ! band_moments Mills moment, raw (no tail correction). Catches a wrong normalisation, sign, vertex or orbital sum.
   q = [-1.0_rp, 1.0_rp, 1.0_rp]/12.0_rp
   eta = 0.02_rp
   w = 5.0_rp
   nw = nint(2.0_rp*w/(eta/4.0_rp)) + 1
   allocate (omega(nw), c_a(1, 1, nw))
   omega = [(-w + (iw - 1)*2.0_rp*w/real(nw - 1, rp), iw=1, nw)]
   call sr%bare_response(q, omega, eta, .false., c_a)
   dev = abs(sum(0.5_rp*(-aimag(c_a(1, 1, 1:nw - 1)) - aimag(c_a(1, 1, 2:)))*(omega(2:) - omega(1:nw - 1)))/pi - mref)/abs(mref)
   call report('C3.3 E1 |int tr L - M_mills|/M_mills, 2001 omega, eta 0.02', dev, e1_fe_rel)
   deallocate (omega, c_a)

   ! C3.4 (S1): chi0(q) = chi0(-q), Juelich amplitudes, at a mesh-commensurate q and at an arbitrary q.
   ! Catches k+q folding or inversion handling errors (it cannot see q being ignored: see C3.5).
   nw = 81
   allocate (omega(nw), c_a(1, 1, nw), c_b(1, 1, nw))
   omega = [(-0.5_rp + (iw - 1)*1.0_rp/real(nw - 1, rp), iw=1, nw)]
   qs(:, 1) = q
   qs(:, 2) = [0.0731_rp, -0.0412_rp, 0.1127_rp]
   do i = 1, 2
      call sr%bare_response(qs(:, i), omega, 0.02_rp, .true., c_a)
      call sr%bare_response(-qs(:, i), omega, 0.02_rp, .true., c_b)
      call report('C3.4 S1 max|chi0(q) - chi0(-q)|/max|chi0(-q)|, q #'//achar(48 + i), &
                  maxval(abs(c_a(1, 1, :) - c_b(1, 1, :)))/maxval(abs(c_b(1, 1, :))), s1_rel)
   end do

   ! C3.5: chi0(q) from bare_response vs accumulate_chi0 fed with mesh eigenpairs, k+q taken at the mesh index.
   ! Catches q ignored or misapplied on the run path.
   call mesh_reference(q, c_b)
   call sr%bare_response(q, omega, 0.02_rp, .false., c_a)
   call report('C3.5 chi0(q) vs mesh-index reference, max|a-b|/max|b|', &
               maxval(abs(c_a(1, 1, :) - c_b(1, 1, :)))/maxval(abs(c_b(1, 1, :))), moment_rel)
   call sr%bare_response([0.0_rp, 0.0_rp, 0.0_rp], omega, 0.02_rp, .false., c_a)
   call check('C3.5 non-vacuity: chi0(q) differs from chi0(0) by more than 1e-3 relative', &
              maxval(abs(c_a(1, 1, :) - c_b(1, 1, :)))/maxval(abs(c_b(1, 1, :))) > 1.0e-3_rp)

   ! C4: U values and the Mills residual on Fe, printed, not graded.
   write (*, '(a,f10.7,a,f10.7,a)') 'C4 U_Juelich ', sr%uj(1), ' Ry, U_Mills ', sr%um(1), ' Ry'
   write (*, '(a,f10.7,a,f10.7,a)') 'C4 Mills residual delta ', sr%delta_mills(1), ', implied gap delta*Delta_d ', &
      sr%delta_mills(1)*sr%gap_d(1), ' Ry'

   call driver_files()
   call finish()

contains

   !> chi0(q, omega) at eta = 0.02 from the pre-campaign mesh eigenpairs (rec_ref), Mills amplitudes, no spin swap.
   subroutine mesh_reference(qq, chi0)
      real(rp), intent(in) :: qq(3)
      complex(rp), intent(out) :: chi0(:, :, :)

      integer :: nk, ik, n, m
      real(rp), allocatable :: e(:, :), f(:, :), e_kq(:, :), f_kq(:, :)
      complex(rp), allocatable :: a(:, :, :, :, :), a_kq(:, :, :, :, :)
      real(rp) :: kt

      kt = rec_ref%temperature*kB_Ry_per_K
      nk = rec_ref%k_workset%nk_local
      call check('C3.5 reference assumes no spin relabelling', .not. sr%spin_swapped)
      e = transpose(rec_ref%eigenvalues)
      call d_amplitudes_coefficient(rec_ref%eigenvectors, 1, a)
      allocate (f(nk, size(e, 2)), e_kq(nk, size(e, 2)), f_kq(nk, size(e, 2)), a_kq(nk, size(e, 2), size(a, 3), 1, 2))
      do ik = 1, nk
         do n = 1, size(e, 2)
            f(ik, n) = fermi_dirac_occupation(e(ik, n), rec_ref%fermi_level, kt)
         end do
      end do
      do ik = 1, nk
         m = mesh_index(rec_ref%k_workset%points, rec_ref%k_workset%points(:, ik) + qq)
         e_kq(ik, :) = e(m, :)
         f_kq(ik, :) = f(m, :)
         a_kq(ik, :, :, :, :) = a(m, :, :, :, :)
      end do
      chi0 = (0.0_rp, 0.0_rp)
      call accumulate_chi0(kt, rec_ref%k_workset%weights, e, f, a, e_kq, f_kq, a_kq, omega, 0.02_rp, chi0)
   end subroutine mesh_reference

   integer function mesh_index(points, p) result(m)
      real(rp), intent(in) :: points(:, :), p(3)

      real(rp) :: d(3)

      do m = 1, size(points, 2)
         d = points(:, m) - p
         if (all(abs(d - nint(d)) < 1.0e-8_rp)) return
      end do
      error stop 'mesh_index: k+q is not a mesh point'
   end function mesh_index

   !> C5/C6: the driver's q files and dispersion file against the in-process state.
   subroutine driver_files()
      integer :: u, iq, ios, nq
      real(rp) :: row(7), val(2), ref_dev
      character(len=24) :: tpeak, tcross
      character(len=256) :: fname
      logical :: found
      complex(rp), allocatable :: chi0(:, :, :), from_file(:)
      real(rp), allocatable :: om(:)

      nq = size(sr%q_direct, 2)
      open (newunit=u, file='spin_response_dispersion.dat', status='old', action='read', iostat=ios)
      call check('driver: dispersion file exists', ios == 0)
      read (u, *)
      do iq = 1, nq
         write (fname, '(a,i3.3,a)') 'spin_response_q', iq, '.dat'
         inquire (file=trim(fname), exist=found)
         call check('driver: '//trim(fname)//' exists', found)
         read (u, *) row, tpeak, tcross
         write (*, '(a,i2,a,f9.5,a,a,a,a,a)') 'dispersion q #', iq, ': |q| (1/A) ', row(7), ', peak (Ry) ', trim(tpeak), &
            ', crossing (Ry) ', trim(tcross), '   (printed, not graded)'
         if (iq == 1) then
            ! Static/dynamic consistency check, not physics: U_Juelich is defined by chi0(0,0) U = -1 at omega = 0, eta = 0,
            ! so the q = 0 pole must sit at omega = 0 up to eta; this compares the static route with the dynamic one.
            call check('C6 q = 0 peak and crossing found', trim(tpeak) /= 'n/a' .and. trim(tcross) /= 'n/a')
            read (tpeak, *) val(1)
            read (tcross, *) val(2)
            call report('C6 q = 0 |peak| (Ry), static/dynamic consistency', abs(val(1)), q0_pole_abs)
            call report('C6 q = 0 |crossing| (Ry), static/dynamic consistency', abs(val(2)), q0_pole_abs)
         end if
      end do
      close (u)

      ! Driver's tr chi0 for the second q vs the in-process run path (same q, omega, eta, method).
      open (newunit=u, file='spin_response_q002.dat', status='old', action='read')
      read (u, *)
      allocate (om(size(sr%omega)), from_file(size(sr%omega)), chi0(1, 1, size(sr%omega)))
      do iq = 1, size(sr%omega)
         read (u, *) om(iq), row(1:2)
         from_file(iq) = cmplx(row(1), row(2), rp)
      end do
      close (u)
      call sr%bare_response(sr%q_direct(:, 2), sr%omega, sr%eta, trim(sr%method) == 'juelich', chi0)
      ref_dev = maxval(abs(from_file - chi0(1, 1, :)))/maxval(abs(chi0(1, 1, :)))
      call report('driver: tr chi0 file vs in-process, q #2, max|a-b|/max|b|', ref_dev, moment_rel)
   end subroutine driver_files

end program test_spin_response_c3
