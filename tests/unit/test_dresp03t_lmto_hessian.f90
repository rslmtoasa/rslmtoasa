!------------------------------------------------------------------------------
! DRESP-03T -- complete local-torque and mixed-contact Hessian audit.
!
! The fixture is made from one underlying LMTO channel set.  The live
! predls-equivalent transformation creates cx/wx/obx/enim channels; the H,
! T, and C paths then consume those channels through the same bond and HOH
! product algebra used by the reciprocal ham_only representation.
!------------------------------------------------------------------------------
program test_dresp03t_lmto_hessian
   use precision_mod, only: rp
   use lr_kl_hessian_mod, only: lmto_live_hamiltonian_fixture, lmto_fixture_init, &
      assemble_lmto_hamiltonian, assemble_lmto_torque, assemble_lmto_mixed_derivative, &
      force_theorem_hessian_from_eigenbasis, mixed_second_difference, grand_potential_from_eigenvalues
   use lmto_magnetic_tangent_mod, only: lmto_transform_potential
   implicit none

   type(lmto_live_hamiltonian_fixture) :: fixture
   real(rp) :: max_first_error, mixed_error, theorem_error, contact_gap
   logical :: failed

   failed = .false.
   call make_legitimate_fixture(fixture)
   call test_first_derivative_oracle(fixture, max_first_error)
   call test_mixed_derivative_oracle(fixture, mixed_error)
   call test_grand_potential_hessian(fixture, theorem_error, contact_gap)
   write (*, '(a,es12.4)') 'DRESP-03T complete first-derivative max error = ', max_first_error
   write (*, '(a,es12.4)') 'DRESP-03T complete mixed-derivative max error = ', mixed_error
   write (*, '(a,es12.4)') 'DRESP-03T grand-potential complete-Hessian error = ', theorem_error
   write (*, '(a,es12.4)') 'DRESP-03T TT-only/contact separation          = ', contact_gap
   if (max_first_error > 3.0e-9_rp) failed = .true.
   if (mixed_error > 3.0e-8_rp) failed = .true.
   if (theorem_error > 3.0e-5_rp) failed = .true.
   if (contact_gap < 1.0e-7_rp) failed = .true.
   call fixture%clear()
   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine make_legitimate_fixture(fix)
      type(lmto_live_hamiltonian_fixture), intent(out) :: fix
      integer, parameter :: nsite = 2, norb = 3, nbond = 4, nchan = 2
      real(rp) :: c(nchan,nchan,nsite), enu(nchan,nchan,nsite), srdel(nchan,nchan,nsite), &
         qpar(nchan,nchan,nsite), center(nchan,nchan), width(nchan,nchan), shifted(nchan,nchan), &
         obar(nchan,nchan), dele(nchan,nchan), qm(nchan)
      integer :: site, orb, l

      call lmto_fixture_init(fix, nsite, norb, nbond, .true.)
      ! All channels below are the underlying LMTO input.  No transformed wx,
      ! cx, or obar values are invented independently for the Hessian test.
      c(:,:,1) = reshape([-0.72_rp,-0.53_rp,-0.24_rp,-0.39_rp,-0.18_rp,-0.07_rp], [nchan,nchan])
      c(:,:,2) = reshape([-0.61_rp,-0.43_rp,-0.28_rp,-0.31_rp,-0.11_rp,-0.02_rp], [nchan,nchan])
      enu(:,:,1) = reshape([-0.21_rp,-0.17_rp,-0.13_rp,-0.16_rp,-0.12_rp,-0.08_rp], [nchan,nchan])
      enu(:,:,2) = reshape([-0.19_rp,-0.15_rp,-0.11_rp,-0.14_rp,-0.10_rp,-0.06_rp], [nchan,nchan])
      srdel(:,:,1) = reshape([0.74_rp,0.69_rp,0.62_rp,0.81_rp,0.73_rp,0.67_rp], [nchan,nchan])
      srdel(:,:,2) = reshape([0.71_rp,0.66_rp,0.59_rp,0.78_rp,0.70_rp,0.64_rp], [nchan,nchan])
      qpar(:,:,1) = reshape([0.11_rp,0.08_rp,0.05_rp,0.13_rp,0.09_rp,0.06_rp], [nchan,nchan])
      qpar(:,:,2) = reshape([0.10_rp,0.07_rp,0.045_rp,0.12_rp,0.085_rp,0.055_rp], [nchan,nchan])
      qm = [0.021_rp, -0.014_rp]
      do site = 1, nsite
         call lmto_transform_potential(c(:,:,site), enu(:,:,site), srdel(:,:,site), qpar(:,:,site), &
            0.17_rp, 1.41_rp, 1.0_rp, qm, center, width, shifted, obar, dele)
         do orb = 1, norb
            l = merge(0, 1, orb == 1)
            ! For HOH, ham0m_nc uses cex0/cex1 (the transformed shifted
            ! centers) in the live bond builder; center_band still enters the
            ! same underlying potential transformation and the e_nu term.
            fix%c0(orb,site) = 0.5_rp*(shifted(l+1,1)+shifted(l+1,2))
            fix%c1(orb,site) = 0.5_rp*(shifted(l+1,1)-shifted(l+1,2))
            fix%wx0(orb,site) = 0.5_rp*(width(l+1,1)+width(l+1,2))
            fix%wx1(orb,site) = 0.5_rp*(width(l+1,1)-width(l+1,2))
            fix%obar0(orb,site) = 0.5_rp*(obar(l+1,1)+obar(l+1,2))
            fix%obar1(orb,site) = 0.5_rp*(obar(l+1,1)-obar(l+1,2))
            ! cx-cex is exactly the live e_nu coefficient assembled by
            ! hamiltonian_build%build_enim.
            fix%enu0(orb,site) = 0.5_rp*((center(l+1,1)-shifted(l+1,1)) + &
               (center(l+1,2)-shifted(l+1,2)))
            fix%enu1(orb,site) = 0.5_rp*((center(l+1,1)-shifted(l+1,1)) - &
               (center(l+1,2)-shifted(l+1,2)))
         end do
      end do
      fix%moments(:,1) = [0.0_rp, 0.0_rp, 1.0_rp]
      fix%moments(:,2) = [sin(0.37_rp), 0.0_rp, cos(0.37_rp)]
      fix%bond_source = [1, 2, 1, 2]
      fix%bond_target = [1, 2, 2, 1]
      fix%onsite = [.true., .true., .false., .false.]
      fix%bond_vector(:,1) = 0.0_rp; fix%bond_vector(:,2) = 0.0_rp
      fix%bond_vector(:,3) = [0.19_rp, -0.07_rp, 0.11_rp]
      fix%bond_vector(:,4) = -fix%bond_vector(:,3)
      fix%hhh(:,:,1) = 0.0_rp; fix%hhh(:,:,2) = 0.0_rp
      fix%hhh(1,1,1) = 0.38_rp; fix%hhh(2,2,1) = -0.24_rp; fix%hhh(3,3,1) = 0.17_rp
      fix%hhh(1,1,2) = -0.31_rp; fix%hhh(2,2,2) = 0.29_rp; fix%hhh(3,3,2) = -0.12_rp
      fix%hhh(:,:,3) = reshape([cmplx(0.18_rp,0.021_rp,rp), cmplx(0.041_rp,-0.013_rp,rp), cmplx(-0.027_rp,0.018_rp,rp), &
         cmplx(-0.033_rp,0.009_rp,rp), cmplx(-0.11_rp,0.017_rp,rp), cmplx(0.052_rp,0.014_rp,rp), &
         cmplx(0.024_rp,-0.016_rp,rp), cmplx(-0.045_rp,0.008_rp,rp), cmplx(0.09_rp,-0.019_rp,rp)], [norb,norb])
      fix%hhh(:,:,4) = transpose(conjg(fix%hhh(:,:,3)))
   end subroutine make_legitimate_fixture

   subroutine test_first_derivative_oracle(fix, max_error)
      type(lmto_live_hamiltonian_fixture), intent(inout) :: fix
      real(rp), intent(out) :: max_error
      real(rp), parameter :: delta = 1.0e-5_rp
      real(rp), parameter :: axes(3,2) = reshape([0.0_rp,1.0_rp,0.0_rp, 1.0_rp,0.0_rp,0.0_rp], [3,2])
      real(rp), parameter :: kpoints(3,3) = reshape([0.0_rp,0.0_rp,0.0_rp, 0.17_rp,-0.09_rp,0.13_rp, -0.23_rp,0.11_rp,0.07_rp], [3,3])
      complex(rp) :: h0(12,12), hp(12,12), hm(12,12), analytic(12,12)
      real(rp) :: saved(3,2)
      integer :: ik, site

      saved = fix%moments; max_error = 0.0_rp
      do ik = 1, size(kpoints,2)
         do site = 1, 2
            call assemble_lmto_hamiltonian(fix, kpoints(:,ik), h0)
            fix%moments(:,site) = rotate_moment(saved(:,site), axes(:,site), delta)
            call assemble_lmto_hamiltonian(fix, kpoints(:,ik), hp)
            fix%moments(:,site) = rotate_moment(saved(:,site), axes(:,site), -delta)
            call assemble_lmto_hamiltonian(fix, kpoints(:,ik), hm)
            fix%moments = saved
            call assemble_lmto_torque(fix, kpoints(:,ik), site, axes(:,site), analytic)
            max_error = max(max_error, maxval(abs((hp-hm)/(2.0_rp*delta)-analytic)))
         end do
      end do
      fix%moments = saved
   end subroutine test_first_derivative_oracle

   subroutine test_mixed_derivative_oracle(fix, max_error)
      type(lmto_live_hamiltonian_fixture), intent(inout) :: fix
      real(rp), intent(out) :: max_error
      real(rp), parameter :: delta = 2.0e-4_rp
      real(rp), parameter :: axis_i(3) = [0.0_rp,1.0_rp,0.0_rp]
      real(rp), parameter :: axis_j(3) = [1.0_rp,0.0_rp,0.0_rp]
      real(rp), parameter :: kpoints(3,2) = reshape([0.13_rp,-0.08_rp,0.06_rp, -0.21_rp,0.12_rp,0.09_rp], [3,2])
      complex(rp) :: hpp(12,12), hpm(12,12), hmp(12,12), hmm(12,12), analytic(12,12), fd(12,12)
      real(rp) :: saved(3,2)
      integer :: ik

      saved = fix%moments; max_error = 0.0_rp
      do ik = 1, size(kpoints,2)
         fix%moments(:,1) = rotate_moment(saved(:,1), axis_i, delta)
         fix%moments(:,2) = rotate_moment(saved(:,2), axis_j, delta)
         call assemble_lmto_hamiltonian(fix, kpoints(:,ik), hpp)
         fix%moments(:,2) = rotate_moment(saved(:,2), axis_j, -delta)
         call assemble_lmto_hamiltonian(fix, kpoints(:,ik), hpm)
         fix%moments(:,1) = rotate_moment(saved(:,1), axis_i, -delta)
         fix%moments(:,2) = rotate_moment(saved(:,2), axis_j, delta)
         call assemble_lmto_hamiltonian(fix, kpoints(:,ik), hmp)
         fix%moments(:,2) = rotate_moment(saved(:,2), axis_j, -delta)
         call assemble_lmto_hamiltonian(fix, kpoints(:,ik), hmm)
         fix%moments = saved
         call assemble_lmto_mixed_derivative(fix, kpoints(:,ik), 1, axis_i, 2, axis_j, analytic)
         fd = (hpp-hpm-hmp+hmm)/(4.0_rp*delta*delta)
         max_error = max(max_error, maxval(abs(fd-analytic)))
      end do
      fix%moments = saved
   end subroutine test_mixed_derivative_oracle

   subroutine test_grand_potential_hessian(fix, max_error, contact_gap)
      type(lmto_live_hamiltonian_fixture), intent(inout) :: fix
      real(rp), intent(out) :: max_error, contact_gap
      real(rp), parameter :: delta = 2.0e-4_rp
      real(rp), parameter :: axis_i(3) = [0.0_rp,1.0_rp,0.0_rp]
      real(rp), parameter :: axis_j(3) = [1.0_rp,0.0_rp,0.0_rp]
      real(rp), parameter :: kpoint(3) = [0.13_rp,-0.08_rp,0.06_rp]
      complex(rp) :: h0(12,12), hpp(12,12), hpm(12,12), hmp(12,12), hmm(12,12), ti(12,12), tj(12,12), cij(12,12)
      real(rp) :: eval0(12), eval(12), fermi, tt, cc, complete
      complex(rp) :: vec0(12,12), vec(12,12)
      real(rp) :: omega_pp, omega_pm, omega_mp, omega_mm, hessian_fd, saved(3,2)

      saved = fix%moments
      call assemble_lmto_hamiltonian(fix, kpoint, h0)
      call assemble_lmto_torque(fix, kpoint, 1, axis_i, ti)
      call assemble_lmto_torque(fix, kpoint, 2, axis_j, tj)
      call assemble_lmto_mixed_derivative(fix, kpoint, 1, axis_i, 2, axis_j, cij)
      call diagonalize(h0, eval0, vec0)
      fermi = 0.5_rp*(eval0(4)+eval0(5))
      call force_theorem_hessian_from_eigenbasis(eval0, vec0, fermi, ti, tj, cij, tt, cc, complete)

      fix%moments(:,1) = rotate_moment(saved(:,1), axis_i, delta)
      fix%moments(:,2) = rotate_moment(saved(:,2), axis_j, delta)
      call assemble_lmto_hamiltonian(fix, kpoint, hpp); call diagonalize(hpp, eval, vec)
      omega_pp = grand_potential_from_eigenvalues(eval, fermi)
      fix%moments(:,2) = rotate_moment(saved(:,2), axis_j, -delta)
      call assemble_lmto_hamiltonian(fix, kpoint, hpm); call diagonalize(hpm, eval, vec)
      omega_pm = grand_potential_from_eigenvalues(eval, fermi)
      fix%moments(:,1) = rotate_moment(saved(:,1), axis_i, -delta)
      fix%moments(:,2) = rotate_moment(saved(:,2), axis_j, delta)
      call assemble_lmto_hamiltonian(fix, kpoint, hmp); call diagonalize(hmp, eval, vec)
      omega_mp = grand_potential_from_eigenvalues(eval, fermi)
      fix%moments(:,2) = rotate_moment(saved(:,2), axis_j, -delta)
      call assemble_lmto_hamiltonian(fix, kpoint, hmm); call diagonalize(hmm, eval, vec)
      omega_mm = grand_potential_from_eigenvalues(eval, fermi)
      fix%moments = saved
      hessian_fd = mixed_second_difference(omega_pp, omega_pm, omega_mp, omega_mm, delta, delta)
      max_error = abs(hessian_fd-complete)
      contact_gap = abs(hessian_fd-tt)-abs(hessian_fd-complete)
   end subroutine test_grand_potential_hessian

   subroutine diagonalize(matrix, values, vectors)
      complex(rp), intent(in) :: matrix(:, :)
      real(rp), intent(out) :: values(:)
      complex(rp), intent(out) :: vectors(:, :)
      complex(rp), allocatable :: work(:), query(:)
      real(rp), allocatable :: rwork(:)
      integer :: n, lwork, info
      n = size(values); vectors = matrix
      allocate(query(1), rwork(max(1,3*n-2)))
      call zheev('V','U',n,vectors,n,values,query,-1,rwork,info)
      lwork = max(1, int(real(query(1),rp))); allocate(work(lwork))
      call zheev('V','U',n,vectors,n,values,work,lwork,rwork,info)
      if (info /= 0) error stop 'DRESP-03T diagonalization failed'
      deallocate(work, query, rwork)
   end subroutine diagonalize

   pure function rotate_moment(moment, axis, angle) result(rotated)
      real(rp), intent(in) :: moment(3), axis(3), angle
      real(rp) :: rotated(3)
      rotated = moment*cos(angle) + cross3(axis,moment)*sin(angle) + axis*dot_product(axis,moment)*(1.0_rp-cos(angle))
   end function rotate_moment

   pure function cross3(a,b) result(c)
      real(rp), intent(in) :: a(3), b(3)
      real(rp) :: c(3)
      c = [a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), a(1)*b(2)-a(2)*b(1)]
   end function cross3

end program test_dresp03t_lmto_hessian
