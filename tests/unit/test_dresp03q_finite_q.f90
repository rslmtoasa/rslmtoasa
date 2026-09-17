!------------------------------------------------------------------------------
! DRESP-03Q -- finite-q torque/Hessian and independent supercell oracle.
!
! The primitive finite-q path supplies vertices and contacts.  The oracle below
! builds a period-two Hamiltonian directly from the LMTO bond value and never
! calls a torque or Hessian routine while evaluating the finite differences.
!------------------------------------------------------------------------------
program test_dresp03q_finite_q
   use precision_mod, only: rp
   use lmto_magnetic_tangent_mod, only: lmto_bond_value, lmto_hhmag_to_spinor
   use lr_kl_hessian_mod, only: lmto_live_hamiltonian_fixture, lmto_fixture_init, &
      assemble_lmto_hamiltonian, assemble_lmto_finite_q_torques, assemble_lmto_finite_q_mixed_derivative, &
      assemble_lmto_torque, assemble_lmto_mixed_derivative, force_theorem_finite_q_hessian_from_eigenbasis, &
      force_theorem_finite_q_hessian_from_eigenbasis_batch, force_theorem_hessian_from_eigenbasis
   implicit none

   type(lmto_live_hamiltonian_fixture) :: fixture
   integer, parameter :: ncell = 2, nmat = 2, nk = 2
   real(rp), parameter :: q(3) = [0.5_rp, 0.0_rp, 0.0_rp]
   real(rp), parameter :: axes(3,1) = reshape([0.0_rp, 1.0_rp, 0.0_rp], [3,1])
   real(rp), parameter :: kpoints(3,nk) = reshape([0.0_rp,0.0_rp,0.0_rp, 0.5_rp,0.0_rp,0.0_rp], [3,nk])
   real(rp), parameter :: weights(nk) = [0.5_rp, 0.5_rp]
   real(rp), parameter :: fermi = -0.1_rp, delta = 2.0e-5_rp
   complex(rp) :: torques_q(nmat,nmat,1,nk), torques_minus_q(nmat,nmat,1,nk)
   complex(rp) :: contacts(nmat,nmat,1,1,nk)
   real(rp) :: evals(nmat,nk), endpoint_evals(nmat,nk)
   complex(rp) :: hessian(1,1), tt(1,1), contact(1,1), complete(1,1)
   complex(rp) :: vecs(nmat,nmat,nk), endpoint_vecs(nmat,nmat,nk)
   real(rp) :: fd, contact_gap, q0_error, hermitian_error, reverse_error, periodic_error
   real(rp) :: multisite_hermitian_error
   integer :: ik

   call make_fixture(fixture)
   call prepare_primitive_data(fixture, evals, vecs, endpoint_evals, endpoint_vecs, torques_q, torques_minus_q, contacts)
   call force_theorem_finite_q_hessian_from_eigenbasis_batch(evals, vecs, endpoint_evals, endpoint_vecs, fermi, weights, &
      torques_q, torques_minus_q, contacts, hessian, tt, contact, complete)

   fd = supercell_mixed_difference(fixture, delta, fermi)
   contact_gap = abs(fd - tt(1,1)) - abs(fd - complete(1,1))
   q0_error = qzero_reduction(fixture)
   hermitian_error = abs(finite_hessian_for_q(fixture,q) - finite_hessian_for_q(fixture,-q))
   reverse_error = q_reversal_check(fixture)
   periodic_error = q_periodicity_check(fixture)
   multisite_hermitian_error = multisite_hessian_hermiticity_check()

   write (*, '(a,es12.4)') 'DRESP-03Q period-two finite-q complete-vs-FD error = ', abs(complete(1,1)-fd)
   write (*, '(a,es12.4)') 'DRESP-03Q finite-difference Hessian             = ', fd
   write (*, '(a,es12.4)') 'DRESP-03Q spectral TT                         = ', real(tt(1,1),rp)
   write (*, '(a,es12.4)') 'DRESP-03Q spectral contact                    = ', real(contact(1,1),rp)
   write (*, '(a,es12.4)') 'DRESP-03Q spectral complete                   = ', real(complete(1,1),rp)
   write (*, '(a,es12.4)') 'DRESP-03Q period-two TT-only-vs-FD error       = ', abs(tt(1,1)-fd)
   write (*, '(a,es12.4)') 'DRESP-03Q contact contribution                 = ', abs(contact(1,1))
   write (*, '(a,es12.4)') 'DRESP-03Q q=0 reduction error                  = ', q0_error
   write (*, '(a,es12.4)') 'DRESP-03Q Hermitian Hessian error               = ', hermitian_error
   write (*, '(a,es12.4)') 'DRESP-03Q q/-q covariance error                 = ', reverse_error
   write (*, '(a,es12.4)') 'DRESP-03Q reciprocal-period covariance error   = ', periodic_error
   write (*, '(a,es12.4)') 'DRESP-03Q multi-site Hessian Hermitian error   = ', multisite_hermitian_error

   if (abs(complete(1,1)-fd) > 3.0e-7_rp .or. abs(contact(1,1)) < 1.0e-7_rp .or. contact_gap < 1.0e-7_rp .or. &
       q0_error > 3.0e-12_rp .or. hermitian_error > 3.0e-12_rp .or. reverse_error > 3.0e-10_rp .or. &
       periodic_error > 3.0e-10_rp .or. multisite_hermitian_error > 3.0e-10_rp) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'
   call fixture%clear()

contains

   subroutine make_fixture(fix)
      type(lmto_live_hamiltonian_fixture), intent(out) :: fix
      call lmto_fixture_init(fix, 1, 1, 3, .false.)
      fix%moments(:,1) = [0.0_rp, 0.0_rp, 1.0_rp]
      fix%site_position(:,1) = 0.0_rp
      fix%wx0(:,1) = [1.0_rp]
      fix%wx1(:,1) = [0.08_rp]
      fix%c0(:,1) = [-0.2_rp]
      fix%c1(:,1) = [-0.5_rp]
      fix%bond_source = [1, 1, 1]
      fix%bond_target = [1, 1, 1]
      fix%onsite = [.true., .false., .false.]
      fix%bond_vector(:,1) = 0.0_rp
      fix%bond_vector(:,2) = [1.0_rp, 0.0_rp, 0.0_rp]
      fix%bond_vector(:,3) = [-1.0_rp, 0.0_rp, 0.0_rp]
      fix%hhh(:,:,1) = cmplx(0.0_rp,0.0_rp,rp)
      fix%hhh(:,:,2) = cmplx(0.10_rp,0.0_rp,rp)
      fix%hhh(:,:,3) = cmplx(0.10_rp,0.0_rp,rp)
   end subroutine make_fixture

   subroutine prepare_primitive_data(fix, values, vectors, endpoint_values, endpoint_vectors, tq, tm, c)
      type(lmto_live_hamiltonian_fixture), intent(in) :: fix
      real(rp), intent(out) :: values(:, :), endpoint_values(:, :)
      complex(rp), intent(out) :: vectors(:, :, :), endpoint_vectors(:, :, :), tq(:, :, :, :), tm(:, :, :, :), c(:, :, :, :, :)
      complex(rp) :: h(nmat,nmat), hq(nmat,nmat), vq(nmat,nmat), vm(nmat,nmat), cm(nmat,nmat)
      integer :: ik

      do ik = 1, nk
         call assemble_lmto_hamiltonian(fix, kpoints(:,ik), h)
         call assemble_lmto_hamiltonian(fix, kpoints(:,ik)+q, hq)
         call diagonalize(h, values(:,ik), vectors(:,:,ik))
         call diagonalize(hq, endpoint_values(:,ik), endpoint_vectors(:,:,ik))
         call assemble_lmto_finite_q_torques(fix, kpoints(:,ik), q, axes, tq(:,:,:,ik))
         call assemble_lmto_finite_q_torques(fix, kpoints(:,ik)+q, -q, axes, tm(:,:,:,ik))
         call assemble_lmto_finite_q_mixed_derivative(fix, kpoints(:,ik), q, 1, axes(:,1), 1, axes(:,1), cm)
         c(:,:,:,:,ik) = 0.0_rp
         c(:,:,1,1,ik) = cm
         ! Keep the local variables visible in the test: these are the exact
         ! k -> k+q and k+q -> k vertices used by the spectral contraction.
         vq = tq(:,:,1,ik); vm = tm(:,:,1,ik)
      end do
   end subroutine prepare_primitive_data

   function supercell_mixed_difference(fix, step, ef) result(value)
      type(lmto_live_hamiltonian_fixture), intent(in) :: fix
      real(rp), intent(in) :: step, ef
      real(rp) :: value
      real(rp) :: omega(2,2)
      integer :: sx, sy
      do sx = 1, 2
         do sy = 1, 2
            omega(sx,sy) = supercell_grand_potential(fix, merge(step,-step,sx == 1), &
               merge(step,-step,sy == 1), ef)
         end do
      end do
      value = (omega(1,1)-omega(1,2)-omega(2,1)+omega(2,2))/(4.0_rp*step*step)
   end function supercell_mixed_difference

   function supercell_grand_potential(fix, amplitude_q, amplitude_minus_q, ef) result(omega)
      type(lmto_live_hamiltonian_fixture), intent(in) :: fix
      real(rp), intent(in) :: amplitude_q, amplitude_minus_q, ef
      real(rp) :: omega
      complex(rp) :: h(2*ncell,2*ncell), block(2,2), hmag(1,1,4), vectors(2*ncell,2*ncell)
      real(rp) :: values(2*ncell), moment(3), angle
      integer :: cell, target_cell, ibond, source_index, target_index

      h = cmplx(0.0_rp,0.0_rp,rp)
      do cell = 0, ncell-1
         angle = (amplitude_q + amplitude_minus_q)*merge(1.0_rp,-1.0_rp,mod(cell,2) == 0)
         moment = rotate_moment(fix%moments(:,1), axes(:,1), angle)
         do ibond = 1, fix%nbond
            target_cell = modulo(cell + nint(fix%bond_vector(1,ibond)), ncell)
            call lmto_bond_value(fix%hhh(:,:,ibond), fix%wx0(:,1), fix%wx1(:,1), fix%wx0(:,1), fix%wx1(:,1), &
               fix%c0(:,1), fix%c1(:,1), moment, rotate_moment(fix%moments(:,1), axes(:,1), &
               (amplitude_q + amplitude_minus_q)*merge(1.0_rp,-1.0_rp,mod(target_cell,2) == 0)), fix%onsite(ibond), hmag)
            call lmto_hhmag_to_spinor(hmag, block)
            source_index = 2*cell+1; target_index = 2*target_cell+1
            h(source_index:source_index+1,target_index:target_index+1) = &
               h(source_index:source_index+1,target_index:target_index+1) + block
         end do
      end do
      call diagonalize(h, values, vectors)
      omega = 0.0_rp
      do cell = 1, size(values)
         if (values(cell) < ef) omega = omega + values(cell) - ef
      end do
      omega = omega/real(ncell,rp)
   end function supercell_grand_potential

   function qzero_reduction(fix) result(error)
      type(lmto_live_hamiltonian_fixture), intent(in) :: fix
      real(rp) :: error
      complex(rp) :: old_t(2,2), new_t(2,2), new_stack(2,2,1), minus_stack(2,2,1), cstack(2,2,1,1)
      complex(rp) :: h(2,2), vectors(2,2)
      real(rp) :: values(2), old_tt, old_contact, old_complete
      complex(rp) :: new_h(1,1), new_tt(1,1), new_contact(1,1), new_complete(1,1)
      real(rp) :: k(3)
      k = [0.17_rp,0.0_rp,0.0_rp]
      call assemble_lmto_hamiltonian(fix,k,h)
      call diagonalize(h,values,vectors)
      call assemble_lmto_torque(fix,k,1,axes(:,1),old_t)
      call assemble_lmto_finite_q_torques(fix,k,[0.0_rp,0.0_rp,0.0_rp],axes,new_stack)
      call assemble_lmto_finite_q_torques(fix,k,[0.0_rp,0.0_rp,0.0_rp],axes,minus_stack)
      call assemble_lmto_finite_q_mixed_derivative(fix,k,[0.0_rp,0.0_rp,0.0_rp],1,axes(:,1),1,axes(:,1),cstack(:,:,1,1))
      new_t=new_stack(:,:,1)
      call force_theorem_hessian_from_eigenbasis(values,vectors,fermi,old_t,old_t,cstack(:,:,1,1), &
         old_tt,old_contact,old_complete)
      call force_theorem_finite_q_hessian_from_eigenbasis(values,vectors,values,vectors,fermi,new_stack,minus_stack,cstack, &
         new_h,new_tt,new_contact,new_complete)
      error = max(maxval(abs(old_t-new_t)),abs(old_tt-new_tt(1,1)),abs(old_contact-new_contact(1,1)), &
         abs(old_complete-new_complete(1,1)))
   end function qzero_reduction

   function q_reversal_check(fix) result(error)
      type(lmto_live_hamiltonian_fixture), intent(in) :: fix
      real(rp) :: error
      complex(rp) :: plus(2,2,1), minus(2,2,1), expected(2,2)
      call assemble_lmto_finite_q_torques(fix,[0.13_rp,0.0_rp,0.0_rp],q,axes,plus)
      call assemble_lmto_finite_q_torques(fix,[0.13_rp,0.0_rp,0.0_rp]+q,-q,axes,minus)
      expected = transpose(conjg(minus(:,:,1)))
      error = maxval(abs(plus(:,:,1)-expected))
   end function q_reversal_check

   function finite_hessian_for_q(fix, q_in) result(value)
      type(lmto_live_hamiltonian_fixture), intent(in) :: fix
      real(rp), intent(in) :: q_in(3)
      complex(rp) :: value
      real(rp) :: values(nmat,nk), endpoint_values(nmat,nk)
      complex(rp) :: hessian(1,1), tt(1,1), cc(1,1), complete(1,1)
      complex(rp) :: vectors(nmat,nmat,nk), endpoint_vectors(nmat,nmat,nk)
      complex(rp) :: tq(nmat,nmat,1,nk), tm(nmat,nmat,1,nk), c(nmat,nmat,1,1,nk)
      complex(rp) :: h(nmat,nmat), hq(nmat,nmat)
      integer :: local_ik

      do local_ik = 1, nk
         call assemble_lmto_hamiltonian(fix, kpoints(:,local_ik), h)
         call assemble_lmto_hamiltonian(fix, kpoints(:,local_ik)+q_in, hq)
         call diagonalize(h, values(:,local_ik), vectors(:,:,local_ik))
         call diagonalize(hq, endpoint_values(:,local_ik), endpoint_vectors(:,:,local_ik))
         call assemble_lmto_finite_q_torques(fix, kpoints(:,local_ik), q_in, axes, tq(:,:,:,local_ik))
         call assemble_lmto_finite_q_torques(fix, kpoints(:,local_ik)+q_in, -q_in, axes, tm(:,:,:,local_ik))
         call assemble_lmto_finite_q_mixed_derivative(fix, kpoints(:,local_ik), q_in, 1, axes(:,1), 1, axes(:,1), h)
         c(:,:,1,1,local_ik) = h
      end do
      call force_theorem_finite_q_hessian_from_eigenbasis_batch(values, vectors, endpoint_values, endpoint_vectors, fermi, weights, &
         tq, tm, c, hessian, tt, cc, complete)
      value = complete(1,1)
   end function finite_hessian_for_q

   function multisite_hessian_hermiticity_check() result(error)
      type(lmto_live_hamiltonian_fixture) :: fix
      complex(rp) :: hessian(2,2), hessian_minus(2,2)
      real(rp) :: qpoint(3)
      real(rp) :: error, vertex_error
      qpoint = [0.50_rp, 0.0_rp, 0.0_rp]
      call make_multisite_fixture(fix)
      call evaluate_multisite_hessian(fix, qpoint, hessian, vertex_error)
      call evaluate_multisite_hessian(fix, -qpoint, hessian_minus)
      error = max(maxval(abs(hessian-transpose(conjg(hessian)))), vertex_error, &
         maxval(abs(hessian_minus-conjg(hessian))))
      call fix%clear()
   end function multisite_hessian_hermiticity_check

   subroutine make_multisite_fixture(fix)
      type(lmto_live_hamiltonian_fixture), intent(out) :: fix
      call lmto_fixture_init(fix, 2, 1, 4, .false.)
      fix%moments(:,1) = [0.0_rp, 0.0_rp, 1.0_rp]
      fix%moments(:,2) = [0.0_rp, 0.0_rp, 1.0_rp]
      fix%site_position(:,1) = [0.0_rp, 0.0_rp, 0.0_rp]
      fix%site_position(:,2) = [0.25_rp, 0.0_rp, 0.0_rp]
      fix%wx0 = 1.0_rp
      fix%wx1 = 0.08_rp
      fix%c0(:,1) = [-0.20_rp]
      fix%c0(:,2) = [-0.15_rp]
      fix%c1(:,1) = [-0.50_rp]
      fix%c1(:,2) = [-0.42_rp]
      fix%bond_source = [1, 2, 1, 2]
      fix%bond_target = [1, 2, 2, 1]
      fix%onsite = [.true., .true., .false., .false.]
      fix%bond_vector(:,1) = 0.0_rp
      fix%bond_vector(:,2) = 0.0_rp
      fix%bond_vector(:,3) = [0.25_rp, 0.0_rp, 0.0_rp]
      fix%bond_vector(:,4) = [-0.25_rp, 0.0_rp, 0.0_rp]
      fix%hhh(:,:,1) = 0.0_rp
      fix%hhh(:,:,2) = 0.0_rp
      fix%hhh(:,:,3) = cmplx(0.10_rp, 0.0_rp, rp)
      fix%hhh(:,:,4) = cmplx(0.10_rp, 0.0_rp, rp)
   end subroutine make_multisite_fixture

   subroutine evaluate_multisite_hessian(fix, qpoint, result, vertex_adjoint_error)
      type(lmto_live_hamiltonian_fixture), intent(in) :: fix
      real(rp), intent(in) :: qpoint(3)
      complex(rp), intent(out) :: result(2,2)
      real(rp), intent(out), optional :: vertex_adjoint_error
      integer, parameter :: nmat2 = 4, nk2 = 2
      real(rp), parameter :: kpoints2(3,nk2) = reshape([0.0_rp,0.0_rp,0.0_rp, 0.50_rp,0.0_rp,0.0_rp], [3,nk2])
      real(rp), parameter :: weights2(nk2) = [0.5_rp, 0.5_rp], fermi2 = -0.10_rp
      real(rp) :: values(nmat2,nk2), endpoint_values(nmat2,nk2)
      complex(rp) :: vectors(nmat2,nmat2,nk2), endpoint_vectors(nmat2,nmat2,nk2)
      complex(rp) :: tq(nmat2,nmat2,2,nk2), tm(nmat2,nmat2,2,nk2), c(nmat2,nmat2,2,2,nk2)
      complex(rp) :: h(nmat2,nmat2)
      complex(rp) :: tt(2,2), cc(2,2), complete(2,2)
      real(rp), parameter :: axes2(3,2) = reshape([0.0_rp,1.0_rp,0.0_rp, 0.0_rp,1.0_rp,0.0_rp], [3,2])
      real(rp) :: vertex_error
      integer :: ik, site, site_other

      do ik = 1, nk2
         call assemble_lmto_hamiltonian(fix, kpoints2(:,ik), h)
         call assemble_lmto_hamiltonian(fix, kpoints2(:,ik)+qpoint, c(:,:,1,1,ik))
         call diagonalize(h, values(:,ik), vectors(:,:,ik))
         call diagonalize(c(:,:,1,1,ik), endpoint_values(:,ik), endpoint_vectors(:,:,ik))
         call assemble_lmto_finite_q_torques(fix, kpoints2(:,ik), qpoint, axes2, tq(:,:,:,ik))
         call assemble_lmto_finite_q_torques(fix, kpoints2(:,ik)+qpoint, -qpoint, axes2, tm(:,:,:,ik))
         do site = 1, 2
            do site_other = 1, 2
               call assemble_lmto_finite_q_mixed_derivative(fix, kpoints2(:,ik), qpoint, site, axes2(:,site), &
                  site_other, axes2(:,site_other), c(:,:,site,site_other,ik))
            end do
         end do
      end do
      vertex_error = 0.0_rp
      do ik = 1, nk2
         do site = 1, 2
            vertex_error = max(vertex_error, maxval(abs(tq(:,:,site,ik)-transpose(conjg(tm(:,:,site,ik))))))
         end do
      end do
      if (present(vertex_adjoint_error)) vertex_adjoint_error = vertex_error
      call force_theorem_finite_q_hessian_from_eigenbasis_batch(values, vectors, endpoint_values, endpoint_vectors, fermi2, weights2, &
         tq, tm, c, result, tt, cc, complete)
   end subroutine evaluate_multisite_hessian

   function q_periodicity_check(fix) result(error)
      type(lmto_live_hamiltonian_fixture), intent(in) :: fix
      real(rp) :: error
      complex(rp) :: a(2,2,1), b(2,2,1)
      call assemble_lmto_finite_q_torques(fix,[0.13_rp,0.0_rp,0.0_rp],q,axes,a)
      call assemble_lmto_finite_q_torques(fix,[0.13_rp,0.0_rp,0.0_rp],q+[1.0_rp,0.0_rp,0.0_rp],axes,b)
      error = maxval(abs(a-b))
   end function q_periodicity_check

   subroutine diagonalize(matrix, values, vectors)
      complex(rp), intent(in) :: matrix(:, :)
      real(rp), intent(out) :: values(:)
      complex(rp), intent(out) :: vectors(:, :)
      complex(rp), allocatable :: work(:)
      real(rp), allocatable :: rwork(:)
      complex(rp) :: query(1)
      integer :: n, lwork, info
      n=size(values); vectors=matrix; allocate(rwork(max(1,3*n-2)))
      call zheev('V','U',n,vectors,n,values,query,-1,rwork,info)
      lwork=max(1,nint(real(query(1),rp))); allocate(work(lwork))
      call zheev('V','U',n,vectors,n,values,work,lwork,rwork,info)
      if (info /= 0) error stop 'DRESP-03Q diagonalization failed'
      deallocate(work,rwork)
   end subroutine diagonalize

   subroutine supercell_matrix(fix, amplitude_q, amplitude_minus_q, h)
      type(lmto_live_hamiltonian_fixture), intent(in) :: fix
      real(rp), intent(in) :: amplitude_q, amplitude_minus_q
      complex(rp), intent(out) :: h(:, :)
      complex(rp) :: block(2,2), hmag(1,1,4)
      real(rp) :: moment(3), angle
      integer :: cell, target_cell, ibond, source_index, target_index

      h = cmplx(0.0_rp,0.0_rp,rp)
      do cell = 0, ncell-1
         angle = (amplitude_q + amplitude_minus_q)*merge(1.0_rp,-1.0_rp,mod(cell,2) == 0)
         moment = rotate_moment(fix%moments(:,1), axes(:,1), angle)
         do ibond = 1, fix%nbond
            target_cell = modulo(cell + nint(fix%bond_vector(1,ibond)), ncell)
            call lmto_bond_value(fix%hhh(:,:,ibond), fix%wx0(:,1), fix%wx1(:,1), fix%wx0(:,1), fix%wx1(:,1), &
               fix%c0(:,1), fix%c1(:,1), moment, rotate_moment(fix%moments(:,1), axes(:,1), &
               (amplitude_q + amplitude_minus_q)*merge(1.0_rp,-1.0_rp,mod(target_cell,2) == 0)), fix%onsite(ibond), hmag)
            call lmto_hhmag_to_spinor(hmag, block)
            source_index = 2*cell+1; target_index = 2*target_cell+1
            h(source_index:source_index+1,target_index:target_index+1) = &
               h(source_index:source_index+1,target_index:target_index+1) + block
         end do
      end do
   end subroutine supercell_matrix

   pure function rotate_moment(moment, axis, angle) result(rotated)
      real(rp), intent(in) :: moment(3), axis(3), angle
      real(rp) :: rotated(3)
      rotated=moment*cos(angle)+cross3(axis,moment)*sin(angle)+axis*dot_product(axis,moment)*(1.0_rp-cos(angle))
   end function rotate_moment

   pure function cross3(a,b) result(c)
      real(rp), intent(in) :: a(3), b(3)
      real(rp) :: c(3)
      c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
   end function cross3

end program test_dresp03q_finite_q
