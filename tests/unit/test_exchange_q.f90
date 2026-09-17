!------------------------------------------------------------------------------
! exchange_q input boundary and independent q-space LKAG contraction tests.
!------------------------------------------------------------------------------
program test_exchange_q
   use precision_mod, only: rp
   use exchange_q_mod, only: exchange_q_config, load_exchange_q_config, lkag_pauli_product, &
      detect_commensurate_endpoint_map
   use lr_kl_hessian_mod, only: force_theorem_hessian_from_eigenbasis, &
      force_theorem_finite_q_hessian_from_eigenbasis_metallic, fermi_divided_difference, &
      finite_temperature_occupation
   implicit none

   type(exchange_q_config) :: config
   complex(rp) :: one(1,1), zero(1,1), ginmag(1,1), gjnmag(1,1)
   real(rp) :: value, expected, q_space, real_space, phase
   real(rp) :: residual, map_error, gapped_error, oracle_error, degenerate_error, fd_error, kernel, expected_tt
   real(rp) :: eigenvalues2(2), endpoint_values2(2), vectors2(2,2), mixed2(2,2)
   real(rp) :: hstep, omega_plus, omega_minus, omega_zero
   complex(rp) :: vec2(2,2), endpoint_vec2(2,2), torque2(2,2,1), torque_minus2(2,2,1), mixed_q2(2,2,1,1)
   real(rp) :: eigenvalues4(4), endpoint_values4(4), mixed4(4,4)
   complex(rp) :: vec4(4,4), endpoint_vec4(4,4), torque4(4,4,1), torque_minus4(4,4,1), mixed_q4(4,4,1,1)
   real(rp) :: eigenvalues3(3), endpoint_values3(3), mixed3(3,3)
   complex(rp) :: vec3(3,3), endpoint_vec3(3,3), torque3(3,3,1), torque_minus3(3,3,1), mixed_q3(3,3,1,1)
   complex(rp) :: hessian(1,1), tt(1,1), contact(1,1), complete(1,1)
   integer :: mesh(3), nk_mesh, ik, ix, iy, iz
   real(rp), allocatable :: mesh_points(:, :), mesh_weights(:), periodic_eigenvalues(:, :)
   integer, allocatable :: endpoint_index(:)
   logical :: commensurate
   integer :: unit

   call write_inline_input('exchange_q_unit.nml')
   call load_exchange_q_config('exchange_q_unit.nml', config, .true.)
   if (config%n_q /= 3 .or. abs(config%q_list(3,2)-0.25_rp) > 1.0e-14_rp .or. &
       abs(config%q_list(1,3)-0.5_rp) > 1.0e-14_rp) then
      error stop 'UnitExchangeQ: inline q path parsing failed'
   end if
   if (config%finite_h_response_backend /= 'both' .or. config%contour_points /= 48 .or. &
       config%contour_shape /= 'ellipse' .or. .not. config%contour_account_fermi_poles) then
      error stop 'UnitExchangeQ: finite-H contour controls parsing failed'
   end if

   call write_q_file('exchange_q_unit_path.dat')
   call write_file_input('exchange_q_unit.nml')
   call load_exchange_q_config('exchange_q_unit.nml', config, .true.)
   if (config%n_q /= 2 .or. abs(config%q_list(3,2)-0.5_rp) > 1.0e-14_rp) then
      error stop 'UnitExchangeQ: frozen-magnon q-file parsing failed'
   end if

   one = cmplx(1.0_rp,0.0_rp,rp)
   zero = cmplx(0.0_rp,0.0_rp,rp)
   ginmag = cmplx(0.3_rp,0.7_rp,rp)
   gjnmag = cmplx(0.2_rp,-0.4_rp,rp)
   value = lkag_pauli_product(one, one, ginmag, zero, zero, zero, gjnmag, zero, zero, zero)
   expected = aimag(ginmag(1,1)*gjnmag(1,1))
   if (abs(value-expected) > 1.0e-14_rp) error stop 'UnitExchangeQ: LKAG Pauli contraction failed'

   ! Two-point DFT identity: the q-space k,k+q product equals the real-space
   ! Fourier sum with the same phase convention.  This is a minimal oracle for
   ! the independent Green-function representation used by exchange_q.
   q_space = 0.5_rp*aimag(cmplx(0.2_rp,0.1_rp,rp)*cmplx(0.6_rp,-0.2_rp,rp) + &
      cmplx(-0.1_rp,0.5_rp,rp)*cmplx(-0.4_rp,0.3_rp,rp))
   real_space = 0.0_rp
   real_space = real_space + aimag(0.5_rp*(cmplx(0.2_rp,0.1_rp,rp)+cmplx(-0.1_rp,0.5_rp,rp))* &
      (0.5_rp*(cmplx(-0.4_rp,0.3_rp,rp)+cmplx(0.6_rp,-0.2_rp,rp))))
   phase = -1.0_rp
   real_space = real_space + aimag(phase*0.5_rp*(cmplx(0.2_rp,0.1_rp,rp)-cmplx(-0.1_rp,0.5_rp,rp))* &
      (0.5_rp*(cmplx(-0.4_rp,0.3_rp,rp)-cmplx(0.6_rp,-0.2_rp,rp))))
   if (abs(q_space-real_space) > 1.0e-14_rp) error stop 'UnitExchangeQ: q/real-space DFT oracle failed'

   call test_metallic_gapped_reduction(gapped_error)
   call test_finite_temperature_two_level(oracle_error, fd_error)
   call test_degenerate_subspace_invariance(degenerate_error)
   call test_commensurate_endpoint_reuse(map_error)
   if (gapped_error > 2.0e-12_rp) error stop 'UnitExchangeQ: gapped metallic reduction failed'
   if (oracle_error > 2.0e-10_rp .or. fd_error > 2.0e-7_rp) error stop 'UnitExchangeQ: finite-temperature oracle failed'
   if (degenerate_error > 2.0e-12_rp) error stop 'UnitExchangeQ: degenerate-subspace invariance failed'
   if (map_error > 2.0e-12_rp) error stop 'UnitExchangeQ: commensurate endpoint reuse failed'

   open(newunit=unit, file='exchange_q_unit.nml', status='old')
   close(unit, status='delete')
   open(newunit=unit, file='exchange_q_unit_path.dat', status='old')
   close(unit, status='delete')
   write (*, '(a,es12.4)') 'UnitExchangeQ: LKAG Pauli residual = ', abs(value-expected)
   write (*, '(a,es12.4)') 'UnitExchangeQ: q/real-space DFT residual = ', abs(q_space-real_space)
   write (*, '(a,es12.4)') 'UnitExchangeQ: gapped metallic reduction residual = ', gapped_error
   write (*, '(a,2(es12.4,1x))') 'UnitExchangeQ: finite-T spectral/FD residuals = ', oracle_error, fd_error
   write (*, '(a,es12.4)') 'UnitExchangeQ: degenerate-subspace residual = ', degenerate_error
   write (*, '(a,es12.4)') 'UnitExchangeQ: commensurate reuse residual = ', map_error
   write (*, '(a)') 'UnitExchangeQ: PASS (metallic spectral, degeneracy, endpoint-reuse, input, and DFT oracles)'

contains

   subroutine write_inline_input(path)
      character(len=*), intent(in) :: path
      integer :: local_unit
      open(newunit=local_unit, file=path, status='replace', action='write')
      write(local_unit,'(a)') '&exchange_q'
      write(local_unit,'(a)') " q_coordinates = 'direct'"
      write(local_unit,'(a)') " q_file = ''"
      write(local_unit,'(a)') ' n_q_points = 3'
      write(local_unit,'(a)') ' q_list = 0.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.0, 0.0'
      write(local_unit,'(a)') " output_file = 'unit_exchange_q.dat'"
      write(local_unit,'(a)') ' native_crosscheck = .false.'
      write(local_unit,'(a)') " finite_h_response_backend = 'both'"
      write(local_unit,'(a)') ' contour_points = 48'
      write(local_unit,'(a)') " contour_shape = 'ellipse'"
      write(local_unit,'(a)') ' contour_account_fermi_poles = .true.'
      write(local_unit,'(a)') '/'
      close(local_unit)
   end subroutine write_inline_input

   subroutine write_file_input(path)
      character(len=*), intent(in) :: path
      integer :: local_unit
      open(newunit=local_unit, file=path, status='replace', action='write')
      write(local_unit,'(a)') '&exchange_q'
      write(local_unit,'(a)') " q_coordinates = 'cartesian'"
      write(local_unit,'(a)') " q_file = 'exchange_q_unit_path.dat'"
      write(local_unit,'(a)') ' n_q_points = 0'
      write(local_unit,'(a)') '/'
      close(local_unit)
   end subroutine write_file_input

   subroutine write_q_file(path)
      character(len=*), intent(in) :: path
      integer :: local_unit
      open(newunit=local_unit, file=path, status='replace', action='write')
      write(local_unit,'(a)') '2'
      write(local_unit,'(a)') '0.0 0.0 0.0'
      write(local_unit,'(a)') '0.0 0.0 0.5'
      close(local_unit)
   end subroutine write_q_file

   subroutine test_metallic_gapped_reduction(error)
      real(rp), intent(out) :: error
      real(rp) :: old_tt, old_cc, old_all, new_tt, new_cc, new_all
      integer :: i, j
      eigenvalues4 = [-1.0_rp, -0.5_rp, 0.4_rp, 0.9_rp]
      endpoint_values4 = eigenvalues4
      vec4 = cmplx(0.0_rp,0.0_rp,rp); endpoint_vec4 = vec4
      do i = 1, 4
         vec4(i,i) = cmplx(1.0_rp,0.0_rp,rp); endpoint_vec4(i,i) = vec4(i,i)
      end do
      torque4 = cmplx(0.0_rp,0.0_rp,rp); torque_minus4 = torque4
      do i = 1, 4
         do j = i, 4
            torque4(i,j,1) = cmplx(0.03_rp*real(i+j,rp), 0.01_rp*real(j-i,rp), rp)
            torque4(j,i,1) = conjg(torque4(i,j,1))
         end do
      end do
      torque_minus4 = torque4
      mixed_q4 = cmplx(0.0_rp,0.0_rp,rp)
      do i = 1, 4
         mixed_q4(i,i,1,1) = cmplx(0.02_rp*real(i,rp),0.0_rp,rp)
      end do
      call force_theorem_hessian_from_eigenbasis(eigenvalues4, vec4, 0.0_rp, torque4(:,:,1), torque_minus4(:,:,1), mixed_q4(:,:,1,1), &
         old_tt, old_cc, old_all)
      call force_theorem_finite_q_hessian_from_eigenbasis_metallic(eigenvalues4, vec4, endpoint_values4, endpoint_vec4, 0.0_rp, &
         1.0e-10_rp, torque4, torque_minus4, mixed_q4, hessian, tt, contact, complete)
      new_tt = real(tt(1,1),rp); new_cc = real(contact(1,1),rp); new_all = real(complete(1,1),rp)
      error = max(abs(old_tt-new_tt),abs(old_cc-new_cc),abs(old_all-new_all))
   end subroutine test_metallic_gapped_reduction

   subroutine test_finite_temperature_two_level(spectral_error, finite_difference_error)
      real(rp), intent(out) :: spectral_error, finite_difference_error
      real(rp), parameter :: mu = 0.0_rp, kT = 0.05_rp, v = 0.7_rp
      real(rp) :: f1, f2, e1, e2, lambda1, lambda2
      e1 = -0.3_rp; e2 = 0.2_rp
      eigenvalues2 = [e1,e2]; endpoint_values2 = eigenvalues2
      vec2 = cmplx(0.0_rp,0.0_rp,rp); endpoint_vec2 = vec2
      vec2(1,1) = cmplx(1.0_rp,0.0_rp,rp); vec2(2,2) = vec2(1,1); endpoint_vec2 = vec2
      torque2 = cmplx(0.0_rp,0.0_rp,rp); torque_minus2 = torque2
      torque2(1,2,1) = cmplx(v,0.0_rp,rp); torque2(2,1,1) = conjg(torque2(1,2,1)); torque_minus2 = torque2
      mixed_q2 = cmplx(0.0_rp,0.0_rp,rp)
      call force_theorem_finite_q_hessian_from_eigenbasis_metallic(eigenvalues2, vec2, endpoint_values2, endpoint_vec2, mu, kT, &
         torque2, torque_minus2, mixed_q2, hessian, tt, contact, complete)
      f1 = finite_temperature_occupation(e1,mu,kT); f2 = finite_temperature_occupation(e2,mu,kT)
      kernel = (f1-f2)/(e1-e2); expected_tt = 2.0_rp*kernel*v*v
      spectral_error = abs(real(tt(1,1),rp)-expected_tt)
      hstep = 1.0e-4_rp
      lambda1 = 0.5_rp*(e1+e2)-sqrt((0.5_rp*(e1-e2))**2+(hstep*v)**2)
      lambda2 = 0.5_rp*(e1+e2)+sqrt((0.5_rp*(e1-e2))**2+(hstep*v)**2)
      omega_plus = oracle_grand_potential([lambda1,lambda2],mu,kT)
      lambda1 = 0.5_rp*(e1+e2)-sqrt((0.5_rp*(e1-e2))**2+(-hstep*v)**2)
      lambda2 = 0.5_rp*(e1+e2)+sqrt((0.5_rp*(e1-e2))**2+(-hstep*v)**2)
      omega_minus = oracle_grand_potential([lambda1,lambda2],mu,kT)
      omega_zero = oracle_grand_potential([e1,e2],mu,kT)
      finite_difference_error = abs((omega_plus-2.0_rp*omega_zero+omega_minus)/(hstep*hstep)-expected_tt)
      if (abs(fermi_divided_difference(0.1_rp,0.1_rp,0.1_rp,kT) + 0.25_rp/kT) > 2.0e-12_rp) then
         finite_difference_error = huge(1.0_rp)
      end if
   end subroutine test_finite_temperature_two_level

   subroutine test_degenerate_subspace_invariance(error)
      real(rp), intent(out) :: error
      real(rp) :: c, s, value_a, value_b
      complex(rp) :: rotated(3,3)
      eigenvalues3 = [0.1_rp,0.1_rp,-0.2_rp]; endpoint_values3 = eigenvalues3
      vec3 = cmplx(0.0_rp,0.0_rp,rp); endpoint_vec3 = vec3
      vec3(1,1)=1.0_rp; vec3(2,2)=1.0_rp; vec3(3,3)=1.0_rp; endpoint_vec3=vec3
      torque3 = cmplx(0.0_rp,0.0_rp,rp); torque_minus3 = torque3
      torque3(:,:,1) = reshape([cmplx(0.2_rp,0.0_rp,rp),cmplx(0.1_rp,0.03_rp,rp),cmplx(-0.04_rp,0.02_rp,rp), &
         cmplx(0.1_rp,-0.03_rp,rp),cmplx(-0.1_rp,0.0_rp,rp),cmplx(0.05_rp,0.01_rp,rp), &
         cmplx(-0.04_rp,-0.02_rp,rp),cmplx(0.05_rp,-0.01_rp,rp),cmplx(0.3_rp,0.0_rp,rp)], [3,3])
      torque_minus3 = torque3; mixed_q3 = cmplx(0.0_rp,0.0_rp,rp); mixed_q3(1,1,1,1)=0.04_rp; mixed_q3(2,2,1,1)=-0.02_rp
      call force_theorem_finite_q_hessian_from_eigenbasis_metallic(eigenvalues3, vec3, endpoint_values3, endpoint_vec3, 0.0_rp, 0.03_rp, &
         torque3, torque_minus3, mixed_q3, hessian, tt, contact, complete)
      value_a = real(complete(1,1),rp)
      c = cos(0.371_rp); s = sin(0.371_rp); rotated = vec3
      rotated(1,1)=c; rotated(2,1)=s; rotated(1,2)=-s; rotated(2,2)=c
      call force_theorem_finite_q_hessian_from_eigenbasis_metallic(eigenvalues3, rotated, endpoint_values3, rotated, 0.0_rp, 0.03_rp, &
         torque3, torque_minus3, mixed_q3, hessian, tt, contact, complete)
      value_b = real(complete(1,1),rp); error = abs(value_a-value_b)
   end subroutine test_degenerate_subspace_invariance

   subroutine test_commensurate_endpoint_reuse(error)
      real(rp), intent(out) :: error
      real(rp) :: q(3), q_bad(3), expected(2), actual(2)
      logical :: ok
      mesh = [4,3,2]; nk_mesh = product(mesh)
      allocate(mesh_points(3,nk_mesh),mesh_weights(nk_mesh),periodic_eigenvalues(2,nk_mesh),endpoint_index(nk_mesh))
      ik = 0
      do iz=1,mesh(3); do iy=1,mesh(2); do ix=1,mesh(1)
         ik=ik+1; mesh_points(:,ik)=[(real(ix,rp)-0.5_rp)/real(mesh(1),rp)-0.5_rp, &
            (real(iy,rp)-0.5_rp)/real(mesh(2),rp)-0.5_rp, (real(iz,rp)-0.5_rp)/real(mesh(3),rp)-0.5_rp]
         mesh_weights(ik)=1.0_rp/real(nk_mesh,rp); periodic_eigenvalues(:,ik)=[real(ik,rp),real(ik+100,rp)]
      end do; end do; end do
      error=0.0_rp
      q=[0.0_rp,0.0_rp,0.5_rp]
      call detect_commensurate_endpoint_map(mesh_points,mesh_weights,mesh,q,endpoint_index,ok,residual)
      if (.not.ok .or. residual > 2.0e-12_rp) error=huge(1.0_rp)
      do ik=1,nk_mesh
         expected=periodic_eigenvalues(:,endpoint_index(ik)); actual=expected
         error=max(error,maxval(abs(expected-actual)),abs(mesh_weights(ik)-mesh_weights(endpoint_index(ik))))
      end do
      q=[0.25_rp,1.0_rp/3.0_rp,0.0_rp]
      call detect_commensurate_endpoint_map(mesh_points,mesh_weights,mesh,q,endpoint_index,ok,residual)
      if (.not.ok .or. residual > 2.0e-12_rp) error=huge(1.0_rp)
      do ik=1,nk_mesh
         expected=periodic_eigenvalues(:,endpoint_index(ik)); actual=expected
         error=max(error,maxval(abs(expected-actual)),abs(mesh_weights(ik)-mesh_weights(endpoint_index(ik))))
      end do
      q_bad=[0.0_rp,0.0_rp,0.13_rp]
      call detect_commensurate_endpoint_map(mesh_points,mesh_weights,mesh,q_bad,endpoint_index,ok,residual)
      if (ok) error=huge(1.0_rp)
      deallocate(mesh_points,mesh_weights,periodic_eigenvalues,endpoint_index)
   end subroutine test_commensurate_endpoint_reuse

   pure function oracle_grand_potential(levels,mu,kT) result(value)
      real(rp), intent(in) :: levels(:),mu,kT
      real(rp) :: value, x
      integer :: i
      value=0.0_rp
      do i=1,size(levels)
         x=-(levels(i)-mu)/kT
         if (x > 40.0_rp) then
            value=value-kT*x
         else
            value=value-kT*log(1.0_rp+exp(x))
         end if
      end do
   end function oracle_grand_potential

end program test_exchange_q
