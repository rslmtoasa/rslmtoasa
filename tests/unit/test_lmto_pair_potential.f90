program test_lmto_pair_potential
   use precision_mod, only: rp
   use math_mod, only: i_unit, init_math_operators, hcpx
   use lmto_magnetic_tangent_mod, only: lmto_bond_value, lmto_bond_tangent, lmto_hhmag_to_spinor
   use lmto_pair_potential_mod, only: lmto_circular_pair_potential, lmto_bloch_phase, lmto_endpoint_phases, &
      lmto_circular_pair_potential_from_reverse, lmto_pair_transition_metadata, lmto_transition_metadata, lmto_unfold_site_spinors
   use reciprocal_mod, only: reciprocal
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use charge_mod, only: charge
   use control_mod, only: control
   use basis_mod, only: norb
   use logger_mod, only: g_logger
   implicit none
   real(rp), parameter :: tol = 2.0e-9_rp
   real(rp) :: maximum_error
   logical :: failed

   maximum_error = 0.0_rp; failed = .false.
   call g_logger%init()
   call init_math_operators()
   call test_rotation_oracle_and_adjoint()
   call test_unequal_orbital_negative_control()
   call test_signed_moment_and_bloch_phase()
   call test_finite_q_endpoint_phases_and_gauge()
   call test_commensurate_supercell_oracle()
   call test_reciprocal_service_fixture()
   call test_reciprocal_two_site_identity()
   call test_full_normal_builder_rotation_oracle()
   call test_two_sublattice_service_supercell_oracle()
   write(*,'(a,es12.4)') 'LMTO pair-potential maximum error: ', maximum_error
   if (failed .or. maximum_error > tol) error stop 1
   write(*,'(a)') 'RESULT: PASS'

contains
   subroutine check(error,label)
      real(rp),intent(in)::error; character(len=*),intent(in)::label
      maximum_error=max(maximum_error,error)
      if(error>tol) then; failed=.true.; write(*,'(a,a,a,es12.4)') 'FAIL ',trim(label),': ',error; end if
   end subroutine check

   subroutine test_rotation_oracle_and_adjoint()
      complex(rp)::h(2,2),w0(2),w1(2),c0(2),c1(2),dxh(2,2,4),dyh(2,2,4),vp(2,2,4),vm(2,2,4)
      complex(rp)::dx(4,4),dy(4,4),qminus(4,4),qplus(4,4),fd(4,4),total(4,4),txh(2,2,4),tyh(2,2,4)
      complex(rp)::dx_reverse(4,4),dy_reverse(4,4)
      real(rp)::ez(3),ex(3),ey(3),theta
      logical::supported
      character(len=120)::reason
      call sample(h,w0,w1,c0,c1); ez=[0.0_rp,0.0_rp,1.0_rp]; ex=[1.0_rp,0.0_rp,0.0_rp]; ey=[0.0_rp,1.0_rp,0.0_rp]
      call lmto_bond_tangent(h,w0,w1,w0,w1,c1,ez,ez,ex,[0.0_rp,0.0_rp,0.0_rp],.true.,dxh)
      call lmto_bond_tangent(h,w0,w1,w0,w1,c1,ez,ez,ey,[0.0_rp,0.0_rp,0.0_rp],.true.,dyh)
      call lmto_hhmag_to_spinor(dxh,dx); call lmto_hhmag_to_spinor(dyh,dy)
      call lmto_circular_pair_potential(dx,dy,2.0_rp,qminus,qplus,supported,reason)
      if(.not.supported) then; failed=.true.; return; end if
      dx_reverse=transpose(conjg(dx)); dy_reverse=transpose(conjg(dy))
      call lmto_circular_pair_potential_from_reverse(dx_reverse,dy_reverse,2.0_rp,qplus,supported,reason)
      if(.not.supported) then; failed=.true.; return; end if
      call check(maxval(abs(qplus-transpose(conjg(qminus)))), 'Qplus equals adjoint Qminus')
      theta=5.0e-5_rp
      call lmto_bond_value(h,w0,w1,w0,w1,c0,c1,[sin(theta),0.0_rp,cos(theta)],ez,.true.,vp)
      call lmto_bond_value(h,w0,w1,w0,w1,c0,c1,[sin(-theta),0.0_rp,cos(theta)],ez,.true.,vm)
      call lmto_hhmag_to_spinor((vp-vm)/(2.0_rp*theta),fd)
      call check(maxval(abs(fd-dx)), 'one-site x finite-rotation oracle')
      call lmto_bond_value(h,w0,w1,w0,w1,c0,c1,[0.0_rp,sin(theta),cos(theta)],ez,.true.,vp)
      call lmto_bond_value(h,w0,w1,w0,w1,c0,c1,[0.0_rp,sin(-theta),cos(theta)],ez,.true.,vm)
      call lmto_hhmag_to_spinor((vp-vm)/(2.0_rp*theta),fd)
      call check(maxval(abs(fd-dy)), 'one-site y finite-rotation oracle')
      ! A rigid rotation of both endpoints is the moment-weighted sum of the
      ! separately normalised site pair potentials (M_i=M_j=2 here).
      call lmto_bond_tangent(h,w0,w1,w0,w1,c1,ez,ez,ex,ex,.true.,txh)
      call lmto_bond_tangent(h,w0,w1,w0,w1,c1,ez,ez,ey,ey,.true.,tyh)
      call lmto_hhmag_to_spinor(txh,dx); call lmto_hhmag_to_spinor(tyh,dy)
      call lmto_circular_pair_potential(dx,dy,2.0_rp,qminus,qplus,supported,reason)
      total = 2.0_rp*qminus
      call lmto_bond_value(h,w0,w1,w0,w1,c0,c1,[sin(theta),0.0_rp,cos(theta)], &
         [sin(theta),0.0_rp,cos(theta)],.true.,vp)
      call lmto_bond_value(h,w0,w1,w0,w1,c0,c1,[sin(-theta),0.0_rp,cos(theta)], &
         [sin(-theta),0.0_rp,cos(theta)],.true.,vm)
      call lmto_hhmag_to_spinor((vp-vm)/(2.0_rp*theta),fd)
      ! This is D_x; compare it to the circular relation after rebuilding D_y.
      call check(maxval(abs((total+transpose(conjg(total)))-fd)), 'rigid endpoint moment-weighted identity')
   end subroutine test_rotation_oracle_and_adjoint

   subroutine test_unequal_orbital_negative_control()
      complex(rp)::h(2,2),w0(2),w1(2),c0(2),c1(2),dxh(2,2,4),dyh(2,2,4),dx(4,4),dy(4,4),q(4,4),qp(4,4),scalar(4,4)
      real(rp)::ez(3),ex(3),ey(3); logical::supported; character(len=80)::reason
      call sample(h,w0,w1,c0,c1); ez=[0.0_rp,0.0_rp,1.0_rp]; ex=[1.0_rp,0.0_rp,0.0_rp]; ey=[0.0_rp,1.0_rp,0.0_rp]
      call lmto_bond_tangent(h,w0,w1,w0,w1,c1,ez,ez,ex,[0.0_rp,0.0_rp,0.0_rp],.true.,dxh)
      call lmto_bond_tangent(h,w0,w1,w0,w1,c1,ez,ez,ey,[0.0_rp,0.0_rp,0.0_rp],.true.,dyh)
      call lmto_hhmag_to_spinor(dxh,dx); call lmto_hhmag_to_spinor(dyh,dy)
      call lmto_circular_pair_potential(dx,dy,1.0_rp,q,qp,supported,reason)
      scalar=cmplx(0.0_rp,0.0_rp,rp); scalar(3,1)=cmplx(0.25_rp,0.0_rp,rp); scalar(4,2)=scalar(3,1)
      if(maxval(abs(q-scalar))<1.0e-3_rp) then; failed=.true.; write(*,'(a)') 'FAIL unequal-orbital operator collapsed to site scalar'; end if
   end subroutine test_unequal_orbital_negative_control

   subroutine test_signed_moment_and_bloch_phase()
      complex(rp)::dx(2,2),dy(2,2),qm(2,2),qp(2,2),phase
      logical::supported; character(len=80)::reason
      dx=cmplx(0.0_rp,0.0_rp,rp); dy=dx; dx(1,2)=cmplx(2.0_rp,0.0_rp,rp)
      call lmto_circular_pair_potential(dx,dy,-2.0_rp,qm,qp,supported,reason)
      call check(abs(qm(1,2)+0.5_rp), 'signed reversed-sublattice normalization')
      phase=lmto_bloch_phase([0.0_rp,0.0_rp,0.0_rp],[2.0_rp,0.0_rp,0.0_rp])
      call check(abs(phase-cmplx(1.0_rp,0.0_rp,rp)), 'q=0 Bloch phase convention')
      phase=lmto_bloch_phase([0.25_rp,0.0_rp,0.0_rp],[2.0_rp,0.0_rp,0.0_rp])
      call check(abs(phase-cmplx(-1.0_rp,0.0_rp,rp)), 'q=0-compatible Bloch phase convention')
      call lmto_circular_pair_potential(dx,dy,0.0_rp,qm,qp,supported,reason)
      if(supported) then; failed=.true.; write(*,'(a)') 'FAIL zero signed moment was accepted'; end if
   end subroutine test_signed_moment_and_bloch_phase

   subroutine test_finite_q_endpoint_phases_and_gauge()
      real(rp)::k(3),q(3),d(3),tau(3,2)
      complex(rp)::left,right,common,folded(4,1),unfolded(4,1),expected(4,1),ket(4,1),qmat(4,4),qfold(4,4),vertex_unfold,vertex_fold,ug(4)
      type(lmto_pair_transition_metadata)::meta
      k=[0.17_rp,0.0_rp,0.0_rp]; q=[0.5_rp,0.0_rp,0.0_rp]; d=[1.0_rp,0.0_rp,0.0_rp]
      call lmto_endpoint_phases(k,q,d,left,right)
      call check(abs(left-lmto_bloch_phase(k,d)), 'left endpoint carries k phase')
      call check(abs(right-lmto_bloch_phase(k+q,d)), 'right endpoint carries k+q phase')
      common=left
      if(abs(right-common)<1.0e-3_rp) then; failed=.true.; write(*,'(a)') 'FAIL common left/right phase negative control'; end if
      q=[1.0_rp/3.0_rp,0.0_rp,0.0_rp]
      call lmto_endpoint_phases(k,q,d,left,right)
      if(abs(right-left)<1.0e-3_rp) then; failed=.true.; write(*,'(a)') 'FAIL q=1/3 endpoint phase is a no-op'; end if
      meta=lmto_transition_metadata([0.45_rp,0.0_rp,0.0_rp],[0.2_rp,0.0_rp,0.0_rp])
      if(.not.meta%unfolded_gauge_required .or. meta%reciprocal_shift(1)/=1) then; failed=.true.; write(*,'(a)') 'FAIL folded/unfolded metadata'; end if
      tau(:,1)=[0.0_rp,0.0_rp,0.0_rp]; tau(:,2)=[0.25_rp,0.0_rp,0.0_rp]
      folded=cmplx(1.0_rp,0.0_rp,rp); call lmto_unfold_site_spinors(meta,tau,folded,unfolded)
      expected=folded; expected(3:4,:)=exp(cmplx(0.0_rp,-0.5_rp*acos(-1.0_rp),rp))*expected(3:4,:)
      call check(maxval(abs(unfolded-expected)), 'multisublattice unfolded site gauge')
      ! The physical contraction is invariant only after transforming the
      ! folded target state or, equivalently, the vertex.  This is deliberately
      ! a vertex check rather than a false raw-matrix periodicity check.
      ket(:,1)=[cmplx(0.3_rp,0.1_rp,rp),cmplx(-0.2_rp,0.4_rp,rp),cmplx(0.5_rp,-0.1_rp,rp),cmplx(0.1_rp,0.2_rp,rp)]
      qmat=cmplx(0.0_rp,0.0_rp,rp); qmat(1,1)=cmplx(0.2_rp,-0.1_rp,rp); qmat(3,2)=cmplx(-0.3_rp,0.4_rp,rp); qmat(4,4)=cmplx(0.1_rp,0.2_rp,rp)
      ug=[cmplx(1.0_rp,0.0_rp,rp),cmplx(1.0_rp,0.0_rp,rp),exp(cmplx(0.0_rp,0.5_rp*acos(-1.0_rp),rp)),exp(cmplx(0.0_rp,0.5_rp*acos(-1.0_rp),rp))]
      qfold=qmat; qfold(1,:)=ug(1)*qmat(1,:); qfold(2,:)=ug(2)*qmat(2,:); qfold(3,:)=ug(3)*qmat(3,:); qfold(4,:)=ug(4)*qmat(4,:)
      vertex_unfold=sum(conjg(unfolded(:,1))*matmul(qmat,ket(:,1)))
      vertex_fold=sum(conjg(folded(:,1))*matmul(qfold,ket(:,1)))
      call check(abs(vertex_unfold-vertex_fold), 'folded/unfolded physical transition vertex')
   end subroutine test_finite_q_endpoint_phases_and_gauge

   ! Exercise the production assembler with a completed minimal ordinary-LMTO
   ! object.  Geometry is supplied through the same direct-displacement cache
   ! used by normal reciprocal Fourier assembly.
   subroutine test_reciprocal_service_fixture()
      type(reciprocal) :: recip
      type(hamiltonian), target :: ham
      type(lattice), target :: lat
      type(charge), target :: chg
      type(control), target :: ctl
      integer :: nmat
      real(rp) :: k(3), q0(3), qhalf(3), qthird(3)
      complex(rp), allocatable :: q_default(:,:), q_zero(:,:), q_half(:,:), q_third(:,:), qplus(:,:)
      logical :: supported
      character(len=160) :: reason

      call setup_reciprocal_service_fixture(recip,ham,lat,chg,ctl)
      nmat=2*norb
      allocate(q_default(nmat,nmat),q_zero(nmat,nmat),q_half(nmat,nmat),q_third(nmat,nmat),qplus(nmat,nmat))
      k=[0.17_rp,0.0_rp,0.0_rp]; q0=0.0_rp; qhalf=[0.5_rp,0.0_rp,0.0_rp]; qthird=[1.0_rp/3.0_rp,0.0_rp,0.0_rp]
      call recip%build_lmto_pair_potential_at_kpoint(1,k,2.0_rp,q_default,qplus,supported,reason)
      if(.not.supported) then; failed=.true.; write(*,'(a,a)') 'FAIL reciprocal service default: ',trim(reason); return; end if
      call recip%build_lmto_pair_potential_at_kpoint(1,k,2.0_rp,q_zero,qplus,supported,reason,q0)
      if(.not.supported) then; failed=.true.; write(*,'(a,a)') 'FAIL reciprocal service q=0: ',trim(reason); return; end if
      call check(maxval(abs(q_default-q_zero)), 'reciprocal service q=0 reduction')
      call recip%build_lmto_pair_potential_at_kpoint(1,k,2.0_rp,q_half,qplus,supported,reason,qhalf)
      if(.not.supported) then; failed=.true.; write(*,'(a,a)') 'FAIL reciprocal service q=1/2: ',trim(reason); return; end if
      call check(maxval(abs(qplus-transpose(conjg(q_half)))), 'reciprocal service reverse Qplus q=1/2')
      call expected_reciprocal_fixture_operator(k,qhalf,q_default)
      call check(maxval(abs(q_half-q_default)), 'reciprocal service analytic q=1/2 operator')
      call recip%build_lmto_pair_potential_at_kpoint(1,k,2.0_rp,q_third,qplus,supported,reason,qthird)
      if(.not.supported) then; failed=.true.; write(*,'(a,a)') 'FAIL reciprocal service q=1/3: ',trim(reason); return; end if
      call check(maxval(abs(qplus-transpose(conjg(q_third)))), 'reciprocal service reverse Qplus q=1/3')
      call expected_reciprocal_fixture_operator(k,qthird,q_default)
      call check(maxval(abs(q_third-q_default)), 'reciprocal service analytic q=1/3 operator')
      if(maxval(abs(q_half-q_zero)) < 1.0e-6_rp .or. maxval(abs(q_third-q_zero)) < 1.0e-6_rp) then
         failed=.true.; write(*,'(a)') 'FAIL finite-q reciprocal service did not retain endpoint phases'
      end if
      call recip%build_lmto_pair_potential_at_kpoint(1,k,0.0_rp,q_half,qplus,supported,reason,qhalf)
      if(supported) then; failed=.true.; write(*,'(a)') 'FAIL reciprocal service accepted zero signed moment'; end if
      recip%reciprocal_mode='generalized_overlap_proxy'
      call recip%build_lmto_pair_potential_at_kpoint(1,k,2.0_rp,q_half,qplus,supported,reason,qhalf)
      if(supported) then; failed=.true.; write(*,'(a)') 'FAIL reciprocal service accepted overlap proxy'; end if
      recip%reciprocal_mode='ham_only'; ham%hoh=.true.
      call recip%build_lmto_pair_potential_at_kpoint(1,k,2.0_rp,q_half,qplus,supported,reason,qhalf)
      if(supported) then; failed=.true.; write(*,'(a)') 'FAIL reciprocal service accepted HOH'; end if
      ham%hoh=.false.; ham%magnetic_representation='gbt_single_q'
      call recip%build_lmto_pair_potential_at_kpoint(1,k,2.0_rp,q_half,qplus,supported,reason,qhalf)
      if(supported) then; failed=.true.; write(*,'(a)') 'FAIL reciprocal service accepted GBT'; end if
      ham%magnetic_representation='periodic_nc'; ham%ccor_2c=.true.
      call recip%build_lmto_pair_potential_at_kpoint(1,k,2.0_rp,q_half,qplus,supported,reason,qhalf)
      if(supported) then; failed=.true.; write(*,'(a)') 'FAIL reciprocal service accepted CCOR'; end if
      ham%ccor_2c=.false.; allocate(ham%lsham(norb,norb,1)); ham%lsham=cmplx(0.0_rp,0.0_rp,rp); ham%lsham(1,1,1)=cmplx(1.0_rp,0.0_rp,rp)
      call recip%build_lmto_pair_potential_at_kpoint(1,k,2.0_rp,q_half,qplus,supported,reason,qhalf)
      if(supported) then; failed=.true.; write(*,'(a)') 'FAIL reciprocal service accepted SOC'; end if
      deallocate(ham%lsham)
   end subroutine test_reciprocal_service_fixture

   subroutine setup_reciprocal_service_fixture(recip,ham,lat,chg,ctl)
      type(reciprocal), intent(out) :: recip
      type(hamiltonian), target, intent(out) :: ham
      type(lattice), target, intent(out) :: lat
      type(charge), target, intent(out) :: chg
      type(control), target, intent(out) :: ctl
      integer :: i
      call ctl%restore_to_default()
      lat%nrec=1; lat%ntype=1; lat%nn_max=3; lat%kk=1; lat%nmax=1
      allocate(lat%ib(1),lat%atlist(1),lat%iz(1),lat%num(1),lat%nn(1,3),lat%sbar(norb,norb,3,1),lat%symbolic_atoms(1))
      lat%ib=1; lat%atlist=1; lat%iz=1; lat%num=1; lat%nn(1,:)=[3,1,1]; lat%sbar=cmplx(0.0_rp,0.0_rp,rp)
      call lat%symbolic_atoms(1)%restore_to_default()
      lat%symbolic_atoms(1)%potential%mom=[0.0_rp,0.0_rp,1.0_rp]
      lat%symbolic_atoms(1)%potential%wx0=cmplx(1.0_rp,0.0_rp,rp)
      lat%symbolic_atoms(1)%potential%wx1=cmplx(0.31_rp,0.0_rp,rp)
      lat%symbolic_atoms(1)%potential%cx1=cmplx(0.23_rp,0.0_rp,rp)
      do i=1,norb
         lat%sbar(i,i,1,1)=cmplx(0.12_rp,0.0_rp,rp); lat%sbar(i,i,2,1)=cmplx(0.37_rp,0.0_rp,rp); lat%sbar(i,i,3,1)=cmplx(0.37_rp,0.0_rp,rp)
      end do
      chg%lattice=>lat; chg%symbolic_atom=>lat%symbolic_atoms
      ham%charge=>chg; ham%lattice=>lat; ham%control=>ctl; ham%hoh=.false.; ham%ccor_2c=.false.; ham%hubbard_u_general_check=.false.; ham%hubbard_v_check=.false.; ham%hubbard_u_impurity_check=.false.; ham%local_axis=.false.; ham%magnetic_representation='periodic_nc'; ham%operator_generation=1
      recip%hamiltonian=>ham; recip%lattice=>lat; recip%control=>ctl; recip%reciprocal_mode='ham_only'; recip%max_orbs=2*norb
      allocate(recip%ham_vec_type(3,3,1),recip%ham_vec_type_direct(3,3,1)); recip%ham_vec_type=0.0_rp; recip%ham_vec_type_direct=0.0_rp
      recip%ham_vec_type(1,2,1)=1.0_rp; recip%ham_vec_type_direct(1,2,1)=1.0_rp; recip%ham_vec_type(1,3,1)=-1.0_rp; recip%ham_vec_type_direct(1,3,1)=-1.0_rp
   end subroutine setup_reciprocal_service_fixture

   subroutine test_reciprocal_two_site_identity()
      type(reciprocal) :: recip
      type(hamiltonian), target :: ham
      type(lattice), target :: lat
      type(charge), target :: chg
      type(control), target :: ctl
      integer :: nblock,nmat
      complex(rp), allocatable :: qa(:,:),qb(:,:),qplus(:,:)
      logical :: supported
      character(len=160) :: reason
      call setup_two_site_reciprocal_fixture(recip,ham,lat,chg,ctl)
      nblock=2*norb; nmat=2*nblock; allocate(qa(nmat,nmat),qb(nmat,nmat),qplus(nmat,nmat))
      call recip%build_lmto_pair_potential_at_kpoint(1,[0.17_rp,0.0_rp,0.0_rp],2.0_rp,qa,qplus,supported,reason,[0.5_rp,0.0_rp,0.0_rp])
      if(.not.supported) then; failed=.true.; write(*,'(a,a)') 'FAIL two-site response a: ',trim(reason); return; end if
      call check(maxval(abs(qplus-transpose(conjg(qa)))), 'two-site response a reverse Qplus')
      call recip%build_lmto_pair_potential_at_kpoint(2,[0.17_rp,0.0_rp,0.0_rp],-2.0_rp,qb,qplus,supported,reason,[0.5_rp,0.0_rp,0.0_rp])
      if(.not.supported) then; failed=.true.; write(*,'(a,a)') 'FAIL two-site response b: ',trim(reason); return; end if
      call check(maxval(abs(qa(nblock+1:nmat,nblock+1:nmat))), 'same-type site a does not leak to site b onsite block')
      call check(maxval(abs(qb(1:nblock,1:nblock))), 'same-type site b does not leak to site a onsite block')
   end subroutine test_reciprocal_two_site_identity

   subroutine setup_two_site_reciprocal_fixture(recip,ham,lat,chg,ctl)
      type(reciprocal), intent(out) :: recip
      type(hamiltonian), target, intent(out) :: ham
      type(lattice), target, intent(out) :: lat
      type(charge), target, intent(out) :: chg
      type(control), target, intent(out) :: ctl
      integer :: isite,i
      call ctl%restore_to_default()
      lat%nrec=2; lat%ntype=2; lat%nn_max=3; lat%kk=2; lat%nmax=2
      allocate(lat%ib(2),lat%atlist(2),lat%iz(2),lat%num(2),lat%nn(2,3),lat%sbar(norb,norb,3,2),lat%symbolic_atoms(2))
      lat%ib=[1,2]; lat%atlist=[1,2]; lat%iz=[1,2]; lat%num=[1,2]; lat%nn(1,:)=[3,2,2]; lat%nn(2,:)=[3,1,1]; lat%sbar=cmplx(0.0_rp,0.0_rp,rp)
      do isite=1,2
         call lat%symbolic_atoms(isite)%restore_to_default(); lat%symbolic_atoms(isite)%potential%mom=[0.0_rp,0.0_rp,1.0_rp]
         lat%symbolic_atoms(isite)%potential%wx0=cmplx(1.0_rp,0.0_rp,rp); lat%symbolic_atoms(isite)%potential%wx1=cmplx(0.31_rp,0.0_rp,rp); lat%symbolic_atoms(isite)%potential%cx1=cmplx(0.23_rp,0.0_rp,rp)
         do i=1,norb
            lat%sbar(i,i,1,isite)=cmplx(0.12_rp,0.0_rp,rp); lat%sbar(i,i,2,isite)=cmplx(0.37_rp,0.0_rp,rp); lat%sbar(i,i,3,isite)=cmplx(0.37_rp,0.0_rp,rp)
         end do
      end do
      chg%lattice=>lat; chg%symbolic_atom=>lat%symbolic_atoms
      ham%charge=>chg; ham%lattice=>lat; ham%control=>ctl; ham%hoh=.false.; ham%ccor_2c=.false.; ham%hubbard_u_general_check=.false.; ham%hubbard_v_check=.false.; ham%hubbard_u_impurity_check=.false.; ham%local_axis=.false.; ham%magnetic_representation='periodic_nc'; ham%operator_generation=1
      recip%hamiltonian=>ham; recip%lattice=>lat; recip%control=>ctl; recip%reciprocal_mode='ham_only'; recip%max_orbs=2*norb
      allocate(recip%ham_vec_type(3,3,2),recip%ham_vec_type_direct(3,3,2)); recip%ham_vec_type=0.0_rp; recip%ham_vec_type_direct=0.0_rp
      ! Full directed displacements for tau_A=0, tau_B=1/4.  The second
      ! neighbour of each type is in the adjacent primitive cell, which makes
      ! this a genuine two-sublattice periodic fixture rather than a pair of
      ! duplicated intra-cell bonds.
      recip%ham_vec_type(1,2,1)=0.25_rp; recip%ham_vec_type_direct(1,2,1)=0.25_rp
      recip%ham_vec_type(1,3,1)=-0.75_rp; recip%ham_vec_type_direct(1,3,1)=-0.75_rp
      recip%ham_vec_type(1,2,2)=0.75_rp; recip%ham_vec_type_direct(1,2,2)=0.75_rp
      recip%ham_vec_type(1,3,2)=-0.25_rp; recip%ham_vec_type_direct(1,3,2)=-0.25_rp
   end subroutine setup_two_site_reciprocal_fixture

   ! WR-02 q=0 closure: the finite rotations rebuild every directed block via
   ! the production normal Hamiltonian builder, ham0m_nc.  This deliberately
   ! does not use lmto_bond_value, endpoint tangents, or a bond-only oracle.
   subroutine test_full_normal_builder_rotation_oracle()
      type(reciprocal) :: recip
      type(hamiltonian), target :: ham
      type(lattice), target :: lat
      type(charge), target :: chg
      type(control), target :: ctl
      integer :: nmat
      real(rp) :: theta, signed_moment(2), moments_plus(3,2), moments_minus(3,2)
      complex(rp), allocatable :: qminus(:, :, :), qplus(:, :, :), hplus(:,:), hminus(:,:), fdx(:,:), fdy(:,:), qfd(:,:), &
         hzero(:,:), spin_rotation(:,:), predicted(:,:)
      logical :: supported
      character(len=160) :: reason
      integer :: isite

      call setup_two_site_reciprocal_fixture(recip,ham,lat,chg,ctl)
      nmat=4*norb; theta=1.0e-8_rp; signed_moment=[2.0_rp,-2.0_rp]
      allocate(qminus(nmat,nmat,2),qplus(nmat,nmat,2),hplus(nmat,nmat),hminus(nmat,nmat),fdx(nmat,nmat),fdy(nmat,nmat), &
         qfd(nmat,nmat),hzero(nmat,nmat),spin_rotation(nmat,nmat),predicted(nmat,nmat))
      do isite=1,2
         call recip%build_lmto_pair_potential_at_kpoint(isite,[0.0_rp,0.0_rp,0.0_rp],signed_moment(isite), &
            qminus(:,:,isite),qplus(:,:,isite),supported,reason,[0.0_rp,0.0_rp,0.0_rp])
         if(.not.supported) then; failed=.true.; write(*,'(a,a)') 'FAIL full-builder service: ',trim(reason); return; end if
         moments_plus=0.0_rp; moments_plus(3,:)=1.0_rp; moments_minus=moments_plus
         moments_plus(1,isite)=sin(theta); moments_plus(3,isite)=cos(theta)
         moments_minus(1,isite)=-sin(theta); moments_minus(3,isite)=cos(theta)
         call build_two_site_normal_hamiltonian(ham,moments_plus,hplus)
         call build_two_site_normal_hamiltonian(ham,moments_minus,hminus)
         fdx=(hplus-hminus)/(2.0_rp*theta)
         moments_plus=0.0_rp; moments_plus(3,:)=1.0_rp; moments_minus=moments_plus
         moments_plus(2,isite)=sin(theta); moments_plus(3,isite)=cos(theta)
         moments_minus(2,isite)=-sin(theta); moments_minus(3,isite)=cos(theta)
         call build_two_site_normal_hamiltonian(ham,moments_plus,hplus)
         call build_two_site_normal_hamiltonian(ham,moments_minus,hminus)
         fdy=(hplus-hminus)/(2.0_rp*theta)
         qfd=(fdx-i_unit*fdy)/(2.0_rp*signed_moment(isite))
         call check(maxval(abs(qfd-qminus(:,:,isite))), 'q=0 full normal-builder finite rotation')
      end do
      moments_plus=0.0_rp; moments_plus(3,:)=cos(theta); moments_plus(1,:)=sin(theta)
      moments_minus=0.0_rp; moments_minus(3,:)=1.0_rp
      call build_two_site_normal_hamiltonian(ham,moments_plus,hplus)
      call build_two_site_normal_hamiltonian(ham,moments_minus,hzero)
      call spin_rotation_about_y(theta,2*norb,spin_rotation)
      predicted=matmul(spin_rotation,matmul(hzero,transpose(conjg(spin_rotation))))
      call check(maxval(abs(hplus-predicted)), 'q=0 full normal-builder rigid rotation covariance')
      call ham%clear_texture_moments()
      ham%magnetic_representation='periodic_nc'
   end subroutine test_full_normal_builder_rotation_oracle

   subroutine spin_rotation_about_y(theta,nblock,rotation)
      real(rp), intent(in) :: theta
      integer, intent(in) :: nblock
      complex(rp), intent(out) :: rotation(:,:)
      integer :: isite,ibeg,iend,iorb
      real(rp) :: c,s
      c=cos(0.5_rp*theta); s=sin(0.5_rp*theta); rotation=cmplx(0.0_rp,0.0_rp,rp)
      do isite=1,2
         ibeg=(isite-1)*nblock+1; iend=ibeg+nblock/2-1
         do iorb=0,nblock/2-1
            rotation(ibeg+iorb,ibeg+iorb)=cmplx(c,0.0_rp,rp)
            rotation(iend+1+iorb,iend+1+iorb)=cmplx(c,0.0_rp,rp)
            rotation(ibeg+iorb,iend+1+iorb)=cmplx(-s,0.0_rp,rp)
            rotation(iend+1+iorb,ibeg+iorb)=cmplx(s,0.0_rp,rp)
         end do
      end do
   end subroutine spin_rotation_about_y

   ! WR-01c closure: compare the actual reciprocal service directly with an
   ! explicit N-cell, two-sublattice finite-rotation supercell.  The oracle
   ! invokes only the normal ham0m_nc builder and its own absolute-position
   ! Fourier projection; it cannot consume the analytic tangent or Q service.
   subroutine test_two_sublattice_service_supercell_oracle()
      call test_two_sublattice_service_supercell_q(2,1)
      call test_two_sublattice_service_supercell_q(2,2)
      call test_two_sublattice_service_supercell_q(3,1)
      call test_two_sublattice_service_supercell_q(3,2)
   end subroutine test_two_sublattice_service_supercell_oracle

   subroutine test_two_sublattice_service_supercell_q(ncell,response_site)
      integer, intent(in) :: ncell,response_site
      type(reciprocal) :: recip
      type(hamiltonian), target :: ham
      type(lattice), target :: lat
      type(charge), target :: chg
      type(control), target :: ctl
      integer :: nmat,itheta
      real(rp) :: q(3),theta(3),signed_moment,error(3),relative_error
      complex(rp), allocatable :: qminus(:,:),qplus(:,:),oracle(:,:)
      logical :: supported
      character(len=160) :: reason

      call setup_two_site_reciprocal_fixture(recip,ham,lat,chg,ctl)
      nmat=4*norb; q=[1.0_rp/real(ncell,rp),0.0_rp,0.0_rp]
      signed_moment=merge(2.0_rp,-2.0_rp,response_site==1)
      allocate(qminus(nmat,nmat),qplus(nmat,nmat),oracle(nmat,nmat))
      call recip%build_lmto_pair_potential_at_kpoint(response_site,[0.0_rp,0.0_rp,0.0_rp],signed_moment, &
         qminus,qplus,supported,reason,q)
      if(.not.supported) then; failed=.true.; write(*,'(a,a)') 'FAIL two-sublattice reciprocal service: ',trim(reason); return; end if
      theta=[2.5e-4_rp,1.25e-4_rp,6.25e-5_rp]
      do itheta=1,3
         call two_sublattice_supercell_qminus(ham,ncell,q(1),response_site,signed_moment,theta(itheta),oracle)
         error(itheta)=maxval(abs(oracle-qminus))
      end do
      relative_error=error(3)/max(maxval(abs(qminus)),tiny(1.0_rp))
      write(*,'(a,i0,a,i0,a,3(1x,es10.3),a,1x,es10.3)') 'Two-sublattice service/supercell q=1/',ncell, &
         ' site=',response_site,' Q errors:',error,' relative:',relative_error
      call check(error(3), 'two-sublattice reciprocal-service supercell oracle')
      call check(maxval(abs(qplus-transpose(conjg(oracle)))), 'two-sublattice independently accumulated reverse Qplus')
      if (error(1) < 3.5_rp*error(2) .or. error(2) < 3.0_rp*error(3)) then
         failed=.true.; write(*,'(a,i0,a,i0)') 'FAIL two-sublattice central difference is not O(theta^2), q denominator ',ncell,' site ',response_site
      end if
      call ham%clear_texture_moments()
      ham%magnetic_representation='periodic_nc'
   end subroutine test_two_sublattice_service_supercell_q

   subroutine build_two_site_normal_hamiltonian(ham,moments,hk)
      type(hamiltonian), intent(inout) :: ham
      real(rp), intent(in) :: moments(3,2)
      complex(rp), intent(out) :: hk(:,:)
      integer :: isite,ineigh,ja,it,jt,ibeg,iend,jbeg,jend,ilm,jlm
      real(rp) :: vet(3),hhh(norb,norb)
      complex(rp) :: cart(norb,norb,4),spinor(2*norb,2*norb)

      hk=cmplx(0.0_rp,0.0_rp,rp)
      if (.not. allocated(ham%hhmag)) allocate(ham%hhmag(norb,norb,4))
      ham%magnetic_representation='explicit_texture'
      call ham%set_texture_moments(moments)
      do isite=1,2
         it=isite; ibeg=(isite-1)*2*norb+1; iend=isite*2*norb
         do ineigh=1,3
            ja=isite; if (ineigh/=1) ja=ham%lattice%nn(isite,ineigh)
            jt=ham%lattice%iz(ja)
            jbeg=(ja-1)*2*norb+1; jend=ja*2*norb
            call two_site_bond_displacement(isite,ineigh,vet)
            do jlm=1,norb; do ilm=1,norb
               hhh(ilm,jlm)=real(ham%lattice%sbar(jlm,ilm,ineigh,isite),rp)
            end do; end do
            call ham%ham0m_nc(isite,ja,it,jt,vet,hhh)
            cart=ham%hhmag
            call hcpx(cart(:,:,1),'cart2sph'); call hcpx(cart(:,:,2),'cart2sph')
            call hcpx(cart(:,:,3),'cart2sph'); call hcpx(cart(:,:,4),'cart2sph')
            call lmto_hhmag_to_spinor(cart,spinor)
            hk(ibeg:iend,jbeg:jend)=hk(ibeg:iend,jbeg:jend)+spinor
         end do
      end do
   end subroutine build_two_site_normal_hamiltonian

   subroutine two_sublattice_supercell_qminus(ham,ncell,q,response_site,signed_moment,theta,qminus)
      type(hamiltonian), intent(inout) :: ham
      integer, intent(in) :: ncell,response_site
      real(rp), intent(in) :: q,signed_moment,theta
      complex(rp), intent(out) :: qminus(:,:)
      complex(rp) :: dxcos(size(qminus,1),size(qminus,2)),dxsin(size(qminus,1),size(qminus,2))
      complex(rp) :: dycos(size(qminus,1),size(qminus,2)),dysin(size(qminus,1),size(qminus,2))
      complex(rp) :: dx(size(qminus,1),size(qminus,2)),dy(size(qminus,1),size(qminus,2))
      call two_sublattice_supercell_derivative(ham,ncell,q,response_site,1,.true.,theta,dxcos)
      call two_sublattice_supercell_derivative(ham,ncell,q,response_site,1,.false.,theta,dxsin)
      call two_sublattice_supercell_derivative(ham,ncell,q,response_site,2,.true.,theta,dycos)
      call two_sublattice_supercell_derivative(ham,ncell,q,response_site,2,.false.,theta,dysin)
      dx=dxcos+i_unit*dxsin; dy=dycos+i_unit*dysin
      qminus=(dx-i_unit*dy)/(2.0_rp*signed_moment)
   end subroutine two_sublattice_supercell_qminus

   subroutine two_sublattice_supercell_derivative(ham,ncell,q,response_site,component,is_cosine,theta,derivative)
      type(hamiltonian), intent(inout) :: ham
      integer, intent(in) :: ncell,response_site,component
      real(rp), intent(in) :: q,theta
      logical, intent(in) :: is_cosine
      complex(rp), intent(out) :: derivative(:,:)
      integer :: nsuper
      real(rp), allocatable :: moments_plus(:,:),moments_minus(:,:)
      complex(rp), allocatable :: hplus(:,:),hminus(:,:)
      nsuper=2*ncell
      allocate(moments_plus(3,nsuper),moments_minus(3,nsuper),hplus(2*norb*nsuper,2*norb*nsuper),hminus(2*norb*nsuper,2*norb*nsuper))
      call make_two_sublattice_texture(ncell,q,response_site,component,is_cosine,theta,moments_plus)
      call make_two_sublattice_texture(ncell,q,response_site,component,is_cosine,-theta,moments_minus)
      call build_two_sublattice_supercell_hamiltonian(ham,ncell,moments_plus,hplus)
      call build_two_sublattice_supercell_hamiltonian(ham,ncell,moments_minus,hminus)
      call project_two_sublattice_supercell(ncell,q,(hplus-hminus)/(2.0_rp*theta),derivative)
   end subroutine two_sublattice_supercell_derivative

   subroutine make_two_sublattice_texture(ncell,q,response_site,component,is_cosine,theta,moments)
      integer, intent(in) :: ncell,response_site,component
      real(rp), intent(in) :: q,theta
      logical, intent(in) :: is_cosine
      real(rp), intent(out) :: moments(:,:)
      integer :: icell,isite
      real(rp) :: tau(2),phase,amplitude
      tau=[0.0_rp,0.25_rp]; moments=0.0_rp; moments(3,:)=1.0_rp
      do icell=1,ncell
         isite=2*(icell-1)+response_site
         phase=2.0_rp*acos(-1.0_rp)*q*(real(icell-1,rp)+tau(response_site))
         amplitude=theta*merge(cos(phase),sin(phase),is_cosine)
         moments(component,isite)=amplitude
         moments(:,isite)=moments(:,isite)/norm2(moments(:,isite))
      end do
   end subroutine make_two_sublattice_texture

   subroutine build_two_sublattice_supercell_hamiltonian(ham,ncell,moments,hk)
      type(hamiltonian), intent(inout) :: ham
      integer, intent(in) :: ncell
      real(rp), intent(in) :: moments(:,:)
      complex(rp), intent(out) :: hk(:,:)
      integer :: icell,isite,ineigh,jcell,ja,it,jt,ip,jp,ibeg,iend,jbeg,jend,ilm,jlm
      real(rp) :: vet(3),hhh(norb,norb)
      complex(rp) :: cart(norb,norb,4),spinor(2*norb,2*norb)
      hk=cmplx(0.0_rp,0.0_rp,rp)
      if (.not. allocated(ham%hhmag)) allocate(ham%hhmag(norb,norb,4))
      ham%magnetic_representation='explicit_texture'; call ham%set_texture_moments(moments)
      do icell=1,ncell
         do isite=1,2
            ip=2*(icell-1)+isite; it=isite; ibeg=(ip-1)*2*norb+1; iend=ip*2*norb
            do ineigh=1,3
               call two_site_supercell_target(ncell,icell,isite,ineigh,jcell,ja)
               jp=2*(jcell-1)+ja; jt=ja; jbeg=(jp-1)*2*norb+1; jend=jp*2*norb
               call two_site_bond_displacement(isite,ineigh,vet)
               do jlm=1,norb; do ilm=1,norb
                  hhh(ilm,jlm)=real(ham%lattice%sbar(jlm,ilm,ineigh,isite),rp)
               end do; end do
               call ham%ham0m_nc(ip,jp,it,jt,vet,hhh)
               cart=ham%hhmag
               call hcpx(cart(:,:,1),'cart2sph'); call hcpx(cart(:,:,2),'cart2sph')
               call hcpx(cart(:,:,3),'cart2sph'); call hcpx(cart(:,:,4),'cart2sph')
               call lmto_hhmag_to_spinor(cart,spinor)
               hk(ibeg:iend,jbeg:jend)=hk(ibeg:iend,jbeg:jend)+spinor
            end do
         end do
      end do
   end subroutine build_two_sublattice_supercell_hamiltonian

   subroutine two_site_supercell_target(ncell,icell,isite,ineigh,jcell,ja)
      integer, intent(in) :: ncell,icell,isite,ineigh
      integer, intent(out) :: jcell,ja
      ja=isite; jcell=icell
      if (ineigh==1) return
      if (isite==1) then
         ja=2; if(ineigh==3) jcell=modulo(icell-2,ncell)+1
      else
         ja=1; if(ineigh==2) jcell=modulo(icell,ncell)+1
      end if
   end subroutine two_site_supercell_target

   subroutine two_site_bond_displacement(isite,ineigh,vet)
      integer, intent(in) :: isite,ineigh
      real(rp), intent(out) :: vet(3)
      vet=0.0_rp
      if(isite==1 .and. ineigh==2) vet(1)=0.25_rp
      if(isite==1 .and. ineigh==3) vet(1)=-0.75_rp
      if(isite==2 .and. ineigh==2) vet(1)=0.75_rp
      if(isite==2 .and. ineigh==3) vet(1)=-0.25_rp
   end subroutine two_site_bond_displacement

   subroutine project_two_sublattice_supercell(ncell,q,delta_h,derivative)
      integer, intent(in) :: ncell
      real(rp), intent(in) :: q
      complex(rp), intent(in) :: delta_h(:,:)
      complex(rp), intent(out) :: derivative(:,:)
      integer :: icell,jcell,isite,jsite,ip,jp,ibeg,iend,jbeg,jend
      real(rp) :: tau(2),ri,rj
      complex(rp) :: bra,ket
      tau=[0.0_rp,0.25_rp]; derivative=cmplx(0.0_rp,0.0_rp,rp)
      do jcell=1,ncell; do jsite=1,2
         jp=2*(jcell-1)+jsite; jbeg=(jp-1)*2*norb+1; jend=jp*2*norb
         rj=real(jcell-1,rp)+tau(jsite); ket=cmplx(1.0_rp,0.0_rp,rp)
         do icell=1,ncell; do isite=1,2
            ip=2*(icell-1)+isite; ibeg=(ip-1)*2*norb+1; iend=ip*2*norb
            ri=real(icell-1,rp)+tau(isite)
            bra=exp(cmplx(0.0_rp,-2.0_rp*acos(-1.0_rp)*q*ri,rp))
            derivative((isite-1)*2*norb+1:isite*2*norb,(jsite-1)*2*norb+1:jsite*2*norb)= &
               derivative((isite-1)*2*norb+1:isite*2*norb,(jsite-1)*2*norb+1:jsite*2*norb)+ &
               bra*delta_h(ibeg:iend,jbeg:jend)*ket/real(ncell,rp)
         end do; end do
      end do; end do
   end subroutine project_two_sublattice_supercell

   ! This is an explicit bond-by-bond oracle for the completed fixture.  It
   ! does not call the reciprocal service or consume endpoint records.
   subroutine expected_reciprocal_fixture_operator(k,q,qminus)
      real(rp), intent(in) :: k(3),q(3)
      complex(rp), intent(out) :: qminus(:,:)
      complex(rp) :: h(norb,norb), w0(norb), w1(norb), c1(norb), tangent(norb,norb,4), cart(norb,norb,4), spinor(2*norb,2*norb)
      complex(rp) :: dx(2*norb,2*norb), dy(2*norb,2*norb), left_phase,right_phase,throwaway(2*norb,2*norb)
      real(rp) :: ez(3),delta(3),d(3), value
      integer :: ineigh, idir
      logical :: supported
      character(len=80) :: reason
      dx=cmplx(0.0_rp,0.0_rp,rp); dy=dx; ez=[0.0_rp,0.0_rp,1.0_rp]
      w0=cmplx(1.0_rp,0.0_rp,rp); w1=cmplx(0.31_rp,0.0_rp,rp); c1=cmplx(0.23_rp,0.0_rp,rp)
      do ineigh=1,3
         value=merge(0.12_rp,0.37_rp,ineigh==1); h=cmplx(0.0_rp,0.0_rp,rp)
         do idir=1,norb; h(idir,idir)=cmplx(value,0.0_rp,rp); end do
         d=0.0_rp; if(ineigh==2)d(1)=1.0_rp; if(ineigh==3)d(1)=-1.0_rp
         call lmto_endpoint_phases(k,q,d,left_phase,right_phase)
         do idir=1,2
            delta=0.0_rp; delta(idir)=1.0_rp
            call lmto_bond_tangent(h,w0,w1,w0,w1,c1,ez,ez,delta,[0.0_rp,0.0_rp,0.0_rp],ineigh==1,tangent)
            cart=tangent; call hcpx(cart(:,:,1),'cart2sph'); call hcpx(cart(:,:,2),'cart2sph'); call hcpx(cart(:,:,3),'cart2sph'); call hcpx(cart(:,:,4),'cart2sph')
            call lmto_hhmag_to_spinor(cart,spinor)
            if(idir==1) dx=dx+spinor*left_phase
            if(idir==2) dy=dy+spinor*left_phase
            call lmto_bond_tangent(h,w0,w1,w0,w1,c1,ez,ez,[0.0_rp,0.0_rp,0.0_rp],delta,ineigh==1,tangent)
            cart=tangent; call hcpx(cart(:,:,1),'cart2sph'); call hcpx(cart(:,:,2),'cart2sph'); call hcpx(cart(:,:,3),'cart2sph'); call hcpx(cart(:,:,4),'cart2sph')
            call lmto_hhmag_to_spinor(cart,spinor)
            if(idir==1) dx=dx+spinor*right_phase
            if(idir==2) dy=dy+spinor*right_phase
         end do
      end do
      call lmto_circular_pair_potential(dx,dy,2.0_rp,qminus,throwaway,supported,reason)
      if(.not.supported) error stop 'expected_reciprocal_fixture_operator failed'
   end subroutine expected_reciprocal_fixture_operator

   ! The finite-q oracle deliberately rebuilds real-space Hamiltonians from
   ! explicit cosine/sine moment textures.  It does not call an endpoint-tangent
   ! record, the finite-q pair-potential assembler, or the phase assembler used
   ! by the analytic path below.
   subroutine test_commensurate_supercell_oracle()
      call test_supercell_q(2)
      call test_supercell_q(3)
   end subroutine test_commensurate_supercell_oracle

   subroutine test_supercell_q(ncell)
      integer, intent(in) :: ncell
      integer :: itheta
      real(rp) :: k, q, theta(3), error(3), relative_error
      complex(rp) :: analytic_x(2,2), analytic_y(2,2), analytic_q(2,2)
      complex(rp) :: cosine_x(2,2), sine_x(2,2), cosine_y(2,2), sine_y(2,2), supercell_q(2,2)

      ! Both primitive endpoints must be commensurate with the explicit
      ! Gamma-supercell; k=0 therefore makes K=k+q an allowed supercell mode.
      k = 0.0_rp
      q = 1.0_rp/real(ncell,rp)
      theta = [2.5e-4_rp, 1.25e-4_rp, 6.25e-5_rp]
      call analytic_primitive_tangent(k,q,1,analytic_x)
      call analytic_primitive_tangent(k,q,2,analytic_y)
      analytic_q = (analytic_x-i_unit*analytic_y)/2.0_rp
      do itheta=1,3
         call supercell_texture_derivative(ncell,k,q,1,.true., theta(itheta),cosine_x)
         call supercell_texture_derivative(ncell,k,q,1,.false.,theta(itheta),sine_x)
         call supercell_texture_derivative(ncell,k,q,2,.true., theta(itheta),cosine_y)
         call supercell_texture_derivative(ncell,k,q,2,.false.,theta(itheta),sine_y)
         ! cos(qR)+i sin(qR)=exp(iqR), so this isolates the +q component
         ! with unit normalization.  No fitted scale factor is applied.
         supercell_q = ((cosine_x+i_unit*sine_x)-i_unit*(cosine_y+i_unit*sine_y))/2.0_rp
         error(itheta) = maxval(abs(supercell_q-analytic_q))
      end do
      relative_error = error(3)/max(maxval(abs(analytic_q)),tiny(1.0_rp))
      write(*,'(a,i0,a,3(1x,es10.3),a,1x,es10.3)') 'Supercell q=1/',ncell,' Q errors:',error,' relative:',relative_error
      call check(error(3), 'commensurate-supercell Q oracle')
      if (error(1) < 3.5_rp*error(2) .or. error(2) < 3.0_rp*error(3)) then
         failed=.true.; write(*,'(a,i0)') 'FAIL supercell central difference is not O(theta^2), q denominator ',ncell
      end if
   end subroutine test_supercell_q

   subroutine analytic_primitive_tangent(k,q,component,derivative)
      real(rp), intent(in) :: k,q
      integer, intent(in) :: component
      complex(rp), intent(out) :: derivative(2,2)
      complex(rp) :: hhop(1,1), hzero(1,1), w0(1), w1(1), c0(1), c1(1), tangent(1,1,4), spinor(2,2)
      real(rp) :: ez(3), direction(3), zero(3), displacement(3)
      integer :: sign

      hhop(1,1)=cmplx(-0.37_rp,0.0_rp,rp); hzero=cmplx(0.0_rp,0.0_rp,rp)
      w0(1)=cmplx(1.0_rp,0.0_rp,rp); w1(1)=cmplx(0.23_rp,0.0_rp,rp)
      c0=cmplx(0.0_rp,0.0_rp,rp); c1(1)=cmplx(0.61_rp,0.0_rp,rp)
      ez=[0.0_rp,0.0_rp,1.0_rp]; zero=0.0_rp; direction=zero; direction(component)=1.0_rp
      ! On-site ownership is one physical moment, so both endpoint variations
      ! enter once while the onsite c1 contribution is retained once.
      call lmto_bond_tangent(hzero,w0,w1,w0,w1,c1,ez,ez,direction,direction,.true.,tangent)
      call lmto_hhmag_to_spinor(tangent,derivative)
      do sign=-1,1,2
         displacement=[real(sign,rp),0.0_rp,0.0_rp]
         call lmto_bond_tangent(hhop,w0,w1,w0,w1,c0,ez,ez,direction,zero,.false.,tangent)
         call lmto_hhmag_to_spinor(tangent,spinor)
         derivative = derivative + spinor*lmto_bloch_phase([k,0.0_rp,0.0_rp],displacement)
         call lmto_bond_tangent(hhop,w0,w1,w0,w1,c0,ez,ez,zero,direction,.false.,tangent)
         call lmto_hhmag_to_spinor(tangent,spinor)
         derivative = derivative + spinor*lmto_bloch_phase([k+q,0.0_rp,0.0_rp],displacement)
      end do
   end subroutine analytic_primitive_tangent

   subroutine supercell_texture_derivative(ncell,k,q,component,is_cosine,theta,derivative)
      integer, intent(in) :: ncell,component
      real(rp), intent(in) :: k,q,theta
      logical, intent(in) :: is_cosine
      complex(rp), intent(out) :: derivative(2,2)
      real(rp) :: moments_plus(3,ncell), moments_minus(3,ncell)
      complex(rp) :: hplus(2*ncell,2*ncell), hminus(2*ncell,2*ncell)

      call make_texture(ncell,q,component,is_cosine, theta,moments_plus)
      call make_texture(ncell,q,component,is_cosine,-theta,moments_minus)
      call build_explicit_supercell(moments_plus,hplus)
      call build_explicit_supercell(moments_minus,hminus)
      call project_supercell_derivative(ncell,k,q,(hplus-hminus)/(2.0_rp*theta),derivative)
   end subroutine supercell_texture_derivative

   subroutine make_texture(ncell,q,component,is_cosine,theta,moments)
      integer, intent(in) :: ncell,component
      real(rp), intent(in) :: q,theta
      logical, intent(in) :: is_cosine
      real(rp), intent(out) :: moments(3,ncell)
      integer :: icell
      real(rp) :: amplitude, phase
      do icell=1,ncell
         phase=2.0_rp*acos(-1.0_rp)*q*real(icell-1,rp)
         amplitude=theta*merge(cos(phase),sin(phase),is_cosine)
         moments(:,icell)=[0.0_rp,0.0_rp,1.0_rp]
         moments(component,icell)=amplitude
         moments(:,icell)=moments(:,icell)/norm2(moments(:,icell))
      end do
   end subroutine make_texture

   subroutine build_explicit_supercell(moments,hamiltonian)
      real(rp), intent(in) :: moments(:,:)
      complex(rp), intent(out) :: hamiltonian(:,:)
      complex(rp) :: hhop(1,1), hzero(1,1), w0(1), w1(1), c0(1), c1(1), hhmag(1,1,4), spinor(2,2)
      integer :: ncell,icell,jcell,idir,ibeg,iend,jbeg,jend

      ncell=size(moments,2); hamiltonian=cmplx(0.0_rp,0.0_rp,rp)
      hhop(1,1)=cmplx(-0.37_rp,0.0_rp,rp); hzero=cmplx(0.0_rp,0.0_rp,rp)
      w0(1)=cmplx(1.0_rp,0.0_rp,rp); w1(1)=cmplx(0.23_rp,0.0_rp,rp)
      c0=cmplx(0.0_rp,0.0_rp,rp); c1(1)=cmplx(0.61_rp,0.0_rp,rp)
      do icell=1,ncell
         ibeg=2*(icell-1)+1; iend=ibeg+1
         call lmto_bond_value(hzero,w0,w1,w0,w1,c0,c1,moments(:,icell),moments(:,icell),.true.,hhmag)
         call lmto_hhmag_to_spinor(hhmag,spinor)
         hamiltonian(ibeg:iend,ibeg:iend)=hamiltonian(ibeg:iend,ibeg:iend)+spinor
         do idir=1,2
            jcell=merge(modulo(icell,ncell)+1,modulo(icell-2,ncell)+1,idir==1)
            jbeg=2*(jcell-1)+1; jend=jbeg+1
            call lmto_bond_value(hhop,w0,w1,w0,w1,c0,c0,moments(:,icell),moments(:,jcell),.false.,hhmag)
            call lmto_hhmag_to_spinor(hhmag,spinor)
            hamiltonian(ibeg:iend,jbeg:jend)=hamiltonian(ibeg:iend,jbeg:jend)+spinor
         end do
      end do
   end subroutine build_explicit_supercell

   subroutine project_supercell_derivative(ncell,k,q,delta_h,derivative)
      integer, intent(in) :: ncell
      real(rp), intent(in) :: k,q
      complex(rp), intent(in) :: delta_h(:,:)
      complex(rp), intent(out) :: derivative(2,2)
      integer :: i,j,ibeg,iend,jbeg,jend
      complex(rp) :: bra,ket
      derivative=cmplx(0.0_rp,0.0_rp,rp)
      do j=1,ncell
         jbeg=2*(j-1)+1; jend=jbeg+1
         ket=exp(cmplx(0.0_rp,2.0_rp*acos(-1.0_rp)*k*real(j-1,rp),rp))
         do i=1,ncell
            ibeg=2*(i-1)+1; iend=ibeg+1
            bra=exp(cmplx(0.0_rp,-2.0_rp*acos(-1.0_rp)*(k+q)*real(i-1,rp),rp))
            derivative=derivative+bra*delta_h(ibeg:iend,jbeg:jend)*ket/real(ncell,rp)
         end do
      end do
   end subroutine project_supercell_derivative

   subroutine sample(h,w0,w1,c0,c1)
      complex(rp),intent(out)::h(2,2),w0(2),w1(2),c0(2),c1(2)
      h=reshape([cmplx(0.8_rp,0.0_rp,rp),cmplx(0.0_rp,0.0_rp,rp),cmplx(0.0_rp,0.0_rp,rp),cmplx(1.2_rp,0.0_rp,rp)],[2,2])
      w0=[cmplx(1.0_rp,0.0_rp,rp),cmplx(0.7_rp,0.0_rp,rp)]; w1=[cmplx(0.2_rp,0.0_rp,rp),cmplx(-0.45_rp,0.0_rp,rp)]
      c0=cmplx(0.0_rp,0.0_rp,rp); c1=[cmplx(0.3_rp,0.0_rp,rp),cmplx(-0.6_rp,0.0_rp,rp)]
   end subroutine sample
end program test_lmto_pair_potential
