program test_lmto_pair_potential
   use precision_mod, only: rp
   use math_mod, only: i_unit
   use lmto_magnetic_tangent_mod, only: lmto_bond_value, lmto_bond_tangent, lmto_hhmag_to_spinor
   use lmto_pair_potential_mod, only: lmto_circular_pair_potential, lmto_bloch_phase
   implicit none
   real(rp), parameter :: tol = 2.0e-9_rp
   real(rp) :: maximum_error
   logical :: failed

   maximum_error = 0.0_rp; failed = .false.
   call test_rotation_oracle_and_adjoint()
   call test_unequal_orbital_negative_control()
   call test_signed_moment_and_bloch_phase()
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
      real(rp)::ez(3),ex(3),ey(3),theta
      logical::supported
      character(len=120)::reason
      call sample(h,w0,w1,c0,c1); ez=[0.0_rp,0.0_rp,1.0_rp]; ex=[1.0_rp,0.0_rp,0.0_rp]; ey=[0.0_rp,1.0_rp,0.0_rp]
      call lmto_bond_tangent(h,w0,w1,w0,w1,c1,ez,ez,ex,[0.0_rp,0.0_rp,0.0_rp],.true.,dxh)
      call lmto_bond_tangent(h,w0,w1,w0,w1,c1,ez,ez,ey,[0.0_rp,0.0_rp,0.0_rp],.true.,dyh)
      call lmto_hhmag_to_spinor(dxh,dx); call lmto_hhmag_to_spinor(dyh,dy)
      call lmto_circular_pair_potential(dx,dy,2.0_rp,qminus,qplus,supported,reason)
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

   subroutine sample(h,w0,w1,c0,c1)
      complex(rp),intent(out)::h(2,2),w0(2),w1(2),c0(2),c1(2)
      h=reshape([cmplx(0.8_rp,0.0_rp,rp),cmplx(0.0_rp,0.0_rp,rp),cmplx(0.0_rp,0.0_rp,rp),cmplx(1.2_rp,0.0_rp,rp)],[2,2])
      w0=[cmplx(1.0_rp,0.0_rp,rp),cmplx(0.7_rp,0.0_rp,rp)]; w1=[cmplx(0.2_rp,0.0_rp,rp),cmplx(-0.45_rp,0.0_rp,rp)]
      c0=cmplx(0.0_rp,0.0_rp,rp); c1=[cmplx(0.3_rp,0.0_rp,rp),cmplx(-0.6_rp,0.0_rp,rp)]
   end subroutine sample
end program test_lmto_pair_potential
