! Independent finite-q oracle for the second-order local-rotation adapter.
program test_dresp03q_second_order_rotation
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use precision_mod, only: rp
   use math_mod, only: pi
   use lr_kl_hessian_mod, only: lmto_live_hamiltonian_fixture, lmto_fixture_init, &
      assemble_lmto_hamiltonian, assemble_lmto_finite_q_torque, assemble_lmto_finite_q_mixed_derivative
   implicit none
   integer, parameter :: nc=4, ns=2
   real(rp), parameter :: q0(3)=[0.0_rp,0.0_rp,0.0_rp], qf(3)=[0.25_rp,0.0_rp,0.0_rp]
   real(rp), parameter :: k(3)=[0.0_rp,0.0_rp,0.0_rp], axis(3)=[1.0_rp,0.0_rp,0.0_rp]
   real(rp), parameter :: steps(3)=[0.04_rp,0.01_rp,0.005_rp]
   type(lmto_live_hamiltonian_fixture) :: fix
   complex(rp) :: h(2,2), h_direct(2,2), torque(2,2), contact(2,2), fd(2,2)
   real(rp) :: base_error, torque_err(2,3), contact_err(2,3), conv_t, conv_c
   integer :: iq, is, fail

   call make_fixture(fix)
   call assemble_lmto_hamiltonian(fix,k,h)
   call direct_projection(fix,k,k,cmplx(0.0_rp,0.0_rp,rp),cmplx(0.0_rp,0.0_rp,rp),h_direct)
   base_error=relative_matrix(h_direct,h)
   write(*,'(a,es14.6)') 'second-order H fixture/direct baseline residual = ',base_error
   if (base_error>5.0e-13_rp) error stop 'second-order independent H construction does not match the fixture'

   fail=0
   do iq=1,2
      if (iq==1) then
         call assemble_lmto_finite_q_torque(fix,k,q0,1,axis,torque)
      else
         call assemble_lmto_finite_q_torque(fix,k,qf,1,axis,torque)
      end if
      do is=1,size(steps)
         fd=torque_fd(fix,k,merge(q0,qf,iq==1),steps(is))
         torque_err(iq,is)=relative_matrix(fd,torque)
      end do
      write(*,'(a,3(es14.6,1x))') 'second-order torque FD errors q0/finite-q steps = ',torque_err(iq,:)
      if (torque_err(iq,3)>2.0e-5_rp .or. torque_err(iq,2)>0.45_rp*torque_err(iq,1)) fail=fail+1
   end do

   do iq=1,2
      if (iq==1) then
         call assemble_lmto_finite_q_mixed_derivative(fix,k,q0,1,axis,1,axis,contact)
      else
         call assemble_lmto_finite_q_mixed_derivative(fix,k,qf,1,axis,1,axis,contact)
      end if
      do is=1,size(steps)
         fd=contact_fd(fix,k,merge(q0,qf,iq==1),steps(is))
         contact_err(iq,is)=relative_matrix(fd,contact)
      end do
      write(*,'(a,3(es14.6,1x))') 'second-order contact four-point errors q0/finite-q steps = ',contact_err(iq,:)
      if (contact_err(iq,3)>3.0e-5_rp .or. contact_err(iq,2)>0.45_rp*contact_err(iq,1)) fail=fail+1
   end do
   conv_t=max(torque_err(1,3)/max(torque_err(1,1),tiny(1.0_rp)),torque_err(2,3)/max(torque_err(2,1),tiny(1.0_rp)))
   conv_c=max(contact_err(1,3)/max(contact_err(1,1),tiny(1.0_rp)),contact_err(2,3)/max(contact_err(2,1),tiny(1.0_rp)))
   write(*,'(a,es14.6)') 'second-order torque step error ratio small/large = ',conv_t
   write(*,'(a,es14.6)') 'second-order contact step error ratio small/large = ',conv_c
   if (.not.ieee_is_finite(conv_t) .or. .not.ieee_is_finite(conv_c)) fail=fail+1
   call fix%clear()
   if (fail/=0) error stop 'second-order local-rotation adapter oracle failed'
   write(*,'(a)') 'second-order local-rotation adapter oracle: PASS'

contains

   subroutine make_fixture(f)
      type(lmto_live_hamiltonian_fixture), intent(out) :: f
      call lmto_fixture_init(f,1,1,3,.true.)
      f%include_enu=.true.
      f%moments(:,1)=[0.0_rp,0.0_rp,1.0_rp]
      f%wx0(:,1)=cmplx(1.0_rp,0.0_rp,rp); f%wx1(:,1)=cmplx(0.13_rp,0.0_rp,rp)
      f%c0(:,1)=cmplx(-0.20_rp,0.0_rp,rp); f%c1(:,1)=cmplx(-0.47_rp,0.0_rp,rp)
      f%obar0(:,1)=cmplx(0.035_rp,0.0_rp,rp); f%obar1(:,1)=cmplx(0.082_rp,0.0_rp,rp)
      f%enu0(:,1)=cmplx(-0.11_rp,0.0_rp,rp); f%enu1(:,1)=cmplx(-0.027_rp,0.0_rp,rp)
      f%bond_source=[1,1,1]; f%bond_target=[1,1,1]; f%onsite=[.true.,.false.,.false.]
      f%bond_vector=0.0_rp; f%bond_vector(1,2)=1.0_rp; f%bond_vector(1,3)=-1.0_rp
      f%hhh=cmplx(0.0_rp,0.0_rp,rp); f%hhh(1,1,2)=cmplx(0.12_rp,0.0_rp,rp); f%hhh(1,1,3)=cmplx(0.12_rp,0.0_rp,rp)
   end subroutine make_fixture

   function torque_fd(f,kpoint,qpoint,step) result(value)
      type(lmto_live_hamiltonian_fixture), intent(in) :: f
      real(rp), intent(in) :: kpoint(3),qpoint(3),step
      complex(rp) :: value(2,2), hp(2,2),hm(2,2)
      call direct_projection(f,kpoint+qpoint,kpoint,cmplx(step,0.0_rp,rp),cmplx(0.0_rp,0.0_rp,rp),hp,qpoint,axis)
      call direct_projection(f,kpoint+qpoint,kpoint,cmplx(-step,0.0_rp,rp),cmplx(0.0_rp,0.0_rp,rp),hm,qpoint,axis)
      value=(hp-hm)/(2.0_rp*step)
   end function torque_fd

   function contact_fd(f,kpoint,qpoint,step) result(value)
      type(lmto_live_hamiltonian_fixture), intent(in) :: f
      real(rp), intent(in) :: kpoint(3),qpoint(3),step
      complex(rp) :: value(2,2),hpp(2,2),hpm(2,2),hmp(2,2),hmm(2,2)
      call direct_projection(f,kpoint,kpoint,cmplx(step,0.0_rp,rp),cmplx(step,0.0_rp,rp),hpp,qpoint,axis)
      call direct_projection(f,kpoint,kpoint,cmplx(step,0.0_rp,rp),cmplx(-step,0.0_rp,rp),hpm,qpoint,axis)
      call direct_projection(f,kpoint,kpoint,cmplx(-step,0.0_rp,rp),cmplx(step,0.0_rp,rp),hmp,qpoint,axis)
      call direct_projection(f,kpoint,kpoint,cmplx(-step,0.0_rp,rp),cmplx(-step,0.0_rp,rp),hmm,qpoint,axis)
      value=(hpp-hpm-hmp+hmm)/(4.0_rp*step*step)
   end function contact_fd

   subroutine direct_projection(f,krow,kcol,xamp,yamp,result,qpoint,rot_axis)
      type(lmto_live_hamiltonian_fixture), intent(in) :: f
      real(rp), intent(in) :: krow(3),kcol(3)
      complex(rp), intent(in) :: xamp,yamp
      complex(rp), intent(out) :: result(2,2)
      real(rp), intent(in), optional :: qpoint(3),rot_axis(3)
      complex(rp) :: hs(2*nc,2*nc)
      if (present(qpoint) .and. present(rot_axis)) then
         call direct_h2_supercell(f,xamp,yamp,qpoint,rot_axis,hs)
      else
         call direct_h2_supercell(f,cmplx(0.0_rp,0.0_rp,rp),cmplx(0.0_rp,0.0_rp,rp),q0,axis,hs)
      end if
      call fourier_block(hs,krow,kcol,result)
   end subroutine direct_projection

   subroutine direct_h2_supercell(f,xamp,yamp,qpoint,rot_axis,h)
      type(lmto_live_hamiltonian_fixture), intent(in) :: f
      complex(rp), intent(in) :: xamp,yamp
      real(rp), intent(in) :: qpoint(3),rot_axis(3)
      complex(rp), intent(out) :: h(2*nc,2*nc)
      complex(rp) :: b(2*nc,2*nc),qmat(2*nc,2*nc),enu(2*nc,2*nc),bond(2,2),oblock(2,2),eblock(2,2)
      complex(rp) :: msrc(3),mtgt(3),hmag(1,1,4),phase
      integer :: cell,ibond,target_cell,disp,first,last
      b=0.0_rp; qmat=0.0_rp; enu=0.0_rp
      do cell=0,nc-1
         do ibond=1,f%nbond
            disp=nint(f%bond_vector(1,ibond)); target_cell=modulo(cell+disp,nc)
            msrc=rotated_moment(f%moments(:,1),rot_axis,cell_angle(cell,xamp,yamp,qpoint))
            mtgt=rotated_moment(f%moments(:,1),rot_axis,cell_angle(target_cell,xamp,yamp,qpoint))
            call direct_bond_value(f,ibond,msrc,mtgt,hmag); call spinor(hmag,bond)
            first=2*cell+1; last=2*target_cell+1
            b(first:first+1,last:last+1)=b(first:first+1,last:last+1)+bond
            call direct_coefficient(f%obar0(:,1),f%obar1(:,1),mtgt,hmag); call spinor(hmag,oblock)
            qmat(first:first+1,last:last+1)=qmat(first:first+1,last:last+1)+matmul(bond,oblock)
         end do
         msrc=rotated_moment(f%moments(:,1),rot_axis,cell_angle(cell,xamp,yamp,qpoint))
         call direct_coefficient(f%enu0(:,1),f%enu1(:,1),msrc,hmag); call spinor(hmag,eblock)
         first=2*cell+1; enu(first:first+1,first:first+1)=eblock
      end do
      h=b+enu
      if (f%hoh) h=h-matmul(qmat,b)
   end subroutine direct_h2_supercell

   function cell_angle(cell,xamp,yamp,qpoint) result(angle)
      integer, intent(in) :: cell
      complex(rp), intent(in) :: xamp,yamp
      real(rp), intent(in) :: qpoint(3)
      complex(rp) :: angle,phase
      phase=exp(cmplx(0.0_rp,2.0_rp*pi*qpoint(1)*real(cell,rp),rp))
      angle=xamp*phase+yamp*conjg(phase)
   end function cell_angle

   function rotated_moment(moment,rot_axis,angle) result(rotated)
      real(rp), intent(in) :: moment(3),rot_axis(3)
      complex(rp), intent(in) :: angle
      complex(rp) :: rotated(3),m(3),a(3),cross(3),dot
      m=cmplx(moment,0.0_rp,rp); a=cmplx(rot_axis,0.0_rp,rp)
      cross=[a(2)*m(3)-a(3)*m(2),a(3)*m(1)-a(1)*m(3),a(1)*m(2)-a(2)*m(1)]
      dot=sum(a*m)
      rotated=m*cos(angle)+cross*sin(angle)+a*dot*(1.0_rp-cos(angle))
   end function rotated_moment

   subroutine direct_bond_value(f,ibond,mi,mj,hmag)
      type(lmto_live_hamiltonian_fixture), intent(in) :: f
      integer, intent(in) :: ibond
      complex(rp), intent(in) :: mi(3),mj(3)
      complex(rp), intent(out) :: hmag(1,1,4)
      complex(rp) :: dot,cross(3)
      integer :: i
      dot=sum(mi*mj)
      cross=[mi(2)*mj(3)-mi(3)*mj(2),mi(3)*mj(1)-mi(1)*mj(3),mi(1)*mj(2)-mi(2)*mj(1)]
      hmag=0.0_rp
      hmag(1,1,4)=f%wx0(1,1)*f%hhh(1,1,ibond)*f%wx0(1,1)+ &
         f%wx1(1,1)*f%hhh(1,1,ibond)*f%wx1(1,1)*dot
      do i=1,3
         hmag(1,1,i)=f%wx1(1,1)*f%hhh(1,1,ibond)*f%wx0(1,1)*mi(i)+ &
            f%wx0(1,1)*f%hhh(1,1,ibond)*f%wx1(1,1)*mj(i)+ &
            cmplx(0.0_rp,1.0_rp,rp)*f%wx1(1,1)*f%hhh(1,1,ibond)*f%wx1(1,1)*cross(i)
      end do
      if (f%onsite(ibond)) then
         hmag(1,1,4)=hmag(1,1,4)+f%c0(1,1)
         hmag(1,1,1:3)=hmag(1,1,1:3)+f%c1(1,1)*mi
      end if
   end subroutine direct_bond_value

   subroutine direct_coefficient(c0,c1,m,hmag)
      complex(rp), intent(in) :: c0(:),c1(:),m(3)
      complex(rp), intent(out) :: hmag(1,1,4)
      hmag=0.0_rp; hmag(1,1,4)=c0(1); hmag(1,1,1:3)=c1(1)*m
   end subroutine direct_coefficient

   subroutine spinor(hmag,block)
      complex(rp), intent(in) :: hmag(1,1,4)
      complex(rp), intent(out) :: block(2,2)
      block(1,1)=hmag(1,1,4)+hmag(1,1,3)
      block(2,2)=hmag(1,1,4)-hmag(1,1,3)
      block(1,2)=hmag(1,1,1)-cmplx(0.0_rp,1.0_rp,rp)*hmag(1,1,2)
      block(2,1)=hmag(1,1,1)+cmplx(0.0_rp,1.0_rp,rp)*hmag(1,1,2)
   end subroutine spinor

   subroutine fourier_block(h,krow,kcol,block)
      complex(rp), intent(in) :: h(2*nc,2*nc)
      real(rp), intent(in) :: krow(3),kcol(3)
      complex(rp), intent(out) :: block(2,2)
      complex(rp) :: phase
      integer :: r,s
      block=0.0_rp
      do r=0,nc-1
         do s=0,nc-1
            phase=exp(cmplx(0.0_rp,2.0_rp*pi*(kcol(1)*real(s,rp)-krow(1)*real(r,rp)),rp))
            block=block+phase*h(2*r+1:2*r+2,2*s+1:2*s+2)/real(nc,rp)
         end do
      end do
   end subroutine fourier_block

   function relative_matrix(a,b) result(value)
      complex(rp), intent(in) :: a(:,:),b(:,:)
      real(rp) :: value
      value=sqrt(sum(abs(a-b)**2))/max(sqrt(sum(abs(b)**2)),tiny(1.0_rp))
   end function relative_matrix

end program test_dresp03q_second_order_rotation
