program test_lmto_magnetic_tangents
   use precision_mod, only: rp
   use math_mod, only: i_unit
   use lmto_magnetic_tangent_mod, only: lmto_bond_value, lmto_bond_tangent, lmto_hhmag_to_spinor, &
      lmto_endpoint_tangent_record, lmto_make_endpoint_record, lmto_ordinary_tangent_supported
   implicit none

   real(rp), parameter :: tol = 1.0e-9_rp
   real(rp) :: maximum_error
   logical :: failed

   maximum_error = 0.0_rp
   failed = .false.
   call test_value_preservation_and_spinor_map()
   call test_endpoint_finite_differences()
   call test_dot_and_cross_negative_controls()
   call test_reverse_directed_bond_covariance()
   call test_endpoint_identity_and_capability_guard()
   write (*, '(a,es12.4)') 'LMTO magnetic-tangent maximum error: ', maximum_error
   if (failed .or. maximum_error > tol) error stop 1
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine check(error, label)
      real(rp), intent(in) :: error
      character(len=*), intent(in) :: label
      maximum_error = max(maximum_error, error)
      if (error > tol) then
         write (*, '(a,a,a,es12.4)') 'FAIL ', trim(label), ': ', error
         failed = .true.
      end if
   end subroutine check

   subroutine sample_data(h, w0i, w1i, w0j, w1j, c0, c1, mi, mj)
      complex(rp), intent(out) :: h(2,2), w0i(2), w1i(2), w0j(2), w1j(2), c0(2), c1(2)
      real(rp), intent(out) :: mi(3), mj(3)
      h(1,1) = cmplx(0.7_rp, 0.2_rp, rp); h(1,2) = cmplx(-0.3_rp, 0.4_rp, rp)
      h(2,1) = cmplx(0.1_rp, -0.5_rp, rp); h(2,2) = cmplx(0.8_rp, -0.1_rp, rp)
      w0i = [cmplx(1.1_rp,0.0_rp,rp), cmplx(0.7_rp,0.0_rp,rp)]
      w1i = [cmplx(0.2_rp,0.0_rp,rp), cmplx(-0.4_rp,0.0_rp,rp)]
      w0j = [cmplx(0.9_rp,0.0_rp,rp), cmplx(1.2_rp,0.0_rp,rp)]
      w1j = [cmplx(-0.3_rp,0.0_rp,rp), cmplx(0.25_rp,0.0_rp,rp)]
      c0 = [cmplx(-0.2_rp,0.1_rp,rp), cmplx(0.3_rp,-0.2_rp,rp)]
      c1 = [cmplx(0.45_rp,-0.1_rp,rp), cmplx(-0.15_rp,0.2_rp,rp)]
      mi = normalized([0.3_rp, -0.4_rp, 0.8_rp])
      mj = normalized([-0.5_rp, 0.6_rp, 0.4_rp])
   end subroutine sample_data

   subroutine test_value_preservation_and_spinor_map()
      complex(rp) :: h(2,2), w0i(2), w1i(2), w0j(2), w1j(2), c0(2), c1(2), value(2,2,4), reference(2,2,4)
      complex(rp) :: spinor(4,4), spinor_reference(4,4)
      real(rp) :: mi(3), mj(3)
      call sample_data(h,w0i,w1i,w0j,w1j,c0,c1,mi,mj)
      call lmto_bond_value(h,w0i,w1i,w0j,w1j,c0,c1,mi,mj,.true.,value)
      call old_value(h,w0i,w1i,w0j,w1j,c0,c1,mi,mj,.true.,reference)
      call check(maxval(abs(value-reference)), 'shared core preserves old onsite value')
      call lmto_bond_value(h,w0i,w1i,w0j,w1j,c0,c1,mi,mj,.false.,value)
      call old_value(h,w0i,w1i,w0j,w1j,c0,c1,mi,mj,.false.,reference)
      call check(maxval(abs(value-reference)), 'shared core preserves finite-bond value')
      call lmto_hhmag_to_spinor(value, spinor)
      call old_spinor(value, spinor_reference)
      call check(maxval(abs(spinor-spinor_reference)), 'complex H0+Hvec sigma spinor map')
   end subroutine test_value_preservation_and_spinor_map

   subroutine test_endpoint_finite_differences()
      complex(rp) :: h(2,2), w0i(2), w1i(2), w0j(2), w1j(2), c0(2), c1(2), tangent(2,2,4), vp(2,2,4), vm(2,2,4)
      real(rp) :: mi(3), mj(3), di(3), dj(3), theta(3), err(3)
      integer :: n
      call sample_data(h,w0i,w1i,w0j,w1j,c0,c1,mi,mj)
      di = cross3([0.2_rp,-0.7_rp,0.6_rp],mi)
      dj = cross3([-0.5_rp,0.1_rp,0.8_rp],mj)
      theta = [1.0e-3_rp, 5.0e-4_rp, 2.5e-4_rp]
      call lmto_bond_tangent(h,w0i,w1i,w0j,w1j,c1,mi,mj,di,[0.0_rp,0.0_rp,0.0_rp],.false.,tangent)
      do n=1,3
         call lmto_bond_value(h,w0i,w1i,w0j,w1j,c0,c1,rotate(mi,di,theta(n)),mj,.false.,vp)
         call lmto_bond_value(h,w0i,w1i,w0j,w1j,c0,c1,rotate(mi,di,-theta(n)),mj,.false.,vm)
         err(n)=maxval(abs((vp-vm)/(2.0_rp*theta(n))-tangent))
      end do
      call fd_report('left endpoint',theta,err)
      call check(err(3), 'left endpoint central tangent')
      if (err(1) < 3.5_rp*err(2) .or. err(2) < 3.5_rp*err(3)) then; failed=.true.; write(*,'(a)') 'FAIL left FD is not O(theta^2)'; end if
      call lmto_bond_tangent(h,w0i,w1i,w0j,w1j,c1,mi,mj,[0.0_rp,0.0_rp,0.0_rp],dj,.false.,tangent)
      do n=1,3
         call lmto_bond_value(h,w0i,w1i,w0j,w1j,c0,c1,mi,rotate(mj,dj,theta(n)),.false.,vp)
         call lmto_bond_value(h,w0i,w1i,w0j,w1j,c0,c1,mi,rotate(mj,dj,-theta(n)),.false.,vm)
         err(n)=maxval(abs((vp-vm)/(2.0_rp*theta(n))-tangent))
      end do
      call fd_report('right endpoint',theta,err); call check(err(3), 'right endpoint central tangent')
      ! A common rigid rotation has both endpoint contributions; neither is lost.
      di = cross3([0.4_rp,0.3_rp,-0.2_rp],mi); dj = cross3([0.4_rp,0.3_rp,-0.2_rp],mj)
      call lmto_bond_tangent(h,w0i,w1i,w0j,w1j,c1,mi,mj,di,dj,.false.,tangent)
      call lmto_bond_value(h,w0i,w1i,w0j,w1j,c0,c1,rotate(mi,di,theta(3)),rotate(mj,dj,theta(3)),.false.,vp)
      call lmto_bond_value(h,w0i,w1i,w0j,w1j,c0,c1,rotate(mi,di,-theta(3)),rotate(mj,dj,-theta(3)),.false.,vm)
      call check(maxval(abs((vp-vm)/(2.0_rp*theta(3))-tangent)), 'common rigid rotation tangent')
   end subroutine test_endpoint_finite_differences

   subroutine test_dot_and_cross_negative_controls()
      complex(rp) :: h(2,2), w0i(2), w1i(2), w0j(2), w1j(2), c0(2), c1(2), tangent(2,2,4), vp(2,2,4), vm(2,2,4)
      real(rp) :: mi(3), mj(3), di(3), dcross(3), exact_error, bad_error
      integer :: k
      call sample_data(h,w0i,w1i,w0j,w1j,c0,c1,mi,mj)
      di = cross3([0.6_rp,-0.1_rp,0.5_rp],mi)
      call lmto_bond_tangent(h,w0i,w1i,w0j,w1j,c1,mi,mj,di,[0.0_rp,0.0_rp,0.0_rp],.false.,tangent)
      call lmto_bond_value(h,w0i,w1i,w0j,w1j,c0,c1,rotate(mi,di,1.0e-5_rp),mj,.false.,vp)
      call lmto_bond_value(h,w0i,w1i,w0j,w1j,c0,c1,rotate(mi,di,-1.0e-5_rp),mj,.false.,vm)
      exact_error=maxval(abs((vp-vm)/2.0e-5_rp-tangent))
      dcross = cross3(di,mj)
      do k=1,3
         tangent(:,:,k)=tangent(:,:,k)-i_unit*outer(w1i,w1j,h)*cmplx(dcross(k),0.0_rp,rp)
      end do
      bad_error=maxval(abs((vp-vm)/2.0e-5_rp-tangent))
      if (bad_error < 1.0e3_rp*max(exact_error,1.0e-14_rp)) then; failed=.true.; write(*,'(a)') 'FAIL cross-term negative control'; end if
      call lmto_bond_tangent(h,w0i,w1i,w0j,w1j,c1,mi,mj,di,[0.0_rp,0.0_rp,0.0_rp],.false.,tangent)
      tangent(:,:,4)=tangent(:,:,4)-outer(w1i,w1j,h)*cmplx(dot_product(di,mj),0.0_rp,rp)
      bad_error=maxval(abs((vp-vm)/2.0e-5_rp-tangent))
      if (bad_error < 1.0e3_rp*max(exact_error,1.0e-14_rp)) then; failed=.true.; write(*,'(a)') 'FAIL dot-term negative control'; end if
   end subroutine test_dot_and_cross_negative_controls

   subroutine test_reverse_directed_bond_covariance()
      complex(rp) :: h(2,2), w0i(2), w1i(2), w0j(2), w1j(2), c0(2), c1(2), forward(2,2,4), reverse(2,2,4), df(2,2,4), dr(2,2,4)
      real(rp) :: mi(3), mj(3), di(3)
      call sample_data(h,w0i,w1i,w0j,w1j,c0,c1,mi,mj); di=cross3([0.1_rp,0.9_rp,-0.2_rp],mi)
      call lmto_bond_value(h,w0i,w1i,w0j,w1j,c0,c1,mi,mj,.false.,forward)
      call lmto_bond_value(transpose(conjg(h)),w0j,w1j,w0i,w1i,c0,c1,mj,mi,.false.,reverse)
      call check(maxval(abs(reverse-adjoint4(forward))), 'reverse directed bond Hermiticity')
      call lmto_bond_tangent(h,w0i,w1i,w0j,w1j,c1,mi,mj,di,[0.0_rp,0.0_rp,0.0_rp],.false.,df)
      call lmto_bond_tangent(transpose(conjg(h)),w0j,w1j,w0i,w1i,c1,mj,mi,[0.0_rp,0.0_rp,0.0_rp],di,.false.,dr)
      call check(maxval(abs(dr-adjoint4(df))), 'reverse endpoint tangent covariance')
   end subroutine test_reverse_directed_bond_covariance

   subroutine test_endpoint_identity_and_capability_guard()
      type(lmto_endpoint_tangent_record) :: record
      record = lmto_make_endpoint_record(4, 9, 1, 1, 7, 23, [1.2_rp, -0.5_rp, 0.0_rp], .true.)
      if (record%source_site == record%neighbor_site .or. record%source_type /= record%neighbor_type .or. &
          record%directed_bond /= 7 .or. .not. record%supported) then
         failed = .true.; write(*,'(a)') 'FAIL same-type endpoint identities aliased'
      end if
      if (lmto_ordinary_tangent_supported(.true.,.false.,.false.,.false.,.false.,.false.,.false.) .or. &
          lmto_ordinary_tangent_supported(.false.,.true.,.false.,.false.,.false.,.false.,.false.) .or. &
          lmto_ordinary_tangent_supported(.false.,.false.,.true.,.false.,.false.,.false.,.false.) .or. &
          lmto_ordinary_tangent_supported(.false.,.false.,.false.,.true.,.false.,.false.,.false.) .or. &
          lmto_ordinary_tangent_supported(.false.,.false.,.false.,.false.,.true.,.false.,.false.) .or. &
          lmto_ordinary_tangent_supported(.false.,.false.,.false.,.false.,.false.,.true.,.false.) .or. &
          lmto_ordinary_tangent_supported(.false.,.false.,.false.,.false.,.false.,.false.,.true.)) then
         failed = .true.; write(*,'(a)') 'FAIL unsupported tangent capability was accepted'
      end if
   end subroutine test_endpoint_identity_and_capability_guard

   subroutine old_value(h,w0i,w1i,w0j,w1j,c0,c1,mi,mj,onsite,value)
      complex(rp),intent(in)::h(2,2),w0i(2),w1i(2),w0j(2),w1j(2),c0(2),c1(2); real(rp),intent(in)::mi(3),mj(3); logical,intent(in)::onsite; complex(rp),intent(out)::value(2,2,4)
      integer::i,j,k; complex(rp)::cross(3); value=cmplx(0.0_rp,0.0_rp,rp)
      cross=cmplx(cross3(mi,mj),0.0_rp,rp)
      do j=1,2; do i=1,2
         value(i,j,4)=w0i(i)*h(i,j)*w0j(j)+w1i(i)*h(i,j)*w1j(j)*dot_product(mi,mj)
         do k=1,3; value(i,j,k)=w1i(i)*h(i,j)*w0j(j)*mi(k)+w0i(i)*h(i,j)*w1j(j)*mj(k)+i_unit*w1i(i)*h(i,j)*w1j(j)*cross(k); end do
      end do; end do
      if(onsite) then; do i=1,2; value(i,i,4)=value(i,i,4)+c0(i); do k=1,3; value(i,i,k)=value(i,i,k)+c1(i)*mi(k); end do; end do; end if
   end subroutine old_value

   subroutine old_spinor(h, s)
      complex(rp),intent(in)::h(2,2,4); complex(rp),intent(out)::s(4,4); integer::i,j
      s=cmplx(0.0_rp,0.0_rp,rp); do j=1,2; do i=1,2; s(i,j)=h(i,j,4)+h(i,j,3); s(i+2,j+2)=h(i,j,4)-h(i,j,3); s(i,j+2)=h(i,j,1)-i_unit*h(i,j,2); s(i+2,j)=h(i,j,1)+i_unit*h(i,j,2); end do; end do
   end subroutine old_spinor

   pure function outer(a,b,h) result(x)
      complex(rp),intent(in)::a(2),b(2),h(2,2); complex(rp)::x(2,2); integer::i,j
      do j=1,2; do i=1,2; x(i,j)=a(i)*h(i,j)*b(j); end do; end do
   end function outer
   pure function adjoint4(h) result(a)
      complex(rp),intent(in)::h(2,2,4); complex(rp)::a(2,2,4); integer::k
      do k=1,4; a(:,:,k)=transpose(conjg(h(:,:,k))); end do
   end function adjoint4
   pure function normalized(x) result(y); real(rp),intent(in)::x(3); real(rp)::y(3); y=x/norm2(x); end function normalized
   pure function rotate(m,dm,theta) result(r); real(rp),intent(in)::m(3),dm(3),theta; real(rp)::r(3),axis(3),speed; speed=norm2(dm); axis=cross3(m,dm)/(speed+tiny(1.0_rp)); r=cos(theta*speed)*m+sin(theta*speed)*cross3(axis,m); end function rotate
   pure function cross3(a,b) result(c); real(rp),intent(in)::a(3),b(3); real(rp)::c(3); c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]; end function cross3
   subroutine fd_report(label,theta,error); character(len=*),intent(in)::label; real(rp),intent(in)::theta(3),error(3); write(*,'(a,3(1x,es10.3))') trim(label)//' FD errors:',error; end subroutine fd_report
end program test_lmto_magnetic_tangents
