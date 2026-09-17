!------------------------------------------------------------------------------
! DRESP-03G -- independent finite-H spectral/contour oracle.
!------------------------------------------------------------------------------
program test_dresp03g_contour
   use precision_mod, only: rp
   use lr_kl_hessian_mod, only: force_theorem_hessian_from_eigenbasis, &
      force_theorem_finite_q_hessian_from_eigenbasis_metallic, force_theorem_finite_q_hessian_from_eigenbasis
   use lr_kl_contour_mod, only: finite_h_contour_options, finite_h_contour_report, &
      force_theorem_finite_q_hessian_from_resolvent, force_theorem_finite_q_hessian_from_zero_temperature_contour
   implicit none

   complex(rp) :: h(4,4), ti(4,4), tj(4,4), cij(4,4), tq(4,4,1), tm(4,4,1), cstack(4,4,1,1), vectors(4,4)
   complex(rp) :: hs(4,4), he(4,4), tqs(4,4,1), tms(4,4,1), cs(4,4,1,1)
   complex(rp) :: gf_h(1,1), gf_tt(1,1), gf_cc(1,1), gf_all(1,1), gf_reverse(1,1)
   complex(rp) :: sp_h(1,1), sp_tt(1,1), sp_cc(1,1), sp_all(1,1)
   real(rp) :: values(4), fermi, kT, fd, residual, spectral_residual, fd_residual
   real(rp) :: zero_residual, degenerate_residual, q_covariance_residual, convergence(4), old_tt, old_cc, old_all
   integer :: contour_sequence(4)
   type(finite_h_contour_options) :: options
   type(finite_h_contour_report) :: finite_report, zero_report
   integer :: i
   contour_sequence = [8,16,32,64]

   call make_noncommuting_fixture(h,ti,tj,cij)
   call diagonalize(h,values,vectors)
   fermi = 0.12_rp; kT = 0.08_rp
   tq(:,:,1) = ti; tm(:,:,1) = tj; cstack(:,:,1,1) = cij
   options%contour_margin = 0.30_rp
   options%contour_height_fraction = 0.35_rp
   options%account_fermi_poles = .true.
   options%contour_points = 128
   options%account_fermi_poles = .true.
   call force_theorem_finite_q_hessian_from_eigenbasis_metallic(values, vectors, values, vectors, fermi, kT, &
      tq,tm,cstack,sp_h,sp_tt,sp_cc,sp_all)
   call force_theorem_finite_q_hessian_from_resolvent(h,h,fermi,kT,tq,tm,cstack,options,gf_h,gf_tt,gf_cc,gf_all,finite_report)
   spectral_residual = maxval(abs([real(gf_tt(1,1)-sp_tt(1,1),rp), real(gf_cc(1,1)-sp_cc(1,1),rp), &
      real(gf_all(1,1)-sp_all(1,1),rp)]))
   fd = finite_difference_hessian(h,ti,tj,cij,fermi,kT,2.0e-4_rp)
   fd_residual = abs(fd-real(gf_all(1,1),rp))
   write(*,'(a,es14.6)') 'DRESP-03G noncommuting spectral/contour residual = ', spectral_residual
   write(*,'(a,3(es14.6,1x))') 'DRESP-03G spectral TT/contact/total = ', real(sp_tt(1,1),rp), real(sp_cc(1,1),rp), real(sp_all(1,1),rp)
   write(*,'(a,es14.6)') 'DRESP-03G direct finite-difference/contour residual = ', fd_residual
   write(*,'(a,i0,a,i0)') 'DRESP-03G first contour nodes/poles = ', finite_report%contour_points, '/', finite_report%fermi_poles
   write(*,'(a,3(es14.6,1x))') 'DRESP-03G finite-T contour TT/contact/total = ', real(gf_tt(1,1),rp), &
      real(gf_cc(1,1),rp), real(gf_all(1,1),rp)
   if (spectral_residual > 2.0e-9_rp .or. fd_residual > 3.0e-6_rp) error stop 'DRESP-03G finite-T oracle failed'

   call make_gapped_fixture(h,ti,tj,cij)
   call diagonalize(h,values,vectors)
   tq(:,:,1) = ti; tm(:,:,1) = tj; cstack(:,:,1,1) = cij
   options%contour_points = 128
   call force_theorem_hessian_from_eigenbasis(values,vectors,0.0_rp,ti,tj,cij,old_tt,old_cc,old_all)
   call force_theorem_finite_q_hessian_from_zero_temperature_contour(h,h,[-1.2_rp,-0.25_rp],tq,tm,cstack,options, &
      gf_h,gf_tt,gf_cc,gf_all,zero_report)
   zero_residual = maxval(abs([real(gf_tt(1,1),rp)-old_tt, real(gf_cc(1,1),rp)-old_cc, &
      real(gf_all(1,1),rp)-old_all]))
   write(*,'(a,es14.6)') 'DRESP-03G zero-T occupied-contour residual = ', zero_residual
   if (zero_residual > 2.0e-9_rp) error stop 'DRESP-03G zero-T oracle failed'

   call degeneracy_similarity_test(degenerate_residual, h,ti,tj,cij)
   write(*,'(a,es14.6)') 'DRESP-03G degenerate similarity residual = ', degenerate_residual
   if (degenerate_residual > 2.0e-9_rp) error stop 'DRESP-03G degeneracy invariance failed'

   call make_noncommuting_fixture(h,ti,tj,cij)
   tq(:,:,1) = ti; tm(:,:,1) = tj; cstack(:,:,1,1) = cij
   do i = 1, 4
      options%contour_points = contour_sequence(i)
      call force_theorem_finite_q_hessian_from_resolvent(h,h,fermi,kT,tq,tm,cstack,options,gf_h,gf_tt,gf_cc,gf_all)
      convergence(i) = abs(real(gf_all(1,1),rp)-real(sp_all(1,1),rp))
      write(*,'(a,i0,a,3(es14.6,1x))') 'DRESP-03G Ncontour=',options%contour_points,' TT/contact/total residual = ', &
         abs(real(gf_tt(1,1)-sp_tt(1,1),rp)), abs(real(gf_cc(1,1)-sp_cc(1,1),rp)), convergence(i)
   end do
   if (convergence(4) > 2.0e-9_rp .or. convergence(4) > convergence(1)) error stop 'DRESP-03G contour convergence failed'

   ! A genuinely finite-q ordering check: H(k) and H(k+q) differ, and the
   ! two directional vertices are not Hermitian copies of one another.
   hs = h; he = h; he(1,1) = he(1,1) + 0.17_rp; he(2,3) = he(2,3) + cmplx(0.03_rp,-0.04_rp,rp); he(3,2) = conjg(he(2,3))
   tqs(:,:,1) = ti + cmplx(0.01_rp,-0.02_rp,rp); tms(:,:,1) = transpose(conjg(tqs(:,:,1)))
   ! The contact block is k-dependent in the production object.  Use the
   ! zero-contact special case here so the covariance check isolates the
   ! ordered finite-q torque product while the noncommuting contact block is
   ! covered by the preceding oracle and finite-difference gates.
   cs(:,:,1,1) = 0.0_rp
   call force_theorem_finite_q_hessian_from_resolvent(hs,he,fermi,kT,tqs,tms,cs,options,gf_h,gf_tt,gf_cc,gf_all)
   call force_theorem_finite_q_hessian_from_resolvent(he,hs,fermi,kT,tms,tqs,cs,options,gf_h,gf_tt,gf_cc,gf_reverse)
   q_covariance_residual = abs(gf_all(1,1)-gf_reverse(1,1))
   write(*,'(a,es14.6)') 'DRESP-03G q/-q covariance residual = ', q_covariance_residual
   if (q_covariance_residual > 2.0e-9_rp) error stop 'DRESP-03G finite-q covariance test failed'

   write(*,'(a,i0,a,i0)') 'DRESP-03G finite-T contour nodes/pole residues = ', finite_report%contour_points, '/', finite_report%fermi_poles
   write(*,'(a,i0)') 'DRESP-03G zero-T contour nodes = ', zero_report%contour_points
   write(*,'(a)') 'DRESP-03G: PASS'

contains

   subroutine make_noncommuting_fixture(h,ti,tj,c)
      complex(rp), intent(out) :: h(4,4),ti(4,4),tj(4,4),c(4,4)
      h = 0.0_rp; h(1,1)=-0.91_rp; h(2,2)=-0.19_rp; h(3,3)=0.37_rp; h(4,4)=0.83_rp
      h(1,2)=cmplx(0.12_rp,0.07_rp,rp); h(2,1)=conjg(h(1,2)); h(2,3)=cmplx(-0.08_rp,0.11_rp,rp); h(3,2)=conjg(h(2,3))
      h(3,4)=cmplx(0.09_rp,-0.06_rp,rp); h(4,3)=conjg(h(3,4)); h(1,4)=cmplx(-0.04_rp,0.03_rp,rp); h(4,1)=conjg(h(1,4))
      ti = 0.0_rp; ti(1,1)=0.21_rp; ti(2,2)=-0.13_rp; ti(3,3)=0.17_rp; ti(4,4)=-0.09_rp
      ti(1,2)=cmplx(0.16_rp,0.05_rp,rp); ti(2,1)=conjg(ti(1,2)); ti(2,4)=cmplx(-0.07_rp,0.12_rp,rp); ti(4,2)=conjg(ti(2,4)); ti(1,3)=cmplx(0.03_rp,-0.09_rp,rp); ti(3,1)=conjg(ti(1,3))
      tj = 0.0_rp; tj(1,1)=-0.18_rp; tj(2,2)=0.22_rp; tj(3,3)=0.08_rp; tj(4,4)=0.14_rp
      tj(1,3)=cmplx(-0.11_rp,0.04_rp,rp); tj(3,1)=conjg(tj(1,3)); tj(2,3)=cmplx(0.06_rp,0.15_rp,rp); tj(3,2)=conjg(tj(2,3)); tj(2,4)=cmplx(0.04_rp,-0.08_rp,rp); tj(4,2)=conjg(tj(2,4))
      c = 0.0_rp; c(1,1)=0.05_rp; c(2,2)=-0.03_rp; c(3,3)=0.04_rp; c(4,4)=0.02_rp
      c(1,4)=cmplx(0.07_rp,0.02_rp,rp); c(4,1)=conjg(c(1,4)); c(2,3)=cmplx(-0.05_rp,0.06_rp,rp); c(3,2)=conjg(c(2,3))
   end subroutine make_noncommuting_fixture

   subroutine make_gapped_fixture(h,ti,tj,c)
      complex(rp), intent(out) :: h(4,4),ti(4,4),tj(4,4),c(4,4)
      call make_noncommuting_fixture(h,ti,tj,c)
      h = h*0.35_rp
      h(1,1) = -1.10_rp; h(2,2) = -0.70_rp; h(3,3) = 0.60_rp; h(4,4) = 0.95_rp
   end subroutine make_gapped_fixture

   subroutine diagonalize(matrix,values,vectors)
      complex(rp), intent(in) :: matrix(:,:)
      real(rp), intent(out) :: values(:)
      complex(rp), intent(out) :: vectors(:,:)
      complex(rp), allocatable :: work(:)
      real(rp), allocatable :: rwork(:)
      complex(rp) :: copy(size(matrix,1),size(matrix,2)), query(1)
      integer :: n,lwork,info
      external :: zheev
      n=size(values); copy=matrix; allocate(rwork(max(1,3*n-2)))
      call zheev('V','U',n,copy,n,values,query,-1,rwork,info); lwork=max(1,int(real(query(1),rp))); allocate(work(lwork))
      call zheev('V','U',n,copy,n,values,work,lwork,rwork,info); if(info/=0) error stop 'DRESP-03G diagonalization failed'
      vectors=copy
      deallocate(work,rwork)
   end subroutine diagonalize

   function finite_difference_hessian(h,ti,tj,c,mu,kT,step) result(value)
      complex(rp), intent(in) :: h(:,:),ti(:,:),tj(:,:),c(:,:)
      real(rp), intent(in) :: mu,kT,step
      real(rp) :: value,opp,opm,omp,omm
      opp=perturbed_omega(h,ti,tj,c,step,step,mu,kT); opm=perturbed_omega(h,ti,tj,c,step,-step,mu,kT)
      omp=perturbed_omega(h,ti,tj,c,-step,step,mu,kT); omm=perturbed_omega(h,ti,tj,c,-step,-step,mu,kT)
      value=(opp-opm-omp+omm)/(4.0_rp*step*step)
   end function finite_difference_hessian

   function perturbed_omega(h,ti,tj,c,x,y,mu,kT) result(value)
      complex(rp), intent(in) :: h(:,:),ti(:,:),tj(:,:),c(:,:)
      real(rp), intent(in) :: x,y,mu,kT
      real(rp) :: value,levels(size(h,1)),arg
      complex(rp) :: hp(size(h,1),size(h,2)),vec(size(h,1),size(h,2)),query(1),work(32)
      real(rp) :: rwork(32)
      integer :: n,info
      external :: zheev
      n=size(h,1); hp=h+x*ti+y*tj+x*y*c; vec=hp
      call zheev('N','U',n,vec,n,levels,query,-1,rwork,info); call zheev('N','U',n,vec,n,levels,work,32,rwork,info)
      value=0.0_rp
      do info=1,n
         arg=-(levels(info)-mu)/kT
         if(arg>40.0_rp) then; value=value-kT*arg; else; value=value-kT*log(1.0_rp+exp(arg)); end if
      end do
   end function perturbed_omega

   subroutine degeneracy_similarity_test(error,h,ti,tj,c)
      real(rp), intent(out) :: error
      complex(rp), intent(out) :: h(4,4),ti(4,4),tj(4,4),c(4,4)
      complex(rp) :: u(4,4),hp(4,4),tip(4,4),tjp(4,4),cp(4,4),a(1,1),b(1,1),d(1,1),e(1,1), &
         tqa(4,4,1),tma(4,4,1),ca(4,4,1,1),tqb(4,4,1),tmb(4,4,1),cb(4,4,1,1)
      type(finite_h_contour_options) :: local_options
      integer :: j
      call make_noncommuting_fixture(h,ti,tj,c)
      h=0.0_rp; h(1,1)=0.10_rp; h(2,2)=0.10_rp; h(3,3)=-0.20_rp; h(4,4)=0.70_rp
      u=0.0_rp; do j=1,4; u(j,j)=1.0_rp; end do
      u(1,1)=cos(0.37_rp); u(2,1)=sin(0.37_rp); u(1,2)=-sin(0.37_rp); u(2,2)=cos(0.37_rp)
      hp=matmul(u,matmul(h,conjg(transpose(u)))); tip=matmul(u,matmul(ti,conjg(transpose(u)))); tjp=matmul(u,matmul(tj,conjg(transpose(u)))); cp=matmul(u,matmul(c,conjg(transpose(u))))
      tqa(:,:,1)=ti; tma(:,:,1)=tj; ca(:,:,1,1)=c; tqb(:,:,1)=tip; tmb(:,:,1)=tjp; cb(:,:,1,1)=cp
      local_options%contour_points=48; local_options%contour_margin=0.30_rp; local_options%contour_height_fraction=0.35_rp
      call force_theorem_finite_q_hessian_from_resolvent(h,h,0.0_rp,0.03_rp,tqa,tma,ca,local_options,a,b,d,e)
      call force_theorem_finite_q_hessian_from_resolvent(hp,hp,0.0_rp,0.03_rp,tqb,tmb,cb,local_options,a,b,d,e)
      error=0.0_rp
      ! Re-evaluate into distinct variables to compare the invariant result.
      call force_theorem_finite_q_hessian_from_resolvent(h,h,0.0_rp,0.03_rp,tqa,tma,ca,local_options,a,b,d,e)
      call force_theorem_finite_q_hessian_from_resolvent(hp,hp,0.0_rp,0.03_rp,tqb,tmb,cb,local_options, a,b,d,e)
      ! The two calls above deliberately share result slots; use a fresh direct
      ! comparison through the complete scalar saved before the second call.
      error=0.0_rp
      call compare_similarity(h,ti,tj,c,hp,tip,tjp,cp,error)
   end subroutine degeneracy_similarity_test

   subroutine compare_similarity(h,ti,tj,c,hp,tip,tjp,cp,error)
      complex(rp), intent(in) :: h(:,:),ti(:,:),tj(:,:),c(:,:),hp(:,:),tip(:,:),tjp(:,:),cp(:,:)
      real(rp), intent(out) :: error
      complex(rp) :: x(4,4,1),y(4,4,1),z(4,4,1,1),a(1,1),b(1,1),d(1,1),e(1,1),aa(1,1),bb(1,1),dd(1,1),ee(1,1)
      type(finite_h_contour_options) :: o
      x(:,:,1)=ti; y(:,:,1)=tj; z(:,:,1,1)=c; o%contour_points=48; o%contour_margin=0.30_rp; o%contour_height_fraction=0.35_rp
      call force_theorem_finite_q_hessian_from_resolvent(h,h,0.0_rp,0.03_rp,x,y,z,o,a,b,d,e)
      x(:,:,1)=tip; y(:,:,1)=tjp; z(:,:,1,1)=cp
      call force_theorem_finite_q_hessian_from_resolvent(hp,hp,0.0_rp,0.03_rp,x,y,z,o,aa,bb,dd,ee)
      error=maxval(abs([real(e-ee,rp),real(b-bb,rp),real(d-dd,rp)]))
   end subroutine compare_similarity

end program test_dresp03g_contour
