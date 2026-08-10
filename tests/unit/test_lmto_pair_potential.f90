program test_lmto_pair_potential
   use precision_mod, only: rp
   use math_mod, only: i_unit
   use lmto_magnetic_tangent_mod, only: lmto_bond_value, lmto_bond_tangent, lmto_hhmag_to_spinor
   use lmto_pair_potential_mod, only: lmto_circular_pair_potential, lmto_bloch_phase, lmto_endpoint_phases, &
      lmto_pair_transition_metadata, lmto_transition_metadata, lmto_unfold_site_spinors
   implicit none
   real(rp), parameter :: tol = 2.0e-9_rp
   real(rp) :: maximum_error
   logical :: failed

   maximum_error = 0.0_rp; failed = .false.
   call test_rotation_oracle_and_adjoint()
   call test_unequal_orbital_negative_control()
   call test_signed_moment_and_bloch_phase()
   call test_finite_q_endpoint_phases_and_gauge()
   call test_commensurate_supercell_oracle()
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

   subroutine test_finite_q_endpoint_phases_and_gauge()
      real(rp)::k(3),q(3),d(3),tau(3,2)
      complex(rp)::left,right,common,folded(4,1),unfolded(4,1),expected(4,1)
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
   end subroutine test_finite_q_endpoint_phases_and_gauge

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
