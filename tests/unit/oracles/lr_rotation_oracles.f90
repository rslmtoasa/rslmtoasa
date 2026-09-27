!------------------------------------------------------------------------------
! Rotation and force-theorem helpers retained only for independent unit oracles.
!------------------------------------------------------------------------------
module lr_rotation_oracles_mod
   use precision_mod, only: rp
   use math_mod, only: hcpx, pi
   use lmto_magnetic_tangent_mod, only: lmto_bond_value, lmto_bond_derivative, lmto_hhmag_to_spinor
   use linear_response_mod, only: force_theorem_pi, lmto_live_hamiltonian_fixture, &
      assemble_lmto_finite_q_torque, assemble_lmto_finite_q_mixed_derivative
   implicit none
   private

   public :: assemble_lmto_torque
   public :: assemble_lmto_rotation_terms
   public :: assemble_lmto_mixed_derivative
   public :: force_theorem_integrand
   public :: force_theorem_hessian_from_eigenbasis
   public :: mixed_second_difference
   public :: grand_potential_from_eigenvalues

contains

   subroutine assemble_lmto_torque(this, k_point, site, axis, torque)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: k_point(3), axis(3)
      integer, intent(in) :: site
      complex(rp), intent(out) :: torque(:,:)
      real(rp) :: zero_q(3)

      zero_q = 0.0_rp
      call assemble_lmto_finite_q_torque(this, k_point, zero_q, site, axis, torque)
   end subroutine assemble_lmto_torque

   subroutine assemble_lmto_rotation_terms(this, k_point, site, axis, b, q, enu, bi, qi, enui, torque)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: k_point(3), axis(3)
      integer, intent(in) :: site
      complex(rp), intent(out) :: b(:,:), q(:,:), enu(:,:), bi(:,:), qi(:,:), enui(:,:), torque(:,:)
      integer :: nmat
      real(rp) :: zero_q(3)

      call validate_site_axis(this, site, axis)
      nmat = 2*this%norb*this%nsite
      if (any(shape(b) /= [nmat,nmat]) .or. any(shape(q) /= [nmat,nmat]) .or. &
          any(shape(enu) /= [nmat,nmat]) .or. any(shape(bi) /= [nmat,nmat]) .or. &
          any(shape(qi) /= [nmat,nmat]) .or. any(shape(enui) /= [nmat,nmat]) .or. &
          any(shape(torque) /= [nmat,nmat])) then
         error stop 'assemble_lmto_rotation_terms: output shape mismatch'
      end if

      call assemble_base_terms(this, k_point, b, q)
      enu = cmplx(0.0_rp, 0.0_rp, rp)
      if (this%include_enu) call assemble_onsite_coefficient(this, this%enu0, this%enu1, enu)
      zero_q = 0.0_rp
      call assemble_finite_q_directional_terms(this, k_point, zero_q, site, axis, bi, qi)
      enui = cmplx(0.0_rp, 0.0_rp, rp)
      if (this%include_enu) call assemble_finite_q_onsite_derivative(this, zero_q, site, axis, enui)

      torque = bi + enui
      if (this%hoh) torque = torque - matmul(qi, b) - matmul(q, bi)
   end subroutine assemble_lmto_rotation_terms

   subroutine assemble_lmto_mixed_derivative(this, k_point, site_i, axis_i, site_j, axis_j, mixed)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: k_point(3), axis_i(3), axis_j(3)
      integer, intent(in) :: site_i, site_j
      complex(rp), intent(out) :: mixed(:,:)
      real(rp) :: zero_q(3)

      if (site_i == site_j) error stop 'assemble_lmto_mixed_derivative: sites must be distinct'
      zero_q = 0.0_rp
      call assemble_lmto_finite_q_mixed_derivative(this, k_point, zero_q, site_i, axis_i, site_j, axis_j, mixed)
   end subroutine assemble_lmto_mixed_derivative

   pure subroutine force_theorem_integrand(torque_i, green, torque_j, mixed, torque_torque, mixed_contact, complete, include_contact)
      complex(rp), intent(in) :: torque_i(:,:), green(:,:), torque_j(:,:), mixed(:,:)
      real(rp), intent(out) :: torque_torque, mixed_contact, complete
      logical, intent(in), optional :: include_contact
      complex(rp) :: tt_trace, contact_trace
      logical :: with_contact

      if (any(shape(torque_i) /= shape(green)) .or. any(shape(torque_j) /= shape(green)) .or. &
          any(shape(mixed) /= shape(green))) error stop 'force_theorem_integrand: shape mismatch'
      with_contact = .true.
      if (present(include_contact)) with_contact = include_contact
      tt_trace = trace_product(torque_i, green, torque_j, green)
      contact_trace = trace_product(mixed, green)
      torque_torque = -aimag(tt_trace)/force_theorem_pi
      mixed_contact = -aimag(contact_trace)/force_theorem_pi
      complete = torque_torque
      if (with_contact) complete = complete + mixed_contact
   end subroutine force_theorem_integrand

   subroutine force_theorem_hessian_from_eigenbasis(eigenvalues, eigenvectors, fermi, torque_i, torque_j, mixed, &
                                                    torque_torque, mixed_contact, complete)
      real(rp), intent(in) :: eigenvalues(:), fermi
      complex(rp), intent(in) :: eigenvectors(:,:), torque_i(:,:), torque_j(:,:), mixed(:,:)
      real(rp), intent(out) :: torque_torque, mixed_contact, complete
      complex(rp) :: tim, tmj, tjm, tmi, term
      integer :: n, m, nbands

      nbands = size(eigenvalues)
      if (size(eigenvectors,1) /= size(eigenvectors,2) .or. size(eigenvectors,2) /= nbands .or. &
          any(shape(torque_i) /= [nbands,nbands]) .or. any(shape(torque_j) /= [nbands,nbands]) .or. &
          any(shape(mixed) /= [nbands,nbands])) error stop 'force_theorem_hessian_from_eigenbasis: shape mismatch'
      torque_torque = 0.0_rp; mixed_contact = 0.0_rp
      do n = 1, nbands
         if (eigenvalues(n) >= fermi) cycle
         mixed_contact = mixed_contact + real(dot_product(eigenvectors(:,n), matmul(mixed,eigenvectors(:,n))), rp)
         do m = 1, nbands
            if (m == n) cycle
            if (abs(eigenvalues(n)-eigenvalues(m)) <= 100.0_rp*epsilon(1.0_rp)) then
               error stop 'force_theorem_hessian_from_eigenbasis: degenerate occupied response'
            end if
            tim = dot_product(eigenvectors(:,n), matmul(torque_i,eigenvectors(:,m)))
            tmj = dot_product(eigenvectors(:,m), matmul(torque_j,eigenvectors(:,n)))
            tjm = dot_product(eigenvectors(:,n), matmul(torque_j,eigenvectors(:,m)))
            tmi = dot_product(eigenvectors(:,m), matmul(torque_i,eigenvectors(:,n)))
            term = (tim*tmj + tjm*tmi)/cmplx(eigenvalues(n)-eigenvalues(m),0.0_rp,rp)
            torque_torque = torque_torque + real(term, rp)
         end do
      end do
      complete = torque_torque + mixed_contact
   end subroutine force_theorem_hessian_from_eigenbasis

   pure function mixed_second_difference(omega_pp, omega_pm, omega_mp, omega_mm, delta_i, delta_j) result(hessian)
      real(rp), intent(in) :: omega_pp, omega_pm, omega_mp, omega_mm, delta_i, delta_j
      real(rp) :: hessian
      if (delta_i == 0.0_rp .or. delta_j == 0.0_rp) error stop 'mixed_second_difference: zero step'
      hessian = (omega_pp - omega_pm - omega_mp + omega_mm)/(4.0_rp*delta_i*delta_j)
   end function mixed_second_difference

   pure function grand_potential_from_eigenvalues(eigenvalues, fermi) result(omega)
      real(rp), intent(in) :: eigenvalues(:), fermi
      real(rp) :: omega
      integer :: i
      omega = 0.0_rp
      do i = 1, size(eigenvalues)
         if (eigenvalues(i) < fermi) omega = omega + eigenvalues(i) - fermi
      end do
   end function grand_potential_from_eigenvalues

   subroutine assemble_finite_q_directional_terms(this, k_point, q_point, site, axis, db, dq)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: k_point(3), q_point(3), axis(3)
      integer, intent(in) :: site
      complex(rp), intent(out) :: db(:,:), dq(:,:)
      complex(rp) :: bond(2*this%norb,2*this%norb), db_source(2*this%norb,2*this%norb), &
         db_target(2*this%norb,2*this%norb), obar(2*this%norb,2*this%norb), dobar(2*this%norb,2*this%norb), &
         hmag(this%norb,this%norb,4), dhmag(this%norb,this%norb,4), phase_k, phase_source, phase_target
      real(rp) :: dm_source(3), dm_target(3), zero(3)
      integer :: ibond, source, target

      db = cmplx(0.0_rp,0.0_rp,rp); dq = db; zero = 0.0_rp
      do ibond = 1, this%nbond
         source=this%bond_source(ibond); target=this%bond_target(ibond)
         dm_source=zero; dm_target=zero
         if (source == site) dm_source=cross3(axis,this%moments(:,source))
         if (target == site) dm_target=cross3(axis,this%moments(:,target))
         if (maxval(abs(dm_source))+maxval(abs(dm_target)) == 0.0_rp) cycle
         call build_bond(this,ibond,hmag); call fixture_hhmag_to_spinor(this,hmag,bond)
         call build_bond_derivative(this,ibond,dm_source,zero,dhmag); call fixture_hhmag_to_spinor(this,dhmag,db_source)
         call build_bond_derivative(this,ibond,zero,dm_target,dhmag); call fixture_hhmag_to_spinor(this,dhmag,db_target)
         phase_k=bond_phase(this,k_point,ibond)
         phase_source=endpoint_phase(this,q_point,ibond,.false.)
         phase_target=endpoint_phase(this,q_point,ibond,.true.)
         call add_site_block(db,phase_source*db_source+phase_target*db_target,source,target,phase_k,this%norb)
         if (this%hoh) then
            call build_onsite_coefficient(this,this%obar0(:,target),this%obar1(:,target),this%moments(:,target),obar)
            call build_onsite_derivative(this,this%obar1(:,target),this%moments(:,target),dm_target,dobar)
            call add_site_block(dq, phase_source*matmul(db_source,obar) + phase_target*(matmul(db_target,obar)+ &
               matmul(bond,dobar)), source,target,phase_k,this%norb)
         end if
      end do
   end subroutine assemble_finite_q_directional_terms

   subroutine assemble_finite_q_onsite_derivative(this, q_point, site, axis, matrix)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: q_point(3), axis(3)
      integer, intent(in) :: site
      complex(rp), intent(out) :: matrix(:,:)
      complex(rp) :: block(2*this%norb,2*this%norb), phase
      real(rp) :: dm(3)
      matrix=cmplx(0.0_rp,0.0_rp,rp); dm=cross3(axis,this%moments(:,site))
      call build_onsite_derivative(this,this%enu1(:,site),this%moments(:,site),dm,block)
      phase=site_phase(this,q_point,site)
      call add_site_block(matrix,phase*block,site,site,cmplx(1.0_rp,0.0_rp,rp),this%norb)
   end subroutine assemble_finite_q_onsite_derivative

   subroutine assemble_base_terms(this, k_point, b, q)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: k_point(3)
      complex(rp), intent(out) :: b(:,:), q(:,:)
      complex(rp) :: bond(2*this%norb,2*this%norb), hmag(this%norb,this%norb,4), &
         obar(2*this%norb,2*this%norb), phase
      integer :: ibond, source, target

      b=cmplx(0.0_rp,0.0_rp,rp); q=b
      do ibond=1,this%nbond
         source=this%bond_source(ibond); target=this%bond_target(ibond)
         call build_bond(this,ibond,hmag); call fixture_hhmag_to_spinor(this,hmag,bond)
         phase=bond_phase(this,k_point,ibond)
         call add_site_block(b,bond,source,target,phase,this%norb)
         if (this%hoh) then
            call build_onsite_coefficient(this,this%obar0(:,target),this%obar1(:,target),this%moments(:,target),obar)
            call add_site_block(q,matmul(bond,obar),source,target,phase,this%norb)
         end if
      end do
   end subroutine assemble_base_terms

   subroutine build_bond(this, ibond, hmag)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      integer, intent(in) :: ibond
      complex(rp), intent(out) :: hmag(:,:,:)
      integer :: source, target
      source=this%bond_source(ibond); target=this%bond_target(ibond)
      call lmto_bond_value(this%hhh(:,:,ibond),this%wx0(:,source),this%wx1(:,source), &
         this%wx0(:,target),this%wx1(:,target),this%c0(:,source),this%c1(:,source), &
         this%moments(:,source),this%moments(:,target),this%onsite(ibond),hmag)
   end subroutine build_bond

   subroutine fixture_hhmag_to_spinor(this, hhmag, spinor)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      complex(rp), intent(in) :: hhmag(:,:,:)
      complex(rp), intent(out) :: spinor(:,:)
      complex(rp), allocatable :: work(:,:,:)
      integer :: idir
      allocate(work,source=hhmag)
      if (this%cartesian_to_spherical) then
         do idir=1,size(work,3)
            call hcpx(work(:,:,idir),'cart2sph')
         end do
      end if
      call lmto_hhmag_to_spinor(work,spinor)
      deallocate(work)
   end subroutine fixture_hhmag_to_spinor

   subroutine build_bond_derivative(this, ibond, dm_source, dm_target, dhmag)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      integer, intent(in) :: ibond
      real(rp), intent(in) :: dm_source(3), dm_target(3)
      complex(rp), intent(out) :: dhmag(:,:,:)
      integer :: source, target
      source=this%bond_source(ibond); target=this%bond_target(ibond)
      call lmto_bond_derivative(this%hhh(:,:,ibond),this%wx0(:,source),this%wx1(:,source), &
         this%wx0(:,target),this%wx1(:,target),this%c1(:,source),this%moments(:,source), &
         this%moments(:,target),dm_source,dm_target,this%onsite(ibond),dhmag)
   end subroutine build_bond_derivative

   subroutine build_onsite_coefficient(this, c0, c1, moment, spinor)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      complex(rp), intent(in) :: c0(:), c1(:)
      real(rp), intent(in) :: moment(3)
      complex(rp), intent(out) :: spinor(:,:)
      complex(rp) :: zero_hhh(size(c0),size(c0)), zero(size(c0)), hmag(size(c0),size(c0),4)
      zero_hhh=cmplx(0.0_rp,0.0_rp,rp); zero=cmplx(0.0_rp,0.0_rp,rp)
      call lmto_bond_value(zero_hhh,zero,zero,zero,zero,c0,c1,moment,moment,.true.,hmag)
      call fixture_hhmag_to_spinor(this,hmag,spinor)
   end subroutine build_onsite_coefficient

   subroutine build_onsite_derivative(this, c1, moment, dm, spinor)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      complex(rp), intent(in) :: c1(:)
      real(rp), intent(in) :: moment(3), dm(3)
      complex(rp), intent(out) :: spinor(:,:)
      complex(rp) :: zero_hhh(size(c1),size(c1)), zero(size(c1)), dhmag(size(c1),size(c1),4)
      zero_hhh=cmplx(0.0_rp,0.0_rp,rp); zero=cmplx(0.0_rp,0.0_rp,rp)
      call lmto_bond_derivative(zero_hhh,zero,zero,zero,zero,c1,moment,moment,dm,dm,.true.,dhmag)
      call fixture_hhmag_to_spinor(this,dhmag,spinor)
   end subroutine build_onsite_derivative

   subroutine assemble_onsite_coefficient(this, c0, c1, matrix)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      complex(rp), intent(in) :: c0(:,:), c1(:,:)
      complex(rp), intent(out) :: matrix(:,:)
      complex(rp) :: block(2*this%norb,2*this%norb)
      integer :: site
      matrix=cmplx(0.0_rp,0.0_rp,rp)
      do site=1,this%nsite
         call build_onsite_coefficient(this,c0(:,site),c1(:,site),this%moments(:,site),block)
         call add_site_block(matrix,block,site,site,cmplx(1.0_rp,0.0_rp,rp),this%norb)
      end do
   end subroutine assemble_onsite_coefficient

   pure function endpoint_phase(this,q_point,ibond,target_endpoint) result(phase)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: q_point(3)
      integer, intent(in) :: ibond
      logical, intent(in) :: target_endpoint
      complex(rp) :: phase
      real(rp) :: position(3)
      position=0.0_rp
      if (target_endpoint) position=this%bond_vector(:,ibond)
      phase=cmplx(cos(2.0_rp*pi*dot_product(q_point,position)),sin(2.0_rp*pi*dot_product(q_point,position)),rp)
   end function endpoint_phase

   pure function site_phase(this,q_point,site) result(phase)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: q_point(3)
      integer, intent(in) :: site
      complex(rp) :: phase
      phase=cmplx(1.0_rp,0.0_rp,rp)
   end function site_phase

   pure function bond_phase(this,k_point,ibond) result(phase)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: k_point(3)
      integer, intent(in) :: ibond
      complex(rp) :: phase
      real(rp) :: angle
      angle=2.0_rp*pi*dot_product(k_point,this%bond_vector(:,ibond))
      phase=cmplx(cos(angle),sin(angle),rp)
   end function bond_phase

   subroutine add_site_block(matrix, block, source, target, factor, norb)
      complex(rp), intent(inout) :: matrix(:,:)
      complex(rp), intent(in) :: block(:,:), factor
      integer, intent(in) :: source, target, norb
      integer :: i1, i2, j1, j2
      i1=(source-1)*2*norb+1; i2=source*2*norb
      j1=(target-1)*2*norb+1; j2=target*2*norb
      matrix(i1:i2,j1:j2)=matrix(i1:i2,j1:j2)+factor*block
   end subroutine add_site_block

   subroutine validate_site_axis(this, site, axis)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      integer, intent(in) :: site
      real(rp), intent(in) :: axis(3)
      if (this%nsite < 1 .or. this%norb < 1 .or. this%nbond < 1 .or. .not. allocated(this%hhh)) &
         error stop 'DRESP-03T fixture is not initialized'
      if (site < 1 .or. site > this%nsite .or. abs(sqrt(dot_product(axis,axis))-1.0_rp) > 1.0e-10_rp) &
         error stop 'DRESP-03T invalid local rotation site or non-unit axis'
   end subroutine validate_site_axis

   pure function trace_product(a, b, c, d) result(value)
      complex(rp), intent(in) :: a(:,:), b(:,:)
      complex(rp), intent(in), optional :: c(:,:), d(:,:)
      complex(rp) :: value
      complex(rp) :: product(size(a,1),size(a,2))
      integer :: n, i
      n=size(a,1)
      if (present(c)) then
         product=matmul(a,matmul(b,matmul(c,d)))
      else
         product=matmul(a,b)
      end if
      value=sum([(product(i,i),i=1,n)])
   end function trace_product

   pure function cross3(a,b) result(c)
      real(rp), intent(in) :: a(3), b(3)
      real(rp) :: c(3)
      c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
   end function cross3

end module lr_rotation_oracles_mod
