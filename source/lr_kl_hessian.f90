!------------------------------------------------------------------------------
! DRESP-03T -- local-rotation derivatives of the live orthogonal LMTO H.
!------------------------------------------------------------------------------
module lr_kl_hessian_mod
   use, intrinsic :: ieee_arithmetic
   use precision_mod, only: rp
   use math_mod, only: i_unit, pi, inverse_3x3
   use hamiltonian_mod, only: hamiltonian
   use lmto_magnetic_tangent_mod, only: lmto_bond_value, lmto_bond_derivative, &
                                        lmto_bond_mixed_derivative, lmto_hhmag_to_spinor
   implicit none
   private

   real(rp), parameter, public :: force_theorem_pi = pi

   ! A small, explicit representation of the live production inputs.  It is
   ! intentionally bond based: every directed production bond is represented
   ! once, including its Fourier displacement and target-site obarm factor.
   ! The fixture is also useful as a read-only audit seam for production code.
   type, public :: lmto_live_hamiltonian_fixture
      integer :: nsite = 0
      integer :: norb = 0
      integer :: nbond = 0
      logical :: hoh = .false.
      complex(rp), allocatable :: hhh(:, :, :)       ! orbital directed bonds
      integer, allocatable :: bond_source(:), bond_target(:)
      real(rp), allocatable :: bond_vector(:, :)     ! Cartesian lattice vector
      logical, allocatable :: onsite(:)
      real(rp), allocatable :: moments(:, :)         ! (3,nsite), unit moments
      complex(rp), allocatable :: wx0(:, :), wx1(:, :) ! (norb,nsite)
      complex(rp), allocatable :: c0(:, :), c1(:, :)   ! live onsite c channels
      complex(rp), allocatable :: obar0(:, :), obar1(:, :) ! live O channels
      complex(rp), allocatable :: enu0(:, :), enu1(:, :)   ! live e_nu channels
   contains
      procedure :: clear => lmto_fixture_clear
   end type lmto_live_hamiltonian_fixture

   public :: lmto_fixture_init
   public :: lmto_fixture_from_hamiltonian
   public :: assemble_lmto_hamiltonian
   public :: assemble_lmto_torque
   public :: assemble_lmto_mixed_derivative
   public :: force_theorem_integrand
   public :: force_theorem_hessian_from_green
   public :: force_theorem_hessian_from_eigenbasis
   public :: mixed_second_difference
   public :: grand_potential_from_eigenvalues

contains

   subroutine lmto_fixture_init(this, nsite, norb, nbond, hoh)
      class(lmto_live_hamiltonian_fixture), intent(out) :: this
      integer, intent(in) :: nsite, norb, nbond
      logical, intent(in), optional :: hoh

      if (nsite < 1 .or. norb < 1 .or. nbond < 1) error stop 'lmto_fixture_init: invalid dimensions'
      this%nsite = nsite
      this%norb = norb
      this%nbond = nbond
      this%hoh = .false.
      if (present(hoh)) this%hoh = hoh
      allocate(this%hhh(norb, norb, nbond), this%bond_source(nbond), this%bond_target(nbond), &
               this%bond_vector(3, nbond), this%onsite(nbond), this%moments(3, nsite), &
               this%wx0(norb, nsite), this%wx1(norb, nsite), this%c0(norb, nsite), this%c1(norb, nsite), &
               this%obar0(norb, nsite), this%obar1(norb, nsite), this%enu0(norb, nsite), this%enu1(norb, nsite))
      this%hhh = cmplx(0.0_rp, 0.0_rp, rp)
      this%bond_source = 0; this%bond_target = 0; this%bond_vector = 0.0_rp; this%onsite = .false.
      this%moments = 0.0_rp; this%moments(3, :) = 1.0_rp
      this%wx0 = cmplx(0.0_rp, 0.0_rp, rp); this%wx1 = cmplx(0.0_rp, 0.0_rp, rp)
      this%c0 = cmplx(0.0_rp, 0.0_rp, rp); this%c1 = cmplx(0.0_rp, 0.0_rp, rp)
      this%obar0 = cmplx(0.0_rp, 0.0_rp, rp); this%obar1 = cmplx(0.0_rp, 0.0_rp, rp)
      this%enu0 = cmplx(0.0_rp, 0.0_rp, rp); this%enu1 = cmplx(0.0_rp, 0.0_rp, rp)
   end subroutine lmto_fixture_init

   !> Copy the live production Hamiltonian inputs into the audit seam.  The
   !> source object is read only.  In particular, the directed orbital block
   !> and fractional Fourier vector are taken from the same lattice neighbor
   !> tables used by reciprocal_fourier, while c/w/o/e_nu channels are copied
   !> from the transformed symbolic-atom potential.
   subroutine lmto_fixture_from_hamiltonian(source, fixture)
      type(hamiltonian), intent(in) :: source
      type(lmto_live_hamiltonian_fixture), intent(out) :: fixture
      integer :: site, ntype, ia, nr, m, ja, target, it, count, ibond
      real(rp), allocatable :: cralat(:, :), ham_vec(:, :)
      integer :: nn_max_loc

      if (.not. associated(source%lattice) .or. .not. associated(source%charge)) then
         error stop 'lmto_fixture_from_hamiltonian: production object is not wired'
      end if
      count = 0
      do site = 1, source%lattice%nrec
         ntype = source%lattice%ib(site); ia = source%lattice%atlist(ntype); nr = source%lattice%nn(ia,1)
         do m = 1, nr
            if (m == 1) then
               target = site
            else
               ja = source%lattice%nn(ia,m)
               if (ja < 1 .or. ja > source%lattice%kk) cycle
               target = source%lattice%iz(ja)
               if (target < 1 .or. target > source%lattice%nrec) cycle
            end if
            count = count + 1
         end do
      end do
      call lmto_fixture_init(fixture, source%lattice%nrec, size(source%charge%lattice%sbar,1), count, source%hoh)
      do site = 1, fixture%nsite
         ntype = source%lattice%ib(site); ia = source%lattice%atlist(ntype); it = source%lattice%iz(ia)
         fixture%moments(:,site) = source%charge%lattice%symbolic_atoms(it)%potential%mom
         fixture%wx0(:,site) = source%charge%lattice%symbolic_atoms(it)%potential%wx0(1:fixture%norb)
         fixture%wx1(:,site) = source%charge%lattice%symbolic_atoms(it)%potential%wx1(1:fixture%norb)
         if (source%hoh) then
            fixture%c0(:,site) = source%charge%lattice%symbolic_atoms(it)%potential%cex0(1:fixture%norb)
            fixture%c1(:,site) = source%charge%lattice%symbolic_atoms(it)%potential%cex1(1:fixture%norb)
         else
            fixture%c0(:,site) = source%charge%lattice%symbolic_atoms(it)%potential%cx0(1:fixture%norb)
            fixture%c1(:,site) = source%charge%lattice%symbolic_atoms(it)%potential%cx1(1:fixture%norb)
         end if
         fixture%obar0(:,site) = source%charge%lattice%symbolic_atoms(it)%potential%obx0(1:fixture%norb)
         fixture%obar1(:,site) = source%charge%lattice%symbolic_atoms(it)%potential%obx1(1:fixture%norb)
         fixture%enu0(:,site) = 0.5_rp*(source%charge%lattice%symbolic_atoms(it)%potential%cx0(1:fixture%norb) - &
                    source%charge%lattice%symbolic_atoms(it)%potential%cex0(1:fixture%norb))
         fixture%enu1(:,site) = 0.5_rp*(source%charge%lattice%symbolic_atoms(it)%potential%cx1(1:fixture%norb) - &
                                        source%charge%lattice%symbolic_atoms(it)%potential%cex1(1:fixture%norb))
      end do
      allocate(cralat(3,source%lattice%kk)); cralat = source%lattice%cr(:,1:source%lattice%kk)*source%lattice%alat
      ibond = 0
      do site = 1, source%lattice%nrec
         ntype = source%lattice%ib(site); ia = source%lattice%atlist(ntype); nr = source%lattice%nn(ia,1)
         allocate(ham_vec(3,nr)); nn_max_loc = nr
         call source%lattice%clusba(source%lattice%r2, cralat, ia, source%lattice%kk, source%lattice%kk, nn_max_loc, ham_vec)
         it = source%lattice%iz(ia)
         do m = 1, nr
            if (m == 1) then
               target = site
            else
               ja = source%lattice%nn(ia,m)
               if (ja < 1 .or. ja > source%lattice%kk) cycle
               target = source%lattice%iz(ja)
               if (target < 1 .or. target > source%lattice%nrec) cycle
            end if
            ibond = ibond + 1
            fixture%bond_source(ibond) = site; fixture%bond_target(ibond) = target; fixture%onsite(ibond) = (m == 1)
            fixture%hhh(:,:,ibond) = transpose(cmplx(real(source%charge%lattice%sbar(1:fixture%norb,1:fixture%norb,m,source%lattice%num(ia)),rp),0.0_rp,rp))
            if (m == 1) then
               fixture%bond_vector(:,ibond) = 0.0_rp
            else if (source%lattice%a_cart_inv_ready) then
               fixture%bond_vector(:,ibond) = matmul(source%lattice%a_cart_inv, ham_vec(:,m)/source%lattice%alat)
            else
               fixture%bond_vector(:,ibond) = matmul(inverse_3x3(source%lattice%a), ham_vec(:,m)/source%lattice%alat)
            end if
         end do
         deallocate(ham_vec)
      end do
      deallocate(cralat)
      if (ibond /= fixture%nbond) error stop 'lmto_fixture_from_hamiltonian: bond count changed during extraction'
   end subroutine lmto_fixture_from_hamiltonian

   subroutine lmto_fixture_clear(this)
      class(lmto_live_hamiltonian_fixture), intent(inout) :: this
      if (allocated(this%hhh)) deallocate(this%hhh)
      if (allocated(this%bond_source)) deallocate(this%bond_source)
      if (allocated(this%bond_target)) deallocate(this%bond_target)
      if (allocated(this%bond_vector)) deallocate(this%bond_vector)
      if (allocated(this%onsite)) deallocate(this%onsite)
      if (allocated(this%moments)) deallocate(this%moments)
      if (allocated(this%wx0)) deallocate(this%wx0)
      if (allocated(this%wx1)) deallocate(this%wx1)
      if (allocated(this%c0)) deallocate(this%c0)
      if (allocated(this%c1)) deallocate(this%c1)
      if (allocated(this%obar0)) deallocate(this%obar0)
      if (allocated(this%obar1)) deallocate(this%obar1)
      if (allocated(this%enu0)) deallocate(this%enu0)
      if (allocated(this%enu1)) deallocate(this%enu1)
      this%nsite = 0; this%norb = 0; this%nbond = 0; this%hoh = .false.
   end subroutine lmto_fixture_clear

   !> Assemble the same first/second-order orthogonal H used by reciprocal
   !> `ham_only`: H=B-QB+E_nu, with Q=ee*obar in reciprocal space.
   subroutine assemble_lmto_hamiltonian(this, k_point, hamiltonian)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: k_point(3)
      complex(rp), intent(out) :: hamiltonian(:, :)
      complex(rp), allocatable :: b(:, :), q(:, :), enu(:, :)
      integer :: nmat

      call validate_fixture(this)
      nmat = 2*this%norb*this%nsite
      if (any(shape(hamiltonian) /= [nmat, nmat])) error stop 'assemble_lmto_hamiltonian: output shape mismatch'
      allocate(b(nmat,nmat), q(nmat,nmat), enu(nmat,nmat))
      call assemble_base_terms(this, k_point, b, q)
      call assemble_onsite_coefficient(this, this%enu0, this%enu1, enu)
      hamiltonian = b + enu
      if (this%hoh) hamiltonian = hamiltonian - matmul(q, b)
      deallocate(b, q, enu)
   end subroutine assemble_lmto_hamiltonian

   !> Complete T_i(k)=dH_live/dtheta_i at theta=0.  The result is global in
   !> site/orbital/spin space and therefore contains offsite blocks.
   subroutine assemble_lmto_torque(this, k_point, site, axis, torque)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: k_point(3), axis(3)
      integer, intent(in) :: site
      complex(rp), intent(out) :: torque(:, :)
      complex(rp), allocatable :: b(:, :), q(:, :), db(:, :), dq(:, :), denu(:, :)
      integer :: nmat

      call validate_site_axis(this, site, axis)
      nmat = 2*this%norb*this%nsite
      if (any(shape(torque) /= [nmat, nmat])) error stop 'assemble_lmto_torque: output shape mismatch'
      allocate(b(nmat,nmat), q(nmat,nmat), db(nmat,nmat), dq(nmat,nmat), denu(nmat,nmat))
      call assemble_base_terms(this, k_point, b, q)
      call assemble_directional_terms(this, k_point, site, axis, db, dq)
      call assemble_onsite_derivative(this, site, axis, denu)
      torque = db + denu
      if (this%hoh) torque = torque - matmul(dq,b) - matmul(q,db)
      deallocate(b, q, db, dq, denu)
   end subroutine assemble_lmto_torque

   !> Complete C_ij=d2H_live/(dtheta_i dtheta_j), i/=j.  This includes the
   !> product rule for the global Q*B HOH term and is independently testable
   !> against four evaluations of assemble_lmto_hamiltonian.
   subroutine assemble_lmto_mixed_derivative(this, k_point, site_i, axis_i, site_j, axis_j, mixed)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: k_point(3), axis_i(3), axis_j(3)
      integer, intent(in) :: site_i, site_j
      complex(rp), intent(out) :: mixed(:, :)
      complex(rp), allocatable :: b(:, :), q(:, :), dbi(:, :), dbj(:, :), dqi(:, :), dqj(:, :), &
                                  d2b(:, :), d2q(:, :), d2enu(:, :)
      integer :: nmat

      if (site_i == site_j) error stop 'assemble_lmto_mixed_derivative: sites must be distinct'
      call validate_site_axis(this, site_i, axis_i)
      call validate_site_axis(this, site_j, axis_j)
      nmat = 2*this%norb*this%nsite
      if (any(shape(mixed) /= [nmat, nmat])) error stop 'assemble_lmto_mixed_derivative: output shape mismatch'
      allocate(b(nmat,nmat), q(nmat,nmat), dbi(nmat,nmat), dbj(nmat,nmat), dqi(nmat,nmat), dqj(nmat,nmat), &
               d2b(nmat,nmat), d2q(nmat,nmat), d2enu(nmat,nmat))
      call assemble_base_terms(this, k_point, b, q)
      call assemble_directional_terms(this, k_point, site_i, axis_i, dbi, dqi)
      call assemble_directional_terms(this, k_point, site_j, axis_j, dbj, dqj)
      call assemble_mixed_terms(this, k_point, site_i, axis_i, site_j, axis_j, d2b, d2q)
      d2enu = cmplx(0.0_rp, 0.0_rp, rp)
      mixed = d2b + d2enu
      if (this%hoh) mixed = mixed - matmul(d2q,b) - matmul(dqi,dbj) - matmul(dqj,dbi) - matmul(q,d2b)
      deallocate(b, q, dbi, dbj, dqi, dqj, d2b, d2q, d2enu)
   end subroutine assemble_lmto_mixed_derivative

   !> Integrand of the exact zero-temperature grand-potential Hessian in the
   !> resolvent convention G=(E-H)^-1.  The contact term is C_ij G; it is not
   !> optional physics, only an optional diagnostic split.
   pure subroutine force_theorem_integrand(torque_i, green, torque_j, mixed, torque_torque, mixed_contact, complete, include_contact)
      complex(rp), intent(in) :: torque_i(:, :), green(:, :), torque_j(:, :), mixed(:, :)
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

   !> Integrate the force-theorem integrand using caller-supplied contour
   !> weights.  The same routine returns A (TT only), B (C contact), and A+B.
   subroutine force_theorem_hessian_from_green(greens, weights, torque_i, torque_j, mixed, &
                                                torque_torque, mixed_contact, complete, include_contact)
      complex(rp), intent(in) :: greens(:, :, :), torque_i(:, :), torque_j(:, :), mixed(:, :)
      real(rp), intent(in) :: weights(:)
      real(rp), intent(out) :: torque_torque, mixed_contact, complete
      logical, intent(in), optional :: include_contact
      integer :: ie
      real(rp) :: tt, cc, all_terms

      if (size(greens,3) /= size(weights) .or. size(greens,1) /= size(greens,2) .or. &
          any(shape(torque_i) /= shape(greens(:,:,1))) .or. any(shape(torque_j) /= shape(greens(:,:,1))) .or. &
          any(shape(mixed) /= shape(greens(:,:,1)))) error stop 'force_theorem_hessian_from_green: shape mismatch'
      torque_torque = 0.0_rp; mixed_contact = 0.0_rp; complete = 0.0_rp
      do ie = 1, size(weights)
         call force_theorem_integrand(torque_i, greens(:,:,ie), torque_j, mixed, tt, cc, all_terms, include_contact)
         torque_torque = torque_torque + weights(ie)*tt
         mixed_contact = mixed_contact + weights(ie)*cc
         complete = complete + weights(ie)*all_terms
      end do
   end subroutine force_theorem_hessian_from_green

   !> Spectral form of the same theorem, useful for a finite fixture with a
   !> Fermi level in a gap.  It avoids any quadrature or finite-difference
   !> reuse: the TT term is the ordinary first-order eigenvector response and
   !> the contact term is the diagonal expectation of C_ij.
   subroutine force_theorem_hessian_from_eigenbasis(eigenvalues, eigenvectors, fermi, torque_i, torque_j, mixed, &
                                                    torque_torque, mixed_contact, complete)
      real(rp), intent(in) :: eigenvalues(:), fermi
      complex(rp), intent(in) :: eigenvectors(:, :), torque_i(:, :), torque_j(:, :), mixed(:, :)
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

   subroutine assemble_base_terms(this, k_point, b, q)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: k_point(3)
      complex(rp), intent(out) :: b(:, :), q(:, :)
      complex(rp) :: bond(2*this%norb,2*this%norb), hmag(this%norb,this%norb,4), &
                     obar(2*this%norb,2*this%norb), phase
      integer :: ibond, source, target

      b = cmplx(0.0_rp, 0.0_rp, rp); q = b
      do ibond = 1, this%nbond
         source = this%bond_source(ibond); target = this%bond_target(ibond)
         call build_bond(this, ibond, hmag)
         call lmto_hhmag_to_spinor(hmag, bond)
         phase = cmplx(cos(2.0_rp*pi*dot_product(k_point,this%bond_vector(:,ibond))), &
                       sin(2.0_rp*pi*dot_product(k_point,this%bond_vector(:,ibond))), rp)
         call add_site_block(b, bond, source, target, phase, this%norb)
         if (this%hoh) then
            call build_onsite_coefficient(this%obar0(:,target), this%obar1(:,target), this%moments(:,target), obar)
            call add_site_block(q, matmul(bond, obar), source, target, phase, this%norb)
         end if
      end do
   end subroutine assemble_base_terms

   subroutine assemble_directional_terms(this, k_point, site, axis, db, dq)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: k_point(3), axis(3)
      integer, intent(in) :: site
      complex(rp), intent(out) :: db(:, :), dq(:, :)
      complex(rp) :: bond(2*this%norb,2*this%norb), dbond(2*this%norb,2*this%norb), &
                     hmag(this%norb,this%norb,4), dhmag(this%norb,this%norb,4), &
                     obar(2*this%norb,2*this%norb), dobar(2*this%norb,2*this%norb), phase
      real(rp) :: dm_source(3), dm_target(3)
      integer :: ibond, source, target

      db = cmplx(0.0_rp, 0.0_rp, rp); dq = db
      do ibond = 1, this%nbond
         source = this%bond_source(ibond); target = this%bond_target(ibond)
         call endpoint_variation(this, source, target, site, axis, dm_source, dm_target)
         if (maxval(abs(dm_source)) == 0.0_rp .and. maxval(abs(dm_target)) == 0.0_rp) cycle
         call build_bond(this, ibond, hmag)
         call lmto_hhmag_to_spinor(hmag, bond)
         call build_bond_derivative(this, ibond, dm_source, dm_target, dhmag)
         call lmto_hhmag_to_spinor(dhmag, dbond)
         phase = bond_phase(this, k_point, ibond)
         call add_site_block(db, dbond, source, target, phase, this%norb)
         if (this%hoh) then
            call build_onsite_coefficient(this%obar0(:,target), this%obar1(:,target), this%moments(:,target), obar)
            call build_onsite_derivative(this%obar1(:,target), this%moments(:,target), dm_target, dobar)
            call add_site_block(dq, matmul(dbond,obar) + matmul(bond,dobar), source,target,phase,this%norb)
         end if
      end do
   end subroutine assemble_directional_terms

   subroutine assemble_mixed_terms(this, k_point, site_i, axis_i, site_j, axis_j, d2b, d2q)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: k_point(3), axis_i(3), axis_j(3)
      integer, intent(in) :: site_i, site_j
      complex(rp), intent(out) :: d2b(:, :), d2q(:, :)
      complex(rp) :: d2bond(2*this%norb,2*this%norb), bond(2*this%norb,2*this%norb), &
                     dbond_i(2*this%norb,2*this%norb), dbond_j(2*this%norb,2*this%norb), &
                     hmag(this%norb,this%norb,4), dhmag(this%norb,this%norb,4), &
                     d2hmag(this%norb,this%norb,4), obar(2*this%norb,2*this%norb), &
                     dobar_i(2*this%norb,2*this%norb), dobar_j(2*this%norb,2*this%norb), phase
      real(rp) :: dm_si(3), dm_ti(3), dm_sj(3), dm_tj(3)
      integer :: ibond, source, target

      d2b = cmplx(0.0_rp, 0.0_rp, rp); d2q = d2b
      do ibond = 1, this%nbond
         source = this%bond_source(ibond); target = this%bond_target(ibond)
         call endpoint_variation(this, source, target, site_i, axis_i, dm_si, dm_ti)
         call endpoint_variation(this, source, target, site_j, axis_j, dm_sj, dm_tj)
         if (maxval(abs(dm_si))+maxval(abs(dm_ti))+maxval(abs(dm_sj))+maxval(abs(dm_tj)) == 0.0_rp) cycle
         call build_bond(this, ibond, hmag)
         call lmto_hhmag_to_spinor(hmag, bond)
         call build_bond_derivative(this, ibond, dm_si, dm_ti, dhmag)
         call lmto_hhmag_to_spinor(dhmag, dbond_i)
         call build_bond_derivative(this, ibond, dm_sj, dm_tj, dhmag)
         call lmto_hhmag_to_spinor(dhmag, dbond_j)
         call lmto_bond_mixed_derivative(this%hhh(:,:,ibond), this%wx0(:,source), this%wx1(:,source), &
            this%wx0(:,target), this%wx1(:,target), this%moments(:,source), this%moments(:,target), &
            dm_si, dm_ti, dm_sj, dm_tj, d2hmag)
         call lmto_hhmag_to_spinor(d2hmag, d2bond)
         phase = bond_phase(this, k_point, ibond)
         call add_site_block(d2b, d2bond, source, target, phase, this%norb)
         if (this%hoh) then
            call build_onsite_coefficient(this%obar0(:,target), this%obar1(:,target), this%moments(:,target), obar)
            call build_onsite_derivative(this%obar1(:,target), this%moments(:,target), dm_ti, dobar_i)
            call build_onsite_derivative(this%obar1(:,target), this%moments(:,target), dm_tj, dobar_j)
            call add_site_block(d2q, matmul(d2bond,obar) + matmul(dbond_i,dobar_j) + &
               matmul(dbond_j,dobar_i), source, target, phase, this%norb)
         end if
      end do
   end subroutine assemble_mixed_terms

   subroutine build_bond(this, ibond, hmag)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      integer, intent(in) :: ibond
      complex(rp), intent(out) :: hmag(:, :, :)
      complex(rp) :: zero(this%norb)
      integer :: source, target
      source = this%bond_source(ibond); target = this%bond_target(ibond)
      zero = cmplx(0.0_rp, 0.0_rp, rp)
      call lmto_bond_value(this%hhh(:,:,ibond), this%wx0(:,source), this%wx1(:,source), &
         this%wx0(:,target), this%wx1(:,target), this%c0(:,source), this%c1(:,source), &
         this%moments(:,source), this%moments(:,target), this%onsite(ibond), hmag)
   end subroutine build_bond

   subroutine build_bond_derivative(this, ibond, dm_source, dm_target, dhmag)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      integer, intent(in) :: ibond
      real(rp), intent(in) :: dm_source(3), dm_target(3)
      complex(rp), intent(out) :: dhmag(:, :, :)
      integer :: source, target
      source = this%bond_source(ibond); target = this%bond_target(ibond)
      call lmto_bond_derivative(this%hhh(:,:,ibond), this%wx0(:,source), this%wx1(:,source), &
         this%wx0(:,target), this%wx1(:,target), this%c1(:,source), this%moments(:,source), &
         this%moments(:,target), dm_source, dm_target, this%onsite(ibond), dhmag)
   end subroutine build_bond_derivative

   subroutine build_onsite_coefficient(c0, c1, moment, spinor)
      complex(rp), intent(in) :: c0(:), c1(:)
      real(rp), intent(in) :: moment(3)
      complex(rp), intent(out) :: spinor(:, :)
      complex(rp) :: zero_hhh(size(c0),size(c0)), zero(size(c0)), hmag(size(c0),size(c0),4)
      zero_hhh = cmplx(0.0_rp, 0.0_rp, rp); zero = cmplx(0.0_rp, 0.0_rp, rp)
      call lmto_bond_value(zero_hhh, zero, zero, zero, zero, c0, c1, moment, moment, .true., hmag)
      call lmto_hhmag_to_spinor(hmag, spinor)
   end subroutine build_onsite_coefficient

   subroutine build_onsite_derivative(c1, moment, dm, spinor)
      complex(rp), intent(in) :: c1(:)
      real(rp), intent(in) :: moment(3), dm(3)
      complex(rp), intent(out) :: spinor(:, :)
      complex(rp) :: zero_hhh(size(c1),size(c1)), zero(size(c1)), dhmag(size(c1),size(c1),4)
      zero_hhh = cmplx(0.0_rp, 0.0_rp, rp); zero = cmplx(0.0_rp, 0.0_rp, rp)
      call lmto_bond_derivative(zero_hhh, zero, zero, zero, zero, c1, moment, moment, dm, dm, .true., dhmag)
      call lmto_hhmag_to_spinor(dhmag, spinor)
   end subroutine build_onsite_derivative

   subroutine assemble_onsite_coefficient(this, c0, c1, matrix)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      complex(rp), intent(in) :: c0(:, :), c1(:, :)
      complex(rp), intent(out) :: matrix(:, :)
      complex(rp) :: block(2*this%norb,2*this%norb)
      integer :: site
      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, this%nsite
         call build_onsite_coefficient(c0(:,site), c1(:,site), this%moments(:,site), block)
         call add_site_block(matrix, block, site, site, cmplx(1.0_rp,0.0_rp,rp), this%norb)
      end do
   end subroutine assemble_onsite_coefficient

   subroutine assemble_onsite_derivative(this, site, axis, matrix)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      integer, intent(in) :: site
      real(rp), intent(in) :: axis(3)
      complex(rp), intent(out) :: matrix(:, :)
      complex(rp) :: block(2*this%norb,2*this%norb)
      real(rp) :: dm(3)
      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      dm = cross3(axis, this%moments(:,site))
      call build_onsite_derivative(this%enu1(:,site), this%moments(:,site), dm, block)
      call add_site_block(matrix, block, site, site, cmplx(1.0_rp,0.0_rp,rp), this%norb)
   end subroutine assemble_onsite_derivative

   subroutine endpoint_variation(this, source, target, site, axis, dm_source, dm_target)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      integer, intent(in) :: source, target, site
      real(rp), intent(in) :: axis(3)
      real(rp), intent(out) :: dm_source(3), dm_target(3)
      dm_source = 0.0_rp; dm_target = 0.0_rp
      if (source == site) dm_source = cross3(axis, this%moments(:,source))
      if (target == site) dm_target = cross3(axis, this%moments(:,target))
   end subroutine endpoint_variation

   pure function bond_phase(this, k_point, ibond) result(phase)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: k_point(3)
      integer, intent(in) :: ibond
      complex(rp) :: phase
      real(rp) :: angle
      angle = 2.0_rp*pi*dot_product(k_point, this%bond_vector(:,ibond))
      phase = cmplx(cos(angle), sin(angle), rp)
   end function bond_phase

   subroutine add_site_block(matrix, block, source, target, factor, norb)
      complex(rp), intent(inout) :: matrix(:, :)
      complex(rp), intent(in) :: block(:, :), factor
      integer, intent(in) :: source, target, norb
      integer :: i1, i2, j1, j2
      i1 = (source-1)*2*norb+1; i2 = source*2*norb
      j1 = (target-1)*2*norb+1; j2 = target*2*norb
      matrix(i1:i2,j1:j2) = matrix(i1:i2,j1:j2) + factor*block
   end subroutine add_site_block

   subroutine validate_fixture(this)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      integer :: ibond
      if (this%nsite < 1 .or. this%norb < 1 .or. this%nbond < 1 .or. .not. allocated(this%hhh)) &
         error stop 'DRESP-03T fixture is not initialized'
      if (any(shape(this%moments) /= [3,this%nsite]) .or. any(shape(this%wx0) /= [this%norb,this%nsite]) .or. &
          any(shape(this%wx1) /= [this%norb,this%nsite])) error stop 'DRESP-03T fixture channel shape mismatch'
      do ibond = 1, this%nbond
         if (this%bond_source(ibond) < 1 .or. this%bond_source(ibond) > this%nsite .or. &
             this%bond_target(ibond) < 1 .or. this%bond_target(ibond) > this%nsite) error stop 'DRESP-03T invalid bond endpoint'
      end do
   end subroutine validate_fixture

   subroutine validate_site_axis(this, site, axis)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      integer, intent(in) :: site
      real(rp), intent(in) :: axis(3)
      call validate_fixture(this)
      if (site < 1 .or. site > this%nsite .or. abs(sqrt(dot_product(axis,axis))-1.0_rp) > 1.0e-10_rp) &
         error stop 'DRESP-03T invalid local rotation site or non-unit axis'
   end subroutine validate_site_axis

   pure function trace_product(a, b, c, d) result(value)
      complex(rp), intent(in) :: a(:, :), b(:, :)
      complex(rp), intent(in), optional :: c(:, :), d(:, :)
      complex(rp) :: value
      complex(rp) :: product(size(a,1),size(a,2))
      integer :: n, i
      n = size(a,1)
      if (present(c)) then
         product = matmul(a,matmul(b,matmul(c,d)))
         value = sum([(product(i,i), i=1,n)])
      else
         value = sum([(a(i,i), i=1,n)])
      end if
   end function trace_product

   pure function cross3(a, b) result(c)
      real(rp), intent(in) :: a(3), b(3)
      real(rp) :: c(3)
      c = [a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), a(1)*b(2)-a(2)*b(1)]
   end function cross3

end module lr_kl_hessian_mod
