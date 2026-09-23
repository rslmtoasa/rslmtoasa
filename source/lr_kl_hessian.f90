!------------------------------------------------------------------------------
! DRESP-03T -- local-rotation derivatives of the live orthogonal LMTO H.
!------------------------------------------------------------------------------
module lr_kl_hessian_mod
   use, intrinsic :: ieee_arithmetic
   use precision_mod, only: rp
   use math_mod, only: i_unit, pi, inverse_3x3, hcpx
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
      logical :: include_enu = .true.
      logical :: cartesian_to_spherical = .false.
      complex(rp), allocatable :: hhh(:, :, :)       ! orbital directed bonds
      integer, allocatable :: bond_source(:), bond_target(:)
      real(rp), allocatable :: bond_vector(:, :)     ! Cartesian lattice vector
      logical, allocatable :: onsite(:)
      real(rp), allocatable :: moments(:, :)         ! (3,nsite), unit moments
      ! Fractional basis positions.  A bond vector is the full target-source
      ! displacement R+tau_target-tau_source.  The endpoint phases below are
      ! relative to the source endpoint because the Bloch basis already carries
      ! the absolute basis position.
      real(rp), allocatable :: site_position(:, :)   ! (3,nsite), fractional
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
   public :: assemble_lmto_rotation_terms
   public :: assemble_lmto_mixed_derivative
   public :: assemble_lmto_finite_q_torque
   public :: assemble_lmto_finite_q_torques
   public :: assemble_lmto_finite_q_mixed_derivative
   public :: force_theorem_finite_q_hessian_from_eigenbasis
   public :: force_theorem_finite_q_hessian_from_eigenbasis_batch
   public :: force_theorem_finite_q_hessian_from_eigenbasis_metallic
   public :: force_theorem_finite_q_hessian_from_eigenbasis_metallic_batch
   public :: finite_temperature_occupation
   public :: fermi_divided_difference
   public :: lmto_fixture_adapter_residual
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
      this%include_enu = .true.
      this%cartesian_to_spherical = .false.
      allocate(this%hhh(norb, norb, nbond), this%bond_source(nbond), this%bond_target(nbond), &
               this%bond_vector(3, nbond), this%onsite(nbond), this%moments(3, nsite), &
               this%site_position(3, nsite), &
               this%wx0(norb, nsite), this%wx1(norb, nsite), this%c0(norb, nsite), this%c1(norb, nsite), &
               this%obar0(norb, nsite), this%obar1(norb, nsite), this%enu0(norb, nsite), this%enu1(norb, nsite))
      this%hhh = cmplx(0.0_rp, 0.0_rp, rp)
      this%bond_source = 0; this%bond_target = 0; this%bond_vector = 0.0_rp; this%onsite = .false.
      this%moments = 0.0_rp; this%moments(3, :) = 1.0_rp
      this%site_position = 0.0_rp
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
      ! The reciprocal auto mode is first order for an orthogonal Hamiltonian
      ! and second order only when HOH is active.  Keep the extracted fixture
      ! on that same ham_only convention; synthetic DRESP fixtures retain the
      ! historical include_enu=.true. default from lmto_fixture_init.
      fixture%include_enu = source%hoh
      fixture%cartesian_to_spherical = .true.
      do site = 1, fixture%nsite
         ntype = source%lattice%ib(site); ia = source%lattice%atlist(ntype); it = source%lattice%iz(ia)
         fixture%moments(:,site) = source%charge%lattice%symbolic_atoms(it)%potential%mom
         fixture%site_position(:,site) = source%lattice%cr(:,ia)
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
         ! cx0/cx1 and cex0/cex1 already are the spin average/difference;
         ! their direct difference is E_nu (do not halve it again).
         fixture%enu0(:,site) = source%charge%lattice%symbolic_atoms(it)%potential%cx0(1:fixture%norb) - &
                    source%charge%lattice%symbolic_atoms(it)%potential%cex0(1:fixture%norb)
         fixture%enu1(:,site) = source%charge%lattice%symbolic_atoms(it)%potential%cx1(1:fixture%norb) - &
                                        source%charge%lattice%symbolic_atoms(it)%potential%cex1(1:fixture%norb)
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
      if (allocated(this%site_position)) deallocate(this%site_position)
      if (allocated(this%wx0)) deallocate(this%wx0)
      if (allocated(this%wx1)) deallocate(this%wx1)
      if (allocated(this%c0)) deallocate(this%c0)
      if (allocated(this%c1)) deallocate(this%c1)
      if (allocated(this%obar0)) deallocate(this%obar0)
      if (allocated(this%obar1)) deallocate(this%obar1)
      if (allocated(this%enu0)) deallocate(this%enu0)
      if (allocated(this%enu1)) deallocate(this%enu1)
      this%nsite = 0; this%norb = 0; this%nbond = 0; this%hoh = .false.; this%include_enu = .true.; &
         this%cartesian_to_spherical = .false.
   end subroutine lmto_fixture_clear

   !> Assemble the same first/second-order orthogonal H used by reciprocal
   !> `ham_only`: H=B in first order, and H=B-QB+E_nu in second order,
   !> with Q=ee*obar in reciprocal space.  The production fixture records
   !> which convention is active so the first-order path does not acquire an
   !> E_nu term that reciprocal_fourier deliberately omits.
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
      enu = cmplx(0.0_rp, 0.0_rp, rp)
      if (this%include_enu) call assemble_onsite_coefficient(this, this%enu0, this%enu1, enu)
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
      real(rp) :: zero_q(3)
      integer :: nmat

      zero_q = 0.0_rp
      call assemble_lmto_finite_q_torque(this, k_point, zero_q, site, axis, torque)
   end subroutine assemble_lmto_torque

   !> Expose the live fixed-q product-rule terms for a rotation audit.
   !>
   !> The returned matrices are exactly the objects used by the certified
   !> finite-H torque path:
   !>   H = B - Q B + E_nu,
   !>   T_i = B_i - Q_i B - Q B_i + E_nu,i.
   !>
   !> Here Q is the assembled reciprocal `eeo` operator, i.e. the Fourier sum
   !> of `ee(R)*obarm_target`; it is not potential%qpar.  This routine is a
   !> read-only diagnostic seam and deliberately does not expose or mutate any
   !> production Hamiltonian state.
   subroutine assemble_lmto_rotation_terms(this, k_point, site, axis, b, q, enu, bi, qi, enui, torque)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: k_point(3), axis(3)
      integer, intent(in) :: site
      complex(rp), intent(out) :: b(:, :), q(:, :), enu(:, :), bi(:, :), qi(:, :), enui(:, :), torque(:, :)
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

   !> Complete C_ij=d2H_live/(dtheta_i dtheta_j), i/=j.  This includes the
   !> product rule for the global Q*B HOH term and is independently testable
   !> against four evaluations of assemble_lmto_hamiltonian.
   subroutine assemble_lmto_mixed_derivative(this, k_point, site_i, axis_i, site_j, axis_j, mixed)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: k_point(3), axis_i(3), axis_j(3)
      integer, intent(in) :: site_i, site_j
      complex(rp), intent(out) :: mixed(:, :)
      real(rp) :: zero_q(3)
      if (site_i == site_j) error stop 'assemble_lmto_mixed_derivative: sites must be distinct'
      zero_q = 0.0_rp
      call assemble_lmto_finite_q_mixed_derivative(this, k_point, zero_q, site_i, axis_i, site_j, axis_j, mixed)
   end subroutine assemble_lmto_mixed_derivative

   !> First finite-q vertex for the convention
   !> theta_(aR)=theta_a(q) exp(+i 2*pi*q.(R+tau_a)).  After factoring the
   !> production Bloch basis, the returned matrix is the k -> k+q block and
   !> the target endpoint carries the relative phase exp(+i q.D), with
   !> D=R+tau_target-tau_source.  Its second argument is deliberately not folded:
   !> reciprocal-vector changes are retained as the corresponding basis gauge.
   subroutine assemble_lmto_finite_q_torque(this, k_point, q_point, site, axis, torque)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: k_point(3), q_point(3), axis(3)
      integer, intent(in) :: site
      complex(rp), intent(out) :: torque(:, :)
      complex(rp), allocatable :: b(:, :), b_right(:, :), qmat(:, :), db(:, :), dq(:, :), denu(:, :)
      integer :: nmat

      call validate_site_axis(this, site, axis)
      nmat = 2*this%norb*this%nsite
      if (any(shape(torque) /= [nmat,nmat])) error stop 'assemble_lmto_finite_q_torque: output shape mismatch'
      allocate(b(nmat,nmat), b_right(nmat,nmat), qmat(nmat,nmat), db(nmat,nmat), dq(nmat,nmat), denu(nmat,nmat))
      call assemble_base_terms(this, k_point, b, qmat)
      call assemble_base_terms(this, k_point + q_point, b_right, qmat)
      call assemble_finite_q_directional_terms(this, k_point, q_point, site, axis, db, dq)
      if (this%include_enu) then
         call assemble_finite_q_onsite_derivative(this, q_point, site, axis, denu)
      else
         denu = cmplx(0.0_rp, 0.0_rp, rp)
      end if
      torque = db + denu
      ! For the k -> k+q matrix element of QB, the differentiated Q acts
      ! before the source-side B(k); the undifferentiated Q is evaluated at
      ! the row endpoint and multiplies the finite-q derivative of B.
      if (this%hoh) then
         torque = torque - matmul(dq,b) - matmul(qmat,db)
      end if
      deallocate(b,b_right,qmat,db,dq,denu)
   end subroutine assemble_lmto_finite_q_torque

   !> Assemble all site vertices at once.  The site dimension is the final
   !> dimension and is convenient for the q-Hessian contraction.
   subroutine assemble_lmto_finite_q_torques(this, k_point, q_point, axes, torques)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: k_point(3), q_point(3), axes(:, :)
      complex(rp), intent(out) :: torques(:, :, :)
      integer :: site, nmat
      nmat = 2*this%norb*this%nsite
      if (size(axes,1) /= 3 .or. size(axes,2) /= this%nsite .or. any(shape(torques) /= [nmat,nmat,this%nsite])) then
         error stop 'assemble_lmto_finite_q_torques: shape mismatch'
      end if
      do site = 1, this%nsite
         call assemble_lmto_finite_q_torque(this, k_point, q_point, site, axes(:,site), torques(:,:,site))
      end do
   end subroutine assemble_lmto_finite_q_torques

   !> Complete C_ab(k;q,-q), including same-site rotation curvature and the
   !> two distinct momentum orderings in the HOH product rule.
   subroutine assemble_lmto_finite_q_mixed_derivative(this, k_point, q_point, site_i, axis_i, site_j, axis_j, mixed)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: k_point(3), q_point(3), axis_i(3), axis_j(3)
      integer, intent(in) :: site_i, site_j
      complex(rp), intent(out) :: mixed(:, :)
      complex(rp), allocatable :: b(:, :), qmat(:, :), dbi(:, :), dqi(:, :), dbj_minus(:, :), dqj_minus(:, :), &
         dbi_left(:, :), dqi_left(:, :), dbj_at_k(:, :), dqj_at_k(:, :), d2b(:, :), d2q(:, :), d2enu(:, :)
      integer :: nmat

      call validate_site_axis(this, site_i, axis_i)
      call validate_site_axis(this, site_j, axis_j)
      nmat = 2*this%norb*this%nsite
      if (any(shape(mixed) /= [nmat,nmat])) error stop 'assemble_lmto_finite_q_mixed_derivative: output shape mismatch'
      allocate(b(nmat,nmat), qmat(nmat,nmat), dbi(nmat,nmat), dqi(nmat,nmat), dbj_minus(nmat,nmat), dqj_minus(nmat,nmat), &
         dbi_left(nmat,nmat), dqi_left(nmat,nmat), dbj_at_k(nmat,nmat), dqj_at_k(nmat,nmat), &
         d2b(nmat,nmat), d2q(nmat,nmat), d2enu(nmat,nmat))
      call assemble_base_terms(this, k_point, b, qmat)
      call assemble_finite_q_directional_terms(this, k_point, q_point, site_i, axis_i, dbi, dqi)
      call assemble_finite_q_directional_terms(this, k_point + q_point, -q_point, site_j, axis_j, dbj_minus, dqj_minus)
      call assemble_finite_q_directional_terms(this, k_point, -q_point, site_j, axis_j, dbj_at_k, dqj_at_k)
      call assemble_finite_q_directional_terms(this, k_point - q_point, q_point, site_i, axis_i, dbi_left, dqi_left)
      call assemble_finite_q_mixed_base_terms(this, k_point, q_point, site_i, axis_i, site_j, axis_j, d2b, d2q, d2enu)
      mixed = d2b + d2enu
      if (this%hoh) mixed = mixed - matmul(d2q,b) - matmul(dqi_left,dbj_at_k) - matmul(dqj_minus,dbi) - matmul(qmat,d2b)
      deallocate(b,qmat,dbi,dqi,dbj_minus,dqj_minus,dbi_left,dqi_left,dbj_at_k,dqj_at_k,d2b,d2q,d2enu)
   end subroutine assemble_lmto_finite_q_mixed_derivative

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

   !> Finite-q spectral Hessian for one k -> k+q endpoint pair.  The first
   !> vertex stack is evaluated at k for +q and the second at k+q for -q.
   !> The second ordering is retained explicitly; it is needed for a Hermitian
   !> sublattice Hessian when a and b are different.
   subroutine force_theorem_finite_q_hessian_from_eigenbasis(eigenvalues, eigenvectors, endpoint_values, endpoint_vectors, &
      fermi, torques_q, torques_minus_q, mixed, hessian, torque_torque, mixed_contact, complete)
      real(rp), intent(in) :: eigenvalues(:), endpoint_values(:), fermi
      complex(rp), intent(in) :: eigenvectors(:, :), endpoint_vectors(:, :)
      complex(rp), intent(in) :: torques_q(:, :, :), torques_minus_q(:, :, :), mixed(:, :, :, :)
      complex(rp), intent(out) :: hessian(:, :), torque_torque(:, :), mixed_contact(:, :), complete(:, :)
      complex(rp) :: ta_nm, tb_mn, tb_nm, ta_mn, term
      integer :: n, m, a, b, nbands, nsite, nmat

      nbands = size(eigenvalues); nmat = size(eigenvectors,1); nsite = size(torques_q,3)
      if (size(endpoint_values) /= nbands .or. size(eigenvectors,2) /= nbands .or. &
          any(shape(endpoint_vectors) /= [nmat,nbands]) .or. any(shape(torques_q) /= [nmat,nmat,nsite]) .or. &
          any(shape(torques_minus_q) /= [nmat,nmat,nsite]) .or. any(shape(mixed) /= [nmat,nmat,nsite,nsite]) .or. &
          any(shape(hessian) /= [nsite,nsite]) .or. any(shape(torque_torque) /= [nsite,nsite]) .or. &
          any(shape(mixed_contact) /= [nsite,nsite]) .or. any(shape(complete) /= [nsite,nsite])) then
         error stop 'force_theorem_finite_q_hessian_from_eigenbasis: shape mismatch'
      end if
      hessian = 0.0_rp; torque_torque = 0.0_rp; mixed_contact = 0.0_rp; complete = 0.0_rp
      do a = 1, nsite
         do b = 1, nsite
            do n = 1, nbands
               if (eigenvalues(n) >= fermi) cycle
               mixed_contact(a,b) = mixed_contact(a,b) + real(dot_product(eigenvectors(:,n), &
                  matmul(mixed(:,:,a,b),eigenvectors(:,n))),rp)
               do m = 1, nbands
                  if (abs(eigenvalues(n)-endpoint_values(m)) <= 100.0_rp*epsilon(1.0_rp)) then
                     ! A degenerate subspace wholly below or wholly above the
                     ! Fermi level contributes no resolvent denominator to the
                     ! trace of the occupied projector; its basis-invariant
                     ! contribution is obtained by the contact term.  A
                     ! degeneracy crossing the zero-temperature occupation
                     ! edge is genuinely non-analytic and remains an error.
                     if ((eigenvalues(n) < fermi) .eqv. (endpoint_values(m) < fermi)) cycle
                     error stop 'force_theorem_finite_q_hessian_from_eigenbasis: Fermi-edge degeneracy'
                  end if
                  ! T_q(k) is stored in the production direction: rows at
                  ! k+q and columns at k.  The -q partner at k+q has rows at
                  ! k and columns at k+q.  Contract those directions directly
                  ! rather than silently taking their reverse matrix elements.
                  ta_nm = dot_product(endpoint_vectors(:,m), matmul(torques_q(:,:,a), eigenvectors(:,n)))
                  tb_mn = dot_product(eigenvectors(:,n), matmul(torques_minus_q(:,:,b), endpoint_vectors(:,m)))
                  tb_nm = dot_product(endpoint_vectors(:,m), matmul(torques_q(:,:,b), eigenvectors(:,n)))
                  ta_mn = dot_product(eigenvectors(:,n), matmul(torques_minus_q(:,:,a), endpoint_vectors(:,m)))
                  term = (ta_nm*tb_mn + tb_nm*ta_mn) / cmplx(eigenvalues(n)-endpoint_values(m),0.0_rp,rp)
                  torque_torque(a,b) = torque_torque(a,b) + term
               end do
            end do
         end do
      end do
      hessian = torque_torque + mixed_contact
      complete = hessian
   end subroutine force_theorem_finite_q_hessian_from_eigenbasis

   !> Brillouin-zone average of the finite-q spectral Hessian.  Endpoint
   !> arrays and every vertex/contact stack use the same k ordering.
   subroutine force_theorem_finite_q_hessian_from_eigenbasis_batch(eigenvalues, eigenvectors, endpoint_values, endpoint_vectors, &
      fermi, weights, torques_q, torques_minus_q, mixed, hessian, torque_torque, mixed_contact, complete)
      real(rp), intent(in) :: eigenvalues(:, :), endpoint_values(:, :), fermi, weights(:)
      complex(rp), intent(in) :: eigenvectors(:, :, :), endpoint_vectors(:, :, :)
      complex(rp), intent(in) :: torques_q(:, :, :, :), torques_minus_q(:, :, :, :), mixed(:, :, :, :, :)
      complex(rp), intent(out) :: hessian(:, :), torque_torque(:, :), mixed_contact(:, :), complete(:, :)
      complex(rp), allocatable :: h(:, :), tt(:, :), cc(:, :), allh(:, :)
      integer :: ik, nk, nsite
      real(rp) :: weight_sum

      nk = size(eigenvalues,2); nsite = size(torques_q,3)
      if (size(weights) /= nk .or. size(endpoint_values,2) /= nk .or. size(eigenvectors,3) /= nk .or. &
          size(endpoint_vectors,3) /= nk .or. size(torques_q,4) /= nk .or. size(torques_minus_q,4) /= nk .or. &
          size(mixed,5) /= nk .or. any(shape(hessian) /= [nsite,nsite]) .or. &
          any(shape(torque_torque) /= [nsite,nsite]) .or. any(shape(mixed_contact) /= [nsite,nsite]) .or. &
          any(shape(complete) /= [nsite,nsite])) error stop 'force_theorem_finite_q_hessian_from_eigenbasis_batch: shape mismatch'
      weight_sum = sum(weights)
      if (abs(weight_sum) <= tiny(1.0_rp)) error stop 'force_theorem_finite_q_hessian_from_eigenbasis_batch: zero weight sum'
      allocate(h(nsite,nsite),tt(nsite,nsite),cc(nsite,nsite),allh(nsite,nsite))
      hessian = 0.0_rp; torque_torque = 0.0_rp; mixed_contact = 0.0_rp; complete = 0.0_rp
      do ik = 1, nk
         call force_theorem_finite_q_hessian_from_eigenbasis(eigenvalues(:,ik), eigenvectors(:,:,ik), endpoint_values(:,ik), &
            endpoint_vectors(:,:,ik), fermi, torques_q(:,:,:,ik), torques_minus_q(:,:,:,ik), mixed(:,:,:,:,ik), h, tt, cc, allh)
         hessian = hessian + weights(ik)*h/weight_sum
         torque_torque = torque_torque + weights(ik)*tt/weight_sum
         mixed_contact = mixed_contact + weights(ik)*cc/weight_sum
         complete = complete + weights(ik)*allh/weight_sum
      end do
      deallocate(h,tt,cc,allh)
   end subroutine force_theorem_finite_q_hessian_from_eigenbasis_batch

   !> Finite-temperature finite-q grand-potential Hessian for one k -> k+q
   !> endpoint pair.  The TT term is written with the symmetric divided
   !> difference
   !>
   !>   K_nm = (f(e_n)-f(e_m))/(e_n-e_m),
   !>
   !> including m=n through f'(e_n).  The factor one half is essential when
   !> all ordered band pairs are retained.  It makes this form exactly equal
   !> to the occupied-state expression in a gapped zero-temperature limit,
   !> while avoiding cancellation between two nearly degenerate equally
   !> occupied (or equally empty) states in a metal.
   subroutine force_theorem_finite_q_hessian_from_eigenbasis_metallic(eigenvalues, eigenvectors, endpoint_values, endpoint_vectors, &
      fermi, kT, torques_q, torques_minus_q, mixed, hessian, torque_torque, mixed_contact, complete)
      real(rp), intent(in) :: eigenvalues(:), endpoint_values(:), fermi, kT
      complex(rp), intent(in) :: eigenvectors(:, :), endpoint_vectors(:, :)
      complex(rp), intent(in) :: torques_q(:, :, :), torques_minus_q(:, :, :), mixed(:, :, :, :)
      complex(rp), intent(out) :: hessian(:, :), torque_torque(:, :), mixed_contact(:, :), complete(:, :)
      complex(rp) :: ta_nm, tb_mn, tb_nm, ta_mn
      complex(rp), allocatable :: torque_q_band(:, :, :), torque_minus_band(:, :, :), mixed_band(:, :, :, :)
      real(rp) :: fn, kernel
      integer :: n, m, a, b, nbands, nsite, nmat

      nbands = size(eigenvalues); nmat = size(eigenvectors,1); nsite = size(torques_q,3)
      if (size(endpoint_values) /= nbands .or. size(eigenvectors,2) /= nbands .or. &
          any(shape(endpoint_vectors) /= [nmat,nbands]) .or. any(shape(torques_q) /= [nmat,nmat,nsite]) .or. &
          any(shape(torques_minus_q) /= [nmat,nmat,nsite]) .or. any(shape(mixed) /= [nmat,nmat,nsite,nsite]) .or. &
          any(shape(hessian) /= [nsite,nsite]) .or. any(shape(torque_torque) /= [nsite,nsite]) .or. &
          any(shape(mixed_contact) /= [nsite,nsite]) .or. any(shape(complete) /= [nsite,nsite])) then
         error stop 'force_theorem_finite_q_hessian_from_eigenbasis_metallic: shape mismatch'
      end if
      hessian = 0.0_rp; torque_torque = 0.0_rp; mixed_contact = 0.0_rp; complete = 0.0_rp
      ! Transform each vertex/contact once into the two endpoint eigenbases.
      ! The band-pair contraction below then contains only scalar products;
      ! this avoids repeating orbital-space matrix-vector products for every
      ! (n,m) pair.
      allocate(torque_q_band(nbands,nbands,nsite), torque_minus_band(nbands,nbands,nsite), &
         mixed_band(nbands,nbands,nsite,nsite))
      do a = 1, nsite
         torque_q_band(:,:,a) = matmul(conjg(transpose(endpoint_vectors)), &
            matmul(torques_q(:,:,a),eigenvectors))
         torque_minus_band(:,:,a) = matmul(conjg(transpose(eigenvectors)), &
            matmul(torques_minus_q(:,:,a),endpoint_vectors))
         do b = 1, nsite
            mixed_band(:,:,a,b) = matmul(conjg(transpose(eigenvectors)), &
               matmul(mixed(:,:,a,b),eigenvectors))
         end do
      end do
      do a = 1, nsite
         do b = 1, nsite
            do n = 1, nbands
               fn = finite_temperature_occupation(eigenvalues(n), fermi, kT)
               mixed_contact(a,b) = mixed_contact(a,b) + fn*mixed_band(n,n,a,b)
               do m = 1, nbands
                  kernel = fermi_divided_difference(eigenvalues(n), endpoint_values(m), fermi, kT)
                  ! T_q(k) has rows at k+q and columns at k.  Its -q
                  ! partner at k+q has the reverse endpoint ordering.
                  ta_nm = torque_q_band(m,n,a)
                  tb_mn = torque_minus_band(n,m,b)
                  tb_nm = torque_q_band(m,n,b)
                  ta_mn = torque_minus_band(n,m,a)
                  torque_torque(a,b) = torque_torque(a,b) + 0.5_rp*kernel*(ta_nm*tb_mn + tb_nm*ta_mn)
               end do
            end do
         end do
      end do
      hessian = torque_torque + mixed_contact
      complete = hessian
      deallocate(torque_q_band, torque_minus_band, mixed_band)
   end subroutine force_theorem_finite_q_hessian_from_eigenbasis_metallic

   !> Brillouin-zone average of the finite-temperature metallic Hessian.
   subroutine force_theorem_finite_q_hessian_from_eigenbasis_metallic_batch(eigenvalues, eigenvectors, endpoint_values, endpoint_vectors, &
      fermi, kT, weights, torques_q, torques_minus_q, mixed, hessian, torque_torque, mixed_contact, complete)
      real(rp), intent(in) :: eigenvalues(:, :), endpoint_values(:, :), fermi, kT, weights(:)
      complex(rp), intent(in) :: eigenvectors(:, :, :), endpoint_vectors(:, :, :)
      complex(rp), intent(in) :: torques_q(:, :, :, :), torques_minus_q(:, :, :, :), mixed(:, :, :, :, :)
      complex(rp), intent(out) :: hessian(:, :), torque_torque(:, :), mixed_contact(:, :), complete(:, :)
      complex(rp), allocatable :: h(:, :), tt(:, :), cc(:, :), allh(:, :)
      integer :: ik, nk, nsite
      real(rp) :: weight_sum

      nk = size(eigenvalues,2); nsite = size(torques_q,3)
      if (size(weights) /= nk .or. size(endpoint_values,2) /= nk .or. size(eigenvectors,3) /= nk .or. &
          size(endpoint_vectors,3) /= nk .or. size(torques_q,4) /= nk .or. size(torques_minus_q,4) /= nk .or. &
          size(mixed,5) /= nk .or. any(shape(hessian) /= [nsite,nsite]) .or. &
          any(shape(torque_torque) /= [nsite,nsite]) .or. any(shape(mixed_contact) /= [nsite,nsite]) .or. &
          any(shape(complete) /= [nsite,nsite])) error stop 'force_theorem_finite_q_hessian_from_eigenbasis_metallic_batch: shape mismatch'
      weight_sum = sum(weights)
      if (abs(weight_sum) <= tiny(1.0_rp)) error stop 'force_theorem_finite_q_hessian_from_eigenbasis_metallic_batch: zero weight sum'
      allocate(h(nsite,nsite),tt(nsite,nsite),cc(nsite,nsite),allh(nsite,nsite))
      hessian = 0.0_rp; torque_torque = 0.0_rp; mixed_contact = 0.0_rp; complete = 0.0_rp
      do ik = 1, nk
         call force_theorem_finite_q_hessian_from_eigenbasis_metallic(eigenvalues(:,ik), eigenvectors(:,:,ik), endpoint_values(:,ik), &
            endpoint_vectors(:,:,ik), fermi, kT, torques_q(:,:,:,ik), torques_minus_q(:,:,:,ik), mixed(:,:,:,:,ik), h, tt, cc, allh)
         hessian = hessian + weights(ik)*h/weight_sum
         torque_torque = torque_torque + weights(ik)*tt/weight_sum
         mixed_contact = mixed_contact + weights(ik)*cc/weight_sum
         complete = complete + weights(ik)*allh/weight_sum
      end do
      deallocate(h,tt,cc,allh)
   end subroutine force_theorem_finite_q_hessian_from_eigenbasis_metallic_batch

   !> Fermi-Dirac occupation with the same positive kT floor used by the
   !> reciprocal SCF occupation solver.
   pure real(rp) function finite_temperature_occupation(eigenvalue, fermi, kT) result(occupation)
      real(rp), intent(in) :: eigenvalue, fermi, kT
      real(rp) :: argument, effective_kT
      effective_kT = max(kT, 1.0e-10_rp)
      argument = (eigenvalue-fermi)/effective_kT
      if (argument >= 50.0_rp) then
         occupation = 0.0_rp
      else if (argument <= -50.0_rp) then
         occupation = 1.0_rp
      else
         occupation = 1.0_rp/(exp(argument)+1.0_rp)
      end if
   end function finite_temperature_occupation

   !> Stable divided difference of the Fermi occupation.  It is a response
   !> kernel, so the coincident-energy value is f'(e), not zero and not a
   !> dropped small denominator.  A midpoint Taylor value is used only when
   !> subtraction would lose floating-point digits; it is the smooth analytic
   !> continuation of the same kernel.
   pure real(rp) function fermi_divided_difference(e1, e2, fermi, kT) result(kernel)
      real(rp), intent(in) :: e1, e2, fermi, kT
      real(rp) :: effective_kT, delta, midpoint, fm, third_derivative, f1, f2
      effective_kT = max(kT, 1.0e-10_rp)
      delta = e1-e2
      midpoint = 0.5_rp*(e1+e2)
      fm = finite_temperature_occupation(midpoint, fermi, effective_kT)
      if (delta == 0.0_rp .or. abs(delta) <= sqrt(epsilon(1.0_rp))*max(1.0_rp,abs(e1),abs(e2),effective_kT)) then
         ! f'''(e) = -p(1-p)(1-6p+6p^2)/kT^3.
         third_derivative = -fm*(1.0_rp-fm)*(1.0_rp-6.0_rp*fm+6.0_rp*fm*fm)/(effective_kT**3)
         kernel = -fm*(1.0_rp-fm)/effective_kT + third_derivative*delta*delta/24.0_rp
      else
         f1 = finite_temperature_occupation(e1, fermi, effective_kT)
         f2 = finite_temperature_occupation(e2, fermi, effective_kT)
         kernel = (f1-f2)/delta
      end if
   end function fermi_divided_difference

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

   !> Compare the reconstructed fixture with the production reciprocal
   !> assembler real-space blocks at caller-selected k points.  This routine
   !> deliberately reads the production object only; it never writes ee, eeo,
   !> enim, or any native exchange state.
   subroutine lmto_fixture_adapter_residual(source, fixture, k_points, max_error)
      type(hamiltonian), intent(in) :: source
      type(lmto_live_hamiltonian_fixture), intent(in) :: fixture
      real(rp), intent(in) :: k_points(:, :)
      real(rp), intent(out) :: max_error
      complex(rp), allocatable :: expected(:, :), b(:, :), qmat(:, :), fixture_h(:, :), block(:, :), phase
      integer :: ik, ibond, site, m, nmat, ntype, ia, i_start, i_end, j_start, j_end

      call validate_fixture(fixture)
      if (size(k_points,1) /= 3 .or. .not. associated(source%lattice) .or. &
          .not. allocated(source%ee)) error stop 'lmto_fixture_adapter_residual: incomplete production state'
      if (source%ccor_2c) error stop 'lmto_fixture_adapter_residual: CCOR is outside the clean adapter gate'
      nmat = 2*fixture%norb*fixture%nsite
      if (size(source%ee,1) /= nmat .or. size(source%ee,2) /= nmat) then
         error stop 'lmto_fixture_adapter_residual: production basis mismatch'
      end if
      allocate(expected(nmat,nmat), b(nmat,nmat), qmat(nmat,nmat), fixture_h(nmat,nmat), block(nmat,nmat))
      max_error = 0.0_rp
      do ik = 1, size(k_points,2)
         expected = cmplx(0.0_rp,0.0_rp,rp)
         b = cmplx(0.0_rp,0.0_rp,rp); qmat = b
         ibond = 0
         do site = 1, fixture%nsite
            ntype = source%lattice%ib(site); ia = source%lattice%atlist(ntype)
            do m = 1, source%lattice%nn(ia,1)
               if (m > 1) then
                  if (source%lattice%nn(ia,m) < 1 .or. source%lattice%nn(ia,m) > source%lattice%kk) cycle
                  if (source%lattice%iz(source%lattice%nn(ia,m)) < 1 .or. &
                      source%lattice%iz(source%lattice%nn(ia,m)) > fixture%nsite) cycle
               end if
               ibond = ibond + 1
               block = source%ee(:,:,m,ntype)
               phase = bond_phase(fixture,k_points(:,ik),ibond)
               i_start=(site-1)*nmat/fixture%nsite+1; i_end=site*nmat/fixture%nsite
               j_start=(fixture%bond_target(ibond)-1)*nmat/fixture%nsite+1
               j_end=fixture%bond_target(ibond)*nmat/fixture%nsite
               b(i_start:i_end,j_start:j_end) = b(i_start:i_end,j_start:j_end) + phase*block
               if (fixture%hoh) qmat(i_start:i_end,j_start:j_end) = qmat(i_start:i_end,j_start:j_end) + &
                  phase*source%eeo(:,:,m,ntype)
            end do
         end do
         expected = b
         if (fixture%hoh) expected = expected - matmul(qmat,b)
         if (fixture%hoh .and. allocated(source%enim)) then
            do site = 1, fixture%nsite
               ntype = source%lattice%ib(site)
               i_start=(site-1)*nmat/fixture%nsite+1; i_end=site*nmat/fixture%nsite
               expected(i_start:i_end,i_start:i_end) = expected(i_start:i_end,i_start:i_end) + source%enim(:,:,ntype)
               if (allocated(source%lsham)) expected(i_start:i_end,i_start:i_end) = &
                  expected(i_start:i_end,i_start:i_end) + source%lsham(:,:,ntype)
            end do
         end if
         call assemble_lmto_hamiltonian(fixture,k_points(:,ik),fixture_h)
         max_error = max(max_error,maxval(abs(fixture_h-expected)))
      end do
      deallocate(expected,b,qmat,fixture_h,block)
   end subroutine lmto_fixture_adapter_residual

   subroutine assemble_finite_q_directional_terms(this, k_point, q_point, site, axis, db, dq)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: k_point(3), q_point(3), axis(3)
      integer, intent(in) :: site
      complex(rp), intent(out) :: db(:, :), dq(:, :)
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
            call build_onsite_derivative(this,this%obar1(:,target),this%moments(:,target),dm_target, dobar)
            call add_site_block(dq, phase_source*matmul(db_source,obar) + phase_target*(matmul(db_target,obar)+ &
               matmul(bond,dobar)), source,target,phase_k,this%norb)
         end if
      end do
   end subroutine assemble_finite_q_directional_terms

   subroutine assemble_finite_q_onsite_derivative(this, q_point, site, axis, matrix)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: q_point(3), axis(3)
      integer, intent(in) :: site
      complex(rp), intent(out) :: matrix(:, :)
      complex(rp) :: block(2*this%norb,2*this%norb), phase
      real(rp) :: dm(3)
      matrix=cmplx(0.0_rp,0.0_rp,rp); dm=cross3(axis,this%moments(:,site))
      call build_onsite_derivative(this,this%enu1(:,site),this%moments(:,site),dm,block)
      phase=site_phase(this,q_point,site)
      call add_site_block(matrix,phase*block,site,site,cmplx(1.0_rp,0.0_rp,rp),this%norb)
   end subroutine assemble_finite_q_onsite_derivative

   subroutine assemble_finite_q_mixed_base_terms(this, k_point, q_point, site_i, axis_i, site_j, axis_j, d2b, d2q, d2enu)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: k_point(3), q_point(3), axis_i(3), axis_j(3)
      integer, intent(in) :: site_i, site_j
      complex(rp), intent(out) :: d2b(:, :), d2q(:, :), d2enu(:, :)
      complex(rp) :: bond(2*this%norb,2*this%norb), dsi(2*this%norb,2*this%norb), dti(2*this%norb,2*this%norb), &
         dsj(2*this%norb,2*this%norb), dtj(2*this%norb,2*this%norb), d2bond(2*this%norb,2*this%norb), &
         obar(2*this%norb,2*this%norb), dobar_i(2*this%norb,2*this%norb), dobar_j(2*this%norb,2*this%norb), &
         d2obar(2*this%norb,2*this%norb), hmag(this%norb,this%norb,4), dhmag(this%norb,this%norb,4), &
         d2hmag(this%norb,this%norb,4), phase_k, pis, pit, pjs, pjt
      complex(rp) :: dbi(2*this%norb,2*this%norb), dbj(2*this%norb,2*this%norb), &
         d2phase(2*this%norb,2*this%norb)
      real(rp) :: dmi_s(3), dmi_t(3), dmj_s(3), dmj_t(3), d2mi(3), d2mj(3), zero(3)
      integer :: ibond, source, target

      d2b=cmplx(0.0_rp,0.0_rp,rp); d2q=d2b; d2enu=d2b; zero=0.0_rp
      d2mi=zero; d2mj=zero
      if (site_i == site_j) then
         d2mi=second_moment_variation(axis_i,this%moments(:,site_i))
         d2mj=second_moment_variation(axis_j,this%moments(:,site_j))
      end if
      do ibond=1,this%nbond
         source=this%bond_source(ibond); target=this%bond_target(ibond)
         dmi_s=zero; dmi_t=zero; dmj_s=zero; dmj_t=zero
         if (source == site_i) dmi_s=cross3(axis_i,this%moments(:,source))
         if (target == site_i) dmi_t=cross3(axis_i,this%moments(:,target))
         if (source == site_j) dmj_s=cross3(axis_j,this%moments(:,source))
         if (target == site_j) dmj_t=cross3(axis_j,this%moments(:,target))
         if (maxval(abs(dmi_s))+maxval(abs(dmi_t))+maxval(abs(dmj_s))+maxval(abs(dmj_t)) == 0.0_rp .and. &
             site_i /= site_j) cycle
         call build_bond(this,ibond,hmag); call fixture_hhmag_to_spinor(this,hmag,bond)
         call build_bond_derivative(this,ibond,dmi_s,zero,dhmag); call fixture_hhmag_to_spinor(this,dhmag,dsi)
         call build_bond_derivative(this,ibond,zero,dmi_t,dhmag); call fixture_hhmag_to_spinor(this,dhmag,dti)
         call build_bond_derivative(this,ibond,dmj_s,zero,dhmag); call fixture_hhmag_to_spinor(this,dhmag,dsj)
         call build_bond_derivative(this,ibond,zero,dmj_t,dhmag); call fixture_hhmag_to_spinor(this,dhmag,dtj)
         d2bond=cmplx(0.0_rp,0.0_rp,rp)
         call add_mixed_bond_component(this,ibond,dmi_s,zero,dmj_s,zero,d2bond,p_endpoint(this,q_point,ibond,.false.)* &
            p_endpoint(this,-q_point,ibond,.false.))
         call add_mixed_bond_component(this,ibond,dmi_s,zero,zero,dmj_t,d2bond,p_endpoint(this,q_point,ibond,.false.)* &
            p_endpoint(this,-q_point,ibond,.true.))
         call add_mixed_bond_component(this,ibond,zero,dmi_t,dmj_s,zero,d2bond,p_endpoint(this,q_point,ibond,.true.)* &
            p_endpoint(this,-q_point,ibond,.false.))
         call add_mixed_bond_component(this,ibond,zero,dmi_t,zero,dmj_t,d2bond,p_endpoint(this,q_point,ibond,.true.)* &
            p_endpoint(this,-q_point,ibond,.true.))
         if (site_i == site_j) then
            if (source == site_i) then
               call build_bond_derivative(this,ibond,d2mi,zero,dhmag); call fixture_hhmag_to_spinor(this,dhmag,d2phase)
               d2bond=d2bond+d2phase
            end if
            if (target == site_i) then
               call build_bond_derivative(this,ibond,zero,d2mi,dhmag); call fixture_hhmag_to_spinor(this,dhmag,d2phase)
               d2bond=d2bond+d2phase
            end if
         end if
         phase_k=bond_phase(this,k_point,ibond); call add_site_block(d2b,d2bond,source,target,phase_k,this%norb)
         if (this%hoh) then
            call build_onsite_coefficient(this,this%obar0(:,target),this%obar1(:,target),this%moments(:,target),obar)
            dbi=cmplx(0.0_rp,0.0_rp,rp); dbj=dbi
            pis=endpoint_phase(this,q_point,ibond,.false.); pit=endpoint_phase(this,q_point,ibond,.true.)
            pjs=endpoint_phase(this,-q_point,ibond,.false.); pjt=endpoint_phase(this,-q_point,ibond,.true.)
            dbi=pis*dsi+pit*dti; dbj=pjs*dsj+pjt*dtj
            call build_onsite_derivative(this,this%obar1(:,target),this%moments(:,target),dmj_t,dobar_j)
            call build_onsite_derivative(this,this%obar1(:,target),this%moments(:,target),dmi_t,dobar_i)
            d2obar=cmplx(0.0_rp,0.0_rp,rp)
            if (site_i == site_j .and. target == site_i) then
               call build_onsite_derivative(this,this%obar1(:,target),this%moments(:,target),d2mi,d2obar)
            end if
            ! The first term is d2B*O; the remaining terms are the
            ! endpoint/product rules.
            call add_site_block(d2q, matmul(d2bond,obar) + matmul(dsi*pis+dti*pit, pjt*dobar_j) + &
               matmul(dsj*pjs+dtj*pjt,pit*dobar_i) + matmul(bond,d2obar), source,target,phase_k,this%norb)
         end if
      end do
      if (this%include_enu .and. site_i == site_j) then
         call build_onsite_derivative(this,this%enu1(:,site_i),this%moments(:,site_i),d2mi,d2phase)
         call add_site_block(d2enu,d2phase,site_i,site_i,cmplx(1.0_rp,0.0_rp,rp),this%norb)
      end if
   end subroutine assemble_finite_q_mixed_base_terms

   subroutine add_mixed_bond_component(this, ibond, dm_i_s, dm_i_t, dm_j_s, dm_j_t, accum, factor)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      integer, intent(in) :: ibond
      real(rp), intent(in) :: dm_i_s(3), dm_i_t(3), dm_j_s(3), dm_j_t(3)
      complex(rp), intent(inout) :: accum(:, :)
      complex(rp), intent(in) :: factor
      complex(rp) :: hmag(this%norb,this%norb,4), block(2*this%norb,2*this%norb)
      real(rp) :: zero(3)
      zero=0.0_rp
      if (maxval(abs(dm_i_s))+maxval(abs(dm_i_t))+maxval(abs(dm_j_s))+maxval(abs(dm_j_t)) == 0.0_rp) return
      call lmto_bond_mixed_derivative(this%hhh(:,:,ibond),this%wx0(:,this%bond_source(ibond)),this%wx1(:,this%bond_source(ibond)), &
         this%wx0(:,this%bond_target(ibond)),this%wx1(:,this%bond_target(ibond)),this%moments(:,this%bond_source(ibond)), &
         this%moments(:,this%bond_target(ibond)),dm_i_s,dm_i_t,dm_j_s,dm_j_t,hmag)
      call fixture_hhmag_to_spinor(this,hmag,block); accum=accum+factor*block
   end subroutine add_mixed_bond_component

   pure function endpoint_phase(this,q_point,ibond,target_endpoint) result(phase)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: q_point(3)
      integer, intent(in) :: ibond
      logical, intent(in) :: target_endpoint
      complex(rp) :: phase
      real(rp) :: position(3)
      position=0.0_rp
      if (target_endpoint) position=this%bond_vector(:,ibond)
      phase=cmplx(cos(2.0_rp*pi*dot_product(q_point,position)), &
         sin(2.0_rp*pi*dot_product(q_point,position)),rp)
   end function endpoint_phase

   pure function p_endpoint(this,q_point,ibond,target_endpoint) result(phase)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: q_point(3)
      integer, intent(in) :: ibond
      logical, intent(in) :: target_endpoint
      complex(rp) :: phase
      phase=endpoint_phase(this,q_point,ibond,target_endpoint)
   end function p_endpoint

   pure function site_phase(this,q_point,site) result(phase)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: q_point(3)
      integer, intent(in) :: site
      complex(rp) :: phase
      ! The onsite endpoint phase cancels against the absolute basis phase in
      ! <a,k+q|delta H|a,k>; retain the site argument for a single convention
      ! point and to make that cancellation explicit.
      phase=cmplx(1.0_rp,0.0_rp,rp)
   end function site_phase

   pure function second_moment_variation(axis,moment) result(dm2)
      real(rp), intent(in) :: axis(3), moment(3)
      real(rp) :: dm2(3)
      dm2=cross3(axis,cross3(axis,moment))
   end function second_moment_variation

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
         call fixture_hhmag_to_spinor(this,hmag,bond)
         phase = cmplx(cos(2.0_rp*pi*dot_product(k_point,this%bond_vector(:,ibond))), &
                       sin(2.0_rp*pi*dot_product(k_point,this%bond_vector(:,ibond))), rp)
         call add_site_block(b, bond, source, target, phase, this%norb)
         if (this%hoh) then
            call build_onsite_coefficient(this,this%obar0(:,target), this%obar1(:,target), this%moments(:,target), obar)
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
         call fixture_hhmag_to_spinor(this,hmag,bond)
         call build_bond_derivative(this, ibond, dm_source, dm_target, dhmag)
         call fixture_hhmag_to_spinor(this,dhmag,dbond)
         phase = bond_phase(this, k_point, ibond)
         call add_site_block(db, dbond, source, target, phase, this%norb)
         if (this%hoh) then
            call build_onsite_coefficient(this,this%obar0(:,target), this%obar1(:,target), this%moments(:,target), obar)
            call build_onsite_derivative(this,this%obar1(:,target), this%moments(:,target), dm_target, dobar)
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
         call fixture_hhmag_to_spinor(this,hmag,bond)
         call build_bond_derivative(this, ibond, dm_si, dm_ti, dhmag)
         call fixture_hhmag_to_spinor(this,dhmag,dbond_i)
         call build_bond_derivative(this, ibond, dm_sj, dm_tj, dhmag)
         call fixture_hhmag_to_spinor(this,dhmag,dbond_j)
         call lmto_bond_mixed_derivative(this%hhh(:,:,ibond), this%wx0(:,source), this%wx1(:,source), &
            this%wx0(:,target), this%wx1(:,target), this%moments(:,source), this%moments(:,target), &
            dm_si, dm_ti, dm_sj, dm_tj, d2hmag)
         call fixture_hhmag_to_spinor(this,d2hmag,d2bond)
         phase = bond_phase(this, k_point, ibond)
         call add_site_block(d2b, d2bond, source, target, phase, this%norb)
         if (this%hoh) then
            call build_onsite_coefficient(this,this%obar0(:,target), this%obar1(:,target), this%moments(:,target), obar)
            call build_onsite_derivative(this,this%obar1(:,target), this%moments(:,target), dm_ti, dobar_i)
            call build_onsite_derivative(this,this%obar1(:,target), this%moments(:,target), dm_tj, dobar_j)
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

   subroutine fixture_hhmag_to_spinor(this, hhmag, spinor)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      complex(rp), intent(in) :: hhmag(:, :, :)
      complex(rp), intent(out) :: spinor(:, :)
      complex(rp), allocatable :: work(:, :, :)
      integer :: idir

      allocate(work, source=hhmag)
      if (this%cartesian_to_spherical) then
         do idir = 1, size(work, 3)
            call hcpx(work(:, :, idir), 'cart2sph')
         end do
      end if
      call lmto_hhmag_to_spinor(work, spinor)
      deallocate(work)
   end subroutine fixture_hhmag_to_spinor

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

   subroutine build_onsite_coefficient(this, c0, c1, moment, spinor)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      complex(rp), intent(in) :: c0(:), c1(:)
      real(rp), intent(in) :: moment(3)
      complex(rp), intent(out) :: spinor(:, :)
      complex(rp) :: zero_hhh(size(c0),size(c0)), zero(size(c0)), hmag(size(c0),size(c0),4)
      zero_hhh = cmplx(0.0_rp, 0.0_rp, rp); zero = cmplx(0.0_rp, 0.0_rp, rp)
      call lmto_bond_value(zero_hhh, zero, zero, zero, zero, c0, c1, moment, moment, .true., hmag)
      call fixture_hhmag_to_spinor(this,hmag,spinor)
   end subroutine build_onsite_coefficient

   subroutine build_onsite_derivative(this, c1, moment, dm, spinor)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      complex(rp), intent(in) :: c1(:)
      real(rp), intent(in) :: moment(3), dm(3)
      complex(rp), intent(out) :: spinor(:, :)
      complex(rp) :: zero_hhh(size(c1),size(c1)), zero(size(c1)), dhmag(size(c1),size(c1),4)
      zero_hhh = cmplx(0.0_rp, 0.0_rp, rp); zero = cmplx(0.0_rp, 0.0_rp, rp)
      call lmto_bond_derivative(zero_hhh, zero, zero, zero, zero, c1, moment, moment, dm, dm, .true., dhmag)
      call fixture_hhmag_to_spinor(this,dhmag,spinor)
   end subroutine build_onsite_derivative

   subroutine assemble_onsite_coefficient(this, c0, c1, matrix)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      complex(rp), intent(in) :: c0(:, :), c1(:, :)
      complex(rp), intent(out) :: matrix(:, :)
      complex(rp) :: block(2*this%norb,2*this%norb)
      integer :: site
      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, this%nsite
         call build_onsite_coefficient(this,c0(:,site), c1(:,site), this%moments(:,site), block)
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
      call build_onsite_derivative(this,this%enu1(:,site), this%moments(:,site), dm, block)
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
      if (size(a,2) /= n .or. size(b,1) /= n .or. size(b,2) /= n) then
         error stop 'trace_product: invalid two-factor matrix shape'
      end if
      if (present(c)) then
         if (.not. present(d)) error stop 'trace_product: d is required when c is present'
         if (size(c,1) /= n .or. size(c,2) /= n .or. size(d,1) /= n .or. size(d,2) /= n) then
            error stop 'trace_product: invalid four-factor matrix shape'
         end if
         product = matmul(a,matmul(b,matmul(c,d)))
      else
         if (present(d)) error stop 'trace_product: c is required when d is present'
         product = matmul(a,b)
      end if
      value = sum([(product(i,i), i=1,n)])
   end function trace_product

   pure function cross3(a, b) result(c)
      real(rp), intent(in) :: a(3), b(3)
      real(rp) :: c(3)
      c = [a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), a(1)*b(2)-a(2)*b(1)]
   end function cross3

end module lr_kl_hessian_mod
