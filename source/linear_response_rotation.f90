!------------------------------------------------------------------------------
! Linear-response rotation, finite-H Hessian, and finite-H contour services.
!------------------------------------------------------------------------------
submodule (linear_response_mod) linear_response_rotation
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use precision_mod, only: rp
   use math_mod, only: i_unit, pi, inverse_3x3, hcpx, ang2au
   use control_mod, only: control
   use energy_mod, only: energy
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use lmto_magnetic_tangent_mod, only: lmto_bond_value, lmto_bond_derivative, &
                                        lmto_bond_mixed_derivative, lmto_hhmag_to_spinor
   use logger_mod, only: g_logger
   use lmto_path_operator_mod, only: native_turek_contour_options, native_turek_contour_report, &
      native_turek_static_reference
   use reciprocal_mod, only: reciprocal
   use self_mod, only: self
   use string_mod, only: int2str, lower
   implicit none

   real(rp), parameter :: kB_ry_per_k = 6.3336814e-6_rp
   real(rp), parameter :: q_tolerance = 1.0e-12_rp

   interface
      subroutine zgetrf(m, n, a, lda, ipiv, info)
         import :: rp
         integer, intent(in) :: m, n, lda
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
         complex(rp), intent(inout) :: a(lda, *)
      end subroutine zgetrf
      subroutine zgetri(n, a, lda, ipiv, work, lwork, info)
         import :: rp
         integer, intent(in) :: n, lda, lwork
         integer, intent(in) :: ipiv(*)
         integer, intent(out) :: info
         complex(rp), intent(inout) :: a(lda, *), work(*)
      end subroutine zgetri
      subroutine zgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: rp
         character(len=1), intent(in) :: trans
         integer, intent(in) :: n, nrhs, lda, ldb
         complex(rp), intent(in) :: a(lda, *)
         integer, intent(in) :: ipiv(*)
         complex(rp), intent(inout) :: b(ldb, *)
         integer, intent(out) :: info
      end subroutine zgetrs
   end interface

contains

   ! --- from lr_kl_hessian ---

   module subroutine lmto_fixture_init(this, nsite, norb, nbond, hoh)
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
   module subroutine lmto_fixture_from_hamiltonian(source, fixture)
      type(hamiltonian), intent(in) :: source
      type(lmto_live_hamiltonian_fixture), intent(out) :: fixture
      integer :: site, ntype, ia, nr, m, ja, target, it, count, ibond
      real(rp), allocatable :: cralat(:, :), ham_vec(:, :)
      integer :: nn_max_loc

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

   module subroutine lmto_fixture_clear(this)
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
   module subroutine assemble_lmto_hamiltonian(this, k_point, hamiltonian)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: k_point(3)
      complex(rp), intent(out) :: hamiltonian(:, :)
      complex(rp), allocatable :: b(:, :), q(:, :), enu(:, :)
      integer :: nmat

      call validate_fixture(this)
      nmat = 2*this%norb*this%nsite
      allocate(b(nmat,nmat), q(nmat,nmat), enu(nmat,nmat))
      call assemble_base_terms(this, k_point, b, q)
      enu = cmplx(0.0_rp, 0.0_rp, rp)
      if (this%include_enu) call assemble_onsite_coefficient(this, this%enu0, this%enu1, enu)
      hamiltonian = b + enu
      if (this%hoh) hamiltonian = hamiltonian - matmul(q, b)
      deallocate(b, q, enu)
   end subroutine assemble_lmto_hamiltonian

   !> First finite-q vertex for the convention
   !> theta_(aR)=theta_a(q) exp(+i 2*pi*q.(R+tau_a)).  After factoring the
   !> production Bloch basis, the returned matrix is the k -> k+q block and
   !> the target endpoint carries the relative phase exp(+i q.D), with
   !> D=R+tau_target-tau_source.  Its second argument is deliberately not folded:
   !> reciprocal-vector changes are retained as the corresponding basis gauge.
   module subroutine assemble_lmto_finite_q_torque(this, k_point, q_point, site, axis, torque)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: k_point(3), q_point(3), axis(3)
      integer, intent(in) :: site
      complex(rp), intent(out) :: torque(:, :)
      complex(rp), allocatable :: b(:, :), b_right(:, :), qmat(:, :), db(:, :), dq(:, :), denu(:, :)
      integer :: nmat

      call validate_site_axis(this, site, axis)
      nmat = 2*this%norb*this%nsite
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
   module subroutine assemble_lmto_finite_q_torques(this, k_point, q_point, axes, torques)
      type(lmto_live_hamiltonian_fixture), intent(in) :: this
      real(rp), intent(in) :: k_point(3), q_point(3), axes(:, :)
      complex(rp), intent(out) :: torques(:, :, :)
      integer :: site, nmat
      nmat = 2*this%norb*this%nsite
      do site = 1, this%nsite
         call assemble_lmto_finite_q_torque(this, k_point, q_point, site, axes(:,site), torques(:,:,site))
      end do
   end subroutine assemble_lmto_finite_q_torques

   !> Complete C_ab(k;q,-q), including same-site rotation curvature and the
   !> two distinct momentum orderings in the HOH product rule.
   module subroutine assemble_lmto_finite_q_mixed_derivative(this, k_point, q_point, site_i, axis_i, site_j, axis_j, mixed)
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

   !> Finite-q spectral Hessian for one k -> k+q endpoint pair.  The first
   !> vertex stack is evaluated at k for +q and the second at k+q for -q.
   !> The second ordering is retained explicitly; it is needed for a Hermitian
   !> sublattice Hessian when a and b are different.
   module subroutine force_theorem_finite_q_hessian_from_eigenbasis(eigenvalues, eigenvectors, endpoint_values, endpoint_vectors, &
      fermi, torques_q, torques_minus_q, mixed, hessian, torque_torque, mixed_contact, complete)
      real(rp), intent(in) :: eigenvalues(:), endpoint_values(:), fermi
      complex(rp), intent(in) :: eigenvectors(:, :), endpoint_vectors(:, :)
      complex(rp), intent(in) :: torques_q(:, :, :), torques_minus_q(:, :, :), mixed(:, :, :, :)
      complex(rp), intent(out) :: hessian(:, :), torque_torque(:, :), mixed_contact(:, :), complete(:, :)
      complex(rp) :: ta_nm, tb_mn, tb_nm, ta_mn, term
      integer :: n, m, a, b, nbands, nsite, nmat

      nbands = size(eigenvalues); nmat = size(eigenvectors,1); nsite = size(torques_q,3)
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
   module subroutine force_theorem_finite_q_hessian_from_eigenbasis_batch(eigenvalues, eigenvectors, endpoint_values, endpoint_vectors, &
      fermi, weights, torques_q, torques_minus_q, mixed, hessian, torque_torque, mixed_contact, complete)
      real(rp), intent(in) :: eigenvalues(:, :), endpoint_values(:, :), fermi, weights(:)
      complex(rp), intent(in) :: eigenvectors(:, :, :), endpoint_vectors(:, :, :)
      complex(rp), intent(in) :: torques_q(:, :, :, :), torques_minus_q(:, :, :, :), mixed(:, :, :, :, :)
      complex(rp), intent(out) :: hessian(:, :), torque_torque(:, :), mixed_contact(:, :), complete(:, :)
      complex(rp), allocatable :: h(:, :), tt(:, :), cc(:, :), allh(:, :)
      integer :: ik, nk, nsite
      real(rp) :: weight_sum

      nk = size(eigenvalues,2); nsite = size(torques_q,3)
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
   module subroutine force_theorem_finite_q_hessian_from_eigenbasis_metallic(eigenvalues, eigenvectors, endpoint_values, endpoint_vectors, &
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
   module subroutine force_theorem_finite_q_hessian_from_eigenbasis_metallic_batch(eigenvalues, eigenvectors, endpoint_values, endpoint_vectors, &
      fermi, kT, weights, torques_q, torques_minus_q, mixed, hessian, torque_torque, mixed_contact, complete)
      real(rp), intent(in) :: eigenvalues(:, :), endpoint_values(:, :), fermi, kT, weights(:)
      complex(rp), intent(in) :: eigenvectors(:, :, :), endpoint_vectors(:, :, :)
      complex(rp), intent(in) :: torques_q(:, :, :, :), torques_minus_q(:, :, :, :), mixed(:, :, :, :, :)
      complex(rp), intent(out) :: hessian(:, :), torque_torque(:, :), mixed_contact(:, :), complete(:, :)
      complex(rp), allocatable :: h(:, :), tt(:, :), cc(:, :), allh(:, :)
      integer :: ik, nk, nsite
      real(rp) :: weight_sum

      nk = size(eigenvalues,2); nsite = size(torques_q,3)
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
   module pure real(rp) function finite_temperature_occupation(eigenvalue, fermi, kT) result(occupation)
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
   module pure real(rp) function fermi_divided_difference(e1, e2, fermi, kT) result(kernel)
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

   !> Compare the reconstructed fixture with the production reciprocal
   !> assembler real-space blocks at caller-selected k points.  This routine
   !> deliberately reads the production object only; it never writes ee, eeo,
   !> enim, or any native exchange state.
   module subroutine lmto_fixture_adapter_residual(source, fixture, k_points, max_error)
      type(hamiltonian), intent(in) :: source
      type(lmto_live_hamiltonian_fixture), intent(in) :: fixture
      real(rp), intent(in) :: k_points(:, :)
      real(rp), intent(out) :: max_error
      complex(rp), allocatable :: expected(:, :), b(:, :), qmat(:, :), fixture_h(:, :), block(:, :), phase
      integer :: ik, ibond, site, m, nmat, ntype, ia, i_start, i_end, j_start, j_end

      call validate_fixture(fixture)
      if (source%ccor_2c) error stop 'lmto_fixture_adapter_residual: CCOR is outside the clean adapter gate'
      nmat = 2*fixture%norb*fixture%nsite
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
      real(rp) :: dmi_s(3), dmi_t(3), dmj_s(3), dmj_t(3), d2mi(3), zero(3)
      integer :: ibond, source, target

      d2b=cmplx(0.0_rp,0.0_rp,rp); d2q=d2b; d2enu=d2b; zero=0.0_rp
      d2mi=zero
      if (site_i == site_j) then
         ! The same-site contact is a mixed derivative in two independent
         ! rotation coordinates.  Use the symmetric Hessian of the unit-vector
         ! rotation map; for equal axes this reduces to the established second
         ! directional variation.
         d2mi=mixed_moment_variation(axis_i,axis_j,this%moments(:,site_i))
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

   pure function mixed_moment_variation(axis_a,axis_b,moment) result(dm2)
      real(rp), intent(in) :: axis_a(3),axis_b(3),moment(3)
      real(rp) :: dm2(3)
      dm2=0.5_rp*(cross3(axis_a,cross3(axis_b,moment))+cross3(axis_b,cross3(axis_a,moment)))
   end function mixed_moment_variation

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

   pure function cross3(a, b) result(c)
      real(rp), intent(in) :: a(3), b(3)
      real(rp) :: c(3)
      c = [a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), a(1)*b(2)-a(2)*b(1)]
   end function cross3


   ! --- from lr_kl_contour ---

   !> Construct the counter-clockwise ellipse used by the finite-T bridge.
   !>
   !> The ellipse encloses the real-axis spectra of both supplied matrices.
   !> Its semiminor axis is intentionally allowed to cross Matsubara poles.
   !> Those poles are returned explicitly so that their residues can be added
   !> back.  This permits a well-conditioned, broad contour at low temperature
   !> without silently changing the Fermi operator.
   module subroutine build_finite_temperature_contour(h_source, h_endpoint, fermi, kT, options, nodes, weights, fermi_poles)
      complex(rp), intent(in) :: h_source(:, :), h_endpoint(:, :)
      real(rp), intent(in) :: fermi, kT
      type(finite_h_contour_options), intent(in) :: options
      complex(rp), allocatable, intent(out) :: nodes(:), weights(:), fermi_poles(:)
      real(rp) :: lower, upper, center, semimajor, semiminor, theta, dtheta
      real(rp) :: desired_height, pole_y, ellipse_value
      integer :: i, l, max_l, npoint, npole

      if (kT <= 0.0_rp) error stop 'build_finite_temperature_contour: kT must be positive'
      call validate_options(options)
      call hermitian_bounds(h_source, h_endpoint, lower, upper)
      center = 0.5_rp*(lower + upper)
      semimajor = max(0.5_rp*(upper-lower) + options%contour_margin, options%contour_margin)
      desired_height = max(options%contour_height_fraction*semimajor, 1.0e-8_rp)
      if (.not. options%account_fermi_poles) desired_height = min(desired_height, 0.45_rp*pi*kT)
      semiminor = desired_height

      npoint = options%contour_points
      allocate(nodes(npoint), weights(npoint))
      dtheta = 2.0_rp*pi/real(npoint,rp)
      do i = 1, npoint
         theta = dtheta*real(i-1,rp)
         nodes(i) = cmplx(center + semimajor*cos(theta), semiminor*sin(theta), rp)
         ! z(theta) is counter-clockwise for theta increasing from zero.
         weights(i) = ((-semimajor*sin(theta) + i_unit*semiminor*cos(theta))*dtheta)/(2.0_rp*pi*i_unit)
      end do

      if (options%account_fermi_poles) then
         max_l = max(2, int(ceiling(semiminor/(pi*kT))) + 2)
         npole = 0
         do l = -max_l, max_l
            pole_y = pi*kT*real(2*l+1,rp)
            ellipse_value = ((fermi-center)/semimajor)**2 + (pole_y/semiminor)**2
            if (ellipse_value < 1.0_rp-1.0e-12_rp) npole = npole + 1
         end do
         allocate(fermi_poles(npole))
         npole = 0
         do l = -max_l, max_l
            pole_y = pi*kT*real(2*l+1,rp)
            ellipse_value = ((fermi-center)/semimajor)**2 + (pole_y/semiminor)**2
            if (ellipse_value < 1.0_rp-1.0e-12_rp) then
               npole = npole + 1
               fermi_poles(npole) = cmplx(fermi,pole_y,rp)
            end if
         end do
      else
         allocate(fermi_poles(0))
      end if
   end subroutine build_finite_temperature_contour

   !> Construct an occupied-state contour for the T=0 limit.
   !>
   !> `occupied_bounds=[emin,emax]` must lie strictly inside the insulating
   !> gap, with emax below the first unoccupied eigenvalue.  The weight is
   !> f(z)=1, and no finite-T Fermi poles are present.  The routine does not
   !> diagonalize either matrix; the bounds are an explicit fixture/driver
   !> contract for this classical limit.
   module subroutine build_zero_temperature_occupied_contour(occupied_bounds, options, nodes, weights)
      real(rp), intent(in) :: occupied_bounds(2)
      type(finite_h_contour_options), intent(in) :: options
      complex(rp), allocatable, intent(out) :: nodes(:), weights(:)
      real(rp) :: center, semimajor, semiminor, theta, dtheta
      integer :: i, npoint

      call validate_options(options)
      if (occupied_bounds(2) <= occupied_bounds(1)) then
         error stop 'build_zero_temperature_occupied_contour: invalid occupied bounds'
      end if
      center = 0.5_rp*sum(occupied_bounds)
      semimajor = 0.5_rp*(occupied_bounds(2)-occupied_bounds(1)) + options%contour_margin
      semiminor = max(options%contour_height_fraction*semimajor, 1.0e-8_rp)
      npoint = options%contour_points
      allocate(nodes(npoint), weights(npoint))
      dtheta = 2.0_rp*pi/real(npoint,rp)
      do i = 1, npoint
         theta = dtheta*real(i-1,rp)
         nodes(i) = cmplx(center + semimajor*cos(theta), semiminor*sin(theta), rp)
         weights(i) = ((-semimajor*sin(theta) + i_unit*semiminor*cos(theta))*dtheta)/(2.0_rp*pi*i_unit)
      end do
   end subroutine build_zero_temperature_occupied_contour

   !> Stable complex Fermi function on a contour node.
   module pure complex(rp) function finite_temperature_complex_fermi(z, fermi, kT) result(value)
      complex(rp), intent(in) :: z
      real(rp), intent(in) :: fermi, kT
      complex(rp) :: argument, reduced

      if (kT <= 0.0_rp) error stop 'finite_temperature_complex_fermi: kT must be positive'
      argument = (z-fermi)/kT
      if (real(argument,rp) > 40.0_rp) then
         reduced = exp(-argument)
         value = reduced/(1.0_rp+reduced)
      else if (real(argument,rp) < -40.0_rp) then
         reduced = exp(argument)
         value = 1.0_rp/(1.0_rp+reduced)
      else
         value = 1.0_rp/(exp(argument)+1.0_rp)
      end if
   end function finite_temperature_complex_fermi

   !> Fermi function with all poles enclosed by the selected contour removed.
   !>
   !> The residue of f at z_l=mu+i*pi*(2l+1)kT is -kT.  Adding kT/(z-z_l)
   !> cancels that pole.  Cauchy theorem then gives
   !>
   !>   integral_C f_reg(z) R(z) dz/(2*pi*i)
   !>     = integral_C f(z) R(z) dz/(2*pi*i) + kT sum_l R(z_l),
   !>
   !> which is precisely the physical-spectrum contour after the Fermi-pole
   !> residues have been removed.  This regularization is numerically much
   !> better than integrating a meromorphic weight close to its last enclosed
   !> pole, while remaining an exact finite-dimensional identity.
   module pure complex(rp) function finite_temperature_regularized_fermi(z, fermi, kT, fermi_poles) result(value)
      complex(rp), intent(in) :: z, fermi_poles(:)
      real(rp), intent(in) :: fermi, kT
      integer :: i
      value = finite_temperature_complex_fermi(z,fermi,kT)
      do i = 1, size(fermi_poles)
         value = value + kT/(z-fermi_poles(i))
      end do
   end function finite_temperature_regularized_fermi

   !> Finite-T direct-resolvent Hessian for one k -> k+q pair.
   !>
   !> H_source is H(k), H_endpoint is H(k+q).  T_q has rows at k+q and
   !> columns at k; T_minus_q has rows at k and columns at k+q.  The returned
   !> TT term is the explicit ordered-pair average, while contact is
   !> Tr[f(H_k) C(q,-q;k)].
   module subroutine force_theorem_finite_q_hessian_from_resolvent(h_source, h_endpoint, fermi, kT, torques_q, torques_minus_q, mixed, &
      options, hessian, torque_torque, mixed_contact, complete, report)
      complex(rp), intent(in) :: h_source(:, :), h_endpoint(:, :)
      real(rp), intent(in) :: fermi, kT
      complex(rp), intent(in) :: torques_q(:, :, :), torques_minus_q(:, :, :), mixed(:, :, :, :)
      type(finite_h_contour_options), intent(in) :: options
      complex(rp), intent(out) :: hessian(:, :), torque_torque(:, :), mixed_contact(:, :), complete(:, :)
      type(finite_h_contour_report), intent(out), optional :: report
      complex(rp), allocatable :: nodes(:), weights(:), poles(:)

      call build_finite_temperature_contour(h_source, h_endpoint, fermi, kT, options, nodes, weights, poles)
      call resolvent_hessian_core(h_source, h_endpoint, fermi, kT, .true., nodes, weights, poles, torques_q, torques_minus_q, mixed, &
         hessian, torque_torque, mixed_contact, complete, report)
      deallocate(nodes, weights, poles)
   end subroutine force_theorem_finite_q_hessian_from_resolvent

   !> Brillouin-zone average of the direct-resolvent Hessian.
   module subroutine force_theorem_finite_q_hessian_from_resolvent_batch(h_source, h_endpoint, fermi, kT, k_weights, torques_q, &
      torques_minus_q, mixed, options, hessian, torque_torque, mixed_contact, complete, report)
      complex(rp), intent(in) :: h_source(:, :, :), h_endpoint(:, :, :)
      real(rp), intent(in) :: fermi, kT, k_weights(:)
      complex(rp), intent(in) :: torques_q(:, :, :, :), torques_minus_q(:, :, :, :), mixed(:, :, :, :, :)
      type(finite_h_contour_options), intent(in) :: options
      complex(rp), intent(out) :: hessian(:, :), torque_torque(:, :), mixed_contact(:, :), complete(:, :)
      type(finite_h_contour_report), intent(out), optional :: report
      complex(rp), allocatable :: h(:, :), tt(:, :), cc(:, :), allh(:, :)
      type(finite_h_contour_report) :: one_report
      real(rp) :: weight_sum
      integer :: ik, nk, nsite

      nk = size(h_source,3); nsite = size(torques_q,3)
      weight_sum = sum(k_weights)
      if (weight_sum <= tiny(1.0_rp)) error stop 'resolvent batch: zero k-weight sum'
      allocate(h(nsite,nsite), tt(nsite,nsite), cc(nsite,nsite), allh(nsite,nsite))
      hessian = 0.0_rp; torque_torque = 0.0_rp; mixed_contact = 0.0_rp; complete = 0.0_rp
      if (present(report)) report = finite_h_contour_report()
      do ik = 1, nk
         call force_theorem_finite_q_hessian_from_resolvent(h_source(:,:,ik), h_endpoint(:,:,ik), fermi, kT, torques_q(:,:,:,ik), &
            torques_minus_q(:,:,:,ik), mixed(:,:,:,:,ik), options, h, tt, cc, allh, one_report)
         hessian = hessian + k_weights(ik)*h/weight_sum
         torque_torque = torque_torque + k_weights(ik)*tt/weight_sum
         mixed_contact = mixed_contact + k_weights(ik)*cc/weight_sum
         complete = complete + k_weights(ik)*allh/weight_sum
         if (present(report)) then
            report%contour_points = report%contour_points + one_report%contour_points
            report%fermi_poles = report%fermi_poles + one_report%fermi_poles
            report%solve_seconds = report%solve_seconds + one_report%solve_seconds
            report%contour_seconds = report%contour_seconds + one_report%contour_seconds
            report%pole_seconds = report%pole_seconds + one_report%pole_seconds
         end if
      end do
      deallocate(h, tt, cc, allh)
   end subroutine force_theorem_finite_q_hessian_from_resolvent_batch

   !> Zero-temperature occupied-contour version for a gapped fixture.
   module subroutine force_theorem_finite_q_hessian_from_zero_temperature_contour(h_source, h_endpoint, occupied_bounds, torques_q, &
      torques_minus_q, mixed, options, hessian, torque_torque, mixed_contact, complete, report)
      complex(rp), intent(in) :: h_source(:, :), h_endpoint(:, :)
      real(rp), intent(in) :: occupied_bounds(2)
      complex(rp), intent(in) :: torques_q(:, :, :), torques_minus_q(:, :, :), mixed(:, :, :, :)
      type(finite_h_contour_options), intent(in) :: options
      complex(rp), intent(out) :: hessian(:, :), torque_torque(:, :), mixed_contact(:, :), complete(:, :)
      type(finite_h_contour_report), intent(out), optional :: report
      complex(rp), allocatable :: nodes(:), weights(:), poles(:)

      call build_zero_temperature_occupied_contour(occupied_bounds, options, nodes, weights)
      allocate(poles(0))
      call resolvent_hessian_core(h_source, h_endpoint, 0.0_rp, 0.0_rp, .false., nodes, weights, poles, torques_q, torques_minus_q, mixed, &
         hessian, torque_torque, mixed_contact, complete, report)
      deallocate(nodes, weights, poles)
   end subroutine force_theorem_finite_q_hessian_from_zero_temperature_contour

   ! The common contraction is intentionally expressed in orbital matrix space
   ! rather than in an eigenbasis.  Each energy node performs direct LU
   ! factorization of z I-H and solves for the vertex/contact right-hand sides.
   subroutine resolvent_hessian_core(h_source, h_endpoint, fermi, kT, finite_temperature, nodes, weights, poles, torques_q, &
      torques_minus_q, mixed, hessian, torque_torque, mixed_contact, complete, report)
      complex(rp), intent(in) :: h_source(:, :), h_endpoint(:, :), nodes(:), weights(:), poles(:)
      real(rp), intent(in) :: fermi, kT
      logical, intent(in) :: finite_temperature
      complex(rp), intent(in) :: torques_q(:, :, :), torques_minus_q(:, :, :), mixed(:, :, :, :)
      complex(rp), intent(out) :: hessian(:, :), torque_torque(:, :), mixed_contact(:, :), complete(:, :)
      type(finite_h_contour_report), intent(out), optional :: report
      complex(rp), allocatable :: a0(:, :), a1(:, :), rhs0(:, :), rhs1(:, :)
      complex(rp) :: coefficient, trace_first, trace_second, trace_contact
      complex(rp), allocatable :: contour_tt(:,:), contour_contact(:,:), pole_tt(:,:), pole_contact(:,:)
      integer, allocatable :: piv0(:), piv1(:)
      integer :: nmat, nsite, nnode, npole, a, b, info, clock_start, clock_stop, clock_rate
      integer :: contour_start, contour_stop, pole_start, pole_stop
      real(rp) :: solve_seconds, contour_seconds, pole_seconds

      nmat = size(h_source,1); nsite = size(torques_q,3); nnode = size(nodes); npole = size(poles)

      allocate(a0(nmat,nmat), a1(nmat,nmat), rhs0(nmat,nmat), rhs1(nmat,nmat), piv0(nmat), piv1(nmat), &
         contour_tt(nsite,nsite), contour_contact(nsite,nsite), pole_tt(nsite,nsite), pole_contact(nsite,nsite))
      contour_tt = 0.0_rp; contour_contact = 0.0_rp; pole_tt = 0.0_rp; pole_contact = 0.0_rp
      solve_seconds = 0.0_rp; contour_seconds = 0.0_rp; pole_seconds = 0.0_rp
      call system_clock(count_rate=clock_rate)

      call system_clock(contour_start)
      do nnode = 1, size(nodes)
         a0 = -h_source; a1 = -h_endpoint
         do a = 1, nmat
            a0(a,a) = a0(a,a) + nodes(nnode)
            a1(a,a) = a1(a,a) + nodes(nnode)
         end do
         call system_clock(clock_start)
         call zgetrf(nmat, nmat, a0, nmat, piv0, info)
         if (info /= 0) error stop 'resolvent Hessian: source LU factorization failed'
         call zgetrf(nmat, nmat, a1, nmat, piv1, info)
         if (info /= 0) error stop 'resolvent Hessian: endpoint LU factorization failed'
         call system_clock(clock_stop)
         solve_seconds = solve_seconds + elapsed_seconds(clock_start,clock_stop,clock_rate)
         coefficient = weights(nnode)
         if (finite_temperature) coefficient = coefficient*finite_temperature_regularized_fermi(nodes(nnode),fermi,kT,poles)
         call accumulate_at_factorized_node(a0,a1,piv0,piv1,coefficient,contour_tt,contour_contact)
      end do
      call system_clock(contour_stop)
      contour_seconds = elapsed_seconds(contour_start,contour_stop,clock_rate)

      ! The regularized contour integrand has no enclosed Fermi poles, but the
      ! physical-spectrum contour is recovered by adding their residues.  The
      ! explicit solve here is part of the exact finite-T identity, not a
      ! quadrature correction or a spectral reconstruction.
      call system_clock(pole_start)
      do nnode = 1, size(poles)
         a0 = -h_source; a1 = -h_endpoint
         do a = 1, nmat
            a0(a,a) = a0(a,a) + poles(nnode)
            a1(a,a) = a1(a,a) + poles(nnode)
         end do
         call system_clock(clock_start)
         call zgetrf(nmat, nmat, a0, nmat, piv0, info)
         if (info /= 0) error stop 'resolvent Hessian: source pole LU factorization failed'
         call zgetrf(nmat, nmat, a1, nmat, piv1, info)
         if (info /= 0) error stop 'resolvent Hessian: endpoint pole LU factorization failed'
         call system_clock(clock_stop)
         solve_seconds = solve_seconds + elapsed_seconds(clock_start,clock_stop,clock_rate)
         call accumulate_at_factorized_node(a0,a1,piv0,piv1,cmplx(kT,0.0_rp,rp),pole_tt,pole_contact)
      end do
      call system_clock(pole_stop)
      pole_seconds = elapsed_seconds(pole_start,pole_stop,clock_rate)

      torque_torque = contour_tt + pole_tt
      mixed_contact = contour_contact + pole_contact
      hessian = torque_torque + mixed_contact
      complete = hessian
      if (present(report)) then
         report%contour_points = size(nodes)
         report%fermi_poles = size(poles)
         report%solve_seconds = solve_seconds
         report%contour_seconds = contour_seconds
         report%pole_seconds = pole_seconds
      end if
      deallocate(a0,a1,rhs0,rhs1,piv0,piv1,contour_tt,contour_contact,pole_tt,pole_contact)

   contains

      subroutine accumulate_at_factorized_node(lu0,lu1,ip0,ip1,weight,tt_sum,contact_sum)
         complex(rp), intent(in) :: lu0(:, :), lu1(:, :), weight
         integer, intent(in) :: ip0(:), ip1(:)
         complex(rp), intent(inout) :: tt_sum(:,:), contact_sum(:,:)
         integer :: ia, ib, solve_start, solve_stop, local_info

         do ia = 1, nsite
            do ib = 1, nsite
               rhs0 = torques_minus_q(:,:,ib)
               call system_clock(solve_start)
               call zgetrs('N',nmat,nmat,lu0,nmat,ip0,rhs0,nmat,local_info)
               if (local_info /= 0) error stop 'resolvent Hessian: source vertex solve failed'
               rhs1 = matmul(torques_q(:,:,ia),rhs0)
               call zgetrs('N',nmat,nmat,lu1,nmat,ip1,rhs1,nmat,local_info)
               if (local_info /= 0) error stop 'resolvent Hessian: endpoint vertex solve failed'
               trace_first = trace_matrix(rhs1)

               rhs0 = torques_minus_q(:,:,ia)
               call zgetrs('N',nmat,nmat,lu0,nmat,ip0,rhs0,nmat,local_info)
               if (local_info /= 0) error stop 'resolvent Hessian: reverse source vertex solve failed'
               rhs1 = matmul(torques_q(:,:,ib),rhs0)
               call zgetrs('N',nmat,nmat,lu1,nmat,ip1,rhs1,nmat,local_info)
               if (local_info /= 0) error stop 'resolvent Hessian: reverse endpoint vertex solve failed'
               trace_second = trace_matrix(rhs1)

               rhs0 = mixed(:,:,ia,ib)
               call zgetrs('N',nmat,nmat,lu0,nmat,ip0,rhs0,nmat,local_info)
               if (local_info /= 0) error stop 'resolvent Hessian: contact solve failed'
               trace_contact = trace_matrix(rhs0)
               tt_sum(ia,ib) = tt_sum(ia,ib) + weight*0.5_rp*(trace_first+trace_second)
               contact_sum(ia,ib) = contact_sum(ia,ib) + weight*trace_contact
               call system_clock(solve_stop)
               solve_seconds = solve_seconds + elapsed_seconds(solve_start,solve_stop,clock_rate)
            end do
         end do
      end subroutine accumulate_at_factorized_node
   end subroutine resolvent_hessian_core

   subroutine validate_options(options)
      type(finite_h_contour_options), intent(in) :: options
      if (options%contour_points < 8) error stop 'finite-H contour: contour_points must be at least 8'
      if (trim(options%contour_shape) /= 'ellipse') error stop 'finite-H contour: only contour_shape=ellipse is implemented'
      if (options%contour_margin <= 0.0_rp) error stop 'finite-H contour: contour_margin must be positive'
      if (options%contour_height_fraction <= 0.0_rp) error stop 'finite-H contour: height fraction must be positive'
   end subroutine validate_options

   subroutine hermitian_bounds(h_source, h_endpoint, lower, upper)
      complex(rp), intent(in) :: h_source(:, :), h_endpoint(:, :)
      real(rp), intent(out) :: lower, upper
      real(rp) :: lo0, hi0, lo1, hi1
      call matrix_gershgorin_bounds(h_source,lo0,hi0)
      call matrix_gershgorin_bounds(h_endpoint,lo1,hi1)
      lower = min(lo0,lo1); upper = max(hi0,hi1)
   end subroutine hermitian_bounds

   subroutine matrix_gershgorin_bounds(matrix, lower, upper)
      complex(rp), intent(in) :: matrix(:, :)
      real(rp), intent(out) :: lower, upper
      real(rp) :: radius
      integer :: i
      lower = huge(1.0_rp); upper = -huge(1.0_rp)
      do i = 1, size(matrix,1)
         radius = sum(abs(matrix(i,:))) - abs(matrix(i,i))
         lower = min(lower,real(matrix(i,i),rp)-radius)
         upper = max(upper,real(matrix(i,i),rp)+radius)
      end do
   end subroutine matrix_gershgorin_bounds

   pure function trace_matrix(matrix) result(value)
      complex(rp), intent(in) :: matrix(:, :)
      complex(rp) :: value
      integer :: i
      value = 0.0_rp
      do i = 1, min(size(matrix,1),size(matrix,2))
         value = value + matrix(i,i)
      end do
   end function trace_matrix

   pure function elapsed_seconds(start_count, end_count, rate) result(seconds)
      integer, intent(in) :: start_count, end_count, rate
      real(rp) :: seconds
      if (rate > 0) then
         seconds = real(end_count-start_count,rp)/real(rate,rp)
      else
         seconds = 0.0_rp
      end if
   end function elapsed_seconds


   ! --- from lr_rotation_response ---

   !> Freeze the accepted second-order state and pretransform both endpoints'
   !> torque vertices once.  The band-pair response then preserves Pi_AB and
   !> Pi_BA independently at every requested frequency.
   module subroutine prepare_rotation_response(fixture, recip, q, state)
      type(lmto_live_hamiltonian_fixture), target, intent(in) :: fixture
      type(reciprocal), target, intent(inout) :: recip
      real(rp), intent(in) :: q(3)
      type(rotation_state), intent(inout) :: state
      real(rp), allocatable :: points(:, :), folded(:, :), axes(:, :)
      complex(rp), allocatable :: vectors(:, :, :), endpoint_vectors(:, :, :), vertex(:, :), mixed(:, :), band(:, :)
      integer :: ik, ia, ib, n, m, isite_a, isite_b, axis_a, axis_b, norb_site, i0, info
      real(rp) :: f, spin_z, norm_m, weight

      call state%clear()
      if (fixture%nsite < 1 .or. fixture%norb < 1 .or. .not. fixture%hoh .or. .not. fixture%include_enu) &
         error stop 'STATE_PROVENANCE_OPEN: response requires accepted H=B-QB+E_nu fixture'
      if (trim(recip%reciprocal_mode) /= 'ham_only' .or. trim(recip%kspace_ham_order) /= 'second') &
         error stop 'STATE_PROVENANCE_OPEN: response requires second-order ham_only reciprocal state'
      if (recip%include_so) error stop 'STATE_PROVENANCE_OPEN: native rotation response currently requires SOC off'
      if (.not. allocated(recip%k_points) .or. .not. allocated(recip%k_weights)) &
         error stop 'STATE_PROVENANCE_OPEN: reciprocal mesh is unavailable'
      if (fixture%nsite /= recip%lattice%nrec) error stop 'STATE_PROVENANCE_OPEN: fixture/state site counts differ'
      call recip%require_replicated_k_workset('prepare_rotation_response')

      state%fixture => fixture
      state%reciprocal_state => recip
      state%q = q
      state%nsite = fixture%nsite
      state%ncoord = 2*fixture%nsite
      state%nmat = 2*fixture%norb*fixture%nsite
      state%nk = size(recip%k_points,2)
      state%nbands = state%nmat
      state%fermi = recip%fermi_level
      state%kT = max(recip%temperature*kB_ry_per_k, 1.0e-10_rp)
      state%weight_sum = sum(recip%k_weights)
      if (state%nk < 1 .or. size(recip%k_weights) /= state%nk .or. state%weight_sum <= tiny(1.0_rp)) &
         error stop 'STATE_PROVENANCE_OPEN: invalid reciprocal k weights'
      if (abs(state%weight_sum-1.0_rp) > 1.0e-8_rp) &
         error stop 'STATE_PROVENANCE_OPEN: reciprocal k weights are not normalized'

      allocate(state%values(state%nbands,state%nk), state%endpoint_values(state%nbands,state%nk), &
         state%weights(state%nk), state%occupations(state%nbands,state%nk), &
         state%endpoint_occupations(state%nbands,state%nk), &
         state%torque_band(state%nbands,state%nbands,state%ncoord,state%nk), &
         state%partner_band(state%nbands,state%nbands,state%ncoord,state%nk), &
         state%contact(state%ncoord,state%ncoord))
      state%weights = recip%k_weights/state%weight_sum

      ! The normal endpoint comes from the accepted SCF eigensystem whenever
      ! it is present; the arbitrary-point API supplies every k+q endpoint.
      if (allocated(recip%eigenvalues) .and. allocated(recip%eigenvectors)) then
         if (all(shape(recip%eigenvalues) == [state%nbands,state%nk]) .and. &
             all(shape(recip%eigenvectors) == [state%nmat,state%nbands,state%nk])) then
            state%values = recip%eigenvalues
            allocate(vectors(state%nmat,state%nbands,state%nk))
            vectors = recip%eigenvectors
         else
            call recip%calculate_eigenpairs_at_kpoints(recip%k_points,state%values,vectors)
         end if
      else
         call recip%calculate_eigenpairs_at_kpoints(recip%k_points,state%values,vectors)
      end if
      if (maxval(abs(q)) <= q_tolerance) then
         state%endpoint_values = state%values
         allocate(endpoint_vectors(state%nmat,state%nbands,state%nk))
         endpoint_vectors = vectors
      else
         points = recip%k_points + spread(q,2,state%nk)
         call recip%calculate_eigenpairs_at_kpoints(points,state%endpoint_values,endpoint_vectors,folded)
      end if

      do ik = 1, state%nk
         do n = 1, state%nbands
            state%occupations(n,ik) = finite_temperature_occupation(state%values(n,ik),state%fermi,state%kT)
            state%endpoint_occupations(n,ik) = finite_temperature_occupation(state%endpoint_values(n,ik),state%fermi,state%kT)
         end do
      end do
      call rotation_axes(fixture%moments,axes)
      state%contact = cmplx(0.0_rp,0.0_rp,rp)
      state%magnetization = 0.0_rp
      state%electron_count = 0.0_rp
      norb_site = fixture%norb
      allocate(vertex(state%nmat,state%nmat),mixed(state%nmat,state%nmat),band(state%nbands,state%nbands))

      do ik = 1, state%nk
         weight = state%weights(ik)
         do n = 1, state%nbands
            f = state%occupations(n,ik)
            state%electron_count = state%electron_count + weight*f
            spin_z = 0.0_rp
            do isite_a = 1, state%nsite
               i0 = (isite_a-1)*2*norb_site
               spin_z = spin_z + sum(abs(vectors(i0+1:i0+norb_site,n,ik))**2) - &
                  sum(abs(vectors(i0+norb_site+1:i0+2*norb_site,n,ik))**2)
            end do
            state%magnetization = state%magnetization + weight*f*spin_z
         end do

         do ia = 1, state%ncoord
            isite_a = (ia+1)/2
            axis_a = 1+mod(ia-1,2)
            call assemble_lmto_finite_q_torque(fixture,recip%k_points(:,ik),q,isite_a,axes(:,2*(isite_a-1)+axis_a),vertex)
            band = matmul(conjg(transpose(endpoint_vectors(:,:,ik))),matmul(vertex,vectors(:,:,ik)))
            state%torque_band(:,:,ia,ik) = band
            call assemble_lmto_finite_q_torque(fixture,recip%k_points(:,ik)+q,-q,isite_a, &
               axes(:,2*(isite_a-1)+axis_a),vertex)
            band = matmul(conjg(transpose(vectors(:,:,ik))),matmul(vertex,endpoint_vectors(:,:,ik)))
            state%partner_band(:,:,ia,ik) = band
         end do

         do ia = 1, state%ncoord
            isite_a = (ia+1)/2
            axis_a = 1+mod(ia-1,2)
            do ib = 1, state%ncoord
               isite_b = (ib+1)/2
               axis_b = 1+mod(ib-1,2)
               call assemble_lmto_finite_q_mixed_derivative(fixture,recip%k_points(:,ik),q,isite_a, &
                  axes(:,2*(isite_a-1)+axis_a),isite_b,axes(:,2*(isite_b-1)+axis_b),mixed)
               band = matmul(conjg(transpose(vectors(:,:,ik))),matmul(mixed,vectors(:,:,ik)))
               do n = 1, state%nbands
                  state%contact(ia,ib) = state%contact(ia,ib) + weight*state%occupations(n,ik)*band(n,n)
               end do
            end do
         end do
      end do
      call global_spin_berry(vectors,state%occupations,state%weights,state%nsite,norb_site,state%berry)
      state%berry_residual=abs(state%berry-0.5_rp*state%magnetization)/ &
         max(abs(state%berry),abs(state%magnetization),1.0e-12_rp)
      state%electron_residual = state%electron_count-recip%total_electrons
      state%prepared = .true.
      deallocate(vectors,endpoint_vectors,axes,vertex,mixed,band)
      if (allocated(points)) deallocate(points)
      if (allocated(folded)) deallocate(folded)
   end subroutine prepare_rotation_response

   !> Evaluate the literal unsymmetrized retarded band sum.  At omega=eta=0,
   !> exact_static selects the finite-temperature divided-difference branch.
   module subroutine evaluate_rotation_response(request,result)
      type(rotation_request), intent(in) :: request
      type(rotation_result), intent(out) :: result
      type(rotation_state), pointer :: state
      complex(rp), allocatable :: unitary(:, :)
      complex(rp) :: denominator, product
      real(rp) :: fn, fm, divided, energy_difference, weight
      integer :: ia, ib, ik, n, m, info, ncoord

      state => request%state
      if (.not. state%prepared) error stop 'rotation response state is not prepared'
      if (request%exact_static) then
         if (abs(request%omega) > q_tolerance .or. abs(request%eta) > q_tolerance) &
            error stop 'exact static response requires omega=eta=0'
      else if (request%eta <= 0.0_rp .or. .not.ieee_is_finite(request%eta)) then
         error stop 'retarded rotation response requires finite positive eta'
      end if
      ncoord = state%ncoord
      result%q = state%q
      result%omega = request%omega
      result%eta = request%eta
      result%berry = state%berry
      result%magnetization = state%magnetization
      result%berry_residual = state%berry_residual
      result%electron_count = state%electron_count
      result%electron_residual = state%electron_residual
      allocate(result%bubble(ncoord,ncoord),result%contact(ncoord,ncoord), &
         result%kernel(ncoord,ncoord),result%kernel_pm(ncoord,ncoord),unitary(ncoord,ncoord))
      result%bubble = cmplx(0.0_rp,0.0_rp,rp)
      do ik = 1, state%nk
         weight = state%weights(ik)
         do ia = 1, ncoord
            do ib = 1, ncoord
               do n = 1, state%nbands
                  fn = state%occupations(n,ik)
                  do m = 1, state%nbands
                     fm = state%endpoint_occupations(m,ik)
                     product = state%torque_band(m,n,ia,ik)*state%partner_band(n,m,ib,ik)
                     energy_difference = state%values(n,ik)-state%endpoint_values(m,ik)
                     if (request%exact_static) then
                        divided = fermi_divided_difference(state%values(n,ik),state%endpoint_values(m,ik), &
                           state%fermi,state%kT)
                        result%bubble(ia,ib) = result%bubble(ia,ib)+weight*divided*product
                     else
                        denominator = cmplx(request%omega+energy_difference,request%eta,rp)
                        result%bubble(ia,ib) = result%bubble(ia,ib)+weight*(fn-fm)*product/denominator
                     end if
                  end do
               end do
            end do
         end do
      end do
      result%contact = state%contact
      result%kernel = result%contact+result%bubble
      call rotation_circular_unitary(state%nsite,unitary)
      result%kernel_pm = matmul(unitary,matmul(result%kernel,conjg(transpose(unitary))))
      if (request%want_inverse) then
         ! The exact q=0 static Goldstone kernel is singular by construction.
         if (request%exact_static .and. maxval(abs(state%q)) <= q_tolerance) then
            result%inverse_available = .false.
         else
            allocate(result%inverse_kernel(ncoord,ncoord),result%inverse_kernel_pm(ncoord,ncoord))
            result%inverse_kernel = result%kernel
            call invert_complex_matrix(result%inverse_kernel,info)
            if (info == 0) then
               result%inverse_available = .true.
               result%inverse_kernel_pm = matmul(unitary, &
                  matmul(result%inverse_kernel,conjg(transpose(unitary))))
            else
               deallocate(result%inverse_kernel,result%inverse_kernel_pm)
               result%inverse_available = .false.
            end if
         end if
      end if
      deallocate(unitary)
   end subroutine evaluate_rotation_response

   !> Independent literal band-loop oracle.  It assembles orbital-space
   !> vertices and evaluates <m,k+q|T_A|n,k><n,k|T_B|m,k+q> directly; it does
   !> not call or consume the production accumulator or its transformed bands.
   module subroutine evaluate_rotation_response_oracle(fixture,recip,q,omega,eta,bubble,contact,kernel)
      type(lmto_live_hamiltonian_fixture), intent(in) :: fixture
      type(reciprocal), intent(inout) :: recip
      real(rp), intent(in) :: q(3),omega,eta
      complex(rp), intent(out) :: bubble(:, :),contact(:, :),kernel(:, :)
      real(rp), allocatable :: values(:, :),endpoint_values(:, :),points(:, :),axes(:, :)
      complex(rp), allocatable :: vectors(:, :, :),endpoint_vectors(:, :, :),vq(:, :, :),vm(:, :, :),mixed(:, :)
      complex(rp) :: amp_a,amp_b,amp_c,denominator
      real(rp) :: weight_sum,weight,fermi,kT,fn,fm
      integer :: nk,nband,nmat,ncoord,nsite,ik,ia,ib,n,m,site_a,site_b,axis_a,axis_b
      integer :: norb_site,i0

      if (eta<=0.0_rp) error stop 'dynamic band oracle requires positive eta'
      if (.not.fixture%hoh .or. .not.fixture%include_enu .or. trim(recip%kspace_ham_order)/='second' .or. &
          trim(recip%reciprocal_mode)/='ham_only' .or. recip%include_so) &
         error stop 'STATE_PROVENANCE_OPEN: direct oracle requires second-order ham_only state with SOC off'
      call recip%require_replicated_k_workset('evaluate_rotation_response_oracle')
      nk=size(recip%k_points,2); nsite=fixture%nsite; ncoord=2*nsite; nmat=2*fixture%norb*nsite; nband=nmat
      fermi=recip%fermi_level; kT=max(recip%temperature*kB_ry_per_k,1.0e-10_rp); weight_sum=sum(recip%k_weights)
      allocate(values(nband,nk),endpoint_values(nband,nk),vectors(nmat,nband,nk), &
         endpoint_vectors(nmat,nband,nk),axes(3,ncoord),vq(nmat,nmat,ncoord),vm(nmat,nmat,ncoord),mixed(nmat,nmat))
      if (allocated(recip%eigenvalues) .and. allocated(recip%eigenvectors)) then
         if (all(shape(recip%eigenvalues)==[nband,nk]) .and. all(shape(recip%eigenvectors)==[nmat,nband,nk])) then
            values=recip%eigenvalues; vectors=recip%eigenvectors
         else
            call recip%calculate_eigenpairs_at_kpoints(recip%k_points,values,vectors)
         end if
      else
         call recip%calculate_eigenpairs_at_kpoints(recip%k_points,values,vectors)
      end if
      if (maxval(abs(q))<=q_tolerance) then
         endpoint_values=values; endpoint_vectors=vectors
      else
         points=recip%k_points+spread(q,2,nk)
         call recip%calculate_eigenpairs_at_kpoints(points,endpoint_values,endpoint_vectors)
         deallocate(points)
      end if
      call rotation_axes(fixture%moments,axes)
      bubble=cmplx(0.0_rp,0.0_rp,rp); contact=cmplx(0.0_rp,0.0_rp,rp)
      norb_site=fixture%norb
      do ik=1,nk
         weight=recip%k_weights(ik)/weight_sum
         do ia=1,ncoord
            site_a=(ia+1)/2; axis_a=1+mod(ia-1,2)
            call assemble_lmto_finite_q_torque(fixture,recip%k_points(:,ik),q,site_a, &
               axes(:,2*(site_a-1)+axis_a),vq(:,:,ia))
            call assemble_lmto_finite_q_torque(fixture,recip%k_points(:,ik)+q,-q,site_a, &
               axes(:,2*(site_a-1)+axis_a),vm(:,:,ia))
         end do
         do ia=1,ncoord
            site_a=(ia+1)/2; axis_a=1+mod(ia-1,2)
            do ib=1,ncoord
               site_b=(ib+1)/2; axis_b=1+mod(ib-1,2)
               call assemble_lmto_finite_q_mixed_derivative(fixture,recip%k_points(:,ik),q,site_a, &
                  axes(:,2*(site_a-1)+axis_a),site_b,axes(:,2*(site_b-1)+axis_b),mixed)
               do n=1,nband
                  fn=finite_temperature_occupation(values(n,ik),fermi,kT)
                  amp_c=dot_product(vectors(:,n,ik),matmul(mixed,vectors(:,n,ik)))
                  contact(ia,ib)=contact(ia,ib)+weight*fn*amp_c
                  do m=1,nband
                     fm=finite_temperature_occupation(endpoint_values(m,ik),fermi,kT)
                     denominator=cmplx(omega+values(n,ik)-endpoint_values(m,ik),eta,rp)
                     amp_a=dot_product(endpoint_vectors(:,m,ik),matmul(vq(:,:,ia),vectors(:,n,ik)))
                     amp_b=dot_product(vectors(:,n,ik),matmul(vm(:,:,ib),endpoint_vectors(:,m,ik)))
                     bubble(ia,ib)=bubble(ia,ib)+weight*(fn-fm)*amp_a*amp_b/denominator
                  end do
               end do
            end do
         end do
      end do
      kernel=contact+bubble
      deallocate(values,endpoint_values,vectors,endpoint_vectors,axes,vq,vm,mixed)
   end subroutine evaluate_rotation_response_oracle

   !> The stored torque orientation makes the exact static force-theorem
   !> Hessian the Cartesian-index symmetric part at the same q:
   !> H_AB(q)=0.5*(K_AB^R(q,0)+K_BA^R(q,0)).  The q/-q covariance is a
   !> separate identity; inserting a conjugation here retains a spurious
   !> antisymmetric-in-frequency component at finite q.
   module pure subroutine reduce_static_rotation_kernel(kernel_q,reduced)
      complex(rp), intent(in) :: kernel_q(:, :)
      complex(rp), intent(out) :: reduced(:, :)
      integer :: a,b,n
      n = size(kernel_q,1)
      if (size(kernel_q,2)/=n .or. any(shape(reduced)/=[n,n])) &
         error stop 'reduce_static_rotation_kernel: shape mismatch'
      do a = 1,n
         do b = 1,n
            reduced(a,b)=0.5_rp*(kernel_q(a,b)+kernel_q(b,a))
         end do
      end do
   end subroutine reduce_static_rotation_kernel

   !> Evaluate -i Tr rho [Gx,Gy] directly in the accepted band basis.
   !> The spin generators are global spin-1/2 generators replicated over
   !> sites/orbitals; the result should independently equal M_band/2.
   subroutine global_spin_berry(vectors,occupations,weights,nsite,norb,berry)
      complex(rp), intent(in) :: vectors(:, :, :)
      real(rp), intent(in) :: occupations(:, :),weights(:)
      integer, intent(in) :: nsite,norb
      real(rp), intent(out) :: berry
      complex(rp), allocatable :: gx(:, :),gy(:, :),commutator(:, :)
      complex(rp) :: trace_value
      integer :: nmat,nband,nk,site,io,up,dn,ik,n
      nmat=2*nsite*norb; nband=size(vectors,2); nk=size(vectors,3)
      if (any(shape(vectors)/=[nmat,nband,nk]) .or. any(shape(occupations)/=[nband,nk]) .or. size(weights)/=nk) &
         error stop 'global_spin_berry: shape mismatch'
      allocate(gx(nmat,nmat),gy(nmat,nmat),commutator(nmat,nmat))
      gx=cmplx(0.0_rp,0.0_rp,rp); gy=gx
      do site=1,nsite
         do io=1,norb
            up=(site-1)*2*norb+io; dn=(site-1)*2*norb+norb+io
            gx(up,dn)=0.5_rp; gx(dn,up)=0.5_rp
            gy(up,dn)=-0.5_rp*i_unit; gy(dn,up)=0.5_rp*i_unit
         end do
      end do
      commutator=matmul(gx,gy)-matmul(gy,gx)
      berry=0.0_rp
      do ik=1,nk
         do n=1,nband
            trace_value=dot_product(vectors(:,n,ik),matmul(commutator,vectors(:,n,ik)))
            berry=berry+weights(ik)*occupations(n,ik)*real(-i_unit*trace_value,rp)
         end do
      end do
      deallocate(gx,gy,commutator)
   end subroutine global_spin_berry

   !> Construct local transverse rotation axes from each accepted unit moment.
   module pure subroutine rotation_axes(moments,axes)
      real(rp), intent(in) :: moments(:, :)
      real(rp), allocatable, intent(out) :: axes(:, :)
      real(rp) :: m(3), reference(3), e1(3), e2(3), norm_m, dot_m
      integer :: site
      allocate(axes(3,2*size(moments,2)))
      do site = 1,size(moments,2)
         m = moments(:,site)
         norm_m = sqrt(dot_product(m,m))
         if (norm_m<=tiny(1.0_rp)) error stop 'rotation_axes: zero site moment'
         m = m/norm_m
         if (abs(m(1))<0.8_rp) then
            reference=[1.0_rp,0.0_rp,0.0_rp]
         else
            reference=[0.0_rp,1.0_rp,0.0_rp]
         end if
         dot_m=dot_product(reference,m)
         e1=reference-dot_m*m
         e1=e1/sqrt(dot_product(e1,e1))
         e2=[m(2)*e1(3)-m(3)*e1(2),m(3)*e1(1)-m(1)*e1(3),m(1)*e1(2)-m(2)*e1(1)]
         axes(:,2*site-1)=e1
         axes(:,2*site)=e2
      end do
   end subroutine rotation_axes

   !> theta_+=(theta_x-i theta_y)/sqrt(2), theta_-=(theta_x+i theta_y)/sqrt(2).
   module pure subroutine rotation_circular_unitary(nsite,unitary)
      integer, intent(in) :: nsite
      complex(rp), intent(out) :: unitary(:, :)
      complex(rp), parameter :: inv_sqrt2=cmplx(1.0_rp/sqrt(2.0_rp),0.0_rp,rp)
      integer :: site, i0
      unitary=cmplx(0.0_rp,0.0_rp,rp)
      do site=1,nsite
         i0=2*(site-1)
         unitary(i0+1,i0+1)=inv_sqrt2
         unitary(i0+1,i0+2)=-i_unit*inv_sqrt2
         unitary(i0+2,i0+1)=inv_sqrt2
         unitary(i0+2,i0+2)= i_unit*inv_sqrt2
      end do
   end subroutine rotation_circular_unitary

   subroutine invert_complex_matrix(matrix,info)
      complex(rp), intent(inout) :: matrix(:, :)
      integer, intent(out) :: info
      complex(rp), allocatable :: work(:)
      complex(rp) :: query(1)
      integer, allocatable :: pivots(:)
      integer :: n,lwork
      n=size(matrix,1)
      allocate(pivots(n))
      call zgetrf(n,n,matrix,n,pivots,info)
      if (info/=0) then
         deallocate(pivots)
         return
      end if
      call zgetri(n,matrix,n,pivots,query,-1,info)
      if (info/=0) then
         deallocate(pivots)
         return
      end if
      lwork=max(n,int(real(query(1),rp)))
      allocate(work(lwork))
      call zgetri(n,matrix,n,pivots,work,lwork,info)
      deallocate(pivots,work)
   end subroutine invert_complex_matrix

   module subroutine rotation_state_clear(this)
      class(rotation_state), intent(inout) :: this
      if (allocated(this%values)) deallocate(this%values)
      if (allocated(this%endpoint_values)) deallocate(this%endpoint_values)
      if (allocated(this%weights)) deallocate(this%weights)
      if (allocated(this%occupations)) deallocate(this%occupations)
      if (allocated(this%endpoint_occupations)) deallocate(this%endpoint_occupations)
      if (allocated(this%torque_band)) deallocate(this%torque_band)
      if (allocated(this%partner_band)) deallocate(this%partner_band)
      if (allocated(this%contact)) deallocate(this%contact)
      nullify(this%fixture,this%reciprocal_state)
      this%q=0.0_rp; this%nsite=0; this%ncoord=0; this%nmat=0; this%nbands=0; this%nk=0
      this%fermi=0.0_rp; this%kT=0.0_rp; this%weight_sum=0.0_rp
      this%magnetization=0.0_rp; this%berry=0.0_rp
      this%electron_count=0.0_rp; this%electron_residual=0.0_rp
      this%prepared=.false.
   end subroutine rotation_state_clear


   module subroutine compute_static_rotation_curvature(q_coordinates, q_list, rotation_axis, finite_h_spectral_mode, &
      finite_h_response_backend, contour_points, contour_shape, contour_margin, contour_height_fraction, &
      contour_account_fermi_poles, native_crosscheck, native_turek, native_contour_points, native_contour_margin, &
      native_contour_height_fraction, native_contour_account_fermi_poles, native_contour_target_fermi_poles, &
      lattice_obj, hamiltonian_obj, energy_obj, self_obj, reciprocal_obj, fixture, q_direct, q_cart, finite_total, &
      finite_tt, finite_contact, spectral_total, spectral_tt, spectral_contact, contour_total, contour_tt, contour_contact, &
      native_jq_ud, native_jq_du, native_jq_sym, native_delta_j, native_curvature, native_report, native_ready, &
      endpoint_mode, endpoint_reused, q_commensurate, endpoint_residual, endpoint_seconds, assembly_seconds, &
      contraction_seconds, hamiltonian_seconds, gf_seconds, solve_seconds, contour_seconds, total_response_seconds)
      character(len=*), intent(in) :: q_coordinates, finite_h_spectral_mode, finite_h_response_backend, contour_shape
      real(rp), intent(in) :: q_list(:, :), rotation_axis(3), contour_margin, contour_height_fraction
      logical, intent(in) :: contour_account_fermi_poles, native_crosscheck, native_turek
      integer, intent(in) :: contour_points, native_contour_points, native_contour_target_fermi_poles
      real(rp), intent(in) :: native_contour_margin, native_contour_height_fraction
      logical, intent(in) :: native_contour_account_fermi_poles
      type(lattice), intent(inout) :: lattice_obj
      type(hamiltonian), intent(inout) :: hamiltonian_obj
      type(energy), intent(in) :: energy_obj
      type(self), intent(inout) :: self_obj
      type(reciprocal), intent(inout) :: reciprocal_obj
      type(lmto_live_hamiltonian_fixture), intent(out) :: fixture
      real(rp), allocatable, intent(out) :: q_direct(:, :), q_cart(:, :)
      real(rp), allocatable, intent(out) :: finite_total(:, :, :), finite_tt(:, :, :), finite_contact(:, :, :)
      real(rp), allocatable, intent(out) :: spectral_total(:, :, :), spectral_tt(:, :, :), spectral_contact(:, :, :)
      real(rp), allocatable, intent(out) :: contour_total(:, :, :), contour_tt(:, :, :), contour_contact(:, :, :)
      real(rp), allocatable, intent(out) :: native_jq_ud(:), native_jq_du(:), native_jq_sym(:), native_delta_j(:), native_curvature(:)
      type(native_turek_contour_report), intent(out) :: native_report
      logical, intent(out) :: native_ready
      character(len=24), allocatable, intent(out) :: endpoint_mode(:)
      logical, allocatable, intent(out) :: endpoint_reused(:), q_commensurate(:)
      real(rp), allocatable, intent(out) :: endpoint_residual(:)
      real(rp), intent(out) :: endpoint_seconds, assembly_seconds, contraction_seconds, hamiltonian_seconds
      real(rp), intent(out) :: gf_seconds, solve_seconds, contour_seconds, total_response_seconds
      real(rp), allocatable :: evals(:, :), endpoint_evals(:, :), endpoint_eval_check(:, :), weights(:)
      complex(rp), allocatable :: evecs(:, :, :), endpoint_evecs(:, :, :), endpoint_evec_check(:, :, :)
      complex(rp), allocatable :: torques_q(:, :, :, :), torques_minus_q(:, :, :, :), mixed(:, :, :, :, :)
      complex(rp), allocatable :: torque_minus_check(:, :, :)
      complex(rp), allocatable :: hessian(:, :), torque_torque(:, :), contact(:, :), complete(:, :)
      complex(rp), allocatable :: h_source(:, :, :), h_endpoint(:, :, :)
      real(rp), allocatable :: axes(:, :)
      integer, allocatable :: endpoint_index(:)
      real(rp) :: adapter_before, adapter_after, axis_norm, finite_h_kT, vertex_error
      real(rp) :: endpoint_identity_error
      type(finite_h_contour_options) :: contour_options
      type(finite_h_contour_report) :: contour_report
      type(native_turek_contour_options) :: native_contour_options
      integer :: nq, nsite, nmat, nk, iq, ik, ia, ja, clock_start, clock_end, clock_rate, response_start
      logical :: vertex_identity_checked, endpoint_identity_checked, do_spectral, do_contour

      if (size(q_list, 1) /= 3 .or. size(q_list, 2) < 1) error stop 'compute_static_rotation_curvature: q path shape is invalid'
      nq = size(q_list, 2)
      call exchange_q_convert_points(q_coordinates, q_list, lattice_obj, q_direct, q_cart)
      nsite = hamiltonian_obj%charge%lattice%nrec
      axis_norm = sqrt(sum(rotation_axis**2))
      allocate(axes(3,nsite))
      do ia = 1, nsite
         axes(:,ia) = rotation_axis/axis_norm
      end do

      ! Freeze the finite-H accepted state before representation-specific
      ! native-LKAG preparation.  predls() is allowed to mutate symbolic-atom
      ! P/S channels, but must not alter the production Hamiltonian consumed by
      ! finite-H or reciprocal endpoint solves.
      call lmto_fixture_from_hamiltonian(hamiltonian_obj, fixture)
      call lmto_fixture_adapter_residual(hamiltonian_obj, fixture, reciprocal_obj%k_points(:,1:1), adapter_before)
      do ik = 1, lattice_obj%ntype
         call lattice_obj%symbolic_atoms(ik)%predls(lattice_obj%wav*ang2au)
      end do
      call lmto_fixture_adapter_residual(hamiltonian_obj, fixture, reciprocal_obj%k_points(:,1:1), adapter_after)
      call g_logger%info('[exchange_q]: accepted-H adapter residual before/after='//trim(real_to_string(adapter_before))// &
         '/'//trim(real_to_string(adapter_after)), __FILE__, __LINE__)
      if (adapter_before > 3.0e-12_rp .or. adapter_after > 3.0e-12_rp .or. &
          abs(adapter_after-adapter_before) > 1.0e-13_rp) then
         call g_logger%fatal('[exchange_q]: predls() changed the accepted production-H adapter state', __FILE__, __LINE__)
      end if
      if (abs(self_obj%reciprocal_scf_cache%fermi_level-energy_obj%fermi) > 3.0e-12_rp) then
         call g_logger%fatal('[exchange_q]: accepted Fermi-level provenance mismatch', __FILE__, __LINE__)
      end if

      nsite = fixture%nsite
      nmat = 2*fixture%norb*nsite
      nk = size(reciprocal_obj%k_points,2)
      do_spectral = finite_h_response_backend == 'spectral' .or. finite_h_response_backend == 'both'
      do_contour = finite_h_response_backend == 'contour' .or. finite_h_response_backend == 'both'
      if (nsite < 1 .or. nk < 1 .or. (do_spectral .and. .not. allocated(reciprocal_obj%eigenvectors))) then
         call g_logger%fatal('[exchange_q]: accepted reciprocal eigensystem is incomplete', __FILE__, __LINE__)
      end if
      if (do_spectral .and. size(reciprocal_obj%eigenvalues,1) /= nmat) then
         call g_logger%fatal('[exchange_q]: accepted eigensystem/Hamiltonian dimensions disagree', __FILE__, __LINE__)
      end if
      allocate(weights(nk))
      weights = reciprocal_obj%k_weights
      if (do_spectral) then
         allocate(evals(nmat,nk), evecs(nmat,nmat,nk))
         evals = reciprocal_obj%eigenvalues
         evecs = reciprocal_obj%eigenvectors
      end if
      if (size(weights) /= nk .or. sum(weights) <= tiny(1.0_rp)) then
         call g_logger%fatal('[exchange_q]: accepted k-mesh weights are invalid', __FILE__, __LINE__)
      end if

      allocate(finite_total(nsite,nsite,nq), finite_tt(nsite,nsite,nq), &
         finite_contact(nsite,nsite,nq))
      finite_total = 0.0_rp; finite_tt = 0.0_rp; finite_contact = 0.0_rp
      allocate(spectral_total(nsite,nsite,nq), spectral_tt(nsite,nsite,nq), spectral_contact(nsite,nsite,nq), &
         contour_total(nsite,nsite,nq), contour_tt(nsite,nsite,nq), contour_contact(nsite,nsite,nq))
      spectral_total = 0.0_rp; spectral_tt = 0.0_rp; spectral_contact = 0.0_rp
      contour_total = 0.0_rp; contour_tt = 0.0_rp; contour_contact = 0.0_rp
      allocate(hessian(nsite,nsite), torque_torque(nsite,nsite), contact(nsite,nsite), complete(nsite,nsite))
      allocate(endpoint_index(nk), endpoint_reused(nq), q_commensurate(nq), &
         endpoint_residual(nq), endpoint_mode(nq))
      endpoint_reused = .false.; q_commensurate = .false.; endpoint_residual = huge(1.0_rp)
      endpoint_mode = 'explicit_diagonalization'
      endpoint_seconds = 0.0_rp; assembly_seconds = 0.0_rp; contraction_seconds = 0.0_rp
      hamiltonian_seconds = 0.0_rp; gf_seconds = 0.0_rp; solve_seconds = 0.0_rp; contour_seconds = 0.0_rp; total_response_seconds = 0.0_rp
      call system_clock(count_rate=clock_rate)
      finite_h_kT = max(reciprocal_obj%temperature*kB_ry_per_k, 1.0e-10_rp)
      vertex_identity_checked = .false.
      endpoint_identity_checked = .false.
      native_ready = native_crosscheck .or. native_turek
      allocate(native_jq_ud(nq), native_jq_du(nq), native_jq_sym(nq), &
         native_delta_j(nq), native_curvature(nq))
      native_jq_ud = 0.0_rp; native_jq_du = 0.0_rp; native_jq_sym = 0.0_rp
      native_delta_j = 0.0_rp; native_curvature = 0.0_rp
      contour_options%contour_points = contour_points
      contour_options%contour_shape = contour_shape
      contour_options%contour_margin = contour_margin
      contour_options%contour_height_fraction = contour_height_fraction
      contour_options%account_fermi_poles = contour_account_fermi_poles
      native_contour_options%contour_points = native_contour_points
      native_contour_options%contour_shape = 'ellipse'
      native_contour_options%contour_margin = native_contour_margin
      native_contour_options%contour_height_fraction = native_contour_height_fraction
      native_contour_options%account_fermi_poles = native_contour_account_fermi_poles
      native_contour_options%target_fermi_poles = native_contour_target_fermi_poles

      call system_clock(response_start)

      do iq = 1, nq
         call detect_commensurate_endpoint_map(reciprocal_obj%k_points, reciprocal_obj%k_weights, reciprocal_obj%nk_mesh, &
            q_direct(:,iq), endpoint_index, q_commensurate(iq), endpoint_residual(iq))
         call system_clock(clock_start)
         if (do_spectral .and. q_commensurate(iq)) then
            allocate(endpoint_evals(nmat,nk), endpoint_evecs(nmat,nmat,nk))
            do ik = 1, nk
               endpoint_evals(:,ik) = evals(:,endpoint_index(ik))
               endpoint_evecs(:,:,ik) = evecs(:,:,endpoint_index(ik))
            end do
            endpoint_reused(iq) = .true.
            endpoint_mode(iq) = 'mesh_reuse'
            if (.not. endpoint_identity_checked .and. sqrt(sum(q_direct(:,iq)**2)) > q_tolerance) then
               ! The integer-cell map is a performance optimization, so verify
               ! its spectral identity once against the exact endpoint solve.
               call solve_unfolded_endpoints(reciprocal_obj, reciprocal_obj%k_points + &
                  spread_q(q_direct(:,iq),nk), endpoint_eval_check, endpoint_evec_check)
               endpoint_identity_error = maxval(abs(endpoint_evals-endpoint_eval_check))
               call g_logger%info('[exchange_q]: commensurate endpoint eigenvalue residual='// &
                  trim(real_to_string(endpoint_identity_error)), __FILE__, __LINE__)
               if (endpoint_identity_error > 2.0e-10_rp) then
                  call g_logger%fatal('[exchange_q]: commensurate endpoint spectral identity check failed', __FILE__, __LINE__)
               end if
               deallocate(endpoint_eval_check, endpoint_evec_check)
               endpoint_identity_checked = .true.
            end if
         else if (do_spectral) then
            call solve_unfolded_endpoints(reciprocal_obj, reciprocal_obj%k_points + spread_q(q_direct(:,iq),nk), &
                                          endpoint_evals, endpoint_evecs)
         end if
         call system_clock(clock_end)
         endpoint_seconds = endpoint_seconds + elapsed_clock_seconds(clock_start,clock_end,clock_rate)
         allocate(torques_q(nmat,nmat,nsite,nk), torques_minus_q(nmat,nmat,nsite,nk), &
                  mixed(nmat,nmat,nsite,nsite,nk))
         if (do_contour) allocate(h_source(nmat,nmat,nk), h_endpoint(nmat,nmat,nk))
         call system_clock(clock_start)
         if (do_contour) then
            do ik = 1, nk
               call assemble_lmto_hamiltonian(fixture, reciprocal_obj%k_points(:,ik), h_source(:,:,ik))
               call assemble_lmto_hamiltonian(fixture, reciprocal_obj%k_points(:,ik)+q_direct(:,iq), h_endpoint(:,:,ik))
            end do
         end if
         call system_clock(clock_end)
         hamiltonian_seconds = hamiltonian_seconds + elapsed_clock_seconds(clock_start,clock_end,clock_rate)
         call system_clock(clock_start)
         do ik = 1, nk
            call assemble_lmto_finite_q_torques(fixture, reciprocal_obj%k_points(:,ik), q_direct(:,iq), axes, &
               torques_q(:,:,:,ik))
            do ia = 1, nsite
               torques_minus_q(:,:,ia,ik) = transpose(conjg(torques_q(:,:,ia,ik)))
            end do
            do ia = 1, nsite
               do ja = 1, nsite
                  call assemble_lmto_finite_q_mixed_derivative(fixture, reciprocal_obj%k_points(:,ik), q_direct(:,iq), &
                     ia, axes(:,ia), ja, axes(:,ja), mixed(:,:,ia,ja,ik))
               end do
            end do
         end do
         if (.not. vertex_identity_checked .and. sqrt(sum(q_direct(:,iq)**2)) > q_tolerance) then
            allocate(torque_minus_check(nmat,nmat,nsite)); vertex_error = 0.0_rp
            call assemble_lmto_finite_q_torques(fixture, reciprocal_obj%k_points(:,1)+q_direct(:,iq), &
               -q_direct(:,iq), axes, torque_minus_check)
            do ia = 1, nsite
               vertex_error = max(vertex_error, maxval(abs(torque_minus_check(:,:,ia)-torques_minus_q(:,:,ia,1))))
            end do
            if (vertex_error > 5.0e-10_rp) call g_logger%fatal('[exchange_q]: finite-q vertex adjoint/gauge check failed', &
               __FILE__, __LINE__)
            vertex_identity_checked = .true.
            deallocate(torque_minus_check)
         end if
         call system_clock(clock_end)
         assembly_seconds = assembly_seconds + elapsed_clock_seconds(clock_start,clock_end,clock_rate)
         if (do_spectral) then
            call system_clock(clock_start)
            if (finite_h_spectral_mode == 'legacy_occupied') then
               call force_theorem_finite_q_hessian_from_eigenbasis_batch(evals, evecs, endpoint_evals, endpoint_evecs, &
                  reciprocal_obj%fermi_level, weights, torques_q, torques_minus_q, mixed, hessian, torque_torque, contact, complete)
            else
               call force_theorem_finite_q_hessian_from_eigenbasis_metallic_batch(evals, evecs, endpoint_evals, endpoint_evecs, &
                  reciprocal_obj%fermi_level, finite_h_kT, weights, torques_q, torques_minus_q, mixed, hessian, torque_torque, contact, complete)
            end if
            call system_clock(clock_end)
            contraction_seconds = contraction_seconds + elapsed_clock_seconds(clock_start,clock_end,clock_rate)
            spectral_tt(:,:,iq) = real(torque_torque,rp)
            spectral_contact(:,:,iq) = real(contact,rp)
            spectral_total(:,:,iq) = real(complete,rp)
         end if
         if (do_contour) then
            call system_clock(clock_start)
            call force_theorem_finite_q_hessian_from_resolvent_batch(h_source, h_endpoint, reciprocal_obj%fermi_level, finite_h_kT, &
               weights, torques_q, torques_minus_q, mixed, contour_options, hessian, torque_torque, contact, complete, contour_report)
            call system_clock(clock_end)
            gf_seconds = gf_seconds + elapsed_clock_seconds(clock_start,clock_end,clock_rate)
            solve_seconds = solve_seconds + contour_report%solve_seconds
            contour_seconds = contour_seconds + contour_report%contour_seconds
            contour_tt(:,:,iq) = real(torque_torque,rp)
            contour_contact(:,:,iq) = real(contact,rp)
            contour_total(:,:,iq) = real(complete,rp)
         end if
         if (finite_h_response_backend == 'contour') then
            finite_tt(:,:,iq) = contour_tt(:,:,iq); finite_contact(:,:,iq) = contour_contact(:,:,iq); finite_total(:,:,iq) = contour_total(:,:,iq)
         else
            finite_tt(:,:,iq) = spectral_tt(:,:,iq); finite_contact(:,:,iq) = spectral_contact(:,:,iq); finite_total(:,:,iq) = spectral_total(:,:,iq)
         end if
         if (allocated(h_source)) deallocate(h_source, h_endpoint)
         deallocate(torques_q, torques_minus_q, mixed)
         if (allocated(endpoint_evals)) deallocate(endpoint_evals, endpoint_evecs)
      end do
      call system_clock(clock_end)
      total_response_seconds = elapsed_clock_seconds(response_start,clock_end,clock_rate)

      if (native_ready) then
         call g_logger%info('[exchange_q]: native reference uses accepted reciprocal SCF k mesh, weights, EF, temperature, lattice and potential; order='// &
            trim(reciprocal_obj%kspace_ham_order),__FILE__,__LINE__)
         call native_exchange_q_driver(lattice_obj, reciprocal_obj, q_direct, native_contour_options, native_jq_ud, native_jq_du, &
            native_jq_sym, native_delta_j, native_curvature, native_report)
         call g_logger%info('[exchange_q]: native Turek contour points/poles='//trim(int2str(native_report%contour_points))//'/'// &
            trim(int2str(native_report%fermi_poles))//' solve/contour/pole='// &
            trim(real_to_string(native_report%solve_seconds))//'/'//trim(real_to_string(native_report%contour_seconds))//'/'// &
            trim(real_to_string(native_report%pole_seconds))//' s; native max ellipse='// &
            trim(real_to_string(native_report%native_max_ellipse_value)), __FILE__, __LINE__)
      end if
      call g_logger%info('[exchange_q]: timing endpoint/hamiltonian+T/C/spectral='//trim(real_to_string(endpoint_seconds))//'/'// &
         trim(real_to_string(hamiltonian_seconds))//'/'//trim(real_to_string(assembly_seconds))//'/'//trim(real_to_string(contraction_seconds))//' s', __FILE__, __LINE__)

      if (allocated(evals)) deallocate(evals)
      if (allocated(evecs)) deallocate(evecs)
      deallocate(weights, axes, hessian, torque_torque, contact, complete, endpoint_index)
   end subroutine compute_static_rotation_curvature

   subroutine exchange_q_convert_points(q_coordinates, q_list, lattice_obj, q_direct, q_cart)
      character(len=*), intent(in) :: q_coordinates
      real(rp), intent(in) :: q_list(:, :)
      type(lattice), intent(in) :: lattice_obj
      real(rp), allocatable, intent(out) :: q_direct(:, :), q_cart(:, :)
      real(rp) :: direct_to_cart(3,3)

      allocate(q_direct(3,size(q_list,2)), q_cart(3,size(q_list,2)))
      direct_to_cart = transpose(inverse_3x3(lattice_obj%a))
      if (trim(q_coordinates) == 'direct') then
         q_direct = q_list
         q_cart = matmul(direct_to_cart, q_direct)
      else
         q_cart = q_list
         q_direct = matmul(transpose(lattice_obj%a), q_cart)
      end if
   end subroutine exchange_q_convert_points

   module subroutine detect_commensurate_endpoint_map(k_points, k_weights, nk_mesh, q_point, endpoint_index, commensurate, residual)
      real(rp), intent(in) :: k_points(:, :), k_weights(:), q_point(3)
      integer, intent(in) :: nk_mesh(3)
      integer, intent(out) :: endpoint_index(:)
      logical, intent(out) :: commensurate
      real(rp), intent(out) :: residual
      integer, allocatable :: grid_index(:, :, :)
      integer :: nk, ix, iy, iz, ik, target, qcell(3), cell(3), flat
      integer :: n1, n2, n3, nint_coord
      real(rp) :: origin(3), folded(3), delta, qscaled
      real(rp), parameter :: mesh_tolerance = 2.0e-10_rp

      nk = size(k_points,2); n1 = nk_mesh(1); n2 = nk_mesh(2); n3 = nk_mesh(3)
      commensurate = .false.; residual = huge(1.0_rp)
      if (size(k_points,1) /= 3 .or. size(k_weights) /= nk .or. size(endpoint_index) /= nk .or. &
          n1 < 1 .or. n2 < 1 .or. n3 < 1 .or. nk /= n1*n2*n3) return

      ! Gamma is always a safe identity map, including a symmetry-reduced
      ! mesh.  Nonzero reuse below deliberately requires uniform weights and
      ! a full mesh.
      if (sqrt(sum(q_point**2)) <= q_tolerance) then
         do ik = 1, nk
            endpoint_index(ik) = ik
         end do
         commensurate = .true.; residual = 0.0_rp
         return
      end if
      if (nk > 1 .and. maxval(abs(k_weights-k_weights(1))) > mesh_tolerance*max(1.0_rp,abs(k_weights(1)))) return

      origin = fold_direct_point(k_points(:,1))
      allocate(grid_index(n1,n2,n3)); grid_index = 0
      do ik = 1, nk
         folded = fold_direct_point(k_points(:,ik))
         do ix = 1, 3
            delta = real(nk_mesh(ix),rp)*(folded(ix)-origin(ix))
            nint_coord = nint(delta)
            if (abs(delta-real(nint_coord,rp)) > mesh_tolerance) then
               deallocate(grid_index)
               return
            end if
            cell(ix) = modulo(nint_coord,nk_mesh(ix)) + 1
         end do
         flat = cell(1) + (cell(2)-1)*n1 + (cell(3)-1)*n1*n2
         ix = 1 + modulo(flat-1,n1)
         iy = 1 + modulo((flat-1)/n1,n2)
         iz = 1 + (flat-1)/(n1*n2)
         if (grid_index(ix,iy,iz) /= 0) then
            deallocate(grid_index)
            return
         end if
         grid_index(ix,iy,iz) = ik
      end do
      if (any(grid_index == 0)) then
         deallocate(grid_index)
         return
      end if

      residual = 0.0_rp
      do ix = 1, 3
         qscaled = q_point(ix)*real(nk_mesh(ix),rp)
         qcell(ix) = nint(qscaled)
         residual = max(residual,abs(qscaled-real(qcell(ix),rp))/real(nk_mesh(ix),rp))
      end do
      if (residual > mesh_tolerance) then
         deallocate(grid_index)
         return
      end if

      do ik = 1, nk
         folded = fold_direct_point(k_points(:,ik))
         do ix = 1, 3
            delta = real(nk_mesh(ix),rp)*(folded(ix)-origin(ix))
            cell(ix) = modulo(nint(delta)+qcell(ix),nk_mesh(ix)) + 1
         end do
         target = grid_index(cell(1),cell(2),cell(3))
         endpoint_index(ik) = target
         ! Check the actual folded endpoint, including the boundary crossing.
         folded = fold_direct_point(k_points(:,ik)+q_point)
         residual = max(residual,maxval(abs(folded-fold_direct_point(k_points(:,target)))))
         if (abs(k_weights(ik)-k_weights(target)) > mesh_tolerance*max(1.0_rp,abs(k_weights(ik)))) then
            deallocate(grid_index)
            return
         end if
      end do
      commensurate = residual <= mesh_tolerance
      deallocate(grid_index)
   end subroutine detect_commensurate_endpoint_map

   pure function fold_direct_point(point) result(folded)
      real(rp), intent(in) :: point(3)
      real(rp) :: folded(3)
      folded = point-floor(point+0.5_rp)
   end function fold_direct_point

   subroutine validate_native_turek_capability(native_crosscheck, native_turek, ham)
      logical, intent(in) :: native_crosscheck, native_turek
      type(hamiltonian), intent(in) :: ham
      if ((native_crosscheck .or. native_turek) .and. ham%charge%lattice%nrec/=1) then
         call g_logger%fatal('[linear_response]: native contour production is currently capability-gated to one sublattice', &
            __FILE__,__LINE__)
      end if
   end subroutine validate_native_turek_capability

   subroutine validate_finite_h_capability(ham,recip)
      type(hamiltonian), intent(in) :: ham
      type(reciprocal), intent(in) :: recip
      if (trim(recip%kspace_ham_order)=='second') then
         if (.not.ham%hoh .or. .not.allocated(ham%eeo) .or. .not.allocated(ham%enim)) then
            call g_logger%fatal('[linear_response]: SECOND_ORDER_STATE_NOT_ACTIVE in reciprocal finite-H path',__FILE__,__LINE__)
         end if
      end if
   end subroutine validate_finite_h_capability

   integer function fixture_basis_size(ham) result(size_orb)
      type(hamiltonian), intent(in) :: ham
      if (.not. associated(ham%charge)) then
         size_orb = 0
      else
         size_orb = size(ham%charge%lattice%sbar,1)
      end if
   end function fixture_basis_size

   module subroutine validate_rotation_capability(n_q_points, native_crosscheck, native_turek, control_obj, ham, self_obj, recip)
      integer, intent(in) :: n_q_points
      logical, intent(in) :: native_crosscheck, native_turek
      type(control), intent(in) :: control_obj
      type(hamiltonian), intent(in) :: ham
      type(self), intent(in) :: self_obj
      type(reciprocal), intent(in) :: recip
      logical :: native_requested

      native_requested = native_crosscheck .or. native_turek
      if (n_q_points < 2) call g_logger%fatal('[linear_response]: q path is empty or has no finite-q point',__FILE__,__LINE__)
      if (trim(control_obj%calctype)/='B' .or. control_obj%nsp/=1 .or. control_obj%has_soc()) then
         call g_logger%fatal('[linear_response]: capability gate requires bulk nsp=1 scalar-relativistic collinear state with SOC off', &
            __FILE__,__LINE__)
      end if
      if (.not.self_obj%use_kspace .or. .not.allocated(self_obj%reciprocal_scf_cache)) then
         call g_logger%fatal('[linear_response]: requires self%use_kspace=.true. to consume the accepted reciprocal SCF state', &
            __FILE__,__LINE__)
      end if
      if (ham%ccor_2c .or. ham%hubbard_u_general_check .or. ham%hubbard_v_check .or. control_obj%constraints_enable) then
         call g_logger%fatal('[linear_response]: unsupported CCOR/Hubbard/constraint combination',__FILE__,__LINE__)
      end if
      if (trim(recip%reciprocal_mode)/='ham_only') then
         call g_logger%fatal('[linear_response]: capability gate requires reciprocal_mode=ham_only',__FILE__,__LINE__)
      end if
      if (trim(recip%kspace_ham_order)/='first' .and. trim(recip%kspace_ham_order)/='second') then
         call g_logger%fatal('[linear_response]: reciprocal Hamiltonian order must be first or second',__FILE__,__LINE__)
      end if
      if (fixture_basis_size(ham)/=9) then
         call g_logger%fatal('[linear_response]: capability gate requires the full spd production basis',__FILE__,__LINE__)
      end if
      if (native_requested) call validate_native_turek_capability(native_crosscheck, native_turek, ham)
      call validate_finite_h_capability(ham,recip)
   end subroutine validate_rotation_capability

   function spread_q(q, count) result(points)
      real(rp), intent(in) :: q(3)
      integer, intent(in) :: count
      real(rp) :: points(3,count)
      integer :: i
      do i = 1, count
         points(:,i) = q
      end do
   end function spread_q

   subroutine solve_unfolded_endpoints(recip, points, eigenvalues, eigenvectors)
      type(reciprocal), intent(inout) :: recip
      real(rp), intent(in) :: points(:, :)
      real(rp), allocatable, intent(out) :: eigenvalues(:, :)
      complex(rp), allocatable, intent(out) :: eigenvectors(:, :, :)
      complex(rp), allocatable :: h(:, :), work(:)
      real(rp), allocatable :: rwork(:)
      integer :: nmat, nk, ik, lwork, info
      external :: zheev

      nk = size(points,2)
      nmat = size(recip%eigenvalues,1)
      allocate(eigenvalues(nmat,nk), eigenvectors(nmat,nmat,nk), h(nmat,nmat), rwork(max(1,3*nmat-2)))
      lwork = max(1, 2*nmat-1)
      allocate(work(lwork))
      do ik = 1, nk
         call recip%build_hamiltonian_at_kpoint(points(:,ik), h)
         call zheev('V','U',nmat,h,nmat,eigenvalues(:,ik),work,lwork,rwork,info)
         if (info /= 0) call g_logger%fatal('[exchange_q]: endpoint diagonalization failed, info='//int2str(info), &
            __FILE__, __LINE__)
         eigenvectors(:,:,ik) = h
      end do
      deallocate(h,work,rwork)
   end subroutine solve_unfolded_endpoints


   subroutine native_exchange_q_driver(lat, recip, q_direct, options, jq_ud_out, jq_du_out, jq_sym_out, delta_j, curvature, report)
      type(lattice), intent(inout) :: lat
      type(reciprocal), intent(in) :: recip
      real(rp), intent(in) :: q_direct(:, :)
      type(native_turek_contour_options), intent(in) :: options
      real(rp), intent(out) :: jq_ud_out(:), jq_du_out(:), jq_sym_out(:), delta_j(:), curvature(:)
      type(native_turek_contour_report), intent(out) :: report
      complex(rp), allocatable :: jq_ud(:, :, :), jq_du(:, :, :), jq_sym(:, :, :), native_delta(:, :, :), native_curvature(:, :, :)
      real(rp) :: kT
      integer :: nsite, nq

      nsite=lat%nrec; nq=size(q_direct,2)
      if (nsite/=1) error stop 'native_exchange_q_driver: production scalar reduction requires nsite=1'
      kT=max(recip%temperature*kB_ry_per_k,1.0e-10_rp)
      allocate(jq_ud(nsite,nsite,nq),jq_du(nsite,nsite,nq),jq_sym(nsite,nsite,nq), &
         native_delta(nsite,nsite,nq),native_curvature(nsite,nsite,nq))
      call native_turek_static_reference(lat,recip%k_points,recip%k_weights,q_direct,recip%fermi_level,kT,options, &
         jq_ud,jq_du,jq_sym,native_delta,native_curvature,report)
      jq_ud_out=real(jq_ud(1,1,:),rp); jq_du_out=real(jq_du(1,1,:),rp)
      jq_sym_out=real(jq_sym(1,1,:),rp); delta_j=real(native_delta(1,1,:),rp); curvature=real(native_curvature(1,1,:),rp)
      deallocate(jq_ud,jq_du,jq_sym,native_delta,native_curvature)
   end subroutine native_exchange_q_driver


   pure function elapsed_clock_seconds(start_count, end_count, rate) result(seconds)
      integer, intent(in) :: start_count, end_count, rate
      real(rp) :: seconds
      seconds = real(end_count-start_count,rp)/real(max(rate,1),rp)
   end function elapsed_clock_seconds

   pure function real_to_string(value) result(text)
      real(rp), intent(in) :: value
      character(len=48) :: text
      write(text,'(es24.16)') value
   end function real_to_string

   module subroutine lr_config_restore_to_default(this)
      class(linear_response_config), intent(out) :: this
      if (allocated(this%q_list)) deallocate(this%q_list)
      if (allocated(this%omega_grid)) deallocate(this%omega_grid)
      if (allocated(this%eta_grid)) deallocate(this%eta_grid)
      this%fname = ''
      this%formulation = 'rotation'
      this%representation = 'radial_points'
      this%bare_response = 'lehmann'
      this%realspace_solver = 'auto'
      this%interaction = 'alsda'
      this%projection = 'spd'
      this%diagnostics = 'none'
      this%channel = 'chi_plus'
      this%q_coordinates = 'direct'
      this%q_file = ''
      this%n_q = 1
      this%n_q_points = 0
      this%n_omega = 1
      this%use_omega_grid = .false.
      this%omega_min = 0.0_rp
      this%omega_max = 0.0_rp
      this%eta = 0.01_rp
      this%n_eta = 1
      this%response_lmax = -1
      this%write_full_matrix = .true.
      this%rotation_axis = [1.0_rp, 0.0_rp, 0.0_rp]
      this%finite_h_spectral_mode = 'metallic'
      this%finite_h_response_backend = 'spectral'
      this%contour_points = 32
      this%contour_shape = 'ellipse'
      this%contour_margin = 0.25_rp
      this%contour_height_fraction = 0.35_rp
      this%contour_account_fermi_poles = .true.
      this%native_turek = .false.
      this%native_crosscheck = .false.
      this%native_green_eta = 1.0e-3_rp
      this%native_energy_points = 0
      this%native_contour_points = 64
      this%native_contour_margin = 0.25_rp
      this%native_contour_height_fraction = 0.35_rp
      this%native_contour_account_fermi_poles = .true.
      this%native_contour_target_fermi_poles = 0
      this%native_rsgf_provider = 'auto'
      this%gf_integration_points = 2001
      this%gf_integration_eta = 0.0_rp
      this%gf_energy_margin = 1.0_rp
      this%output_file = 'rotation_dynamics.dat'
      this%rotation_eta_ladder = [1.0e-4_rp, 2.5e-5_rp, 6.25e-6_rp]
      this%rotation_probe_omega = 1.0e-5_rp
      this%rotation_probe_eta = 1.0e-9_rp
      this%rotation_slope_step = 1.0e-5_rp
      this%rotation_pole_window_floor = 1.0e-3_rp
      this%rotation_pole_window_scale = 2.5_rp
      this%rotation_pole_window_max = 2.0e-2_rp
      this%rotation_pole_coarse_points = 61
      this%rotation_pole_fine_points = 41
      this%rotation_pole_refinement_half_width = 2.0_rp
      this%rotation_n_omega = 0
      this%rotation_omega_min = 0.0_rp
      this%rotation_omega_max = 2.0e-2_rp
      this%rotation_eta = 1.0e-4_rp
      this%rotation_grid_file = 'rotation_response_grid.dat'
   end subroutine lr_config_restore_to_default

   module subroutine lr_restore_to_default(this)
      class(linear_response), intent(out) :: this
      call this%config%restore_to_default()
   end subroutine lr_restore_to_default

   module subroutine lr_load_config(this, filename, validate_request)
      class(linear_response), intent(inout) :: this
      character(len=*), intent(in) :: filename
      logical, intent(in), optional :: validate_request
      logical :: validate
      integer :: iostatus, funit, n_q_points, n_q, n_omega, n_eta, iw, iq, scan_unit, scan_status
      character(len=16) :: formulation, representation, bare_response, realspace_solver, interaction, diagnostics, q_coordinates
      character(len=8) :: projection
      character(len=32) :: channel, native_rsgf_provider
      character(len=256) :: q_file, output_file, scan_line
      logical :: native_turek, native_crosscheck, native_contour_account_fermi_poles, contour_account_fermi_poles
      character(len=24) :: finite_h_spectral_mode
      character(len=16) :: finite_h_response_backend, contour_shape
      real(rp) :: rotation_axis(3), native_green_eta
      real(rp) :: contour_margin, contour_height_fraction, native_contour_margin, native_contour_height_fraction
      integer :: native_energy_points, contour_points, native_contour_points, native_contour_target_fermi_poles
      real(rp) :: rotation_eta_ladder(3), rotation_probe_omega, rotation_probe_eta, rotation_slope_step
      real(rp) :: rotation_pole_window_floor, rotation_pole_window_scale, rotation_pole_window_max, rotation_pole_refinement_half_width
      integer :: rotation_pole_coarse_points, rotation_pole_fine_points
      integer :: rotation_n_omega
      real(rp) :: rotation_omega_min, rotation_omega_max, rotation_eta
      character(len=sl) :: rotation_grid_file
      real(rp) :: q_list(3, linear_response_max_points)
      real(rp) :: omega_grid(4096), eta_grid(16), omega_min, omega_max, eta
      logical :: use_omega_grid, write_full_matrix
      integer :: response_lmax, gf_integration_points
      real(rp) :: gf_integration_eta, gf_energy_margin
      namelist /linear_response/ formulation, q_coordinates, q_file, n_q_points, q_list, rotation_axis, &
         representation, bare_response, realspace_solver, interaction, projection, diagnostics, channel, n_q, n_omega, &
         use_omega_grid, omega_grid, omega_min, omega_max, eta, n_eta, eta_grid, response_lmax, write_full_matrix, &
         native_rsgf_provider, gf_integration_points, gf_integration_eta, gf_energy_margin, &
         finite_h_spectral_mode, finite_h_response_backend, contour_points, contour_shape, contour_margin, &
         contour_height_fraction, contour_account_fermi_poles, native_turek, native_crosscheck, native_green_eta, &
         native_energy_points, native_contour_points, native_contour_margin, native_contour_height_fraction, &
         native_contour_account_fermi_poles, native_contour_target_fermi_poles, output_file, rotation_eta_ladder, &
         rotation_probe_omega, rotation_probe_eta, rotation_slope_step, rotation_pole_window_floor, &
         rotation_pole_window_scale, rotation_pole_window_max, rotation_pole_coarse_points, rotation_pole_fine_points, &
         rotation_pole_refinement_half_width, rotation_n_omega, rotation_omega_min, rotation_omega_max, rotation_eta, &
         rotation_grid_file

      call this%config%restore_to_default()
      this%config%fname = filename
      formulation = this%config%formulation
      representation = this%config%representation
      bare_response = this%config%bare_response
      realspace_solver = this%config%realspace_solver
      interaction = this%config%interaction
      projection = this%config%projection
      diagnostics = this%config%diagnostics
      channel = this%config%channel
      q_coordinates = this%config%q_coordinates
      q_file = this%config%q_file
      n_q = this%config%n_q
      n_q_points = 0
      q_list = 0.0_rp
      n_omega = this%config%n_omega
      use_omega_grid = this%config%use_omega_grid
      omega_grid = 0.0_rp
      omega_min = this%config%omega_min
      omega_max = this%config%omega_max
      eta = this%config%eta
      n_eta = this%config%n_eta
      eta_grid = 0.0_rp
      eta_grid(1) = eta
      response_lmax = this%config%response_lmax
      write_full_matrix = this%config%write_full_matrix
      native_rsgf_provider = this%config%native_rsgf_provider
      gf_integration_points = this%config%gf_integration_points
      gf_integration_eta = this%config%gf_integration_eta
      gf_energy_margin = this%config%gf_energy_margin
      rotation_axis = this%config%rotation_axis
      finite_h_spectral_mode = this%config%finite_h_spectral_mode
      finite_h_response_backend = this%config%finite_h_response_backend
      contour_points = this%config%contour_points
      contour_shape = this%config%contour_shape
      contour_margin = this%config%contour_margin
      contour_height_fraction = this%config%contour_height_fraction
      contour_account_fermi_poles = this%config%contour_account_fermi_poles
      native_turek = this%config%native_turek
      native_crosscheck = this%config%native_crosscheck
      native_green_eta = this%config%native_green_eta
      native_energy_points = this%config%native_energy_points
      native_contour_points = this%config%native_contour_points
      native_contour_margin = this%config%native_contour_margin
      native_contour_height_fraction = this%config%native_contour_height_fraction
      native_contour_account_fermi_poles = this%config%native_contour_account_fermi_poles
      native_contour_target_fermi_poles = this%config%native_contour_target_fermi_poles
      output_file = this%config%output_file
      rotation_eta_ladder = this%config%rotation_eta_ladder
      rotation_probe_omega = this%config%rotation_probe_omega
      rotation_probe_eta = this%config%rotation_probe_eta
      rotation_slope_step = this%config%rotation_slope_step
      rotation_pole_window_floor = this%config%rotation_pole_window_floor
      rotation_pole_window_scale = this%config%rotation_pole_window_scale
      rotation_pole_window_max = this%config%rotation_pole_window_max
      rotation_pole_coarse_points = this%config%rotation_pole_coarse_points
      rotation_pole_fine_points = this%config%rotation_pole_fine_points
      rotation_pole_refinement_half_width = this%config%rotation_pole_refinement_half_width
      rotation_n_omega = this%config%rotation_n_omega
      rotation_omega_min = this%config%rotation_omega_min
      rotation_omega_max = this%config%rotation_omega_max
      rotation_eta = this%config%rotation_eta
      rotation_grid_file = this%config%rotation_grid_file

      open(newunit=scan_unit, file=filename, action='read', status='old', iostat=scan_status)
      if (scan_status == 0) then
         do
            read(scan_unit, '(A)', iostat=scan_status) scan_line
            if (scan_status /= 0) exit
            if (index(adjustl(lower(scan_line)), '&'//'tddft') == 1) then
               close(scan_unit)
               call g_logger%fatal('&'//'tddft was replaced by &linear_response (see docs/linear_response/FORMULATION.md)', __FILE__, __LINE__)
            end if
         end do
         close(scan_unit)
      end if

      open(newunit=funit, file=filename, action='read', status='old', iostat=iostatus)
      if (iostatus /= 0) call g_logger%fatal('[linear_response]: input file '//trim(filename)//' not found', __FILE__, __LINE__)
      read(funit, nml=linear_response, iostat=iostatus)
      close(funit)
      if (iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
         call g_logger%fatal('[linear_response]: error reading &linear_response, iostatus='//int2str(iostatus), __FILE__, __LINE__)
      end if

      this%config%formulation = lower(trim(formulation))
      this%config%representation = lower(trim(representation))
      this%config%bare_response = lower(trim(bare_response))
      this%config%realspace_solver = lower(trim(realspace_solver))
      this%config%interaction = lower(trim(interaction))
      this%config%projection = lower(trim(projection))
      this%config%diagnostics = lower(trim(diagnostics))
      this%config%channel = lower(trim(channel))
      this%config%q_coordinates = lower(trim(q_coordinates))
      this%config%q_file = trim(q_file)
      this%config%n_q = n_q
      this%config%n_omega = n_omega
      this%config%use_omega_grid = use_omega_grid
      this%config%omega_min = omega_min
      this%config%omega_max = omega_max
      this%config%eta = eta
      this%config%n_eta = n_eta
      this%config%response_lmax = response_lmax
      this%config%write_full_matrix = write_full_matrix
      this%config%native_rsgf_provider = lower(trim(native_rsgf_provider))
      this%config%gf_integration_points = gf_integration_points
      this%config%gf_integration_eta = gf_integration_eta
      this%config%gf_energy_margin = gf_energy_margin
      this%config%rotation_axis = rotation_axis
      this%config%finite_h_spectral_mode = lower(trim(finite_h_spectral_mode))
      this%config%finite_h_response_backend = lower(trim(finite_h_response_backend))
      this%config%contour_points = contour_points
      this%config%contour_shape = lower(trim(contour_shape))
      this%config%contour_margin = contour_margin
      this%config%contour_height_fraction = contour_height_fraction
      this%config%contour_account_fermi_poles = contour_account_fermi_poles
      this%config%native_turek = native_turek
      this%config%native_crosscheck = native_crosscheck
      this%config%native_green_eta = native_green_eta
      this%config%native_energy_points = native_energy_points
      this%config%native_contour_points = native_contour_points
      this%config%native_contour_margin = native_contour_margin
      this%config%native_contour_height_fraction = native_contour_height_fraction
      this%config%native_contour_account_fermi_poles = native_contour_account_fermi_poles
      this%config%native_contour_target_fermi_poles = native_contour_target_fermi_poles
      this%config%output_file = trim(output_file)
      if (this%config%formulation /= 'rotation' .and. trim(output_file) == 'rotation_dynamics.dat') then
         this%config%output_file = 'tddft_response.dat'
      end if
      this%config%rotation_eta_ladder = rotation_eta_ladder
      this%config%rotation_probe_omega = rotation_probe_omega
      this%config%rotation_probe_eta = rotation_probe_eta
      this%config%rotation_slope_step = rotation_slope_step
      this%config%rotation_pole_window_floor = rotation_pole_window_floor
      this%config%rotation_pole_window_scale = rotation_pole_window_scale
      this%config%rotation_pole_window_max = rotation_pole_window_max
      this%config%rotation_pole_coarse_points = rotation_pole_coarse_points
      this%config%rotation_pole_fine_points = rotation_pole_fine_points
      this%config%rotation_pole_refinement_half_width = rotation_pole_refinement_half_width
      this%config%rotation_n_omega = rotation_n_omega
      this%config%rotation_omega_min = rotation_omega_min
      this%config%rotation_omega_max = rotation_omega_max
      this%config%rotation_eta = rotation_eta
      this%config%rotation_grid_file = trim(rotation_grid_file)

      if (this%config%q_coordinates /= 'direct' .and. this%config%q_coordinates /= 'cartesian') then
         call g_logger%fatal("[linear_response]: q_coordinates must be 'direct' or 'cartesian'", __FILE__, __LINE__)
      end if

      validate = .false.
      if (present(validate_request)) validate = validate_request
      if (len_trim(this%config%q_file) > 0) then
         call read_linear_response_q_file(this, this%config%q_file)
         this%config%n_q = this%config%n_q_points
      else
         if (this%config%formulation == 'rotation') then
            if (n_q_points < 1 .or. n_q_points > linear_response_max_points) then
               if (validate) call g_logger%fatal('[linear_response]: n_q_points must be in [1,'// &
                  int2str(linear_response_max_points)//']', __FILE__, __LINE__)
            else
               this%config%n_q_points = n_q_points
               this%config%n_q = n_q_points
               allocate(this%config%q_list(3, n_q_points))
               this%config%q_list = q_list(:, 1:n_q_points)
            end if
         else if (n_q < 1 .or. n_q > linear_response_max_points) then
            if (validate) call g_logger%fatal('[linear_response]: n_q_points must be in [1,'// &
               int2str(linear_response_max_points)//']', __FILE__, __LINE__)
         else
            this%config%n_q_points = n_q
            allocate(this%config%q_list(3, n_q))
            this%config%q_list = q_list(:, 1:n_q)
         end if
      end if
      if (.not. validate) return

      if (this%config%formulation /= 'rotation') then
         if (n_omega < 1 .or. n_omega > 4096) then
            call g_logger%fatal('[linear_response]: n_omega must be in [1,4096]', __FILE__, __LINE__)
         end if
         if (n_eta < 1 .or. n_eta > 16) then
            call g_logger%fatal('[linear_response]: n_eta must be in [1,16]', __FILE__, __LINE__)
         end if
         allocate(this%config%omega_grid(n_omega), this%config%eta_grid(n_eta))
         if (use_omega_grid) then
            this%config%omega_grid = omega_grid(1:n_omega)
         else if (n_omega == 1) then
            this%config%omega_grid(1) = omega_min
         else
            do iw = 1, n_omega
               this%config%omega_grid(iw) = omega_min + real(iw-1,rp)*(omega_max-omega_min)/real(n_omega-1,rp)
            end do
         end if
         if (n_eta == 1) then
            this%config%eta_grid(1) = eta
         else
            this%config%eta_grid = eta_grid(1:n_eta)
            this%config%eta = this%config%eta_grid(1)
         end if
         call validate_linear_response_config(this%config)
         return
      end if

      if (len_trim(this%config%output_file) == 0) call g_logger%fatal('[linear_response]: output_file must not be blank', __FILE__, __LINE__)
      if (any(this%config%rotation_eta_ladder <= 0.0_rp)) then
         call g_logger%fatal('[linear_response]: rotation_eta_ladder values must be positive', __FILE__, __LINE__)
      end if
      if (this%config%rotation_probe_omega <= 0.0_rp .or. this%config%rotation_probe_eta <= 0.0_rp .or. &
          this%config%rotation_slope_step <= 0.0_rp) then
         call g_logger%fatal('[linear_response]: rotation probe and slope controls must be positive', __FILE__, __LINE__)
      end if
      if (abs(this%config%rotation_slope_step-this%config%rotation_probe_omega) > &
          1.0e-12_rp*max(1.0_rp,abs(this%config%rotation_probe_omega))) then
         call g_logger%info('[linear_response]: rotation_slope_step is deprecated and ignored; '// &
            'rotation_probe_omega is the canonical finite-difference step', __FILE__, __LINE__)
      end if
      if (this%config%rotation_pole_window_floor <= 0.0_rp .or. this%config%rotation_pole_window_scale <= 0.0_rp .or. &
          this%config%rotation_pole_window_max <= 0.0_rp .or. &
          this%config%rotation_pole_window_max < this%config%rotation_pole_window_floor) then
         call g_logger%fatal('[linear_response]: rotation pole-window controls are invalid', __FILE__, __LINE__)
      end if
      if (this%config%rotation_pole_coarse_points < 2 .or. this%config%rotation_pole_fine_points < 2 .or. &
          this%config%rotation_pole_refinement_half_width <= 0.0_rp) then
         call g_logger%fatal('[linear_response]: rotation pole-scan grid controls are invalid', __FILE__, __LINE__)
      end if
      if (this%config%rotation_n_omega < 0) call g_logger%fatal('[linear_response]: rotation_n_omega must be non-negative', __FILE__, __LINE__)
      if (this%config%rotation_n_omega > 0) then
         if (this%config%rotation_omega_max < this%config%rotation_omega_min .or. this%config%rotation_eta <= 0.0_rp) then
            call g_logger%fatal('[linear_response]: rotation frequency-grid bounds or eta are invalid', __FILE__, __LINE__)
         end if
         if (len_trim(this%config%rotation_grid_file) == 0) then
            call g_logger%fatal('[linear_response]: rotation_grid_file must not be blank when rotation_n_omega is enabled', __FILE__, __LINE__)
         end if
      end if
      if (this%config%native_green_eta <= 0.0_rp) call g_logger%fatal('[linear_response]: native_green_eta must be positive', __FILE__, __LINE__)
      if (this%config%native_energy_points < 0) call g_logger%fatal('[linear_response]: native_energy_points must be non-negative', __FILE__, __LINE__)
      if (this%config%native_contour_points < 8) call g_logger%fatal('[linear_response]: native_contour_points must be at least 8', __FILE__, __LINE__)
      if (this%config%native_contour_margin <= 0.0_rp .or. this%config%native_contour_height_fraction <= 0.0_rp) then
         call g_logger%fatal('[linear_response]: native contour margin and height fraction must be positive', __FILE__, __LINE__)
      end if
      if (this%config%native_contour_target_fermi_poles < 0 .or. &
          (this%config%native_contour_target_fermi_poles > 0 .and. mod(this%config%native_contour_target_fermi_poles,2) /= 0)) then
         call g_logger%fatal('[linear_response]: native_contour_target_fermi_poles must be zero or positive even', __FILE__, __LINE__)
      end if
      if (maxval(abs(this%config%rotation_axis)) <= tiny(1.0_rp)) call g_logger%fatal('[linear_response]: rotation_axis must be nonzero', __FILE__, __LINE__)
      if (this%config%finite_h_spectral_mode /= 'metallic' .and. this%config%finite_h_spectral_mode /= 'legacy_occupied') then
         call g_logger%fatal("[linear_response]: finite_h_spectral_mode must be 'metallic' or 'legacy_occupied'", __FILE__, __LINE__)
      end if
      if (this%config%finite_h_response_backend /= 'spectral' .and. this%config%finite_h_response_backend /= 'contour' .and. &
          this%config%finite_h_response_backend /= 'both') then
         call g_logger%fatal("[linear_response]: finite_h_response_backend must be 'spectral', 'contour', or 'both'", __FILE__, __LINE__)
      end if
      if (this%config%contour_points < 8) call g_logger%fatal('[linear_response]: contour_points must be at least 8', __FILE__, __LINE__)
      if (this%config%contour_shape /= 'ellipse') call g_logger%fatal("[linear_response]: contour_shape must be 'ellipse'", __FILE__, __LINE__)
      if (this%config%contour_margin <= 0.0_rp .or. this%config%contour_height_fraction <= 0.0_rp) then
         call g_logger%fatal('[linear_response]: contour margin and height fraction must be positive', __FILE__, __LINE__)
      end if
      if (this%config%n_q_points < 2) call g_logger%fatal('[linear_response]: provide at least q=0 and one finite-q point', __FILE__, __LINE__)
      if (any(.not. ieee_is_finite(this%config%q_list))) call g_logger%fatal('[linear_response]: q path contains a non-finite coordinate', __FILE__, __LINE__)
   end subroutine lr_load_config

   subroutine validate_linear_response_config(config)
      type(linear_response_config), intent(in) :: config
      character(len=96) :: compatibility_key
      integer :: iq, jq
      logical :: gamma_found, covariance_found

      if (trim(config%diagnostics) /= 'none' .and. trim(config%diagnostics) /= 'invariants') then
         call g_logger%fatal("[linear_response]: diagnostics must be 'none' or 'invariants'", __FILE__, __LINE__)
      end if
      if (trim(config%channel) /= 'chi_plus' .and. trim(config%channel) /= 'chi_minus') then
         call g_logger%fatal('[linear_response]: channel must be chi_plus or chi_minus', __FILE__, __LINE__)
      end if
      if (trim(config%bare_response) == 'kspace_resolvent') then
         call g_logger%fatal("[linear_response]: bare_response='kspace_resolvent' is reserved; not implemented", __FILE__, __LINE__)
      end if
      if (trim(config%bare_response) /= 'lehmann' .and. trim(config%bare_response) /= 'realspace_gf') then
         call g_logger%fatal('[linear_response]: unsupported bare_response', __FILE__, __LINE__)
      end if
      if (trim(config%realspace_solver) /= 'auto' .and. trim(config%realspace_solver) /= 'block' .and. &
          trim(config%realspace_solver) /= 'chebyshev') then
         call g_logger%fatal('[linear_response]: realspace_solver must be auto, block or chebyshev', __FILE__, __LINE__)
      end if
      if (trim(config%realspace_solver) /= 'auto' .and. trim(config%bare_response) /= 'realspace_gf') then
         call g_logger%fatal('[linear_response]: realspace_solver is only valid with bare_response=realspace_gf', __FILE__, __LINE__)
      end if
      if (config%eta <= 0.0_rp .or. .not. allocated(config%eta_grid) .or. any(config%eta_grid <= 0.0_rp)) then
         call g_logger%fatal('[linear_response]: eta and every eta_grid value must be positive', __FILE__, __LINE__)
      end if
      if (config%response_lmax < -1 .or. config%response_lmax > 4) then
         call g_logger%fatal('[linear_response]: response_lmax must be -1 through 4', __FILE__, __LINE__)
      end if
      if (.not. allocated(config%q_list) .or. size(config%q_list,2) < 1 .or. any(.not. ieee_is_finite(config%q_list))) then
         call g_logger%fatal('[linear_response]: q_list must contain finite q points', __FILE__, __LINE__)
      end if
      if (.not. allocated(config%omega_grid) .or. size(config%omega_grid) < 1 .or. &
          any(.not. ieee_is_finite(config%omega_grid))) then
         call g_logger%fatal('[linear_response]: omega grid must contain finite values', __FILE__, __LINE__)
      end if
      if (len_trim(config%output_file) == 0) call g_logger%fatal('[linear_response]: output_file must not be blank', __FILE__, __LINE__)

      select case (trim(config%formulation))
      case ('tddft')
         compatibility_key = trim(config%representation)//':'//trim(config%bare_response)//':'//trim(config%interaction)
         select case (compatibility_key)
         case ('radial_points:lehmann:alsda', 'radial_points:lehmann:lcmm', &
               'radial_points:realspace_gf:alsda', 'radial_points:realspace_gf:lcmm', &
               'product_compact:lehmann:alsda', 'product_compact:lehmann:none')
            continue
         case default
            call g_logger%fatal('[linear_response]: formulation/representation/bare_response/interaction is not an allowed compatibility row', &
               __FILE__, __LINE__)
         end select
      case ('projected')
         if (trim(config%bare_response) /= 'lehmann' .or. trim(config%representation) /= 'radial_points') then
            call g_logger%fatal('[linear_response]: projected requires radial_points and lehmann', __FILE__, __LINE__)
         end if
         select case (trim(config%interaction))
         case ('none', 'stoner_fit', 'lcmm')
            continue
         case default
            call g_logger%fatal('[linear_response]: projected interaction is not an allowed compatibility row', __FILE__, __LINE__)
         end select
         if (trim(config%projection) /= 'd' .and. trim(config%projection) /= 'spd' .and. trim(config%projection) /= 'both') then
            call g_logger%fatal('[linear_response]: projection must be d, spd or both', __FILE__, __LINE__)
         end if
      case default
         call g_logger%fatal("[linear_response]: formulation must be 'rotation', 'tddft' or 'projected'", __FILE__, __LINE__)
      end select
      if (trim(config%bare_response) == 'realspace_gf' .and. config%gf_integration_points < 3) then
         call g_logger%fatal('[linear_response]: realspace_gf requires gf_integration_points >= 3', __FILE__, __LINE__)
      end if
      if (trim(config%formulation) == 'tddft' .and. trim(config%representation) == 'product_compact' .and. &
          trim(config%interaction) == 'alsda' .and. trim(config%diagnostics) == 'invariants') then
         gamma_found = any(sum(abs(config%q_list), dim=1) <= 1.0e-12_rp)
         if (.not. gamma_found) call g_logger%fatal('[linear_response]: diagnostics=invariants requires Gamma for compact Dyson', __FILE__, __LINE__)
         covariance_found = .false.
         do iq = 1, size(config%q_list, 2)
            if (sum(abs(config%q_list(:,iq))) <= 1.0e-12_rp) cycle
            do jq = 1, size(config%q_list, 2)
               if (maxval(abs(config%q_list(:,jq) + config%q_list(:,iq))) <= 2.0e-11_rp) then
                  covariance_found = .true.
                  exit
               end if
            end do
            if (covariance_found) exit
         end do
         if (.not. covariance_found) call g_logger%fatal('[linear_response]: diagnostics=invariants requires an exact +q/-q pair', &
            __FILE__, __LINE__)
      end if
   end subroutine validate_linear_response_config

   subroutine read_linear_response_q_file(this, filename)
      class(linear_response), intent(inout) :: this
      character(len=*), intent(in) :: filename
      integer :: funit, ios, iq, n_q

      open(newunit=funit, file=trim(filename), action='read', status='old', iostat=ios)
      if (ios /= 0) call g_logger%fatal('[linear_response]: q_file '//trim(filename)//' not found', __FILE__, __LINE__)
      read(funit, *, iostat=ios) n_q
      if (ios /= 0 .or. n_q < 1 .or. n_q > linear_response_max_points) then
         close(funit)
         call g_logger%fatal('[linear_response]: malformed q_file count in '//trim(filename), __FILE__, __LINE__)
      end if
      allocate(this%config%q_list(3,n_q))
      this%config%n_q_points = n_q
      do iq = 1, n_q
         read(funit, *, iostat=ios) this%config%q_list(:,iq)
         if (ios /= 0) then
            close(funit)
            call g_logger%fatal('[linear_response]: malformed q_file row '//int2str(iq)//' in '//trim(filename), __FILE__, __LINE__)
         end if
      end do
      close(funit)
   end subroutine read_linear_response_q_file

   module subroutine lr_validate_capability(this, control_obj, hamiltonian_obj, self_obj, reciprocal_obj)
      class(linear_response), intent(in) :: this
      type(control), intent(in) :: control_obj
      type(hamiltonian), intent(in) :: hamiltonian_obj
      type(self), intent(in) :: self_obj
      type(reciprocal), intent(in) :: reciprocal_obj
      call validate_rotation_capability(this%config%n_q_points, this%config%native_crosscheck, this%config%native_turek, &
         control_obj, hamiltonian_obj, self_obj, reciprocal_obj)
   end subroutine lr_validate_capability

   module subroutine lr_run_rotation(this, control_obj, lattice_obj, hamiltonian_obj, energy_obj, self_obj, reciprocal_obj)
      class(linear_response), intent(in) :: this
      type(control), intent(in) :: control_obj
      type(lattice), intent(inout) :: lattice_obj
      type(hamiltonian), intent(inout) :: hamiltonian_obj
      type(energy), intent(in) :: energy_obj
      type(self), intent(inout) :: self_obj
      type(reciprocal), intent(inout) :: reciprocal_obj
      type(lmto_live_hamiltonian_fixture) :: fixture
      real(rp), allocatable :: q_direct(:, :), q_cart(:, :)
      real(rp), allocatable :: finite_total(:, :, :), finite_tt(:, :, :), finite_contact(:, :, :)
      real(rp), allocatable :: spectral_total(:, :, :), spectral_tt(:, :, :), spectral_contact(:, :, :)
      real(rp), allocatable :: contour_total(:, :, :), contour_tt(:, :, :), contour_contact(:, :, :)
      real(rp), allocatable :: native_jq_ud(:), native_jq_du(:), native_jq_sym(:), native_delta_j(:), native_curvature(:)
      logical, allocatable :: endpoint_reused(:), q_commensurate(:)
      real(rp), allocatable :: endpoint_residual(:)
      character(len=24), allocatable :: endpoint_mode(:)
      type(native_turek_contour_report) :: native_report
      real(rp) :: endpoint_seconds, assembly_seconds, contraction_seconds, hamiltonian_seconds
      real(rp) :: gf_seconds, solve_seconds, contour_seconds, total_response_seconds
      logical :: native_ready

      call this%validate_capability(control_obj, hamiltonian_obj, self_obj, reciprocal_obj)
      call compute_static_rotation_curvature(this%config%q_coordinates, this%config%q_list, this%config%rotation_axis, &
         this%config%finite_h_spectral_mode, this%config%finite_h_response_backend, this%config%contour_points, this%config%contour_shape, &
         this%config%contour_margin, this%config%contour_height_fraction, this%config%contour_account_fermi_poles, &
         this%config%native_crosscheck, this%config%native_turek, this%config%native_contour_points, this%config%native_contour_margin, &
         this%config%native_contour_height_fraction, this%config%native_contour_account_fermi_poles, &
         this%config%native_contour_target_fermi_poles, lattice_obj, hamiltonian_obj, energy_obj, self_obj, reciprocal_obj, fixture, &
         q_direct, q_cart, finite_total, finite_tt, finite_contact, spectral_total, spectral_tt, spectral_contact, contour_total, &
         contour_tt, contour_contact, native_jq_ud, native_jq_du, native_jq_sym, native_delta_j, native_curvature, native_report, &
         native_ready, endpoint_mode, endpoint_reused, q_commensurate, endpoint_residual, endpoint_seconds, assembly_seconds, &
         contraction_seconds, hamiltonian_seconds, gf_seconds, solve_seconds, contour_seconds, total_response_seconds)
      call run_native_rotation_dynamics_campaign(trim(this%config%output_file), lattice_obj, reciprocal_obj, self_obj, fixture, &
         q_direct, q_cart, finite_total, native_delta_j, native_curvature, native_ready, this%config%rotation_eta_ladder, &
         this%config%rotation_probe_omega, this%config%rotation_probe_eta, &
         this%config%rotation_pole_window_floor, this%config%rotation_pole_window_scale, &
         this%config%rotation_pole_window_max, this%config%rotation_pole_coarse_points, this%config%rotation_pole_fine_points, &
         this%config%rotation_pole_refinement_half_width)
      if (this%config%rotation_n_omega > 0) then
         call write_rotation_response_grid(trim(this%config%rotation_grid_file), this%config%rotation_n_omega, &
            this%config%rotation_omega_min, this%config%rotation_omega_max, this%config%rotation_eta, fixture, &
            reciprocal_obj, q_direct)
      end if

      call fixture%clear()
      deallocate(q_direct, q_cart, finite_total, finite_tt, finite_contact, spectral_total, spectral_tt, spectral_contact, &
         contour_total, contour_tt, contour_contact, native_jq_ud, native_jq_du, native_jq_sym, native_delta_j, native_curvature, &
         endpoint_reused, q_commensurate, endpoint_residual, endpoint_mode)
   end subroutine lr_run_rotation

   subroutine write_rotation_response_grid(filename,n_omega,omega_min,omega_max,eta,fixture,recip,q_direct)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: n_omega
      real(rp), intent(in) :: omega_min,omega_max,eta
      type(lmto_live_hamiltonian_fixture), intent(in) :: fixture
      type(reciprocal), intent(inout) :: recip
      real(rp), intent(in) :: q_direct(:, :)
      type(rotation_state), target :: state
      type(rotation_request) :: request
      type(rotation_result) :: response
      integer :: unit,iq,iomega
      real(rp) :: omega

      if (n_omega<1 .or. omega_max<omega_min .or. eta<=0.0_rp) then
         error stop 'write_rotation_response_grid: invalid frequency-grid controls'
      end if
      open(newunit=unit,file=trim(filename),status='replace',action='write')
      write(unit,'(a)') '# local-rotation response frequency grid; energies in Ry; loss in 1/Ry'
      write(unit,'(a)') '# q_fraction_x q_fraction_y q_fraction_z omega_Ry ReK_plus_Ry ImK_plus_Ry ReK_minus_Ry ImK_minus_Ry loss_plus_invRy loss_minus_invRy'
      do iq=1,size(q_direct,2)
         call prepare_rotation_response(fixture,recip,q_direct(:,iq),state)
         request%state=>state; request%eta=eta; request%exact_static=.false.; request%want_inverse=.true.
         do iomega=1,n_omega
            if (n_omega==1) then
               omega=omega_min
            else
               omega=omega_min+(omega_max-omega_min)*real(iomega-1,rp)/real(n_omega-1,rp)
            end if
            request%omega=omega
            call evaluate_rotation_response(request,response)
            write(unit,'(3(es20.12,1x),7(es20.12,1x))') q_direct(:,iq),omega, &
               real(response%kernel_pm(1,1),rp),aimag(response%kernel_pm(1,1)), &
               real(response%kernel_pm(2,2),rp),aimag(response%kernel_pm(2,2)), &
               -aimag(response%inverse_kernel_pm(1,1)),-aimag(response%inverse_kernel_pm(2,2))
         end do
         call state%clear()
      end do
      close(unit)
   end subroutine write_rotation_response_grid

   subroutine run_native_rotation_dynamics_campaign(rotation_output_file,lat,recip,self_obj,fixture,q_direct,q_cart,finite_h,native_delta,native_curv,native_ready, &
      rotation_eta_ladder,rotation_probe_omega,rotation_probe_eta,rotation_pole_window_floor, &
      rotation_pole_window_scale,rotation_pole_window_max,rotation_pole_coarse_points,rotation_pole_fine_points,rotation_pole_refinement_half_width)
      character(len=*), intent(in) :: rotation_output_file
      type(lattice), intent(in) :: lat
      type(reciprocal), target, intent(inout) :: recip
      type(self), intent(in) :: self_obj
      type(lmto_live_hamiltonian_fixture), target, intent(in) :: fixture
      real(rp), intent(in) :: q_direct(:, :),q_cart(:, :),finite_h(:, :, :),native_delta(:),native_curv(:)
      logical, intent(in) :: native_ready
      real(rp), intent(in) :: rotation_eta_ladder(3), rotation_probe_omega, rotation_probe_eta
      real(rp), intent(in) :: rotation_pole_window_floor, rotation_pole_window_scale, rotation_pole_window_max
      integer, intent(in) :: rotation_pole_coarse_points, rotation_pole_fine_points
      real(rp), intent(in) :: rotation_pole_refinement_half_width
      type(rotation_state), target :: state,minus_state
      type(rotation_request) :: request
      type(rotation_result) :: response,minus_response,plus_probe,minus_probe
      complex(rp), allocatable :: reduced(:, :)
      real(rp), parameter :: ry_to_mev=13605.693122994_rp
      real(rp) :: bplus,bminus,slope_rel,slope_step,qmag,scf_moment,field_max,cov_q(3)
      real(rp) :: static_residual,q0_residual,circular_offdiag,covariance_residual,omega_cov,eta_cov
      real(rp) :: finite_mev,turek_mev,omega_max,pred_plus,pred_minus,expected,df,fit_d,fit_resid
      real(rp) :: pole_re,pole_min,loss_peak,loss_height,fwhm,fwhm_mev,resolution,re_k,im_k,abs_k,pole_slope
      real(rp) :: energy_eta(3),resid_eta(3),eta_move,window_min,window_max
      integer :: unit,iq,igamma,ik,site,channel,ieta,ipole,nfinite,clean_count
      logical :: static_ok,berry_ok,circular_ok,covariance_ok,causal_ok,pole_ok(3),resolved,fit_ok
      character(len=8) :: channel_name
      character(len=24) :: gate_status
      character(len=10) :: pole_status

      if (fixture%nsite/=1) error stop 'first Fe rotation pole driver currently reports the one-site primitive bcc state'
      if (recip%hamiltonian%ccor_2c) error stop 'first Fe rotation pole state requires CCOR off'
      if (abs(recip%temperature-300.0_rp)>1.0e-8_rp) error stop 'first Fe rotation pole state requires T=300 K'
      if (len_trim(rotation_output_file)==0) error stop 'rotation dynamics output path is blank'
      open(newunit=unit,file=trim(rotation_output_file),status='replace',action='write')
      write(unit,'(a)') '# native second-order local-rotation dynamics; K^R = contact + unsymmetrized retarded bubble'
      write(unit,'(a)') '# energies are Ry internally; reported q_Ainv uses 2*pi/alat times the actual Cartesian reciprocal vector'
      write(unit,'(a,l1)') '# native_turek_diagnostic_enabled = ',native_ready
      write(unit,'(a)') '# pole_window_basis = H2 static rotation curvature from finite-H; Turek is diagnostic only'
      write(unit,'(a)') '# missing diagnostic fields are written as -'
      write(unit,'(a)') '# q_fraction q_Ainv adiabatic_finiteH_meV adiabatic_Turek_meV_or_missing pole_ReK_meV pole_minK_meV loss_peak_meV eta_Ry FWHM_meV max_minus_ImG_invRy circular_channel ReK_Ry ImK_Ry absK_Ry frequency_resolution_Ry status'
      write(*,'(a,3(i0,1x))') 'Rotation dynamics accepted k mesh = ',recip%nk_mesh
      write(*,'(a,es16.8)') 'Rotation dynamics EF (Ry) = ',recip%fermi_level
      write(*,'(a,es16.8)') 'Rotation dynamics temperature (K) = ',recip%temperature
      write(*,'(a,es16.8)') 'Rotation dynamics SCF constraining field max (Ry) = ',max_constraint_field(self_obj,lat%nrec,lat%nbulk)

      igamma=0
      do iq=1,size(q_direct,2)
         if (maxval(abs(q_direct(:,iq)))<q_tolerance) igamma=iq
      end do
      if (igamma==0) error stop 'rotation dynamics q list has no Gamma point'
      allocate(reduced(2,2))
      call prepare_rotation_response(fixture,recip,q_direct(:,igamma),state)
      request%state=>state; request%omega=0.0_rp; request%eta=0.0_rp; request%exact_static=.true.; request%want_inverse=.false.
      call evaluate_rotation_response(request,response)
      call reduce_static_rotation_kernel(response%kernel,reduced)
      q0_residual=maxval(abs(reduced))
      ! One canonical h controls both the sampled frequencies and the
      ! symmetric-difference denominator.  rotation_slope_step is retained
      ! only as a deprecated namelist compatibility field.
      slope_step=rotation_probe_omega
      request%exact_static=.false.; request%omega=slope_step; request%eta=rotation_probe_eta
      call evaluate_rotation_response(request,plus_probe)
      request%omega=-slope_step
      call evaluate_rotation_response(request,minus_probe)
      bplus=real((plus_probe%kernel_pm(1,1)-minus_probe%kernel_pm(1,1))/(2.0_rp*slope_step),rp)
      bminus=real((plus_probe%kernel_pm(2,2)-minus_probe%kernel_pm(2,2))/(2.0_rp*slope_step),rp)
      slope_rel=max(abs(abs(bplus)-abs(response%berry)),abs(abs(bminus)-abs(response%berry)),abs(bplus+bminus))/ &
         max(abs(response%berry),1.0e-12_rp)
      circular_offdiag=max(abs(reduced(1,1)-reduced(2,2)),abs(reduced(1,2)),abs(reduced(2,1)))
      scf_moment=0.0_rp
      do site=1,lat%nrec
         scf_moment=scf_moment+self_obj%symbolic_atom(lat%nbulk+site)%potential%mtot
      end do
      write(*,'(a,es16.8)') 'Rotation dynamics electron-number residual = ',state%electron_residual
      write(*,'(a,es16.8)') 'Rotation dynamics M_band = ',state%magnetization
      write(*,'(a,es16.8)') 'Rotation dynamics SCF magnetic moment = ',scf_moment
      write(*,'(a,es16.8)') 'Rotation dynamics Berry commutator = ',response%berry
      write(*,'(a,es12.4)') 'Rotation dynamics direct Berry vs M/2 residual = ',response%berry_residual
      write(*,'(a,2(es16.8,1x),a,es12.4)') 'Rotation dynamics circular slopes (+,-) = ',bplus,bminus, &
         ' relative residual=',slope_rel
      write(*,'(a,es12.4)') 'Rotation dynamics q0 Goldstone residual = ',q0_residual
      write(*,'(a,es12.4)') 'Rotation dynamics circular off-diagonal residual = ',max(abs(reduced(1,2)),abs(reduced(2,1)))
      berry_ok=slope_rel<=5.0e-2_rp .and. abs(abs(bplus)-abs(bminus))<=5.0e-2_rp*max(abs(response%berry),1.0e-12_rp)
      circular_ok=max(abs(reduced(1,2)),abs(reduced(2,1)))<=1.0e-7_rp
      q0_residual=maxval(abs(reduced))
      static_ok=q0_residual<=2.0e-7_rp
      if (abs(response%berry)<=1.0e-8_rp) berry_ok=.false.

      write(unit,'(a)') '# state_provenance = accepted second-order ham_only; SOC off; CCOR off; constraining field zero'
      write(unit,'(a,3(i0,1x))') '# k_mesh = ',recip%nk_mesh
      write(unit,'(a,es20.12)') '# EF_Ry = ',recip%fermi_level
      write(unit,'(a,es20.12)') '# electron_number_residual = ',state%electron_residual
      write(unit,'(a,es20.12)') '# M_band = ',state%magnetization
      write(unit,'(a,es20.12)') '# SCF_magnetic_moment = ',scf_moment
      write(unit,'(a,es20.12)') '# Berry_commutator = ',response%berry
      write(unit,'(a,es20.12)') '# direct_Berry_vs_Mband_over_2_residual = ',response%berry_residual
      write(unit,'(a,es20.12)') '# q0_slope_plus = ',bplus
      write(unit,'(a,es20.12)') '# q0_slope_minus = ',bminus
      write(unit,'(a,es20.12)') '# q0_Goldstone_residual = ',q0_residual
      write(unit,'(a,es20.12)') '# circular_offdiagonal_residual = ',max(abs(reduced(1,2)),abs(reduced(2,1)))

      pole_ok=.false.; energy_eta=-1.0_rp; resid_eta=huge(1.0_rp); nfinite=0; clean_count=0
      covariance_ok=.false.; covariance_residual=huge(1.0_rp); causal_ok=.false.
      do iq=1,size(q_direct,2)
         if (iq==igamma) cycle
         nfinite=nfinite+1
         call prepare_rotation_response(fixture,recip,q_direct(:,iq),state)
         request%state=>state; request%omega=0.0_rp; request%eta=0.0_rp; request%exact_static=.true.; request%want_inverse=.false.
         call evaluate_rotation_response(request,response)
         finite_mev=2.0_rp*real(finite_h(1,1,iq),rp)/max(abs(state%magnetization),1.0e-12_rp)*ry_to_mev
         if (native_ready) then
            turek_mev=4.0_rp*native_delta(iq)/max(abs(state%magnetization),1.0e-12_rp)*ry_to_mev
         else
            ! Keep the internal value harmless, but never serialize it as a
            ! physical zero.  The row writer emits an explicit missing marker.
            turek_mev=0.0_rp
         end if
         qmag=2.0_rp*pi/lat%alat*sqrt(sum(q_cart(:,iq)**2))
         if (iq==2) then
            call prepare_rotation_response(fixture,recip,-q_direct(:,iq),minus_state)
            request%state=>minus_state; request%exact_static=.true.; request%omega=0.0_rp; request%eta=0.0_rp
            call evaluate_rotation_response(request,minus_response)
            call reduce_static_rotation_kernel(response%kernel,reduced)
            static_residual=maxval(abs(reduced(1:1,1:1)-finite_h(1:1,1:1,iq)))
            write(*,'(a,es12.4)') 'Rotation dynamics q,-q static reduction residual = ',static_residual
            static_ok=static_ok .and. static_residual<=max(1.0e-8_rp,1.0e-7_rp*maxval(abs(finite_h(:,:,iq))))
            ! The covariance identity is evaluated at a non-self-inverse q
            ! that translates the sampled mesh exactly; arbitrary q remains
            ! supported by the response API and the separate direct oracle.
            cov_q=[1.0_rp/real(recip%nk_mesh(1),rp),0.0_rp,0.0_rp]
            call state%clear(); call minus_state%clear()
            call prepare_rotation_response(fixture,recip,cov_q,state)
            call prepare_rotation_response(fixture,recip,-cov_q,minus_state)
            request%state=>state; request%exact_static=.false.; request%omega=2.0e-4_rp; request%eta=5.0e-5_rp
            call evaluate_rotation_response(request,plus_probe)
            request%state=>minus_state; request%omega=-2.0e-4_rp
            call evaluate_rotation_response(request,minus_probe)
            ! The live Bloch endpoint convention satisfies K_AB(q,w) =
            ! conjg(K_AB(-q,-w)) for real Cartesian coordinates.
            covariance_residual=maxval(abs(plus_probe%kernel-conjg(minus_probe%kernel)))/ &
               max(1.0_rp,maxval(abs(plus_probe%kernel)))
            write(*,'(a,3(es12.4,1x))') 'Rotation dynamics q/-q/-omega covariance q/residual = ',cov_q,covariance_residual
            covariance_ok=covariance_residual<=1.0e-7_rp
            call state%clear(); call minus_state%clear()
            call prepare_rotation_response(fixture,recip,q_direct(:,iq),state)
         else
            static_residual=abs(real(response%kernel(1,1),rp)-finite_h(1,1,iq))
            if (nfinite==1) static_ok=static_ok .and. static_residual<=max(1.0e-8_rp,1.0e-7_rp*abs(finite_h(1,1,iq)))
         end if
         if (iq==2) then
            pred_plus=-real(response%kernel_pm(1,1),rp)/bplus
            pred_minus=-real(response%kernel_pm(2,2),rp)/bminus
         else
            pred_plus=-real(response%kernel_pm(1,1),rp)/bplus
            pred_minus=-real(response%kernel_pm(2,2),rp)/bminus
         end if
         channel=0
         if (pred_plus>0.0_rp .and. pred_minus<=0.0_rp) channel=1
         if (pred_minus>0.0_rp .and. pred_plus<=0.0_rp) channel=2
         if (channel==0 .and. pred_plus>0.0_rp .and. pred_minus>0.0_rp) then
            if (abs(pred_plus-finite_mev/ry_to_mev)<=abs(pred_minus-finite_mev/ry_to_mev)) then
               channel=1
            else
               channel=2
            end if
         end if
         if (channel==1) then
            channel_name='+ '
            expected=pred_plus
         else if (channel==2) then
            channel_name='- '
            expected=pred_minus
         else
            channel_name='none'
            expected=0.0_rp
         end if
         ! The H2 finite-H curvature is the production static scale.  Turek
         ! is an independent diagnostic and must not influence branch,
         ! window, or pole acceptance decisions.
         omega_max=max(rotation_pole_window_floor,rotation_pole_window_scale*abs(finite_mev)/ry_to_mev)
         omega_max=min(omega_max,rotation_pole_window_max)
         if (native_ready) then
            write(*,'(a,i0,a,3(es12.4,1x),a,es12.4,a,es12.4)') 'Rotation q index ',iq,' finite-H/Turek/meV=', &
               finite_mev,turek_mev,qmag,' A^-1 omega_window_Ry=',omega_max,' predicted_Ry=',expected
         else
            write(*,'(a,i0,a,es12.4,a,es12.4,a,es12.4)') 'Rotation q index ',iq,' finite-H/meV=',finite_mev, &
               ' q_A^-1=',qmag,' omega_window_Ry=',omega_max
         end if
         if (channel>0) then
            if (channel==1) then
               ipole=1
            else
               ipole=2
            end if
            do ieta=1,size(rotation_eta_ladder)
               call scan_rotation_pole(state,ipole,omega_max,rotation_eta_ladder(ieta),rotation_pole_coarse_points, &
                  rotation_pole_fine_points,rotation_pole_refinement_half_width,pole_re,pole_min,loss_peak,loss_height, &
                  fwhm,resolution,re_k,im_k,abs_k,pole_slope,resolved)
               pole_ok(nfinite)=resolved .or. pole_ok(nfinite)
               causal_ok=causal_ok .or. (pole_re>0.0_rp .and. pole_min>0.0_rp .and. pole_slope*im_k>0.0_rp)
               if (ieta==size(rotation_eta_ladder)) then
                  energy_eta(nfinite)=pole_re*ry_to_mev
                  resid_eta(nfinite)=resolution*ry_to_mev
               end if
               if (resolved) then
                  pole_status='RESOLVED'
               else
                  pole_status='UNRESOLVED'
               end if
               fwhm_mev=-1.0_rp
               if (fwhm>=0.0_rp) fwhm_mev=fwhm*ry_to_mev
               call write_rotation_dynamics_row(unit,q_direct(:,iq),qmag,finite_mev,native_ready,turek_mev, &
                  pole_re*ry_to_mev,pole_min*ry_to_mev,loss_peak*ry_to_mev,rotation_eta_ladder(ieta),fwhm_mev,loss_height, &
                  trim(channel_name),re_k,im_k,abs_k,resolution,trim(pole_status))
               write(*,'(a,3(es14.6,1x),a,es12.4,a,a,a,l1)') '  pole ReK/minK/loss (meV) = ', &
                  pole_re*ry_to_mev,pole_min*ry_to_mev,loss_peak*ry_to_mev,' eta=',rotation_eta_ladder(ieta), &
                  ' channel=',trim(channel_name),' resolved=',resolved
            end do
            if (pole_ok(nfinite)) clean_count=clean_count+1
         else
            call write_rotation_dynamics_row(unit,q_direct(:,iq),qmag,finite_mev,native_ready,turek_mev, &
               -1.0_rp,-1.0_rp,-1.0_rp,rotation_eta_ladder(3),-1.0_rp,0.0_rp,'none', &
               0.0_rp,0.0_rp,0.0_rp,0.0_rp,'NO_POSITIVE_CHANNEL')
         end if
         call state%clear()
         call minus_state%clear()
      end do

      ! covariance_ok is set only by the explicit q/-q/-omega comparison above.
      berry_ok=berry_ok .and. slope_rel<=5.0e-2_rp
      static_ok=static_ok .and. q0_residual<=2.0e-7_rp
      fit_ok=.false.; fit_d=0.0_rp; fit_resid=-1.0_rp; window_min=0.0_rp; window_max=0.0_rp
      if (nfinite>=3 .and. clean_count>=3) then
         fit_ok=.true.; expected=0.0_rp; window_min=huge(1.0_rp); window_max=0.0_rp
         do ipole=1,3
            if (energy_eta(ipole)<=0.0_rp) fit_ok=.false.
            if (ipole>1 .and. energy_eta(ipole)<=energy_eta(ipole-1)) fit_ok=.false.
            expected=expected+energy_eta(ipole)*(2.0_rp*pi/lat%alat*sqrt(sum(q_cart(:,ipole+1)**2)))**2
            fit_d=fit_d+(2.0_rp*pi/lat%alat*sqrt(sum(q_cart(:,ipole+1)**2)))**4
            window_min=min(window_min,2.0_rp*pi/lat%alat*sqrt(sum(q_cart(:,ipole+1)**2)))
            window_max=max(window_max,2.0_rp*pi/lat%alat*sqrt(sum(q_cart(:,ipole+1)**2)))
         end do
         if (fit_ok .and. fit_d>tiny(1.0_rp)) then
            fit_d=expected/fit_d
            fit_resid=0.0_rp
            do ipole=1,3
               qmag=2.0_rp*pi/lat%alat*sqrt(sum(q_cart(:,ipole+1)**2))
               fit_resid=fit_resid+(energy_eta(ipole)-fit_d*qmag*qmag)**2
            end do
            fit_resid=sqrt(fit_resid/3.0_rp)
         else
            fit_ok=.false.
         end if
      end if
      if (fit_ok .and. clean_count>=3 .and. static_ok .and. berry_ok .and. circular_ok .and. covariance_ok .and. causal_ok) then
         gate_status='PASS-A'
      else if (static_ok .and. berry_ok .and. circular_ok .and. covariance_ok .and. causal_ok) then
         gate_status='PASS-B'
      else
         gate_status='BLOCKED'
      end if
      write(*,'(a,l1)') 'Rotation dynamics causal positive-frequency pole = ',causal_ok
      write(*,'(a,a)') 'Rotation dynamics implementation = ',trim(gate_status)
      if (fit_ok) then
         write(*,'(a)') 'Fe material convergence = PRELIMINARY — MATERIAL CONVERGENCE OPEN'
      else
         write(*,'(a)') 'Fe material convergence = OPEN'
      end if
      if (fit_ok) write(*,'(a,es16.8,a,es16.8,a,2(es12.4,1x))') 'Preliminary D (meV A2) = ',fit_d, &
         ' fit residual (meV) = ',fit_resid,' q window (A^-1) = ',window_min,window_max
      write(*,'(a)') 'Goldstone correction = OFF; strict-ASA L0 dynamics = NOT RUN'
      close(unit)
      call state%clear(); call minus_state%clear()
      deallocate(reduced)
   end subroutine run_native_rotation_dynamics_campaign

   subroutine write_rotation_dynamics_row(unit,q,qmag,finite_mev,turek_available,turek_mev,pole_re,pole_min,loss_peak,eta, &
      fwhm,loss_height,channel_name,re_k,im_k,abs_k,resolution,status)
      integer, intent(in) :: unit
      real(rp), intent(in) :: q(3),qmag,finite_mev,turek_mev,pole_re,pole_min,loss_peak,eta,fwhm,loss_height
      logical, intent(in) :: turek_available
      character(len=*), intent(in) :: channel_name,status
      real(rp), intent(in) :: re_k,im_k,abs_k,resolution

      if (turek_available) then
         write(unit,'(3(es14.6,1x),9(es14.6,1x),a,1x,4(es14.6,1x),a)') q,qmag,finite_mev,turek_mev, &
            pole_re,pole_min,loss_peak,eta,fwhm,loss_height,trim(channel_name),re_k,im_k,abs_k,resolution,trim(status)
      else
         write(unit,'(3(es14.6,1x),2(es14.6,1x),a,1x,6(es14.6,1x),a,1x,4(es14.6,1x),a)') q,qmag,finite_mev,'-', &
            pole_re,pole_min,loss_peak,eta,fwhm,loss_height,trim(channel_name),re_k,im_k,abs_k,resolution,trim(status)
      end if
   end subroutine write_rotation_dynamics_row

   subroutine scan_rotation_pole(state,channel,omega_max,eta,coarse_points,fine_points,refinement_half_width, &
      pole_re,pole_min,loss_peak,loss_height,fwhm,resolution,re_k,im_k,abs_k,pole_slope,resolved)
      type(rotation_state), target, intent(inout) :: state
      integer, intent(in) :: channel
      real(rp), intent(in) :: omega_max,eta
      integer, intent(in) :: coarse_points,fine_points
      real(rp), intent(in) :: refinement_half_width
      real(rp), intent(out) :: pole_re,pole_min,loss_peak,loss_height,fwhm,resolution,re_k,im_k,abs_k,pole_slope
      logical, intent(out) :: resolved
      type(rotation_request) :: request
      type(rotation_result) :: response
      integer :: ncoarse,nfine
      real(rp), allocatable :: omega_c(:),kabs_c(:),kre_c(:),loss_c(:)
      real(rp), allocatable :: omega_f(:),kabs_f(:),kre_f(:),loss_f(:)
      real(rp) :: center,dw,left,right,half,peak,root_distance,x1,x2,t
      integer :: i,imin,ipeak,iroot
      logical :: have_root,have_left,have_right
      ncoarse=coarse_points; nfine=fine_points
      if (ncoarse<2 .or. nfine<2 .or. refinement_half_width<=0.0_rp) then
         error stop 'scan_rotation_pole: invalid scan controls'
      end if
      allocate(omega_c(ncoarse),kabs_c(ncoarse),kre_c(ncoarse),loss_c(ncoarse), &
         omega_f(nfine),kabs_f(nfine),kre_f(nfine),loss_f(nfine))
      request%state=>state; request%eta=eta; request%exact_static=.false.; request%want_inverse=.true.
      do i=1,ncoarse
         omega_c(i)=omega_max*real(i-1,rp)/real(ncoarse-1,rp)
         request%omega=omega_c(i)
         call evaluate_rotation_response(request,response)
         kabs_c(i)=abs(response%kernel_pm(channel,channel))
         kre_c(i)=real(response%kernel_pm(channel,channel),rp)
         loss_c(i)=-aimag(response%inverse_kernel_pm(channel,channel))
      end do
      imin=minloc(kabs_c,dim=1); center=omega_c(imin); dw=omega_max/real(ncoarse-1,rp)
      left=max(0.0_rp,center-refinement_half_width*dw); right=min(omega_max,center+refinement_half_width*dw)
      if (right<=left) right=min(omega_max,left+dw)
      resolution=(right-left)/real(nfine-1,rp)
      do i=1,nfine
         omega_f(i)=left+resolution*real(i-1,rp)
         request%omega=omega_f(i)
         call evaluate_rotation_response(request,response)
         kabs_f(i)=abs(response%kernel_pm(channel,channel))
         kre_f(i)=real(response%kernel_pm(channel,channel),rp)
         loss_f(i)=-aimag(response%inverse_kernel_pm(channel,channel))
      end do
      imin=minloc(kabs_f,dim=1); pole_min=omega_f(imin); pole_re=-1.0_rp; have_root=.false.; iroot=0
      root_distance=huge(1.0_rp)
      do i=1,nfine-1
         if (kre_f(i)==0.0_rp) then
            t=0.0_rp
         else if (kre_f(i)*kre_f(i+1)<0.0_rp) then
            t=abs(kre_f(i))/(abs(kre_f(i))+abs(kre_f(i+1)))
         else
            cycle
         end if
         x1=omega_f(i)+t*(omega_f(i+1)-omega_f(i))
         if (abs(x1-center)<root_distance) then
            root_distance=abs(x1-center); pole_re=x1; iroot=i; have_root=.true.
         end if
      end do
      ipeak=maxloc(loss_f,dim=1); loss_peak=omega_f(ipeak); peak=loss_f(ipeak); loss_height=peak; fwhm=-1.0_rp
      half=0.5_rp*peak; have_left=.false.; have_right=.false.; x1=left; x2=right
      if (peak>0.0_rp) then
         do i=ipeak-1,1,-1
            if (loss_f(i)<=half .and. loss_f(i+1)>half) then
               t=(half-loss_f(i))/(loss_f(i+1)-loss_f(i))
               x1=omega_f(i)+t*(omega_f(i+1)-omega_f(i)); have_left=.true.; exit
            end if
         end do
         do i=ipeak,nfine-1
            if (loss_f(i)>half .and. loss_f(i+1)<=half) then
               t=(half-loss_f(i))/(loss_f(i+1)-loss_f(i))
               x2=omega_f(i)+t*(omega_f(i+1)-omega_f(i)); have_right=.true.; exit
            end if
         end do
         if (have_left .and. have_right) fwhm=x2-x1
      end if
      pole_slope=0.0_rp
      if (iroot>0) pole_slope=(kre_f(iroot+1)-kre_f(iroot))/(omega_f(iroot+1)-omega_f(iroot))
      ! Report the kernel components at the refined Re K=0 crossing.  The
      ! separate pole_min diagnostic remains the minimum of |K|.
      request%omega=merge(pole_re,pole_min,pole_re>0.0_rp)
      call evaluate_rotation_response(request,response)
      re_k=real(response%kernel_pm(channel,channel),rp)
      im_k=aimag(response%kernel_pm(channel,channel))
      abs_k=abs(response%kernel_pm(channel,channel))
      resolved=have_root .and. peak>0.0_rp .and. abs(pole_re-pole_min)<=max(2.0_rp*eta,2.0_rp*resolution)
      if (peak>0.0_rp) resolved=resolved .and. abs(loss_peak-pole_min)<=max(2.0_rp*eta,2.0_rp*resolution)
      if (iroot==0) resolved=.false.
      deallocate(omega_c,kabs_c,kre_c,loss_c,omega_f,kabs_f,kre_f,loss_f)
   end subroutine scan_rotation_pole

   real(rp) function max_constraint_field(self_obj,nsite,nbulk) result(value)
      type(self), intent(in) :: self_obj
      integer, intent(in) :: nsite,nbulk
      integer :: site
      value=0.0_rp
      do site=1,nsite
         value=max(value,maxval(abs(self_obj%symbolic_atom(nbulk+site)%mag_cfield)))
      end do
   end function max_constraint_field


end submodule linear_response_rotation
