!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Finite-Hamiltonian Katsnelson--Lichtenstein static bridge.
!>
!> This module is deliberately separate from `exchange_mod`.  The latter is
!> the Rung-0 native LKAG implementation and remains unchanged.  The present
!> service evaluates the spectral finite-Hamiltonian contraction using an
!> explicitly supplied, site-local spin-flip exchange vertex.  In particular,
!> it never invents a scalar site splitting from an operator vertex.
!>
!> For the source convention used by `exchange%dGdG_Jnc`, the zero-frequency
!> spectral identity is
!>
!>   J_ij = -1/4 sum_(k,n,m) w_k (f_n-f_m)
!>                   X_i(n,m) X_j(n,m)^* /(e_n-e_m+i eta),
!>
!> where X_i=<n,k|Delta_i|m,k+q>.  The minus sign and 1/4 are not adjustable:
!> they follow from the imaginary part of the one-energy LKAG contour and the
!> source `J = Im Tr[...]` convention.  The corresponding DRESP-02 chi0
!> convention has a factor two, so an exact scalar compression obeys
!>
!>   J = -1/8 Delta_i chi0_ij Delta_j.
!------------------------------------------------------------------------------
module lr_kl_static_bridge_mod

   use, intrinsic :: ieee_arithmetic
   use precision_mod, only: rp
   use lmto_radial_augmentation_mod, only: lmto_orbital_l
   use lr_ks_susceptibility_mod, only: lr_electronic_state
   implicit none
   private

   real(rp), parameter, public :: kl_static_prefactor = -0.25_rp
   real(rp), parameter, public :: kl_scalar_chi_prefactor = -0.125_rp
   real(rp), parameter :: state_tolerance = 2.0e-11_rp

   character(len=*), parameter, public :: kl_vertex_ham_only = &
      'finite ham_only local H_up-H_down spin-flip vertex'
   character(len=*), parameter, public :: kl_vertex_operator = &
      'explicit operator-valued exchange vertex'

   !> Request for the static finite-Hamiltonian K/L contraction.
   type, public :: kl_static_bridge_request
      real(rp) :: q(3) = 0.0_rp
      real(rp) :: eta = 0.0_rp
      character(len=64) :: vertex_source = kl_vertex_operator
      type(lr_electronic_state), pointer :: electronic_state => null()
      type(lr_electronic_state), pointer :: q_endpoint_state => null()
      ! (basis,basis,site), with the spin-flip block in the same global
      ! coefficient ordering as the DRESP-02 eigenvectors.
      complex(rp), pointer :: exchange_vertex(:, :, :) => null()
   end type kl_static_bridge_request

   !> Static bridge output in the native internal LKAG energy normalization.
   type, public :: kl_static_bridge_result
      real(rp) :: q(3) = 0.0_rp
      real(rp) :: eta = 0.0_rp
      character(len=64) :: vertex_source = ''
      character(len=128) :: representation = ''
      integer :: nsite = 0
      integer :: nbasis = 0
      integer :: nk = 0
      integer :: ntransitions_evaluated = 0
      integer :: noccupation_skips = 0
      complex(rp), allocatable :: exchange(:, :) ! (site_i,site_j)
   end type kl_static_bridge_result

   !> Diagnostic for a supplied, never-fitted scalar compression.
   type, public :: kl_scalar_vertex_diagnostic
      logical :: exact = .false.
      integer :: nsite = 0
      real(rp) :: max_absolute_residual = 0.0_rp
      real(rp) :: relative_residual = 0.0_rp
      real(rp) :: operator_norm = 0.0_rp
   end type kl_scalar_vertex_diagnostic

   public :: evaluate_kl_static_bridge
   public :: build_ham_only_exchange_vertex
   public :: assess_scalar_exchange_vertex
   public :: contract_scalar_site_chi0

contains

   !> Evaluate the explicit operator-vertex finite-Hamiltonian bridge.
   !>
   !> The transition amplitudes are formed from the supplied operator itself,
   !> not from a scalar Delta and not from a fitted ratio to chi0.  The same
   !> immutable eigenstate snapshots used by DRESP-02 are consumed directly.
   subroutine evaluate_kl_static_bridge(request, result)
      type(kl_static_bridge_request), intent(in) :: request
      type(kl_static_bridge_result), intent(out) :: result

      type(lr_electronic_state), pointer :: left_state, right_state
      complex(rp), allocatable :: transition(:), action(:)
      complex(rp) :: denominator, pair_factor
      real(rp) :: weight_sum, occupation_difference
      integer :: ik, ib, jb, isite, nsite, nbasis

      call validate_request(request)
      left_state => request%electronic_state
      right_state => request%q_endpoint_state
      nsite = size(request%exchange_vertex, 3)
      nbasis = left_state%nbasis
      weight_sum = sum(left_state%k_weights)

      allocate(result%exchange(nsite, nsite), transition(nsite), action(nbasis))
      result%exchange = cmplx(0.0_rp, 0.0_rp, rp)
      result%q = request%q
      result%eta = request%eta
      result%vertex_source = request%vertex_source
      result%representation = 'orthogonal finite-Hamiltonian operator-vertex spectral contraction'
      result%nsite = nsite
      result%nbasis = nbasis
      result%nk = left_state%nk

      do ik = 1, left_state%nk
         do ib = 1, left_state%nbands
            do jb = 1, right_state%nbands
               occupation_difference = left_state%occupations(ib, ik) - right_state%occupations(jb, ik)
               if (abs(occupation_difference) <= tiny(1.0_rp)) then
                  result%noccupation_skips = result%noccupation_skips + 1
                  cycle
               end if
               denominator = cmplx(left_state%eigenvalues(ib, ik) - right_state%eigenvalues(jb, ik), request%eta, rp)
               if (abs(denominator) <= tiny(1.0_rp)) then
                  error stop 'evaluate_kl_static_bridge: zero static transition denominator'
               end if
               do isite = 1, nsite
                  action = matmul(request%exchange_vertex(:, :, isite), right_state%eigenvectors(:, jb, ik))
                  transition(isite) = dot_product(left_state%eigenvectors(:, ib, ik), action)
               end do
               pair_factor = cmplx(kl_static_prefactor*left_state%k_weights(ik)*occupation_difference/weight_sum, &
                                   0.0_rp, rp)/denominator
               do isite = 1, nsite
                  result%exchange(:, isite) = result%exchange(:, isite) + &
                     pair_factor*transition*conjg(transition(isite))
               end do
               result%ntransitions_evaluated = result%ntransitions_evaluated + 1
            end do
         end do
      end do

      deallocate(transition, action)
   end subroutine evaluate_kl_static_bridge

   !> Build the finite `ham_only` local exchange vertex from the coefficient
   !> Hamiltonian blocks.  `up-down` follows the native `rs2pao` source trace;
   !> the returned matrix is the plus-channel block (up row, down column).
   !>
   !> This routine only embeds a supplied finite Hamiltonian difference.  It
   !> does not claim that this object equals native LMTO `d_matrix(E)`.
   subroutine build_ham_only_exchange_vertex(up_blocks, down_blocks, selected_l, vertices)
      complex(rp), intent(in) :: up_blocks(:, :, :), down_blocks(:, :, :)
      logical, intent(in) :: selected_l(0:)
      complex(rp), allocatable, intent(out) :: vertices(:, :, :)

      integer :: norb, nsite, site, iorb, jorb, l_left, l_right, offset

      if (size(up_blocks, 1) /= size(up_blocks, 2) .or. any(shape(down_blocks) /= shape(up_blocks))) then
         error stop 'build_ham_only_exchange_vertex: up/down block shape mismatch'
      end if
      norb = size(up_blocks, 1)
      nsite = size(up_blocks, 3)
      if (norb < 1 .or. size(selected_l) < 1) then
         error stop 'build_ham_only_exchange_vertex: invalid orbital or selector shape'
      end if
      allocate(vertices(2*norb*nsite, 2*norb*nsite, nsite))
      vertices = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, nsite
         offset = (site - 1)*2*norb
         do iorb = 1, norb
            l_left = lmto_orbital_l(iorb)
            if (l_left < lbound(selected_l, 1) .or. l_left > ubound(selected_l, 1) .or. &
                .not. selected_l(l_left)) cycle
            do jorb = 1, norb
               l_right = lmto_orbital_l(jorb)
               if (l_right < lbound(selected_l, 1) .or. l_right > ubound(selected_l, 1) .or. &
                   .not. selected_l(l_right)) cycle
               vertices(offset + iorb, offset + norb + jorb, site) = &
                  up_blocks(iorb, jorb, site) - down_blocks(iorb, jorb, site)
            end do
         end do
      end do
   end subroutine build_ham_only_exchange_vertex

   !> Assess a predeclared scalar site splitting without fitting it.
   !>
   !> The reference operator is Delta(site) times the selected orbital identity
   !> in the plus spin block.  A nonzero residual is evidence that scalar site
   !> compression is not exact; no best-fit Delta is computed here.
   subroutine assess_scalar_exchange_vertex(vertices, selected_l, scalar_splitting, diagnostic, tolerance)
      complex(rp), intent(in) :: vertices(:, :, :)
      logical, intent(in) :: selected_l(0:)
      real(rp), intent(in) :: scalar_splitting(:)
      type(kl_scalar_vertex_diagnostic), intent(out) :: diagnostic
      real(rp), intent(in), optional :: tolerance

      integer :: nbasis, nsite, norb, site, iorb, jorb, offset, l_left, l_right
      complex(rp) :: expected
      real(rp) :: residual, norm_value, tol

      nsite = size(vertices, 3)
      if (size(vertices, 1) /= size(vertices, 2) .or. nsite < 1 .or. size(scalar_splitting) /= nsite) then
         error stop 'assess_scalar_exchange_vertex: vertex/scalar shape mismatch'
      end if
      if (mod(size(vertices, 1), 2*nsite) /= 0) then
         error stop 'assess_scalar_exchange_vertex: basis is not site-blocked spin basis'
      end if
      norb = size(vertices, 1)/(2*nsite)
      nbasis = size(vertices, 1)
      tol = 100.0_rp*epsilon(1.0_rp)
      if (present(tolerance)) tol = tolerance
      diagnostic%nsite = nsite
      diagnostic%max_absolute_residual = 0.0_rp
      diagnostic%operator_norm = sqrt(sum(abs(vertices)**2))
      do site = 1, nsite
         offset = (site - 1)*2*norb
         do iorb = 1, 2*norb
            do jorb = 1, 2*norb
               expected = cmplx(0.0_rp, 0.0_rp, rp)
               if (iorb <= norb .and. jorb > norb) then
                  l_left = lmto_orbital_l(iorb)
                  l_right = lmto_orbital_l(jorb - norb)
                  if (l_left >= lbound(selected_l, 1) .and. l_left <= ubound(selected_l, 1) .and. &
                      l_right >= lbound(selected_l, 1) .and. l_right <= ubound(selected_l, 1) .and. &
                      l_left == l_right .and. selected_l(l_left) .and. selected_l(l_right) .and. iorb == jorb - norb) then
                     expected = cmplx(scalar_splitting(site), 0.0_rp, rp)
                  end if
               end if
               residual = abs(vertices(offset + iorb, offset + jorb, site) - expected)
               diagnostic%max_absolute_residual = max(diagnostic%max_absolute_residual, residual)
            end do
         end do
      end do
      diagnostic%relative_residual = diagnostic%max_absolute_residual/max(diagnostic%operator_norm, tiny(1.0_rp))
      diagnostic%exact = diagnostic%max_absolute_residual <= tol
      ! Keep NBASIS live in bounds/debug builds and make the intended global
      ! shape explicit to readers of the source.
      if (nbasis < 1) error stop 'assess_scalar_exchange_vertex: empty vertex'
   end subroutine assess_scalar_exchange_vertex

   !> Apply the exact scalar compression identity to a DRESP-02 site matrix.
   !> This is an algebraic helper only; it never infers `scalar_splitting`.
   subroutine contract_scalar_site_chi0(chi0, scalar_splitting, exchange)
      complex(rp), intent(in) :: chi0(:, :)
      real(rp), intent(in) :: scalar_splitting(:)
      complex(rp), intent(out) :: exchange(:, :)
      integer :: i, j, nsite

      nsite = size(chi0, 1)
      if (size(chi0, 2) /= nsite .or. size(scalar_splitting) /= nsite .or. &
          any(.not. ieee_is_finite(real(chi0, rp))) .or. any(.not. ieee_is_finite(aimag(chi0)))) then
         error stop 'contract_scalar_site_chi0: input shape or finiteness failure'
      end if
      if (any(shape(exchange) /= [nsite, nsite])) then
         error stop 'contract_scalar_site_chi0: output shape mismatch'
      end if
      do i = 1, nsite
         do j = 1, nsite
            exchange(i, j) = cmplx(kl_scalar_chi_prefactor*scalar_splitting(i)*scalar_splitting(j), 0.0_rp, rp)*chi0(i, j)
         end do
      end do
   end subroutine contract_scalar_site_chi0

   subroutine validate_request(request)
      type(kl_static_bridge_request), intent(in) :: request
      type(lr_electronic_state), pointer :: left_state, right_state

      if (.not. associated(request%electronic_state) .or. .not. associated(request%q_endpoint_state) .or. &
          .not. associated(request%exchange_vertex)) then
         error stop 'evaluate_kl_static_bridge: request references are incomplete'
      end if
      if (request%eta < 0.0_rp) then
         error stop 'evaluate_kl_static_bridge: static regulator eta must be nonnegative'
      end if
      left_state => request%electronic_state
      right_state => request%q_endpoint_state
      call left_state%validate('evaluate_kl_static_bridge:left_state')
      call right_state%validate('evaluate_kl_static_bridge:q_endpoint_state')
      if (trim(left_state%reciprocal_mode) /= 'ham_only' .or. trim(right_state%reciprocal_mode) /= 'ham_only' .or. &
          trim(left_state%hamiltonian_order) /= 'second' .or. trim(right_state%hamiltonian_order) /= 'second' .or. &
          .not. left_state%orthogonal .or. .not. right_state%orthogonal .or. .not. left_state%collinear .or. &
          .not. right_state%collinear .or. left_state%has_soc .or. right_state%has_soc .or. &
          left_state%has_extra_operator .or. right_state%has_extra_operator) then
         error stop 'evaluate_kl_static_bridge: only the certified collinear orthogonal ham_only baseline is supported'
      end if
      if (right_state%nbasis /= left_state%nbasis .or. right_state%nbands /= left_state%nbands .or. &
          right_state%nk /= left_state%nk .or. any(shape(request%exchange_vertex) /= &
          [left_state%nbasis, left_state%nbasis, size(request%exchange_vertex, 3)])) then
         error stop 'evaluate_kl_static_bridge: endpoint or exchange-vertex shape mismatch'
      end if
      if (size(request%exchange_vertex, 3) < 1 .or. any(.not. ieee_is_finite(real(request%exchange_vertex, rp))) .or. &
          any(.not. ieee_is_finite(aimag(request%exchange_vertex)))) then
         error stop 'evaluate_kl_static_bridge: exchange vertex is empty or non-finite'
      end if
      if (maxval(abs(left_state%k_weights - right_state%k_weights)) > state_tolerance) then
         error stop 'evaluate_kl_static_bridge: endpoint k weights differ from left state'
      end if
   end subroutine validate_request

end module lr_kl_static_bridge_mod
