!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Reciprocal-GF real-axis construction of the bare transverse response.
!>
!> This module is deliberately independent of the LR-06 band-pair
!> band-pair accumulator.  It evaluates the retarded/advanced Green-function
!> bubble on a real energy mesh and only shares the one-electron eigenpair and
!> Pauli/radial vertex contracts.
!------------------------------------------------------------------------------
module lr_gf_susceptibility_mod

   use precision_mod, only: rp
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, lmto_orbital_l
   use response_angular_basis_mod, only: response_angular_pi, response_gaunt
   use response_basis_mapping_mod, only: response_super_index, response_unflatten_superindex
   use lr_response_space_mod, only: response_space_layout, response_raw_to_canonical
   use lr_pauli_transition_vertex_mod, only: pauli_vertex_capabilities, pauli_sigma_plus_matrix, &
      pauli_sigma_minus_matrix
   use lr_lmto_endpoint_branches_mod, only: lmto_product_nbranch, lmto_product_max_gf_moment, &
      lmto_product_branch_powers, lmto_product_second_order_radial_branch, lmto_product_energy_power
   use lr_ks_susceptibility_mod, only: lr_electronic_state, lr_ks_susceptibility_result, &
      lr_channel_plus, lr_channel_minus, lr_fermi_dirac_occupation
   implicit none
   private

   real(rp), parameter :: state_match_tolerance = 2.0e-11_rp
   real(rp), parameter :: endpoint_match_tolerance = 2.0e-11_rp
   real(rp), parameter :: occupation_match_tolerance = 2.0e-10_rp

   !> Request for one reciprocal-GF q/frequency/channel sweep.
   !>
   !> `integration_eta` is the small broadening used only to represent the
   !> spectral discontinuity A(E)=i(G^R-G^A)/(2 pi).  `eta` remains the
   !> physical response broadening in G^R(E+omega) and G^A(E-omega).
   type, public :: lr_gf_susceptibility_request
      real(rp) :: q(3) = 0.0_rp
      real(rp), allocatable :: frequencies(:)
      real(rp) :: eta = 0.0_rp
      character(len=32) :: channel = lr_channel_plus
      integer :: integration_points = 2001
      real(rp) :: integration_eta = 0.0_rp
      real(rp) :: energy_margin = 1.0_rp
      type(response_space_layout), pointer :: response_space => null()
      type(lmto_radial_basis), pointer :: radial_bases(:) => null()
      type(lr_electronic_state), pointer :: electronic_state => null()
      type(lr_electronic_state), pointer :: q_endpoint_state => null()
   end type lr_gf_susceptibility_request

   public :: evaluate_lr_gf_susceptibility
   public :: build_weighted_resolvent
   public :: energy_power

contains

   subroutine evaluate_lr_gf_susceptibility(request, result)
      type(lr_gf_susceptibility_request), intent(in) :: request
      type(lr_ks_susceptibility_result), intent(out) :: result

      type(response_space_layout), pointer :: space
      type(lmto_radial_basis), pointer :: radial_bases(:)
      type(lr_electronic_state), pointer :: left_state, right_state
      type(pauli_vertex_capabilities) :: capabilities
      complex(rp), allocatable :: vertices(:, :, :, :)
      complex(rp), allocatable :: raw(:, :, :)
      complex(rp), allocatable :: left_gr(:, :, :), left_ga(:, :, :), left_a(:, :, :)
      complex(rp), allocatable :: right_gr(:, :, :), right_ga(:, :, :), right_a(:, :, :)
      real(rp) :: integration_eta, energy_min, energy_max, step, energy, fermi_weight
      real(rp) :: quadrature_weight, weight_sum
      complex(rp) :: operator_matrix(2, 2)
      integer :: channel_kind, nbasis, nfrequency, ne, ie, ik, ifrequency

      if (.not. associated(request%response_space) .or. .not. associated(request%radial_bases) .or. &
          .not. associated(request%electronic_state) .or. .not. associated(request%q_endpoint_state)) then
         error stop 'evaluate_lr_gf_susceptibility: request references are incomplete'
      end if
      space => request%response_space
      radial_bases => request%radial_bases
      left_state => request%electronic_state
      right_state => request%q_endpoint_state

      call validate_gf_inputs(space, radial_bases, left_state, right_state, request, channel_kind)
      call left_state%validate('evaluate_lr_gf_susceptibility:left_state')
      call right_state%validate('evaluate_lr_gf_susceptibility:q_endpoint_state')
      capabilities = pauli_vertex_capabilities()
      if (channel_kind == 1) then
         operator_matrix = pauli_sigma_plus_matrix()
      else
         operator_matrix = pauli_sigma_minus_matrix()
      end if

      integration_eta = request%integration_eta
      if (integration_eta <= 0.0_rp) integration_eta = request%eta/40.0_rp
      if (integration_eta >= request%eta) then
         error stop 'evaluate_lr_gf_susceptibility: integration_eta must be smaller than response eta'
      end if

      call build_gf_vertex_tensor(space, radial_bases, operator_matrix, capabilities, vertices)
      nbasis = left_state%nbasis
      nfrequency = size(request%frequencies)
      ne = request%integration_points
      energy_min = min(minval(left_state%eigenvalues), minval(right_state%eigenvalues)) - request%energy_margin
      energy_max = max(maxval(left_state%eigenvalues), maxval(right_state%eigenvalues)) + request%energy_margin
      if (energy_max <= energy_min) error stop 'evaluate_lr_gf_susceptibility: invalid integration interval'
      step = (energy_max - energy_min)/real(ne - 1, rp)
      weight_sum = sum(left_state%k_weights)

      allocate(raw(space%ndim, space%ndim, nfrequency), &
               left_gr(nbasis, nbasis, lmto_product_max_gf_moment + 1), &
               left_ga(nbasis, nbasis, lmto_product_max_gf_moment + 1), &
               left_a(nbasis, nbasis, lmto_product_max_gf_moment + 1), &
               right_gr(nbasis, nbasis, lmto_product_max_gf_moment + 1), &
               right_ga(nbasis, nbasis, lmto_product_max_gf_moment + 1), &
               right_a(nbasis, nbasis, lmto_product_max_gf_moment + 1))
      raw = cmplx(0.0_rp, 0.0_rp, rp)

      ! Fixed composite Simpson quadrature is intentional.  The GF route is
      ! a correctness/cross-validation backend; its resolution is an explicit
      ! request parameter rather than an inherited DOS/GF mesh setting.
      do ie = 1, ne
         energy = energy_min + real(ie - 1, rp)*step
         if (ie == 1 .or. ie == ne) then
            quadrature_weight = 1.0_rp
         else if (mod(ie, 2) == 0) then
            quadrature_weight = 4.0_rp
         else
            quadrature_weight = 2.0_rp
         end if
         quadrature_weight = quadrature_weight*step/3.0_rp
         fermi_weight = lr_fermi_dirac_occupation(energy, left_state%fermi_level, left_state%temperature)

         do ik = 1, left_state%nk
            call build_weighted_resolvent(left_state, ik, cmplx(energy, integration_eta, rp), left_gr)
            call build_weighted_resolvent(left_state, ik, cmplx(energy, -integration_eta, rp), left_ga)
            call build_weighted_resolvent(right_state, ik, cmplx(energy, integration_eta, rp), right_gr)
            call build_weighted_resolvent(right_state, ik, cmplx(energy, -integration_eta, rp), right_ga)
            left_a = cmplx(0.0_rp, 1.0_rp/(2.0_rp*response_angular_pi), rp)*(left_gr - left_ga)
            right_a = cmplx(0.0_rp, 1.0_rp/(2.0_rp*response_angular_pi), rp)*(right_gr - right_ga)

            do ifrequency = 1, nfrequency
               call build_weighted_resolvent(right_state, ik, &
                  cmplx(energy + request%frequencies(ifrequency), request%eta, rp), right_gr)
               call build_weighted_resolvent(left_state, ik, &
                  cmplx(energy - request%frequencies(ifrequency), -request%eta, rp), left_ga)
               call accumulate_gf_bubble(space, vertices, left_a, right_a, right_gr, left_ga, &
                  fermi_weight*quadrature_weight*2.0_rp*left_state%k_weights(ik)/weight_sum, &
                  raw(:, :, ifrequency))
            end do
         end do
      end do

      allocate(result%susceptibility(space%ndim, space%ndim, nfrequency), result%frequencies(nfrequency))
      do ifrequency = 1, nfrequency
         call response_raw_to_canonical(space, raw(:, :, ifrequency), result%susceptibility(:, :, ifrequency))
      end do
      result%q = request%q
      result%frequencies = request%frequencies
      result%eta = request%eta
      if (channel_kind == 1) then
         result%channel = lr_channel_plus
      else
         result%channel = lr_channel_minus
      end if
      result%response_representation = 'LR-04 canonical right-weighted B=chi_raw*W'
      write (result%response_space_metadata, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,es12.4)') 'ndim=', space%ndim, &
         ' nsite=', space%nsite, ' radial_points=', space%npoint, ' channels=', space%nchannel, &
         ' integration_points=', ne, ' integration_eta=', integration_eta
   end subroutine evaluate_lr_gf_susceptibility

   !> Build G^(p)(z)=sum_n eps_n^p |n><n|/(z-eps_n), through the requested
   !> output order.  The p>0 moments are required because the LMTO transition
   !> vertex is polynomial in each endpoint energy.  This remains a resolvent evaluation;
   !> no LR-06 transition or susceptibility accumulator is called.
   subroutine build_weighted_resolvent(state, ik, z, weighted_green)
      type(lr_electronic_state), intent(in) :: state
      integer, intent(in) :: ik
      complex(rp), intent(in) :: z
      complex(rp), intent(out) :: weighted_green(:, :, :)

      integer :: power, ib, i, j
      complex(rp) :: factor

      if (size(weighted_green, 1) /= state%nbasis .or. size(weighted_green, 2) /= state%nbasis .or. &
          size(weighted_green, 3) < 1 .or. size(weighted_green, 3) > lmto_product_max_gf_moment + 1) then
         error stop 'build_weighted_resolvent: output shape mismatch'
      end if
      weighted_green = cmplx(0.0_rp, 0.0_rp, rp)
      do power = 0, size(weighted_green, 3) - 1
         do ib = 1, state%nbands
            factor = cmplx(energy_power(state%eigenvalues(ib, ik), power), 0.0_rp, rp)/ &
                     (z - state%eigenvalues(ib, ik))
            do j = 1, state%nbasis
               do i = 1, state%nbasis
                  weighted_green(i, j, power + 1) = weighted_green(i, j, power + 1) + &
                     factor*state%eigenvectors(i, ib, ik)*conjg(state%eigenvectors(j, ib, ik))
               end do
            end do
         end do
      end do
   end subroutine build_weighted_resolvent

   subroutine accumulate_gf_bubble(space, vertices, left_a, right_a, right_gr, left_ga, scale, raw)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: vertices(:, :, :, :), left_a(:, :, :), right_a(:, :, :), &
         right_gr(:, :, :), left_ga(:, :, :)
      real(rp), intent(in) :: scale
      complex(rp), intent(inout) :: raw(:, :)

      integer :: component_i, component_j, left_power_i, right_power_i, left_power_j, right_power_j
      integer :: combined_left, combined_right, i, j
      complex(rp) :: temporary(size(left_a, 1), size(left_a, 2))
      complex(rp) :: vertex_adjoint(size(left_a, 1), size(left_a, 2))

      if (size(raw, 1) /= space%ndim .or. size(raw, 2) /= space%ndim) then
         error stop 'accumulate_gf_bubble: response matrix shape mismatch'
      end if

      do component_i = 1, size(vertices, 3)
         call lmto_product_branch_powers(component_i, left_power_i, right_power_i)
         if (maxval(abs(vertices(:, :, component_i, :))) == 0.0_rp) cycle
         do component_j = 1, size(vertices, 3)
            call lmto_product_branch_powers(component_j, left_power_j, right_power_j)
            if (maxval(abs(vertices(:, :, component_j, :))) == 0.0_rp) cycle
            combined_left = left_power_i + left_power_j + 1
            combined_right = right_power_i + right_power_j + 1

            ! First Kubo term: f(E) Tr[A_L(E) V_I G_R^R(E+omega) V_J^dagger].
            do i = 1, space%ndim
               if (maxval(abs(vertices(:, :, component_i, i))) == 0.0_rp) cycle
               temporary = matmul(left_a(:, :, combined_left), vertices(:, :, component_i, i))
               temporary = matmul(temporary, right_gr(:, :, combined_right))
               do j = 1, space%ndim
                  if (maxval(abs(vertices(:, :, component_j, j))) == 0.0_rp) cycle
                  raw(i, j) = raw(i, j) + scale*sum(temporary*conjg(vertices(:, :, component_j, j)))
               end do
            end do

            ! Second Kubo term: f(E) Tr[A_R(E) V_J^dagger G_L^A(E-omega) V_I].
            do j = 1, space%ndim
               if (maxval(abs(vertices(:, :, component_j, j))) == 0.0_rp) cycle
               vertex_adjoint = conjg(transpose(vertices(:, :, component_j, j)))
               temporary = matmul(right_a(:, :, combined_right), vertex_adjoint)
               temporary = matmul(temporary, left_ga(:, :, combined_left))
               do i = 1, space%ndim
                  if (maxval(abs(vertices(:, :, component_i, i))) == 0.0_rp) cycle
                  raw(i, j) = raw(i, j) + scale*sum(temporary*transpose(vertices(:, :, component_i, i)))
               end do
            end do
         end do
      end do
   end subroutine accumulate_gf_bubble

   !> V_I(E_L,E_R)=sum_{branch} E_L^p E_R^q V_I^(p,q), with the
   !> authoritative six-branch map supplying (p,q).
   subroutine build_gf_vertex_tensor(space, radial_bases, operator_matrix, capabilities, vertices)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      complex(rp), intent(in) :: operator_matrix(2, 2)
      type(pauli_vertex_capabilities), intent(in) :: capabilities
      complex(rp), allocatable, intent(out) :: vertices(:, :, :, :)

      type(response_super_index) :: item
      integer :: nbasis, norb, flat, p, q, component, branch, isite, iorb, jorb, ispin, jspin
      integer :: orbital_l, orbital_lp, left_offset, right_offset
      real(rp) :: radial_product

      if (size(radial_bases) /= space%nsite) error stop 'build_gf_vertex_tensor: radial/site mismatch'
      call radial_bases(1)%require_supported(capabilities%reciprocal_mode, capabilities%hamiltonian_order, &
         capabilities%orthogonal, capabilities%collinear, capabilities%has_soc, capabilities%has_extra_operator)
      norb = (radial_bases(1)%lmax + 1)**2
      nbasis = 2*norb*space%nsite
      allocate(vertices(nbasis, nbasis, lmto_product_nbranch, space%ndim))
      vertices = cmplx(0.0_rp, 0.0_rp, rp)

      do flat = 1, space%ndim
         call response_unflatten_superindex(flat, space%nsite, space%response_lmax, space%npoint, &
            space%nchannel, item)
         if (item%channel /= 1) error stop 'build_gf_vertex_tensor: one circular channel is required'
         isite = item%site
         left_offset = (isite - 1)*2*norb
         right_offset = left_offset
         do branch = 1, lmto_product_nbranch
            call lmto_product_branch_powers(branch, p, q)
            component = branch
               do iorb = 1, norb
                  orbital_l = lmto_orbital_l(iorb)
                  do jorb = 1, norb
                     orbital_lp = lmto_orbital_l(jorb)
                     do ispin = 1, 2
                        do jspin = 1, 2
                           call lmto_product_second_order_radial_branch(radial_bases(isite), item%radial_point, &
                              orbital_l, orbital_lp, ispin, jspin, branch, radial_product)
                           vertices(left_offset + (ispin - 1)*norb + iorb, right_offset + (jspin - 1)*norb + jorb, &
                              component, flat) = operator_matrix(ispin, jspin)*radial_product* &
                              response_gaunt(orbital_l, orbital_m(iorb, orbital_l), orbital_lp, &
                              orbital_m(jorb, orbital_lp), item%response_l, item%response_m)
                        end do
                     end do
                  end do
               end do
         end do
      end do
   end subroutine build_gf_vertex_tensor

   pure integer function orbital_m(iorb, l) result(m)
      integer, intent(in) :: iorb, l
      m = iorb - l*l - l - 1
   end function orbital_m

   pure real(rp) function energy_power(energy, power) result(value)
      real(rp), intent(in) :: energy
      integer, intent(in) :: power

      value = lmto_product_energy_power(energy, power)
   end function energy_power

   subroutine validate_gf_inputs(space, radial_bases, left_state, right_state, request, channel_kind)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(lr_electronic_state), intent(in) :: left_state, right_state
      type(lr_gf_susceptibility_request), intent(in) :: request
      integer, intent(out) :: channel_kind

      type(pauli_vertex_capabilities) :: capabilities
      real(rp) :: expected_k(3), scale, expected_occupation
      integer :: ik, ib, expected_nbasis

      if (.not. allocated(request%frequencies) .or. size(request%frequencies) < 1) then
         error stop 'evaluate_lr_gf_susceptibility: at least one frequency is required'
      end if
      if (request%eta <= 0.0_rp) error stop 'evaluate_lr_gf_susceptibility: retarded eta must be positive'
      if (request%integration_points < 3 .or. mod(request%integration_points, 2) == 0) then
         error stop 'evaluate_lr_gf_susceptibility: integration_points must be odd and at least three'
      end if
      if (request%energy_margin <= 0.0_rp) error stop 'evaluate_lr_gf_susceptibility: energy_margin must be positive'
      select case (trim(request%channel))
      case ('plus', 'chi_plus')
         channel_kind = 1
      case ('minus', 'chi_minus')
         channel_kind = 2
      case default
         error stop 'evaluate_lr_gf_susceptibility: channel must be chi_plus or chi_minus'
      end select
      if (space%nchannel /= 1) error stop 'evaluate_lr_gf_susceptibility: one circular channel is required'
      if (size(radial_bases) /= space%nsite) error stop 'evaluate_lr_gf_susceptibility: radial/site mismatch'
      expected_nbasis = 2*(radial_bases(1)%lmax + 1)**2*space%nsite
      if (left_state%nbasis /= expected_nbasis .or. right_state%nbasis /= expected_nbasis) then
         error stop 'evaluate_lr_gf_susceptibility: eigensystem basis is incompatible with radial basis'
      end if
      if (left_state%nk /= right_state%nk .or. left_state%nbands /= right_state%nbands .or. &
          left_state%nbasis /= right_state%nbasis) then
         error stop 'evaluate_lr_gf_susceptibility: endpoint dimensions differ'
      end if
      scale = max(1.0_rp, abs(left_state%fermi_level), abs(right_state%fermi_level), &
         abs(left_state%temperature), abs(right_state%temperature), abs(left_state%energy_zero), &
         abs(right_state%energy_zero))
      if (abs(left_state%fermi_level - right_state%fermi_level) > state_match_tolerance*scale .or. &
          abs(left_state%temperature - right_state%temperature) > state_match_tolerance*scale .or. &
          abs(left_state%energy_zero - right_state%energy_zero) > state_match_tolerance*scale) then
         error stop 'evaluate_lr_gf_susceptibility: EF/temperature/energy-zero provenance differs'
      end if
      if (trim(left_state%reciprocal_mode) /= 'ham_only' .or. trim(right_state%reciprocal_mode) /= 'ham_only' .or. &
          trim(left_state%hamiltonian_order) /= 'second' .or. trim(right_state%hamiltonian_order) /= 'second' .or. &
          .not. left_state%orthogonal .or. .not. right_state%orthogonal .or. .not. left_state%collinear .or. &
          .not. right_state%collinear .or. left_state%has_soc .or. right_state%has_soc .or. &
          left_state%has_extra_operator .or. right_state%has_extra_operator) then
         error stop 'evaluate_lr_gf_susceptibility: representation is outside the LR-GF-02 baseline'
      end if
      do ik = 1, left_state%nk
         expected_k = fold_fractional_kpoint(left_state%k_points(:, ik) + request%q)
         if (maxval(abs(expected_k - right_state%k_points(:, ik))) > endpoint_match_tolerance) then
            error stop 'evaluate_lr_gf_susceptibility: k+q endpoint is not the exact folded endpoint'
         end if
         do ib = 1, left_state%nbands
            expected_occupation = lr_fermi_dirac_occupation(left_state%eigenvalues(ib, ik), &
               left_state%fermi_level, left_state%temperature)
            if (abs(expected_occupation - left_state%occupations(ib, ik)) > occupation_match_tolerance) then
               error stop 'evaluate_lr_gf_susceptibility: occupations are not the certified Fermi snapshot'
            end if
            expected_occupation = lr_fermi_dirac_occupation(right_state%eigenvalues(ib, ik), &
               right_state%fermi_level, right_state%temperature)
            if (abs(expected_occupation - right_state%occupations(ib, ik)) > occupation_match_tolerance) then
               error stop 'evaluate_lr_gf_susceptibility: endpoint occupations are not the certified Fermi snapshot'
            end if
         end do
      end do
      capabilities = pauli_vertex_capabilities()
      call require_gf_response_mesh(space, radial_bases, capabilities)
   end subroutine validate_gf_inputs

   subroutine require_gf_response_mesh(space, radial_bases, capabilities)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(pauli_vertex_capabilities), intent(in) :: capabilities
      integer :: isite
      real(rp), parameter :: mesh_tolerance = 2.0e-13_rp

      if (space%ndim < 1 .or. .not. allocated(space%radius)) then
         error stop 'evaluate_lr_gf_susceptibility: response space is not initialized'
      end if
      do isite = 1, size(radial_bases)
         call radial_bases(isite)%require_supported(capabilities%reciprocal_mode, capabilities%hamiltonian_order, &
            capabilities%orthogonal, capabilities%collinear, capabilities%has_soc, capabilities%has_extra_operator)
         if (radial_bases(isite)%nspin /= 2 .or. radial_bases(isite)%npoint /= space%npoint .or. &
             .not. allocated(radial_bases(isite)%rofi) .or. .not. allocated(radial_bases(isite)%channel_present)) then
            error stop 'evaluate_lr_gf_susceptibility: incomplete radial basis provenance'
         end if
         if (radial_bases(isite)%lmax /= radial_bases(1)%lmax .or. &
             size(radial_bases(isite)%channel_present, 1) /= radial_bases(isite)%lmax + 1 .or. &
             size(radial_bases(isite)%channel_present, 2) /= 2 .or. &
             .not. all(radial_bases(isite)%channel_present)) then
            error stop 'evaluate_lr_gf_susceptibility: radial channel shape or completeness is invalid'
         end if
         if (abs(radial_bases(isite)%mesh_a - space%a) > mesh_tolerance*max(1.0_rp, abs(space%a)) .or. &
             abs(radial_bases(isite)%mesh_b - space%b) > mesh_tolerance*max(1.0_rp, abs(space%b)) .or. &
             maxval(abs(radial_bases(isite)%rofi - space%radius)) > mesh_tolerance* &
             max(1.0_rp, maxval(abs(space%radius)))) then
            error stop 'evaluate_lr_gf_susceptibility: radial mesh provenance differs from response space'
         end if
      end do
   end subroutine require_gf_response_mesh

   pure function fold_fractional_kpoint(k_point) result(folded)
      real(rp), intent(in) :: k_point(3)
      real(rp) :: folded(3)

      folded = k_point - floor(k_point + 0.5_rp)
   end function fold_fractional_kpoint

end module lr_gf_susceptibility_mod
