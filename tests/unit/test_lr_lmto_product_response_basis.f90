!------------------------------------------------------------------------------
! TDVK-02R0 LMTO radial-product response-basis closure audit.
!
! This is an audit-only executable.  It deliberately reconstructs the
! candidate basis independently from the live LR-05 point-vector evaluator;
! no production response path is changed.
!------------------------------------------------------------------------------
program test_lr_lmto_product_response_basis
   use precision_mod, only: rp
   use, intrinsic :: ieee_arithmetic
   use basis_mod, only: basis_init
   use logger_mod, only: g_logger
   use math_mod, only: init_math_operators
   use self_mod, only: legacy_radial_fixture
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use response_basis_mapping_mod, only: response_super_index, response_unflatten_superindex
   use lr_response_space_mod, only: response_space_layout, response_vector_inner_product, response_vector_norm
   use lr_lmto_endpoint_branches_mod, only: lmto_product_nbranch, lmto_product_branch_powers, &
      lmto_product_energy_power
   use lr_pauli_transition_vertex_mod, only: pauli_vertex_capabilities, pauli_endpoint_state, &
      pauli_sigma_plus_matrix, pauli_sigma_minus_matrix, evaluate_pauli_transition_vertex
   implicit none

   integer, parameter :: nr = 51
   integer, parameter :: norb_spd = 9
   real(rp), parameter :: mesh_a = 0.03_rp, mesh_b = 0.10_rp, nuclear_z = 1.0_rp
   real(rp), parameter :: closure_tolerance = 1.0e-10_rp
   real(rp) :: radius(nr)
   type(lmto_radial_basis) :: radial_sp, radial_spd, radial_fe
   type(response_space_layout) :: space_sp, space_spd, space_fe
   type(pauli_vertex_capabilities) :: capabilities
   logical :: failed
   character(len=512) :: fe_file

   type :: candidate_descriptor
      integer :: l = 0
      integer :: lp = 0
      integer :: branch = 0
   end type candidate_descriptor

   call basis_init(2)
   call init_math_operators()
   call g_logger%init()
   call build_mesh(radius)
   call setup_radial_basis(radial_sp, radius, 1)
   call setup_radial_basis(radial_spd, radius, 2)
   capabilities = pauli_vertex_capabilities()
   failed = .false.

   call space_sp%initialize(1, 2, radius, mesh_a, mesh_b, 1)
   call space_spd%initialize(1, 4, radius, mesh_a, mesh_b, 1)

   call inventory_and_spectra(radial_sp, space_sp, 1, failed)
   call inventory_and_spectra(radial_spd, space_spd, 2, failed)
   call scaled_svd_audit(radial_spd, space_spd, failed)
   if (command_argument_count() > 0) then
      call get_command_argument(1, fe_file)
      call load_radial_basis_dump(trim(fe_file), radial_fe)
      call space_fe%initialize(1, 4, radial_fe%rofi, radial_fe%mesh_a, radial_fe%mesh_b, 1)
      call scaled_svd_audit(radial_fe, space_fe, failed, 'accepted Fe')
   end if
   call energy_affine_oracle(radial_spd, failed)
   call gf_component_oracle(radial_spd, failed)
   call lr05_closure_oracle(radial_spd, space_spd, capabilities, failed)

   if (failed) then
      write (*, '(a)') 'UnitLrLmtoProductResponseBasis: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'UnitLrLmtoProductResponseBasis: PASS (inventory, Gram/SVD spectra, LR-05 closure, affine, LR-GF-02)'

contains

   subroutine build_mesh(mesh)
      real(rp), intent(out) :: mesh(:)
      integer :: ir

      mesh(1) = 0.0_rp
      do ir = 2, size(mesh)
         mesh(ir) = mesh_b*(exp(mesh_a*real(ir - 1, rp)) - 1.0_rp)
      end do
   end subroutine build_mesh

   subroutine setup_radial_basis(basis, mesh, basis_lmax)
      type(lmto_radial_basis), intent(out) :: basis
      real(rp), intent(in) :: mesh(:)
      integer, intent(in) :: basis_lmax
      real(rp) :: potential(size(mesh)), energy
      real(rp), allocatable :: g(:, :), gp(:, :), gpp(:, :)
      real(rp), allocatable :: gpack(:), gdotpack(:), gddotpack(:)
      integer :: ispin, l

      potential = 0.0_rp
      call basis%initialize(size(mesh), basis_lmax, 2)
      do ispin = 1, 2
         do l = 0, basis_lmax
            call legacy_radial_fixture(nuclear_z, l, mesh_a, mesh_b, mesh, potential, energy, g, gp, gpp, &
               2.0_rp)
            gpack = reshape(g, [2*size(mesh)])
            gdotpack = reshape(gp, [2*size(mesh)])
            gddotpack = reshape(gpp, [2*size(mesh)])
            call basis%capture_channel(l, ispin, energy, mesh, potential, mesh_a, mesh_b, nuclear_z, &
               gpack, gdotpack, gddotpack, energy)
         end do
      end do
   end subroutine setup_radial_basis

   subroutine enumerate_candidates(lmax, response_l, candidates)
      integer, intent(in) :: lmax, response_l
      type(candidate_descriptor), allocatable, intent(out) :: candidates(:)
      integer :: l, lp, branch, count, k

      count = 0
      do l = 0, lmax
         do lp = 0, lmax
            if (allowed_pair(l, lp, response_l)) count = count + lmto_product_nbranch
         end do
      end do
      allocate(candidates(count))
      k = 0
      do l = 0, lmax
         do lp = 0, lmax
            if (.not. allowed_pair(l, lp, response_l)) cycle
            do branch = 1, lmto_product_nbranch
               k = k + 1
               candidates(k) = candidate_descriptor(l, lp, branch)
            end do
         end do
      end do
   end subroutine enumerate_candidates

   pure logical function allowed_pair(l, lp, response_l) result(is_allowed)
      integer, intent(in) :: l, lp, response_l

      is_allowed = abs(l - lp) <= response_l .and. response_l <= l + lp .and. &
         mod(l + lp + response_l, 2) == 0
   end function allowed_pair

   subroutine inventory_and_spectra(radial, space, lmax, failed)
      type(lmto_radial_basis), intent(in) :: radial
      type(response_space_layout), intent(in) :: space
      integer, intent(in) :: lmax
      logical, intent(inout) :: failed
      type(candidate_descriptor), allocatable :: candidates(:)
      integer :: response_l, response_m, channel, ordered_count, block_dimension
      integer :: total_count, all_m_count
      character(len=8) :: basis_label

      if (lmax == 1) then
         basis_label = 'sp'
      else
         basis_label = 'spd'
      end if
      total_count = 0
      write (*, '(a,a)') 'PRODUCT INVENTORY basis=', trim(basis_label)
      do response_l = 0, 2*lmax
         call enumerate_candidates(lmax, response_l, candidates)
         ordered_count = size(candidates)/lmto_product_nbranch
         block_dimension = size(candidates)
         all_m_count = (2*response_l + 1)*block_dimension
         total_count = total_count + all_m_count
         write (*, '(a,i0,a,i0,a,i0,a,i0)') '  L=', response_l, ' ordered_pairs=', ordered_count, &
            ' block_dimension=', block_dimension, ' all_M_count=', all_m_count
         deallocate(candidates)
      end do
      write (*, '(a,i0)') '  total_candidate_count=', total_count
      if (lmax == 1 .and. total_count /= 78) failed = .true.
      if (lmax == 2 .and. total_count /= 348) then
         write (*, '(a,i0,a)') '  REQUIRED spd count 348, IMPLEMENTED ', total_count, '; STOP'
         failed = .true.
         return
      end if

      do channel = 1, 2
         do response_l = 0, 2*lmax
            do response_m = -response_l, response_l
               call gram_block_audit(radial, space, lmax, response_l, response_m, channel, failed)
            end do
         end do
      end do
   end subroutine inventory_and_spectra

   subroutine scaled_svd_audit(radial, space, failed, audit_label)
      type(lmto_radial_basis), intent(in) :: radial
      type(response_space_layout), intent(in) :: space
      logical, intent(inout) :: failed
      character(len=*), intent(in), optional :: audit_label
      type(candidate_descriptor), allocatable :: candidates(:)
      complex(rp), allocatable :: weighted_basis(:, :), u(:, :), vt(:, :), work(:)
      complex(rp) :: work_query(1)
      real(rp), allocatable :: column_norm(:), singular_values(:), rwork(:)
      real(rp) :: tau1, sigma_min, condition_estimate
      integer :: response_l, channel, spin_left, spin_right, nr_local, ncolumn, nsv, lwork, info
      integer :: ir, k, rank1, rank10, rank100
      character(len=5) :: channel_label
      character(len=64) :: label
      external :: zgesvd

      nr_local = space%npoint
      label = 'spd fixture'
      if (present(audit_label)) label = audit_label
      write (*, '(a,a,a)') 'SCALED SVD audit basis=', trim(label), &
         ' (direct LAPACK zgesvd; no Gram eigensolve)'
      do channel = 1, 2
         if (channel == 1) then
            spin_left = 1
            spin_right = 2
            channel_label = 'plus'
         else
            spin_left = 2
            spin_right = 1
            channel_label = 'minus'
         end if
         do response_l = 0, space%response_lmax
            call enumerate_candidates(radial%lmax, response_l, candidates)
            ncolumn = size(candidates)
            nsv = min(nr_local, ncolumn)
            allocate(weighted_basis(nr_local, ncolumn), column_norm(ncolumn), singular_values(nsv), &
               u(1, 1), vt(1, 1), rwork(max(1, 5*nsv)))
            weighted_basis = cmplx(0.0_rp, 0.0_rp, rp)
            do k = 1, ncolumn
               do ir = 1, nr_local
                  weighted_basis(ir, k) = cmplx(radial_product(radial, ir, candidates(k)%l, &
                     candidates(k)%lp, spin_left, spin_right, candidates(k)%branch), 0.0_rp, rp)
               end do
            end do
            do k = 1, ncolumn
               column_norm(k) = sqrt(sum(space%radial_weights*real(weighted_basis(:, k)*conjg(weighted_basis(:, k)), rp)))
               if (.not. ieee_is_finite(column_norm(k)) .or. column_norm(k) <= tiny(1.0_rp)) then
                  write (*, '(a,i0,a,i0,a,es12.4)') 'SVD invalid column norm L=', response_l, ' column=', k, &
                     ' norm=', column_norm(k)
                  failed = .true.
               else
                  weighted_basis(:, k) = sqrt(space%radial_weights)*weighted_basis(:, k)/column_norm(k)
               end if
            end do
            if (failed) then
               deallocate(candidates, weighted_basis, column_norm, singular_values, u, vt, rwork)
               return
            end if
            call zgesvd('N', 'N', nr_local, ncolumn, weighted_basis, nr_local, singular_values, u, 1, vt, 1, &
               work_query, -1, rwork, info)
            lwork = max(1, nint(real(work_query(1), rp)))
            allocate(work(lwork))
            call zgesvd('N', 'N', nr_local, ncolumn, weighted_basis, nr_local, singular_values, u, 1, vt, 1, &
               work, lwork, rwork, info)
            deallocate(work)
            if (info /= 0) then
               write (*, '(a,i0)') 'SVD zgesvd info=', info
               failed = .true.
               deallocate(candidates, weighted_basis, column_norm, singular_values, u, vt, rwork)
               return
            end if
            tau1 = real(max(nr_local, ncolumn), rp)*epsilon(1.0_rp)*singular_values(1)
            rank1 = count(singular_values > tau1)
            rank10 = count(singular_values > 10.0_rp*tau1)
            rank100 = count(singular_values > 100.0_rp*tau1)
            sigma_min = singular_values(nsv)
            condition_estimate = singular_values(1)/sigma_min
            write (*, '(a,1x,a,1x,a,i0,1x,a,i0,1x,a,es16.8,1x,a,es16.8,1x,a,es16.8,1x,a,es16.8)') &
               'SVD', trim(channel_label), 'L=', response_l, 'ncand=', ncolumn, 'dmin=', minval(column_norm), &
               'dmax=', maxval(column_norm), 'sigma_max=', singular_values(1), 'sigma_min=', sigma_min
            write (*, '(a,es16.8,1x,a,es16.8,1x,a,i0,1x,a,i0,1x,a,i0)') '  condition=', condition_estimate, &
               'tau1=', tau1, 'rank_tau1=', rank1, 'rank_tau10=', rank10, 'rank_tau100=', rank100
            write (*, '(a)', advance='no') '  singular_values:'
            do k = 1, nsv
               write (*, '(1x,es20.12)', advance='no') singular_values(k)
            end do
            write (*, *)
            if (rank1 /= rank10 .or. rank1 /= rank100) then
               write (*, '(a)') '  REVIEW REQUIRED: retained rank changes across sensitivity thresholds.'
            end if
            deallocate(candidates, weighted_basis, column_norm, singular_values, u, vt, rwork)
         end do
      end do
   end subroutine scaled_svd_audit

   subroutine load_radial_basis_dump(filename, basis)
      character(len=*), intent(in) :: filename
      type(lmto_radial_basis), intent(out) :: basis
      character(len=256) :: line
      integer :: unit, ios, npoint, lmax, ir, l, ispin, ir_read, l_read, spin_read
      real(rp) :: mesh_a_read, mesh_b_read, radius_read, phi_read, dot_read, enu_read

      open (newunit=unit, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) error stop 'load_radial_basis_dump: could not open radial extraction'
      read (unit, '(a)', iostat=ios) line
      read (unit, '(a)', iostat=ios) line
      read (line(index(line, '=') + 1:), *) npoint
      read (unit, '(a)', iostat=ios) line
      read (line(index(line, '=') + 1:), *) lmax
      read (unit, '(a)', iostat=ios) line
      read (line(index(line, '=') + 1:), *) mesh_a_read
      read (unit, '(a)', iostat=ios) line
      read (line(index(line, '=') + 1:), *) mesh_b_read
      read (unit, '(a)', iostat=ios) line
      call basis%initialize(npoint, lmax, 2)
      basis%mesh_a = mesh_a_read
      basis%mesh_b = mesh_b_read
      do ispin = 1, 2
         do l = 0, lmax
            do ir = 1, npoint
               read (unit, *, iostat=ios) ir_read, l_read, spin_read, radius_read, phi_read, dot_read, enu_read
               if (ios /= 0 .or. ir_read /= ir .or. l_read /= l .or. spin_read /= ispin) then
                  error stop 'load_radial_basis_dump: malformed radial extraction'
               end if
               basis%rofi(ir) = radius_read
               basis%phi_large(ir, l + 1, ispin) = phi_read
               basis%phidot_large(ir, l + 1, ispin) = dot_read
               basis%enu_work(l + 1, ispin) = enu_read
               basis%enu_radial(l + 1, ispin) = enu_read
            end do
         end do
      end do
      basis%channel_present = .true.
      close (unit)
   end subroutine load_radial_basis_dump

   subroutine build_candidate_block(radial, space, response_l, response_m, channel, candidates, basis)
      type(lmto_radial_basis), intent(in) :: radial
      type(response_space_layout), intent(in) :: space
      integer, intent(in) :: response_l, response_m, channel
      type(candidate_descriptor), intent(in) :: candidates(:)
      complex(rp), intent(out) :: basis(:, :)
      type(response_super_index) :: item
      integer :: flat, k, spin_left, spin_right

      if (size(basis, 1) /= space%ndim .or. size(basis, 2) /= size(candidates)) then
         error stop 'build_candidate_block: shape mismatch'
      end if
      if (channel == 1) then
         spin_left = 1
         spin_right = 2
      else
         spin_left = 2
         spin_right = 1
      end if
      basis = cmplx(0.0_rp, 0.0_rp, rp)
      do flat = 1, space%ndim
         call response_unflatten_superindex(flat, space%nsite, space%response_lmax, space%npoint, &
            space%nchannel, item)
         if (item%site /= 1 .or. item%response_l /= response_l .or. item%response_m /= response_m .or. &
             item%channel /= 1) cycle
         do k = 1, size(candidates)
            basis(flat, k) = cmplx(radial_product(radial, item%radial_point, candidates(k)%l, candidates(k)%lp, &
               spin_left, spin_right, candidates(k)%branch), 0.0_rp, rp)
         end do
      end do
   end subroutine build_candidate_block

   subroutine factor_block(radial, space, lmax, response_l, response_m, channel, candidates, basis, eigenvalues, &
                           eigenvectors, rank, rank_tolerance, normalize_columns)
      type(lmto_radial_basis), intent(in) :: radial
      type(response_space_layout), intent(in) :: space
      integer, intent(in) :: lmax, response_l, response_m, channel
      type(candidate_descriptor), intent(in) :: candidates(:)
      complex(rp), allocatable, intent(out) :: basis(:, :), eigenvectors(:, :)
      real(rp), allocatable, intent(out) :: eigenvalues(:)
      integer, intent(out) :: rank
      real(rp), intent(out) :: rank_tolerance
      logical, intent(in) :: normalize_columns
      complex(rp), allocatable :: gram(:, :)
      real(rp), allocatable :: column_norm(:)
      integer :: i, j

      allocate(basis(space%ndim, size(candidates)), gram(size(candidates), size(candidates)))
      call build_candidate_block(radial, space, response_l, response_m, channel, candidates, basis)
      do i = 1, size(candidates)
         do j = i, size(candidates)
            gram(i, j) = response_vector_inner_product(space, basis(:, i), basis(:, j))
            gram(j, i) = conjg(gram(i, j))
         end do
      end do
      if (normalize_columns) then
         allocate(column_norm(size(candidates)))
         do i = 1, size(candidates)
            column_norm(i) = sqrt(max(real(gram(i, i), rp), tiny(1.0_rp)))
            basis(:, i) = basis(:, i)/column_norm(i)
         end do
         do i = 1, size(candidates)
            do j = i, size(candidates)
               gram(i, j) = response_vector_inner_product(space, basis(:, i), basis(:, j))
               gram(j, i) = conjg(gram(i, j))
            end do
         end do
         deallocate(column_norm)
      end if
      call hermitian_eigenpairs(gram, eigenvalues, eigenvectors)
      rank_tolerance = epsilon(1.0_rp)*real(size(candidates), rp)*maxval(eigenvalues)
      rank = count(eigenvalues > rank_tolerance)
      deallocate(gram)
   end subroutine factor_block

   subroutine gram_block_audit(radial, space, lmax, response_l, response_m, channel, failed)
      type(lmto_radial_basis), intent(in) :: radial
      type(response_space_layout), intent(in) :: space
      integer, intent(in) :: lmax, response_l, response_m, channel
      logical, intent(inout) :: failed
      type(candidate_descriptor), allocatable :: candidates(:)
      complex(rp), allocatable :: basis(:, :), eigenvectors(:, :)
      real(rp), allocatable :: eigenvalues(:)
      real(rp) :: rank_tolerance, minimum_resolved, condition_estimate
      integer :: rank, k
      character(len=5) :: channel_label

      call enumerate_candidates(lmax, response_l, candidates)
      call factor_block(radial, space, lmax, response_l, response_m, channel, candidates, basis, eigenvalues, &
         eigenvectors, rank, rank_tolerance, .false.)
      if (channel == 1) then
         channel_label = 'plus'
      else
         channel_label = 'minus'
      end if
      minimum_resolved = minval(eigenvalues, mask=eigenvalues > rank_tolerance)
      condition_estimate = maxval(eigenvalues)/minimum_resolved
      write (*, '(a,1x,a,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,es12.4,1x,a,es12.4,1x,a,es12.4,1x,a,i0,1x,a,es12.4)') &
         'GRAM', trim(channel_label), 'L=', response_l, 'M=', response_m, 'dim=', size(candidates), &
         'rank=', rank, 'largest=', maxval(eigenvalues), 'smallest_resolved=', minimum_resolved, &
         'condition=', condition_estimate, 'nullity=', size(candidates) - rank, 'rank_tol=', rank_tolerance
      write (*, '(a)', advance='no') '  eigenvalues:'
      do k = 1, size(eigenvalues)
         write (*, '(1x,es14.6)', advance='no') eigenvalues(k)
      end do
      write (*, *)
      if (rank < size(candidates)) then
         write (*, '(a)') '  NOTE numerical nullity is diagnostic only; no production cutoff is selected.'
      end if
      deallocate(candidates, basis, eigenvalues, eigenvectors)
   end subroutine gram_block_audit

   subroutine hermitian_eigenpairs(matrix, eigenvalues, eigenvectors)
      complex(rp), intent(in) :: matrix(:, :)
      real(rp), allocatable, intent(out) :: eigenvalues(:)
      complex(rp), allocatable, intent(out) :: eigenvectors(:, :)
      complex(rp) :: work_query(1)
      complex(rp), allocatable :: work(:)
      real(rp), allocatable :: rwork(:)
      integer :: n, lwork, info
      logical :: halt_divide_by_zero
      external :: zheev

      n = size(matrix, 1)
      if (size(matrix, 2) /= n) error stop 'hermitian_eigenpairs: matrix must be square'
      allocate(eigenvalues(n), eigenvectors(n, n), rwork(max(1, 3*n - 2)))
      eigenvectors = matrix
      call ieee_get_halting_mode(ieee_divide_by_zero, halt_divide_by_zero)
      call ieee_set_halting_mode(ieee_divide_by_zero, .false.)
      call zheev('V', 'U', n, eigenvectors, n, eigenvalues, work_query, -1, rwork, info)
      call ieee_set_flag(ieee_divide_by_zero, .false.)
      lwork = max(1, nint(real(work_query(1), rp)))
      allocate(work(lwork))
      call zheev('V', 'U', n, eigenvectors, n, eigenvalues, work, lwork, rwork, info)
      call ieee_set_flag(ieee_divide_by_zero, .false.)
      call ieee_set_halting_mode(ieee_divide_by_zero, halt_divide_by_zero)
      deallocate(work, rwork)
      if (info /= 0) then
         write (*, '(a,i0)') 'hermitian_eigenpairs: zheev info=', info
         error stop 'hermitian_eigenpairs: zheev failed'
      end if
   end subroutine hermitian_eigenpairs

   pure real(rp) function endpoint_component(radial, ir, l, spin, power) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, spin, power

      if (power == 0) then
         value = radial%phi_large(ir, l + 1, spin) - radial%enu_work(l + 1, spin)* &
            radial%phidot_large(ir, l + 1, spin)
      else
         value = radial%phidot_large(ir, l + 1, spin)
      end if
   end function endpoint_component

   recursive real(rp) function radial_product(radial, ir, l, lp, spin_left, spin_right, branch) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, lp, spin_left, spin_right, branch
      real(rp) :: phi_l, phi_r, dot_l, dot_r, ddot_l, ddot_r, enu_l, enu_r
      real(rp) :: first, second

      if (ir == 1) then
         if (l /= 0 .or. lp /= 0) then
            value = 0.0_rp
            return
         end if
         first = radial_product(radial, 2, l, lp, spin_left, spin_right, branch)
         second = radial_product(radial, 3, l, lp, spin_left, spin_right, branch)
         value = (first*radial%rofi(3)**2 - second*radial%rofi(2)**2)/ &
            (radial%rofi(3)**2 - radial%rofi(2)**2)
         return
      end if
      if (radial%rofi(ir) <= tiny(1.0_rp)) error stop 'radial_product: invalid positive radial point'
      phi_l = radial%phi_large(ir, l + 1, spin_left)
      phi_r = radial%phi_large(ir, lp + 1, spin_right)
      dot_l = radial%phidot_large(ir, l + 1, spin_left)
      dot_r = radial%phidot_large(ir, lp + 1, spin_right)
      ddot_l = radial%phiddot_large(ir, l + 1, spin_left)
      ddot_r = radial%phiddot_large(ir, lp + 1, spin_right)
      enu_l = radial%enu_work(l + 1, spin_left)
      enu_r = radial%enu_work(lp + 1, spin_right)
      select case (branch)
      case (1)
         value = phi_l*phi_r - enu_l*dot_l*phi_r - enu_r*phi_l*dot_r + enu_l*enu_r*dot_l*dot_r + &
            0.5_rp*enu_l**2*ddot_l*phi_r + 0.5_rp*enu_r**2*phi_l*ddot_r
      case (2)
         value = dot_l*phi_r - enu_r*dot_l*dot_r - enu_l*ddot_l*phi_r
      case (3)
         value = phi_l*dot_r - enu_l*dot_l*dot_r - enu_r*phi_l*ddot_r
      case (4)
         value = dot_l*dot_r
      case (5)
         value = 0.5_rp*ddot_l*phi_r
      case (6)
         value = 0.5_rp*phi_l*ddot_r
      case default
         error stop 'radial_product: invalid branch'
      end select
      value = value/radial%rofi(ir)**2
   end function radial_product

   subroutine energy_affine_oracle(radial, failed)
      type(lmto_radial_basis), intent(in) :: radial
      logical, intent(inout) :: failed
      real(rp), parameter :: energies(4) = [-0.37_rp, -0.11_rp, 0.19_rp, 0.43_rp]
      integer :: l, spin, ir, ie
      real(rp) :: live, affine, error, scale, maximum_error, maximum_relative

      maximum_error = 0.0_rp
      maximum_relative = 0.0_rp
      do spin = 1, 2
         do l = 0, radial%lmax
            do ie = 1, size(energies)
               do ir = 1, radial%npoint
                  live = radial%phi_large(ir, l + 1, spin) + (energies(ie) - radial%enu_work(l + 1, spin))* &
                     radial%phidot_large(ir, l + 1, spin)
                  affine = endpoint_component(radial, ir, l, spin, 0) + energies(ie)* &
                     endpoint_component(radial, ir, l, spin, 1)
                  error = abs(live - affine)
                  scale = max(abs(live), abs(affine), epsilon(1.0_rp))
                  maximum_error = max(maximum_error, error)
                  maximum_relative = max(maximum_relative, error/scale)
               end do
            end do
         end do
      end do
      write (*, '(a,es12.4,a,es12.4)') 'AFFINE maximum_abs=', maximum_error, ' maximum_relative=', maximum_relative
      if (maximum_relative > 2.0e-14_rp) failed = .true.
   end subroutine energy_affine_oracle

   subroutine gf_component_oracle(radial, failed)
      type(lmto_radial_basis), intent(in) :: radial
      logical, intent(inout) :: failed
      integer :: component, p, q, l, lp, spin_left, spin_right, ir
      real(rp) :: candidate, gf_component, maximum_error

      maximum_error = 0.0_rp
      do component = 1, lmto_product_nbranch
         call lmto_product_branch_powers(component, p, q)
         if (p > 2 .or. q > 2) failed = .true.
         do spin_left = 1, 2
            do spin_right = 1, 2
               do l = 0, radial%lmax
                  do lp = 0, radial%lmax
                     do ir = 1, radial%npoint
                        candidate = radial_product(radial, ir, l, lp, spin_left, spin_right, component)
                        gf_component = gf_radial_vertex_component(radial, ir, l, lp, spin_left, spin_right, component)
                        maximum_error = max(maximum_error, abs(candidate - gf_component))
                     end do
                  end do
               end do
            end do
         end do
      end do
      write (*, '(a,i0,a,es12.4)') 'GF-02 component mapping count=', lmto_product_nbranch, &
         ' maximum_radial_factor_error=', maximum_error
      if (maximum_error > 0.0_rp) failed = .true.
   end subroutine gf_component_oracle

   real(rp) function gf_radial_vertex_component(radial, ir, l, lp, spin_left, spin_right, branch) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, lp, spin_left, spin_right, branch
      integer :: p, q

      call lmto_product_branch_powers(branch, p, q)
      value = radial_product(radial, ir, l, lp, spin_left, spin_right, branch)
   end function gf_radial_vertex_component

   subroutine lr05_closure_oracle(radial, space, capabilities, failed)
      type(lmto_radial_basis), intent(in) :: radial
      type(response_space_layout), intent(in) :: space
      type(pauli_vertex_capabilities), intent(in) :: capabilities
      logical, intent(inout) :: failed
      integer, parameter :: nstate = 2*norb_spd, npair = 4
      integer, parameter :: left_index(npair) = [1, 3, 5, 9]
      integer, parameter :: right_index(npair) = [2, 8, 12, 18]
      complex(rp) :: hamiltonian(nstate, nstate), eigenvectors(nstate, nstate), operator_matrix(2, 2)
      real(rp) :: eigenvalues(nstate), error, maximum_error
      type(pauli_endpoint_state) :: left_state, right_state
      complex(rp), allocatable :: transition(:)
      integer :: pair, channel, info
      external :: zheev

      hamiltonian = cmplx(0.0_rp, 0.0_rp, rp)
      do pair = 1, nstate
         hamiltonian(pair, pair) = cmplx(-0.85_rp + 0.11_rp*real(pair, rp), 0.0_rp, rp)
         do channel = pair + 1, nstate
            hamiltonian(pair, channel) = cmplx(0.0013_rp*real(pair + channel, rp), &
               0.0007_rp*real(channel - pair, rp), rp)
            hamiltonian(channel, pair) = conjg(hamiltonian(pair, channel))
         end do
      end do
      call diagonalize_fixture(hamiltonian, eigenvalues, eigenvectors, info)
      if (info /= 0) then
         failed = .true.
         return
      end if

      maximum_error = 0.0_rp
      do pair = 1, npair
         call left_state%initialize(eigenvalues(left_index(pair)), eigenvectors(:, left_index(pair)))
         call right_state%initialize(eigenvalues(right_index(pair)), eigenvectors(:, right_index(pair)))
         do channel = 1, 2
            if (channel == 1) then
               operator_matrix = pauli_sigma_plus_matrix()
            else
               operator_matrix = pauli_sigma_minus_matrix()
            end if
            allocate(transition(space%ndim))
            call evaluate_pauli_transition_vertex(space, [radial], left_state, right_state, operator_matrix, &
               capabilities, transition)
            call project_and_report(space, radial, transition, pair, channel, error)
            maximum_error = max(maximum_error, error)
            deallocate(transition)
         end do
      end do
      write (*, '(a,es12.4)') 'LR-05 maximum_relative_closure_error=', maximum_error
      if (maximum_error > closure_tolerance) failed = .true.
   end subroutine lr05_closure_oracle

   subroutine diagonalize_fixture(matrix, eigenvalues, eigenvectors, info)
      complex(rp), intent(in) :: matrix(:, :)
      real(rp), intent(out) :: eigenvalues(:)
      complex(rp), intent(out) :: eigenvectors(:, :)
      integer, intent(out) :: info
      complex(rp) :: work_query(1)
      complex(rp), allocatable :: work(:)
      real(rp), allocatable :: rwork(:)
      integer :: n, lwork
      logical :: halt_divide_by_zero
      external :: zheev

      n = size(matrix, 1)
      eigenvectors = matrix
      allocate(rwork(max(1, 3*n - 2)))
      call ieee_get_halting_mode(ieee_divide_by_zero, halt_divide_by_zero)
      call ieee_set_halting_mode(ieee_divide_by_zero, .false.)
      call zheev('V', 'U', n, eigenvectors, n, eigenvalues, work_query, -1, rwork, info)
      call ieee_set_flag(ieee_divide_by_zero, .false.)
      lwork = max(1, nint(real(work_query(1), rp)))
      allocate(work(lwork))
      call zheev('V', 'U', n, eigenvectors, n, eigenvalues, work, lwork, rwork, info)
      call ieee_set_flag(ieee_divide_by_zero, .false.)
      call ieee_set_halting_mode(ieee_divide_by_zero, halt_divide_by_zero)
      deallocate(work, rwork)
   end subroutine diagonalize_fixture

   subroutine project_and_report(space, radial, transition, pair, channel, relative_error)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial
      complex(rp), intent(in) :: transition(:)
      integer, intent(in) :: pair, channel
      real(rp), intent(out) :: relative_error
      type(candidate_descriptor), allocatable :: candidates(:)
      complex(rp), allocatable :: basis(:, :), eigenvectors(:, :), rhs(:), modal(:), coefficients(:), reconstructed(:)
      real(rp), allocatable :: eigenvalues(:)
      real(rp) :: rank_tolerance, absolute_error
      integer :: response_l, response_m, rank, k

      allocate(reconstructed(space%ndim))
      reconstructed = cmplx(0.0_rp, 0.0_rp, rp)
      do response_l = 0, space%response_lmax
         call enumerate_candidates(radial%lmax, response_l, candidates)
         do response_m = -response_l, response_l
            call factor_block(radial, space, radial%lmax, response_l, response_m, channel, candidates, basis, &
               eigenvalues, eigenvectors, rank, rank_tolerance, .true.)
            allocate(rhs(size(candidates)), modal(size(candidates)), coefficients(size(candidates)))
            do k = 1, size(candidates)
               rhs(k) = response_vector_inner_product(space, basis(:, k), transition)
            end do
            modal = matmul(conjg(transpose(eigenvectors)), rhs)
            coefficients = cmplx(0.0_rp, 0.0_rp, rp)
            do k = 1, size(candidates)
               if (eigenvalues(k) > rank_tolerance) coefficients = coefficients + &
                  eigenvectors(:, k)*modal(k)/eigenvalues(k)
            end do
            reconstructed = reconstructed + matmul(basis, coefficients)
            deallocate(basis, eigenvalues, eigenvectors, rhs, modal, coefficients)
         end do
         deallocate(candidates)
      end do
      absolute_error = response_vector_norm(space, transition - reconstructed)
      relative_error = absolute_error/max(response_vector_norm(space, transition), epsilon(1.0_rp))
      write (*, '(a,i0,a,i0,a,es12.4)') 'LR-05 closure pair=', pair, ' channel=', channel, ' relative_error=', relative_error
      deallocate(reconstructed)
   end subroutine project_and_report

end program test_lr_lmto_product_response_basis
