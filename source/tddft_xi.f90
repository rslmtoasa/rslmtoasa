!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Direct kernel-weighted self-enhancement operator from eigenpairs.
!>
!> This module deliberately accepts weighted electronic-space operators Q_b
!> rather than a response-space kernel.  The k-resolved entry point is used by
!> the production pair-potential workflow: Q_b(k,q) is kept in the active
!> LMTO representation until after its transition matrix elements are formed.
module tddft_xi_mod
   use precision_mod, only: rp
   use response_vertices_mod, only: response_channel, response_transition_vertex, weighted_transition_vertex
   use tddft_chi0_mod, only: tddft_chi0_options, tddft_chi0_metadata, tddft_fermi_occupation, &
      tddft_kB_Ry_per_K, tddft_occupation_kT_floor, tddft_static_divided_difference
   use tddft_transition_engine_mod, only: tddft_transition_engine, pair_operator_vertex_provider, &
      pair_operator_tile_source, make_pair_operator_vertex_provider
   implicit none

   private

   type, public :: tddft_direct_xi_result
      complex(rp), allocatable :: xi(:, :, :)
      type(tddft_chi0_metadata) :: metadata
   end type tddft_direct_xi_result

   public :: build_direct_xi_from_eigenpairs
   public :: build_direct_xi_from_k_dependent_eigenpairs
   public :: build_static_direct_xi_from_k_dependent_eigenpairs
   public :: build_direct_xi_from_operator_source
   public :: build_static_direct_xi_from_operator_source

contains

   !> Dynamic Xi from a producer that supplies exactly one operator tile for
   !> each local k point.  The source remains caller-owned; the transition
   !> provider allocates only its active `(nmat,nmat,nright)` tile.
   subroutine build_direct_xi_from_operator_source(k_weights, eigenvalues_k, eigenvectors_k, eigenvalues_kq, &
      eigenvectors_kq, site_orbital_counts, left_channels, operator_source, omega, options, result)
      real(rp), target, intent(in) :: k_weights(:), eigenvalues_k(:, :), eigenvalues_kq(:, :), omega(:)
      complex(rp), target, intent(in) :: eigenvectors_k(:, :, :), eigenvectors_kq(:, :, :)
      integer, target, intent(in) :: site_orbital_counts(:)
      type(response_channel), target, intent(in) :: left_channels(:)
      class(pair_operator_tile_source), target, intent(inout) :: operator_source
      type(tddft_chi0_options), intent(in) :: options
      type(tddft_direct_xi_result), intent(out) :: result

      integer :: nk, nbands, nspinor, nleft, nright, nw, band_first, band_last, batch_size
      real(rp) :: weight_sum
      type(tddft_transition_engine) :: engine
      type(pair_operator_vertex_provider) :: provider

      nk = size(k_weights); nbands = size(eigenvalues_k, 1); nspinor = 2*sum(site_orbital_counts)
      nleft = size(left_channels); nright = operator_source%channel_dimension(); nw = size(omega)
      call validate_direct_xi_source_inputs(nk, nbands, nspinor, nleft, nright, nw, k_weights, eigenvalues_k, eigenvectors_k, &
         eigenvalues_kq, eigenvectors_kq, options)
      band_first = options%band_first; band_last = options%band_last
      if (band_last == 0) band_last = nbands
      if (band_first < 1 .or. band_last < band_first .or. band_last > nbands) then
         error stop 'build_direct_xi_from_operator_source: invalid selected band window'
      end if
      weight_sum = sum(k_weights)
      allocate(result%xi(nleft, nright, nw)); result%xi = cmplx(0.0_rp, 0.0_rp, rp)
      call set_direct_xi_metadata(result%metadata, options, nk, nbands, band_first, band_last, weight_sum)
      result%metadata%backend = 'direct_xi_operator_source'
      batch_size = min(options%transition_batch_size, (band_last-band_first+1)**2)
      call make_pair_operator_vertex_provider(provider, site_orbital_counts, left_channels, eigenvectors_k, eigenvectors_kq, &
         operator_source)
      call engine%accumulate_dynamic(k_weights, eigenvalues_k, eigenvalues_kq, omega, options%eta, options%fermi_level, &
         options%electronic_temperature, band_first, band_last, options%occupation_prune_tolerance, batch_size, &
         options%use_batched_accumulation, provider, result%xi, result%metadata%vertex_cpu_seconds, &
         result%metadata%transition_preparation_cpu_seconds, result%metadata%denominator_cpu_seconds, &
         result%metadata%accumulation_cpu_seconds)
   end subroutine build_direct_xi_from_operator_source

   !> Static q=0 Xi from the same streamed tile source contract.  As in the
   !> legacy static entry point, no dynamic eta enters its divided difference.
   subroutine build_static_direct_xi_from_operator_source(k_weights, eigenvalues, eigenvectors, site_orbital_counts, &
      left_channels, operator_source, options, result)
      real(rp), target, intent(in) :: k_weights(:), eigenvalues(:, :)
      complex(rp), target, intent(in) :: eigenvectors(:, :, :)
      integer, target, intent(in) :: site_orbital_counts(:)
      type(response_channel), target, intent(in) :: left_channels(:)
      class(pair_operator_tile_source), target, intent(inout) :: operator_source
      type(tddft_chi0_options), intent(in) :: options
      type(tddft_direct_xi_result), intent(out) :: result

      integer :: nk, nbands, nspinor, nleft, nright, band_first, band_last, batch_size
      real(rp) :: weight_sum
      type(tddft_transition_engine) :: engine
      type(pair_operator_vertex_provider) :: provider

      nk = size(k_weights); nbands = size(eigenvalues, 1); nspinor = 2*sum(site_orbital_counts)
      nleft = size(left_channels); nright = operator_source%channel_dimension()
      call validate_direct_xi_source_inputs(nk, nbands, nspinor, nleft, nright, 1, k_weights, eigenvalues, eigenvectors, &
         eigenvalues, eigenvectors, options)
      band_first = options%band_first; band_last = options%band_last
      if (band_last == 0) band_last = nbands
      if (band_first < 1 .or. band_last < band_first .or. band_last > nbands) then
         error stop 'build_static_direct_xi_from_operator_source: invalid selected band window'
      end if
      weight_sum = sum(k_weights)
      allocate(result%xi(nleft, nright, 1)); result%xi = cmplx(0.0_rp, 0.0_rp, rp)
      call set_direct_xi_metadata(result%metadata, options, nk, nbands, band_first, band_last, weight_sum)
      result%metadata%backend = 'static_xi_operator_source'
      result%metadata%frequency_convention = 'static q=0 omega=0 divided difference; no dynamical eta'
      result%metadata%eta = 0.0_rp
      batch_size = min(options%transition_batch_size, (band_last-band_first+1)**2)
      call make_pair_operator_vertex_provider(provider, site_orbital_counts, left_channels, eigenvectors, eigenvectors, &
         operator_source)
      call engine%accumulate_static(k_weights, eigenvalues, options%fermi_level, options%electronic_temperature, &
         band_first, band_last, options%occupation_prune_tolerance, batch_size, provider, result%xi, &
         result%metadata%vertex_cpu_seconds, result%metadata%transition_preparation_cpu_seconds, &
         result%metadata%accumulation_cpu_seconds)
   end subroutine build_static_direct_xi_from_operator_source

   !> Dedicated real static q=0 operator.  This is intentionally separate
   !> from the dynamic builder so a response eta can never enter Xi_static.
   subroutine build_static_direct_xi_from_k_dependent_eigenpairs(k_weights, eigenvalues, eigenvectors, &
      site_orbital_counts, left_channels, weighted_right_operators, options, result)
      real(rp), target, intent(in) :: k_weights(:), eigenvalues(:, :)
      complex(rp), target, intent(in) :: eigenvectors(:, :, :), weighted_right_operators(:, :, :, :)
      integer, target, intent(in) :: site_orbital_counts(:)
      type(response_channel), target, intent(in) :: left_channels(:)
      type(tddft_chi0_options), intent(in) :: options
      type(tddft_direct_xi_result), intent(out) :: result

      integer :: nk, nbands, nspinor, nleft, nright, ik, n, m, band_first, band_last
      real(rp) :: weight_sum, prefactor, factor
      complex(rp), allocatable :: left_vertex(:), right_vertex(:)
      type(tddft_transition_engine) :: engine
      type(pair_operator_vertex_provider) :: provider
      integer :: batch_size

      nk = size(k_weights); nbands = size(eigenvalues, 1); nspinor = 2*sum(site_orbital_counts)
      nleft = size(left_channels); nright = size(weighted_right_operators, 3)
      call validate_direct_xi_inputs(nk, nbands, nspinor, nleft, nright, 1, k_weights, eigenvalues, eigenvectors, &
         eigenvalues, eigenvectors, weighted_right_operators(:, :, :, 1), options)
      if (size(weighted_right_operators, 4) /= nk) then
         error stop 'build_static_direct_xi_from_k_dependent_eigenpairs: operator k dimension is incompatible'
      end if
      band_first = options%band_first; band_last = options%band_last
      if (band_last == 0) band_last = nbands
      if (band_first < 1 .or. band_last < band_first .or. band_last > nbands) then
         error stop 'build_static_direct_xi_from_k_dependent_eigenpairs: invalid selected band window'
      end if
      weight_sum = sum(k_weights)
      allocate(result%xi(nleft, nright, 1)); result%xi = cmplx(0.0_rp, 0.0_rp, rp)
      call set_direct_xi_metadata(result%metadata, options, nk, nbands, band_first, band_last, weight_sum)
      result%metadata%backend = 'static_xi_k_resolved_eigenpairs'
      result%metadata%frequency_convention = 'static q=0 omega=0 divided difference; no dynamical eta'
      result%metadata%eta = 0.0_rp
      batch_size = min(options%transition_batch_size, (band_last-band_first+1)**2)
      call make_pair_operator_vertex_provider(provider, site_orbital_counts, left_channels, eigenvectors, eigenvectors, &
         weighted_right_operators, .true.)
      call engine%accumulate_static(k_weights, eigenvalues, options%fermi_level, options%electronic_temperature, &
         band_first, band_last, options%occupation_prune_tolerance, batch_size, provider, result%xi, &
         result%metadata%vertex_cpu_seconds, result%metadata%transition_preparation_cpu_seconds, &
         result%metadata%accumulation_cpu_seconds)
      if (.false.) then
      allocate(left_vertex(nleft), right_vertex(nright))
      do ik = 1, nk
         prefactor = k_weights(ik)/weight_sum
         do n = band_first, band_last
            do m = band_first, band_last
               factor = tddft_static_divided_difference(eigenvalues(n, ik), eigenvalues(m, ik), options%fermi_level, &
                  options%electronic_temperature)
               if (options%occupation_prune_tolerance > 0.0_rp .and. abs(factor) <= options%occupation_prune_tolerance) cycle
               call build_transition_vertices(left_channels, weighted_right_operators(:, :, :, ik), site_orbital_counts, &
                  eigenvectors(:, n, ik), eigenvectors(:, m, ik), left_vertex, right_vertex)
               result%xi(:, :, 1) = result%xi(:, :, 1) + prefactor*factor*outer_product(left_vertex, right_vertex)
            end do
         end do
      end do
      deallocate(left_vertex, right_vertex)
      end if
   end subroutine build_static_direct_xi_from_k_dependent_eigenpairs

   !> Assemble Xi_ab=sum (f_n-f_m)V_a W_b/(omega+en-em+i eta), where the
   !> weighted right vertex has the required order W_b=<m,k+q|Q_b|n,k>.
   !> `weighted_right_operators(:,:,b)` must be in the same electronic
   !> coefficient representation as the supplied eigenvectors.
   subroutine build_direct_xi_from_eigenpairs(k_weights, eigenvalues_k, eigenvectors_k, eigenvalues_kq, &
      eigenvectors_kq, site_orbital_counts, left_channels, weighted_right_operators, omega, options, result)
      real(rp), target, intent(in) :: k_weights(:), eigenvalues_k(:, :), eigenvalues_kq(:, :), omega(:)
      complex(rp), target, intent(in) :: eigenvectors_k(:, :, :), eigenvectors_kq(:, :, :)
      integer, target, intent(in) :: site_orbital_counts(:)
      type(response_channel), target, intent(in) :: left_channels(:)
      complex(rp), target, intent(in) :: weighted_right_operators(:, :, :)
      type(tddft_chi0_options), intent(in) :: options
      type(tddft_direct_xi_result), intent(out) :: result

      integer :: nk, nbands, nspinor, nleft, nright, nw, ik, n, m, iw, npairs, batch_size
      integer :: band_first, band_last
      real(rp) :: weight_sum, occupation_difference, transition_energy, prefactor, t_start, t_stop
      complex(rp) :: denominator
      complex(rp), allocatable :: left_vertex(:), right_vertex(:)
      complex(rp), allocatable :: left_batch(:, :), right_batch(:, :), weighted_left(:, :), denominator_batch(:)
      real(rp), allocatable :: occupation_batch(:), transition_energy_batch(:)
      type(tddft_transition_engine) :: engine
      type(pair_operator_vertex_provider) :: provider

      nk = size(k_weights); nbands = size(eigenvalues_k, 1); nspinor = 2*sum(site_orbital_counts)
      nleft = size(left_channels); nright = size(weighted_right_operators, 3); nw = size(omega)
      call validate_direct_xi_inputs(nk, nbands, nspinor, nleft, nright, nw, k_weights, eigenvalues_k, eigenvectors_k, &
         eigenvalues_kq, eigenvectors_kq, weighted_right_operators, options)
      band_first = options%band_first; band_last = options%band_last
      if (band_last == 0) band_last = nbands
      if (band_first < 1 .or. band_last < band_first .or. band_last > nbands) then
         error stop 'build_direct_xi_from_eigenpairs: invalid selected band window'
      end if
      weight_sum = sum(k_weights)
      allocate(result%xi(nleft, nright, nw))
      result%xi = cmplx(0.0_rp, 0.0_rp, rp)
      call set_direct_xi_metadata(result%metadata, options, nk, nbands, band_first, band_last, weight_sum)
      batch_size = min(options%transition_batch_size, (band_last-band_first+1)**2)
      call make_pair_operator_vertex_provider(provider, site_orbital_counts, left_channels, eigenvectors_k, eigenvectors_kq, &
         weighted_right_operators)
      call engine%accumulate_dynamic(k_weights, eigenvalues_k, eigenvalues_kq, omega, options%eta, options%fermi_level, &
         options%electronic_temperature, band_first, band_last, options%occupation_prune_tolerance, batch_size, &
         options%use_batched_accumulation, provider, result%xi, result%metadata%vertex_cpu_seconds, &
         result%metadata%transition_preparation_cpu_seconds, result%metadata%denominator_cpu_seconds, &
         result%metadata%accumulation_cpu_seconds)
      if (.false.) then
      allocate(left_vertex(nleft), right_vertex(nright))

      if (options%use_batched_accumulation) then
         batch_size = min(options%transition_batch_size, (band_last-band_first+1)**2)
         result%metadata%transition_batch_size = batch_size
         allocate(left_batch(nleft, batch_size), right_batch(nright, batch_size), weighted_left(nleft, batch_size), &
            denominator_batch(batch_size), occupation_batch(batch_size), transition_energy_batch(batch_size))
         do ik = 1, nk
            prefactor = k_weights(ik)/weight_sum; npairs = 0
            do n = band_first, band_last
               do m = band_first, band_last
                  occupation_difference = tddft_fermi_occupation(eigenvalues_k(n, ik), options%fermi_level, &
                     options%electronic_temperature) - tddft_fermi_occupation(eigenvalues_kq(m, ik), &
                     options%fermi_level, options%electronic_temperature)
                  if (options%occupation_prune_tolerance > 0.0_rp) then
                     if (abs(occupation_difference) <= options%occupation_prune_tolerance) cycle
                  end if
                  npairs = npairs + 1
                  call cpu_time(t_start)
                  call build_transition_vertices(left_channels, weighted_right_operators, site_orbital_counts, &
                     eigenvectors_k(:, n, ik), eigenvectors_kq(:, m, ik), left_batch(:, npairs), right_batch(:, npairs))
                  call cpu_time(t_stop)
                  result%metadata%vertex_cpu_seconds = result%metadata%vertex_cpu_seconds + t_stop-t_start
                  occupation_batch(npairs) = occupation_difference
                  transition_energy_batch(npairs) = eigenvalues_k(n, ik) - eigenvalues_kq(m, ik)
                  if (npairs == batch_size) then
                     call accumulate_direct_xi_batch(result%xi, omega, options%eta, prefactor, left_batch, right_batch, &
                        occupation_batch, transition_energy_batch, npairs, weighted_left, denominator_batch, result%metadata)
                     npairs = 0
                  end if
               end do
            end do
            if (npairs > 0) then
               call accumulate_direct_xi_batch(result%xi, omega, options%eta, prefactor, left_batch, right_batch, &
                  occupation_batch, transition_energy_batch, npairs, weighted_left, denominator_batch, result%metadata)
            end if
         end do
         deallocate(left_batch, right_batch, weighted_left, denominator_batch, occupation_batch, transition_energy_batch)
      else
         do ik = 1, nk
            prefactor = k_weights(ik)/weight_sum
            do n = band_first, band_last
               do m = band_first, band_last
                  occupation_difference = tddft_fermi_occupation(eigenvalues_k(n, ik), options%fermi_level, &
                     options%electronic_temperature) - tddft_fermi_occupation(eigenvalues_kq(m, ik), &
                     options%fermi_level, options%electronic_temperature)
                  if (options%occupation_prune_tolerance > 0.0_rp) then
                     if (abs(occupation_difference) <= options%occupation_prune_tolerance) cycle
                  end if
                  call cpu_time(t_start)
                  call build_transition_vertices(left_channels, weighted_right_operators, site_orbital_counts, &
                     eigenvectors_k(:, n, ik), eigenvectors_kq(:, m, ik), left_vertex, right_vertex)
                  call cpu_time(t_stop)
                  result%metadata%vertex_cpu_seconds = result%metadata%vertex_cpu_seconds + t_stop-t_start
                  transition_energy = eigenvalues_k(n, ik) - eigenvalues_kq(m, ik)
                  do iw = 1, nw
                     call cpu_time(t_start)
                     denominator = cmplx(omega(iw) + transition_energy, options%eta, rp)
                     call cpu_time(t_stop)
                     result%metadata%denominator_cpu_seconds = result%metadata%denominator_cpu_seconds + t_stop-t_start
                     call cpu_time(t_start)
                     result%xi(:, :, iw) = result%xi(:, :, iw) + prefactor*occupation_difference* &
                        outer_product(left_vertex, right_vertex)/denominator
                     call cpu_time(t_stop)
                     result%metadata%accumulation_cpu_seconds = result%metadata%accumulation_cpu_seconds + t_stop-t_start
                  end do
               end do
            end do
         end do
      end if
      deallocate(left_vertex, right_vertex)
      end if
   end subroutine build_direct_xi_from_eigenpairs

   !> K-resolved variant of the direct construction.  The last dimension of
   !> `weighted_right_operators` is the same k index as the eigenpairs, so a
   !> finite-q pair-potential matrix is never replaced by a site scalar or a
   !> k-averaged operator before the Lehmann sum.
   subroutine build_direct_xi_from_k_dependent_eigenpairs(k_weights, eigenvalues_k, eigenvectors_k, eigenvalues_kq, &
      eigenvectors_kq, site_orbital_counts, left_channels, weighted_right_operators, omega, options, result)
      real(rp), target, intent(in) :: k_weights(:), eigenvalues_k(:, :), eigenvalues_kq(:, :), omega(:)
      complex(rp), target, intent(in) :: eigenvectors_k(:, :, :), eigenvectors_kq(:, :, :)
      integer, target, intent(in) :: site_orbital_counts(:)
      type(response_channel), target, intent(in) :: left_channels(:)
      complex(rp), target, intent(in) :: weighted_right_operators(:, :, :, :)
      type(tddft_chi0_options), intent(in) :: options
      type(tddft_direct_xi_result), intent(out) :: result

      integer :: nk, nbands, nspinor, nleft, nright, nw, ik, n, m, iw, npairs, batch_size
      integer :: band_first, band_last
      real(rp) :: weight_sum, occupation_difference, transition_energy, prefactor, t_start, t_stop
      complex(rp) :: denominator
      complex(rp), allocatable :: left_vertex(:), right_vertex(:)
      complex(rp), allocatable :: left_batch(:, :), right_batch(:, :), weighted_left(:, :), denominator_batch(:)
      real(rp), allocatable :: occupation_batch(:), transition_energy_batch(:)
      type(tddft_transition_engine) :: engine
      type(pair_operator_vertex_provider) :: provider

      nk = size(k_weights); nbands = size(eigenvalues_k, 1); nspinor = 2*sum(site_orbital_counts)
      nleft = size(left_channels); nright = size(weighted_right_operators, 3); nw = size(omega)
      call validate_direct_xi_inputs(nk, nbands, nspinor, nleft, nright, nw, k_weights, eigenvalues_k, eigenvectors_k, &
         eigenvalues_kq, eigenvectors_kq, weighted_right_operators(:, :, :, 1), options)
      if (size(weighted_right_operators, 4) /= nk) then
         error stop 'build_direct_xi_from_k_dependent_eigenpairs: operator k dimension is incompatible'
      end if
      band_first = options%band_first; band_last = options%band_last
      if (band_last == 0) band_last = nbands
      if (band_first < 1 .or. band_last < band_first .or. band_last > nbands) then
         error stop 'build_direct_xi_from_k_dependent_eigenpairs: invalid selected band window'
      end if
      weight_sum = sum(k_weights)
      allocate(result%xi(nleft, nright, nw)); result%xi = cmplx(0.0_rp, 0.0_rp, rp)
      call set_direct_xi_metadata(result%metadata, options, nk, nbands, band_first, band_last, weight_sum)
      result%metadata%backend = 'direct_xi_k_resolved_eigenpairs'
      batch_size = min(options%transition_batch_size, (band_last-band_first+1)**2)
      call make_pair_operator_vertex_provider(provider, site_orbital_counts, left_channels, eigenvectors_k, eigenvectors_kq, &
         weighted_right_operators, .true.)
      call engine%accumulate_dynamic(k_weights, eigenvalues_k, eigenvalues_kq, omega, options%eta, options%fermi_level, &
         options%electronic_temperature, band_first, band_last, options%occupation_prune_tolerance, batch_size, &
         options%use_batched_accumulation, provider, result%xi, result%metadata%vertex_cpu_seconds, &
         result%metadata%transition_preparation_cpu_seconds, result%metadata%denominator_cpu_seconds, &
         result%metadata%accumulation_cpu_seconds)
      if (.false.) then
      allocate(left_vertex(nleft), right_vertex(nright))

      if (options%use_batched_accumulation) then
         batch_size = min(options%transition_batch_size, (band_last-band_first+1)**2)
         result%metadata%transition_batch_size = batch_size
         allocate(left_batch(nleft, batch_size), right_batch(nright, batch_size), weighted_left(nleft, batch_size), &
            denominator_batch(batch_size), occupation_batch(batch_size), transition_energy_batch(batch_size))
         do ik = 1, nk
            prefactor = k_weights(ik)/weight_sum; npairs = 0
            do n = band_first, band_last
               do m = band_first, band_last
                  occupation_difference = tddft_fermi_occupation(eigenvalues_k(n, ik), options%fermi_level, &
                     options%electronic_temperature) - tddft_fermi_occupation(eigenvalues_kq(m, ik), &
                     options%fermi_level, options%electronic_temperature)
                  if (options%occupation_prune_tolerance > 0.0_rp) then
                     if (abs(occupation_difference) <= options%occupation_prune_tolerance) cycle
                  end if
                  npairs = npairs + 1
                  call cpu_time(t_start)
                  call build_transition_vertices(left_channels, weighted_right_operators(:, :, :, ik), site_orbital_counts, &
                     eigenvectors_k(:, n, ik), eigenvectors_kq(:, m, ik), left_batch(:, npairs), right_batch(:, npairs))
                  call cpu_time(t_stop)
                  result%metadata%vertex_cpu_seconds = result%metadata%vertex_cpu_seconds + t_stop-t_start
                  occupation_batch(npairs) = occupation_difference
                  transition_energy_batch(npairs) = eigenvalues_k(n, ik) - eigenvalues_kq(m, ik)
                  if (npairs == batch_size) then
                     call accumulate_direct_xi_batch(result%xi, omega, options%eta, prefactor, left_batch, right_batch, &
                        occupation_batch, transition_energy_batch, npairs, weighted_left, denominator_batch, result%metadata)
                     npairs = 0
                  end if
               end do
            end do
            if (npairs > 0) call accumulate_direct_xi_batch(result%xi, omega, options%eta, prefactor, left_batch, right_batch, &
               occupation_batch, transition_energy_batch, npairs, weighted_left, denominator_batch, result%metadata)
         end do
         deallocate(left_batch, right_batch, weighted_left, denominator_batch, occupation_batch, transition_energy_batch)
      else
         do ik = 1, nk
            prefactor = k_weights(ik)/weight_sum
            do n = band_first, band_last
               do m = band_first, band_last
                  occupation_difference = tddft_fermi_occupation(eigenvalues_k(n, ik), options%fermi_level, &
                     options%electronic_temperature) - tddft_fermi_occupation(eigenvalues_kq(m, ik), &
                     options%fermi_level, options%electronic_temperature)
                  if (options%occupation_prune_tolerance > 0.0_rp) then
                     if (abs(occupation_difference) <= options%occupation_prune_tolerance) cycle
                  end if
                  call cpu_time(t_start)
                  call build_transition_vertices(left_channels, weighted_right_operators(:, :, :, ik), site_orbital_counts, &
                     eigenvectors_k(:, n, ik), eigenvectors_kq(:, m, ik), left_vertex, right_vertex)
                  call cpu_time(t_stop)
                  result%metadata%vertex_cpu_seconds = result%metadata%vertex_cpu_seconds + t_stop-t_start
                  transition_energy = eigenvalues_k(n, ik) - eigenvalues_kq(m, ik)
                  do iw = 1, nw
                     denominator = cmplx(omega(iw) + transition_energy, options%eta, rp)
                     result%xi(:, :, iw) = result%xi(:, :, iw) + prefactor*occupation_difference* &
                        outer_product(left_vertex, right_vertex)/denominator
                  end do
               end do
            end do
         end do
      end if
      deallocate(left_vertex, right_vertex)
      end if
   end subroutine build_direct_xi_from_k_dependent_eigenpairs

   subroutine build_transition_vertices(left_channels, operators, site_orbital_counts, n_spinor, m_spinor, left, right)
      type(response_channel), intent(in) :: left_channels(:)
      complex(rp), intent(in) :: operators(:, :, :), n_spinor(:), m_spinor(:)
      integer, intent(in) :: site_orbital_counts(:)
      complex(rp), intent(out) :: left(:), right(:)
      integer :: ia

      do ia = 1, size(left_channels)
         left(ia) = response_transition_vertex(left_channels(ia), site_orbital_counts, n_spinor, m_spinor)
      end do
      do ia = 1, size(right)
         right(ia) = weighted_transition_vertex(operators(:, :, ia), m_spinor, n_spinor)
      end do
   end subroutine build_transition_vertices

   subroutine accumulate_direct_xi_batch(xi, omega, eta, prefactor, left_vertices, right_vertices, occupations, &
      transition_energies, npairs, weighted_left, denominators, metadata)
      complex(rp), intent(inout) :: xi(:, :, :)
      real(rp), intent(in) :: omega(:), eta, prefactor, occupations(:), transition_energies(:)
      complex(rp), intent(in) :: left_vertices(:, :), right_vertices(:, :)
      integer, intent(in) :: npairs
      complex(rp), intent(inout) :: weighted_left(:, :), denominators(:)
      type(tddft_chi0_metadata), intent(inout) :: metadata
      integer :: iw
      real(rp) :: t_start, t_stop

      do iw = 1, size(omega)
         call cpu_time(t_start)
         denominators(1:npairs) = cmplx(omega(iw) + transition_energies(1:npairs), eta, rp)
         call cpu_time(t_stop)
         metadata%denominator_cpu_seconds = metadata%denominator_cpu_seconds + t_stop-t_start
         weighted_left(:, 1:npairs) = left_vertices(:, 1:npairs)*spread( &
            cmplx(prefactor*occupations(1:npairs), 0.0_rp, rp)/denominators(1:npairs), dim=1, ncopies=size(left_vertices, 1))
         call cpu_time(t_start)
         xi(:, :, iw) = xi(:, :, iw) + matmul(weighted_left(:, 1:npairs), transpose(right_vertices(:, 1:npairs)))
         call cpu_time(t_stop)
         metadata%accumulation_cpu_seconds = metadata%accumulation_cpu_seconds + t_stop-t_start
      end do
   end subroutine accumulate_direct_xi_batch

   function outer_product(left, right) result(product)
      complex(rp), intent(in) :: left(:), right(:)
      complex(rp) :: product(size(left), size(right))
      integer :: i, j

      do j = 1, size(right)
         do i = 1, size(left)
            product(i, j) = left(i)*right(j)
         end do
      end do
   end function outer_product

   subroutine set_direct_xi_metadata(metadata, options, nk, nbands, band_first, band_last, weight_sum)
      type(tddft_chi0_metadata), intent(out) :: metadata
      type(tddft_chi0_options), intent(in) :: options
      integer, intent(in) :: nk, nbands, band_first, band_last
      real(rp), intent(in) :: weight_sum

      metadata%backend = 'direct_xi_eigenpairs'
      metadata%energy_integration = 'not applicable'
      metadata%eta = options%eta; metadata%fermi_level = options%fermi_level
      metadata%electronic_temperature = options%electronic_temperature
      metadata%electronic_kT = max(options%electronic_temperature*tddft_kB_Ry_per_K, tddft_occupation_kT_floor)
      metadata%k_weight_sum = weight_sum; metadata%k_mesh_shape = options%k_mesh_shape; metadata%nk = nk
      metadata%available_band_count = nbands; metadata%band_first = band_first; metadata%band_last = band_last
      metadata%occupation_prune_tolerance = options%occupation_prune_tolerance
      metadata%batched_accumulation = options%use_batched_accumulation
      if (options%use_batched_accumulation) metadata%transition_batch_size = options%transition_batch_size
   end subroutine set_direct_xi_metadata

   subroutine validate_direct_xi_inputs(nk, nbands, nspinor, nleft, nright, nw, k_weights, eigenvalues_k, &
      eigenvectors_k, eigenvalues_kq, eigenvectors_kq, operators, options)
      integer, intent(in) :: nk, nbands, nspinor, nleft, nright, nw
      real(rp), intent(in) :: k_weights(:), eigenvalues_k(:, :), eigenvalues_kq(:, :)
      complex(rp), intent(in) :: eigenvectors_k(:, :, :), eigenvectors_kq(:, :, :), operators(:, :, :)
      type(tddft_chi0_options), intent(in) :: options

      if (nk <= 0 .or. nbands <= 0 .or. nspinor <= 0 .or. nleft <= 0 .or. nright <= 0 .or. nw <= 0) then
         error stop 'build_direct_xi_from_eigenpairs: empty input is not valid'
      end if
      if (size(eigenvalues_k, 2) /= nk .or. any(shape(eigenvalues_kq) /= shape(eigenvalues_k)) .or. &
          size(eigenvectors_k, 1) /= nspinor .or. size(eigenvectors_k, 2) /= nbands .or. &
          size(eigenvectors_k, 3) /= nk .or. any(shape(eigenvectors_kq) /= shape(eigenvectors_k)) .or. &
          size(operators, 1) /= nspinor .or. size(operators, 2) /= nspinor) then
         error stop 'build_direct_xi_from_eigenpairs: eigenpair or operator shapes are incompatible'
      end if
      if (any(k_weights < 0.0_rp) .or. sum(k_weights) <= tiny(1.0_rp)) then
         error stop 'build_direct_xi_from_eigenpairs: k weights must have a positive sum'
      end if
      if (options%eta <= 0.0_rp .or. options%electronic_temperature < 0.0_rp .or. &
          options%occupation_prune_tolerance < 0.0_rp .or. options%transition_batch_size < 1) then
         error stop 'build_direct_xi_from_eigenpairs: invalid response options'
      end if
   end subroutine validate_direct_xi_inputs

   subroutine validate_direct_xi_source_inputs(nk, nbands, nspinor, nleft, nright, nw, k_weights, eigenvalues_k, &
      eigenvectors_k, eigenvalues_kq, eigenvectors_kq, options)
      integer, intent(in) :: nk, nbands, nspinor, nleft, nright, nw
      real(rp), intent(in) :: k_weights(:), eigenvalues_k(:, :), eigenvalues_kq(:, :)
      complex(rp), intent(in) :: eigenvectors_k(:, :, :), eigenvectors_kq(:, :, :)
      type(tddft_chi0_options), intent(in) :: options

      if (nk <= 0 .or. nbands <= 0 .or. nspinor <= 0 .or. nleft <= 0 .or. nright <= 0 .or. nw <= 0) then
         error stop 'build_direct_xi_from_operator_source: empty input is not valid'
      end if
      if (size(eigenvalues_k, 2) /= nk .or. any(shape(eigenvalues_kq) /= shape(eigenvalues_k)) .or. &
          size(eigenvectors_k, 1) /= nspinor .or. size(eigenvectors_k, 2) /= nbands .or. &
          size(eigenvectors_k, 3) /= nk .or. any(shape(eigenvectors_kq) /= shape(eigenvectors_k))) then
         error stop 'build_direct_xi_from_operator_source: eigenpair shapes are incompatible'
      end if
      if (any(k_weights < 0.0_rp) .or. sum(k_weights) <= tiny(1.0_rp)) then
         error stop 'build_direct_xi_from_operator_source: k weights must have a positive sum'
      end if
      if (options%eta <= 0.0_rp .or. options%electronic_temperature < 0.0_rp .or. &
          options%occupation_prune_tolerance < 0.0_rp .or. options%transition_batch_size < 1) then
         error stop 'build_direct_xi_from_operator_source: invalid response options'
      end if
   end subroutine validate_direct_xi_source_inputs

end module tddft_xi_mod
