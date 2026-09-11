!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Native real-space coefficient-GF construction of the LR-04 response.
!>
!> This module owns the real-space backend only.  A provider supplies directed
!> coefficient-space G_ij/G_ji blocks at arbitrary complex energies.  The
!> backend applies the certified RSGF-00 endpoint augmentation to both
!> directions, evaluates the same LR-GF-02 retarded/advanced Kubo bubble, and
!> Fourier assembles the complete radial/angular response matrix.
!>
!> The provider boundary is intentional: block recursion and Chebyshev expose
!> different convergence controls, but neither is allowed to change the
!> response-space algebra or the energy-integration contract.
!------------------------------------------------------------------------------
module lr_rs_gf_susceptibility_mod

   use precision_mod, only: rp
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, lmto_orbital_l
   use response_angular_basis_mod, only: response_angular_pi, response_gaunt
   use response_basis_mapping_mod, only: response_super_index, response_unflatten_superindex, &
      response_real_space_phase
   use lr_response_space_mod, only: response_space_layout, response_raw_to_canonical
   use lr_pauli_transition_vertex_mod, only: pauli_vertex_capabilities, pauli_sigma_plus_matrix, &
      pauli_sigma_minus_matrix
   use lr_ks_susceptibility_mod, only: lr_ks_susceptibility_result, lr_channel_plus, lr_channel_minus, &
      lr_fermi_dirac_occupation
   use lr_gf_endpoint_augmentation_mod, only: lr_gf_augmented_block, augment_lr_gf_endpoint_pair
   implicit none
   private

   real(rp), parameter :: mesh_tolerance = 2.0e-13_rp
   real(rp), parameter :: translation_tolerance = 2.0e-13_rp

   !> One directed real-space pair.  The coefficient-GF provider interprets
   !> gij as G_(left,right)(z), with the first index as the row/destination
   !> site and the second index as the column/source site.  gji is the
   !> separately evaluated reverse block G_(right,left)(z); it is not inferred
   !> by transposition.
   type, public :: lr_rs_gf_pair
      integer :: left_site = 0
      integer :: right_site = 0
      real(rp) :: translation(3) = 0.0_rp
   end type lr_rs_gf_pair

   !> Abstract native coefficient-GF provider contract.
   type, abstract, public :: lr_rs_gf_provider
      integer :: nbasis_per_site = 0
      character(len=32) :: provider_kind = ''
   contains
      procedure(lr_rs_get_pair), deferred :: get_pair
      procedure(lr_rs_describe), deferred :: describe
   end type lr_rs_gf_provider

   abstract interface
      subroutine lr_rs_get_pair(this, left_site, right_site, translation, z, gij, gji, hgamma_ij, hgamma_ji)
         import :: lr_rs_gf_provider, rp
         class(lr_rs_gf_provider), intent(inout) :: this
         integer, intent(in) :: left_site, right_site
         real(rp), intent(in) :: translation(3)
         complex(rp), intent(in) :: z
         complex(rp), intent(out) :: gij(:, :), gji(:, :), hgamma_ij(:, :), hgamma_ji(:, :)
      end subroutine lr_rs_get_pair

      function lr_rs_describe(this) result(description)
         import :: lr_rs_gf_provider
         class(lr_rs_gf_provider), intent(in) :: this
         character(len=256) :: description
      end function lr_rs_describe
   end interface

   !> Callback signature for a native block-recursion provider.  The callback
   !> is the seam to the production two-sweep/native block implementation.
   abstract interface
      subroutine lr_rs_native_callback(left_site, right_site, translation, z, gij, gji, hgamma_ij, hgamma_ji)
         import :: rp
         integer, intent(in) :: left_site, right_site
         real(rp), intent(in) :: translation(3)
         complex(rp), intent(in) :: z
         complex(rp), intent(out) :: gij(:, :), gji(:, :), hgamma_ij(:, :), hgamma_ji(:, :)
      end subroutine lr_rs_native_callback
   end interface

   !> Block-recursion native coefficient-GF provider.  It is kept distinct
   !> from the Chebyshev provider so recursion depth and terminator controls
   !> cannot be accidentally reported as polynomial controls.
   type, extends(lr_rs_gf_provider), public :: lr_rs_block_recursion_provider
      procedure(lr_rs_native_callback), pointer, nopass :: callback => null()
      integer :: recursion_depth = 0
      character(len=64) :: terminator = 'certified native block terminator'
   contains
      procedure :: initialize => lr_rs_block_recursion_initialize
      procedure :: get_pair => lr_rs_block_recursion_get_pair
      procedure :: describe => lr_rs_block_recursion_describe
   end type lr_rs_block_recursion_provider

   !> Chebyshev native coefficient-GF provider.  It has an independent
   !> polynomial-order control and the same directed-block contract.
   type, extends(lr_rs_gf_provider), public :: lr_rs_chebyshev_provider
      procedure(lr_rs_native_callback), pointer, nopass :: callback => null()
      integer :: polynomial_order = 0
      character(len=64) :: kernel = 'certified native Chebyshev kernel'
   contains
      procedure :: initialize => lr_rs_chebyshev_initialize
      procedure :: get_pair => lr_rs_chebyshev_get_pair
      procedure :: describe => lr_rs_chebyshev_describe
   end type lr_rs_chebyshev_provider

   !> Dense inverse provider used by the finite exact oracle.  It is also a
   !> useful representation-isolation provider: it never calls reciprocal GF
   !> code and returns both directed blocks from an explicit coefficient-space
   !> Hamiltonian inverse.
   type, extends(lr_rs_gf_provider), public :: lr_rs_dense_gf_provider
      complex(rp), allocatable :: hamiltonian(:, :)
      complex(rp), allocatable :: hgamma(:, :)
      integer :: nsite = 0
      integer :: block_size = 0
      logical :: translation_independent = .false.
   contains
      procedure :: initialize => lr_rs_dense_initialize
      procedure :: get_pair => lr_rs_dense_get_pair
      procedure :: describe => lr_rs_dense_describe
   end type lr_rs_dense_gf_provider

   !> Request for one real-space pair/Fourier/frequency sweep.
   type, public :: lr_rs_gf_susceptibility_request
      real(rp) :: q(3) = 0.0_rp
      real(rp), allocatable :: frequencies(:)
      real(rp) :: eta = 0.0_rp
      character(len=32) :: channel = lr_channel_plus
      integer :: integration_points = 2001
      real(rp) :: integration_eta = 0.0_rp
      real(rp) :: energy_margin = 1.0_rp
      real(rp) :: energy_min = -1.0_rp
      real(rp) :: energy_max = 1.0_rp
      real(rp) :: fermi_level = 0.0_rp
      real(rp) :: temperature = 0.0_rp
      type(response_space_layout), pointer :: response_space => null()
      type(lmto_radial_basis), pointer :: radial_bases(:) => null()
      class(lr_rs_gf_provider), pointer :: provider => null()
      type(lr_rs_gf_pair), allocatable :: pairs(:)
      ! Fractional/direct basis-site positions tau(:,site).  A null pointer
      ! means all endpoint positions are zero (the finite-cluster gauge).
      real(rp), pointer :: site_positions(:, :) => null()
   end type lr_rs_gf_susceptibility_request

   public :: evaluate_lr_rs_gf_susceptibility

contains

   subroutine lr_rs_block_recursion_initialize(this, callback, nbasis_per_site, recursion_depth, terminator)
      class(lr_rs_block_recursion_provider), intent(out) :: this
      procedure(lr_rs_native_callback) :: callback
      integer, intent(in) :: nbasis_per_site, recursion_depth
      character(len=*), intent(in), optional :: terminator

      if (nbasis_per_site < 1 .or. recursion_depth < 1) then
         error stop 'lr_rs_block_recursion_provider%initialize: invalid callback or recursion control'
      end if
      this%callback => callback
      this%nbasis_per_site = nbasis_per_site
      this%recursion_depth = recursion_depth
      this%provider_kind = 'block_recursion'
      this%terminator = 'certified native block terminator'
      if (present(terminator)) this%terminator = trim(terminator)
   end subroutine lr_rs_block_recursion_initialize

   subroutine lr_rs_block_recursion_get_pair(this, left_site, right_site, translation, z, gij, gji, hgamma_ij, hgamma_ji)
      class(lr_rs_block_recursion_provider), intent(inout) :: this
      integer, intent(in) :: left_site, right_site
      real(rp), intent(in) :: translation(3)
      complex(rp), intent(in) :: z
      complex(rp), intent(out) :: gij(:, :), gji(:, :), hgamma_ij(:, :), hgamma_ji(:, :)

      if (.not. associated(this%callback)) error stop 'block-recursion provider: callback is not initialized'
      call this%callback(left_site, right_site, translation, z, gij, gji, hgamma_ij, hgamma_ji)
   end subroutine lr_rs_block_recursion_get_pair

   function lr_rs_block_recursion_describe(this) result(description)
      class(lr_rs_block_recursion_provider), intent(in) :: this
      character(len=256) :: description

      write (description, '(a,i0,a,a)') 'provider=block_recursion recursion_depth=', this%recursion_depth, &
         ' terminator='//trim(this%terminator)
   end function lr_rs_block_recursion_describe

   subroutine lr_rs_chebyshev_initialize(this, callback, nbasis_per_site, polynomial_order, kernel)
      class(lr_rs_chebyshev_provider), intent(out) :: this
      procedure(lr_rs_native_callback) :: callback
      integer, intent(in) :: nbasis_per_site, polynomial_order
      character(len=*), intent(in), optional :: kernel

      if (nbasis_per_site < 1 .or. polynomial_order < 1) then
         error stop 'lr_rs_chebyshev_provider%initialize: invalid callback or polynomial control'
      end if
      this%callback => callback
      this%nbasis_per_site = nbasis_per_site
      this%polynomial_order = polynomial_order
      this%provider_kind = 'chebyshev'
      this%kernel = 'certified native Chebyshev kernel'
      if (present(kernel)) this%kernel = trim(kernel)
   end subroutine lr_rs_chebyshev_initialize

   subroutine lr_rs_chebyshev_get_pair(this, left_site, right_site, translation, z, gij, gji, hgamma_ij, hgamma_ji)
      class(lr_rs_chebyshev_provider), intent(inout) :: this
      integer, intent(in) :: left_site, right_site
      real(rp), intent(in) :: translation(3)
      complex(rp), intent(in) :: z
      complex(rp), intent(out) :: gij(:, :), gji(:, :), hgamma_ij(:, :), hgamma_ji(:, :)

      if (.not. associated(this%callback)) error stop 'Chebyshev provider: callback is not initialized'
      call this%callback(left_site, right_site, translation, z, gij, gji, hgamma_ij, hgamma_ji)
   end subroutine lr_rs_chebyshev_get_pair

   function lr_rs_chebyshev_describe(this) result(description)
      class(lr_rs_chebyshev_provider), intent(in) :: this
      character(len=256) :: description

      write (description, '(a,i0,a,a)') 'provider=chebyshev polynomial_order=', this%polynomial_order, &
         ' kernel='//trim(this%kernel)
   end function lr_rs_chebyshev_describe

   subroutine lr_rs_dense_initialize(this, hamiltonian, hgamma, block_size, translation_independent)
      class(lr_rs_dense_gf_provider), intent(out) :: this
      complex(rp), intent(in) :: hamiltonian(:, :), hgamma(:, :)
      integer, intent(in) :: block_size
      logical, intent(in), optional :: translation_independent

      if (block_size < 1 .or. size(hamiltonian, 1) /= size(hamiltonian, 2) .or. &
          any(shape(hgamma) /= shape(hamiltonian)) .or. mod(size(hamiltonian, 1), block_size) /= 0) then
         error stop 'lr_rs_dense_gf_provider%initialize: inconsistent Hamiltonian/block dimensions'
      end if
      this%hamiltonian = hamiltonian
      this%hgamma = hgamma
      this%block_size = block_size
      this%nsite = size(hamiltonian, 1)/block_size
      this%nbasis_per_site = block_size
      this%translation_independent = .false.
      if (present(translation_independent)) this%translation_independent = translation_independent
      this%provider_kind = 'dense_exact_coefficient_inverse'
   end subroutine lr_rs_dense_initialize

   subroutine lr_rs_dense_get_pair(this, left_site, right_site, translation, z, gij, gji, hgamma_ij, hgamma_ji)
      class(lr_rs_dense_gf_provider), intent(inout) :: this
      integer, intent(in) :: left_site, right_site
      real(rp), intent(in) :: translation(3)
      complex(rp), intent(in) :: z
      complex(rp), intent(out) :: gij(:, :), gji(:, :), hgamma_ij(:, :), hgamma_ji(:, :)
      complex(rp), allocatable :: matrix(:, :), inverse(:, :)
      integer :: left_begin, right_begin, left_end, right_end

      if (.not. allocated(this%hamiltonian)) error stop 'dense coefficient provider: not initialized'
      if (left_site < 1 .or. left_site > this%nsite .or. right_site < 1 .or. right_site > this%nsite) then
         error stop 'dense coefficient provider: site index is outside the finite system'
      end if
      if (maxval(abs(translation)) > translation_tolerance .and. .not. this%translation_independent) then
         error stop 'dense coefficient provider: nonzero translation is unsupported by the finite oracle'
      end if
      if (any(shape(gij) /= [this%block_size, this%block_size]) .or. &
          any(shape(gji) /= [this%block_size, this%block_size]) .or. &
          any(shape(hgamma_ij) /= [this%block_size, this%block_size]) .or. &
          any(shape(hgamma_ji) /= [this%block_size, this%block_size])) then
         error stop 'dense coefficient provider: output block shape mismatch'
      end if

      allocate(matrix(size(this%hamiltonian, 1), size(this%hamiltonian, 2)), &
         inverse(size(this%hamiltonian, 1), size(this%hamiltonian, 2)))
      matrix = -this%hamiltonian
      do left_begin = 1, size(matrix, 1)
         matrix(left_begin, left_begin) = z - this%hamiltonian(left_begin, left_begin)
      end do
      call invert_dense_matrix(matrix, inverse)

      left_begin = (left_site - 1)*this%block_size + 1
      left_end = left_site*this%block_size
      right_begin = (right_site - 1)*this%block_size + 1
      right_end = right_site*this%block_size
      gij = inverse(left_begin:left_end, right_begin:right_end)
      gji = inverse(right_begin:right_end, left_begin:left_end)
      hgamma_ij = this%hgamma(left_begin:left_end, right_begin:right_end)
      hgamma_ji = this%hgamma(right_begin:right_end, left_begin:left_end)
   end subroutine lr_rs_dense_get_pair

   function lr_rs_dense_describe(this) result(description)
      class(lr_rs_dense_gf_provider), intent(in) :: this
      character(len=256) :: description

      write (description, '(a,i0,a,i0,a,l1)') 'provider=dense_exact_coefficient_inverse nsite=', this%nsite, &
         ' block_size=', this%block_size, ' translation_independent=', this%translation_independent
   end function lr_rs_dense_describe

   subroutine evaluate_lr_rs_gf_susceptibility(request, result)
      type(lr_rs_gf_susceptibility_request), intent(in) :: request
      type(lr_ks_susceptibility_result), intent(out) :: result

      type(response_space_layout), pointer :: space
      type(lmto_radial_basis), pointer :: radial_bases(:)
      class(lr_rs_gf_provider), pointer :: provider
      type(pauli_vertex_capabilities) :: capabilities
      complex(rp), allocatable :: vertices(:, :, :, :), raw(:, :, :)
      complex(rp), allocatable :: gij(:, :), gji(:, :), hgamma_ij(:, :), hgamma_ji(:, :)
      complex(rp), allocatable :: gij_m(:, :), gji_m(:, :), hgamma_ij_m(:, :), hgamma_ji_m(:, :)
      complex(rp), allocatable :: gij_p(:, :), gji_p(:, :), hgamma_ij_p(:, :), hgamma_ji_p(:, :)
      complex(rp), allocatable :: gij_w(:, :), gji_w(:, :), hgamma_ij_w(:, :), hgamma_ji_w(:, :)
      complex(rp), allocatable :: gij_v(:, :), gji_v(:, :), hgamma_ij_v(:, :), hgamma_ji_v(:, :)
      type(lr_gf_augmented_block) :: ab_gr, ba_gr, ab_ga, ba_ga, ab_w, ba_w, ab_v, ba_v
      complex(rp) :: operator_matrix(2, 2), phase
      real(rp) :: integration_eta, energy_min, energy_max, step, energy, fermi_weight
      real(rp) :: quadrature_weight, scale
      real(rp), allocatable :: tau(:, :)
      complex(rp), allocatable :: pair_raw(:, :)
      integer :: channel_kind, nlocal, ne, ie, ipair, ifrequency

      call validate_rs_request(request, channel_kind)
      space => request%response_space
      radial_bases => request%radial_bases
      provider => request%provider
      capabilities = pauli_vertex_capabilities()
      if (channel_kind == 1) then
         operator_matrix = pauli_sigma_plus_matrix()
      else
         operator_matrix = pauli_sigma_minus_matrix()
      end if

      integration_eta = request%integration_eta
      if (integration_eta <= 0.0_rp) integration_eta = request%eta/40.0_rp
      if (integration_eta >= request%eta) then
         error stop 'evaluate_lr_rs_gf_susceptibility: integration_eta must be smaller than response eta'
      end if
      energy_min = request%energy_min - request%energy_margin
      energy_max = request%energy_max + request%energy_margin
      step = (energy_max - energy_min)/real(request%integration_points - 1, rp)
      ne = request%integration_points
      nlocal = provider%nbasis_per_site
      call build_native_vertex_tensor(space, radial_bases, operator_matrix, capabilities, vertices)
      allocate(raw(space%ndim, space%ndim, size(request%frequencies)))
      raw = cmplx(0.0_rp, 0.0_rp, rp)

      allocate(tau(3, space%nsite), pair_raw(space%ndim, space%ndim))
      tau = 0.0_rp
      if (associated(request%site_positions)) tau = request%site_positions

      allocate(gij(nlocal,nlocal), gji(nlocal,nlocal), hgamma_ij(nlocal,nlocal), hgamma_ji(nlocal,nlocal), &
         gij_m(nlocal,nlocal), gji_m(nlocal,nlocal), hgamma_ij_m(nlocal,nlocal), hgamma_ji_m(nlocal,nlocal), &
         gij_p(nlocal,nlocal), gji_p(nlocal,nlocal), hgamma_ij_p(nlocal,nlocal), hgamma_ji_p(nlocal,nlocal), &
         gij_w(nlocal,nlocal), gji_w(nlocal,nlocal), hgamma_ij_w(nlocal,nlocal), hgamma_ji_w(nlocal,nlocal), &
         gij_v(nlocal,nlocal), gji_v(nlocal,nlocal), hgamma_ij_v(nlocal,nlocal), hgamma_ji_v(nlocal,nlocal))

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
         fermi_weight = lr_fermi_dirac_occupation(energy, request%fermi_level, request%temperature)
         do ipair = 1, size(request%pairs)
            call provider%get_pair(request%pairs(ipair)%left_site, request%pairs(ipair)%right_site, &
               request%pairs(ipair)%translation, cmplx(energy, integration_eta, rp), gij_p, gji_p, hgamma_ij_p, hgamma_ji_p)
            call provider%get_pair(request%pairs(ipair)%left_site, request%pairs(ipair)%right_site, &
               request%pairs(ipair)%translation, cmplx(energy, -integration_eta, rp), gij_m, gji_m, hgamma_ij_m, hgamma_ji_m)
            call augment_lr_gf_endpoint_pair(radial_bases(request%pairs(ipair)%left_site), &
               radial_bases(request%pairs(ipair)%right_site), request%pairs(ipair)%left_site, &
               request%pairs(ipair)%right_site, cmplx(energy, integration_eta, rp), gij_p, gji_p, ab_gr, ba_gr, &
               hgamma_block=hgamma_ij_p, reverse_hgamma_block=hgamma_ji_p)
            call augment_lr_gf_endpoint_pair(radial_bases(request%pairs(ipair)%left_site), &
               radial_bases(request%pairs(ipair)%right_site), request%pairs(ipair)%left_site, &
               request%pairs(ipair)%right_site, cmplx(energy, -integration_eta, rp), gij_m, gji_m, ab_ga, ba_ga, &
               hgamma_block=hgamma_ij_m, reverse_hgamma_block=hgamma_ji_m)
            do ifrequency = 1, size(request%frequencies)
               call provider%get_pair(request%pairs(ipair)%left_site, request%pairs(ipair)%right_site, &
                  request%pairs(ipair)%translation, cmplx(energy + request%frequencies(ifrequency), request%eta, rp), &
                  gij_w, gji_w, hgamma_ij_w, hgamma_ji_w)
               call augment_lr_gf_endpoint_pair(radial_bases(request%pairs(ipair)%left_site), &
                  radial_bases(request%pairs(ipair)%right_site), request%pairs(ipair)%left_site, &
                  request%pairs(ipair)%right_site, cmplx(energy + request%frequencies(ifrequency), request%eta, rp), &
                  gij_w, gji_w, ab_w, ba_w, hgamma_block=hgamma_ij_w, reverse_hgamma_block=hgamma_ji_w)
               call provider%get_pair(request%pairs(ipair)%left_site, request%pairs(ipair)%right_site, &
                  request%pairs(ipair)%translation, cmplx(energy - request%frequencies(ifrequency), -request%eta, rp), &
                  gij_v, gji_v, hgamma_ij_v, hgamma_ji_v)
               call augment_lr_gf_endpoint_pair(radial_bases(request%pairs(ipair)%left_site), &
                  radial_bases(request%pairs(ipair)%right_site), request%pairs(ipair)%left_site, &
                  request%pairs(ipair)%right_site, cmplx(energy - request%frequencies(ifrequency), -request%eta, rp), &
                  gij_v, gji_v, ab_v, ba_v, hgamma_block=hgamma_ij_v, reverse_hgamma_block=hgamma_ji_v)

               scale = fermi_weight*quadrature_weight*2.0_rp
               pair_raw = cmplx(0.0_rp, 0.0_rp, rp)
               call accumulate_native_pair_bubble(space, vertices, ab_gr, ba_gr, ab_ga, ba_ga, ab_w, ba_v, &
                  request%pairs(ipair)%left_site, request%pairs(ipair)%right_site, scale, pair_raw)
               phase = response_real_space_phase(request%q, request%pairs(ipair)%translation, &
                  tau(:, request%pairs(ipair)%left_site), tau(:, request%pairs(ipair)%right_site))
               raw(:, :, ifrequency) = raw(:, :, ifrequency) + phase*pair_raw
            end do
         end do
      end do

      allocate(result%susceptibility(space%ndim, space%ndim, size(request%frequencies)), &
         result%frequencies(size(request%frequencies)))
      do ifrequency = 1, size(request%frequencies)
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
      write (result%response_space_metadata, '(a,a,a,i0,a,es12.4,a,es12.4,a,es12.4,a,es12.4)') &
         trim(provider%describe()), ' phase=+i2pi*q.(R+tau_right-tau_left)', &
         ' integration_points=', ne, ' integration_eta=', integration_eta, ' eta=', request%eta, &
         ' fermi_level=', request%fermi_level, ' temperature=', request%temperature
   end subroutine evaluate_lr_rs_gf_susceptibility

   subroutine accumulate_native_pair_bubble(space, vertices, ab_gr, ba_gr, ab_ga, ba_ga, ab_w, ba_v, &
                                            left_site, right_site, scale, raw)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: vertices(:, :, :, :)
      type(lr_gf_augmented_block), intent(in) :: ab_gr, ba_gr, ab_ga, ba_ga, ab_w, ba_v
      integer, intent(in) :: left_site, right_site
      real(rp), intent(in) :: scale
      complex(rp), intent(inout) :: raw(:, :)
      type(response_super_index) :: item_i, item_j
      integer :: i, j, p, q, r, s, component, ioffset, joffset, nlocal
      complex(rp), allocatable :: a_ba(:, :), a_ab(:, :), g_ab(:, :), g_ba(:, :), temporary(:, :)
      complex(rp), allocatable :: vi(:, :), vj(:, :), vi_sp(:, :), vj_rq(:, :)

      nlocal = ab_gr%norb*ab_gr%nspin
      allocate(a_ba(nlocal,nlocal), a_ab(nlocal,nlocal), g_ab(nlocal,nlocal), g_ba(nlocal,nlocal), &
         temporary(nlocal,nlocal), vi(nlocal,nlocal), vj(nlocal,nlocal), vi_sp(nlocal,nlocal), vj_rq(nlocal,nlocal))
      ioffset = (left_site - 1)*nlocal
      joffset = (right_site - 1)*nlocal
      do i = 1, space%ndim
         call response_unflatten_superindex(i, space%nsite, space%response_lmax, space%npoint, space%nchannel, item_i)
         if (item_i%site /= left_site) cycle
         do j = 1, space%ndim
            call response_unflatten_superindex(j, space%nsite, space%response_lmax, space%npoint, space%nchannel, item_j)
            if (item_j%site /= right_site) cycle
            do p = 0, 1
               do q = 0, 1
                  call native_branch(ba_gr, p, q, a_ba)
                  call native_branch(ba_ga, p, q, g_ba)
                  a_ba = cmplx(0.0_rp, 1.0_rp/(2.0_rp*response_angular_pi), rp)*(a_ba - g_ba)
                  call native_branch(ab_gr, p, q, g_ab)
                  call native_branch(ab_ga, p, q, a_ab)
                  a_ab = cmplx(0.0_rp, 1.0_rp/(2.0_rp*response_angular_pi), rp)*(g_ab - a_ab)
                  do r = 0, 1
                     do s = 0, 1
                        component = 1 + q + 2*r
                        vi = vertices(ioffset+1:ioffset+nlocal, ioffset+1:ioffset+nlocal, component, i)
                        component = 1 + p + 2*s
                        vj = vertices(joffset+1:joffset+nlocal, joffset+1:joffset+nlocal, component, j)
                        component = 1 + r + 2*s
                        call native_branch(ab_w, r, s, g_ab)
                        temporary = matmul(a_ba, vi)
                        temporary = matmul(temporary, g_ab)
                        raw(i, j) = raw(i, j) + scale*sum(temporary*conjg(vj))

                        component = 1 + q + 2*r
                        vj_rq = vertices(joffset+1:joffset+nlocal, joffset+1:joffset+nlocal, component, j)
                        component = 1 + s + 2*p
                        vi_sp = vertices(ioffset+1:ioffset+nlocal, ioffset+1:ioffset+nlocal, component, i)
                        call native_branch(ba_v, r, s, g_ba)
                        temporary = matmul(a_ab, conjg(transpose(vj_rq)))
                        temporary = matmul(temporary, g_ba)
                        raw(i, j) = raw(i, j) + scale*sum(temporary*transpose(vi_sp))
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine accumulate_native_pair_bubble

   subroutine native_branch(block, p, q, branch)
      type(lr_gf_augmented_block), intent(in) :: block
      integer, intent(in) :: p, q
      complex(rp), intent(out) :: branch(:, :)

      select case (1 + p + 2*q)
      case (1)
         branch = block%coefficient_gf
      case (2)
         branch = block%h_g
      case (3)
         branch = block%g_h
      case (4)
         branch = block%h_g_h
      case default
         error stop 'native_branch: invalid endpoint branch'
      end select
   end subroutine native_branch

   subroutine build_native_vertex_tensor(space, radial_bases, operator_matrix, capabilities, vertices)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      complex(rp), intent(in) :: operator_matrix(2, 2)
      type(pauli_vertex_capabilities), intent(in) :: capabilities
      complex(rp), allocatable, intent(out) :: vertices(:, :, :, :)
      type(response_super_index) :: item
      integer :: nbasis, norb, flat, p, q, component, isite, iorb, jorb, ispin, jspin
      integer :: orbital_l, orbital_lp, offset
      real(rp) :: radial_product

      norb = (radial_bases(1)%lmax + 1)**2
      nbasis = 2*norb*space%nsite
      allocate(vertices(nbasis, nbasis, 4, space%ndim))
      vertices = cmplx(0.0_rp, 0.0_rp, rp)
      do flat = 1, space%ndim
         call response_unflatten_superindex(flat, space%nsite, space%response_lmax, space%npoint, &
            space%nchannel, item)
         isite = item%site
         offset = (isite - 1)*2*norb
         do p = 0, 1
            do q = 0, 1
               component = 1 + p + 2*q
               do iorb = 1, norb
                  orbital_l = lmto_orbital_l(iorb)
                  do jorb = 1, norb
                     orbital_lp = lmto_orbital_l(jorb)
                     do ispin = 1, 2
                        do jspin = 1, 2
                           radial_product = native_radial_vertex_component(radial_bases(isite), item%radial_point, &
                              orbital_l, orbital_lp, ispin, jspin, p, q)
                           vertices(offset+(ispin-1)*norb+iorb, offset+(jspin-1)*norb+jorb, component, flat) = &
                              operator_matrix(ispin,jspin)*radial_product*response_gaunt(orbital_l, orbital_m(iorb, orbital_l), &
                              orbital_lp, orbital_m(jorb, orbital_lp), item%response_l, item%response_m)
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      call radial_bases(1)%require_supported(capabilities%reciprocal_mode, capabilities%hamiltonian_order, &
         capabilities%orthogonal, capabilities%collinear, capabilities%has_soc, capabilities%has_extra_operator)
   end subroutine build_native_vertex_tensor

   function native_radial_vertex_component(basis, radial_point, l, lp, ispin, jspin, p, q) result(value)
      type(lmto_radial_basis), intent(in) :: basis
      integer, intent(in) :: radial_point, l, lp, ispin, jspin, p, q
      real(rp) :: value
      real(rp) :: first, second, rfirst, rsecond

      if (radial_point /= 1) then
         if (basis%rofi(radial_point) <= tiny(1.0_rp)) error stop 'native radial vertex: invalid positive radial point'
         value = native_radial_component(basis, radial_point, l, ispin, p)* &
            native_radial_component(basis, radial_point, lp, jspin, q)/(basis%rofi(radial_point)**2)
         return
      end if
      if (l /= 0 .or. lp /= 0) then
         value = 0.0_rp
         return
      end if
      rfirst = basis%rofi(2)
      rsecond = basis%rofi(3)
      first = native_radial_component(basis, 2, l, ispin, p)*native_radial_component(basis, 2, lp, jspin, q)/(rfirst**2)
      second = native_radial_component(basis, 3, l, ispin, p)*native_radial_component(basis, 3, lp, jspin, q)/(rsecond**2)
      value = (first*rsecond**2 - second*rfirst**2)/(rsecond**2 - rfirst**2)
   end function native_radial_vertex_component

   pure real(rp) function native_radial_component(basis, ir, l, ispin, branch) result(value)
      type(lmto_radial_basis), intent(in) :: basis
      integer, intent(in) :: ir, l, ispin, branch

      if (branch == 0) then
         value = basis%phi_large(ir, l + 1, ispin)
      else
         value = basis%phidot_large(ir, l + 1, ispin)
      end if
   end function native_radial_component

   subroutine validate_rs_request(request, channel_kind)
      type(lr_rs_gf_susceptibility_request), intent(in) :: request
      integer, intent(out) :: channel_kind
      type(response_space_layout), pointer :: space
      type(lmto_radial_basis), pointer :: radial_bases(:)
      type(pauli_vertex_capabilities) :: capabilities
      integer :: isite, ipair, expected_nlocal

      if (.not. associated(request%response_space) .or. .not. associated(request%radial_bases) .or. &
          .not. associated(request%provider)) error stop 'evaluate_lr_rs_gf_susceptibility: request references are incomplete'
      if (.not. allocated(request%frequencies) .or. size(request%frequencies) < 1) then
         error stop 'evaluate_lr_rs_gf_susceptibility: at least one frequency is required'
      end if
      if (request%eta <= 0.0_rp .or. request%integration_points < 3 .or. &
          mod(request%integration_points, 2) == 0) then
         error stop 'evaluate_lr_rs_gf_susceptibility: invalid eta or Simpson integration_points'
      end if
      if (request%energy_margin < 0.0_rp .or. request%energy_max <= request%energy_min) then
         error stop 'evaluate_lr_rs_gf_susceptibility: invalid energy integration window'
      end if
      select case (trim(request%channel))
      case ('plus', 'chi_plus')
         channel_kind = 1
      case ('minus', 'chi_minus')
         channel_kind = 2
      case default
         error stop 'evaluate_lr_rs_gf_susceptibility: channel must be chi_plus or chi_minus'
      end select
      space => request%response_space
      radial_bases => request%radial_bases
      if (space%nchannel /= 1 .or. size(radial_bases) /= space%nsite) then
         error stop 'evaluate_lr_rs_gf_susceptibility: response channel/site mismatch'
      end if
      expected_nlocal = 2*(radial_bases(1)%lmax + 1)**2
      if (request%provider%nbasis_per_site /= expected_nlocal) then
         error stop 'evaluate_lr_rs_gf_susceptibility: provider coefficient block size differs from radial basis'
      end if
      capabilities = pauli_vertex_capabilities()
      do isite = 1, space%nsite
         call radial_bases(isite)%require_supported(capabilities%reciprocal_mode, capabilities%hamiltonian_order, &
            capabilities%orthogonal, capabilities%collinear, capabilities%has_soc, capabilities%has_extra_operator)
         if (radial_bases(isite)%lmax /= radial_bases(1)%lmax .or. radial_bases(isite)%nspin /= 2 .or. &
             radial_bases(isite)%npoint /= space%npoint .or. .not. allocated(radial_bases(isite)%rofi) .or. &
             .not. allocated(radial_bases(isite)%channel_present) .or. .not. all(radial_bases(isite)%channel_present)) then
            error stop 'evaluate_lr_rs_gf_susceptibility: incomplete radial basis provenance'
         end if
         if (abs(radial_bases(isite)%mesh_a - space%a) > mesh_tolerance*max(1.0_rp,abs(space%a)) .or. &
             abs(radial_bases(isite)%mesh_b - space%b) > mesh_tolerance*max(1.0_rp,abs(space%b)) .or. &
             maxval(abs(radial_bases(isite)%rofi - space%radius)) > mesh_tolerance*max(1.0_rp,maxval(abs(space%radius)))) then
            error stop 'evaluate_lr_rs_gf_susceptibility: radial mesh provenance differs from response space'
         end if
      end do
      if (.not. allocated(request%pairs) .or. size(request%pairs) < 1) then
         error stop 'evaluate_lr_rs_gf_susceptibility: at least one real-space pair is required'
      end if
      do ipair = 1, size(request%pairs)
         if (request%pairs(ipair)%left_site < 1 .or. request%pairs(ipair)%left_site > space%nsite .or. &
             request%pairs(ipair)%right_site < 1 .or. request%pairs(ipair)%right_site > space%nsite) then
            error stop 'evaluate_lr_rs_gf_susceptibility: pair site is outside the response space'
         end if
      end do
      if (associated(request%site_positions)) then
         if (any(shape(request%site_positions) /= [3, space%nsite])) then
            error stop 'evaluate_lr_rs_gf_susceptibility: site-position shape mismatch'
         end if
      end if
   end subroutine validate_rs_request

   pure integer function orbital_m(iorb, l) result(m)
      integer, intent(in) :: iorb, l
      m = iorb - l*l - l - 1
   end function orbital_m

   subroutine invert_dense_matrix(matrix, inverse)
      complex(rp), intent(in) :: matrix(:, :)
      complex(rp), intent(out) :: inverse(:, :)
      complex(rp), allocatable :: work(:, :), row(:)
      complex(rp) :: pivot, factor
      integer :: n, i, j, pivot_row

      if (size(matrix, 1) /= size(matrix, 2) .or. any(shape(inverse) /= shape(matrix))) then
         error stop 'invert_dense_matrix: square shape required'
      end if
      n = size(matrix, 1)
      allocate(work(n,n), row(n))
      work = matrix
      inverse = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, n
         inverse(i,i) = cmplx(1.0_rp, 0.0_rp, rp)
      end do
      do i = 1, n
         pivot_row = i
         do j = i + 1, n
            if (abs(work(j,i)) > abs(work(pivot_row,i))) pivot_row = j
         end do
         if (abs(work(pivot_row,i)) <= 100.0_rp*epsilon(1.0_rp)) then
            error stop 'invert_dense_matrix: singular finite fixture Hamiltonian resolvent'
         end if
         if (pivot_row /= i) then
            row = work(i,:); work(i,:) = work(pivot_row,:); work(pivot_row,:) = row
            row = inverse(i,:); inverse(i,:) = inverse(pivot_row,:); inverse(pivot_row,:) = row
         end if
         pivot = work(i,i)
         work(i,:) = work(i,:)/pivot
         inverse(i,:) = inverse(i,:)/pivot
         do j = 1, n
            if (j == i) cycle
            factor = work(j,i)
            work(j,:) = work(j,:) - factor*work(i,:)
            inverse(j,:) = inverse(j,:) - factor*inverse(i,:)
         end do
      end do
   end subroutine invert_dense_matrix

end module lr_rs_gf_susceptibility_mod
