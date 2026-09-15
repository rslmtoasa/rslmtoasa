!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Production adapter from the existing real-space recursion/Green
!>        services to the native RSGF provider contract.
!>
!> This adapter owns orchestration and coefficient-space plumbing only.  The
!> recursion, complex-energy Green reconstruction, and H_eff action remain the
!> existing production implementations; LR-04 augmentation and the native
!> response bubble remain in lr_rs_gf_susceptibility_mod.
!------------------------------------------------------------------------------
module tddft_native_rsgf_provider_mod

   use precision_mod, only: rp
   use basis_mod, only: nb, spin_off
   use math_mod, only: i_unit
   use mpi_mod, only: rank, numprocs, get_mpi_variables
   use string_mod, only: lower
   use control_mod, only: control
   use lattice_mod, only: lattice
   use hamiltonian_mod, only: hamiltonian
   use reciprocal_mod, only: reciprocal
   use recursion_mod, only: recursion
   use green_mod, only: green
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use lr_rs_gf_susceptibility_mod, only: lr_rs_gf_provider, lr_rs_gf_pair
   implicit none
   private

   type, extends(lr_rs_gf_provider), public :: tddft_native_rsgf_provider
      type(green), pointer :: green_obj => null()
      type(recursion), pointer :: recursion_obj => null()
      type(hamiltonian), pointer :: hamiltonian_obj => null()
      type(lattice), pointer :: lattice_obj => null()
      type(reciprocal), pointer :: reciprocal_obj => null()
      type(lmto_radial_basis), pointer :: radial_bases(:) => null()
      character(len=32) :: selection = ''
      integer :: recursion_depth = 0
      integer :: polynomial_order = 0
      type(lr_rs_gf_pair), allocatable :: pairs(:)
      integer, allocatable :: left_atoms(:), right_atoms(:)
      complex(rp), allocatable :: hgamma_ij(:, :, :), hgamma_ji(:, :, :)
      real(rp), allocatable :: a_inf(:, :, :, :), b_inf(:, :, :, :)
      complex(rp), allocatable :: psi(:, :, :), psi_out(:, :, :), identity(:, :)
   contains
      procedure :: initialize => tddft_native_rsgf_initialize
      procedure :: get_pair => tddft_native_rsgf_get_pair
      procedure :: describe => tddft_native_rsgf_describe
   end type tddft_native_rsgf_provider

contains

   subroutine tddft_native_rsgf_initialize(this, provider_name, green_obj, recursion_obj, hamiltonian_obj, &
                                            lattice_obj, reciprocal_obj, radial_bases)
      class(tddft_native_rsgf_provider), intent(out) :: this
      character(len=*), intent(in) :: provider_name
      type(green), target, intent(inout) :: green_obj
      type(recursion), target, intent(inout) :: recursion_obj
      type(hamiltonian), target, intent(in) :: hamiltonian_obj
      type(lattice), target, intent(in) :: lattice_obj
      type(reciprocal), target, intent(in) :: reciprocal_obj
      type(lmto_radial_basis), target, intent(in) :: radial_bases(:)

      integer :: ipair, nsite, left_atom, right_atom, nw, recursion_slot
      real(rp) :: a_inf0(4), b_inf0(4)
      character(len=32) :: selected
      logical :: onsite_only

      if (numprocs /= 1) then
         error stop 'native RSGF production provider: the registered baseline requires a serial replicated pair workset'
      end if
      if (lattice_obj%njij < 1 .or. .not. allocated(lattice_obj%ijpair)) then
         error stop 'native RSGF production provider: a complete lattice%ijpair set is required; site-only fallback is forbidden'
      end if
      if (size(radial_bases) < 1) error stop 'native RSGF production provider: no response-site radial bases were supplied'

      selected = trim(lower(provider_name))
      if (selected == 'auto') selected = trim(lower(recursion_obj%control%recur))
      select case (selected)
      case ('block', 'block_recursion')
         if (trim(lower(recursion_obj%control%recur)) /= 'block') then
            error stop 'native RSGF production provider: block provider requires control%recur=block'
         end if
         this%selection = 'block'
         this%provider_kind = 'block_recursion'
         this%recursion_depth = recursion_obj%control%lld
         if (this%recursion_depth < 1) error stop 'native RSGF production provider: invalid block recursion depth'
      case ('chebyshev')
         if (trim(lower(recursion_obj%control%recur)) /= 'chebyshev') then
            error stop 'native RSGF production provider: Chebyshev provider requires control%recur=chebyshev'
         end if
         this%selection = 'chebyshev'
         this%provider_kind = 'chebyshev'
         this%polynomial_order = recursion_obj%control%lld
         if (this%polynomial_order < 1) error stop 'native RSGF production provider: invalid Chebyshev polynomial order'
      case default
         error stop 'native RSGF production provider: use block or chebyshev'
      end select

      this%green_obj => green_obj
      this%recursion_obj => recursion_obj
      this%hamiltonian_obj => hamiltonian_obj
      this%lattice_obj => lattice_obj
      this%reciprocal_obj => reciprocal_obj
      this%radial_bases => radial_bases
      nsite = size(radial_bases)
      this%nbasis_per_site = 2*(radial_bases(1)%lmax + 1)**2
      if (size(lattice_obj%ijpair, 1) < lattice_obj%njij .or. size(lattice_obj%ijpair, 2) /= 2) then
         error stop 'native RSGF production provider: lattice%ijpair shape is incomplete'
      end if

      allocate(this%pairs(lattice_obj%njij), this%left_atoms(lattice_obj%njij), this%right_atoms(lattice_obj%njij), &
         this%hgamma_ij(this%nbasis_per_site, this%nbasis_per_site, lattice_obj%njij), &
         this%hgamma_ji(this%nbasis_per_site, this%nbasis_per_site, lattice_obj%njij))
      allocate(this%identity(nb, nb))
      this%identity = cmplx(0.0_rp, 0.0_rp, rp)
      do ipair = 1, nb
         this%identity(ipair, ipair) = cmplx(1.0_rp, 0.0_rp, rp)
      end do

      ! Pair recursion uses the same global pair partition as the existing
      ! transport/real-space service.  The serial gate above makes the local
      ! pair index identical to the request pair index.
      call get_mpi_variables(rank, lattice_obj%njij)
      onsite_only = all(lattice_obj%ijpair(1:lattice_obj%njij, 1) == &
                        lattice_obj%ijpair(1:lattice_obj%njij, 2))
      select case (this%selection)
      case ('block')
         if (.not. onsite_only) then
            call recursion_obj%recur_b_ij()
            call recursion_obj%zsqr()
         else if (.not. recursion_obj%b2_b_is_sqrt) then
            ! The accepted SCF state already ran recur_b()/zsqr() for its
            ! response atoms.  Reuse those coefficients for an all-on-site
            ! pair workset instead of repeating the same block recursion.
            call recursion_obj%zsqr()
         end if
         allocate(this%a_inf(nb, nb, 4, lattice_obj%njij), this%b_inf(nb, nb, 4, lattice_obj%njij))
      case ('chebyshev')
         call recursion_obj%chebyshev_recur_ij()
      end select

      do ipair = 1, lattice_obj%njij
         left_atom = lattice_obj%ijpair(ipair, 1)
         right_atom = lattice_obj%ijpair(ipair, 2)
         call validate_atom_index(lattice_obj, left_atom, 'left')
         call validate_atom_index(lattice_obj, right_atom, 'right')
         this%left_atoms(ipair) = left_atom
         this%right_atoms(ipair) = right_atom
         this%pairs(ipair)%left_site = response_site_for_atom(lattice_obj, left_atom, nsite)
         this%pairs(ipair)%right_site = response_site_for_atom(lattice_obj, right_atom, nsite)
         call pair_translation(reciprocal_obj, lattice_obj, left_atom, right_atom, this%pairs(ipair)%translation)
         call native_h_eff_block(this, right_atom, left_atom, this%hgamma_ij(:, :, ipair))
         call native_h_eff_block(this, left_atom, right_atom, this%hgamma_ji(:, :, ipair))
         if (left_atom == right_atom) then
            call remove_radial_enu(this%hgamma_ij(:, :, ipair), this%radial_bases(this%pairs(ipair)%left_site))
            call remove_radial_enu(this%hgamma_ji(:, :, ipair), this%radial_bases(this%pairs(ipair)%right_site))
         end if
         if (this%selection == 'block') then
            nw = 10*this%recursion_depth
            if (onsite_only) then
               recursion_slot = recursion_site_slot(lattice_obj, left_atom)
               if (recursion_slot < 1) then
                  error stop 'native RSGF production provider: on-site pair is not an accepted recursion atom'
               end if
               call recursion_obj%get_terminf(recursion_obj%a_b(:, :, :, recursion_slot), &
                  recursion_obj%b2_b(:, :, :, recursion_slot), 1, this%recursion_depth, nb, nw, &
                  this%a_inf(:, :, 1:1, ipair), this%b_inf(:, :, 1:1, ipair), a_inf0(1:1), b_inf0(1:1))
               this%a_inf(:, :, 2:4, ipair) = 0.0_rp
               this%b_inf(:, :, 2:4, ipair) = 0.0_rp
            else
               call recursion_obj%get_terminf(recursion_obj%a_b(:, :, :, (ipair - 1)*4 + 1), &
                  recursion_obj%b2_b(:, :, :, (ipair - 1)*4 + 1), 4, this%recursion_depth, nb, nw, &
                  this%a_inf(:, :, :, ipair), this%b_inf(:, :, :, ipair), a_inf0, b_inf0)
            end if
         end if
      end do
   end subroutine tddft_native_rsgf_initialize

   subroutine tddft_native_rsgf_get_pair(this, left_site, right_site, translation, z, gij, gji, hgamma_ij, hgamma_ji)
      class(tddft_native_rsgf_provider), intent(inout) :: this
      integer, intent(in) :: left_site, right_site
      real(rp), intent(in) :: translation(3)
      complex(rp), intent(in) :: z
      complex(rp), intent(out) :: gij(:, :), gji(:, :), hgamma_ij(:, :), hgamma_ji(:, :)

      complex(rp) :: gphase(nb, nb, 4), g_one(nb, nb, 1), z_grid(1)
      integer :: ipair, phase

      if (any(shape(gij) /= [nb, nb]) .or. any(shape(gji) /= [nb, nb]) .or. &
          any(shape(hgamma_ij) /= [nb, nb]) .or. any(shape(hgamma_ji) /= [nb, nb])) then
         error stop 'native RSGF production provider: coefficient block shape mismatch'
      end if
      ipair = find_pair(this, left_site, right_site, translation)
      hgamma_ij = this%hgamma_ij(:, :, ipair)
      hgamma_ji = this%hgamma_ji(:, :, ipair)
      z_grid(1) = z

      select case (this%selection)
      case ('block')
         gphase = cmplx(0.0_rp, 0.0_rp, rp)
         if (this%left_atoms(ipair) == this%right_atoms(ipair)) then
            call this%green_obj%bgreen_complex(g_one, (ipair - 1)*4 + 1, z_grid, &
               this%a_inf(:, :, 1, ipair), this%b_inf(:, :, 1, ipair), .false.)
            gphase(:, :, 1) = g_one(:, :, 1)
         else
            do phase = 1, 4
               g_one = cmplx(0.0_rp, 0.0_rp, rp)
               call this%green_obj%bgreen_complex(g_one, (ipair - 1)*4 + phase, z_grid, &
                  this%a_inf(:, :, phase, ipair), this%b_inf(:, :, phase, ipair), .false.)
               gphase(:, :, phase) = g_one(:, :, 1)
            end do
         end if
      case ('chebyshev')
         call this%green_obj%chebyshev_green_complex_ij((ipair - 1)*4 + 1, z, gphase)
      end select

      if (this%left_atoms(ipair) == this%right_atoms(ipair)) then
         gij = gphase(:, :, 1)
         gji = gij
      else
         gij = gphase(:, :, 1) - gphase(:, :, 2) + (gphase(:, :, 3)/i_unit - gphase(:, :, 4)/i_unit)
         gji = gphase(:, :, 1) - gphase(:, :, 2) - (gphase(:, :, 3)/i_unit - gphase(:, :, 4)/i_unit)
         gij = 0.5_rp*gij
         gji = 0.5_rp*gji
      end if
   end subroutine tddft_native_rsgf_get_pair

   function tddft_native_rsgf_describe(this) result(description)
      class(tddft_native_rsgf_provider), intent(in) :: this
      character(len=256) :: description

      if (this%selection == 'block') then
         write(description, '(a,a,a,i0,a,i0,a)') 'provider=native_block_recursion selection=', trim(this%selection), &
            ' recursion_depth=', this%recursion_depth, ' pairs=', size(this%pairs), ' workset=serial_replicated'
      else
         write(description, '(a,a,a,i0,a,i0,a)') 'provider=native_chebyshev selection=', trim(this%selection), &
            ' polynomial_order=', this%polynomial_order, ' pairs=', size(this%pairs), ' workset=serial_replicated'
      end if
   end function tddft_native_rsgf_describe

   subroutine native_h_eff_block(this, source_atom, destination_atom, block)
      class(tddft_native_rsgf_provider), intent(inout) :: this
      integer, intent(in) :: source_atom, destination_atom
      complex(rp), intent(out) :: block(:, :)
      integer :: destination_type, intermediate_atom, intermediate_type, direct_slot, intermediate_slot, ineigh
      integer :: neighbor_count

      if (source_atom < 1 .or. source_atom > this%lattice_obj%kk .or. &
          destination_atom < 1 .or. destination_atom > this%lattice_obj%kk) then
         error stop 'native RSGF production provider: H_eff action atom is outside the cluster'
      end if
      if (.not. allocated(this%hamiltonian_obj%ee) .or. .not. allocated(this%hamiltonian_obj%eeo) .or. &
          .not. allocated(this%hamiltonian_obj%enim)) then
         error stop 'native RSGF production provider: accepted H_eff blocks are incomplete'
      end if
      destination_type = this%lattice_obj%iz(destination_atom)
      block = cmplx(0.0_rp, 0.0_rp, rp)
      direct_slot = direct_hopping_slot(this%lattice_obj, destination_atom, source_atom)
      if (direct_slot > 0) block = this%hamiltonian_obj%ee(:, :, direct_slot, destination_type)

      ! H_eff = h - h*o*h + E_nu, written in the same directed-bond
      ! convention used by recursion%ham_hoh_vec_matmul and H(k) assembly.
      neighbor_count = this%lattice_obj%nn(destination_atom, 1)
      do ineigh = 1, neighbor_count
         if (ineigh == 1) then
            intermediate_atom = destination_atom
         else
            intermediate_atom = this%lattice_obj%nn(destination_atom, ineigh)
         end if
         if (intermediate_atom < 1 .or. intermediate_atom > this%lattice_obj%kk) cycle
         intermediate_type = this%lattice_obj%iz(intermediate_atom)
         intermediate_slot = direct_hopping_slot(this%lattice_obj, intermediate_atom, source_atom)
         if (intermediate_slot > 0) then
            block = block - matmul(this%hamiltonian_obj%eeo(:, :, ineigh, destination_type), &
               this%hamiltonian_obj%ee(:, :, intermediate_slot, intermediate_type))
         end if
      end do
      if (destination_atom == source_atom) block = block + this%hamiltonian_obj%enim(:, :, destination_type)
   end subroutine native_h_eff_block

   integer function direct_hopping_slot(lattice_obj, destination_atom, source_atom) result(slot)
      type(lattice), intent(in) :: lattice_obj
      integer, intent(in) :: destination_atom, source_atom
      integer :: ineigh

      slot = 0
      if (destination_atom == source_atom) then
         slot = 1
         return
      end if
      do ineigh = 2, lattice_obj%nn(destination_atom, 1)
         if (lattice_obj%nn(destination_atom, ineigh) == source_atom) then
            slot = ineigh
            return
         end if
      end do
   end function direct_hopping_slot

   integer function recursion_site_slot(lattice_obj, atom) result(slot)
      type(lattice), intent(in) :: lattice_obj
      integer, intent(in) :: atom
      integer :: i

      slot = 0
      if (.not. allocated(lattice_obj%irec)) return
      do i = 1, size(lattice_obj%irec)
         if (lattice_obj%irec(i) == atom) then
            slot = i
            return
         end if
      end do
   end function recursion_site_slot

   subroutine remove_radial_enu(block, basis)
      complex(rp), intent(inout) :: block(:, :)
      type(lmto_radial_basis), intent(in) :: basis
      integer :: l, orbital, ispin, first_orbital, last_orbital, index

      do ispin = 1, 2
         do l = 0, basis%lmax
            first_orbital = l*l + 1
            last_orbital = (l + 1)*(l + 1)
            do orbital = first_orbital, last_orbital
               index = orbital + merge(0, spin_off, ispin == 1)
               block(index, index) = block(index, index) - basis%enu_work(l + 1, ispin)
            end do
         end do
      end do
   end subroutine remove_radial_enu

   subroutine pair_translation(reciprocal_obj, lattice_obj, left_atom, right_atom, translation)
      type(reciprocal), intent(in) :: reciprocal_obj
      type(lattice), intent(in) :: lattice_obj
      integer, intent(in) :: left_atom, right_atom
      real(rp), intent(out) :: translation(3)
      integer :: ineigh, ntype

      translation = 0.0_rp
      if (left_atom == right_atom) return
      if (.not. allocated(reciprocal_obj%ham_vec_type_direct)) then
         error stop 'native RSGF production provider: reciprocal direct-translation metadata is unavailable'
      end if
      ntype = lattice_obj%iz(left_atom)
      do ineigh = 1, lattice_obj%nn(left_atom, 1)
         if (lattice_obj%nn(left_atom, ineigh) == right_atom) then
            if (ineigh > size(reciprocal_obj%ham_vec_type_direct, 2) .or. &
                ntype > size(reciprocal_obj%ham_vec_type_direct, 3)) exit
            translation = reciprocal_obj%ham_vec_type_direct(:, ineigh, ntype)
            return
         end if
      end do
      error stop 'native RSGF production provider: requested pair is not represented by a direct lattice bond'
   end subroutine pair_translation

   integer function response_site_for_atom(lattice_obj, atom, nsite) result(site)
      type(lattice), intent(in) :: lattice_obj
      integer, intent(in) :: atom, nsite
      integer :: isite

      site = 0
      do isite = 1, min(nsite, lattice_obj%nrec)
         if (lattice_obj%irec(isite) == atom) then
            site = isite
            return
         end if
      end do
      if (allocated(lattice_obj%iz)) then
         do isite = 1, min(nsite, lattice_obj%nrec)
            if (lattice_obj%iz(lattice_obj%irec(isite)) == lattice_obj%iz(atom)) then
               site = isite
               return
            end if
         end do
      end if
      error stop 'native RSGF production provider: pair atom has no accepted response-site representative'
   end function response_site_for_atom

   integer function find_pair(this, left_site, right_site, translation) result(ipair)
      class(tddft_native_rsgf_provider), intent(in) :: this
      integer, intent(in) :: left_site, right_site
      real(rp), intent(in) :: translation(3)
      integer :: candidate

      ipair = 0
      do candidate = 1, size(this%pairs)
         if (this%pairs(candidate)%left_site == left_site .and. this%pairs(candidate)%right_site == right_site .and. &
             maxval(abs(this%pairs(candidate)%translation - translation)) <= 2.0e-10_rp) then
            ipair = candidate
            return
         end if
      end do
      error stop 'native RSGF production provider: response pair is not in the registered complete pair set'
   end function find_pair

   subroutine validate_atom_index(lattice_obj, atom, side)
      type(lattice), intent(in) :: lattice_obj
      integer, intent(in) :: atom
      character(len=*), intent(in) :: side
      if (atom < 1 .or. atom > lattice_obj%kk) then
         error stop 'native RSGF production provider: '//trim(side)//' pair atom is outside the cluster'
      end if
   end subroutine validate_atom_index

end module tddft_native_rsgf_provider_mod
