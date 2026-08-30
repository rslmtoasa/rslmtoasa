!------------------------------------------------------------------------------
! RS-LMTO-ASA
!
! SUBMODULE: Self XC response helpers
!
!> Ground-state XC-response refresh and MPI synchronization helpers for `self`.
!------------------------------------------------------------------------------

submodule(self_mod) self_xc_response

   implicit none

contains

   !> Re-evaluate the atomic SCF XC path for the current ground-state
   !> potential and refresh the provider.  This is intended for a response
   !> post-processing consumer which is constructed after the normal SCF
   !> object has gone out of scope; it does not infer a kernel from Hxc.
   module subroutine refresh_xc_response_kernel(this)
      class(self), intent(inout) :: this
      call this%run_scf()
   end subroutine refresh_xc_response_kernel

   module subroutine synchronize_xc_response_provider(this)
      class(self), intent(inout) :: this
      real(rp), allocatable :: packed(:, :)
      integer :: ia

#ifdef USE_MPI
      allocate(packed(7, this%lattice%nrec))
      packed = 0.0_rp
      do ia = 1, this%lattice%nrec
         if (this%xc_response_provider%site(ia)%has_radial_projection) then
            packed(1, ia) = this%xc_response_provider%site(ia)%radial_spin_population
            packed(2, ia) = this%xc_response_provider%site(ia)%vxc_scalar
            packed(3, ia) = this%xc_response_provider%site(ia)%bxc_spin_moment
            packed(4, ia) = this%xc_response_provider%site(ia)%radial_spin_abs_population
            packed(5, ia) = this%xc_response_provider%site(ia)%radial_vxc_spin_difference
            packed(6, ia) = this%xc_response_provider%site(ia)%radial_vxc_spin_difference_abs
            packed(7, ia) = 1.0_rp
         end if
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE, packed, product(shape(packed)), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      do ia = 1, this%lattice%nrec
         if (packed(7, ia) > 1.5_rp) then
            call g_logger%fatal('SCF XC response provider received duplicate site records across MPI ranks.', __FILE__, __LINE__)
         end if
         if (packed(7, ia) > 0.5_rp) then
            this%xc_response_provider%site(ia)%radial_spin_population = packed(1, ia)
            this%xc_response_provider%site(ia)%vxc_scalar = packed(2, ia)
            this%xc_response_provider%site(ia)%bxc_spin_moment = packed(3, ia)
            this%xc_response_provider%site(ia)%radial_spin_abs_population = packed(4, ia)
            this%xc_response_provider%site(ia)%radial_vxc_spin_difference = packed(5, ia)
            this%xc_response_provider%site(ia)%radial_vxc_spin_difference_abs = packed(6, ia)
            this%xc_response_provider%site(ia)%has_radial_projection = .true.
            this%xc_response_provider%site(ia)%has_k_perp_circular = .false.
         end if
      end do
      deallocate(packed)
#endif

      ! This is the same site population used by P_site sigma response
      ! vertices.  Its setter now normalizes the radial ALSDA numerator to
      ! that projector, rather than replacing it by a radial atomic quantity.
      do ia = 1, this%lattice%nrec
         call this%xc_response_provider%set_site_spin_population(ia, &
            this%symbolic_atom(this%lattice%nbulk + ia)%potential%mtot)
         call this%xc_response_provider%set_site_signed_spin_population(ia, &
            this%symbolic_atom(this%lattice%nbulk + ia)%potential%mtot * &
            this%symbolic_atom(this%lattice%nbulk + ia)%potential%mom(3))
         call this%xc_response_provider%set_site_magnetization_direction(ia, &
            this%symbolic_atom(this%lattice%nbulk + ia)%potential%mom(1:3))
      end do
   end subroutine synchronize_xc_response_provider

end submodule self_xc_response
