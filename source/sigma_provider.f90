!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: sigma_provider
!
!> @author
!> Anders Bergman
!
! DESCRIPTION:
!> Self-energy provider abstraction for the k-space Green’s-function engine
!> (milestone B2, `reciprocal_green`). A `sigma_provider` supplies the
!> retarded self-energy block Sigma(z) on an arbitrary complex energy z for a
!> given unit-cell site, in the same (screened/auxiliary) LMTO representation
!> and orbital ordering as the real-space Green’s-function route.
!>
!> The abstract type is defined here (not in a submodule) so that the Dyson
!> backend in `reciprocal_green` and future concrete providers can share one
!> binding contract:
!>
!>   - `sigma_zero`   : Sigma = 0 (this module) -> makes backend D reduce to
!>                      backend E, the permanent two-backend CI invariant.
!>   - `sigma_cpa`    : single-site coherent potential (B8).
!>   - `sigma_static` : DFT+U static shift (B10 plumbing test).
!>   - `sigma_file`   : external DMFT solver exchange (B10).
!>
!> Keeping this in its own small module means those providers depend only on
!> the interface, never on the whole `reciprocal_mod`, avoiding a dependency
!> cycle (providers `use sigma_provider_mod`; `reciprocal_mod` `use`s it too).
!------------------------------------------------------------------------------
module sigma_provider_mod
   use precision_mod, only: rp
   implicit none

   private
   public :: sigma_provider, sigma_zero

   !> Abstract retarded self-energy provider.
   type, abstract :: sigma_provider
   contains
      !> Return Sigma(z) on the nb x nb block of unit-cell site `isite`.
      procedure(sigma_get_block), deferred :: get_sigma
   end type sigma_provider

   abstract interface
      !> @brief Evaluate the self-energy block for one site at complex energy z.
      !> @param[in]  this        Concrete self-energy provider.
      !> @param[in]  z           Complex energy on the retarded contour (z = E + i*eta).
      !> @param[in]  isite       Unit-cell site index (1..nrec).
      !> @param[out] sigma_block Sigma(z) block, orbital ordering matching the
      !>                         RS route (nb x nb, spin-major Pauli blocks).
      subroutine sigma_get_block(this, z, isite, sigma_block)
         import :: sigma_provider, rp
         class(sigma_provider), intent(in) :: this
         complex(rp), intent(in) :: z
         integer, intent(in) :: isite
         complex(rp), dimension(:, :), intent(out) :: sigma_block
      end subroutine sigma_get_block
   end interface

   !> Trivial provider: Sigma(z) = 0 for every site and energy.
   !> This is the default for backend E (strict Lehmann) and the reference that
   !> makes backend D (Dyson) reproduce backend E to solver tolerance.
   type, extends(sigma_provider) :: sigma_zero
   contains
      procedure :: get_sigma => sigma_zero_get
   end type sigma_zero

contains

   !> @brief sigma_zero implementation: fill the block with zeros.
   !> @param[in]  this        sigma_zero provider (unused).
   !> @param[in]  z           Complex energy (unused).
   !> @param[in]  isite       Site index (unused).
   !> @param[out] sigma_block Zeroed self-energy block.
   subroutine sigma_zero_get(this, z, isite, sigma_block)
      class(sigma_zero), intent(in) :: this
      complex(rp), intent(in) :: z
      integer, intent(in) :: isite
      complex(rp), dimension(:, :), intent(out) :: sigma_block

      sigma_block = cmplx(0.0_rp, 0.0_rp, rp)
   end subroutine sigma_zero_get

end module sigma_provider_mod
