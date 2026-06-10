module hamiltonian_base_mod
   use precision_mod, only: rp
   implicit none

   !> Abstract Hamiltonian type
   type, abstract, public :: hamiltonian_type
      logical :: hoh_enabled ! Flag to trigger double-step
   contains
      procedure(I_apply), deferred, pass :: apply
   end type hamiltonian_type

   abstract interface
      !> Generic interface for H|Psi>
      subroutine I_apply(this, psi_in, psi_out, a, b)
         import :: hamiltonian_type, rp
         class(hamiltonian_type), intent(inout) :: this
         complex(rp), intent(in)                :: psi_in(:,:,:)
         complex(rp), intent(out)               :: psi_out(:,:,:)
         real(rp), intent(in)                   :: a
         real(rp), intent(in)                   :: b
      end subroutine I_apply
   end interface
end module hamiltonian_base_mod