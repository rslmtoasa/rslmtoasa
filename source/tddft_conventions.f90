!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Small executable part of the TDDFT response convention contract.
!>
!> The full derivation is documented in docs/TDDFT_CONVENTIONS.md.  Keeping
!> these factors and the retarded denominator in one module makes a sign or
!> normalization change visible to every response backend and to its tests.
module tddft_conventions_mod
   use precision_mod, only: rp
   implicit none

   private

   ! The measured circular variables are m+/- = mx +/- i my.  Therefore their
   ! matrices are sigma_x +/- i sigma_y = 2 sigma_+/-.
   real(rp), parameter, public :: tddft_circular_operator_factor = 2.0_rp
   ! In the Hamiltonian, B+/- couples to the opposite unhalved circular
   ! operator with one half: dH = (B+ O- + B- O+)/2.
   real(rp), parameter, public :: tddft_circular_source_factor = 0.5_rp

   public :: tddft_retarded_denominator
   public :: tddft_advanced_denominator
   public :: tddft_retarded_green_denominator
   public :: tddft_advanced_green_denominator

contains

   !> Retarded Lehmann denominator used by every dynamic chi_KS backend.
   elemental pure function tddft_retarded_denominator(omega, transition_energy, eta) result(denominator)
      real(rp), intent(in) :: omega, transition_energy, eta
      complex(rp) :: denominator

      denominator = cmplx(omega + transition_energy, eta, rp)
   end function tddft_retarded_denominator

   !> Advanced counterpart used by convention/symmetry fixtures.
   elemental pure function tddft_advanced_denominator(omega, transition_energy, eta) result(denominator)
      real(rp), intent(in) :: omega, transition_energy, eta
      complex(rp) :: denominator

      denominator = cmplx(omega + transition_energy, -eta, rp)
   end function tddft_advanced_denominator

   !> Retarded one-particle Green-function denominator.
   elemental pure function tddft_retarded_green_denominator(energy, eigenvalue, eta) result(denominator)
      real(rp), intent(in) :: energy, eigenvalue, eta
      complex(rp) :: denominator

      denominator = cmplx(energy-eigenvalue, eta, rp)
   end function tddft_retarded_green_denominator

   !> Advanced one-particle Green-function denominator.
   elemental pure function tddft_advanced_green_denominator(energy, eigenvalue, eta) result(denominator)
      real(rp), intent(in) :: energy, eigenvalue, eta
      complex(rp) :: denominator

      denominator = cmplx(energy-eigenvalue, -eta, rp)
   end function tddft_advanced_green_denominator

end module tddft_conventions_mod
