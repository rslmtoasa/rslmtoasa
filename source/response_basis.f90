!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Charge-spin response operators in the RS-LMTO spin-orbital basis.
!>
!> The normal electronic basis is spin-major: (orbital_1 up, ..., orbital_N
!> up, orbital_1 down, ..., orbital_N down).  This module deliberately takes
!> the orbital count as an argument so response code cannot silently inherit a
!> BdG/Nambu extension when it needs an electronic 2N basis.
module response_basis_mod
   use precision_mod, only: rp
   use math_mod, only: i_unit
   use tddft_conventions_mod, only: tddft_circular_operator_factor, tddft_circular_source_factor
   use response_components_mod, only: RESPONSE_CHARGE, RESPONSE_MX, RESPONSE_MY, RESPONSE_MZ, &
      RESPONSE_PLUS, RESPONSE_MINUS, is_response_component
   implicit none

   private

   public :: response_operator, ladder_operator, external_source_operator

contains

   !> Return the operator whose expectation is the named response component.
   !> For PLUS/MINUS this is sigma_x +/- i sigma_y, so its expectation is
   !> exactly m^+/- = m_x +/- i m_y (not a half-normalized ladder operator).
   function response_operator(component, orbital_count) result(op)
      integer, intent(in) :: component, orbital_count
      complex(rp), allocatable :: op(:, :)
      integer :: orbital

      if (.not. is_response_component(component)) error stop 'response_operator: invalid component'
      if (orbital_count <= 0) error stop 'response_operator: orbital_count must be positive'

      allocate(op(2*orbital_count, 2*orbital_count))
      op = cmplx(0.0_rp, 0.0_rp, rp)

      do orbital = 1, orbital_count
         select case (component)
         case (RESPONSE_CHARGE)
            op(orbital, orbital) = cmplx(1.0_rp, 0.0_rp, rp)
            op(orbital + orbital_count, orbital + orbital_count) = cmplx(1.0_rp, 0.0_rp, rp)
         case (RESPONSE_MX)
            op(orbital, orbital + orbital_count) = cmplx(1.0_rp, 0.0_rp, rp)
            op(orbital + orbital_count, orbital) = cmplx(1.0_rp, 0.0_rp, rp)
         case (RESPONSE_MY)
            op(orbital, orbital + orbital_count) = -i_unit
            op(orbital + orbital_count, orbital) = i_unit
         case (RESPONSE_MZ)
            op(orbital, orbital) = cmplx(1.0_rp, 0.0_rp, rp)
            op(orbital + orbital_count, orbital + orbital_count) = cmplx(-1.0_rp, 0.0_rp, rp)
         case (RESPONSE_PLUS)
            op(orbital, orbital + orbital_count) = cmplx(tddft_circular_operator_factor, 0.0_rp, rp)
         case (RESPONSE_MINUS)
            op(orbital + orbital_count, orbital) = cmplx(tddft_circular_operator_factor, 0.0_rp, rp)
         end select
      end do
   end function response_operator

   !> Return conventional sigma^+/- = (sigma_x +/- i sigma_y)/2.
   !> Use response_operator for the physical m^+/- measurement convention.
   function ladder_operator(component, orbital_count) result(op)
      integer, intent(in) :: component, orbital_count
      complex(rp), allocatable :: op(:, :)

      if (component /= RESPONSE_PLUS .and. component /= RESPONSE_MINUS) then
         error stop 'ladder_operator: component must be PLUS or MINUS'
      end if
      op = tddft_circular_source_factor*response_operator(component, orbital_count)
   end function ladder_operator

   !> Return dH/dB for an external magnetic-field coordinate.
   !>
   !> Cartesian field coordinates use sigma_x/y directly.  A circular field
   !> coordinate is paired with the opposite measured circular operator:
   !> dH/dB+ = (sigma_x-i sigma_y)/2 and dH/dB- = (sigma_x+i sigma_y)/2.
   !> This is a source vertex only; it is not the self-consistent XC kernel.
   function external_source_operator(field_component, orbital_count) result(op)
      integer, intent(in) :: field_component, orbital_count
      complex(rp), allocatable :: op(:, :)

      select case (field_component)
      case (RESPONSE_PLUS)
         op = tddft_circular_source_factor*response_operator(RESPONSE_MINUS, orbital_count)
      case (RESPONSE_MINUS)
         op = tddft_circular_source_factor*response_operator(RESPONSE_PLUS, orbital_count)
      case default
         op = response_operator(field_component, orbital_count)
      end select
   end function external_source_operator

end module response_basis_mod
