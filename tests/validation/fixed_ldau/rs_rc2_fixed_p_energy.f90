! Evaluate fixed-P Liechtenstein Hubbard U/J derivatives for C3 RC2.
! Input rows: spin, 1-based d row/column, Re(P), Im(P), exactly 50 rows.
! The same production hamiltonian%evaluate_hubbard_energy implementation used
! by the solver is called at U/J +/- 0.05 and +/- 0.10 eV.
program rs_rc2_fixed_p_energy
   use precision_mod, only: rp
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use logger_mod, only: g_logger
   implicit none

   real(rp), parameter :: ry2ev = 13.605703976_rp
   real(rp), parameter :: h_ev(2) = [0.10_rp, 0.05_rp]
   type(hamiltonian), target :: ham
   type(lattice), target :: lat
   complex(rp) :: density(1,4,2,7,7), potential(1,4,2,7,7)
   real(rp) :: e_int(1,4), e_dc(1,4), e_hub(1,4), tr_vp(1,4)
   real(rp) :: preal, pimag, u_center_ev, j_center_ev, plus_u, minus_u, plus_j, minus_j, center_ehub
   integer :: unit, ios, spin, i, j, n, arg_status
   character(len=1024) :: pfile, arg

   call get_command_argument(1, pfile)
   if (len_trim(pfile) == 0) error stop 'usage: driver P_matrix.txt U_center_eV J_center_eV'
   call get_command_argument(2, arg, status=arg_status)
   if (arg_status /= 0) error stop 'missing U_center_eV'
   read(arg, *, iostat=ios) u_center_ev
   if (ios /= 0) error stop 'invalid U_center_eV'
   call get_command_argument(3, arg, status=arg_status)
   if (arg_status /= 0) error stop 'missing J_center_eV'
   read(arg, *, iostat=ios) j_center_ev
   if (ios /= 0) error stop 'invalid J_center_eV'

   call g_logger%init()
   lat%ntype = 1
   allocate(lat%symbolic_atoms(1))
   call lat%symbolic_atoms(1)%restore_to_default()
   ham%lattice => lat
   ham%hubbard_u_potential_form = 'liechtenstein'
   density = cmplx(0.0_rp, 0.0_rp, rp)
   potential = cmplx(0.0_rp, 0.0_rp, rp)
   open(newunit=unit, file=trim(pfile), status='old', action='read', iostat=ios)
   if (ios /= 0) error stop 'cannot open P_matrix.txt'
   do n = 1, 50
      read(unit, *, iostat=ios) spin, i, j, preal, pimag
      if (ios /= 0 .or. spin < 1 .or. spin > 2 .or. i < 1 .or. i > 5 .or. j < 1 .or. j > 5) &
         error stop 'invalid P matrix row'
      density(1,3,spin,i,j) = cmplx(preal, pimag, rp)
   end do
   close(unit)

   do n = 1, size(h_ev)
      call set_uj(u_center_ev + h_ev(n), j_center_ev)
      call ham%evaluate_hubbard_energy(density, potential, e_int, e_dc, e_hub, tr_vp)
      plus_u = e_hub(1,3)
      call set_uj(u_center_ev - h_ev(n), j_center_ev)
      call ham%evaluate_hubbard_energy(density, potential, e_int, e_dc, e_hub, tr_vp)
      minus_u = e_hub(1,3)
      call set_uj(u_center_ev, j_center_ev + h_ev(n))
      call ham%evaluate_hubbard_energy(density, potential, e_int, e_dc, e_hub, tr_vp)
      plus_j = e_hub(1,3)
      call set_uj(u_center_ev, j_center_ev - h_ev(n))
      call ham%evaluate_hubbard_energy(density, potential, e_int, e_dc, e_hub, tr_vp)
      minus_j = e_hub(1,3)
      call set_uj(u_center_ev, j_center_ev)
      call ham%evaluate_hubbard_energy(density, potential, e_int, e_dc, e_hub, tr_vp)
      center_ehub = e_hub(1,3)
      write(*,'(a,es24.16,a,es24.16,a,es24.16,a,es24.16,a,es24.16,a,es24.16)') &
         'h_eV=', h_ev(n), ' dEHub_dU_Ry_per_eV=', (plus_u-minus_u)/(2.0_rp*h_ev(n)), &
         ' dEHub_dJ_Ry_per_eV=', (plus_j-minus_j)/(2.0_rp*h_ev(n)), &
         ' EHub_Ucenter_Ry=', center_ehub, ' EHub_Uplus_Ry=', plus_u, ' EHub_Uminus_Ry=', minus_u
      write(*,'(a,es24.16,a,es24.16)') ' EHub_Jplus_Ry=', plus_j, ' EHub_Jminus_Ry=', minus_j
   end do

contains
   subroutine set_uj(u, jbar)
      real(rp), intent(in) :: u, jbar
      lat%symbolic_atoms(1)%potential%hubbard_u(3) = u / ry2ev
      lat%symbolic_atoms(1)%potential%hubbard_j(3) = jbar / ry2ev
   end subroutine set_uj
end program rs_rc2_fixed_p_energy
