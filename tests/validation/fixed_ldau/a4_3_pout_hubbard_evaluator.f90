! Re-evaluate Hubbard components from the captured direct P_OUT density using
! the production Liechtenstein evaluate_hubbard_energy routine.
! P file rows: 1-based spin,row,column,real(P),imag(P).
program a4_3_pout_hubbard_evaluator
   use precision_mod, only: rp
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use logger_mod, only: g_logger
   implicit none
   real(rp), parameter :: ry2ev = 13.605703976_rp
   type(hamiltonian), target :: ham
   type(lattice), target :: lat
   complex(rp) :: density(1,4,2,7,7), potential(1,4,2,7,7)
   real(rp) :: e_int(1,4), e_dc(1,4), e_hub(1,4), tr_vp(1,4)
   real(rp) :: u_ev, j_ev, re_p, im_p
   integer :: unit, ios, spin, i, j, n
   character(len=1024) :: pfile, arg

   call get_command_argument(1, pfile)
   if (len_trim(pfile) == 0) error stop 'usage: binary P_matrix.txt U_eV J_eV'
   call get_command_argument(2, arg)
   if (len_trim(arg) == 0) error stop 'missing U_eV'
   read(arg, *, iostat=ios) u_ev
   if (ios /= 0) error stop 'invalid U_eV'
   call get_command_argument(3, arg)
   if (len_trim(arg) == 0) error stop 'missing J_eV'
   read(arg, *, iostat=ios) j_ev
   if (ios /= 0) error stop 'invalid J_eV'

   call g_logger%init()
   lat%ntype = 1
   allocate(lat%symbolic_atoms(1))
   call lat%symbolic_atoms(1)%restore_to_default()
   ham%lattice => lat
   ham%hubbard_u_potential_form = 'liechtenstein'
   lat%symbolic_atoms(1)%potential%hubbard_u(3) = u_ev / ry2ev
   lat%symbolic_atoms(1)%potential%hubbard_j(3) = j_ev / ry2ev
   density = cmplx(0.0_rp, 0.0_rp, rp)
   potential = cmplx(0.0_rp, 0.0_rp, rp)

   open(newunit=unit, file=trim(pfile), status='old', action='read', iostat=ios)
   if (ios /= 0) error stop 'cannot open P_matrix.txt'
   do n = 1, 50
      read(unit, *, iostat=ios) spin, i, j, re_p, im_p
      if (ios /= 0 .or. spin < 1 .or. spin > 2 .or. i < 1 .or. i > 5 .or. j < 1 .or. j > 5) &
         error stop 'invalid P matrix row'
      density(1,3,spin,i,j) = cmplx(re_p, im_p, rp)
   end do
   close(unit)

   call ham%evaluate_hubbard_energy(density, potential, e_int, e_dc, e_hub, tr_vp)
   write(*,'(a,es24.16,a,es24.16,a,es24.16,a,es24.16,a,es24.16,a,es24.16)') &
      'U_eV=', u_ev, ' J_eV=', j_ev, ' E_int_Ry=', e_int(1,3), ' E_DC_Ry=', e_dc(1,3), &
      ' E_Hub_Ry=', e_hub(1,3), ' t_in_out_Ry=', tr_vp(1,3)
end program a4_3_pout_hubbard_evaluator
