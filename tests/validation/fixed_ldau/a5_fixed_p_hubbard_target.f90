! Evaluate fixed-P Hubbard targets through the production Liechtenstein
! evaluator. Input rows are spin, row, column, real(P), imag(P), all 1-based.
program a5_fixed_p_hubbard_target
   use precision_mod, only: rp
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use logger_mod, only: g_logger
   implicit none

   real(rp), parameter :: ry2ev = 13.605703976_rp
   real(rp), parameter :: eps_ev(3) = [1.0e-2_rp, 5.0e-3_rp, 2.5e-3_rp]
   type(hamiltonian), target :: ham
   type(lattice), target :: lat
   complex(rp) :: density(1,4,2,7,7), potential(1,4,2,7,7)
   real(rp) :: e_int(1,4), e_dc(1,4), e_hub(1,4), tr_vp(1,4)
   real(rp) :: preal, pimag, u_ev, j_ev, j_center_ev, plus_u, minus_u, plus_j, minus_j
   integer :: unit, ios, spin, i, j, n
   character(len=1024) :: pfile

   call get_command_argument(1, pfile)
   if (len_trim(pfile) == 0) error stop 'usage: a5_fixed_p_hubbard_target P_matrix.txt [J_center_eV]'
   j_center_ev = 0.0_rp
   call get_command_argument(2, pfile)
   if (len_trim(pfile) > 0) read(pfile, *) j_center_ev
   call get_command_argument(1, pfile)
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

   do n = 1, size(eps_ev)
      call set_uj(2.0_rp + eps_ev(n), 0.0_rp)
      call ham%evaluate_hubbard_energy(density, potential, e_int, e_dc, e_hub, tr_vp)
      plus_u = e_hub(1,3)
      call set_uj(2.0_rp - eps_ev(n), 0.0_rp)
      call ham%evaluate_hubbard_energy(density, potential, e_int, e_dc, e_hub, tr_vp)
      minus_u = e_hub(1,3)
      call set_uj(2.0_rp, j_center_ev + eps_ev(n))
      call ham%evaluate_hubbard_energy(density, potential, e_int, e_dc, e_hub, tr_vp)
      plus_j = e_hub(1,3)
      call set_uj(2.0_rp, j_center_ev - eps_ev(n))
      call ham%evaluate_hubbard_energy(density, potential, e_int, e_dc, e_hub, tr_vp)
      minus_j = e_hub(1,3)
      write(*,'(a,es24.16,a,es24.16,a,es24.16,a,es24.16,a,es24.16)') 'eps_eV=', eps_ev(n), &
         ' dEHub_dU_Ry_per_eV=', (plus_u-minus_u)/(2.0_rp*eps_ev(n)), &
         ' J_center_eV=', j_center_ev, ' dEHub_dJ_Ry_per_eV=', &
         (plus_j-minus_j)/(2.0_rp*eps_ev(n)), ' EHub_Ry=', &
         (plus_u+minus_u)/2.0_rp
   end do

contains

   subroutine set_uj(u, jbar)
      real(rp), intent(in) :: u, jbar
      lat%symbolic_atoms(1)%potential%hubbard_u(3) = u / ry2ev
      lat%symbolic_atoms(1)%potential%hubbard_j(3) = jbar / ry2ev
   end subroutine set_uj

end program a5_fixed_p_hubbard_target
