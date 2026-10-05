!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!> @brief spin_response C1 (ground-state handoff) on the accepted bcc-Fe state.
!> @details Usage, from a scratch copy of tests/spin_response/fe_bcc:
!>          test_spin_response_c1 <oracle directory>. Tolerances come from <dir>/tolerances.nml.
!------------------------------------------------------------------------------
program test_spin_response_c1
   use precision_mod, only: rp
   use control_mod, only: control
   use lattice_mod, only: lattice
   use charge_mod, only: charge
   use energy_mod, only: energy
   use hamiltonian_mod, only: hamiltonian
   use reciprocal_mod, only: reciprocal
   use spin_response_mod, only: spin_response
   use spin_response_fe_state
   implicit none

   type(control), target :: control_obj
   type(lattice), target :: lattice_obj
   type(charge), target :: charge_obj
   type(energy), target :: energy_obj
   type(hamiltonian), target :: ham
   type(reciprocal), target :: rec, rec_ref, rec_fixed
   type(spin_response) :: sr, sr_fixed
   character(len=512) :: dir
   real(rp), allocatable :: e(:, :)
   complex(rp), allocatable :: v(:, :, :)
   real(rp) :: q(3), n_ref, eband, dmax, dshift
   integer :: nk, first, last, ik, m

   call get_command_argument(1, dir)
   if (len_trim(dir) == 0) error stop 'usage: test_spin_response_c1 <oracle directory>'
   call read_tolerances(trim(dir))
   call build_fe_state(.false., control_obj, lattice_obj, charge_obj, energy_obj, ham)

   ! Pre-campaign mesh diagonalization, on its own object, is the reference for tests 1 and 2.
   rec_ref = reciprocal(ham)
   call rec_ref%generate_mp_mesh()
   call rec_ref%build_kspace_hamiltonian()
   call rec_ref%diagonalize_hamiltonian()
   nk = rec_ref%k_workset%nk_local

   rec = reciprocal(ham)
   rec%fermi_level = energy_obj%fermi
   sr = spin_response(rec)
   call sr%prepare()

   ! C1.1: arbitrary-k service on the mesh points vs the mesh diagonalization.
   ! Catches the campaign-era service drifting from the pre-campaign H(k) assembly.
   dmax = 0.0_rp
   dshift = 0.0_rp
   do first = 1, nk, 64
      last = min(nk, first + 63)
      call rec%calculate_eigenpairs_at_kpoints(rec_ref%k_workset%points(:, first:last), e, v)
      dmax = max(dmax, maxval(abs(e - rec_ref%eigenvalues(:, first:last))))
   end do
   call report('C1.1 service vs mesh diagonalization, max |de| (Ry)', dmax, eigen_contract_abs)

   ! C1.2: q = a mesh vector; e(k+q) must be the mesh eigenvalues at the shifted index.
   ! Catches q ignored or applied wrongly; the second line shows the test could fail.
   q = rec_ref%k_workset%points(:, 2) - rec_ref%k_workset%points(:, 1)
   dmax = 0.0_rp
   dshift = 0.0_rp
   do first = 1, nk, 64
      last = min(nk, first + 63)
      call sr%eigenpairs_chunk(q, first, last, e, v)
      do ik = first, last
         m = mesh_index(rec_ref%k_workset%points, rec_ref%k_workset%points(:, ik) + q)
         dmax = max(dmax, maxval(abs(e(:, ik - first + 1) - rec_ref%eigenvalues(:, m))))
         dshift = max(dshift, maxval(abs(rec_ref%eigenvalues(:, ik) - rec_ref%eigenvalues(:, m))))
      end do
   end do
   call report('C1.2 e(k+q) vs mesh eigenvalue at k+q, max |de| (Ry)', dmax, eigen_contract_abs)
   call check('C1.2 shifted point differs from k (q not ignored), max |de| > 1e-3 Ry', dshift > 1.0e-3_rp)

   ! C1.3: auto_find_fermi off, the deck E_F: used E_F equals the input bit for bit.
   ! Catches the Fermi search running anyway. N at the deck value is reported.
   rec_fixed = reciprocal(ham)
   rec_fixed%auto_find_fermi = .false.
   rec_fixed%fermi_level = energy_obj%fermi
   sr_fixed = spin_response(rec_fixed)
   call sr_fixed%prepare()
   call check('C1.3 used E_F == input E_F', sr_fixed%fermi_used == energy_obj%fermi)
   write (*, '(a,f10.7,a,es10.3)') 'C1.3 at deck E_F ', energy_obj%fermi, ' Ry: N - target = ', &
      sr_fixed%electron_count - rec_fixed%total_electrons

   ! C1.4: eps(-k) = eps(k) on the mesh (bcc Fe is centrosymmetric).
   ! Catches a mesh or folding error that breaks inversion.
   dmax = 0.0_rp
   do ik = 1, nk
      m = mesh_index(rec%k_workset%points, -rec%k_workset%points(:, ik))
      dmax = max(dmax, maxval(abs(rec%eigenvalues(:, ik) - rec%eigenvalues(:, m))))
   end do
   call report('C1.4 eps(-k) vs eps(k), max |de| (Ry)', dmax, eigen_contract_abs)

   ! C1.5: N on the response mesh vs the target, and vs the pre-campaign occupation sum at the used E_F.
   ! Catches a wrong E_F or k-weight normalization entering the occupations.
   call rec%evaluate_eigenvalue_occupations(sr%fermi_used, n_ref, eband)
   call report('C1.5 |N - target| (electrons)', abs(sr%electron_count - rec%total_electrons), electron_count_abs)
   call report('C1.5 |N - pre-campaign N(E_F used)| (electrons)', abs(sr%electron_count - n_ref), electron_count_abs)

   call finish()

contains

   !> Index of the mesh point equal to p modulo reciprocal lattice vectors.
   integer function mesh_index(points, p) result(m)
      real(rp), intent(in) :: points(:, :), p(3)

      real(rp) :: d(3)

      do m = 1, size(points, 2)
         d = points(:, m) - p
         if (all(abs(d - nint(d)) < 1.0e-8_rp)) return
      end do
      error stop 'mesh_index: k+q is not a mesh point'
   end function mesh_index

end program test_spin_response_c1
