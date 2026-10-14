! TEST-13: lean contracts for the production Liechtenstein Hubbard correction.
!
! These checks intentionally call hamiltonian%calculate_hubbard_u_potential_general
! directly.  The expected values below are only for the l=0 one-orbital case;
! the Hermiticity check exercises the production Coulomb-matrix path for l=1
! without introducing a second implementation of that algebra.
program test_lda_u_hamiltonian
   use precision_mod, only: rp
   use basis_mod, only: basis_init, nb, spin_off
   use control_mod, only: control
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use reciprocal_mod, only: reciprocal
   use logger_mod, only: g_logger
   use math_mod, only: Coulomb_mat_spherical
   use mix_mod, only: mix
   implicit none

   real(rp), parameter :: tolerance = 1.0e-11_rp
   logical :: failed

   type(hamiltonian), target :: ham
   type(lattice), target :: lat
   type(control), target :: ctl
   type(mix) :: mixer

   failed = .false.
   call g_logger%init()
   call basis_init(1)
   call setup_fixture()

      call test_complex_ldm_flatten_roundtrip()
      call test_legacy_restart_without_imaginary_field()
      call test_complex_ldm_mixing_roundtrip()
   call test_zero_u_j_limit()
   call test_simple_diagonal_occupation()
   call test_collinear_hermiticity()
      call test_complex_hermitian_occupation()
      call test_complex_energy_derivative()
   call test_correlated_trace_and_common_state_fourier()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains

   subroutine test_complex_ldm_flatten_roundtrip()
      complex(rp) :: expected
      integer :: unit, ios
      character(len=*), parameter :: filename = 'UnitLdaUComplexRestart.nml'

      lat%symbolic_atoms(1)%potential%ldm = (0.0_rp, 0.0_rp)
      expected = cmplx(0.21_rp, -0.37_rp, rp)
      lat%symbolic_atoms(1)%potential%ldm(2, 1, 1, 3) = expected
      call lat%symbolic_atoms(1)%potential%flatten_ldm()
      call require(abs(lat%symbolic_atoms(1)%potential%ldm_flatten(2, 1, 3) - 0.21_rp) < tolerance, &
         'legacy real restart field retains real LDM component')
      call require(abs(lat%symbolic_atoms(1)%potential%ldm_flatten_imag(2, 1, 3) + 0.37_rp) < tolerance, &
         'paired imaginary restart field retains complex LDM component')
      lat%symbolic_atoms(1)%potential%ldm = (0.0_rp, 0.0_rp)
      call lat%symbolic_atoms(1)%potential%expand_ldm()
      call require(abs(lat%symbolic_atoms(1)%potential%ldm(2, 1, 1, 3) - expected) < tolerance, &
         'complex LDM flatten/expand round trip')

      ! Exercise the serialized restart path, not only the in-memory helpers.
      call lat%symbolic_atoms(1)%potential%flatten_ldm()
      open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
      call require(ios == 0, 'open temporary complex LDM restart for writing')
      if (ios /= 0) return
      call lat%symbolic_atoms(1)%potential%print_state_formatted(unit=unit)
      close(unit)
      lat%symbolic_atoms(1)%potential%ldm_flatten = 0.0_rp
      lat%symbolic_atoms(1)%potential%ldm_flatten_imag = 0.0_rp
      lat%symbolic_atoms(1)%potential%ldm = (0.0_rp, 0.0_rp)
      call lat%symbolic_atoms(1)%potential%build_from_file(filename)
      call require(abs(lat%symbolic_atoms(1)%potential%ldm(2, 1, 1, 3) - expected) < tolerance, &
         'complex LDM survives serialized write/read restart')
      open(newunit=unit, file=filename, status='old', iostat=ios)
      if (ios == 0) close(unit, status='delete')
   end subroutine test_complex_ldm_flatten_roundtrip

   subroutine test_legacy_restart_without_imaginary_field()
      integer :: unit, ios
      character(len=*), parameter :: filename = 'UnitLdaULegacyRestart.nml'

      open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
      call require(ios == 0, 'open temporary legacy LDM restart')
      if (ios /= 0) return
      write(unit, '(a)') '&par'
      write(unit, '(a)') ' lmax=2,'
      write(unit, '(a)') ' ldm_flatten(2,1,3)=0.21'
      write(unit, '(a)') '/'
      close(unit)

      call lat%symbolic_atoms(1)%potential%build_from_file(filename)
      call require(abs(real(lat%symbolic_atoms(1)%potential%ldm(2,1,1,3), rp) - 0.21_rp) < tolerance, &
         'legacy restart restores real LDM entry')
      call require(abs(aimag(lat%symbolic_atoms(1)%potential%ldm(2,1,1,3))) < tolerance, &
         'legacy restart without imaginary field initializes it to zero')
      open(newunit=unit, file=filename, status='old', iostat=ios)
      if (ios == 0) close(unit, status='delete')
   end subroutine test_legacy_restart_without_imaginary_field

   subroutine test_complex_ldm_mixing_roundtrip()
      complex(rp) :: expected
      integer, parameter :: flat_index = 3

      ! The p-shell (l index 2) flattened entry (m=1,m'=3) carries a
      ! nonzero imaginary component through old/new buffers and current save.
      lat%nrec = 1
      lat%nbulk = 0
      mixer%lattice => lat
      mixer%symbolic_atom => lat%symbolic_atoms
      allocate(mixer%ldm_old(1, 2, 2, 9), mixer%ldm_new(1, 2, 2, 9), mixer%ldm_mix(1, 2, 2, 9))
      mixer%ldm_old = (0.0_rp, 0.0_rp)
      mixer%ldm_new = (0.0_rp, 0.0_rp)
      mixer%ldm_beta = 0.25_rp
      mixer%ldm_old(1, 2, 1, flat_index) = cmplx(0.2_rp, -0.4_rp, rp)
      mixer%ldm_new(1, 2, 1, flat_index) = cmplx(0.6_rp, 0.8_rp, rp)
      expected = (1.0_rp - mixer%ldm_beta)*mixer%ldm_old(1, 2, 1, flat_index) + &
         mixer%ldm_beta*mixer%ldm_new(1, 2, 1, flat_index)

      call mixer%mix_ldm_linear()
      call require(abs(mixer%ldm_mix(1, 2, 1, flat_index) - expected) < tolerance, &
         'linear LDM mixing retains the imaginary component')
      call mixer%save_ldm_to('current')
      call require(abs(lat%symbolic_atoms(1)%potential%ldm(2, 1, 1, 3) - expected) < tolerance, &
         'mixed complex LDM is stored in the current potential')
   end subroutine test_complex_ldm_mixing_roundtrip

   subroutine setup_fixture()
      call ctl%restore_to_default()
      ctl%lmax = 1
      ctl%nsp = 2
      lat%ntype = 1
      allocate (lat%symbolic_atoms(1))
      call lat%symbolic_atoms(1)%restore_to_default()

      ham%lattice => lat
      ham%control => ctl
      allocate (ham%hubbard_u_pot(nb, nb, 1))
      ham%hubbard_u_pot = 0.0_rp
      ham%hubbard_u_potential_form = 'liechtenstein'
   end subroutine setup_fixture

   subroutine reset_fixture(u_s, j_s, u_p, j_p)
      real(rp), intent(in) :: u_s, j_s, u_p, j_p

      lat%symbolic_atoms(1)%potential%hubbard_u = 0.0_rp
      lat%symbolic_atoms(1)%potential%hubbard_j = 0.0_rp
      lat%symbolic_atoms(1)%potential%hubbard_u(1) = u_s
      lat%symbolic_atoms(1)%potential%hubbard_j(1) = j_s
      lat%symbolic_atoms(1)%potential%hubbard_u(2) = u_p
      lat%symbolic_atoms(1)%potential%hubbard_j(2) = j_p
      lat%symbolic_atoms(1)%potential%ldm = 0.0_rp
      ham%hubbard_u_pot = 0.0_rp
   end subroutine reset_fixture

   subroutine test_zero_u_j_limit()
      complex(rp) :: baseline(nb, nb), with_correction(nb, nb)
      integer :: i, j

      call reset_fixture(0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp)
      lat%symbolic_atoms(1)%potential%ldm = 0.23_rp
      do i = 1, nb
         do j = 1, nb
            baseline(i, j) = cmplx(0.013_rp*real(i + 2*j, rp), &
               0.002_rp*real(i - j, rp), rp)
         end do
      end do

      call ham%calculate_hubbard_u_potential_general()
      with_correction = baseline + ham%hubbard_u_pot(:, :, 1)
      call require(maxval(abs(ham%hubbard_u_pot(:, :, 1))) < tolerance, &
         'U=J=0 produces zero Hubbard correction')
      call require(maxval(abs(with_correction - baseline)) < tolerance, &
         'U=J=0 leaves the Hamiltonian contribution unchanged')
   end subroutine test_zero_u_j_limit

   subroutine test_simple_diagonal_occupation()
      real(rp) :: expected_up, expected_down

      ! For one s orbital and J=0, the production expression reduces to
      ! V_up = U*(1/2 - n_up) and V_down = U*(1/2 - n_down).
      call reset_fixture(2.0_rp, 0.0_rp, 0.0_rp, 0.0_rp)
      lat%symbolic_atoms(1)%potential%ldm(1, 1, 1, 1) = 0.25_rp
      lat%symbolic_atoms(1)%potential%ldm(1, 2, 1, 1) = 0.75_rp

      ! Occupation outside the correlated shell must not leak into the s block
      ! or cause the routine to assume a d-sized Hamiltonian.
      lat%symbolic_atoms(1)%potential%ldm(2, 1, 1, 1) = 0.60_rp
      lat%symbolic_atoms(1)%potential%ldm(2, 2, 2, 2) = 0.30_rp

      call ham%calculate_hubbard_u_potential_general()
      expected_up = 2.0_rp*(0.5_rp - 0.25_rp)
      expected_down = 2.0_rp*(0.5_rp - 0.75_rp)
      call require(abs(ham%hubbard_u_pot(1, 1, 1) - expected_up) < tolerance, &
         'diagonal s-shell spin-up correction')
      call require(abs(ham%hubbard_u_pot(1 + spin_off, 1 + spin_off, 1) - expected_down) < tolerance, &
         'diagonal s-shell spin-down correction')
      call require(maxval(abs(ham%hubbard_u_pot(2:spin_off, 2:spin_off, 1))) < tolerance, &
         'uncorrelated p orbitals remain zero')
      call require(maxval(abs(ham%hubbard_u_pot(1:spin_off, spin_off + 1:nb, 1))) < tolerance, &
         'ordinary collinear path has no spin-mixing Hubbard block')
   end subroutine test_simple_diagonal_occupation

   subroutine test_collinear_hermiticity()
      complex(rp) :: v_up(3, 3), v_down(3, 3)

      call reset_fixture(0.0_rp, 0.0_rp, 1.7_rp, 0.3_rp)
      call set_symmetric_p_occupation(1, [0.20_rp, 0.40_rp, 0.10_rp], &
         [0.07_rp, -0.03_rp, 0.05_rp])
      call set_symmetric_p_occupation(2, [0.30_rp, 0.15_rp, 0.25_rp], &
         [0.04_rp, -0.02_rp, 0.06_rp])

      call ham%calculate_hubbard_u_potential_general()
      v_up = ham%hubbard_u_pot(2:4, 2:4, 1)
      v_down = ham%hubbard_u_pot(2 + spin_off:4 + spin_off, 2 + spin_off:4 + spin_off, 1)
      call require(maxval(abs(v_up - conjg(transpose(v_up)))) < tolerance, &
         'spin-up Hubbard potential is Hermitian')
      call require(maxval(abs(v_down - conjg(transpose(v_down)))) < tolerance, &
         'spin-down Hubbard potential is Hermitian')
      call require(maxval(abs(ham%hubbard_u_pot(1:spin_off, spin_off + 1:nb, 1))) < tolerance, &
         'Hermitian collinear correction remains block diagonal in spin')
   end subroutine test_collinear_hermiticity

   subroutine test_complex_hermitian_occupation()
      complex(rp) :: v_up(3, 3), v_down(3, 3)

      call reset_fixture(0.0_rp, 0.0_rp, 1.7_rp, 0.3_rp)
      call set_complex_p_occupation(1, [0.20_rp, 0.40_rp, 0.10_rp], &
         [(0.07_rp, 0.04_rp), (-0.03_rp, 0.02_rp), (0.05_rp, -0.06_rp)])
      call set_complex_p_occupation(2, [0.30_rp, 0.15_rp, 0.25_rp], &
         [(0.04_rp, -0.03_rp), (-0.02_rp, 0.05_rp), (0.06_rp, 0.01_rp)])

      call ham%calculate_hubbard_u_potential_general()
      v_up = ham%hubbard_u_pot(2:4, 2:4, 1)
      v_down = ham%hubbard_u_pot(2 + spin_off:4 + spin_off, 2 + spin_off:4 + spin_off, 1)
      call require(maxval(abs(v_up - conjg(transpose(v_up)))) < tolerance, &
         'complex spin-up Hubbard potential is Hermitian')
      call require(maxval(abs(v_down - conjg(transpose(v_down)))) < tolerance, &
         'complex spin-down Hubbard potential is Hermitian')
      call require(maxval(abs(aimag(v_up))) > 1.0e-8_rp .and. maxval(abs(aimag(v_down))) > 1.0e-8_rp, &
         'complex occupations produce a retained imaginary Hubbard correction')
   end subroutine test_complex_hermitian_occupation

   subroutine test_complex_energy_derivative()
      complex(rp) :: nup(3, 3), ndn(3, 3), dup(3, 3), ddn(3, 3)
      complex(rp) :: vup(3, 3), vdn(3, 3), tval
      real(rp), parameter :: u = 1.7_rp, j = 0.3_rp, h = 1.0e-5_rp
      real(rp) :: analytic, numeric, eplus, eminus
      integer :: a, b

      call reset_fixture(0.0_rp, 0.0_rp, u, j)
      nup = (0.0_rp, 0.0_rp)
      ndn = (0.0_rp, 0.0_rp)
      dup = (0.0_rp, 0.0_rp)
      ddn = (0.0_rp, 0.0_rp)
      do a = 1, 3
         nup(a,a) = 0.23_rp + 0.035_rp*real(a, rp)
         ndn(a,a) = 0.16_rp + 0.025_rp*real(a, rp)
         dup(a,a) = 0.013_rp*real(a - 2, rp)
         ddn(a,a) = -0.009_rp*real(a - 1, rp)
         do b = a + 1, 3
            nup(a,b) = cmplx(0.012_rp*real(a + b, rp), 0.016_rp*real(a - b, rp), rp)
            nup(b,a) = conjg(nup(a,b))
            ndn(a,b) = cmplx(0.008_rp*real(a + b, rp), -0.011_rp*real(a - b, rp), rp)
            ndn(b,a) = conjg(ndn(a,b))
            dup(a,b) = cmplx(0.003_rp*real(a + b, rp), 0.007_rp*real(a - b, rp), rp)
            dup(b,a) = conjg(dup(a,b))
            ddn(a,b) = cmplx(-0.002_rp*real(a + b, rp), 0.005_rp*real(a - b, rp), rp)
            ddn(b,a) = conjg(ddn(a,b))
         end do
      end do
      lat%symbolic_atoms(1)%potential%ldm(2,1,1:3,1:3) = nup
      lat%symbolic_atoms(1)%potential%ldm(2,2,1:3,1:3) = ndn

      call ham%calculate_hubbard_u_potential_general()
      vup = ham%hubbard_u_pot(2:4,2:4,1)
      vdn = ham%hubbard_u_pot(2+spin_off:4+spin_off,2+spin_off:4+spin_off,1)
      analytic = 0.0_rp
      do a = 1, 3
         do b = 1, 3
            analytic = analytic + real(vup(a,b)*dup(b,a) + vdn(a,b)*ddn(b,a), rp)
         end do
      end do
      eplus = complex_p_energy(nup + h*dup, ndn + h*ddn, u, j)
      eminus = complex_p_energy(nup - h*dup, ndn - h*ddn, u, j)
      numeric = (eplus - eminus)/(2.0_rp*h)
      call require(abs(numeric - analytic) < 2.0e-8_rp, &
         'production complex-P potential equals Tr(V delta P) for nonzero J')
   end subroutine test_complex_energy_derivative

   function complex_p_energy(pup, pdn, u, jbar) result(energy)
      complex(rp), intent(in) :: pup(3,3), pdn(3,3)
      real(rp), intent(in) :: u, jbar
      real(rp) :: energy, direct, exchange, nup, ndn
      complex(rp) :: densities(3,3,2), tensor
      integer :: a, b, c, d, s, t

      densities(:,:,1) = pup
      densities(:,:,2) = pdn
      direct = 0.0_rp
      exchange = 0.0_rp
      do s = 1, 2
         do t = 1, 2
            do a = 1, 3
               do b = 1, 3
                  do c = 1, 3
                     do d = 1, 3
                        tensor = Coulomb_mat_spherical(1, a-2, b-2, c-2, d-2, u, 5.0_rp*jbar, 0.0_rp, 0.0_rp)
                        direct = direct + real(tensor*densities(b,a,s)*densities(d,c,t), rp)
                     end do
                  end do
               end do
            end do
         end do
         do a = 1, 3
            do b = 1, 3
               do c = 1, 3
                  do d = 1, 3
                     tensor = Coulomb_mat_spherical(1, a-2, b-2, c-2, d-2, u, 5.0_rp*jbar, 0.0_rp, 0.0_rp)
                     exchange = exchange + real(tensor*densities(d,a,s)*densities(b,c,s), rp)
                  end do
               end do
            end do
         end do
      end do
      nup = sum([(real(pup(a,a), rp), a=1,3)])
      ndn = sum([(real(pdn(a,a), rp), a=1,3)])
      energy = 0.5_rp*(direct - exchange) - 0.5_rp*u*(nup+ndn)*(nup+ndn-1.0_rp) &
         + 0.5_rp*jbar*(nup*(nup-1.0_rp) + ndn*(ndn-1.0_rp))
   end function complex_p_energy

   subroutine test_correlated_trace_and_common_state_fourier()
      type(reciprocal) :: recip
      complex(rp), allocatable :: rs_block(:, :), hk(:, :), blocks(:, :, :, :)
      integer :: i, j
      real(rp) :: trace_up, trace_down

      ! Freeze one production occupation matrix, then compare the exact same
      ! onsite block as inserted in the RS slot and Fourier-transformed at k.
      ! This is a common-state representation check, not two SCF fixed points.
      call reset_fixture(0.0_rp, 0.0_rp, 1.7_rp, 0.3_rp)
      call set_symmetric_p_occupation(1, [0.20_rp, 0.40_rp, 0.10_rp], &
         [0.07_rp, -0.03_rp, 0.05_rp])
      call set_symmetric_p_occupation(2, [0.30_rp, 0.15_rp, 0.25_rp], &
         [0.04_rp, -0.02_rp, 0.06_rp])
      call ham%calculate_hubbard_u_potential_general()

      trace_up = 0.0_rp
      trace_down = 0.0_rp
      do i = 1, 3
         trace_up = trace_up + real(lat%symbolic_atoms(1)%potential%ldm(2, 1, i, i), rp)
         trace_down = trace_down + real(lat%symbolic_atoms(1)%potential%ldm(2, 2, i, i), rp)
      end do
      call require(abs(trace_up - 0.70_rp) < tolerance, 'correlated p-shell spin-up trace')
      call require(abs(trace_down - 0.70_rp) < tolerance, 'correlated p-shell spin-down trace')
      call require(abs(trace_up + trace_down - 1.40_rp) < tolerance, 'correlated p-shell total trace')

      lat%nrec = 1
      lat%nn_max = 1
      lat%kk = 1
      allocate(lat%ib(1), lat%atlist(1), lat%iz(1), lat%nn(1, 1))
      lat%ib = 1
      lat%atlist = 1
      lat%iz = 1
      lat%nn = 1
      recip%lattice => lat
      allocate(recip%ham_vec_type_direct(3, 1, 1))
      recip%ham_vec_type_direct = 0.0_rp
      allocate(rs_block(nb, nb), hk(nb, nb), blocks(nb, nb, 1, 1))
      rs_block = ham%hubbard_u_pot(:, :, 1)
      blocks(:, :, 1, 1) = rs_block
      call recip%fourier_transform_array(blocks, [0.173_rp, 0.071_rp, 0.019_rp], hk)
      call require(maxval(abs(hk - rs_block)) < tolerance, &
         'frozen onsite Hubbard block is identical in RS and reciprocal routes')
   end subroutine test_correlated_trace_and_common_state_fourier

   subroutine set_symmetric_p_occupation(ispin, diagonal, off_diagonal)
      integer, intent(in) :: ispin
      real(rp), intent(in) :: diagonal(3), off_diagonal(3)

      lat%symbolic_atoms(1)%potential%ldm(2, ispin, 1, 1) = diagonal(1)
      lat%symbolic_atoms(1)%potential%ldm(2, ispin, 2, 2) = diagonal(2)
      lat%symbolic_atoms(1)%potential%ldm(2, ispin, 3, 3) = diagonal(3)
      lat%symbolic_atoms(1)%potential%ldm(2, ispin, 1, 2) = off_diagonal(1)
      lat%symbolic_atoms(1)%potential%ldm(2, ispin, 2, 1) = off_diagonal(1)
      lat%symbolic_atoms(1)%potential%ldm(2, ispin, 1, 3) = off_diagonal(2)
      lat%symbolic_atoms(1)%potential%ldm(2, ispin, 3, 1) = off_diagonal(2)
      lat%symbolic_atoms(1)%potential%ldm(2, ispin, 2, 3) = off_diagonal(3)
      lat%symbolic_atoms(1)%potential%ldm(2, ispin, 3, 2) = off_diagonal(3)
   end subroutine set_symmetric_p_occupation

   subroutine set_complex_p_occupation(ispin, diagonal, off_diagonal)
      integer, intent(in) :: ispin
      real(rp), intent(in) :: diagonal(3)
      complex(rp), intent(in) :: off_diagonal(3)

      lat%symbolic_atoms(1)%potential%ldm(2, ispin, 1, 1) = diagonal(1)
      lat%symbolic_atoms(1)%potential%ldm(2, ispin, 2, 2) = diagonal(2)
      lat%symbolic_atoms(1)%potential%ldm(2, ispin, 3, 3) = diagonal(3)
      lat%symbolic_atoms(1)%potential%ldm(2, ispin, 1, 2) = off_diagonal(1)
      lat%symbolic_atoms(1)%potential%ldm(2, ispin, 2, 1) = conjg(off_diagonal(1))
      lat%symbolic_atoms(1)%potential%ldm(2, ispin, 1, 3) = off_diagonal(2)
      lat%symbolic_atoms(1)%potential%ldm(2, ispin, 3, 1) = conjg(off_diagonal(2))
      lat%symbolic_atoms(1)%potential%ldm(2, ispin, 2, 3) = off_diagonal(3)
      lat%symbolic_atoms(1)%potential%ldm(2, ispin, 3, 2) = conjg(off_diagonal(3))
   end subroutine set_complex_p_occupation

   subroutine require(condition, description)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: description

      if (.not. condition) then
         write (*, '(a,a)') 'FAIL: ', trim(description)
         failed = .true.
      end if
   end subroutine require

end program test_lda_u_hamiltonian
