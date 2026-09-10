!------------------------------------------------------------------------------
! LR-BASIS-00 production Hamiltonian and state-wise density oracle
!------------------------------------------------------------------------------
program test_lr_basis_augmentation
   use precision_mod, only: rp
   use basis_mod, only: basis_init, nb, norb
   use logger_mod, only: g_logger
   use math_mod, only: init_math_operators, pi
   use self_mod, only: legacy_radial_fixture, legacy_newrho_fixture
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, lmto_orbital_l
   use reciprocal_mod, only: reciprocal
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use charge_mod, only: charge
   use control_mod, only: control
   implicit none

   integer, parameter :: nr = 101
   real(rp), parameter :: a = 0.03_rp, b = 0.10_rp, z = 1.0_rp
   real(rp), parameter :: matrix_tol = 2.0e-12_rp
   real(rp), parameter :: density_tol = 2.0e-11_rp
   real(rp) :: rofi(nr), potential(nr), potential_spin(nr, 2), energy, work_energy
   real(rp) :: radial_energies(3, 2)
   real(rp), allocatable :: g(:, :), gp(:, :), gpp(:, :)
   real(rp), allocatable :: gpack(:), gdotpack(:), gddotpack(:)
   real(rp), allocatable :: evals(:, :), moments(:, :, :)
   real(rp), allocatable :: density_state(:, :, :), density_sum(:, :, :), density_api_sum(:, :, :), density_poly(:, :, :)
   real(rp), allocatable :: density_newrho(:, :), ql_newrho(:, :, :)
   complex(rp), allocatable :: hk(:, :), h(:, :), eeok(:, :), hgamma(:, :), expected(:, :)
   complex(rp), allocatable :: coeff(:, :)
   complex(rp), allocatable :: evecs(:, :, :)
   type(lmto_radial_basis) :: radial
   type(lmto_radial_basis) :: radial_spdf
   type(reciprocal) :: recip
   type(hamiltonian), target :: ham
   type(lattice), target :: lat
   type(charge), target :: chg
   type(control), target :: ctl
   real(rp) :: kpoint(3), matrix_error, residual, density_error, density_api_error, weight
   real(rp) :: delta, amplitude, cross0, cross1, cross2
   integer :: ir, l, ispin, i, j, ib

   call basis_init(2)
   call init_math_operators()
   call g_logger%init()

   rofi(1) = 0.0_rp
   do ir = 2, nr
      rofi(ir) = b*(exp(a*real(ir - 1, rp)) - 1.0_rp)
   end do
   potential = 0.0_rp
   potential_spin = 0.0_rp
   call radial%initialize(nr, 2, 2)
   call radial_spdf%initialize(nr, 3, 1)
   if (.not. radial%supported('ham_only', 'second', .true., .true., .false., .false.)) then
      error stop 'LR-BASIS supported baseline was rejected'
   end if
   if (radial%supported('ham_only', 'second', .true., .true., .true., .false.) .or. &
       radial%supported('ham_only', 'second', .true., .false., .false., .false.) .or. &
       radial%supported('generalized_overlap_proxy', 'second', .true., .true., .false., .false.) .or. &
       radial%supported('ham_only', 'second', .true., .true., .false., .true.) .or. &
       radial_spdf%supported('ham_only', 'second', .true., .true., .false., .false.)) then
      error stop 'LR-BASIS unsupported mode was not rejected'
   end if
   do ispin = 1, 2
      do l = 0, 2
         call legacy_radial_fixture(z, l, a, b, rofi, potential, energy, g, gp, gpp, &
            tan(pi*(0.5_rp - (real(l + 1, rp) + 0.25_rp))) + 1.0_rp)
         gpack = reshape(g, [2*nr])
         gdotpack = reshape(gp, [2*nr])
         gddotpack = reshape(gpp, [2*nr])
         work_energy = energy
         radial_energies(l + 1, ispin) = energy
         call radial%capture_channel(l, ispin, energy, rofi, potential, a, b, z, gpack, gdotpack, gddotpack, work_energy)
      end do
   end do

   call setup_production_fixture(recip, ham, lat, chg, ctl, radial)
   allocate(hk(nb, nb), h(nb, nb), eeok(nb, nb), hgamma(nb, nb), expected(nb, nb))
   kpoint = [0.137_rp, -0.091_rp, 0.053_rp]
   call recip%build_hamiltonian_at_kpoint(kpoint, hk)
   call recip%fourier_transform_array(ham%ee, kpoint, h)
   call recip%fourier_transform_array(ham%eeo, kpoint, eeok)
   hgamma = h - matmul(eeok, h)
   expected = hgamma + ham%enim(:, :, 1)
   matrix_error = maxval(abs(hk - expected))
   write (*, '(a,es16.8)') 'production H(k) - [h-h*Obar*h+Enu] = ', matrix_error
   if (matrix_error > matrix_tol) error stop 'LR-BASIS HOH production assembly identity failed'

   call recip%calculate_eigenpairs_at_kpoints(reshape(kpoint, [3, 1]), evals, evecs)
   residual = 0.0_rp
   do ib = 1, nb
      residual = max(residual, maxval(abs(matmul(hgamma, evecs(:, ib, 1)) + &
         matmul(ham%enim(:, :, 1), evecs(:, ib, 1)) - evals(ib, 1)*evecs(:, ib, 1))))
   end do
   write (*, '(a,es16.8)') 'production eigenstate baseline residual = ', residual
   if (residual > matrix_tol) error stop 'LR-BASIS eigenstate baseline identity failed'

   allocate(moments(3, 3, 2), density_state(nr, 3, 2), density_sum(nr, 3, 2), &
      density_api_sum(nr, 3, 2), density_poly(nr, 3, 2))
   allocate(density_newrho(nr, 2), ql_newrho(3, 3, 2))
   allocate(coeff(norb, 2))
   moments = 0.0_rp
   density_sum = 0.0_rp
   density_api_sum = 0.0_rp
   do ib = 1, nb
      coeff = cmplx(0.0_rp, 0.0_rp, rp)
      coeff(:, 1) = evecs(1:norb, ib, 1)
      coeff(:, 2) = evecs(norb + 1:nb, ib, 1)
      weight = 0.23_rp + 0.011_rp*real(ib, rp)
      call radial%state_density(evals(ib, 1), coeff, density_state)
      density_api_sum = density_api_sum + weight*density_state
      do ispin = 1, 2
         do i = 1, norb
            l = lmto_orbital_l(i)
            amplitude = real(coeff(i, ispin)*conjg(coeff(i, ispin)), rp)
            delta = evals(ib, 1) - radial%enu_work(l + 1, ispin)
            do ir = 1, nr
               cross0 = radial%gfac(ir, l + 1, ispin)*radial%phi_large(ir, l + 1, ispin)**2 + &
                        radial%phi_small(ir, l + 1, ispin)**2
               cross1 = radial%gfac(ir, l + 1, ispin)*radial%phi_large(ir, l + 1, ispin)* &
                        radial%phidot_large(ir, l + 1, ispin) + radial%phi_small(ir, l + 1, ispin)* &
                        radial%phidot_small(ir, l + 1, ispin)
               cross2 = radial%gfac(ir, l + 1, ispin)*(radial%phidot_large(ir, l + 1, ispin)**2 + &
                        radial%phi_large(ir, l + 1, ispin)*radial%phiddot_large(ir, l + 1, ispin)) + &
                        radial%phidot_small(ir, l + 1, ispin)**2 + radial%phi_small(ir, l + 1, ispin)* &
                        radial%phiddot_small(ir, l + 1, ispin)
               density_sum(ir, l + 1, ispin) = density_sum(ir, l + 1, ispin) + weight*amplitude*(cross0 + &
                  2.0_rp*delta*cross1 + delta**2*cross2)
            end do
            moments(1, l + 1, ispin) = moments(1, l + 1, ispin) + &
               weight*amplitude
            moments(2, l + 1, ispin) = moments(2, l + 1, ispin) + &
               weight*evals(ib, 1)*amplitude
            moments(3, l + 1, ispin) = moments(3, l + 1, ispin) + &
               weight*evals(ib, 1)**2*amplitude
         end do
      end do
   end do
   call radial%density_from_moments(moments, density_poly)
   density_error = maxval(abs(density_sum - density_poly))
   density_api_error = maxval(abs(density_sum - density_api_sum))
   write (*, '(a,es16.8)') 'state-wise API density - independent state density = ', density_api_error
   write (*, '(a,es16.8)') 'state-wise radial density - NEWRHO polynomial = ', density_error
   if (density_api_error > density_tol) error stop 'LR-BASIS state-density API failed'
   if (density_error > density_tol) error stop 'LR-BASIS state-wise density identity failed'

   ql_newrho = 0.0_rp
   do ispin = 1, 2
      do l = 0, 2
         ql_newrho(1, l + 1, ispin) = moments(1, l + 1, ispin)
         ql_newrho(2, l + 1, ispin) = moments(2, l + 1, ispin) - radial%enu_work(l + 1, ispin)*moments(1, l + 1, ispin)
         ql_newrho(3, l + 1, ispin) = moments(3, l + 1, ispin) - &
            2.0_rp*radial%enu_work(l + 1, ispin)*moments(2, l + 1, ispin) + &
            radial%enu_work(l + 1, ispin)**2*moments(1, l + 1, ispin)
      end do
   end do
   call legacy_newrho_fixture(2, z, a, b, rofi, potential_spin, radial_energies, ql_newrho, density_newrho)
   density_error = maxval(abs(density_newrho(2:nr, :) - sum(density_poly(2:nr, :, :), dim=2)))
   write (*, '(a,es16.8)') 'state-wise radial density - live NEWRHO = ', density_error
   if (density_error > density_tol) error stop 'LR-BASIS live NEWRHO density comparison failed'

   write (*, '(a)') 'UnitLrBasisAugmentation: PASS (production H, eigenstates, and radial density)'

contains

   subroutine setup_production_fixture(recip, ham, lat, chg, ctl, radial)
      type(reciprocal), intent(out) :: recip
      type(hamiltonian), target, intent(out) :: ham
      type(lattice), target, intent(out) :: lat
      type(charge), target, intent(out) :: chg
      type(control), target, intent(out) :: ctl
      type(lmto_radial_basis), intent(in) :: radial
      integer :: i, j, l, spin_i

      call recip%restore_to_default()
      call ctl%restore_to_default()
      lat%nrec = 1
      lat%ntype = 1
      lat%nn_max = 1
      lat%kk = 1
      lat%nmax = 1
      lat%alat = 1.0_rp
      lat%r2 = 1.0_rp
      lat%a = 0.0_rp
      lat%a(1, 1) = 1.0_rp; lat%a(2, 2) = 1.0_rp; lat%a(3, 3) = 1.0_rp
      lat%a_cart_inv = lat%a
      lat%a_cart_inv_ready = .true.
      allocate(lat%ib(1), lat%atlist(1), lat%iz(1), lat%num(1), lat%nn(1, 1), lat%sbar(norb, norb, 1, 1), &
         lat%symbolic_atoms(1))
      lat%ib = 1; lat%atlist = 1; lat%iz = 1; lat%num = 1; lat%nn = 1; lat%sbar = cmplx(0.0_rp, 0.0_rp, rp)
      call lat%symbolic_atoms(1)%restore_to_default()

      chg%lattice => lat
      chg%symbolic_atom => lat%symbolic_atoms
      ham%charge => chg
      ham%lattice => lat
      ham%control => ctl
      ham%hoh = .true.
      ham%ccor_2c = .false.
      ham%local_axis = .false.
      ham%magnetic_representation = 'collinear'
      ham%operator_generation = 1

      allocate(ham%ee(nb, nb, 1, 1), ham%eeo(nb, nb, 1, 1), ham%eeoee(nb, nb, 1, 1), &
         ham%lsham(nb, nb, 1), ham%obarm(nb, nb, 1), ham%enim(nb, nb, 1))
      ham%ee = cmplx(0.0_rp, 0.0_rp, rp)
      ham%obarm = cmplx(0.0_rp, 0.0_rp, rp)
      ham%enim = cmplx(0.0_rp, 0.0_rp, rp)
      ham%lsham = cmplx(0.0_rp, 0.0_rp, rp)
      ham%eeoee = cmplx(0.0_rp, 0.0_rp, rp)
      do spin_i = 1, 2
         do i = 1, norb
            ham%obarm(i + (spin_i - 1)*norb, i + (spin_i - 1)*norb, 1) = &
               cmplx(0.12_rp + 0.006_rp*real(i, rp) + 0.02_rp*real(spin_i - 1, rp), 0.0_rp, rp)
            l = lmto_orbital_l(i)
            ham%enim(i + (spin_i - 1)*norb, i + (spin_i - 1)*norb, 1) = &
               cmplx(radial%enu_work(l + 1, spin_i), 0.0_rp, rp)
         end do
         do i = 1, norb
            ham%ee(i + (spin_i - 1)*norb, i + (spin_i - 1)*norb, 1, 1) = &
               cmplx(0.03_rp*real(i, rp) + 0.004_rp*real(spin_i, rp), 0.0_rp, rp)
            do j = i + 1, norb
               ham%ee(i + (spin_i - 1)*norb, j + (spin_i - 1)*norb, 1, 1) = &
                  cmplx(0.0007_rp*real(i + j, rp), 0.0003_rp*real(j - i, rp), rp)
               ham%ee(j + (spin_i - 1)*norb, i + (spin_i - 1)*norb, 1, 1) = &
                  conjg(ham%ee(i + (spin_i - 1)*norb, j + (spin_i - 1)*norb, 1, 1))
            end do
         end do
      end do
      ham%eeo(:, :, 1, 1) = matmul(ham%ee(:, :, 1, 1), ham%obarm(:, :, 1))

      recip%hamiltonian => ham
      recip%lattice => lat
      recip%control => ctl
      recip%reciprocal_mode = 'ham_only'
      recip%kspace_ham_order = 'second'
      recip%max_orbs = nb
      allocate(recip%ham_vec_type(3, 1, 1), recip%ham_vec_type_direct(3, 1, 1))
      recip%ham_vec_type = 0.0_rp
      recip%ham_vec_type_direct = 0.0_rp
   end subroutine setup_production_fixture

end program test_lr_basis_augmentation
