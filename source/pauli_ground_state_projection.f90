!------------------------------------------------------------------------------
! RS-LMTO-ASA
!
! MODULE: pauli_ground_state_projection
!
!> @brief LR-02N scalar-relativistic to Pauli ground-state projection.
!>
!> This module is deliberately a numerical diagnostic.  It does not build a
!> response kernel.  It contracts the occupied reciprocal eigenvectors with
!> the accepted large-component radial augmentation and compares that result
!> with the exact LR-01 radial snapshot.
!------------------------------------------------------------------------------
module pauli_ground_state_projection_mod

   use precision_mod, only: rp
   use reciprocal_mod, only: reciprocal
   use symbolic_atom_mod, only: symbolic_atom
   use radial_ground_state_mod, only: radial_ground_state, RADIAL_PI, radial_simpson_weight
   use lmto_radial_augmentation_mod, only: lmto_orbital_l
   use mpi_mod, only: rank, ierr
#ifdef USE_MPI
   use mpi
#endif
   implicit none
   private

   real(rp), parameter :: kB_RY_PER_K = 6.3336814e-6_rp
   real(rp), parameter :: MAGNETIC_REGION_FRACTION = 0.90_rp
   real(rp), parameter :: OCCUPATION_CUTOFF = 1.0e-14_rp

   public :: pauli_ground_state_projection
   public :: compute_accepted_pauli_magnetization

contains

   !> Return the accepted Pauli number-spin magnetization used by LR-03.
   !>
   !> This is the same occupied reciprocal eigensystem plus accepted POTPAR
   !> large-component/core projection used by the LR-02N diagnostic, but it is
   !> exposed as data for the static TDVK-06 interaction services.  No EF,
   !> occupation, radial basis, or response-space object is rebuilt here.
   subroutine compute_accepted_pauli_magnetization(reciprocal_obj, atoms, nbulk, magnetization)
      type(reciprocal), intent(in) :: reciprocal_obj
      type(symbolic_atom), intent(in) :: atoms(:)
      integer, intent(in) :: nbulk
      real(rp), allocatable, intent(out) :: magnetization(:, :)

      integer :: isite, atom_index, nsite, nmat, norb_site, nr, lmax, ir, ik, ik_global, ib, iorb, ispin, l
      real(rp), allocatable :: valence_weighted(:, :), pauli_weighted(:, :)
      real(rp) :: wk, occ, energy, kT, delta_energy, amplitude, radial_amplitude
      complex(rp) :: coefficient

      if (.not. allocated(reciprocal_obj%eigenvalues) .or. .not. allocated(reciprocal_obj%eigenvectors) .or. &
          .not. allocated(reciprocal_obj%k_weights)) then
         error stop 'TDVK-06 Pauli magnetization: accepted reciprocal eigensystem is incomplete'
      end if
      if (.not. associated(reciprocal_obj%lattice)) then
         error stop 'TDVK-06 Pauli magnetization: reciprocal lattice association is missing'
      end if
      nsite = reciprocal_obj%lattice%nrec
      nmat = size(reciprocal_obj%eigenvectors, 1)
      if (nsite < 1 .or. mod(nmat, nsite) /= 0 .or. mod(nmat/nsite, 2) /= 0) then
         error stop 'TDVK-06 Pauli magnetization: reciprocal basis is not site-major two-spin'
      end if
      norb_site = (nmat/nsite)/2
      kT = reciprocal_obj%temperature*kB_RY_PER_K
      nr = 0
      do isite = 1, nsite
         atom_index = nbulk + isite
         if (atom_index < 1 .or. atom_index > size(atoms)) then
            error stop 'TDVK-06 Pauli magnetization: site-to-atom mapping is outside the accepted state'
         end if
         if (.not. atoms(atom_index)%radial_ground_state%valid .or. &
             .not. atoms(atom_index)%radial_ground_state%accepted .or. &
             .not. atoms(atom_index)%radial_ground_state%pauli_basis_valid .or. &
             .not. atoms(atom_index)%radial_ground_state%core_density_valid) then
            error stop 'TDVK-06 Pauli magnetization: accepted LR-01 Pauli/core projection is incomplete'
         end if
         if (nr == 0) nr = size(atoms(atom_index)%radial_ground_state%r)
         if (size(atoms(atom_index)%radial_ground_state%r) /= nr) then
            error stop 'TDVK-06 Pauli magnetization: response sites do not share one radial mesh'
         end if
      end do
      allocate(magnetization(nsite, nr))
      magnetization = 0.0_rp

      do isite = 1, nsite
         atom_index = nbulk + isite
         lmax = atoms(atom_index)%radial_ground_state%pauli_lmax
         allocate(valence_weighted(nr, 2), pauli_weighted(nr, 2))
         valence_weighted = 0.0_rp
         do ik = 1, size(reciprocal_obj%eigenvalues, 2)
            ik_global = ik
            if (allocated(reciprocal_obj%k_l2g_map)) ik_global = reciprocal_obj%k_l2g_map(ik)
            wk = reciprocal_obj%k_weights(ik_global)
            do ib = 1, size(reciprocal_obj%eigenvalues, 1)
               energy = reciprocal_obj%eigenvalues(ib, ik)
               occ = fermi_dirac(energy, reciprocal_obj%fermi_level, kT)
               if (occ <= OCCUPATION_CUTOFF) cycle
               do ispin = 1, 2
                  do iorb = 1, norb_site
                     l = lmto_orbital_l(iorb)
                     if (l < 0 .or. l > lmax) cycle
                     coefficient = reciprocal_obj%eigenvectors((isite - 1)*2*norb_site + &
                        (ispin - 1)*norb_site + iorb, ib, ik)
                     amplitude = real(coefficient*conjg(coefficient), rp)
                     delta_energy = energy - atoms(atom_index)%radial_ground_state%pauli_enu(l + 1, ispin)
                     do ir = 1, nr
                        radial_amplitude = atoms(atom_index)%radial_ground_state%pauli_large(ir, l + 1, ispin) + &
                           delta_energy*atoms(atom_index)%radial_ground_state%pauli_large_dot(ir, l + 1, ispin)
                        valence_weighted(ir, ispin) = valence_weighted(ir, ispin) + &
                           wk*occ*amplitude*radial_amplitude**2
                     end do
                  end do
               end do
            end do
         end do
         pauli_weighted(:, 1) = valence_weighted(:, 1) + &
            atoms(atom_index)%radial_ground_state%core_pauli_weighted_up
         pauli_weighted(:, 2) = valence_weighted(:, 2) + &
            atoms(atom_index)%radial_ground_state%core_pauli_weighted_down
         do ir = 1, nr
            if (ir == 1) then
               magnetization(isite, ir) = origin_density(pauli_weighted(:, 1), atoms(atom_index)%radial_ground_state%r) - &
                  origin_density(pauli_weighted(:, 2), atoms(atom_index)%radial_ground_state%r)
            else
               magnetization(isite, ir) = pauli_weighted(ir, 1)/(4.0_rp*RADIAL_PI* &
                  atoms(atom_index)%radial_ground_state%r(ir)**2) - &
                  pauli_weighted(ir, 2)/(4.0_rp*RADIAL_PI*atoms(atom_index)%radial_ground_state%r(ir)**2)
            end if
         end do
         deallocate(valence_weighted, pauli_weighted)
      end do
   end subroutine compute_accepted_pauli_magnetization

   !> Evaluate and write the LR-02N diagnostic for every reciprocal site.
   !>
   !> The reciprocal object must already contain the final accepted-potential
   !> eigensystem.  Its k weights, Fermi level, temperature, mode, and mesh are
   !> consumed as-is; no independent occupation reconstruction is permitted.
   subroutine pauli_ground_state_projection(reciprocal_obj, atoms, nbulk, output_prefix)
      type(reciprocal), intent(in) :: reciprocal_obj
      type(symbolic_atom), intent(in) :: atoms(:)
      integer, intent(in) :: nbulk
      character(len=*), intent(in) :: output_prefix

      integer :: isite, atom_index, nsite, nmat, norb_site, nr, lmax, ir, ik, ik_global, ib, iorb, ispin, l
      integer :: unit, ios, region_last
      real(rp), allocatable :: valence_weighted(:, :), pauli_weighted(:, :)
      real(rp), allocatable :: n_sr(:, :), n_pauli(:, :), delta_n(:, :)
      real(rp), allocatable :: m_sr(:), m_pauli(:), delta_m(:), charge_sr(:), charge_pauli(:), delta_charge(:)
      real(rp), allocatable :: abs_m_weight(:)
      real(rp) :: wk, occ, energy, farg, kT, delta_energy, amplitude, radial_amplitude
      real(rp) :: integrated_sr(2), integrated_pauli(2), delta_number(2), integrated_m_sr, integrated_m_pauli, delta_moment
      real(rp) :: norm_m_sr, norm_delta_m, norm_charge_sr, norm_delta_charge
      real(rp) :: region_m_sr, region_delta_m, region_charge_sr, region_delta_charge
      real(rp) :: occupied_weight, k_weight_sum, region_radius, coefficient_weight(2), pauli_valence(2)
      complex(rp) :: coefficient
      character(len=512) :: data_file

      if (.not. allocated(reciprocal_obj%eigenvalues) .or. .not. allocated(reciprocal_obj%eigenvectors)) then
         error stop 'LR-02N requires occupied reciprocal eigenvectors from the accepted eigensystem'
      end if
      if (.not. allocated(reciprocal_obj%k_weights)) then
         error stop 'LR-02N requires the reciprocal k-point weights'
      end if
      if (.not. associated(reciprocal_obj%lattice)) then
         error stop 'LR-02N reciprocal object has no lattice'
      end if

      nsite = reciprocal_obj%lattice%nrec
      nmat = size(reciprocal_obj%eigenvectors, 1)
      if (nsite < 1 .or. mod(nmat, nsite) /= 0 .or. mod(nmat / nsite, 2) /= 0) then
         error stop 'LR-02N requires a site-major, two-spin reciprocal basis'
      end if
      norb_site = (nmat / nsite) / 2
      kT = reciprocal_obj%temperature*kB_RY_PER_K
      k_weight_sum = sum(reciprocal_obj%k_weights)
      occupied_weight = 0.0_rp
      do ik = 1, size(reciprocal_obj%eigenvalues, 2)
         ik_global = ik
         if (allocated(reciprocal_obj%k_l2g_map)) ik_global = reciprocal_obj%k_l2g_map(ik)
         wk = reciprocal_obj%k_weights(ik_global)
         do ib = 1, size(reciprocal_obj%eigenvalues, 1)
            energy = reciprocal_obj%eigenvalues(ib, ik)
            occupied_weight = occupied_weight + wk*fermi_dirac(energy, reciprocal_obj%fermi_level, kT)
         end do
      end do
#ifdef USE_MPI
      if (reciprocal_obj%k_mesh_distributed_active) then
         call MPI_ALLREDUCE(MPI_IN_PLACE, occupied_weight, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      end if
#endif

      do isite = 1, nsite
         atom_index = nbulk + isite
         if (atom_index < 1 .or. atom_index > size(atoms)) then
            error stop 'LR-02N site-to-atom mapping is outside the symbolic-atom array'
         end if
         if (.not. atoms(atom_index)%radial_ground_state%valid .or. &
             .not. atoms(atom_index)%radial_ground_state%accepted) then
            error stop 'LR-02N requires an accepted LR-01 radial snapshot'
         end if
         if (.not. atoms(atom_index)%radial_ground_state%pauli_basis_valid) then
            error stop 'LR-02N requires the accepted POTPAR large-component augmentation'
         end if
         if (.not. atoms(atom_index)%radial_ground_state%core_density_valid) then
            error stop 'LR-02N requires the accepted frozen-core decomposition'
         end if

         nr = size(atoms(atom_index)%radial_ground_state%r)
         lmax = atoms(atom_index)%radial_ground_state%pauli_lmax
         allocate(valence_weighted(nr, 2), pauli_weighted(nr, 2))
         valence_weighted = 0.0_rp
         coefficient_weight = 0.0_rp

         do ik = 1, size(reciprocal_obj%eigenvalues, 2)
            ik_global = ik
            if (allocated(reciprocal_obj%k_l2g_map)) ik_global = reciprocal_obj%k_l2g_map(ik)
            wk = reciprocal_obj%k_weights(ik_global)
            do ib = 1, size(reciprocal_obj%eigenvalues, 1)
               energy = reciprocal_obj%eigenvalues(ib, ik)
               occ = fermi_dirac(energy, reciprocal_obj%fermi_level, kT)
               if (occ <= OCCUPATION_CUTOFF) cycle
               do ispin = 1, 2
                  do iorb = 1, norb_site
                     l = lmto_orbital_l(iorb)
                     if (l < 0 .or. l > lmax) cycle
                     coefficient = reciprocal_obj%eigenvectors((isite - 1)*2*norb_site + &
                        (ispin - 1)*norb_site + iorb, ib, ik)
                     amplitude = real(coefficient*conjg(coefficient), rp)
                     coefficient_weight(ispin) = coefficient_weight(ispin) + wk*occ*amplitude
                     delta_energy = energy - atoms(atom_index)%radial_ground_state%pauli_enu(l + 1, ispin)
                     do ir = 1, nr
                        radial_amplitude = atoms(atom_index)%radial_ground_state%pauli_large(ir, l + 1, ispin) + &
                           delta_energy*atoms(atom_index)%radial_ground_state%pauli_large_dot(ir, l + 1, ispin)
                        valence_weighted(ir, ispin) = valence_weighted(ir, ispin) + wk*occ*amplitude*radial_amplitude**2
                     end do
                  end do
               end do
            end do
         end do

#ifdef USE_MPI
         if (reciprocal_obj%k_mesh_distributed_active) then
            call MPI_ALLREDUCE(MPI_IN_PLACE, valence_weighted, size(valence_weighted), MPI_DOUBLE_PRECISION, MPI_SUM, &
                               MPI_COMM_WORLD, ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE, coefficient_weight, 2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         end if
#endif

         ! The reciprocal LMTO eigensystem is valence-only.  The frozen core
         ! is projected through the same accepted large-component rule, so the
         ! reported total compares all occupied charge without hiding a core
         ! electron-count offset in the LR-02N discrepancy.
         pauli_weighted(:, 1) = valence_weighted(:, 1) + &
            atoms(atom_index)%radial_ground_state%core_pauli_weighted_up
         pauli_weighted(:, 2) = valence_weighted(:, 2) + &
            atoms(atom_index)%radial_ground_state%core_pauli_weighted_down

         allocate(n_sr(nr, 2), n_pauli(nr, 2), delta_n(nr, 2), m_sr(nr), m_pauli(nr), delta_m(nr), &
                  charge_sr(nr), charge_pauli(nr), delta_charge(nr), abs_m_weight(nr))
         do ir = 1, nr
            n_sr(ir, 1) = atoms(atom_index)%radial_ground_state%n_up(ir)
            n_sr(ir, 2) = atoms(atom_index)%radial_ground_state%n_down(ir)
            if (ir == 1) then
               n_pauli(ir, 1) = origin_density(pauli_weighted(:, 1), atoms(atom_index)%radial_ground_state%r)
               n_pauli(ir, 2) = origin_density(pauli_weighted(:, 2), atoms(atom_index)%radial_ground_state%r)
            else
               n_pauli(ir, 1) = pauli_weighted(ir, 1)/(4.0_rp*RADIAL_PI*atoms(atom_index)%radial_ground_state%r(ir)**2)
               n_pauli(ir, 2) = pauli_weighted(ir, 2)/(4.0_rp*RADIAL_PI*atoms(atom_index)%radial_ground_state%r(ir)**2)
            end if
            delta_n(ir, :) = n_sr(ir, :) - n_pauli(ir, :)
            m_sr(ir) = n_sr(ir, 1) - n_sr(ir, 2)
            m_pauli(ir) = n_pauli(ir, 1) - n_pauli(ir, 2)
            delta_m(ir) = m_sr(ir) - m_pauli(ir)
            charge_sr(ir) = n_sr(ir, 1) + n_sr(ir, 2)
            charge_pauli(ir) = n_pauli(ir, 1) + n_pauli(ir, 2)
            delta_charge(ir) = charge_sr(ir) - charge_pauli(ir)
            abs_m_weight(ir) = volume_weight(atoms(atom_index)%radial_ground_state, ir)*abs(m_sr(ir))
         end do

         integrated_sr(1) = atoms(atom_index)%radial_ground_state%charge_integral(1)
         integrated_sr(2) = atoms(atom_index)%radial_ground_state%charge_integral(2)
         integrated_pauli(1) = radial_integral(atoms(atom_index)%radial_ground_state, pauli_weighted(:, 1))
         integrated_pauli(2) = radial_integral(atoms(atom_index)%radial_ground_state, pauli_weighted(:, 2))
         pauli_valence(1) = radial_integral(atoms(atom_index)%radial_ground_state, valence_weighted(:, 1))
         pauli_valence(2) = radial_integral(atoms(atom_index)%radial_ground_state, valence_weighted(:, 2))
         delta_number = integrated_sr - integrated_pauli
         integrated_m_sr = integrated_sr(1) - integrated_sr(2)
         integrated_m_pauli = integrated_pauli(1) - integrated_pauli(2)
         delta_moment = integrated_m_sr - integrated_m_pauli

         norm_m_sr = sqrt(volume_integral(atoms(atom_index)%radial_ground_state, m_sr**2))
         norm_delta_m = sqrt(volume_integral(atoms(atom_index)%radial_ground_state, delta_m**2))
         norm_charge_sr = sqrt(volume_integral(atoms(atom_index)%radial_ground_state, charge_sr**2))
         norm_delta_charge = sqrt(volume_integral(atoms(atom_index)%radial_ground_state, delta_charge**2))
         region_last = magnetic_region_endpoint(atoms(atom_index)%radial_ground_state, abs_m_weight, &
                                                MAGNETIC_REGION_FRACTION)
         region_radius = atoms(atom_index)%radial_ground_state%r(region_last)
         region_m_sr = sqrt(volume_integral_prefix(atoms(atom_index)%radial_ground_state, m_sr**2, region_last))
         region_delta_m = sqrt(volume_integral_prefix(atoms(atom_index)%radial_ground_state, delta_m**2, region_last))
         region_charge_sr = sqrt(volume_integral_prefix(atoms(atom_index)%radial_ground_state, charge_sr**2, region_last))
         region_delta_charge = sqrt(volume_integral_prefix(atoms(atom_index)%radial_ground_state, delta_charge**2, region_last))

         data_file = trim(output_prefix)//'_'//trim(atoms(atom_index)%element%symbol)//'_'//trim(int_string(isite))//'.dat'
         if (rank == 0) call write_data_file(data_file, reciprocal_obj, atoms(atom_index)%radial_ground_state, &
                                             n_sr, n_pauli, delta_n, m_sr, m_pauli, delta_m, &
                                             charge_sr, charge_pauli, delta_charge, region_last, region_radius)
         if (rank == 0) then
            call write_summary_record(reciprocal_obj, atoms(atom_index)%radial_ground_state, isite, data_file, &
               integrated_sr, integrated_pauli, delta_number, integrated_m_sr, integrated_m_pauli, delta_moment, &
               norm_m_sr, norm_delta_m, norm_charge_sr, norm_delta_charge, region_last, region_radius, &
               region_m_sr, region_delta_m, region_charge_sr, region_delta_charge, occupied_weight, k_weight_sum, &
               pauli_valence, coefficient_weight)
         end if

         deallocate(valence_weighted, pauli_weighted, n_sr, n_pauli, delta_n, m_sr, m_pauli, delta_m, &
                    charge_sr, charge_pauli, delta_charge, abs_m_weight)
      end do

      if (rank == 0) then
         open(newunit=unit, file=trim(output_prefix)//'_summary.dat', status='unknown', action='write', position='append', iostat=ios)
         if (ios == 0) then
            write(unit, '(a)') '# LR-02N summary records are written above this line for each site.'
            close(unit)
         end if
      end if
   end subroutine pauli_ground_state_projection

   pure real(rp) function fermi_dirac(energy, fermi_level, kT) result(occupation)
      real(rp), intent(in) :: energy, fermi_level, kT
      real(rp) :: argument

      if (kT <= 1.0e-12_rp) then
         occupation = merge(1.0_rp, 0.0_rp, energy <= fermi_level)
      else
         argument = (energy - fermi_level)/kT
         if (argument > 50.0_rp) then
            occupation = 0.0_rp
         else if (argument < -50.0_rp) then
            occupation = 1.0_rp
         else
            occupation = 1.0_rp/(exp(argument) + 1.0_rp)
         end if
      end if
   end function fermi_dirac

   real(rp) function origin_density(weighted, radius) result(value)
      real(rp), intent(in) :: weighted(:), radius(:)

      if (size(weighted) < 3 .or. size(radius) /= size(weighted)) then
         error stop 'LR-02N origin density requires three radial points'
      end if
      value = (weighted(2)/radius(2)**2*radius(3) - weighted(3)/radius(3)**2*radius(2))/ &
              (4.0_rp*RADIAL_PI*(radius(3) - radius(2)))
   end function origin_density

   pure real(rp) function radial_integral(state, values) result(value)
      type(radial_ground_state), intent(in) :: state
      real(rp), intent(in) :: values(:)
      integer :: ir

      value = 0.0_rp
      do ir = 1, size(values)
         value = value + radial_simpson_weight(ir, size(values))*state%a*(state%r(ir) + state%b)*values(ir)
      end do
   end function radial_integral

   pure real(rp) function volume_weight(state, ir) result(value)
      type(radial_ground_state), intent(in) :: state
      integer, intent(in) :: ir

      value = radial_simpson_weight(ir, size(state%r))*state%a*(state%r(ir) + state%b)* &
              4.0_rp*RADIAL_PI*state%r(ir)**2
   end function volume_weight

   pure real(rp) function volume_integral(state, values) result(value)
      type(radial_ground_state), intent(in) :: state
      real(rp), intent(in) :: values(:)
      integer :: ir

      value = 0.0_rp
      do ir = 1, size(values)
         value = value + volume_weight(state, ir)*values(ir)
      end do
   end function volume_integral

   pure real(rp) function volume_integral_prefix(state, values, last) result(value)
      type(radial_ground_state), intent(in) :: state
      real(rp), intent(in) :: values(:)
      integer, intent(in) :: last
      integer :: ir

      value = 0.0_rp
      do ir = 1, min(last, size(values))
         value = value + volume_weight(state, ir)*values(ir)
      end do
   end function volume_integral_prefix

   pure integer function magnetic_region_endpoint(state, abs_weight, fraction) result(last)
      type(radial_ground_state), intent(in) :: state
      real(rp), intent(in) :: abs_weight(:), fraction
      real(rp) :: total, cumulative
      integer :: ir

      total = sum(abs_weight)
      cumulative = 0.0_rp
      last = size(abs_weight)
      if (total <= tiny(1.0_rp)) return
      do ir = 1, size(abs_weight)
         cumulative = cumulative + abs_weight(ir)
         if (cumulative >= fraction*total) then
            last = ir
            return
         end if
      end do
   end function magnetic_region_endpoint

   subroutine write_data_file(filename, reciprocal_obj, state, n_sr, n_pauli, delta_n, m_sr, m_pauli, delta_m, &
                              charge_sr, charge_pauli, delta_charge, region_last, region_radius)
      character(len=*), intent(in) :: filename
      type(reciprocal), intent(in) :: reciprocal_obj
      type(radial_ground_state), intent(in) :: state
      real(rp), intent(in) :: n_sr(:, :), n_pauli(:, :), delta_n(:, :), m_sr(:), m_pauli(:), delta_m(:)
      real(rp), intent(in) :: charge_sr(:), charge_pauli(:), delta_charge(:), region_radius
      integer, intent(in) :: region_last
      integer :: unit, ir, ios

      open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
      if (ios /= 0) error stop 'LR-02N could not open radial data file'
      write(unit, '(a)') '# LR-02N scalar-relativistic -> Pauli ground-state projection'
      write(unit, '(a)') '# reference = exact accepted LR-01 radial ground-state snapshot'
      write(unit, '(a)') '# pauli_state = large component only; no GFAC, lower component, or effective amplitude'
      write(unit, '(a)') '# core_treatment = occupied frozen core projected with its accepted large component'
      write(unit, '(a,3(i0,1x))') '# k_mesh = ', reciprocal_obj%nk_mesh
      write(unit, '(a,es24.16)') '# k_weight_sum = ', sum(reciprocal_obj%k_weights)
      write(unit, '(a,es24.16)') '# fermi_level_ry = ', reciprocal_obj%fermi_level
      write(unit, '(a,es24.16)') '# temperature_k = ', reciprocal_obj%temperature
      write(unit, '(a,a)') '# reciprocal_mode = ', trim(reciprocal_obj%reciprocal_mode)
      write(unit, '(a,a)') '# kspace_ham_order = ', trim(reciprocal_obj%kspace_ham_order)
      write(unit, '(a,es24.16)') '# magnetic_region_fraction = ', MAGNETIC_REGION_FRACTION
      write(unit, '(a,i0)') '# magnetic_region_last_point = ', region_last
      write(unit, '(a,es24.16)') '# magnetic_region_rmax_bohr = ', region_radius
      write(unit, '(a)') '# columns: r n_up_SR n_down_SR n_up_P n_down_P delta_n_up delta_n_down m_SR m_P delta_m charge_SR charge_P delta_charge'
      do ir = 1, size(state%r)
         write(unit, '(13(es24.16,1x))') state%r(ir), n_sr(ir, 1), n_sr(ir, 2), n_pauli(ir, 1), n_pauli(ir, 2), &
            delta_n(ir, 1), delta_n(ir, 2), m_sr(ir), m_pauli(ir), delta_m(ir), charge_sr(ir), charge_pauli(ir), delta_charge(ir)
      end do
      close(unit)
   end subroutine write_data_file

   subroutine write_summary_record(reciprocal_obj, state, isite, data_file, integrated_sr, integrated_pauli, delta_number, &
                                   integrated_m_sr, integrated_m_pauli, delta_moment, norm_m_sr, norm_delta_m, &
                                   norm_charge_sr, norm_delta_charge, region_last, region_radius, region_m_sr, region_delta_m, &
                                   region_charge_sr, region_delta_charge, occupied_weight, k_weight_sum, pauli_valence, &
                                   coefficient_weight)
      type(reciprocal), intent(in) :: reciprocal_obj
      type(radial_ground_state), intent(in) :: state
      integer, intent(in) :: isite, region_last
      character(len=*), intent(in) :: data_file
      real(rp), intent(in) :: integrated_sr(2), integrated_pauli(2), delta_number(2)
      real(rp), intent(in) :: integrated_m_sr, integrated_m_pauli, delta_moment
      real(rp), intent(in) :: norm_m_sr, norm_delta_m, norm_charge_sr, norm_delta_charge
      real(rp), intent(in) :: region_radius, region_m_sr, region_delta_m, region_charge_sr, region_delta_charge
      real(rp), intent(in) :: occupied_weight, k_weight_sum, pauli_valence(2), coefficient_weight(2)
      integer :: unit, l, ispin, ios
      real(rp) :: correction_norm, large_norm

      open(newunit=unit, file=trim(data_file)//'.summary', status='replace', action='write', iostat=ios)
      if (ios /= 0) error stop 'LR-02N could not open summary file'
      write(unit, '(a)') '# LR-02N scalar-relativistic -> Pauli numerical closure'
      write(unit, '(a,i0)') 'site = ', isite
      write(unit, '(a,es24.16)') 'fermi_level_ry = ', reciprocal_obj%fermi_level
      write(unit, '(a,es24.16)') 'temperature_k = ', reciprocal_obj%temperature
      write(unit, '(a,3(i0,1x))') 'k_mesh = ', reciprocal_obj%nk_mesh
      write(unit, '(a,es24.16)') 'k_weight_sum = ', k_weight_sum
      write(unit, '(a,es24.16)') 'occupied_valence_weight = ', occupied_weight
      write(unit, '(a,es24.16)') 'direct_coefficient_weight_up = ', coefficient_weight(1)
      write(unit, '(a,es24.16)') 'direct_coefficient_weight_down = ', coefficient_weight(2)
      write(unit, '(a,es24.16)') 'pauli_valence_integral_up = ', pauli_valence(1)
      write(unit, '(a,es24.16)') 'pauli_valence_integral_down = ', pauli_valence(2)
      write(unit, '(a,es24.16)') 'pauli_valence_minus_direct_up = ', pauli_valence(1) - coefficient_weight(1)
      write(unit, '(a,es24.16)') 'pauli_valence_minus_direct_down = ', pauli_valence(2) - coefficient_weight(2)
      write(unit, '(a,a)') 'reciprocal_mode = ', trim(reciprocal_obj%reciprocal_mode)
      write(unit, '(a,a)') 'kspace_ham_order = ', trim(reciprocal_obj%kspace_ham_order)
      if (reciprocal_obj%temperature > 0.0_rp) then
         write(unit, '(a)') 'occupation_rule = Fermi-Dirac at reciprocal temperature'
      else
         write(unit, '(a)') 'occupation_rule = zero-temperature step'
      end if
      write(unit, '(a,a)') 'data_file = ', trim(data_file)
      write(unit, '(a)') 'core_treatment = accepted frozen-core large-component projection; reciprocal states are valence-only'
      write(unit, '(a,es24.16)') 'N_up_SR = ', integrated_sr(1)
      write(unit, '(a,es24.16)') 'N_down_SR = ', integrated_sr(2)
      write(unit, '(a,es24.16)') 'N_up_P = ', integrated_pauli(1)
      write(unit, '(a,es24.16)') 'N_down_P = ', integrated_pauli(2)
      write(unit, '(a,es24.16)') 'Delta_N_up = ', delta_number(1)
      write(unit, '(a,es24.16)') 'Delta_N_down = ', delta_number(2)
      write(unit, '(a,es24.16)') 'M_SR = ', integrated_m_sr
      write(unit, '(a,es24.16)') 'M_P = ', integrated_m_pauli
      write(unit, '(a,es24.16)') 'Delta_M = ', delta_moment
      write(unit, '(a,es24.16)') 'relative_L2_magnetization = ', safe_ratio(norm_delta_m, norm_m_sr)
      write(unit, '(a,es24.16)') 'relative_L2_total_charge = ', safe_ratio(norm_delta_charge, norm_charge_sr)
      write(unit, '(a,es24.16)') 'magnetic_region_fraction = ', MAGNETIC_REGION_FRACTION
      write(unit, '(a,i0)') 'magnetic_region_last_point = ', region_last
      write(unit, '(a,es24.16)') 'magnetic_region_rmax_bohr = ', region_radius
      write(unit, '(a,es24.16)') 'magnetic_region_relative_L2_magnetization = ', safe_ratio(region_delta_m, region_m_sr)
      write(unit, '(a,es24.16)') 'magnetic_region_relative_L2_total_charge = ', safe_ratio(region_delta_charge, region_charge_sr)
      write(unit, '(a)') 'omitted_scalar_relativistic_norm = integral[ l(l+1)/(TMC*r)^2 G1^2 + G2^2 ] dr'
      do ispin = 1, 2
         do l = 0, state%pauli_lmax
            correction_norm = radial_integral(state, state%pauli_sr_correction(:, l + 1, ispin))
            large_norm = radial_integral(state, state%pauli_large(:, l + 1, ispin)**2)
            write(unit, '(a,i0,a,i0,a,es24.16)') 'omitted_sr_norm_l', l, '_spin', ispin, ' = ', correction_norm
            write(unit, '(a,i0,a,i0,a,es24.16)') 'large_component_norm_l', l, '_spin', ispin, ' = ', large_norm
            write(unit, '(a,i0,a,i0,a,es24.16)') 'omitted_sr_ratio_l', l, '_spin', ispin, ' = ', safe_ratio(correction_norm, large_norm)
         end do
      end do
      close(unit)

      open(newunit=unit, file=trim(replace_suffix(data_file, '.dat', '_summary.dat')), status='unknown', action='write', &
           position='append', iostat=ios)
      if (ios == 0) then
         write(unit, '(a,i0,1x,4(es16.8,1x),2(es16.8,1x))') 'site', isite, integrated_sr(1), integrated_sr(2), &
            integrated_pauli(1), integrated_pauli(2), delta_moment, safe_ratio(norm_delta_m, norm_m_sr)
         close(unit)
      end if
   end subroutine write_summary_record

   pure real(rp) function safe_ratio(numerator, denominator) result(value)
      real(rp), intent(in) :: numerator, denominator
      value = 0.0_rp
      if (abs(denominator) > tiny(1.0_rp)) value = numerator/denominator
   end function safe_ratio

   pure function int_string(value) result(text)
      integer, intent(in) :: value
      character(len=32) :: text
      write(text, '(i0)') value
   end function int_string

   pure function replace_suffix(filename, old_suffix, new_suffix) result(output)
      character(len=*), intent(in) :: filename, old_suffix, new_suffix
      character(len=512) :: output
      integer :: position

      output = filename
      position = len_trim(filename) - len_trim(old_suffix) + 1
      if (position >= 1 .and. trim(filename(position:)) == trim(old_suffix)) then
         output = filename(:position - 1)//trim(new_suffix)
      else
         output = trim(filename)//trim(new_suffix)
      end if
   end function replace_suffix

end module pauli_ground_state_projection_mod
