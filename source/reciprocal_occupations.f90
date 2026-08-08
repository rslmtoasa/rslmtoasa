submodule (reciprocal_mod) reciprocal_occupations
#ifdef USE_MPI
   use mpi
#endif
   implicit none

   real(rp), parameter :: kB_Ry_per_K = 6.3336814e-6_rp
   real(rp), parameter :: occupation_kT_floor = 1.0e-10_rp

contains

   !> Numerically stable Fermi-Dirac occupation.  A tiny positive kT floor
   !> keeps the electron-number equation continuous for nominal T=0 inputs;
   !> the exact same effective kT is used for EF, electron count, and EBAND.
   pure real(rp) function fermi_dirac_occupation(eigenvalue, fermi_level, kT) result(occupation)
      real(rp), intent(in) :: eigenvalue, fermi_level, kT
      real(rp) :: argument

      argument = (eigenvalue - fermi_level)/max(kT, occupation_kT_floor)
      if (argument >= 50.0_rp) then
         occupation = 0.0_rp
      else if (argument <= -50.0_rp) then
         occupation = 1.0_rp
      else
         occupation = 1.0_rp/(exp(argument) + 1.0_rp)
      end if
   end function fermi_dirac_occupation

   !> @brief Projection-free occupation and energy accumulation.
   !> @details k_weights may be normalized probabilities or raw multiplicities.
   !>          Both observables are divided by the global raw weight sum.  A
   !>          distributed mesh owns each eigenvalue exactly once and performs
   !>          one three-scalar all-reduction; a replicated mesh performs none.
   module subroutine evaluate_eigenvalue_occupations(this, fermi_level, electron_count, band_energy, weight_sum)
      class(reciprocal), intent(in) :: this
      real(rp), intent(in) :: fermi_level
      real(rp), intent(out) :: electron_count, band_energy
      real(rp), intent(out), optional :: weight_sum

      integer :: ik, ib, ik_global
      real(rp) :: wk, occupation, kT
      real(rp) :: sums(3)

      if (.not. allocated(this%eigenvalues)) then
         call g_logger%fatal('evaluate_eigenvalue_occupations: eigenvalues are not available', __FILE__, __LINE__)
      end if
      if (.not. allocated(this%k_weights)) then
         call g_logger%fatal('evaluate_eigenvalue_occupations: k-point weights are not available', __FILE__, __LINE__)
      end if

      kT = max(this%temperature*kB_Ry_per_K, occupation_kT_floor)
      sums = 0.0_rp

      do ik = 1, size(this%eigenvalues, 2)
         ik_global = local_k_index_to_global(this, ik)
         if (ik_global < 1 .or. ik_global > size(this%k_weights)) then
            call g_logger%fatal('evaluate_eigenvalue_occupations: local/global k-point map is invalid', __FILE__, __LINE__)
         end if
         wk = this%k_weights(ik_global)
         if (wk < 0.0_rp) then
            call g_logger%fatal('evaluate_eigenvalue_occupations: negative k-point weight', __FILE__, __LINE__)
         end if
         sums(1) = sums(1) + wk
         do ib = 1, size(this%eigenvalues, 1)
            occupation = fermi_dirac_occupation(this%eigenvalues(ib, ik), fermi_level, kT)
            sums(2) = sums(2) + wk*occupation
            sums(3) = sums(3) + wk*occupation*this%eigenvalues(ib, ik)
         end do
      end do

#ifdef USE_MPI
      if (this%k_mesh_distributed_active) then
         call MPI_ALLREDUCE(MPI_IN_PLACE, sums, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      end if
#endif

      if (sums(1) <= tiny(1.0_rp)) then
         call g_logger%fatal('evaluate_eigenvalue_occupations: global k-point weight sum is zero', __FILE__, __LINE__)
      end if

      electron_count = sums(2)/sums(1)
      band_energy = sums(3)/sums(1)
      if (present(weight_sum)) weight_sum = sums(1)
   end subroutine evaluate_eigenvalue_occupations

   !> @brief Solve the projection-free Fermi-level equation N(EF)=N_target.
   module function find_fermi_level_from_eigenvalues(this, total_electrons) result(fermi_level)
      class(reciprocal), intent(in) :: this
      real(rp), intent(in) :: total_electrons
      real(rp) :: fermi_level

      integer :: iteration, nbands
      real(rp) :: e_min, e_max, e_mid, electron_count, band_energy
      real(rp) :: extrema(2), kT, margin, count_tolerance

      if (.not. allocated(this%eigenvalues)) then
         call g_logger%fatal('find_fermi_level_from_eigenvalues: eigenvalues are not available', __FILE__, __LINE__)
      end if
      nbands = size(this%eigenvalues, 1)
      if (total_electrons < -1.0e-12_rp .or. total_electrons > real(nbands, rp) + 1.0e-12_rp) then
         call g_logger%fatal('find_fermi_level_from_eigenvalues: target electron count is outside [0, nbands]', &
                             __FILE__, __LINE__)
      end if

      if (size(this%eigenvalues, 2) > 0) then
         extrema = [minval(this%eigenvalues), maxval(this%eigenvalues)]
      else
         extrema = [huge(1.0_rp), -huge(1.0_rp)]
      end if
#ifdef USE_MPI
      if (this%k_mesh_distributed_active) then
         call MPI_ALLREDUCE(MPI_IN_PLACE, extrema(1), 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, extrema(2), 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      end if
#endif
      if (extrema(1) > extrema(2)) then
         call g_logger%fatal('find_fermi_level_from_eigenvalues: no eigenvalues are owned by the global mesh', &
                             __FILE__, __LINE__)
      end if

      kT = max(this%temperature*kB_Ry_per_K, occupation_kT_floor)
      margin = max(1.0_rp, 64.0_rp*kT)
      e_min = extrema(1) - margin
      e_max = extrema(2) + margin
      count_tolerance = max(5.0e-12_rp, 5.0e-12_rp*max(1.0_rp, total_electrons))

      do iteration = 1, 256
         e_mid = 0.5_rp*(e_min + e_max)
         call this%evaluate_eigenvalue_occupations(e_mid, electron_count, band_energy)
         if (abs(electron_count - total_electrons) <= count_tolerance) exit
         if (electron_count < total_electrons) then
            e_min = e_mid
         else
            e_max = e_mid
         end if
      end do
      fermi_level = e_mid

      call this%evaluate_eigenvalue_occupations(fermi_level, electron_count, band_energy)
      if (abs(electron_count - total_electrons) > 1.0e-9_rp*max(1.0_rp, total_electrons)) then
         call g_logger%fatal('find_fermi_level_from_eigenvalues: electron-number bisection did not converge', &
                             __FILE__, __LINE__)
      end if
   end function find_fermi_level_from_eigenvalues

   !> @brief Canonical k-space band energy from the current eigensystem.
   module function calculate_canonical_band_energy(this, find_fermi, electron_count) result(eband)
      class(reciprocal), intent(inout) :: this
      logical, intent(in), optional :: find_fermi
      real(rp), intent(out), optional :: electron_count
      real(rp) :: eband

      logical :: solve_fermi
      real(rp) :: nelect, weight_sum, electron_error
      character(len=320) :: message

      solve_fermi = this%auto_find_fermi
      if (present(find_fermi)) solve_fermi = find_fermi
      if (solve_fermi) then
         if (this%total_electrons <= 0.0_rp) then
            call g_logger%fatal('calculate_canonical_band_energy: positive total_electrons required to solve EF', &
                                __FILE__, __LINE__)
         end if
         this%fermi_level = this%find_fermi_level_from_eigenvalues(this%total_electrons)
      end if

      call this%evaluate_eigenvalue_occupations(this%fermi_level, nelect, eband, weight_sum)
      this%canonical_electron_count = nelect
      this%canonical_band_energy = eband
      this%canonical_weight_sum = weight_sum
      this%canonical_energy_valid = .true.
      if (present(electron_count)) electron_count = nelect

      electron_error = nelect - this%total_electrons
      if (solve_fermi) then
         if (abs(electron_error) > 1.0e-9_rp*max(1.0_rp, this%total_electrons)) then
            call g_logger%fatal('calculate_canonical_band_energy: canonical electron conservation failed', &
                                __FILE__, __LINE__)
         end if
      else if (this%total_electrons > 0.0_rp .and. &
               abs(electron_error) > 1.0e-6_rp*max(1.0_rp, this%total_electrons)) then
         call g_logger%warning('calculate_canonical_band_energy: fixed EF does not match requested electron count', &
                               __FILE__, __LINE__)
      end if

      write(message, '(A,ES16.8,A,ES16.8,A,ES12.4,A,ES16.8,A,ES16.8)') &
         'Canonical k-space occupations: EF=', this%fermi_level, ' Ry, N=', nelect, &
         ', dN=', electron_error, ', weight_sum(raw)=', weight_sum, ', EBAND=', eband
      call root_info(trim(message), __FILE__, __LINE__)
   end function calculate_canonical_band_energy

end submodule reciprocal_occupations
