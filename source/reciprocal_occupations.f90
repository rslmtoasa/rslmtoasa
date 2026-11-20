submodule (reciprocal_mod) reciprocal_occupations
#ifdef USE_MPI
   use mpi
#endif
   implicit none

   real(rp), parameter :: occupation_kT_floor = 1.0e-10_rp

contains

   !> Numerically stable Fermi-Dirac occupation.  A tiny positive kT floor
   !> keeps the electron-number equation continuous for nominal T=0 inputs;
   !> the exact same effective kT is used for EF, electron count, and EBAND.
   pure module function fermi_dirac_occupation(eigenvalue, fermi_level, kT) result(occupation)
      real(rp), intent(in) :: eigenvalue, fermi_level, kT
      real(rp) :: occupation
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

   pure module function fermi_dirac_entropy_term(occupation) result(entropy)
      real(rp), intent(in) :: occupation
      real(rp) :: entropy

      if (occupation == 0.0_rp .or. occupation == 1.0_rp) then
         entropy = 0.0_rp
      else
         entropy = -occupation*log(occupation) - (1.0_rp - occupation)*log(1.0_rp - occupation)
      end if
   end function fermi_dirac_entropy_term

   pure module function fermi_dirac_entropy(eigenvalues, weights, fermi_level, kT) result(entropy)
      real(rp), intent(in) :: eigenvalues(:, :), weights(:), fermi_level, kT
      real(rp) :: entropy
      integer :: ik, ib
      real(rp) :: occupation

      entropy = 0.0_rp
      do ik = 1, size(eigenvalues, 2)
         do ib = 1, size(eigenvalues, 1)
            occupation = fermi_dirac_occupation(eigenvalues(ib, ik), fermi_level, kT)
            entropy = entropy + weights(ik)*fermi_dirac_entropy_term(occupation)
         end do
      end do
   end function fermi_dirac_entropy

   !> @brief Gather eigenvalues and k weights of the whole mesh on every rank, one collective per solve.
   !> @details packed(1:nbands, ik) are the eigenvalues and packed(nbands+1, ik) the weight of global k-point ik.
   !>          A replicated mesh copies its local data.  Every later sum runs over the local full copy in global
   !>          k order, so results do not depend on the rank count and need no further communication.
   subroutine gather_occupation_inputs(this, packed)
      class(reciprocal), intent(in) :: this
      real(rp), allocatable, intent(out) :: packed(:, :)

      integer :: nbands, nk_loc
      real(rp), allocatable :: packed_local(:, :)

      if (.not. allocated(this%eigenvalues)) then
         call g_logger%fatal('gather_occupation_inputs: eigenvalues are not available', __FILE__, __LINE__)
      end if
      if (.not. allocated(this%k_workset%weights)) then
         call g_logger%fatal('gather_occupation_inputs: k-point weights are not available', __FILE__, __LINE__)
      end if
      nbands = size(this%eigenvalues, 1)
      nk_loc = size(this%eigenvalues, 2)
      allocate(packed_local(nbands + 1, nk_loc), packed(nbands + 1, this%k_workset%nk_global))
      packed_local(1:nbands, :) = this%eigenvalues
      packed_local(nbands + 1, :) = this%k_workset%weights(1:nk_loc)
      call this%gather_k_columns(nbands + 1, packed_local, packed)
   end subroutine gather_occupation_inputs

   !> @brief Irreducible tetrahedra of the mesh when dos_method='blochl'; left unallocated otherwise (Fermi-Dirac).
   !> @details Built once per mesh and symmetry map and kept on the reciprocal object.
   subroutine blochl_tetrahedra(this, tet, mult)
      class(reciprocal), intent(inout) :: this
      integer, allocatable, intent(out) :: tet(:, :), mult(:)
      integer :: n_tet
      integer, allocatable :: key(:)

      if (trim(this%dos_method) /= 'blochl') return
      if (allocated(this%k_path)) then
         call g_logger%fatal("dos_method='blochl' occupations need a uniform k mesh, not a k-path", __FILE__, __LINE__)
      end if
      key = [this%nk_mesh, this%full_to_irred_k]
      if (allocated(this%occ_tetra_key)) then
         if (size(this%occ_tetra_key) /= size(key)) deallocate(this%occ_tetra_key)
      end if
      if (allocated(this%occ_tetra_key)) then
         if (any(this%occ_tetra_key /= key)) deallocate(this%occ_tetra_key)
      end if
      if (.not. allocated(this%occ_tetra_key)) then
         call this%build_irreducible_tetrahedra(this%occ_tetra, this%occ_tetra_mult, n_tet)
         call move_alloc(key, this%occ_tetra_key)
      end if
      tet = this%occ_tetra
      mult = this%occ_tetra_mult
   end subroutine blochl_tetrahedra

   !> @brief Bloechl-corrected tetrahedron occupations of the packed global mesh at T=0 (PRB 49, 16223 Eqs. 4, 22).
   !> @details occ(ib, ik) is defined so that k_weight(ik)*occ(ib, ik)/sum(k_weight) is the integration weight of band ib
   !>          at irreducible point ik: the corner weights of every full-mesh tetrahedron (irreducible class with
   !>          multiplicity) folded onto the irreducible corners. Corrected weights are not clipped. Points of zero weight get 0.
   subroutine tetra_occupations(this, packed, fermi_level, tet, mult, occ)
      class(reciprocal), intent(in) :: this
      real(rp), intent(in) :: packed(:, :), fermi_level
      integer, intent(in) :: tet(:, :), mult(:)
      real(rp), allocatable, intent(out) :: occ(:, :)

      integer :: it, ib, ic, nbands
      real(rp) :: w(4), dw(4), vol, wsum

      nbands = size(packed, 1) - 1
      allocate(occ(nbands, size(packed, 2)))
      occ = 0.0_rp
      wsum = sum(packed(nbands + 1, :))
      do it = 1, size(tet, 2)
         vol = wsum*real(mult(it), rp)/(6.0_rp*real(product(this%nk_mesh), rp))
         do ib = 1, nbands
            call tetra_integrated_weights(packed(ib, tet(:, it)), fermi_level, w, dw)
            do ic = 1, 4
               occ(ib, tet(ic, it)) = occ(ib, tet(ic, it)) + vol*(w(ic) + dw(ic))
            end do
         end do
      end do
      do ib = 1, nbands
         where (packed(nbands + 1, :) > 0.0_rp)
            occ(ib, :) = occ(ib, :)/packed(nbands + 1, :)
         elsewhere
            occ(ib, :) = 0.0_rp
         end where
      end do
   end subroutine tetra_occupations

   !> @brief Occupations of the current eigensystem at this%fermi_level for the density consumers (dos_method='blochl').
   !> @details Every rank gathers the whole mesh and computes the same weights; a consumer reads its own k columns.
   module subroutine blochl_occupations(this, occ)
      class(reciprocal), intent(inout) :: this
      real(rp), allocatable, intent(out) :: occ(:, :)
      real(rp), allocatable :: packed(:, :)
      integer, allocatable :: tet(:, :), mult(:)

      call gather_occupation_inputs(this, packed)
      call blochl_tetrahedra(this, tet, mult)
      call tetra_occupations(this, packed, this%fermi_level, tet, mult, occ)
   end subroutine blochl_occupations

   !> @brief Sums of weight, weight*occupation and weight*occupation*eigenvalue over the packed global mesh.
   !> @details With tet/mult present (dos_method='blochl') the occupations are the tetrahedron weights, otherwise Fermi-Dirac.
   subroutine occupation_sums(this, packed, fermi_level, kT, sums, tet, mult)
      class(reciprocal), intent(in) :: this
      real(rp), intent(in) :: packed(:, :), fermi_level, kT
      real(rp), intent(out) :: sums(3)
      integer, intent(in), optional :: tet(:, :), mult(:)

      integer :: ik, ib, nbands
      real(rp) :: wk, occupation
      real(rp), allocatable :: occ_t(:, :)

      if (present(tet)) call tetra_occupations(this, packed, fermi_level, tet, mult, occ_t)

      nbands = size(packed, 1) - 1
      sums = 0.0_rp
      do ik = 1, size(packed, 2)
         wk = packed(nbands + 1, ik)
         if (wk < 0.0_rp) then
            call g_logger%fatal('occupation_sums: negative k-point weight', __FILE__, __LINE__)
         end if
         sums(1) = sums(1) + wk
         do ib = 1, nbands
            if (present(tet)) then
               occupation = occ_t(ib, ik)
            else
               occupation = fermi_dirac_occupation(packed(ib, ik), fermi_level, kT)
            end if
            sums(2) = sums(2) + wk*occupation
            sums(3) = sums(3) + wk*occupation*packed(ib, ik)
         end do
      end do

      if (sums(1) <= tiny(1.0_rp)) then
         call g_logger%fatal('occupation_sums: global k-point weight sum is zero', __FILE__, __LINE__)
      end if
   end subroutine occupation_sums

   !> @brief Projection-free occupation and energy accumulation.
   !> @details k_weights may be normalized probabilities or raw multiplicities.
   !>          Both observables are divided by the global raw weight sum.  The sums run over the gathered
   !>          global mesh in global k order (one all-gather on a distributed mesh).
   module subroutine evaluate_eigenvalue_occupations(this, fermi_level, electron_count, band_energy, weight_sum)
      class(reciprocal), intent(inout) :: this
      real(rp), intent(in) :: fermi_level
      real(rp), intent(out) :: electron_count, band_energy
      real(rp), intent(out), optional :: weight_sum

      real(rp) :: kT, sums(3)
      real(rp), allocatable :: packed(:, :)
      integer, allocatable :: tet(:, :), mult(:)

      call gather_occupation_inputs(this, packed)
      call blochl_tetrahedra(this, tet, mult)
      kT = max(this%temperature*kB_Ry_per_K, occupation_kT_floor)
      call occupation_sums(this, packed, fermi_level, kT, sums, tet, mult)

      electron_count = sums(2)/sums(1)
      band_energy = sums(3)/sums(1)
      if (present(weight_sum)) weight_sum = sums(1)
   end subroutine evaluate_eigenvalue_occupations

   !> @brief Solve the projection-free Fermi-level equation N(EF)=N_target.
   !> @details The eigenvalues are gathered once; bisection and the final sums use the local full copy.
   module function find_fermi_level_from_eigenvalues(this, total_electrons) result(fermi_level)
      class(reciprocal), intent(inout) :: this
      real(rp), intent(in) :: total_electrons
      real(rp) :: fermi_level

      integer :: iteration, nbands
      real(rp) :: e_min, e_max, e_mid, electron_count, kT, margin, count_tolerance, snap_tolerance
      real(rp) :: sums(3), extrema(2), gap_below, gap_above
      real(rp), allocatable :: packed(:, :)
      integer, allocatable :: tet(:, :), mult(:)

      call gather_occupation_inputs(this, packed)
      call blochl_tetrahedra(this, tet, mult)
      nbands = size(packed, 1) - 1
      if (total_electrons < -1.0e-12_rp .or. total_electrons > real(nbands, rp) + 1.0e-12_rp) then
         call g_logger%fatal('find_fermi_level_from_eigenvalues: target electron count is outside [0, nbands]', &
                             __FILE__, __LINE__)
      end if

      if (size(packed, 2) > 0) then
         extrema = [minval(packed(1:nbands, :)), maxval(packed(1:nbands, :))]
      else
         extrema = [huge(1.0_rp), -huge(1.0_rp)]
      end if
      if (extrema(1) > extrema(2)) then
         call g_logger%fatal('find_fermi_level_from_eigenvalues: no eigenvalues are owned by the global mesh', &
                             __FILE__, __LINE__)
      end if

      kT = max(this%temperature*kB_Ry_per_K, occupation_kT_floor)
      margin = max(1.0_rp, 64.0_rp*kT)
      e_min = extrema(1) - margin
      e_max = extrema(2) + margin
      count_tolerance = 64*spacing(1.0_rp)*max(1.0_rp, total_electrons)

      snap_tolerance = 1.0e-10_rp*max(1.0_rp, total_electrons)

      do iteration = 1, 256
         e_mid = 0.5_rp*(e_min + e_max)
         call occupation_sums(this, packed, e_mid, kT, sums, tet, mult)
         electron_count = sums(2)/sums(1)
         if (abs(electron_count - total_electrons) <= count_tolerance) exit
         if (e_mid <= e_min .or. e_mid >= e_max) exit
         if (electron_count < total_electrons) then
            e_min = e_mid
         else
            e_max = e_mid
         end if
      end do
      fermi_level = e_mid

      ! Convention for a gapped spectrum: EF is the middle of the gap, independent of where the bisection stopped in it.
      ! It is used when the count there is within 1e-10*N of N (physically negligible, below the 1e-9 convergence check).
      gap_below = maxval(packed(1:nbands, :), mask=packed(1:nbands, :) <= fermi_level)
      gap_above = minval(packed(1:nbands, :), mask=packed(1:nbands, :) > fermi_level)
      if (gap_below > -huge(1.0_rp) .and. gap_above < huge(1.0_rp) .and. gap_above > gap_below) then
         call occupation_sums(this, packed, 0.5_rp*(gap_below + gap_above), kT, sums, tet, mult)
         electron_count = sums(2)/sums(1)
         if (abs(electron_count - total_electrons) <= snap_tolerance) fermi_level = 0.5_rp*(gap_below + gap_above)
      end if

      call occupation_sums(this, packed, fermi_level, kT, sums, tet, mult)
      electron_count = sums(2)/sums(1)
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
