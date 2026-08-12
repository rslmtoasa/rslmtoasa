! TDDFT-11 CPU profiling fixture.  This is intentionally informational: it
! runs small, deterministic reciprocal-space models and prints phase timings,
! but applies no performance threshold in CI.
program test_tddft_cpu_profile
   use precision_mod, only: rp
   use basis_mod, only: nb, norb
   use math_mod, only: init_math_operators
   use response_components_mod, only: RESPONSE_CHARGE, RESPONSE_MZ, RESPONSE_PLUS, RESPONSE_MINUS
   use response_vertices_mod, only: response_channel
   use tddft_chi0_mod, only: tddft_chi0_options, tddft_chi0_result, build_chi_ks_from_eigenpairs
   use tddft_chi0_green_mod, only: green_chi0_options, eigenpair_green_function_provider, &
      build_chi_ks_from_green_functions
   use tddft_dyson_mod, only: tddft_dyson_options, tddft_dyson_result, enhance_tddft_susceptibility
   use tddft_modes_mod, only: tddft_mode_options, tddft_mode_result, analyze_tddft_modes
   use reciprocal_mod, only: reciprocal
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use charge_mod, only: charge
   use control_mod, only: control
   use logger_mod, only: g_logger
   implicit none

   call g_logger%init()
   call init_math_operators()
   call profile_fixture('bccFe_one_site', 1, [4, 4, 1], 96)
   call profile_fixture('fccNi_two_site', 2, [4, 4, 2], 192)
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine profile_fixture(label, nsite, mesh_shape, nw)
      character(len=*), intent(in) :: label
      integer, intent(in) :: nsite, mesh_shape(3), nw
      type(reciprocal) :: recip
      type(hamiltonian), target :: ham
      type(lattice), target :: lat
      type(charge), target :: chg
      type(control), target :: ctl
      real(rp), allocatable :: weights(:), omega(:), kq(:, :), evalq(:, :), trace_loss(:, :)
      complex(rp), allocatable :: evecq(:, :, :), kernel(:, :), xi(:, :, :, :), pair_ops(:, :, :, :), qplus(:, :)
      type(response_channel) :: left(2), right(2)
      type(tddft_chi0_options) :: chi_options
      type(green_chi0_options) :: green_options
      type(tddft_dyson_options) :: dyson_options
      type(tddft_mode_options) :: mode_options
      type(tddft_chi0_result) :: chi, green_chi
      type(eigenpair_green_function_provider) :: green_source
      type(tddft_dyson_result) :: dyson
      type(tddft_mode_result) :: modes
      integer, allocatable :: site_orbital_counts(:)
      integer :: nk, nmat, ik, isite, iw
      logical :: supported
      character(len=160) :: reason
      real(rp) :: q_point(3), t_start, t_stop
      real(rp) :: fourier_seconds, eigensolve_seconds, arbitrary_kq_seconds, pair_operator_seconds

      call setup_profile_model(recip, ham, lat, chg, ctl, nsite)
      nk = product(mesh_shape)
      nmat = nb*nsite
      allocate(weights(nk), omega(nw), kq(3, nk), site_orbital_counts(nsite))
      allocate(pair_ops(nmat, nmat, nsite, nk), qplus(nmat, nmat))
      recip%nk_total = nk
      allocate(recip%k_points(3, nk), recip%k_weights(nk))
      call fill_k_mesh(recip%k_points, mesh_shape)
      recip%k_weights = 1.0_rp/real(nk, rp)
      weights = recip%k_weights
      site_orbital_counts = norb
      q_point = [0.125_rp, 0.0_rp, 0.0_rp]
      kq = recip%k_points + spread(q_point, dim=2, ncopies=nk)
      do iw = 1, nw
         omega(iw) = 0.005_rp + 0.20_rp*real(iw-1, rp)/real(nw-1, rp)
      end do

      call cpu_time(t_start)
      call recip%build_kspace_hamiltonian()
      call cpu_time(t_stop)
      fourier_seconds = t_stop-t_start

      call cpu_time(t_start)
      call recip%diagonalize_hamiltonian()
      call cpu_time(t_stop)
      eigensolve_seconds = t_stop-t_start

      call cpu_time(t_start)
      call recip%calculate_eigenpairs_at_kpoints(kq, evalq, evecq)
      call cpu_time(t_stop)
      arbitrary_kq_seconds = t_stop-t_start

      call cpu_time(t_start)
      do ik = 1, nk
         do isite = 1, nsite
            call recip%build_lmto_pair_potential_at_kpoint(isite, recip%k_points(:, ik), &
               merge(2.0_rp, -2.0_rp, isite == 1), pair_ops(:, :, isite, ik), qplus, supported, reason, q_point)
            if (.not. supported) error stop 'test_tddft_cpu_profile: pair-potential profile fixture unsupported: '//trim(reason)
         end do
      end do
      call cpu_time(t_stop)
      pair_operator_seconds = t_stop-t_start

      left = [response_channel(1, RESPONSE_PLUS), response_channel(1, RESPONSE_MZ)]
      right = [response_channel(1, RESPONSE_MINUS), response_channel(1, RESPONSE_CHARGE)]
      chi_options%eta = 0.006_rp
      chi_options%fermi_level = 0.0_rp
      chi_options%electronic_temperature = 300.0_rp
      chi_options%k_mesh_shape = mesh_shape
      chi_options%transition_batch_size = 128
      call build_chi_ks_from_eigenpairs(weights, recip%eigenvalues, recip%eigenvectors, evalq, evecq, &
         site_orbital_counts, left, right, omega, chi_options, chi)

      call green_source%initialize(recip%eigenvalues, recip%eigenvectors, evalq, evecq)
      green_options%eta = chi_options%eta
      green_options%fermi_level = chi_options%fermi_level
      green_options%electronic_temperature = chi_options%electronic_temperature
      green_options%energy_min = -0.50_rp
      green_options%energy_max = 0.50_rp
      green_options%energy_points = 13
      call build_chi_ks_from_green_functions(green_source, weights(:min(2, nk)), site_orbital_counts, left, right, &
         omega(:min(8, nw)), green_options, green_chi)

      allocate(kernel(2, 2))
      kernel = cmplx(0.0_rp, 0.0_rp, rp)
      kernel(1, 1) = cmplx(0.4_rp, 0.0_rp, rp)
      kernel(2, 2) = cmplx(0.2_rp, 0.0_rp, rp)
      call enhance_tddft_susceptibility(chi%chi, kernel, chi_options%eta, dyson_options, dyson)
      allocate(xi(2, 2, nw, 1), trace_loss(nw, 1))
      xi(:, :, :, 1) = dyson%xi
      trace_loss(:, 1) = dyson%trace_spectral_weight
      call analyze_tddft_modes(omega, xi, trace_loss, chi_options%eta, mode_options, modes)

      call report_profile(label, nsite, nmat, nk, nw, mesh_shape, fourier_seconds, eigensolve_seconds, &
         arbitrary_kq_seconds, pair_operator_seconds, chi, green_chi, dyson, modes)
   end subroutine profile_fixture

   subroutine setup_profile_model(recip, ham, lat, chg, ctl, nsite)
      type(reciprocal), intent(out) :: recip
      type(hamiltonian), target, intent(out) :: ham
      type(lattice), target, intent(out) :: lat
      type(charge), target, intent(out) :: chg
      type(control), target, intent(out) :: ctl
      integer, intent(in) :: nsite
      integer :: isite, i

      call ctl%restore_to_default()
      lat%nrec = nsite; lat%ntype = nsite; lat%nn_max = 3; lat%kk = nsite; lat%nmax = nsite
      allocate(lat%ib(nsite), lat%atlist(nsite), lat%iz(nsite), lat%num(nsite), lat%nn(nsite, 3), &
         lat%sbar(norb, norb, 3, nsite), lat%symbolic_atoms(nsite))
      lat%ib = [(isite, isite=1, nsite)]
      lat%atlist = [(isite, isite=1, nsite)]
      lat%iz = [(isite, isite=1, nsite)]
      lat%num = [(isite, isite=1, nsite)]
      lat%sbar = cmplx(0.0_rp, 0.0_rp, rp)
      do isite = 1, nsite
         if (nsite == 1) then
            lat%nn(isite, :) = [3, 1, 1]
         else if (isite == 1) then
            lat%nn(isite, :) = [3, 2, 2]
         else
            lat%nn(isite, :) = [3, 1, 1]
         end if
         call lat%symbolic_atoms(isite)%restore_to_default()
         lat%symbolic_atoms(isite)%potential%mom = [0.0_rp, 0.0_rp, 1.0_rp]
         lat%symbolic_atoms(isite)%potential%wx0 = cmplx(1.0_rp, 0.0_rp, rp)
         lat%symbolic_atoms(isite)%potential%wx1 = cmplx(0.31_rp, 0.0_rp, rp)
         lat%symbolic_atoms(isite)%potential%cx1 = cmplx(0.23_rp, 0.0_rp, rp)
         do i = 1, norb
            lat%sbar(i, i, 1, isite) = cmplx(0.12_rp, 0.0_rp, rp)
            lat%sbar(i, i, 2, isite) = cmplx(0.37_rp, 0.0_rp, rp)
            lat%sbar(i, i, 3, isite) = cmplx(0.37_rp, 0.0_rp, rp)
         end do
      end do

      chg%lattice => lat
      chg%symbolic_atom => lat%symbolic_atoms
      ham%charge => chg
      ham%lattice => lat
      ham%control => ctl
      ham%hoh = .false.
      ham%ccor_2c = .false.
      ham%hubbard_u_general_check = .false.
      ham%hubbard_v_check = .false.
      ham%hubbard_u_impurity_check = .false.
      ham%local_axis = .false.
      ham%magnetic_representation = 'periodic_nc'
      ham%operator_generation = 1
      allocate(ham%ee(nb, nb, 3, nsite))
      ham%ee = cmplx(0.0_rp, 0.0_rp, rp)
      do isite = 1, nsite
         do i = 1, nb
            ham%ee(i, i, 1, isite) = cmplx(0.05_rp*real(i, rp) + 0.001_rp*real(isite-1, rp), 0.0_rp, rp)
            ham%ee(i, i, 2, isite) = cmplx(0.01_rp*real(i, rp), 0.0_rp, rp)
            ham%ee(i, i, 3, isite) = cmplx(0.01_rp*real(i, rp), 0.0_rp, rp)
         end do
      end do

      recip%hamiltonian => ham
      recip%lattice => lat
      recip%control => ctl
      recip%reciprocal_mode = 'ham_only'
      recip%kspace_ham_order = 'first'
      recip%max_orbs = nb
      recip%dos_method = 'tetrahedron'
      recip%cached_operator_generation = 0
      allocate(recip%ham_vec_type(3, 3, nsite), recip%ham_vec_type_direct(3, 3, nsite))
      recip%ham_vec_type = 0.0_rp
      recip%ham_vec_type_direct = 0.0_rp
      if (nsite == 1) then
         recip%ham_vec_type(1, 2, 1) = 1.0_rp
         recip%ham_vec_type_direct(1, 2, 1) = 1.0_rp
         recip%ham_vec_type(1, 3, 1) = -1.0_rp
         recip%ham_vec_type_direct(1, 3, 1) = -1.0_rp
      else
         recip%ham_vec_type(1, 2, 1) = 0.25_rp
         recip%ham_vec_type_direct(1, 2, 1) = 0.25_rp
         recip%ham_vec_type(1, 3, 1) = -0.75_rp
         recip%ham_vec_type_direct(1, 3, 1) = -0.75_rp
         recip%ham_vec_type(1, 2, 2) = 0.75_rp
         recip%ham_vec_type_direct(1, 2, 2) = 0.75_rp
         recip%ham_vec_type(1, 3, 2) = -0.25_rp
         recip%ham_vec_type_direct(1, 3, 2) = -0.25_rp
      end if
   end subroutine setup_profile_model

   subroutine fill_k_mesh(k_points, mesh_shape)
      real(rp), intent(out) :: k_points(:, :)
      integer, intent(in) :: mesh_shape(3)
      integer :: ix, iy, iz, ik

      ik = 0
      do iz = 1, mesh_shape(3)
         do iy = 1, mesh_shape(2)
            do ix = 1, mesh_shape(1)
               ik = ik + 1
               k_points(:, ik) = [real(ix-1, rp)/real(mesh_shape(1), rp)-0.5_rp, &
                  real(iy-1, rp)/real(mesh_shape(2), rp)-0.5_rp, &
                  real(iz-1, rp)/real(mesh_shape(3), rp)-0.5_rp]
            end do
         end do
      end do
   end subroutine fill_k_mesh

   subroutine report_profile(label, nsite, nmat, nk, nw, mesh_shape, fourier_seconds, eigensolve_seconds, &
      arbitrary_kq_seconds, pair_operator_seconds, chi, green_chi, dyson, modes)
      character(len=*), intent(in) :: label
      integer, intent(in) :: nsite, nmat, nk, nw, mesh_shape(3)
      real(rp), intent(in) :: fourier_seconds, eigensolve_seconds, arbitrary_kq_seconds, pair_operator_seconds
      type(tddft_chi0_result), intent(in) :: chi, green_chi
      type(tddft_dyson_result), intent(in) :: dyson
      type(tddft_mode_result), intent(in) :: modes
      real(rp) :: hk_mib, eig_mib, kq_mib, pair_mib, response_mib, payload_mib

      hk_mib = complex_mib(nmat, nmat, nk)
      eig_mib = real_mib(nmat, nk) + complex_mib(nmat, nmat, nk)
      kq_mib = real_mib(nmat, nk) + complex_mib(nmat, nmat, nk)
      pair_mib = complex_mib(nmat, nmat, nsite, nk) + complex_mib(nmat, nmat)
      response_mib = complex_mib(2, 2, nw) + 2.0_rp*real_mib(2, 2, nw)
      payload_mib = hk_mib + eig_mib + kq_mib + pair_mib + response_mib

      write (*, '(a,1x,a,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,a,i0,a,i0,1x,a,i0)') 'PROFILE_DIMENSIONS', trim(label), &
         'sites=', nsite, 'spinor_basis=', nmat, 'nk=', nk, 'mesh=', mesh_shape(1), 'x', mesh_shape(2), 'x', &
         mesh_shape(3), 'nw=', nw
      write (*, '(a,1x,a,1x,a,es12.4,1x,a,es12.4,1x,a,es12.4,1x,a,es12.4)') 'PROFILE_RECIPROCAL', trim(label), &
         'fourier_assembly=', fourier_seconds, 'k_eigensolution=', eigensolve_seconds, &
         'arbitrary_kq_assembly_eigensolution=', arbitrary_kq_seconds, 'pair_operator_construction=', pair_operator_seconds
      write (*, '(a,1x,a,1x,a,es12.4,1x,a,es12.4,1x,a,es12.4,1x,a,es12.4,1x,a,es12.4,1x,a,es12.4,1x,a,es12.4)') &
         'PROFILE_TDDFT', trim(label), 'vertex_construction=', chi%metadata%vertex_cpu_seconds, &
         'denominator_generation=', chi%metadata%denominator_cpu_seconds, 'response_accumulation=', &
         chi%metadata%accumulation_cpu_seconds, 'green_energy_integration=', &
         green_chi%metadata%green_energy_integration_cpu_seconds, 'dyson_solve=', dyson%metadata%solve_cpu_seconds, &
         'dyson_diagonalization=', dyson%metadata%diagonalization_cpu_seconds, 'mode_analysis=', modes%analysis_cpu_seconds
      write (*, '(a,1x,a,1x,a,es12.4,1x,a,es12.4,1x,a,es12.4,1x,a,es12.4,1x,a,es12.4,1x,a,es12.4)') &
         'PROFILE_MEMORY_MIB', trim(label), 'hk=', hk_mib, 'normal_eigenpairs=', eig_mib, 'arbitrary_kq_eigenpairs=', kq_mib, &
         'pair_operators_and_workspace=', pair_mib, 'response=', response_mib, 'principal_payload=', payload_mib
   end subroutine report_profile

   pure real(rp) function complex_mib(d1, d2, d3, d4) result(value)
      integer, intent(in) :: d1, d2
      integer, intent(in), optional :: d3, d4
      integer :: elements

      elements = d1*d2
      if (present(d3)) elements = elements*d3
      if (present(d4)) elements = elements*d4
      value = 16.0_rp*real(elements, rp)/(1024.0_rp**2)
   end function complex_mib

   pure real(rp) function real_mib(d1, d2, d3) result(value)
      integer, intent(in) :: d1, d2
      integer, intent(in), optional :: d3
      integer :: elements

      elements = d1*d2
      if (present(d3)) elements = elements*d3
      value = 8.0_rp*real(elements, rp)/(1024.0_rp**2)
   end function real_mib

end program test_tddft_cpu_profile
