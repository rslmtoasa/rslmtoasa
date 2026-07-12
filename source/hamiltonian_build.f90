submodule(hamiltonian_mod) hamiltonian_build

contains

   !> @brief Construct a Hamiltonian object from an initialized charge object.
   !> @details Wires the Hamiltonian to charge, lattice, and control state, resets
   !>          all owned arrays/flags, then reads the &hamiltonian namelist. This
   !>          object supplies real-space hopping blocks to recursion and k-space
   !>          Fourier paths.
   !> @param[in] charge_obj Charge object containing lattice and potential state.
   !> @return Initialized Hamiltonian object.
   module function constructor(charge_obj) result(obj)
      type(hamiltonian) :: obj
      type(charge), target, intent(in) :: charge_obj

      obj%charge => charge_obj
      obj%lattice => charge_obj%lattice
      obj%control => charge_obj%lattice%control

      call obj%restore_to_default()
      call obj%build_from_file()
   end function constructor

   !> @brief Finalize a Hamiltonian object.
   !> @details Releases spin-orbit, hopping, overlap, velocity, Hubbard, CCOR,
   !>          export, and cache arrays owned by the object.
   !> @param[inout] this Hamiltonian object being finalized.
   module subroutine destructor(this)
      type(hamiltonian) :: this
#ifdef USE_SAFE_ALLOC
      if (allocated(this%lsham)) call g_safe_alloc%deallocate('hamiltonian.lsham', this%lsham)
      if (allocated(this%tmat)) call g_safe_alloc%deallocate('hamiltonian.tmat', this%tmat)
      if (allocated(this%ee)) call g_safe_alloc%deallocate('hamiltonian.ee', this%ee)
      if (allocated(this%eecc)) call g_safe_alloc%deallocate('hamiltonian.eecc', this%eecc)
      if (allocated(this%hmag)) call g_safe_alloc%deallocate('hamiltonian.hmag', this%hmag)
      if (allocated(this%hhmag)) call g_safe_alloc%deallocate('hamiltonian.hhmag', this%hhmag)
      if (allocated(this%hall)) call g_safe_alloc%deallocate('hamiltonian.hall', this%hall)
      if (allocated(this%hallcc)) call g_safe_alloc%deallocate('hamiltonian.hallcc', this%hallcc)
      if (allocated(this%eeo)) call g_safe_alloc%deallocate('hamiltonian.eeo', this%eeo)
      if (allocated(this%eeoee)) call g_safe_alloc%deallocate('hamiltonian.eeoee', this%eeoee)
      if (allocated(this%hallo)) call g_safe_alloc%deallocate('hamiltonian.hallo', this%hallo)
      if (allocated(this%obarm)) call g_safe_alloc%deallocate('hamiltonian.obarm', this%obarm)
      if (allocated(this%enim)) call g_safe_alloc%deallocate('hamiltonian.enim', this%enim)
      if (allocated(this%ee_glob)) call g_safe_alloc%deallocate('hamiltonian.ee_glob', this%ee_glob)
      if (allocated(this%eecc_glob)) call g_safe_alloc%deallocate('hamiltonian.eecc_glob', this%eecc_glob)
      if (allocated(this%eeo_glob)) call g_safe_alloc%deallocate('hamiltonian.eeo_glob', this%eeo_glob)
      if (allocated(this%hallcc_glob)) call g_safe_alloc%deallocate('hamiltonian.hallcc_glob', this%hallcc_glob)
      if (allocated(this%enim_glob)) call g_safe_alloc%deallocate('hamiltonian.enim_glob', this%enim_glob)
      if (allocated(this%v_a)) call g_safe_alloc%deallocate('hamiltonian.v_a', this%v_a)
      if (allocated(this%v_b)) call g_safe_alloc%deallocate('hamiltonian.v_b', this%v_b)
      if (allocated(this%vo_a)) call g_safe_alloc%deallocate('hamiltonian.vo_a', this%vo_a)
      if (allocated(this%vo_b)) call g_safe_alloc%deallocate('hamiltonian.vo_b', this%vo_b)
      if (allocated(this%jso_a)) call g_safe_alloc%deallocate('hamiltonian.jso_a', this%jso_a)
      if (allocated(this%jlo_a)) call g_safe_alloc%deallocate('hamiltonian.jlo_a', this%jlo_a)
      if (allocated(this%h_sparse)) call g_safe_alloc%deallocate('hamiltonian.h_sparse', this%h_sparse)
      if (allocated(this%velocity_scale)) call g_safe_alloc%deallocate('hamiltonian.velocity_scale', this%velocity_scale)
      if (allocated(this%hxc)) call g_safe_alloc%deallocate('hamiltonian.hxc', this%hxc)
#else
      if (allocated(this%lsham)) deallocate (this%lsham)
      if (allocated(this%tmat)) deallocate (this%tmat)
      if (allocated(this%ee)) deallocate (this%ee)
      if (allocated(this%eecc)) deallocate (this%eecc)
      if (allocated(this%eeo)) deallocate (this%eeo)
      if (allocated(this%eeoee)) deallocate (this%eeoee)
      if (allocated(this%hmag)) deallocate (this%hmag)
      if (allocated(this%hhmag)) deallocate (this%hhmag)
      if (allocated(this%hall)) deallocate (this%hall)
      if (allocated(this%hallcc)) deallocate (this%hallcc)
      if (allocated(this%hallo)) deallocate (this%hallo)
      if (allocated(this%obarm)) deallocate (this%obarm)
      if (allocated(this%enim)) deallocate (this%enim)
      if (allocated(this%ee_glob)) deallocate (this%ee_glob)
      if (allocated(this%eecc_glob)) deallocate (this%eecc_glob)
      if (allocated(this%eeo_glob)) deallocate (this%eeo_glob)
      if (allocated(this%hallcc_glob)) deallocate (this%hallcc_glob)
      if (allocated(this%enim_glob)) deallocate (this%enim_glob)
      if (allocated(this%v_a)) deallocate(this%v_a)
      if (allocated(this%v_b)) deallocate(this%v_b)
      if (allocated(this%vo_a)) deallocate(this%vo_a)
      if (allocated(this%vo_b)) deallocate(this%vo_b)
      if (allocated(this%jso_a)) deallocate(this%jso_a)
      if (allocated(this%jlo_a)) deallocate(this%jlo_a)
      if (allocated(this%h_sparse)) deallocate(this%h_sparse)
      if (allocated(this%velocity_scale)) deallocate(this%velocity_scale)
      if (allocated(this%hxc)) deallocate(this%hxc)
      if (allocated(this%hubbard_u_pot)) deallocate(this%hubbard_u_pot)
      if (allocated(this%hubbard_u_impurity)) deallocate(this%hubbard_u_impurity)
      if (allocated(this%hubbard_j_impurity)) deallocate(this%hubbard_j_impurity)
      if (allocated(this%hubbard_u_sc)) deallocate(this%hubbard_u_sc)
      if (allocated(this%hubbard_v)) deallocate(this%hubbard_v)
      if (allocated(this%hubbard_v_pot)) deallocate(this%hubbard_v_pot)
#endif
   end subroutine destructor

   !> @brief Read the &hamiltonian namelist and install Hamiltonian options.
   !> @details Parses HOH, local-axis rotation, CCOR, velocity directions, spectral
   !>          bounds, export mode, and Hubbard U/J/V inputs. Hubbard inputs are
   !>          accepted in eV and converted to internal Ry units.
   !> @param[inout] this Hamiltonian object whose control%fname selects the input file.
   !> @note This is a true input boundary and may raise fatal diagnostics for invalid options.
   module subroutine build_from_file(this)
      class(hamiltonian), intent(inout) :: this

      ! variables associated with the reading processes
      integer :: iostatus, funit, i, j, l, li, lj, n_l_in
      integer :: nimp_in
      integer :: na
      integer :: ntype_nml, nrec_nml
      logical :: legacy_uj_present
      logical :: has_hubbard_general, has_legacy_uj, has_impurity_uj, has_hubbard_v, has_hubbard_u_sc
      character(len=1) :: orbch
      real(rp), parameter :: hubbard_sentinel = -9.87654321e30_rp
      integer, parameter :: hubbard_sc_sentinel = -999

      include 'include_codes/namelists/hamiltonian.f90'

      hoh = this%hoh
      local_axis = this%local_axis
      orb_pol = this%orb_pol
      ccor_2c = this%ccor_2c
      ccor_elin = this%ccor_elin
      ccor_vmt_mode = this%ccor_vmt_mode
      ccor_debug = this%ccor_debug
      ccor_strict = this%ccor_strict
      v_alpha(:) = this%v_alpha(:)
      v_beta(:) = this%v_beta(:)
      q_ss(:) = this%q_ss(:)
      theta_ss = this%theta_ss
      js_alpha = this%js_alpha
      jl_alpha = this%jl_alpha
      call move_alloc(this%velocity_scale, velocity_scale)
      hubbard_u_potential_form = this%hubbard_u_potential_form
      bounds_algorithm = this%bounds%algorithm
      bounds_scaling = this%bounds%scaling
      export = this%export

      ! Robust namelist buffers: pre-allocate allocatable inputs so namelist
      ! read does not depend on compiler-specific auto-allocation behavior.
      ntype_nml = max(1, this%lattice%ntype)
      nrec_nml = max(1, this%lattice%nrec)

      if (.not. allocated(hubbard_u_general)) allocate(hubbard_u_general(ntype_nml, 4))
      if (.not. allocated(hubbard_j_general)) allocate(hubbard_j_general(ntype_nml, 4))
      if (.not. allocated(hubbard_u)) allocate(hubbard_u(nrec_nml, 4))
      if (.not. allocated(hubbard_j)) allocate(hubbard_j(nrec_nml, 4))
      if (.not. allocated(hubbard_u_impurity)) allocate(hubbard_u_impurity(nrec_nml, 4))
      if (.not. allocated(hubbard_j_impurity)) allocate(hubbard_j_impurity(nrec_nml, 4))
      if (.not. allocated(hubbard_u_sc)) allocate(hubbard_u_sc(ntype_nml, 4))
      if (.not. allocated(hubbard_v)) allocate(hubbard_v(ntype_nml, ntype_nml, 4, 4))
      if (.not. allocated(uj_orb)) allocate(uj_orb(nrec_nml))

      hubbard_u_general(:, :) = hubbard_sentinel
      hubbard_j_general(:, :) = hubbard_sentinel
      hubbard_u(:, :) = hubbard_sentinel
      hubbard_j(:, :) = hubbard_sentinel
      hubbard_u_impurity(:, :) = hubbard_sentinel
      hubbard_j_impurity(:, :) = hubbard_sentinel
      hubbard_v(:, :, :, :) = hubbard_sentinel
      hubbard_u_sc(:, :) = hubbard_sc_sentinel
      uj_orb(:) = ''

      ! Reading
      open (newunit=funit, file=this%control%fname, action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) then
         call g_logger%fatal('file '//trim(this%control%fname)//' not found', __FILE__, __LINE__)
      end if

      read (funit, nml=hamiltonian, iostat=iostatus)
      if (iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
         call g_logger%error('Error while reading namelist', __FILE__, __LINE__)
         call g_logger%error('iostatus = '//int2str(iostatus), __FILE__, __LINE__)
      end if
      close (funit)

      this%hoh = hoh
      this%local_axis = local_axis
      this%orb_pol = orb_pol
         this%ccor_2c = ccor_2c
         this%ccor_elin = ccor_elin
         this%ccor_vmt_mode = lower(trim(ccor_vmt_mode))
         if (len_trim(this%ccor_vmt_mode) == 0) this%ccor_vmt_mode = 'surface_scalar'
         this%ccor_debug = ccor_debug
         this%ccor_strict = ccor_strict
         if (this%ccor_vmt_mode /= 'surface_scalar' .and. this%ccor_vmt_mode /= 'vmad_scalar' .and. &
             this%ccor_vmt_mode /= 'pair_surface') then
            call g_logger%fatal("Invalid ccor_vmt_mode. Use 'surface_scalar', 'vmad_scalar', or 'pair_surface'.", __FILE__, __LINE__)
         end if
         this%v_alpha(:) = v_alpha(:)
      this%v_beta(:) = v_beta(:)
      ! q_ss is given in the &hamiltonian namelist in Cartesian units of 2*pi/alat
      ! (q_ss=0.5 along a simple cubic direction is the zone boundary at pi/alat);
      ! theta_ss is given in degrees. The Hamiltonian phase convention is
      ! phase = 2*pi*q_ss.(bond/alat), with theta_ss stored internally in radians.
      this%q_ss(:) = q_ss(:)
      this%theta_ss = theta_ss * pi / 180.0_rp
      this%js_alpha = js_alpha
      this%jl_alpha = jl_alpha
      this%hubbard_u_potential_form = lower(trim(hubbard_u_potential_form))
      if (this%hubbard_u_potential_form /= 'liechtenstein' .and. this%hubbard_u_potential_form /= 'acbn0') then
         call g_logger%fatal("Invalid hubbard_u_potential_form. Use 'liechtenstein' or 'acbn0'.", __FILE__, __LINE__)
      end if
      this%bounds%algorithm = lower(trim(bounds_algorithm))
      if (len_trim(this%bounds%algorithm) == 0) this%bounds%algorithm = 'none'
      this%export = lower(trim(export))
      if (len_trim(this%export) == 0) this%export = 'none'
      if (bounds_scaling > 0.0_rp) then
         this%bounds%scaling = bounds_scaling
      else
         this%bounds%scaling = 1.05_rp
      end if
      call move_alloc(velocity_scale, this%velocity_scale)

      has_hubbard_general = any(hubbard_u_general /= hubbard_sentinel) .or. any(hubbard_j_general /= hubbard_sentinel)
      has_legacy_uj = any(uj_orb /= '') .and. (any(hubbard_u /= hubbard_sentinel) .or. any(hubbard_j /= hubbard_sentinel))
      has_impurity_uj = any(hubbard_u_impurity /= hubbard_sentinel) .or. any(hubbard_j_impurity /= hubbard_sentinel)
      has_hubbard_v = any(hubbard_v /= hubbard_sentinel)
      has_hubbard_u_sc = any(hubbard_u_sc /= hubbard_sc_sentinel)

      ! Optional Hubbard inputs from &hamiltonian are provided in eV.
      ! Converted below to internal Ry units:
      ! - hubbard_u_general / hubbard_j_general
      ! - legacy hubbard_u / hubbard_j via uj_orb mapping
      ! - hubbard_u_impurity / hubbard_j_impurity
      ! - hubbard_v (intersite)
      if (has_hubbard_general) then
         n_l_in = min(size(hubbard_u_general, 2), size(hubbard_j_general, 2))
         do i = 1, min(this%lattice%ntype, size(hubbard_u_general, 1), size(hubbard_j_general, 1))
            do l = 1, min(n_l_in, size(this%lattice%symbolic_atoms(i)%potential%hubbard_u))
               if (hubbard_u_general(i, l) == hubbard_sentinel .or. hubbard_j_general(i, l) == hubbard_sentinel) cycle
               this%lattice%symbolic_atoms(i)%potential%hubbard_u(l) = hubbard_u_general(i, l)/ry2ev
               this%lattice%symbolic_atoms(i)%potential%hubbard_j(l) = hubbard_j_general(i, l)/ry2ev
            end do
         end do
      end if

      ! Legacy LDA+U(+J) input path from lda_u branch:
      ! map (hubbard_u/hubbard_j + uj_orb) into per-atom per-l potential arrays.
      legacy_uj_present = .false.
      if (has_legacy_uj) then
         do i = 1, min(this%lattice%nrec, size(uj_orb), size(hubbard_u, 1), size(hubbard_j, 1))
            na = this%lattice%nbulk + i
            if (na < 1 .or. na > this%lattice%ntype) cycle
            do j = 1, len_trim(uj_orb(i))
               if (j > size(hubbard_u, 2) .or. j > size(hubbard_j, 2)) exit
               orbch = uj_orb(i)(j:j)
               select case (orbch)
               case ('s', 'S')
                  l = 1
               case ('p', 'P')
                  l = 2
               case ('d', 'D')
                  l = 3
               case ('f', 'F')
                  l = 4
               case default
                  cycle
               end select
               if (l <= size(this%lattice%symbolic_atoms(na)%potential%hubbard_u)) then
                  if (hubbard_u(i, j) == hubbard_sentinel .or. hubbard_j(i, j) == hubbard_sentinel) cycle
                  this%lattice%symbolic_atoms(na)%potential%hubbard_u(l) = hubbard_u(i, j)/ry2ev
                  this%lattice%symbolic_atoms(na)%potential%hubbard_j(l) = hubbard_j(i, j)/ry2ev
                  if (abs(hubbard_u(i, j)) > 1.0e-10_rp .or. abs(hubbard_j(i, j)) > 1.0e-10_rp) legacy_uj_present = .true.
               end if
            end do
         end do
      end if
      if (legacy_uj_present) then
         call g_logger%info('Legacy uj_orb + hubbard_u/hubbard_j input detected and mapped to symbolic-atom Hubbard channels.', __FILE__, __LINE__)
      end if

      this%hubbard_u_impurity_check = .false.
      if (has_impurity_uj) then
         nimp_in = min(size(hubbard_u_impurity, 1), size(hubbard_j_impurity, 1))
         if (nimp_in > 0) then
            this%hubbard_u_impurity(:, :) = 0.0_rp
            this%hubbard_j_impurity(:, :) = 0.0_rp
            do i = 1, min(this%lattice%nrec, nimp_in)
               do l = 1, min(4, size(hubbard_u_impurity, 2), size(hubbard_j_impurity, 2))
                  if (hubbard_u_impurity(i, l) == hubbard_sentinel .or. hubbard_j_impurity(i, l) == hubbard_sentinel) cycle
                  this%hubbard_u_impurity(i, l) = hubbard_u_impurity(i, l)/ry2ev
                  this%hubbard_j_impurity(i, l) = hubbard_j_impurity(i, l)/ry2ev
               end do
            end do
         end if
      end if

      ! Apply impurity-specific U/J to recursive atoms (nbulk+1:ntype)
      if (this%lattice%nrec > 0) then
         do i = 1, this%lattice%nrec
            if (maxval(abs(this%hubbard_u_impurity(i, :))) > 1.0e-10_rp .or. &
                maxval(abs(this%hubbard_j_impurity(i, :))) > 1.0e-10_rp) then
               this%hubbard_u_impurity_check = .true.
               do l = 1, min(4, size(this%lattice%symbolic_atoms(this%lattice%nbulk + i)%potential%hubbard_u))
                  this%lattice%symbolic_atoms(this%lattice%nbulk + i)%potential%hubbard_u(l) = this%hubbard_u_impurity(i, l)
                  this%lattice%symbolic_atoms(this%lattice%nbulk + i)%potential%hubbard_j(l) = this%hubbard_j_impurity(i, l)
               end do
            end if
         end do
      end if

      ! Optional self-consistent U selector (0/1 mask by atom type and l-channel).
      this%hubbard_u_sc_check = .false.
      this%hubbard_u_sc(:, :) = 0
      if (has_hubbard_u_sc) then
         do i = 1, min(size(this%hubbard_u_sc, 1), size(hubbard_u_sc, 1))
            do l = 1, min(size(this%hubbard_u_sc, 2), size(hubbard_u_sc, 2))
               if (hubbard_u_sc(i, l) == hubbard_sc_sentinel) cycle
               if (hubbard_u_sc(i, l) < 0 .or. hubbard_u_sc(i, l) > 1) then
                  call g_logger%fatal('Invalid hubbard_u_sc value. Only 0 or 1 is allowed.', __FILE__, __LINE__)
               end if
               this%hubbard_u_sc(i, l) = hubbard_u_sc(i, l)
               if (this%hubbard_u_sc(i, l) == 1) this%hubbard_u_sc_check = .true.
            end do
         end do
      end if

      this%hubbard_v_check = .false.
      if (has_hubbard_v) then
         this%hubbard_v(:, :, :, :) = 0.0_rp
         do i = 1, min(size(this%hubbard_v, 1), size(hubbard_v, 1))
            do j = 1, min(size(this%hubbard_v, 2), size(hubbard_v, 2))
               do li = 1, min(size(this%hubbard_v, 3), size(hubbard_v, 3))
                  do lj = 1, min(size(this%hubbard_v, 4), size(hubbard_v, 4))
                     if (hubbard_v(i, j, li, lj) == hubbard_sentinel) cycle
                     this%hubbard_v(i, j, li, lj) = hubbard_v(i, j, li, lj)/ry2ev
                  end do
               end do
            end do
         end do

      ! Enforce V_ji(lj,li) = V_ij(li,lj) when only one side is provided.
      ! this%hubbard_v is already in Ry at this point.
      do i = 1, size(this%hubbard_v, 1)
            do j = 1, size(this%hubbard_v, 2)
               do li = 1, size(this%hubbard_v, 3)
                  do lj = 1, size(this%hubbard_v, 4)
                     if (abs(this%hubbard_v(i, j, li, lj)) > 1.0e-10_rp) then
                        this%hubbard_v_check = .true.
                        if (j <= size(this%hubbard_v, 1) .and. i <= size(this%hubbard_v, 2)) then
                           if (abs(this%hubbard_v(j, i, lj, li)) <= 1.0e-10_rp) then
                              this%hubbard_v(j, i, lj, li) = this%hubbard_v(i, j, li, lj)
                           end if
                        end if
                     end if
                  end do
               end do
            end do
         end do
      end if

      this%hubbard_u_general_check = .false.
      do i = 1, this%lattice%ntype
         if (maxval(abs(this%lattice%symbolic_atoms(i)%potential%hubbard_u(:))) > 1.0e-10_rp) then
            this%hubbard_u_general_check = .true.
         end if
         if (maxval(abs(this%lattice%symbolic_atoms(i)%potential%hubbard_j(:))) > 1.0e-10_rp) then
            this%hubbard_u_general_check = .true.
         end if
      end do

      ! Self-consistent U and fixed U/J are mutually exclusive in this implementation.
      if (this%hubbard_u_sc_check .and. this%hubbard_u_general_check) then
         call g_logger%fatal('Both hubbard_u_sc and explicit hubbard_u/hubbard_j are set. Use only one mode.', __FILE__, __LINE__)
      end if

      if (rank == 0) then
         if (this%hubbard_u_general_check .or. this%hubbard_u_sc_check .or. this%hubbard_v_check) then
         call g_logger%info('HUBBARD summary: form='//trim(this%hubbard_u_potential_form)// &
                            ' fixed_UJ='//merge('T', 'F', this%hubbard_u_general_check)// &
                            ' sc_U='//merge('T', 'F', this%hubbard_u_sc_check)// &
                            ' V='//merge('T', 'F', this%hubbard_v_check), __FILE__, __LINE__)
         end if
      end if
      if (this%hubbard_u_general_check) then
         block
            integer :: lch, nch
            character(len=512) :: u_msg, j_msg
            do i = 1, this%lattice%ntype
               u_msg = ''
               j_msg = ''
               nch = size(this%lattice%symbolic_atoms(i)%potential%hubbard_u)
               do lch = 1, nch
                  u_msg = trim(u_msg)//' '//fmt('f10.6', this%lattice%symbolic_atoms(i)%potential%hubbard_u(lch))
               end do
               nch = size(this%lattice%symbolic_atoms(i)%potential%hubbard_j)
               do lch = 1, nch
                  j_msg = trim(j_msg)//' '//fmt('f10.6', this%lattice%symbolic_atoms(i)%potential%hubbard_j(lch))
               end do
               if (rank == 0) then
                  call g_logger%info('HUBBARD fixed U/J type='//fmt('i4', i)//' [Ry] U='//trim(u_msg)//' J='//trim(j_msg), __FILE__, __LINE__)
               end if
            end do
         end block
      end if
      if (this%hubbard_u_sc_check) then
         do i = 1, size(this%hubbard_u_sc, 1)
            if (rank == 0) then
               call g_logger%info('HUBBARD sc_U mask type='//fmt('i4', i)//' [s p d f]='// &
                                  fmt('i2', this%hubbard_u_sc(i, 1))//' '//fmt('i2', this%hubbard_u_sc(i, 2))//' '// &
                                  fmt('i2', this%hubbard_u_sc(i, 3))//' '//fmt('i2', this%hubbard_u_sc(i, 4)), __FILE__, __LINE__)
            end if
         end do
      end if

   end subroutine build_from_file

   !> @brief Reset Hamiltonian flags, arrays, and cached bounds to defaults.
   !> @details Restores namelist-controlled options to their baseline values and
   !>          clears allocatable storage so a constructor can rebuild the object
   !>          for the current lattice/control state.
   !> @param[inout] this Hamiltonian object to reset.
   module subroutine restore_to_default(this)
      class(hamiltonian) :: this

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('hamiltonian.lsham', this%lsham, (/nb, nb, this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.tmat', this%tmat, (/nb, nb, 3, this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.hhmag', this%hhmag, (/norb, norb, 4/))
      call g_safe_alloc%allocate('hamiltonian.hmag', this%hmag, (/norb, norb, this%charge%lattice%kk, 4/))
      call g_safe_alloc%allocate('hamiltonian.ee', this%ee, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.eecc', this%eecc, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.hall', this%hall, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax/))
      call g_safe_alloc%allocate('hamiltonian.hallcc', this%hallcc, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax/))
      call g_safe_alloc%allocate('hamiltonian.hall_glob', this%hall_glob, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax/))
      call g_safe_alloc%allocate('hamiltonian.hallcc_glob', this%hallcc_glob, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax/))
      call g_safe_alloc%allocate('hamiltonian.ee_glob', this%ee_glob, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.eecc_glob', this%eecc_glob, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.eeo', this%eeo, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.eeoee', this%eeoee, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.hallo', this%hallo, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax/))
      call g_safe_alloc%allocate('hamiltonian.obarm', this%obarm, (/nb, nb, this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.enim', this%enim, (/nb, nb, this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.hall_glob', this%hall_glob, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax/))
      call g_safe_alloc%allocate('hamiltonian.ee_glob', this%ee_glob, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.ee0_glob', this%eeo_glob, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.hallo_glob', this%hallo_glob, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax/))
      call g_safe_alloc%allocate('hamiltonian.enim_glob', this%enim_glob, (/nb, nb, this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.v_a', this%v_a, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.v_b', this%v_b, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.vo_a', this%vo_a, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.vo_b', this%vo_b, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.js_a', this%js_a, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.jl_a', this%jl_a, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.jso_a', this%jso_a, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.jlo_a', this%jlo_a, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.velocity_scale', this%velocity_scale, (/this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.hxc', this%hxc, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.hubbard_u_pot', this%hubbard_u_pot, (/nb, nb, this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.hubbard_u_impurity', this%hubbard_u_impurity, (/max(1, this%charge%lattice%nrec), 4/))
      call g_safe_alloc%allocate('hamiltonian.hubbard_j_impurity', this%hubbard_j_impurity, (/max(1, this%charge%lattice%nrec), 4/))
      call g_safe_alloc%allocate('hamiltonian.hubbard_u_sc', this%hubbard_u_sc, (/max(1, this%charge%lattice%ntype), 4/))
      ! Hubbard-V is indexed by symbolic atom type pair (itype,jtype), not by impurity index.
      call g_safe_alloc%allocate('hamiltonian.hubbard_v', this%hubbard_v, (/max(1, this%charge%lattice%ntype), max(1, this%charge%lattice%ntype), 4, 4/))
      call g_safe_alloc%allocate('hamiltonian.hubbard_v_pot', this%hubbard_v_pot, (/nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype/))
#else
      allocate (this%lsham(nb, nb, this%charge%lattice%ntype))
      allocate (this%tmat(nb, nb, 3, this%charge%lattice%ntype))
      allocate (this%hhmag(norb, norb, 4), this%hmag(norb, norb, this%charge%lattice%kk, 4))
      allocate (this%ee(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%eecc(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%hall(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%nmax))
      allocate (this%hallcc(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%nmax))
      allocate (this%hxc(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%eeo(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%eeoee(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%hallo(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%nmax))
      allocate (this%obarm(nb, nb, this%charge%lattice%ntype))
      allocate (this%enim(nb, nb, this%charge%lattice%ntype))
      allocate (this%ee_glob(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%eecc_glob(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%hall_glob(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%nmax))
      allocate (this%hallcc_glob(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%nmax))
      allocate (this%eeo_glob(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%hallo_glob(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%nmax))
      allocate (this%enim_glob(nb, nb, this%charge%lattice%ntype))
      ! Velocity operators
      allocate (this%v_a(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%v_b(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%vo_a(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%vo_b(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%js_a(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%jl_a(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%jso_a(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%jlo_a(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%velocity_scale(this%charge%lattice%ntype))
      allocate (this%hubbard_u_pot(nb, nb, this%charge%lattice%ntype))
      allocate (this%hubbard_u_impurity(max(1, this%charge%lattice%nrec), 4))
      allocate (this%hubbard_j_impurity(max(1, this%charge%lattice%nrec), 4))
      allocate (this%hubbard_u_sc(max(1, this%charge%lattice%ntype), 4))
      ! Hubbard-V is indexed by symbolic atom type pair (itype,jtype), not by impurity index.
      allocate (this%hubbard_v(max(1, this%charge%lattice%ntype), max(1, this%charge%lattice%ntype), 4, 4))
      allocate (this%hubbard_v_pot(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
#endif

      this%lsham(:, :, :) = 0.0d0
      this%tmat(:, :, :, :) = 0.0d0
      this%hhmag(:, :, :) = 0.0d0
      this%hmag(:, :, :, :) = 0.0d0
      this%hall(:, :, :, :) = 0.0d0
      this%hallcc(:, :, :, :) = 0.0d0
      this%hxc(:, :, :, :) = 0.0d0
      this%ee(:, :, :, :) = 0.0d0
      this%eecc(:, :, :, :) = 0.0d0
      this%hallo(:, :, :, :) = 0.0d0
      this%eeo(:, :, :, :) = 0.0d0
      this%eeoee(:, :, :, :) = 0.0d0
      this%obarm(:, :, :) = 0.0d0
      this%enim(:, :, :) = 0.0d0
      this%hall_glob(:, :, :, :) = 0.0d0
      this%hallcc_glob(:, :, :, :) = 0.0d0
      this%ee_glob(:, :, :, :) = 0.0d0
      this%eecc_glob(:, :, :, :) = 0.0d0
      this%hallo_glob(:, :, :, :) = 0.0d0
      this%eeo_glob(:, :, :, :) = 0.0d0
      this%enim_glob(:, :, :) = 0.0d0
      this%v_a(:, :, :, :) = 0.0d0
      this%v_b(:, :, :, :) = 0.0d0
      this%vo_a(:, :, :, :) = 0.0d0
      this%vo_b(:, :, :, :) = 0.0d0
      this%js_a(:, :, :, :) = 0.0d0
      this%jl_a(:, :, :, :) = 0.0d0
      this%jso_a(:, :, :, :) = 0.0d0
      this%jlo_a(:, :, :, :) = 0.0d0
      this%velocity_scale(:) = 1.0d0
      this%hubbard_u_pot(:, :, :) = 0.0d0
      this%bounds%algorithm = 'none'
      this%bounds%scaling = 1.05_rp
      this%bounds%e_min = 0.0_rp
      this%bounds%e_max = 0.0_rp
      this%hubbard_u_impurity(:, :) = 0.0d0
      this%hubbard_j_impurity(:, :) = 0.0d0
      this%hubbard_u_sc(:, :) = 0
      this%hubbard_v(:, :, :, :) = 0.0d0
      this%hubbard_v_pot(:, :, :, :) = 0.0d0
      this%hoh = .false.
      this%local_axis = .false.
      this%orb_pol = .false.
         this%ccor_2c = .false.
         this%ccor_elin = 0.0_rp
         this%ccor_vmt_mode = 'surface_scalar'
      this%ccor_debug = .false.
      this%ccor_strict = .false.
      this%v_alpha(:) = [1, 0, 0]
      this%v_beta(:) = [1, 0, 0]
      this%q_ss(:) = [0.0_rp, 0.0_rp, 0.0_rp]
      this%theta_ss = 0.0_rp
      if (allocated(this%theta_ss_sublattice)) deallocate (this%theta_ss_sublattice)
      if (allocated(this%phi_ss_sublattice)) deallocate (this%phi_ss_sublattice)
      this%js_alpha = 'z'
      this%jl_alpha = 'z'
      this%hubbard_u_potential_form = 'liechtenstein'
      this%hubbard_u_general_check = .false.
      this%hubbard_u_impurity_check = .false.
      this%hubbard_u_sc_check = .false.
      this%hubbard_v_check = .false.
      this%export = 'none'
   end subroutine restore_to_default

   !> @brief Build real-space orbital-current operator blocks.
   !> @details Forms orbital velocity operators from angular-momentum matrices and
   !>          real-space hopping geometry for orbital transport and orbital torque
   !>          workflows.
   !> @param[inout] this Hamiltonian object; fills jl_a/jl_b-style operator arrays.
   !> @note Uses the same site/neighbor layout as ee/hall velocity blocks.
   module subroutine build_realspace_orbital_velocity_operators(this)
      class(hamiltonian), intent(inout) :: this
   
      integer :: ntype, ia, nr, m
      integer :: hblocksize
      complex(rp), allocatable :: tmp1(:,:), tmp2(:,:), L_op(:,:)  ! Temp matrices
      complex(rp), dimension(norb, norb) :: mLx, mLy, mLz

      hblocksize = size(this%v_a, 1)
   
      ! Allocate local matrices for partial products
      allocate(tmp1(hblocksize,hblocksize), tmp2(hblocksize,hblocksize), L_op(hblocksize,hblocksize))
   
      ! Initialize the orbital–velocity array
      this%jl_a(:, :, :, :) = (0.0_rp, 0.0_rp)
      this%jlo_a(:, :, :, :) = (0.0_rp, 0.0_rp) 

      !  Getting the angular momentum operators from the math_mod that are in cartesian coordinates
      mLx(:, :) = L_x(:, :)
      mLy(:, :) = L_y(:, :)
      mLz(:, :) = L_z(:, :)
      
      ! Transforming them into the spherical harmonics coordinates
      call hcpx(mLx, 'cart2sph')
      call hcpx(mLy, 'cart2sph')
      call hcpx(mLz, 'cart2sph')
   
      ! Pick which orbital operator L_x, L_y, or L_z based on some user choice
      select case (this%jl_alpha)   ! or whichever variable holds 'x','y','z'
      case ('x')
         L_op(1:norb, 1:norb) = mLx(:, :)
         L_op(norb+1:nb, norb+1:nb) = mLx(:, :)
      case ('y')
         L_op(1:norb, 1:norb) = mLy(:, :)
         L_op(norb+1:nb, norb+1:nb) = mLy(:, :)
      case ('z')
         L_op(1:norb, 1:norb) = mLz(:, :)
         L_op(norb+1:nb, norb+1:nb) = mLz(:, :)
      end select
   
      ! Loop over each atom type
      do ntype = 1, this%charge%lattice%ntype
         ia = this%charge%lattice%atlist(ntype)
         nr = this%charge%lattice%nn(ia, 1)
   
         ! For each neighbor block
         do m = 2, nr
            ! tmp1 = L_op * v_a(:,:,m,ntype)
            tmp1 = matmul(L_op, this%v_a(:,:,m,ntype))
   
            ! tmp2 = v_a(:,:,m,ntype) * L_op
            tmp2 = matmul(this%v_a(:,:,m,ntype), L_op)
   
            ! jl_a(:,:,m,ntype) = 0.5 * ( tmp1 + tmp2 )
            this%jl_a(:,:,m,ntype) = 0.5_rp * ( tmp1 + tmp2 )
   
            if (this%hoh) then
               tmp1 = 0.0d0; tmp2 = 0.0d0
               ! tmp1 = js_a * v_a(:,:,m,ntype)
               tmp1 = matmul(L_op, this%vo_a(:, :, m, ntype))

               ! tmp2 = v_a(:,:,m,ntype) * js_a
               tmp2 = matmul(this%vo_a(:, :, m, ntype), L_op)

               ! v_sza(:,:,m,ntype) = 0.5 * ( tmp1 + tmp2 )
               this%jlo_a(:,:,m,ntype) = 0.5_rp * ( tmp1 + tmp2 )
            end if

            ! Optional debugging output:
            ! write(*,*) 'm=', m
            ! write(*,'(18f10.6)') real(this%jo_a(:,:,m,ntype))
            ! write(*,*)
            ! write(*,'(18f10.6)') aimag(this%jo_a(:,:,m,ntype))
         end do
      end do
   
      deallocate(tmp1, tmp2, L_op)
   end subroutine build_realspace_orbital_velocity_operators

   !> @brief Build real-space spin-current operator blocks.
   !> @details Combines spin matrices with the velocity-operator layout so
   !>          stochastic conductivity can evaluate spin-current correlations.
   !> @param[inout] this Hamiltonian object; fills js_a and related spin-current arrays.
   module subroutine build_realspace_spin_operators(this)
      class(hamiltonian), intent(inout) :: this
   
      integer :: ntype, ia, nr, m, ji, ja, atom_neighbor
      integer :: hblocksize
      complex(rp), allocatable :: tmp1(:, :), tmp2(:, :), S_op(:, :)  ! Temp matrices for partial products
   
      ! Derive dimension from your velocity array:
      hblocksize = size(this%v_a, 1)  ! e.g. first dimension of v_a
   
      ! Allocate temporary matrices for local block multiplication
      allocate(tmp1(hblocksize, hblocksize), tmp2(hblocksize, hblocksize), S_op(hblocksize, hblocksize))
   
      ! Initialize the spin–velocity array to zero
      this%js_a(:, :, :, :) = (0.0_rp, 0.0_rp)
      this%jso_a(:, :, :, :) = (0.0_rp, 0.0_rp)
      
      select case(this%js_alpha)
      
      case('z')
         S_op = S_z
      case('x')
         S_op = S_x
      case('y')
         S_op = S_y
      end select

      !write(*,'(18f10.6)') real(S_op)
      ! Loop over each atom type
      do ntype = 1, this%charge%lattice%ntype
         ia = this%charge%lattice%atlist(ntype)
         nr = this%charge%lattice%nn(ia, 1)
   
         ! For each neighbor block 
         do m = 2, nr

            tmp1 = 0.0d0; tmp2 = 0.0d0
            ! tmp1 = js_a * v_a(:,:,m,ntype)
            tmp1 = matmul(S_op, this%v_a(:, :, m, ntype))
   
            ! tmp2 = v_a(:,:,m,ntype) * js_a
            tmp2 = matmul(this%v_a(:, :, m, ntype), S_op)

            ! v_sza(:,:,m,ntype) = 0.5 * ( tmp1 + tmp2 )
            this%js_a(:,:,m,ntype) = 0.5_rp * ( tmp1 + tmp2 )
            !write(*,*) 'm=', m
            !write(*,'(18f10.6)') real(this%js_a(:,:,m,ntype))
            !write(*,*)
            !write(*,'(18f10.6)') aimag(this%js_a(:,:,m,ntype)) 
            if (this%hoh) then
               tmp1 = 0.0d0; tmp2 = 0.0d0
               ! tmp1 = js_a * v_a(:,:,m,ntype)
               tmp1 = matmul(S_op, this%vo_a(:, :, m, ntype))
   
               ! tmp2 = v_a(:,:,m,ntype) * js_a
               tmp2 = matmul(this%vo_a(:, :, m, ntype), S_op)
   
               ! v_sza(:,:,m,ntype) = 0.5 * ( tmp1 + tmp2 )
               this%jso_a(:,:,m,ntype) = 0.5_rp * ( tmp1 + tmp2 )
            end if

         end do  ! m
      end do  ! ntype
   
      deallocate(tmp1, tmp2)
   end subroutine build_realspace_spin_operators

   !> @brief Build real-space spin-torque operator blocks.
   !> @details Forms torque-current operators from spin operators and local
   !>          Hamiltonian/SOC blocks for spin-torque response calculations.
   !> @param[inout] this Hamiltonian object; fills js_a/jso_a torque arrays.
   module subroutine build_realspace_spin_torque_operators(this)
      class(hamiltonian), intent(inout) :: this

      integer :: ntype, ia, nr, m, ji, ja, atom_neighbor, ino
      integer :: hblocksize
      complex(rp), allocatable :: tmp1(:, :), tmp2(:, :), S_op(:, :)  ! Temp matrices for partial products
      complex(rp), dimension(nb, nb) :: locham

      ! Derive dimension from your velocity array:
      hblocksize = size(this%v_a, 1)  ! e.g. first dimension of v_a

      ! Allocate temporary matrices for local block multiplication
      allocate(tmp1(hblocksize, hblocksize), tmp2(hblocksize, hblocksize), S_op(hblocksize, hblocksize))

      ! Initialize the spin–velocity array to zero
      this%js_a(:, :, :, :) = (0.0_rp, 0.0_rp)
      this%jso_a(:, :, :, :) = (0.0_rp, 0.0_rp)

      select case(this%js_alpha)

      case('z')
         S_op = S_z
      case('x')
         S_op = S_x
      case('y')
         S_op = S_y
      end select

      !write(*,'(18f10.6)') real(S_op)
      ! Loop over each atom type
      do ntype = 1, this%charge%lattice%ntype
         ia = this%charge%lattice%atlist(ntype)
         nr = this%charge%lattice%nn(ia, 1)
          
         ! For each neighbor block 
         do m = 1, nr

            !if (m==1) then
            !  locham(:,:) = this%ee(:, :, m, ntype) + this%lsham(:, :, ntype) 
            !else
            !  locham(:,:) = this%ee(:, :, m, ntype)
            !end if

            locham(:,:) = this%hxc(:, :, m, ntype)

            tmp1 = 0.0d0; tmp2 = 0.0d0
            ! tmp1 = js_a * hxc(:,:,m,ntype)
            tmp1 = matmul(S_op, locham(:, :))

            ! tmp2 = hxc(:,:,m,ntype) * js_a
            tmp2 = matmul(locham(:, :), S_op)

            ! v_sza(:,:,m,ntype) = 0.5 * ( tmp1 - tmp2 )
            this%js_a(:,:,m,ntype) = (1 / i_unit) * ( tmp1 - tmp2 )
            !write(*,*) 'm=', m
            !write(*,'(18f10.6)') real(this%js_a(:,:,m,ntype))
            !write(*,*)
            !write(*,'(18f10.6)') aimag(this%js_a(:,:,m,ntype)) 
            if (this%hoh) then

               if (m==1) then
                 locham(:,:) = this%eeo(:, :, 1, ntype) + this%lsham(:, :, ntype)
               else
                 locham(:,:) = this%eeo(:, :, m, ntype)
               end if

               tmp1 = 0.0d0; tmp2 = 0.0d0
               ! tmp1 = js_a * eeo(:,:,m,ntype)
               tmp1 = matmul(S_op, locham(:, :))

               ! tmp2 = ee(:,:,m,ntype) * js_a
               tmp2 = matmul(locham(:, :), S_op)

               ! v_sza(:,:,m,ntype) = 0.5 * ( tmp1 + tmp2 )
               this%jso_a(:,:,m,ntype) = (1 / i_unit) * ( tmp1 - tmp2 )
            end if

         end do  ! m
      end do  ! ntype

      deallocate(tmp1, tmp2)
   end subroutine build_realspace_spin_torque_operators

   !> @brief Build real-space orbital-torque operator blocks.
   !> @details Forms orbital torque-current operators from angular-momentum
   !>          matrices and local Hamiltonian/SOC blocks for orbital response.
   !> @param[inout] this Hamiltonian object; fills jl_a/jlo_a torque arrays.
   module subroutine build_realspace_orbital_torque_operators(this)
      class(hamiltonian), intent(inout) :: this

      integer :: ntype, ia, nr, m, ji, ja, atom_neighbor, ino
      integer :: hblocksize
      complex(rp), allocatable :: tmp1(:, :), tmp2(:, :), L_op(:, :)  ! Temp matrices for partial products
      complex(rp), dimension(norb, norb) :: mLx, mLy, mLz
      complex(rp), dimension(nb, nb) :: locham

      ! Derive dimension from your velocity array:
      hblocksize = size(this%v_a, 1)  ! e.g. first dimension of v_a

      ! Allocate temporary matrices for local block multiplication
      allocate(tmp1(hblocksize, hblocksize), tmp2(hblocksize, hblocksize), L_op(hblocksize, hblocksize))

      ! Initialize the spin–velocity array to zero
      this%jl_a(:, :, :, :) = (0.0_rp, 0.0_rp)
      this%jlo_a(:, :, :, :) = (0.0_rp, 0.0_rp)

      !  Getting the angular momentum operators from the math_mod that are in cartesian coordinates
      mLx(:, :) = L_x(:, :)
      mLy(:, :) = L_y(:, :)
      mLz(:, :) = L_z(:, :)

      ! Transforming them into the spherical harmonics coordinates
      call hcpx(mLx, 'cart2sph')
      call hcpx(mLy, 'cart2sph')
      call hcpx(mLz, 'cart2sph')

      ! Pick which orbital operator L_x, L_y, or L_z based on some user choice
      select case (this%jl_alpha)   ! or whichever variable holds 'x','y','z'
      case ('x')
         L_op(1:norb, 1:norb) = mLx(:, :)
         L_op(norb+1:nb, norb+1:nb) = mLx(:, :)
      case ('y')
         L_op(1:norb, 1:norb) = mLy(:, :)
         L_op(norb+1:nb, norb+1:nb) = mLy(:, :)
      case ('z')
         L_op(1:norb, 1:norb) = mLz(:, :)
         L_op(norb+1:nb, norb+1:nb) = mLz(:, :)
      end select

      ! Loop over each atom type
      do ntype = 1, this%charge%lattice%ntype
         ia = this%charge%lattice%atlist(ntype)
         nr = this%charge%lattice%nn(ia, 1)

         ! For each neighbor block 
         do m = 1, nr

            if (m==1) then
              locham(:,:) = this%ee(:, :, m, ntype) + this%lsham(:, :, ntype)
            else
              locham(:,:) = this%ee(:, :, m, ntype)
            end if

            tmp1 = 0.0d0; tmp2 = 0.0d0
            ! tmp1 = jl_a * ee(:,:,m,ntype)
            tmp1 = matmul(L_op, locham(:, :))

            ! tmp2 = ee(:,:,m,ntype) * jl_a
            tmp2 = matmul(locham(:, :), L_op)

            ! v_lza(:,:,m,ntype) = 0.5 * ( tmp1 + tmp2 )
            this%jl_a(:,:,m,ntype) = (1 / i_unit) * ( tmp1 - tmp2 )
            if (this%hoh) then

               if (m==1) then
                 locham(:,:) = this%eeo(:, :, 1, ntype) + this%lsham(:, :, ntype)
               else
                 locham(:,:) = this%eeo(:, :, m, ntype)
               end if

               tmp1 = 0.0d0; tmp2 = 0.0d0
               ! tmp1 = jl_a * eeo(:,:,m,ntype)
               tmp1 = matmul(L_op, locham(:, :))

               ! tmp2 = ee(:,:,m,ntype) * jl_a
               tmp2 = matmul(locham(:, :), L_op)

               ! v_lza(:,:,m,ntype) = 0.5 * ( tmp1 + tmp2 )
               this%jlo_a(:,:,m,ntype) = (1 / i_unit) * ( tmp1 - tmp2 )
            end if

         end do  ! m
      end do  ! ntype

      deallocate(tmp1, tmp2)
   end subroutine build_realspace_orbital_torque_operators

   !> @brief Build charge velocity operators in real-space hopping layout.
   !> @details Projects bond displacement vectors onto v_alpha/v_beta directions
   !>          and weights hopping blocks to form v_a/v_b for Kubo conductivity.
   !> @param[inout] this Hamiltonian object; fills v_a, v_b, vo_a, and vo_b.
   !> @note Honors velocity_scale and the HOH companion-operator layout.
   module subroutine build_realspace_velocity_operators(this)
      ! Arguments
      class(hamiltonian), intent(inout) :: this
   
      ! Local variables
      integer :: ia, ntype, nr, m, i, j, velotype, ja, ji    ! Atom and neighbor indices
      integer :: atom_neighbor                               ! Neighbor atom index
      real(rp) :: veloscale
      real(rp), dimension(3) :: rij                          ! Displacement vector (x, y, z components)
      real(rp), dimension(3) :: dir_a, dir_b                 ! Velocity operator directions
      real(rp) :: norm_a, norm_b, dot_a, dot_b
      ! Initialize velocity operators to zero
      this%v_a(:, :, :, :) = 0.0_rp
      this%v_b(:, :, :, :) = 0.0_rp
   
      norm_a = norm2(this%v_alpha)
      norm_b = norm2(this%v_beta)

      dir_a(:) = this%v_alpha(:) / norm_a
      dir_b(:) = this%v_beta(:) / norm_b

      ! Loop over atom types
      do ntype = 1, this%charge%lattice%ntype

         ia = this%charge%lattice%atlist(ntype)  ! Atom number in the cluster
         nr = this%charge%lattice%nn(ia, 1)     ! Number of neighbors for this atom type
    
         ! Loop over neighbors
         do m = 2, nr   ! Start from 2 to exclude the onsite term

            atom_neighbor = this%charge%lattice%nn(ia, m)  ! Neighbor atom number
            if (atom_neighbor /= 0) then
               ! Compute displacement vector rij = r_i - r_j
               rij(:) = (this%charge%lattice%cr(:, ia) - this%charge%lattice%cr(:, atom_neighbor)) * this%charge%lattice%alat
   
               dot_a = dot_product(dir_a, rij); dot_b = dot_product(dir_b, rij)

               ! Compute velocity operator blocks
               this%v_a(:, :, m, ntype) = ((1 / i_unit) * dot_a * this%ee(:, :, m, ntype)) 
            
               velotype = this%charge%lattice%iz(atom_neighbor)
               veloscale = max(this%velocity_scale(ntype), this%velocity_scale(velotype))
               write(*,*) ntype, velotype, veloscale
               this%v_b(:, :, m, ntype) = ((1 / i_unit) * dot_b * this%ee(:, :, m, ntype)) * veloscale 
               ! If hoh is true, multiply the velocity operator by the overlap matrix, similarly to whats done to the Hamiltonian
               if (this%hoh) then
                  ji = this%charge%lattice%iz(atom_neighbor) 
                  call zgemm('n', 'n', nb, nb, nb, cone, this%v_a(:, :, m, ntype), nb, this%obarm(:, :, ji), nb, czero, this%vo_a(:, :, m, ntype), nb)
                  call zgemm('n', 'n', nb, nb, nb, cone, this%v_b(:, :, m, ntype), nb, this%obarm(:, :, ji), nb, czero, this%vo_b(:, :, m, ntype), nb)
               end if
            end if
         end do
      end do
   end subroutine build_realspace_velocity_operators

   !> @brief Build onsite spin-orbit coupling Hamiltonian blocks.
   !> @details Converts angular-momentum operators to the LMTO basis and combines
   !>          them with per-type SOC strengths to populate lsham for real-space
   !>          and reciprocal Hamiltonian construction.
   !> @param[inout] this Hamiltonian object; fills lsham(:,:,itype).
   module subroutine build_lsham(this)
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      integer :: i, j, k
      complex(rp) :: prefac, sg
      real(rp) :: soc_p, soc_d
      complex(rp), dimension(2) :: rac
      complex(rp), dimension(norb, norb) :: Lx, Ly, Lz
      real(rp) :: lz_loc
      !  Getting the angular momentum operators from the math_mod that are in cartesian coordinates
      Lx(:, :) = L_x(:, :)
      Ly(:, :) = L_y(:, :)
      Lz(:, :) = L_z(:, :)

      ! Transforming them into the spherical harmonics coordinates
      call hcpx(Lx, 'cart2sph')
      call hcpx(Ly, 'cart2sph')
      call hcpx(Lz, 'cart2sph')

      ! Writing the L.S hamiltonian
      this%lsham(:, :, :) = cmplx(0.0d0, 0.0d0)
      do k = 1, this%charge%lattice%ntype
         sg = cmplx(0.5d0, 0.0d0)
         soc_p = sqrt(this%charge%lattice%symbolic_atoms(k)%potential%xi_p(1)*this%charge%lattice%symbolic_atoms(k)%potential%xi_p(2))
         soc_d = sqrt(this%charge%lattice%symbolic_atoms(k)%potential%xi_d(1)*this%charge%lattice%symbolic_atoms(k)%potential%xi_d(2))
         ! For f-orbitals, use a scaling based on d if not explicitly available
         ! Set soc_f to small value (f-orbital spin-orbit is typically weak)
         ! real soc_f = 0.0_rp  ! f-orbital s-o coupling (for future enhancement)

         ! Check if orbital polarization is enabled
         if (this%orb_pol) then
            rac = sqrt(this%charge%lattice%symbolic_atoms(k)%potential%xi_d(1)*this%charge%lattice%symbolic_atoms(k)%potential%rac)
            lz_loc = sqrt(this%charge%lattice%symbolic_atoms(k)%potential%xi_d(1)*this%charge%lattice%symbolic_atoms(k)%potential%lmom(3))
         else
            rac = 0.0_rp
            lz_loc = 0.0_rp
         end if

         prefac = 0.0_rp
         do i = 1, norb
            do j = 1, norb
               ! p-orbitals (indices 2-4)
               if (i >= 2 .and. i <= 4 .and. j >= 2 .and. j <= 4) prefac = sg*soc_p
               ! d-orbitals (indices 5-9)
               if (i >= 5 .and. i <= 9 .and. j >= 5 .and. j <= 9) prefac = sg*soc_d
               ! f-orbitals (indices 10-16) - currently set to zero (no f-orbital s.o. coupling yet)
               if (i >= 10 .and. i <= 16 .and. j >= 10 .and. j <= 16) prefac = cmplx(0.0_rp, 0.0_rp)

               this%lsham(j, i, k) = this%lsham(j, i, k) + prefac*Lz(j, i) + Lz(j, i)*rac(1)*lz_loc ! H11
               this%lsham(j, i +spin_off, k) = this%lsham(j, i +spin_off, k) + prefac*(Lx(j, i) - i_unit*Ly(j, i)) ! H12
               this%lsham(j +spin_off, i, k) = this%lsham(j +spin_off, i, k) + prefac*(Lx(j, i) + i_unit*Ly(j, i)) ! H21
               this%lsham(j +spin_off, i +spin_off, k) = this%lsham(j +spin_off, i +spin_off, k) - prefac*Lz(j, i) - Lz(j, i)*rac(2)*lz_loc ! H22
            end do
         end do

         ! Debug output: Print on-site Hamiltonian for lmax=3
         if (norb == 16) then
            open(unit=999, file='debug_hamiltonian_lsham.txt', action='write', status='replace')
            write(999, '(A)') 'On-site Hamiltonian (lsham) for lmax=3 (SPDF basis)'
            write(999, '(A, I0)') 'Atom type: ', k
            write(999, '(A)') 'Real part:'
            do i = 1, norb+spin_off
               write(999, '(16F12.6)') (real(this%lsham(i, j, k)), j=1, norb+spin_off)
            end do
            write(999, '(A)') 'Imaginary part:'
            do i = 1, norb+spin_off
               write(999, '(16F12.6)') (aimag(this%lsham(i, j, k)), j=1, norb+spin_off)
            end do
            close(999)
         end if
      end do
   end subroutine build_lsham

   !> @brief Build the collinear spin-orbit torque operator.
   !> @details Evaluates the commutator-style torque operator T=[o,H_so] in the
   !>          spin-orbital basis for workflows that need magnetic torques.
   !> @param[inout] this Hamiltonian object; fills tmat.
   module subroutine torque_operator_collinear(this)
      !
      class(hamiltonian), intent(inout) :: this
      !
      ! Local variables
      integer :: i, j, k
      complex(rp) :: prefac, sg, soc_p, soc_d
      complex(rp), dimension(norb, norb) :: Lx, Ly, Lz
      !  Getting the angular momentum operators from the math_mod that are in cartesian coordinates
      Lx(:, :) = L_x(:, :)
      Ly(:, :) = L_y(:, :)
      Lz(:, :) = L_z(:, :)

      ! Transforming them into the spherical harmonics coordinates
      call hcpx(Lx, 'cart2sph')
      call hcpx(Ly, 'cart2sph')
      call hcpx(Lz, 'cart2sph')

      ! Now write the torque operator matrix
      this%tmat(:, :, :, :) = cmplx(0.0_rp, 0.0_rp)
      do k = 1, this%charge%lattice%ntype
         sg = cmplx(0.5_rp, 0.0_rp)
         soc_p = sqrt(this%charge%lattice%symbolic_atoms(k)%potential%xi_p(1)*this%charge%lattice%symbolic_atoms(k)%potential%xi_p(2))
         soc_d = sqrt(this%charge%lattice%symbolic_atoms(k)%potential%xi_d(1)*this%charge%lattice%symbolic_atoms(k)%potential%xi_d(2))
         prefac = 0.0_rp
         do i = 1, norb
            do j = 1, norb
               if (i >= 2 .and. i <= 4 .and. j >= 2 .and. j <= 4) prefac = sg*soc_p
               if (i >= 5 .and. i <= 9 .and. j >= 5 .and. j <= 9) prefac = sg*soc_d
               ! build Tx
               this%tmat(j, i, 1, k) = this%tmat(j, i, 1, k) + prefac*i_unit*Ly(j, i)*2.0_rp ! Tx_11
               this%tmat(j, i +spin_off, 1, k) = this%tmat(j, i +spin_off, 1, k) - prefac*Lz(j, i)*2.0_rp*cone ! Tx_12
               this%tmat(j +spin_off, i, 1, k) = this%tmat(j +spin_off, i, 1, k) + prefac*Lz(j, i)*2.0_rp*cone ! Tx_21
               this%tmat(j +spin_off, i +spin_off, 1, k) = this%tmat(j +spin_off, i +spin_off, 1, k) - prefac*i_unit*Ly(j, i)*2.0_rp ! Tx_22
               ! build Ty
               this%tmat(j, i, 2, k) = this%tmat(j, i, 2, k) - prefac*i_unit*Lx(j, i)*2.0_rp ! Ty_11
               this%tmat(j, i +spin_off, 2, k) = this%tmat(j, i +spin_off, 2, k) + prefac*i_unit*Lz(j, i)*2.0_rp ! Ty_12
               this%tmat(j +spin_off, i, 2, k) = this%tmat(j +spin_off, i, 2, k) + prefac*i_unit*Lz(j, i)*2.0_rp ! Ty_21
               this%tmat(j +spin_off, i +spin_off, 2, k) = this%tmat(j +spin_off, i +spin_off, 2, k) + prefac*i_unit*Lx(j, i)*2.0_rp ! Ty_22
               ! build Tz
               this%tmat(j, i +spin_off, 3, k) = this%tmat(j, i +spin_off, 3, k) + prefac*(Lx(j, i) - i_unit*Ly(j, i))*2.0_rp*cone ! Tz_12
               this%tmat(j +spin_off, i, 3, k) = this%tmat(j +spin_off, i, 3, k) + prefac*(Lx(j, i) + i_unit*Ly(j, i))*2.0_rp*cmone ! Tz_21
            end do
         end do
      end do

   end subroutine torque_operator_collinear

   !> @brief Build overlap-bar matrices for the orthogonal representation.
   !> @details Converts potential/overlap information from symbolic atoms into
   !>          onsite obarm blocks used by HOH and representation transforms.
   !> @param[inout] this Hamiltonian object; fills obarm(:,:,itype).
   module subroutine build_obarm(this)
      implicit none
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      complex(rp), dimension(norb, norb) :: obm0, obm1
      complex(rp), dimension(3) :: mom
      integer :: ntype ! Atom type index
      integer :: l, m ! Orbital index

      this%obarm = 0.d00

      do ntype = 1, this%lattice%ntype
         obm0 = cmplx(0.0d0); obm1 = cmplx(0.0d0)
         do m = 1, norb
            obm0(m, m) = this%lattice%symbolic_atoms(ntype)%potential%obx0(m)
            obm1(m, m) = this%lattice%symbolic_atoms(ntype)%potential%obx1(m)
         end do
         mom(:) = cmplx(this%lattice%symbolic_atoms(ntype)%potential%mom(:), 0.0d0)
         do m = 1, norb
            do l = 1, norb
               this%obarm(m, l, ntype) = obm0(m, l) + obm1(m, l)*mom(3)
               this%obarm(m +spin_off, l +spin_off, ntype) = obm0(m, l) - obm1(m, l)*mom(3)
               this%obarm(l, m +spin_off, ntype) = obm1(m, l)*mom(1) - i_unit*obm1(m, l)*mom(2)
               this%obarm(l +spin_off, m, ntype) = obm1(m, l)*mom(1) + i_unit*obm1(m, l)*mom(2)
            end do
         end do
         call hcpx(this%obarm(1:norb, 1:norb, ntype), 'cart2sph')
         call hcpx(this%obarm(norb+1:nb, norb+1:nb, ntype), 'cart2sph')
         call hcpx(this%obarm(1:norb, norb+1:nb, ntype), 'cart2sph')
         call hcpx(this%obarm(norb+1:nb, 1:norb, ntype), 'cart2sph')
      end do
   end subroutine build_obarm

   !> @brief Build onsite e_nu center matrices.
   !> @details Assembles spin-resolved gravity-center/e_nu blocks from symbolic
   !>          atom potentials for the orthogonalized Hamiltonian correction.
   !> @param[inout] this Hamiltonian object; fills enim(:,:,itype).
   module subroutine build_enim(this)
      implicit none
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      complex(rp), dimension(norb, norb) :: em0, em1
      complex(rp), dimension(norb) :: ex0, ex1
      complex(rp), dimension(3) :: mom
      complex(rp) :: eu, ed
      integer :: ntype ! Atom type index
      integer :: l, m ! Orbital index

      this%enim = 0.0d0

      do ntype = 1, this%lattice%ntype
         em0 = cmplx(0.0d0); em1 = cmplx(0.0d0)
         do m = 1, norb
            eu = this%lattice%symbolic_atoms(ntype)%potential%cx(m, 1) - this%lattice%symbolic_atoms(ntype)%potential%cex(m, 1)
            ed = this%lattice%symbolic_atoms(ntype)%potential%cx(m, 2) - this%lattice%symbolic_atoms(ntype)%potential%cex(m, 2)
            ex0(m) = 0.5*(eu + ed)
            ex1(m) = 0.5*(eu - ed)
            em0(m, m) = ex0(m)
            em1(m, m) = ex1(m)
         end do
         mom(:) = cmplx(this%lattice%symbolic_atoms(ntype)%potential%mom(:), 0.0d0)
         do m = 1, norb
            do l = 1, norb
               this%enim(m, l, ntype) = em0(m, l) + em1(m, l)*mom(3)
               this%enim(m +spin_off, l +spin_off, ntype) = em0(m, l) - em1(m, l)*mom(3)
               this%enim(l, m +spin_off, ntype) = em1(m, l)*mom(1) - i_unit*em1(m, l)*mom(2)
               this%enim(l +spin_off, m, ntype) = em1(m, l)*mom(1) + i_unit*em1(m, l)*mom(2)
            end do
         end do
         call hcpx(this%enim(1:norb, 1:norb, ntype), 'cart2sph')
         call hcpx(this%enim(norb+1:nb, norb+1:nb, ntype), 'cart2sph')
         call hcpx(this%enim(1:norb, norb+1:nb, ntype), 'cart2sph')
         call hcpx(this%enim(norb+1:nb, 1:norb, ntype), 'cart2sph')

         if (this%local_axis) then
            this%enim_glob = this%enim
         end if
      end do
   end subroutine build_enim

   !> @brief Build bulk real-space Hamiltonian hopping blocks.
   !> @details Assembles ee, eeo, and related per-type neighbor blocks from the
   !>          lattice structure constants and symbolic-atom potentials. Real-space
   !>          recursion and reciprocal Fourier transforms consume these blocks.
   !> @param[inout] this Hamiltonian object; fills bulk hopping arrays.
   !> @note Also prepares optional Hubbard/CCOR contributions when enabled.
   module subroutine build_bulkham(this)
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m, n, itype, ino, ja, jo, ji, nr, ia
      integer :: ntype

      if (this%hubbard_u_general_check) then
         if (rank == 0) call g_logger%info('HUBBARD applying on-site +U correction to bulk Hamiltonian', __FILE__, __LINE__)
         call this%calculate_hubbard_u_potential_general()
      end if
      if (this%hubbard_v_check) then
         if (rank == 0) call g_logger%info('HUBBARD applying inter-site +V correction to bulk Hamiltonian', __FILE__, __LINE__)
         call this%calculate_hubbard_v_potential()
      end if

      do ntype = 1, this%charge%lattice%ntype
         ia = this%charge%lattice%atlist(ntype) ! Atom number in clust
         ino = this%charge%lattice%num(ia) ! Atom bravais type of ia
         nr = this%charge%lattice%nn(ia, 1) ! Number of neighbours considered
         !write(123, *)´bulkham´
         ! call g_logger%info('Building Hamiltonian for atom type '//fmt('i5', ntype)//' with '//fmt('i5', nr)//' neighbours', __FILE__, __LINE__)
         call this%chbar_nc(ia, nr, ino, ntype)
         do m = 1, nr
            do i = 1, norb
               do j = 1, norb
                  this%ee(j, i, m, ntype) = this%hmag(j, i, m, 4) + this%hmag(j, i, m, 3)        ! H0+Hz
                  this%ee(j +spin_off, i +spin_off, m, ntype) = this%hmag(j, i, m, 4) - this%hmag(j, i, m, 3)        ! H0-Hz
                  this%ee(j, i +spin_off, m, ntype) = this%hmag(j, i, m, 1) - i_unit*this%hmag(j, i, m, 2) ! Hx-iHy
                  this%ee(j +spin_off, i, m, ntype) = this%hmag(j, i, m, 1) + i_unit*this%hmag(j, i, m, 2) ! Hx+iHy
                  ! Builds the magnetic part of the Hamiltonian only
                  this%hxc(j, i, m, ntype) = this%hmag(j, i, m, 3)        ! Hz
                  this%hxc(j +spin_off, i +spin_off, m, ntype) = - this%hmag(j, i, m, 3)        ! - Hz
                  this%hxc(j, i +spin_off, m, ntype) = this%hmag(j, i, m, 1) - i_unit*this%hmag(j, i, m, 2) ! Hx-iHy
                  this%hxc(j +spin_off, i, m, ntype) = this%hmag(j, i, m, 1) + i_unit*this%hmag(j, i, m, 2) ! Hx+iHy
               end do ! end of orbital j loop
            end do ! end of orbital i loop
            ! write(131,*) 'm=', m
            ! write(131,'(18f10.6)') real(this%ee(:,:,m,ntype))
            ! write(132,*) 'm=', m
            ! write(132,'(18f10.6)') aimag(this%ee(:,:,m,ntype))
         end do ! end of neighbour number
         if (this%hubbard_u_general_check) then
            do i = 1, nb
               do j = 1, nb
                  this%ee(i, j, 1, ntype) = this%ee(i, j, 1, ntype) + cmplx(this%hubbard_u_pot(i, j, ntype), 0.0_rp, kind=rp)
               end do
            end do
         end if
         if (this%hubbard_v_check) then
            do m = 1, nr
               do i = 1, nb
                  do j = 1, nb
                     this%ee(i, j, m, ntype) = this%ee(i, j, m, ntype) + cmplx(this%hubbard_v_pot(i, j, m, ntype), 0.0_rp, kind=rp)
                  end do
               end do
            end do
         end if
         if (this%hoh) then
            call this%build_obarm()
            call this%build_enim()
            do m = 1, nr
               ji = 0
               if (m > 1) then
                  ja = this%charge%lattice%nn(ia, m)
                  if (ja .ne. 0) then
                     ji = this%charge%lattice%iz(ja)
                  end if
               else
                  ji = this%charge%lattice%iz(ia)
               end if
               ! Check if neighbour ´m´ exists for atom ´ntype´, otherwise fill HoH Hamiltonian with zeros.
               if (ji > 0) then
                  call zgemm('n', 'n', nb, nb, nb, cone, this%ee(:, :, m, ntype), nb, this%obarm(:, :, ji), nb, czero, this%eeo(:, :, m, ntype), nb)
                  call zgemm('n', 'c', nb, nb, nb, cone, this%eeo(:, :, m, ntype), nb, this%ee(:, :, m, ntype), nb, czero, this%eeoee(:, :, m, ntype), nb)
               else
                  this%eeo(:, :, m, ntype) = 0.0d0
               end if
               !write(*,*) ´m=´, m
               !write(*,´(18f10.6)´) real(this%eeo(:,:,m,ntype))
               !write(*,*) ´ee´, m
               !write(*,´(18f10.6)´) real(this%ee(:,:,m,ntype))
            end do
         end if
      end do ! end of atom type number
      call this%build_ccor_bulk()
      if (this%local_axis) then
         this%ee_glob = this%ee
         this%eecc_glob = this%eecc
         if (this%hoh) this%eeo_glob = this%eeo
      end if
      if (trim(this%control%recur) == 'chebyshev') then
         call this%compute_hamiltonian_bounds(verbose=.false.)
      end if
      close (128)
   end subroutine build_bulkham

   !> @brief Build local-cluster Hamiltonian hopping blocks.
   !> @details Assembles hall/hallo blocks for impurity or surface local regions
   !>          where atom-specific local geometry replaces bulk type blocks.
   !> @param[inout] this Hamiltonian object; fills local hopping arrays.
   module subroutine build_locham(this)
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      integer :: it, ino, nr, nlim, m, i, j, ja, ji

      if (this%hubbard_u_general_check) then
         if (rank == 0) call g_logger%info('HUBBARD applying on-site +U correction to local Hamiltonian', __FILE__, __LINE__)
         call this%calculate_hubbard_u_potential_general()
      end if
      if (this%hubbard_v_check) then
         if (rank == 0) call g_logger%info('HUBBARD applying inter-site +V correction to local Hamiltonian', __FILE__, __LINE__)
         call this%calculate_hubbard_v_potential()
      end if

      call g_timer%start('build local hamiltonian')
    !!$omp parallel do private(nlim, nr, ino, m, i, j, ji, ja, this)
      do nlim = 1, this%charge%lattice%nmax
         nr = this%charge%lattice%nn(nlim, 1) ! Number of neighbours considered
         ino = this%charge%lattice%num(nlim)
         call this%chbar_nc(nlim, nr, ino, nlim)
         do m = 1, nr
            do i = 1, norb
               do j = 1, norb
                  this%hall(j, i, m, nlim) = this%hmag(j, i, m, 4) + this%hmag(j, i, m, 3) ! H0+Hz
                  this%hall(j +spin_off, i +spin_off, m, nlim) = this%hmag(j, i, m, 4) - this%hmag(j, i, m, 3) ! H0-Hz
                  this%hall(j, i +spin_off, m, nlim) = this%hmag(j, i, m, 1) - i_unit*this%hmag(j, i, m, 2) ! Hx-iHy
                  this%hall(j +spin_off, i, m, nlim) = this%hmag(j, i, m, 1) + i_unit*this%hmag(j, i, m, 2) ! Hx+iHy
               end do
            end do
         end do
         if (this%hubbard_u_general_check) then
            do i = 1, nb
               do j = 1, nb
                  this%hall(i, j, 1, nlim) = this%hall(i, j, 1, nlim) + cmplx(this%hubbard_u_pot(i, j, ino), 0.0_rp, kind=rp)
               end do
            end do
         end if
         if (this%hubbard_v_check) then
            do m = 1, nr
               do i = 1, nb
                  do j = 1, nb
                     this%hall(i, j, m, nlim) = this%hall(i, j, m, nlim) + cmplx(this%hubbard_v_pot(i, j, m, ino), 0.0_rp, kind=rp)
                  end do
               end do
            end do
         end if
         if (this%hoh) then
            call this%build_obarm()
            call this%build_enim()
            do m = 1, nr
               ji = 0
               if (m > 1) then
                  ja = this%charge%lattice%nn(nlim, m)
                  if (ja .ne. 0) then
                     ji = this%charge%lattice%iz(ja)
                  end if
               else
                  ji = this%charge%lattice%iz(nlim)
               end if
               ! Check if neighbour ´m´ exists for atom ´nlim´, otherwise fill HoH Hamiltonian with zeros.
               if (ji > 0) then
                  call zgemm('n', 'n', nb, nb, nb, cone, this%hall(1, 1, m, nlim), nb, this%obarm(1, 1, ji), nb, czero, this%hallo(1, 1, m, nlim), nb)
               else
                  this%hallo(:, :, m, nlim) = 0.0d0
               end if
            end do
         end if
      end do
    !!$omp end parallel do
      call this%build_ccor_local()
      if (this%local_axis) then
         this%hall_glob = this%hall
         this%hallcc_glob = this%hallcc
         if (this%hoh) this%hallo_glob = this%hallo
      end if
      call g_timer%stop('build local hamiltonian')
   end subroutine build_locham

   !> @brief Build a noncollinear spin-orbital hopping block from orbital data.
   !> @details Lifts an orbital hopping/structure block hhh into the spin-orbital
   !>          basis using the local moments of atom sites ia and ja, including
   !>          spin-spiral phase handling where enabled.
   !> @param[inout] this Hamiltonian object containing magnetic moment state.
   !> @param[in] ia Site index of the source atom.
   !> @param[in] ja Site index of the neighbor atom.
   !> @param[in] it Atom type of ia.
   !> @param[in] jt Atom type of ja.
   !> @param[in] vet Bond vector between the sites.
   !> @param[in] hhh Orbital hopping/structure block.
   module subroutine ham0m_nc(this, ia, ja, it, jt, vet, hhh)
      class(hamiltonian), intent(inout) :: this
      ! Input
      integer, intent(in) :: ia, ja ! Atom sites i and j
      integer, intent(in) :: it, jt ! Type of atom i and j
      real(rp), dimension(3), intent(in) :: vet
      real(rp), dimension(norb, norb), intent(in) :: hhh
      ! Local Variables
      integer :: ilm, jlm, m
      real(rp), dimension(3) :: mom_ia, mom_ja
      complex(rp), dimension(3) :: cross
      complex(rp), dimension(norb, norb) :: hhhc
      complex(rp), dimension(3) :: momc_i, momc_j
      complex(rp) :: dot
      real(rp) :: vv

      this%hhmag(:, :, :) = 0.0d0

      vv = norm2(vet)

      ! General non-collinear pair block: H_ij = W(m_i) . hhh . W(m_j) with
      ! W(m) = wx0 + wx1 (sigma.m). The moment directions come from the atoms'
      ! self-consistent moments, or from per-sublattice cone angles when set.
      ! No spin-spiral (q_ss) phase is applied in this routine: for k-space GBT the
      ! spiral is a post-assembly twist in the reciprocal module; for real-space
      ! spirals the site moments are set explicitly (supercell). q_ss = 0 paths are
      ! therefore bit-identical to the collinear / noncollinear-FM code.
      mom_ia = this%charge%lattice%symbolic_atoms(it)%potential%mom(:)
      mom_ja = this%charge%lattice%symbolic_atoms(jt)%potential%mom(:)
      if (allocated(this%theta_ss_sublattice) .and. allocated(this%phi_ss_sublattice)) then
         if (it <= size(this%theta_ss_sublattice) .and. jt <= size(this%theta_ss_sublattice)) then
            mom_ia(1) = sin(this%theta_ss_sublattice(it))*cos(this%phi_ss_sublattice(it))
            mom_ia(2) = sin(this%theta_ss_sublattice(it))*sin(this%phi_ss_sublattice(it))
            mom_ia(3) = cos(this%theta_ss_sublattice(it))
            mom_ja(1) = sin(this%theta_ss_sublattice(jt))*cos(this%phi_ss_sublattice(jt))
            mom_ja(2) = sin(this%theta_ss_sublattice(jt))*sin(this%phi_ss_sublattice(jt))
            mom_ja(3) = cos(this%theta_ss_sublattice(jt))
         end if
      end if

      ! Real to complex. The two pair moments are passed as explicit local vectors
      ! (never through a type-indexed array): for an elemental crystal it == jt and
      ! a shared array would alias m_i and m_j onto the same slot (audit item E3).
      dot = cmplx(dot_product(mom_ia, mom_ja), kind=kind(0.0d0))
      momc_i = cmplx(mom_ia, kind=kind(0.0d0))
      momc_j = cmplx(mom_ja, kind=kind(0.0d0))
      cross = cmplx(cross_product(mom_ia, mom_ja), kind=kind(0.0d0))
      hhhc(:, :) = cmplx(hhh(:, :), kind=kind(0.0d0))

      if (size(this%charge%lattice%symbolic_atoms(it)%potential%wx0) < norb) then
         call g_logger%fatal('wx0(it) too small: type='//int2str(it)// &
                             ' size='//int2str(size(this%charge%lattice%symbolic_atoms(it)%potential%wx0))// &
                             ' norb='//int2str(norb)// &
                             ' lmax='//int2str(this%charge%lattice%symbolic_atoms(it)%potential%lmax), __FILE__, __LINE__)
      end if
      if (size(this%charge%lattice%symbolic_atoms(jt)%potential%wx0) < norb) then
         call g_logger%fatal('wx0(jt) too small: type='//int2str(jt)// &
                             ' size='//int2str(size(this%charge%lattice%symbolic_atoms(jt)%potential%wx0))// &
                             ' norb='//int2str(norb)// &
                             ' lmax='//int2str(this%charge%lattice%symbolic_atoms(jt)%potential%lmax), __FILE__, __LINE__)
      end if

      do ilm = 1, norb
         do jlm = 1, norb
            this%hhmag(ilm, jlm, 4) = &
               this%charge%lattice%symbolic_atoms(it)%potential%wx0(ilm)*hhhc(ilm, jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx0(jlm) + &
               this%charge%lattice%symbolic_atoms(it)%potential%wx1(ilm)*hhhc(ilm, jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx1(jlm)*dot
         end do
      end do


      if (vv <= 0.01d0) then
         do ilm = 1, norb
            if (this%hoh) then
               this%hhmag(ilm, ilm, 4) = this%hhmag(ilm, ilm, 4) + this%charge%lattice%symbolic_atoms(it)%potential%cex0(ilm)
            else
               this%hhmag(ilm, ilm, 4) = this%hhmag(ilm, ilm, 4) + this%charge%lattice%symbolic_atoms(it)%potential%cx0(ilm)
            end if
         end do
      end if

      do m = 1, 3
         do jlm = 1, norb
            do ilm = 1, norb
               this%hhmag(ilm, jlm, m) = &
                  (this%charge%lattice%symbolic_atoms(it)%potential%wx1(ilm)*hhhc(ilm, jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx0(jlm))*momc_i(m) + &
                  (this%charge%lattice%symbolic_atoms(it)%potential%wx0(ilm)*hhhc(ilm, jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx1(jlm))*momc_j(m) + &
                  i_unit*this%charge%lattice%symbolic_atoms(it)%potential%wx1(ilm)*hhhc(ilm, jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx1(jlm)*cross(m)
            end do
         end do
      end do

      if (vv > 0.01d0) return
      do m = 1, 3
         do ilm = 1, norb
            if (this%hoh) then
               this%hhmag(ilm, ilm, m) = this%hhmag(ilm, ilm, m) + this%charge%lattice%symbolic_atoms(it)%potential%cex1(ilm)*momc_i(m)
            else
               this%hhmag(ilm, ilm, m) = this%hhmag(ilm, ilm, m) + this%charge%lattice%symbolic_atoms(it)%potential%cx1(ilm)*momc_i(m)
            end if
         end do
      end do

   end subroutine ham0m_nc

   !> @brief Build noncollinear local hopping data for one cluster atom.
   !> @details Walks the neighbor list around atom ia, obtains orbital hopping
   !>          blocks, converts them with ham0m_nc, and stores the resulting local
   !>          Hamiltonian/field data.
   !> @param[inout] this Hamiltonian object; updates local noncollinear arrays.
   !> @param[in] ia Cluster atom index.
   !> @param[in] nr Number of neighbors considered.
   !> @param[in] ino Bravais/type index for ia.
   !> @param[in] ntype Atom type index.
   module subroutine chbar_nc(this, ia, nr, ino, ntype)
      class(hamiltonian), intent(inout) :: this
      ! Input
      integer, intent(in) :: ia ! Atom number in clust
      integer, intent(in) :: nr ! Number of neighbours considered
      integer, intent(in) :: ino ! Atom bravais type of ia
      integer, intent(in) :: ntype ! Atom type
      ! Local variables
      real(rp) :: r2
      real(rp), dimension(3, size(this%charge%lattice%cr(1, :))) :: cralat ! Clust position times the lattice constant
      real(rp), dimension(3) :: vet
      real(rp), dimension(norb, norb) :: hhh
      integer :: i, j, k, l, m, n, it, jt, jj, nn_max_loc
      integer :: ni, mdir
      integer :: kk ! Clust size number
      real(rp), dimension(:, :), allocatable :: ham_vec

      this%hmag(:, :, :, :) = 0.0d0

      r2 = this%charge%lattice%r2
      cralat(:, :) = this%charge%lattice%cr(:, :)*this%charge%lattice%alat
      kk = this%charge%lattice%kk

      allocate(ham_vec(3, nr))
      nn_max_loc = nr
      call this%charge%lattice%clusba(r2, cralat, ia, kk, kk, nn_max_loc, ham_vec)

      it = this%charge%lattice%iz(ia)
      do m = 1, nr
         jj = this%charge%lattice%nn(ia, m)
         if (m == 1) then
            jj = ia
         end if
         if (jj /= 0) then
            jt = this%charge%lattice%iz(jj)
            if (this%lattice%pbc) then
               call this%lattice%f_wrap_coord_diff(this%lattice%kk,this%lattice%cr*this%lattice%alat,ia,jj,vet)
            else
               vet(:) = (this%charge%lattice%cr(:, jj) - this%charge%lattice%cr(:, ia))*this%charge%lattice%alat
            end if
            call this%hmfind(vet, nr, hhh, m, ia, m, ni, ham_vec)
            if (ni == 0) then
               this%charge%lattice%nn(ia, m) = 0
            end if
            call this%ham0m_nc(ia, jj, it, jt, vet, hhh)
            do mdir = 1, 4
               call hcpx(this%hhmag(:, :, mdir), 'cart2sph')
               this%hmag(:, :, m, mdir) = this%hhmag(:, :, mdir)
            end do
         end if
      end do
      if (allocated(ham_vec)) deallocate(ham_vec)
   end subroutine chbar_nc

   !> @brief Find the orbital hopping block matching a neighbor vector.
   !> @details Searches the precomputed neighbor/hopping table for the vector vet
   !>          and returns the matching orbital block and neighbor index.
   !> @param[inout] this Hamiltonian object containing lattice tolerances.
   !> @param[in] vet Candidate neighbor vector.
   !> @param[in] nr Number of neighbors in the search list.
   !> @param[inout] hhh Orbital hopping block output.
   !> @param[in] m Neighbor counter being resolved.
   !> @param[in] ia Cluster atom index.
   !> @param[in] jn Neighbor-list index.
   !> @param[out] ni Matched neighbor index.
   !> @param[in] ham_vec Precomputed neighbor vectors.
   module subroutine hmfind(this, vet, nr, hhh, m, ia, jn, ni, ham_vec)
      class(hamiltonian), intent(inout) :: this
      ! Input
      integer, intent(in) :: m ! Number of the given neighbour
      integer, intent(in) :: ia ! Atom number in clust
      integer, intent(in) :: jn ! ?
      integer, intent(in) :: nr ! Number of neighbours
      real(rp), dimension(3), intent(in) :: vet
      ! Output
      integer, intent(out) :: ni
      real(rp), dimension(norb, norb), intent(inout) :: hhh
      real(rp), dimension(3, this%lattice%nn_max), intent(in) :: ham_vec
      ! Local variables
      real(rp) :: a1, a2, a3, aaa, eps
      integer :: i, ilm, jlm

      eps = 0.0001d0
      ni = 1
      a1 = 0.0d0
      a2 = 0.0d0
      a3 = 0.0d0
      aaa = 0.0d0
      do i = 1, nr
         !write(123, ´(a, i4, 3f10.4)´)´i´, i, this%charge%lattice%sbarvec(:, i)
         a1 = (vet(1) - ham_vec(1, i))
         a2 = (vet(2) - ham_vec(2, i))
         a3 = (vet(3) - ham_vec(3, i))
         aaa = a1**2 + a2**2 + a3**2
         if (aaa < eps) goto 1000
      end do
      write (*, '(1x, a, i4, a, i4, a, 3f10.6)') ' Error in hamiltonian%hmfind: Neighbour vector not found for atom', ia, &
         ' neighbour', jn, 'vector', vet(:)

      ni = 0
1000  continue
      do ilm = 1, norb
         do jlm = 1, norb
            hhh(ilm, jlm) = real(this%charge%lattice%sbar(jlm, ilm, m, this%charge%lattice%num(ia)))
         end do
      end do

      !do ilm = 1, norb
      !    write(123, ´(9f10.6)´)(hhh(ilm, jlm), jlm=1, 9)
      !end do
   end subroutine hmfind

   !> @brief Convert a global orbital index to site and local-orbital indices.
   !> @details Used by PAOFLOW import/export helpers to translate flat orbital
   !>          numbering into RS-LMTO site-major layout.
   !> @param[in] orb Global orbital index.
   !> @param[out] i_out Site index.
   !> @param[out] ia_out Local orbital index on the site.
   !> @param[in] n_atoms Number of atoms in the exported/imported cell.
   !> @param[in] max_orbital Number of orbitals per atom in the flat layout.
   module subroutine orb2site(orb, i_out, ia_out, n_atoms, max_orbital)
      integer, intent(in) :: orb, n_atoms, max_orbital
      integer, intent(out) :: i_out, ia_out

      if (orb <= n_atoms*max_orbital) then
         i_out = modulo(orb - 1, max_orbital) + 1
         ia_out = int((orb - 1)/max_orbital) + 1
      else
         i_out = modulo(orb - 1, max_orbital) + max_orbital + 1
         ia_out = int((orb - 1 - n_atoms*max_orbital)/max_orbital) + 1
      end if
   end subroutine orb2site

   !> @brief Convert site and local-orbital indices to a global orbital index.
   !> @details Used by PAOFLOW import/export helpers to write flat orbital
   !>          numbering from RS-LMTO site-major data.
   !> @param[in] i_in Site index.
   !> @param[in] ia_in Local orbital index on the site.
   !> @param[out] orb_out Global orbital index.
   !> @param[in] n_atoms Number of atoms in the exported/imported cell.
   !> @param[in] max_orbital Number of orbitals per atom in the flat layout.
   module subroutine site2orb(i_in, ia_in, orb_out, n_atoms, max_orbital)
      integer, intent(in) :: i_in, ia_in, n_atoms, max_orbital
      integer, intent(out) :: orb_out

      if (i_in <= max_orbital) then
         orb_out = (ia_in - 1)*max_orbital + i_in
      else
         orb_out = (ia_in - 1)*max_orbital + i_in + (n_atoms - 1)*max_orbital
      end if
   end subroutine site2orb

   !> @brief Rotate Hamiltonian blocks into a local magnetic axis.
   !> @details Saves global hopping/SOC data and rotates spin blocks so a local
   !>          moment direction can be treated as the quantization axis in legacy
   !>          recursion paths.
   !> @param[inout] this Hamiltonian object; updates rotated block arrays.
   !> @param[in] m_loc Local magnetic moment direction.
   !> @note Call rotate_from_local_axis to restore the global representation.
   module subroutine rotate_to_local_axis(this, m_loc)
      use math_mod, only: rotmag_loc
      class(hamiltonian), intent(inout) :: this
      real(rp), dimension(3), intent(in) :: m_loc

      ! Local variables
      integer :: sdim
      ! Rotate Hamiltonian to local axis if wanted
      if (this%local_axis) then
         sdim = product(shape(this%hall))/nb/nb
         call rotmag_loc(this%hall, this%hall_glob, sdim, m_loc)
         sdim = product(shape(this%ee))/nb/nb
         call rotmag_loc(this%ee, this%ee_glob, sdim, m_loc)
         if (this%ccor_2c) then
            sdim = product(shape(this%hallcc))/nb/nb
            call rotmag_loc(this%hallcc, this%hallcc_glob, sdim, m_loc)
            sdim = product(shape(this%eecc))/nb/nb
            call rotmag_loc(this%eecc, this%eecc_glob, sdim, m_loc)
         end if
         if (this%hoh) then
            sdim = product(shape(this%eeo))/nb/nb
            call rotmag_loc(this%eeo, this%eeo_glob, sdim, m_loc)
            sdim = product(shape(this%hallo))/nb/nb
            call rotmag_loc(this%hallo, this%hallo_glob, sdim, m_loc)
            sdim = product(shape(this%enim))/nb/nb
            call rotmag_loc(this%enim, this%enim_glob, sdim, m_loc)
         end if
      end if
   end subroutine rotate_to_local_axis

   !> @brief Restore Hamiltonian blocks from a local-axis rotation.
   !> @details Rotates local-axis Hamiltonian data back to the global spin frame
   !>          after an atom-specific recursion calculation.
   !> @param[inout] this Hamiltonian object; restores global block arrays.
   !> @param[in] m_loc Local magnetic moment direction used for the rotation.
   module subroutine rotate_from_local_axis(this, m_loc)
      use math_mod, only: rotmag_loc
      class(hamiltonian), intent(inout) :: this
      real(rp), dimension(3), intent(in) :: m_loc

      ! Local variables
      integer :: sdim
      ! Rotate Hamiltonian to local axis if wanted
      if (this%local_axis) then
         this%hall = this%hall_glob
         this%ee = this%ee_glob
         if (this%ccor_2c) then
            this%hallcc = this%hallcc_glob
            this%eecc = this%eecc_glob
         end if
         if (this%hoh) then
            this%eeo = this%eeo_glob
            this%hallo = this%hallo_glob
            this%enim = this%enim_glob
         end if
      end if
   end subroutine rotate_from_local_axis

   !> @brief Estimate spectral bounds for Chebyshev scaling.
   !> @details Computes or selects Hamiltonian energy bounds according to the
   !>          configured bounds algorithm and scaling factor. Chebyshev recursion
   !>          uses these bounds to map H into [-1,1].
   !> @param[inout] this Hamiltonian object; updates this%bounds.
   !> @param[in] verbose Optional flag enabling diagnostic logging.
   module subroutine compute_hamiltonian_bounds(this, verbose)
      class(hamiltonian), intent(inout) :: this
      logical, intent(in), optional :: verbose

      integer :: ntype, ia, nr, i, j, m, n_orb, n_sites
      integer :: isite, jsite, i_start, i_end, j_start, j_end, ineigh, ntype_i, ia_loc, ja
      real(rp) :: g_min, g_max, center, radius
      real(rp) :: hgamma_min, hgamma_max
      logical :: verb, have_gamma
      character(len=16) :: algo
      character(len=256) :: msg
      type(bounds) :: gamma_bounds
      complex(rp), allocatable :: h_gamma(:, :)

      verb = .false.
      if (present(verbose)) verb = verbose

      if (.not. allocated(this%ee)) then
         call g_logger%warning('compute_hamiltonian_bounds: ee is not allocated; skipping bounds', __FILE__, __LINE__)
         return
      end if

      call normalize_bounds_algorithm(this%bounds%algorithm, algo)
      if (algo == 'none') return

      g_min = huge(1.0_rp)
      g_max = -huge(1.0_rp)

      do ntype = 1, this%charge%lattice%ntype
         ia = this%charge%lattice%atlist(ntype)
         nr = this%charge%lattice%nn(ia, 1)
         do i = 1, nb
            center = real(this%ee(i, i, 1, ntype))
            if (this%ccor_2c) center = center + real(this%eecc(i, i, 1, ntype))
            if (allocated(this%lsham)) center = center + real(this%lsham(i, i, ntype))
            radius = 0.0_rp
            do m = 1, nr
               do j = 1, nb
                  if (m == 1 .and. i == j) cycle
                  radius = radius + abs(this%ee(i, j, m, ntype))
                  if (this%ccor_2c) radius = radius + abs(this%eecc(i, j, m, ntype))
               end do
            end do
            g_min = min(g_min, center - radius)
            g_max = max(g_max, center + radius)
         end do
      end do

      have_gamma = .false.
      hgamma_min = g_min
      hgamma_max = g_max
      n_orb = nb
      n_sites = this%charge%lattice%nrec
      if (n_sites > 0) then
         allocate(h_gamma(n_orb*n_sites, n_orb*n_sites))
         h_gamma(:, :) = cmplx(0.0_rp, 0.0_rp, kind=rp)
         do isite = 1, n_sites
            ntype_i = this%charge%lattice%ib(isite)
            ia_loc = this%charge%lattice%atlist(ntype_i)
            i_start = (isite - 1)*n_orb + 1
            i_end = isite*n_orb
            do ineigh = 1, this%charge%lattice%nn(ia_loc, 1)
               if (ineigh == 1) then
                  jsite = isite
               else
                  ja = this%charge%lattice%nn(ia_loc, ineigh)
                  jsite = this%charge%lattice%iz(ja)
                  if (jsite < 1 .or. jsite > n_sites) cycle
               end if
               j_start = (jsite - 1)*n_orb + 1
               j_end = jsite*n_orb
               h_gamma(i_start:i_end, j_start:j_end) = h_gamma(i_start:i_end, j_start:j_end) + this%ee(:, :, ineigh, ntype_i)
               if (this%ccor_2c) then
                  h_gamma(i_start:i_end, j_start:j_end) = h_gamma(i_start:i_end, j_start:j_end) + this%eecc(:, :, ineigh, ntype_i)
               end if
            end do
            if (allocated(this%lsham)) then
               h_gamma(i_start:i_end, i_start:i_end) = h_gamma(i_start:i_end, i_start:i_end) + this%lsham(:, :, ntype_i)
            end if
         end do
         call compute_spectrum_bounds(h_gamma, gamma_bounds, 'sturm', verbose=.false.)
         if (gamma_bounds%sturm_available) then
            have_gamma = .true.
            hgamma_min = gamma_bounds%e_min_sturm
            hgamma_max = gamma_bounds%e_max_sturm
         else
            have_gamma = .true.
            hgamma_min = gamma_bounds%e_min_gershgorin
            hgamma_max = gamma_bounds%e_max_gershgorin
         end if
         deallocate(h_gamma)
      end if

      call select_bounds_interval(algo, g_min, g_max, have_gamma, hgamma_min, hgamma_max, this%bounds%e_min, this%bounds%e_max)
      call apply_bounds_scaling(this%bounds%e_min, this%bounds%e_max, this%bounds%scaling)

      msg = 'Chebyshev bounds algorithm='//trim(algo)//' Gershgorin=['// &
            trim(fmt('f12.6', g_min))//','//trim(fmt('f12.6', g_max))//']'
      call g_logger%info(msg, __FILE__, __LINE__)
      if (have_gamma) then
         msg = 'Chebyshev bounds H(Gamma)=['//trim(fmt('f12.6', hgamma_min))//','// &
               trim(fmt('f12.6', hgamma_max))//']'
      else
         msg = 'Chebyshev bounds H(Gamma)=unavailable, fallback to Gershgorin'
      end if
      call g_logger%info(msg, __FILE__, __LINE__)
      msg = 'Chebyshev final scaled bounds=['//trim(fmt('f12.6', this%bounds%e_min))//','// &
            trim(fmt('f12.6', this%bounds%e_max))//'] scale='//trim(fmt('f8.4', this%bounds%scaling))
      call g_logger%info(msg, __FILE__, __LINE__)

      if (verb) then
         write(msg, '(A,L1)') 'Chebyshev bounds detailed mode, H(Gamma) available=', have_gamma
         call g_logger%info(trim(msg), __FILE__, __LINE__)
      end if
   end subroutine compute_hamiltonian_bounds

end submodule hamiltonian_build
