submodule(hamiltonian_mod) hamiltonian_effective_io
   use iso_fortran_env, only: int64, iostat_end
   use ieee_arithmetic, only: ieee_is_finite
   use hamiltonian_effective_mod, only: resolve_translation, hopping_checksum, publish_hamiltonian_file
   implicit none
   real(rp), parameter :: geometry_tol = 1.e-6_rp

   type :: effective_metadata
      logical :: checked_payload = .false.
      character(len=8) :: canonical_checksum = "", legacy_checksum = ""
      integer :: sites = 0, block_size = 0, nsp = 0, orbitals = 0
      integer(int64) :: records = -1, legacy_records = -1
      logical :: periodic = .false., boundary(3) = .false.
      real(rp) :: cell(3,3) = 0.0_rp, fermi = 0.0_rp, bounds(2) = 0.0_rp
      real(rp) :: mesh_fermi = 0.0_rp
      logical :: have_mesh_fermi = .false.
      logical :: had_hoh = .false., had_ccor = .false.
      real(rp) :: min_velocity_scale = huge(1.0_rp), max_velocity_scale = -huge(1.0_rp)
      complex(rp), allocatable :: onsite_soc(:,:,:)
      real(rp), allocatable :: position(:,:)
      integer, allocatable :: types(:), orbital_map(:)
      character(len=32), allocatable :: labels(:)
   end type effective_metadata

contains

   !> @brief Install a checked effective Hamiltonian before allocating recursion.
   !> @details Replace lattice neighbor slots and ee/hall together, retain the
   !>          actual endpoint image displacements, and return the stored E_F.
   module subroutine import_effective(this, fermi, mesh_fermi, external_fermi)
      class(hamiltonian), intent(inout) :: this
      real(rp), intent(out) :: fermi
      real(rp), optional, intent(out) :: mesh_fermi
      real(rp), optional, intent(in) :: external_fermi
      type(effective_metadata) :: meta
      type(effective_row), allocatable :: rows(:)
      class(lattice), pointer :: lat
      integer :: max_slots, i, j, k, t, idx, slot, expected, nperiod(3)
      integer, allocatable :: order(:), cell_index(:,:)
      type(effective_row) :: lattice_index
      real(rp) :: cell(3,3), inverse(3,3), delta(3), coord(3)
      integer :: r(3), key(4), target_cell(3)
      logical :: legacy, found, canonical_exists, external
      complex(rp) :: empty(nb,nb)

      lat => this%lattice
      external = this%effective_import == 'external'
      if (external) then
         if (.not. present(external_fermi)) &
            call g_logger%fatal('External TB import needs a Fermi energy from the input deck.', __FILE__, __LINE__)
         call external_metadata(lat, external_fermi, meta)
      else
         call read_metadata(meta)
      end if
      if (meta%block_size /= nb .or. meta%orbitals /= norb .or. meta%nsp /= this%control%nsp) &
         call g_logger%fatal('Effective Hamiltonian orbital or spin basis differs from input.', __FILE__, __LINE__)
      inquire(file='rs_tb.tb', exist=canonical_exists)
      legacy = external .or. this%effective_import == 'legacy' .or. (this%effective_import == 'auto' .and. .not. canonical_exists)
      if (legacy .and. nb /= 2*norb) &
         call g_logger%fatal('Effective seven-column import needs a two-spin-block basis.', __FILE__, __LINE__)
      cell = lat%a*lat%alat
      nperiod = [lat%n1,lat%n2,lat%n3]
      if (any(meta%boundary .neqv. (lat%pbc .and. [lat%b1,lat%b2,lat%b3]))) &
         call g_logger%fatal('Effective Hamiltonian boundary flags differ from input.', __FILE__, __LINE__)
      if (meta%periodic) then
         if (lat%ntype /= meta%sites .or. lat%nrec /= meta%sites .or. lat%nmax /= 0 .or. &
             .not. lat%pbc .or. any(nperiod <= 0)) &
            call g_logger%fatal('Compact effective import needs a matching bulk replica with pbc enabled.', __FILE__, __LINE__)
         if (any((.not. meta%boundary) .and. nperiod /= 1)) &
            call g_logger%fatal('Compact import needs one cell along each open direction; use cluster scope otherwise.', &
                                __FILE__, __LINE__)
      else
         if (lat%kk /= meta%sites) &
            call g_logger%fatal('Effective Hamiltonian cluster site count differs from input.', __FILE__, __LINE__)
         do i = 1, 3
            cell(:,i) = cell(:,i)*max(1,nperiod(i))
         end do
      end if
      if (maxval(abs(meta%cell-cell)) > geometry_tol) &
         call g_logger%fatal('Effective Hamiltonian cell differs from input.', __FILE__, __LINE__)

      do i = 1, meta%sites
         j = i
         if (meta%periodic) j = lat%atlist(i)
         t = lat%iz(j)
         if (t /= meta%types(i) .or. lat%symbolic_atoms(t)%element%symbol /= meta%labels(i) .or. &
             norm2(meta%position(:,i)-lat%cr(:,j)*lat%alat) > geometry_tol) &
            call g_logger%fatal('Effective Hamiltonian site map differs from input.', __FILE__, __LINE__)
      end do

      call read_records(meta, legacy, rows)
      if (external .and. rank == 0) call write_external_metadata(meta)
      this%lsham = czero
      if (allocated(meta%onsite_soc)) then
         if (size(meta%onsite_soc,3) /= lat%ntype) &
            call g_logger%fatal('Effective Hamiltonian SOC type count differs from input.', __FILE__, __LINE__)
         this%lsham = meta%onsite_soc
         this%effective_has_soc_split = .true.
      end if
      empty = czero
      do i = 1, meta%sites
         ! Slot 1 is always the zero-image onsite block, even if exactly zero.
         if (rows(i)%lookup([i,0,0,0]) == 0) call rows(i)%add([i,0,0,0],empty)
      end do
      max_slots = maxval(rows%count)

      if (meta%periodic) then
         if (abs(determinant(cell)) < spacing(1.0_rp)) &
            call g_logger%fatal('Singular cell in effective Hamiltonian import.', __FILE__, __LINE__)
         inverse = inverse_3x3(cell)
         allocate(cell_index(3,lat%kk))
         call lattice_index%reset(1)
         do i = 1, lat%kk
            t = lat%iz(i)
            delta = lat%cr(:,i)*lat%alat - meta%position(:,t)
            call resolve_translation(cell,inverse,delta,geometry_tol,r,found)
            if (.not. found) call g_logger%fatal('Cannot map replica site to effective cell.', __FILE__, __LINE__)
            cell_index(:,i) = merge(modulo(r,nperiod),r,meta%boundary)
            key = [t,cell_index(:,i)]
            if (lattice_index%lookup(key) /= 0) &
               call g_logger%fatal('Duplicate site in effective Hamiltonian periodic replica.', __FILE__, __LINE__)
            call lattice_index%add(key,empty(:1,:1))
         end do
      end if

      ! The generic kernels consume the same nn/ee/hall layout as native H.
      ! Keep image-resolved records separate even when the columns wrap to
      ! the same site. Only the position difference is overridden in Kubo.
      deallocate(lat%nn)
      lat%nn_max = max_slots
      allocate(lat%nn(lat%kk,max_slots), source=0)
      allocate(this%effective_bond(3,max_slots,lat%kk), source=0.0_rp)
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%deallocate('hamiltonian.ee',this%ee)
      call g_safe_alloc%deallocate('hamiltonian.hall',this%hall)
      if (meta%periodic) then
         call g_safe_alloc%allocate('hamiltonian.ee',this%ee,[nb,nb,max_slots,lat%ntype])
         call g_safe_alloc%allocate('hamiltonian.hall',this%hall,[nb,nb,max_slots,0])
      else
         call g_safe_alloc%allocate('hamiltonian.ee',this%ee,[nb,nb,0,lat%ntype])
         call g_safe_alloc%allocate('hamiltonian.hall',this%hall,[nb,nb,max_slots,lat%kk])
      end if
      this%ee = czero
      this%hall = czero
      if (.not. meta%periodic) lat%nmax = lat%kk
#else
      deallocate(this%ee,this%hall)
      if (meta%periodic) then
         allocate(this%ee(nb,nb,max_slots,lat%ntype), source=czero)
         allocate(this%hall(nb,nb,max_slots,0))
      else
         lat%nmax = lat%kk
         allocate(this%ee(nb,nb,0,lat%ntype))
         allocate(this%hall(nb,nb,max_slots,lat%kk), source=czero)
      end if
#endif
      do i = 1, lat%kk
         t = i
         if (meta%periodic) t = lat%iz(i)
         order = rows(t)%order()
         lat%nn(i,1) = rows(t)%count
         slot = 1
         do k = 1, rows(t)%count
            idx = order(k)
            j = rows(t)%key(1,idx)
            r = rows(t)%key(2:4,idx)
            expected = 1
            if (j /= t .or. any(r /= 0)) then
               slot = slot + 1
               expected = slot
            end if
            coord = meta%position(:,j) - meta%position(:,t) + matmul(cell,real(r,rp))
            if (meta%periodic) then
               target_cell = cell_index(:,i)+r
               key = [j,merge(modulo(target_cell,nperiod),target_cell,meta%boundary)]
               j = lattice_index%lookup(key)
               if (j == 0) call g_logger%fatal('Effective hopping target missing from periodic replica.', __FILE__, __LINE__)
               this%ee(:,:,expected,t) = rows(t)%block(:,:,idx)
            else
               this%hall(:,:,expected,i) = rows(t)%block(:,:,idx)
            end if
            if (expected == 1) then
               ! Recursion adds lsham onsite. The payload already includes it.
               if (meta%periodic) then
                  this%ee(:,:,1,t) = this%ee(:,:,1,t) - this%lsham(:,:,lat%iz(i))
               else
                  this%hall(:,:,1,i) = this%hall(:,:,1,i) - this%lsham(:,:,lat%iz(i))
               end if
               cycle
            end if
            lat%nn(i,expected) = j
            this%effective_bond(:,expected,i) = -coord
         end do
      end do

      this%hoh = .false.
      this%ccor_2c = .false.
      this%local_axis = .false.
      this%imported_effective = .true.
      if (meta%bounds(2) > meta%bounds(1)) then
         this%bounds%e_min = meta%bounds(1)
         this%bounds%e_max = meta%bounds(2)
      end if
      this%operator_generation = 1
      fermi = meta%fermi
      if (present(mesh_fermi)) then
         mesh_fermi = fermi
         if (meta%have_mesh_fermi) mesh_fermi = meta%mesh_fermi
      end if
      if (rank == 0) then
         if (meta%had_ccor) call g_logger%info('Effective current includes offsite CCOR through the complete hopping.', &
                                              __FILE__, __LINE__)
         if (meta%had_hoh .and. meta%max_velocity_scale > meta%min_velocity_scale) &
            call g_logger%info('Nonuniform native HOH velocity weights are not represented by effective hopping weights.', &
                               __FILE__, __LINE__)
         call g_logger%info('Imported effective Hamiltonian as an ordinary tight-binding graph.', __FILE__, __LINE__)
      end if
   end subroutine import_effective

   !> Derive the site map from the explicit cluster built by the lattice driver.
   subroutine external_metadata(lat, fermi, meta)
      class(lattice), intent(in) :: lat
      real(rp), intent(in) :: fermi
      type(effective_metadata), intent(out) :: meta
      integer :: i, itype, u, ios, nperiod(3)
      logical :: exists

      if (.not. ieee_is_finite(fermi)) &
         call g_logger%fatal('External TB requires a finite Fermi energy in Ry.', __FILE__, __LINE__)
      meta%fermi = fermi
      meta%sites = lat%kk
      meta%orbitals = norb
      meta%block_size = nb
      meta%nsp = lat%control%nsp
      meta%boundary = lat%pbc .and. [lat%b1,lat%b2,lat%b3]
      meta%cell = lat%a*lat%alat
      nperiod = [lat%n1,lat%n2,lat%n3]
      do i = 1, 3
         meta%cell(:,i) = meta%cell(:,i)*max(1,nperiod(i))
      end do
      allocate(meta%position(3,lat%kk),meta%types(lat%kk),meta%labels(lat%kk))
      do i = 1, lat%kk
         itype = lat%iz(i)
         meta%position(:,i) = lat%cr(:,i)*lat%alat
         meta%types(i) = itype
         meta%labels(i) = lat%symbolic_atoms(itype)%element%symbol
      end do

      inquire(file='rs_tb_orbital_map.dat',exist=exists)
      if (.not. exists) return
      open(newunit=u,file='rs_tb_orbital_map.dat',status='old',action='read',iostat=ios)
      if (ios /= 0) call g_logger%fatal('Cannot open rs_tb_orbital_map.dat.', __FILE__, __LINE__)
      call read_orbital_map(u,meta)
      read(u,*,iostat=ios) i
      if (ios /= iostat_end) call g_logger%fatal('Unexpected extra records in rs_tb_orbital_map.dat.', __FILE__, __LINE__)
      close(u)
   end subroutine external_metadata

   !> One row per external index, giving its one-based RS flat index.
   subroutine read_orbital_map(u, meta)
      integer, intent(in) :: u
      type(effective_metadata), intent(inout) :: meta
      integer :: i, source, target, ios, n
      logical, allocatable :: seen(:)

      if (meta%sites < 1 .or. meta%block_size < 1 .or. allocated(meta%orbital_map)) goto 900
      if (meta%sites > huge(n)/meta%block_size) goto 900
      n = meta%sites*meta%block_size
      allocate(meta%orbital_map(n),seen(n),stat=ios)
      if (ios /= 0) goto 900
      seen = .false.
      do i = 1, n
         read(u,*,iostat=ios) source,target
         if (ios /= 0) goto 900
         if (source /= i .or. target < 1 .or. target > n) goto 900
         if (seen(target)) goto 900
         seen(target) = .true.
         meta%orbital_map(i) = target
      end do
      return
900   call g_logger%fatal('Invalid orbital map. Expected a permutation of all RS flat indices.', __FILE__, __LINE__)
   end subroutine read_orbital_map

   function legacy_filename() result(name)
      character(len=32) :: name
      logical :: exists

      name = 'paoham.dat'
      inquire(file=name,exist=exists)
      if (exists) return
      name = 'rs2paoham.dat'
      inquire(file=name,exist=exists)
      if (.not. exists) name = 'rs_tb_paoham.dat'
   end function legacy_filename

   subroutine write_external_metadata(meta)
      type(effective_metadata), intent(in) :: meta
      integer :: u, ios, i
      logical :: exists
      character(len=8) :: checksum

      inquire(file='rs_tb.meta',exist=exists)
      if (exists) then
         call g_logger%info('External TB used the input lattice and Fermi energy. Existing rs_tb.meta was retained.', &
                            __FILE__, __LINE__)
         return
      end if
      checksum = hopping_checksum(trim(legacy_filename()))
      if (len_trim(checksum) == 0) goto 900
      open(newunit=u,file='rs_tb.meta.tmp',status='replace',action='write',iostat=ios)
      if (ios /= 0) goto 900
      write(u,'(a)',err=900) '# rs-tb-meta-v3', 'representation effective', 'scope cluster', &
         'basis native_spherical', 'energy_unit eV', 'length_unit Angstrom', &
         'energy_reference unscaled_unshifted', 'spin_frame global', 'translation column_relative_to_row', &
         'orbital_order spin_then_l_then_m_decreasing', 'flat_index spin_then_site_then_orbital', 'index_base 1', &
         'fourier_phase positive_k_dot_column_minus_row', 'legacy_fnv1a32 '//checksum
      write(u,'(a,1x,i0)',err=900) 'n_atoms',meta%sites,'norb_scalar',meta%orbitals, &
         'block_size',meta%block_size,'nsp',meta%nsp
      write(u,'(a,3(1x,l1))',err=900) 'periodic',meta%boundary
      write(u,'(a,1x,es25.16)',err=900) 'ry_to_eV',ry2ev,'fermi_Ry',meta%fermi
      write(u,'(a)',err=900) 'lattice_vectors_cart'
      do i = 1, 3
         write(u,'(3(es25.16,1x))',err=900) meta%cell(:,i)
      end do
      write(u,'(a)',err=900) 'positions_cart'
      do i = 1, meta%sites
         write(u,'(i0,1x,3(es25.16,1x),2(i0,1x),a)',err=900) &
            i,meta%position(:,i),i,meta%types(i),trim(meta%labels(i))
      end do
      if (allocated(meta%orbital_map)) then
         write(u,'(a)',err=900) 'orbital_index_map'
         do i = 1, size(meta%orbital_map)
            write(u,'(i0,1x,i0)',err=900) i,meta%orbital_map(i)
         end do
      end if
      write(u,'(a,1x,i0)',err=900) 'legacy_records',meta%legacy_records,'end_records',meta%legacy_records
      close(u,iostat=ios)
      if (ios /= 0) goto 900
      ios = publish_hamiltonian_file('rs_tb.meta.tmp','rs_tb.meta')
      if (ios /= 0) goto 900
      return
900   call g_logger%fatal('Cannot write external TB metadata rs_tb.meta.', __FILE__, __LINE__)
   end subroutine write_external_metadata

   !> @brief Read the file conventions and geometry required by the native reader.
   subroutine read_metadata(meta)
      type(effective_metadata), intent(out) :: meta
      character(len=512) :: line
      character(len=48) :: keyword, value
      integer :: u, ios, i, idx, site_id, nt, p, q
      real(rp) :: velocity_scale, re, im
      logical :: representation, scope, basis, energy, length, reference, frame, have_fermi, convention
      real(rp) :: conversion
      logical :: orbital_order, flat_index, index_base, fourier_phase, have_boundary

      orbital_order = .false.
      flat_index = .false.
      index_base = .false.
      fourier_phase = .false.
      have_boundary = .false.
      representation = .false.
      scope = .false.
      basis = .false.
      energy = .false.
      length = .false.
      reference = .false.
      frame = .false.
      have_fermi = .false.
      convention = .false.
      conversion = 0.0_rp
      open(newunit=u,file='rs_tb.meta',status='old',action='read',iostat=ios)
      if (ios /= 0) goto 900
      read(u,'(a)',iostat=ios) line
      if (ios /= 0) goto 900
      select case (line)
      case ('# rs-tb-meta-v3')
         meta%checked_payload = .true.
      case ('# rs-tb-meta-v2')
         if (rank == 0) call g_logger%warning('Reading v2 Hamiltonian metadata without a payload checksum.', __FILE__, __LINE__)
      case default
         goto 900
      end select
      do
         read(u,'(a)',iostat=ios) line
         if (ios == iostat_end) exit
         if (ios /= 0) goto 900
         read(line,*,iostat=ios) keyword
         if (ios /= 0) goto 900
         select case (keyword)
         case ('representation','scope','basis','energy_unit','length_unit','energy_reference','spin_frame','translation', &
               'orbital_order','flat_index','index_base','fourier_phase','canonical_fnv1a32','legacy_fnv1a32')
            read(line,*,iostat=ios) keyword,value
            if (ios /= 0) goto 900
            select case (keyword)
            case ('representation'); representation = value == 'effective'
            case ('scope')
               scope = value == 'periodic' .or. value == 'cluster'
               meta%periodic = value == 'periodic'
            case ('basis'); basis = value == 'native_spherical'
            case ('energy_unit'); energy = value == 'eV'
            case ('length_unit'); length = value == 'Angstrom'
            case ('energy_reference'); reference = value == 'unscaled_unshifted'
            case ('spin_frame'); frame = value == 'global'
            case ('translation'); convention = value == 'column_relative_to_row'
            ! files written before the rename carry '..._increasing' for the same layout (m descending)
            case ('orbital_order'); orbital_order = value == 'spin_then_l_then_m_decreasing' .or. value == 'spin_then_l_then_m_increasing'
            case ('flat_index'); flat_index = value == 'spin_then_site_then_orbital'
            case ('index_base'); index_base = value == '1'
            case ('fourier_phase'); fourier_phase = value == 'positive_k_dot_column_minus_row'
            case ('canonical_fnv1a32'); meta%canonical_checksum = value
            case ('legacy_fnv1a32'); meta%legacy_checksum = value
            end select
         case ('n_atoms'); read(line,*,iostat=ios) keyword,meta%sites
         case ('norb_scalar'); read(line,*,iostat=ios) keyword,meta%orbitals
         case ('block_size'); read(line,*,iostat=ios) keyword,meta%block_size
         case ('nsp'); read(line,*,iostat=ios) keyword,meta%nsp
         case ('periodic')
            read(line,*,iostat=ios) keyword,meta%boundary
            have_boundary = ios == 0
         case ('ry_to_eV'); read(line,*,iostat=ios) keyword,conversion
         case ('chebyshev_bounds_Ry'); read(line,*,iostat=ios) keyword,meta%bounds
         case ('hoh_original'); read(line,*,iostat=ios) keyword,meta%had_hoh
         case ('ccor_original'); read(line,*,iostat=ios) keyword,meta%had_ccor
         case ('type_velocity_scale')
            read(line,*,iostat=ios) keyword,idx,velocity_scale
            if (ios /= 0) goto 900
            if (.not. ieee_is_finite(velocity_scale)) goto 900
            meta%min_velocity_scale = min(meta%min_velocity_scale,velocity_scale)
            meta%max_velocity_scale = max(meta%max_velocity_scale,velocity_scale)
         case ('fermi_Ry')
            read(line,*,iostat=ios) keyword,meta%fermi
            have_fermi = ios == 0 .and. ieee_is_finite(meta%fermi)
         case ('mesh_fermi_Ry')
            read(line,*,iostat=ios) keyword,meta%mesh_fermi
            if (ios /= 0 .or. .not. ieee_is_finite(meta%mesh_fermi)) goto 900
            meta%have_mesh_fermi = .true.
         case ('lattice_vectors_cart')
            do i = 1, 3
               read(u,*,iostat=ios) meta%cell(:,i)
               if (ios /= 0) goto 900
            end do
         case ('onsite_soc_types')
            read(line,*,iostat=ios) keyword,nt
            if (ios /= 0 .or. nt <= 0 .or. allocated(meta%onsite_soc)) goto 900
            allocate(meta%onsite_soc(nb,nb,nt), source=czero, stat=ios)
            if (ios /= 0) goto 900
         case ('onsite_soc_Ry')
            read(line,*,iostat=ios) keyword,idx,p,q,re,im
            if (ios /= 0 .or. .not. allocated(meta%onsite_soc)) goto 900
            if (idx < 1 .or. idx > size(meta%onsite_soc,3) .or. min(p,q) < 1 .or. max(p,q) > nb) goto 900
            if (.not. ieee_is_finite(re) .or. .not. ieee_is_finite(im)) goto 900
            if (meta%onsite_soc(p,q,idx) /= czero) goto 900
            meta%onsite_soc(p,q,idx) = cmplx(re,im,rp)
         case ('positions_cart')
            if (meta%sites <= 0 .or. allocated(meta%position)) goto 900
            allocate(meta%position(3,meta%sites),meta%types(meta%sites),meta%labels(meta%sites),stat=ios)
            if (ios /= 0) goto 900
            do i = 1, meta%sites
               read(u,*,iostat=ios) idx,meta%position(:,i),site_id,meta%types(i),meta%labels(i)
               if (ios /= 0 .or. idx /= i .or. site_id <= 0) goto 900
               if (.not. meta%periodic .and. site_id /= i) goto 900
            end do
         case ('orbital_index_map'); call read_orbital_map(u,meta)
         case ('legacy_records'); read(line,*,iostat=ios) keyword,meta%legacy_records
         case ('end_records')
            read(line,*,iostat=ios) keyword,meta%records
            if (ios /= 0) goto 900
            exit
         end select
         if (ios /= 0) goto 900
      end do
      close(u,iostat=ios)
      if (ios /= 0 .or. meta%records < 1 .or. .not. allocated(meta%position)) goto 900
      if (.not. (representation .and. scope .and. basis .and. energy .and. length .and. reference .and. frame .and. &
          have_fermi .and. convention .and. orbital_order .and. flat_index .and. index_base .and. fourier_phase .and. &
          have_boundary)) goto 900
      if (meta%checked_payload .and. &
          len_trim(meta%canonical_checksum) == 0 .and. len_trim(meta%legacy_checksum) == 0) goto 900
      if (.not. ieee_is_finite(conversion) .or. conversion /= ry2ev) goto 900
      if (.not. all(ieee_is_finite(meta%bounds))) goto 900
      if (.not. all(ieee_is_finite(meta%cell)) .or. .not. all(ieee_is_finite(meta%position))) goto 900
      return
900   call g_logger%fatal('Invalid or incomplete rs_tb.meta for effective real-space import.', __FILE__, __LINE__)
   end subroutine read_metadata

   !> @brief Load scalar entries into sparse spherical-basis hopping rows.
   subroutine read_records(meta, legacy, rows)
      type(effective_metadata), intent(inout) :: meta
      logical, intent(in) :: legacy
      type(effective_row), allocatable, intent(out) :: rows(:)
      character(len=512) :: line
      character(len=32) :: comment, label, unit_name
      character(len=8) :: checksum, expected_checksum
      integer :: u, ios, i, j, p, q, k, rec, ip, jp, r(3), a, b
      integer(int64) :: read_count, end_count
      real(rp) :: re, im
      complex(rp) :: empty(nb,nb), value, block(nb,nb)

      allocate(rows(meta%sites))
      empty = czero
      do i = 1, meta%sites
         call rows(i)%reset(nb)
      end do
      end_count = -1
      read_count = 0
      unit_name = 'rs_tb.tb'
      if (legacy) unit_name = legacy_filename()
      if (meta%checked_payload) then
         expected_checksum = meta%canonical_checksum
         if (legacy) expected_checksum = meta%legacy_checksum
         checksum = hopping_checksum(trim(unit_name))
         if (len_trim(expected_checksum) == 0 .or. checksum /= expected_checksum) &
            call g_logger%fatal('Effective Hamiltonian payload does not match rs_tb.meta: '//trim(unit_name), __FILE__, __LINE__)
      end if
      open(newunit=u,file=unit_name,status='old',action='read',iostat=ios)
      if (ios /= 0) goto 900
      if (.not. legacy) then
         read(u,'(a)',iostat=ios) line
         if (ios /= 0 .or. line /= '# rs-tb-hoppings-v2') goto 900
      end if
      do
         read(u,'(a)',iostat=ios) line
         if (ios == iostat_end) exit
         if (ios /= 0) goto 900
         if (len_trim(line) == 0) cycle
         if (index(line,'# end records ') == 1) then
            read(line(15:),*,iostat=ios) end_count
            if (ios /= 0) goto 900
            exit
         end if
         if (line(1:1) == '#') cycle
         if (legacy) then
            read(line,*,iostat=ios) r,ip,jp,re,im
            if (ios /= 0) goto 900
            if (min(ip,jp) < 1 .or. max(ip,jp) > int(nb,int64)*meta%sites) goto 900
            if (allocated(meta%orbital_map)) then
               ip = meta%orbital_map(ip)
               jp = meta%orbital_map(jp)
            end if
            call orb2site(ip,p,i,meta%sites,norb)
            call orb2site(jp,q,j,meta%sites,norb)
         else
            read(line,*,iostat=ios) r,i,j,rec,p,q,ip,jp,re,im
            if (ios /= 0 .or. rec <= 0) goto 900
         end if
         if (min(i,j) < 1 .or. max(i,j) > meta%sites .or. min(p,q) < 1 .or. max(p,q) > nb) goto 900
         if (any(abs(int(r,int64)) > huge(0)/2) .or. any(r /= 0 .and. .not. meta%boundary)) goto 900
         if (.not. ieee_is_finite(re) .or. .not. ieee_is_finite(im)) goto 900
         if (.not. legacy) then
            if (ip /= int(((p-1)/norb),int64)*meta%sites*norb+int(i-1,int64)*norb+modulo(p-1,norb)+1 .or. &
                jp /= int(((q-1)/norb),int64)*meta%sites*norb+int(j-1,int64)*norb+modulo(q-1,norb)+1) goto 900
         end if
         read_count = read_count+1
         value = cmplx(re,im,rp)/ry2ev
         if (value == czero) then
            if (.not. legacy) goto 900
            cycle
         end if
         k = rows(i)%lookup([j,r])
         if (k == 0) then
            call rows(i)%add([j,r],empty)
            k = rows(i)%count
         end if
         if (rows(i)%block(p,q,k) /= czero) goto 900
         rows(i)%block(p,q,k) = value
      end do
      close(u,iostat=ios)
      if (ios /= 0 .or. read_count == 0) goto 900
      if (.not. legacy) then
         if (read_count /= end_count .or. read_count /= meta%records) goto 900
      else
         if (meta%legacy_records >= 0 .and. read_count /= meta%legacy_records) goto 900
         meta%legacy_records = read_count
         do i = 1, meta%sites
            do k = 1, rows(i)%count
               block = rows(i)%block(:,:,k)
               do b = 1, nb, norb
                  do a = 1, nb, norb
                     call hcpx(block(a:a+norb-1,b:b+norb-1),'cart2sph')
                  end do
               end do
               rows(i)%block(:,:,k) = block
            end do
         end do
      end if
      return
900   call g_logger%fatal('Invalid, duplicate or truncated effective Hamiltonian records in '//trim(unit_name), __FILE__, __LINE__)
   end subroutine read_records
end submodule hamiltonian_effective_io
