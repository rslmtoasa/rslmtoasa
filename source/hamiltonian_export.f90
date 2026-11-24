submodule(hamiltonian_mod) hamiltonian_export
   use iso_fortran_env, only: int64
   use ieee_arithmetic, only: ieee_is_finite
   use hamiltonian_effective_mod, only: resolve_translation, hopping_checksum
   use magnetic_representation_mod, only: gbt_single_q, explicit_texture
   implicit none
   real(rp), parameter :: geometry_tol = 1.e-6_rp

contains

   !> @brief Map each directed native neighbor slot to an effective column key.
   !> @details Keys are [column, image(1:3)]. The image counts primitive cells
   !>          for compact output and whole replicas for finite-cluster output.
   module subroutine effective_graph(this, periodic, sites, keys)
      class(hamiltonian), intent(in) :: this
      logical, intent(in) :: periodic
      integer, allocatable, intent(out) :: sites(:), keys(:,:,:)
      class(lattice), pointer :: lat
      real(rp), allocatable :: position(:,:)
      real(rp) :: cell(3,3), inverse(3,3), delta(3), displacement(3)
      integer :: i, j, k, a, b, r(3), n, slots
      logical :: found, ambiguous, wrapped

      lat => this%lattice
      wrapped = lat%pbc .and. .not. periodic
      if (periodic) then
         if (lat%control%calctype /= 'B' .or. lat%nmax /= 0 .or. &
             this%magnetic_representation == explicit_texture) &
            call g_logger%fatal('Compact periodic export requires a repeating bulk basis; use cluster scope.', __FILE__, __LINE__)
         if (lat%nbulk /= 0 .or. lat%nrec /= lat%ntype .or. lat%nrec /= lat%ntot) &
            call g_logger%fatal('Compact periodic export requires one recursion class per basis site.', __FILE__, __LINE__)
         sites = lat%atlist
      else
         sites = [(i, i=1,lat%kk)]
      end if
      cell = lat%a*lat%alat
      if (wrapped) then
         cell(:,1) = cell(:,1)*max(1,lat%n1)
         cell(:,2) = cell(:,2)*max(1,lat%n2)
         cell(:,3) = cell(:,3)*max(1,lat%n3)
         position = lat%cr(:,:lat%kk)*lat%alat
      end if
      if (periodic .or. wrapped) then
         if (.not. all(ieee_is_finite(cell)) .or. abs(determinant(cell)) < spacing(1.0_rp)) &
            call g_logger%fatal('Invalid cell for Hamiltonian export.', __FILE__, __LINE__)
         inverse = inverse_3x3(cell)
      end if

      n = size(sites)
      slots = maxval(lat%nn(sites,1))
      allocate(keys(4,slots,n), source=0)
      do i = 1, n
         a = sites(i)
         keys(:,1,i) = [i,0,0,0]
         do k = 2, lat%nn(a,1)
            b = lat%nn(a,k)
            if (b == 0) cycle
            if (b < 1 .or. b > lat%kk) &
               call g_logger%fatal('Invalid native neighbor in Hamiltonian export.', __FILE__, __LINE__)
            j = b
            r = 0
            if (periodic) then
               j = lat%iz(b)
               if (j < 1 .or. j > n) &
                  call g_logger%fatal('Invalid periodic basis index in Hamiltonian export.', __FILE__, __LINE__)
               delta = (lat%cr(:,b) - lat%cr(:,sites(j)))*lat%alat
               ! A closed representative row can itself cross a supercell face.
               if (lat%pbc) then
                  if (.not. allocated(position)) position = lat%cr(:,:lat%kk)*lat%alat
                  call lat%f_wrap_coord_diff(lat%kk, position, a, b, displacement, ambiguous)
                  if (ambiguous) call g_logger%fatal('Ambiguous periodic Hamiltonian bond: enlarge the supercell.', __FILE__, __LINE__)
                  delta = displacement + (lat%cr(:,a) - lat%cr(:,sites(j)))*lat%alat
               end if
               call resolve_translation(cell, inverse, delta, geometry_tol, r, found)
               if (.not. found) call g_logger%fatal('Native bond does not resolve to a periodic lattice index.', __FILE__, __LINE__)
            else if (wrapped) then
               call lat%f_wrap_coord_diff(lat%kk, position, a, b, displacement, ambiguous)
               if (ambiguous) call g_logger%fatal('Ambiguous cluster Hamiltonian bond: enlarge the supercell.', __FILE__, __LINE__)
               delta = displacement - (position(:,b) - position(:,a))
               call resolve_translation(cell, inverse, delta, geometry_tol, r, found)
               if (.not. found) call g_logger%fatal('Native bond does not resolve to a cluster image.', __FILE__, __LINE__)
            end if
            keys(:,k,i) = [j,r]
         end do
      end do
   end subroutine effective_graph

   !> @brief Read the native bare block in the global spin frame.
   function bare_block(this, site, slot) result(block)
      class(hamiltonian), intent(in) :: this
      integer, intent(in) :: site, slot
      complex(rp) :: block(nb,nb)
      integer :: t
      t = this%lattice%iz(site)
      if (site <= this%lattice%nmax) then
         if (this%local_axis) then
            block = this%hall_glob(:,:,slot,site)
         else
            block = this%hall(:,:,slot,site)
         end if
      else
         if (this%local_axis) then
            block = this%ee_glob(:,:,slot,t)
         else
            block = this%ee(:,:,slot,t)
         end if
      end if
   end function bare_block

   !> @brief Assemble one row of the unscaled effective Hamiltonian.
   !> @details Add native CCOR and onsite terms once. Only bare hoppings enter
   !>          the HOH product. The two directed image indices add along a path.
   module subroutine assemble_effective_row(this, sites, keys, i, row)
      class(hamiltonian), intent(in) :: this
      integer, intent(in) :: sites(:), keys(:,:,:), i
      type(effective_row), intent(inout) :: row
      integer :: a, t, k, l, j, jtype, key(4)
      complex(rp) :: hop(nb,nb), right(nb,nb), weighted(nb,nb), term(nb,nb)
      class(lattice), pointer :: lat

      lat => this%lattice
      a = sites(i)
      t = lat%iz(a)
      call row%reset(nb)
      do k = 1, lat%nn(a,1)
         j = keys(1,k,i)
         if (j == 0) cycle
         hop = bare_block(this, a, k)
         term = hop
         if (this%ccor_2c) then
            if (a <= lat%nmax) then
               if (this%local_axis) then
                  term = term + this%hallcc_glob(:,:,k,a)
               else
                  term = term + this%hallcc(:,:,k,a)
               end if
            else
               if (this%local_axis) then
                  term = term + this%eecc_glob(:,:,k,t)
               else
                  term = term + this%eecc(:,:,k,t)
               end if
            end if
         end if
         call row%add(keys(:,k,i), term)
         if (.not. this%hoh) cycle
         jtype = lat%iz(sites(j))
         call zgemm('N','N',nb,nb,nb,cone,hop,nb,this%obarm(:,:,jtype),nb,czero,weighted,nb)
         do l = 1, lat%nn(sites(j),1)
            key = keys(:,l,j)
            if (key(1) == 0) cycle
            right = bare_block(this, sites(j), l)
            call zgemm('N','N',nb,nb,nb,-cone,weighted,nb,right,nb,czero,term,nb)
            key(2:4) = key(2:4) + keys(2:4,k,i)
            call row%add(key, term)
         end do
      end do
      term = this%lsham(:,:,t)
      if (this%hoh) then
         if (this%local_axis) then
            term = term + this%enim_glob(:,:,t)
         else
            term = term + this%enim(:,:,t)
         end if
      end if
      call row%add([i,0,0,0], term)
   end subroutine assemble_effective_row

   !> @brief Write both real-space views of the assembled effective matrix.
   !> @details Native HOH, SOC and Hubbard operands are read-only. The metadata
   !>          file is published last and fingerprints each completed payload.
   !>          In the seven-column cluster view, indices refer to actual sites
   !>          and translations count whole replica cells, as in the metadata.
   module subroutine export_rs_tb_all(this, fermi, fermi_provenance, mesh_fermi)
      class(hamiltonian), intent(in) :: this
      real(rp), optional, intent(in) :: fermi
      character(len=*), optional, intent(in) :: fermi_provenance
      real(rp), optional, intent(in) :: mesh_fermi
      type(effective_row) :: row
      integer, allocatable :: sites(:), keys(:,:,:), order(:)
      integer :: i, j, l, m, rec, k, r(3), ipao, jpao, u, meta, legacy, ios, a, b, onsite
      integer(int64) :: count, paths, legacy_count, unique_blocks, original_bonds
      logical :: periodic, flags(3), mixed_boundaries
      character(len=8) :: canonical_checksum, legacy_checksum
      complex(rp) :: block(nb,nb)
      real(rp) :: cell(3,3)
      character(len=32) :: scope, legacy_name
      class(lattice), pointer :: lat

      if (rank /= 0 .or. this%export == 'none') return
      lat => this%lattice
      if (this%export_scope == 'kspace' .or. this%magnetic_representation == gbt_single_q) &
         call g_logger%fatal('This operator requires reciprocal export.', __FILE__, __LINE__)
      if (.not. allocated(this%ee) .or. .not. allocated(this%lsham) .or. this%operator_generation == 0) &
         call g_logger%fatal('Hamiltonian export requires a built native operator.', __FILE__, __LINE__)
      if (this%hoh .and. (.not. allocated(this%obarm) .or. .not. allocated(this%enim))) &
         call g_logger%fatal('HOH export requires the built obarm and enim.', __FILE__, __LINE__)
      if (lat%nmax > 0 .and. .not. allocated(this%hall)) &
         call g_logger%fatal('Cluster export requires the built local Hamiltonian.', __FILE__, __LINE__)
      if (this%ccor_2c) then
         if (.not. allocated(this%eecc) .or. (lat%nmax > 0 .and. .not. allocated(this%hallcc))) &
            call g_logger%fatal('CCOR export requires the built CCOR blocks.', __FILE__, __LINE__)
      end if
      flags = lat%pbc .and. [lat%b1,lat%b2,lat%b3]
      mixed_boundaries = any(flags) .and. .not. all(flags)
      periodic = this%export_scope == 'periodic'
      if (this%export_scope == 'auto') periodic = lat%control%calctype == 'B' .and. lat%nmax == 0 .and. &
                                                 this%magnetic_representation /= explicit_texture .and. .not. mixed_boundaries
      scope = 'cluster'
      if (periodic) scope = 'periodic'
      if (nb /= 2*norb) &
         call g_logger%fatal('Seven-column export requires a two-spin-block basis.', __FILE__, __LINE__)
      call this%effective_graph(periodic, sites, keys)
      if (nb > huge(0)/size(sites)) call g_logger%fatal('Hamiltonian export exceeds orbital index capacity.', __FILE__, __LINE__)
      count = 0_int64
      paths = 0_int64
      legacy_count = 0_int64
      unique_blocks = 0_int64
      original_bonds = count_nonzero_graph(keys)
      open(newunit=u,file='rs_tb.tb.tmp',status='replace',action='write',iostat=ios)
      if (ios /= 0) goto 900
      open(newunit=meta,file='rs_tb.meta.tmp',status='replace',action='write',iostat=ios)
      if (ios /= 0) goto 900
      legacy_name = 'rs_tb_paoham.dat'
      if (this%export == 'rs2pao') legacy_name = 'rs2paoham.dat'
      open(newunit=legacy,file=trim(legacy_name)//'.tmp',status='replace',action='write',iostat=ios)
      if (ios /= 0) goto 900
      write(u,'(a)',err=900) '# rs-tb-hoppings-v2', &
         '# R1 R2 R3 ia jbasis effective_record i_local j_local ipao jpao real_eV imag_eV'
      do i = 1, size(sites)
         call this%assemble_effective_row(sites, keys, i, row)
         order = row%order()
         ! Preserve the legacy first record: the zero-image onsite block.
         onsite = findloc(order,row%lookup([i,0,0,0]),dim=1)
         rec = order(onsite)
         order(2:onsite) = order(1:onsite-1)
         order(1) = rec
         if (this%hoh) then
            do k = 1, lat%nn(sites(i),1)
               j = keys(1,k,i)
               if (j /= 0) paths = paths + count_nonzero(keys(1,:,j))
            end do
         end if
         do k = 1, row%count
            rec = order(k)
            j = row%key(1,rec)
            r = row%key(2:4,rec)
            block = row%block(:,:,rec)
            if (all(block == czero)) cycle
            unique_blocks = unique_blocks + 1
            if (.not. all(ieee_is_finite(real(block,rp))) .or. .not. all(ieee_is_finite(aimag(block)))) &
               call g_logger%fatal('Nonfinite effective Hamiltonian block.', __FILE__, __LINE__)
            do m = 1, nb
               do l = 1, nb
                  if (block(l,m) == czero) cycle
                  ipao = flat_orbital(l,i,size(sites))
                  jpao = flat_orbital(m,j,size(sites))
                  write(u,'(10(i0,1x),2(es26.17e3,1x))',err=900) r,i,j,k,l,m,ipao,jpao, &
                     real(block(l,m),rp)*ry2ev, aimag(block(l,m))*ry2ev
                  count = count + 1
               end do
            end do
            do b = 1, nb, norb
               do a = 1, nb, norb
                  call hcpx(block(a:a+norb-1,b:b+norb-1),'sph2cart')
               end do
            end do
            do m = 1, nb
               do l = 1, nb
                  if (block(l,m) == czero) cycle
                  call site2orb(l,i,ipao,size(sites),norb)
                  call site2orb(m,j,jpao,size(sites),norb)
                  write(legacy,'(5(i0,1x),2(es26.17e3,1x))',err=900) r,ipao,jpao, &
                     real(block(l,m),rp)*ry2ev, aimag(block(l,m))*ry2ev
                  legacy_count = legacy_count + 1
               end do
            end do
         end do
      end do
      write(u,'(a,i0)',err=900) '# end records ', count
      close(u,iostat=ios)
      if (ios /= 0) goto 900
      close(legacy,iostat=ios)
      if (ios /= 0) goto 900
      canonical_checksum = hopping_checksum('rs_tb.tb.tmp')
      if (len_trim(canonical_checksum) == 0) goto 900
      legacy_checksum = hopping_checksum(trim(legacy_name)//'.tmp')
      if (len_trim(legacy_checksum) == 0) goto 900
      flags = [.true.,.true.,.true.]
      cell = lat%a*lat%alat
      if (.not. periodic) then
         flags = lat%pbc .and. [lat%b1,lat%b2,lat%b3]
         cell(:,1) = cell(:,1)*max(1,lat%n1)
         cell(:,2) = cell(:,2)*max(1,lat%n2)
         cell(:,3) = cell(:,3)*max(1,lat%n3)
      end if
      write(meta,'(a)',err=900) '# rs-tb-meta-v3', 'representation effective', 'scope '//trim(scope), &
         'energy_unit eV', 'length_unit Angstrom', 'basis native_spherical', 'spin_frame global', &
         'orbital_order spin_then_l_then_m_increasing', 'flat_index spin_then_site_then_orbital', &
         'energy_reference unscaled_unshifted', 'translation column_relative_to_row', &
         'fourier_phase positive_k_dot_column_minus_row', 'index_base 1', 'pruning_threshold_eV 0', &
         'source_commit '//RSLMTO_SOURCE_COMMIT, 'magnetic_representation '//trim(this%magnetic_representation), &
         'hubbard_radial_moments '//trim(this%hubbard_radial_moments)
      write(meta,'(a,i0)',err=900) 'n_atoms ',size(sites), 'norb_scalar ',norb, 'block_size ',nb, 'nsp ',lat%control%nsp
      write(meta,'(a,l1)',err=900) 'hoh_original ',this%hoh, 'soc_original ',lat%control%has_soc(), &
         'ccor_original ',this%ccor_2c, 'hubbard_u_original ',this%hubbard_u_general_check .or. &
         this%hubbard_u_impurity_check .or. this%hubbard_u_sc_check, 'hubbard_v_original ',this%hubbard_v_check
      write(meta,'(a,3(1x,l1))',err=900) 'periodic',flags
      write(meta,'(a,3(1x,l1))',err=900) 'source_boundary',lat%pbc .and. [lat%b1,lat%b2,lat%b3]
      write(meta,'(a,es26.17e3)',err=900) 'alat_Angstrom ',lat%alat, 'ry_to_eV ',ry2ev
      write(meta,'(a,2(1x,es26.17e3))',err=900) 'chebyshev_bounds_Ry', this%bounds%e_min, this%bounds%e_max
      if (present(fermi)) then
         write(meta,'(a,es26.17e3)',err=900) 'fermi_Ry ',fermi
         if (present(fermi_provenance)) then
            write(meta,'(a)',err=900) 'fermi_provenance '//fermi_provenance
         else
            write(meta,'(a)',err=900) 'fermi_provenance native_energy_object'
         end if
      else
         write(meta,'(a)',err=900) 'fermi_provenance unavailable'
      end if
      if (present(mesh_fermi)) write(meta,'(a,es26.17e3)',err=900) 'mesh_fermi_Ry ',mesh_fermi
      write(meta,'(a)',err=900) 'lattice_vectors_cart'
      do i = 1, 3
         write(meta,'(3(es26.17e3,1x))',err=900) cell(:,i)
      end do
      do i = 1, lat%ntype
         write(meta,'(a,i0,1x,*(es26.17e3,1x))',err=900) 'type_hubbard_U_Ry ',i,lat%symbolic_atoms(i)%potential%hubbard_u
         write(meta,'(a,i0,1x,*(es26.17e3,1x))',err=900) 'type_hubbard_J_Ry ',i,lat%symbolic_atoms(i)%potential%hubbard_j
         write(meta,'(a,i0,1x,es26.17e3)',err=900) 'type_velocity_scale ',i,this%velocity_scale(i)
      end do
      write(meta,'(a)',err=900) 'positions_cart'
      do i = 1, size(sites)
         j = sites(i)
         write(meta,'(i0,1x,3(es26.17e3,1x),2(i0,1x),a)',err=900) &
            i,lat%cr(:,j)*lat%alat,j,lat%iz(j),trim(lat%symbolic_atoms(lat%iz(j))%element%symbol)
      end do
      write(meta,'(a)',err=900) 'canonical_fnv1a32 '//canonical_checksum
      write(meta,'(a)',err=900) 'legacy_fnv1a32 '//legacy_checksum
      write(meta,'(a,i0)',err=900) 'original_bond_records ',original_bonds,'effective_blocks ',unique_blocks, &
         'two_hop_paths ',paths,'legacy_records ',legacy_count,'end_records ',count
      close(meta,iostat=ios)
      if (ios /= 0) goto 900
      call rename('rs_tb.tb.tmp','rs_tb.tb',ios)
      if (ios /= 0) goto 900
      call rename(trim(legacy_name)//'.tmp',trim(legacy_name),ios)
      if (ios /= 0) goto 900
      call rename('rs_tb.meta.tmp','rs_tb.meta',ios)
      if (ios == 0) return
900   call g_logger%fatal('I/O error writing effective Hamiltonian export.', __FILE__, __LINE__)
   end subroutine export_rs_tb_all

   integer function flat_orbital(local, site, nsites) result(flat)
      integer, intent(in) :: local, site, nsites
      flat = ((local-1)/norb*nsites + site-1)*norb + modulo(local-1,norb)+1
   end function flat_orbital

   integer(int64) function count_nonzero_graph(keys) result(n)
      integer, intent(in) :: keys(:,:,:)
      n = count(keys(1,:,:) /= 0, kind=int64)
   end function count_nonzero_graph

   integer function count_nonzero(columns)
      integer, intent(in) :: columns(:)
      count_nonzero = count(columns /= 0)
   end function count_nonzero
end submodule hamiltonian_export
