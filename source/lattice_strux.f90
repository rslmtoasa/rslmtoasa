submodule (lattice_mod) lattice_strux
   implicit none

contains

   module subroutine structb(this, do_str)
      class(lattice), intent(inout) :: this
      logical, intent(in) :: do_str
      ! Local variables
      integer :: i, ia, nr, ii, j, nm, np, nlim, nomx, ncut, kk, nnmx
      integer :: sbar_dim, nm_store, nt_tmp
      integer, dimension(:, :), allocatable :: nn
      integer, dimension(:), allocatable :: idnn
      logical :: do_str_
      real(rp), dimension(:, :, :), allocatable :: set
      real(rp), dimension(3) :: ret
      real(rp) :: t_structb_start, t_cluster_ready, t_nnmap_start, t_nnmap_end
      real(rp) :: t_remd_start, t_remd_end, t_nm_store_start, t_nm_store_end
      real(rp) :: t_outmap_start, t_outmap_end, t_str_stage_end

      call g_timer%start('structb')

      ! Open files
      open (12, file='map', form='unformatted')
      open (13, file="sbar", FORM="unformatted")
      open (16, file="view.sbar")
      if (this%strux_want_sdot) then
         open (14, file="sdot", form="unformatted")
         open (15, file="view.sdot")
      end if
      open (17, file='str.out')

      ! Clust parameters
      ncut = 9
      nnmx = 100 ! 5250
      nomx = this%ntot
      kk = this%kk
      allocate (set(3, nomx, nnmx)); set = 0.0d0
      allocate (nn(kk, nnmx))
      allocate (idnn(nnmx))
      nm = nnmx
      write (17, *) 'irec', this%nrec, this%irec
      write (17, *) 'irec type', this%iz(this%irec(:))
      write (17, *) 'ndi=', kk
      write (17, 10000) kk
      write (17, 10001)
      write (17, 10002) (i, (this%cr(j, i)*this%alat, j=1, 3), i=1, max(this%nmax, this%ntype))
      ! write (*, '(a,10f10.4)') 'NNCAL ct=', this%ct(1:min(this%ntype, size(this%ct)))
      call this%nncal(this%ct, this%cr*this%alat, 3, kk, this%iz, nn, kk, nm, mapa, this%ntype)

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('lattice.nn', this%nn, (/this%kk, nm + 1/))
#else
      allocate (this%nn(this%kk, nm + 1))
#endif
      do ii = 1, nm + 1
         this%nn(:, ii) = nn(:, ii)
      end do
      write (17, *) 'ndi=', kk
      write (17, *) 'remd'
      call this%remd(this%cr*this%alat, this%num, this%iu, this%nn, kk, this%ntot, nomx, kk, nnmx, set, idnn, ret)
      nm_store = max(1, maxval(this%nn(:, 1)))
      do ii = 1, this%ntot
         ia = this%iu(ii)
         if (ia <= 0) cycle
         nt_tmp = this%kk
         call this%clusba(this%r2, this%cr*this%alat, ia, kk, kk, nt_tmp)
         nm_store = max(nm_store, nt_tmp)
      end do
      this%nn_max = nm_store
      sbar_dim = max(norb, this%control%npold)
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('lattice.sbar', this%sbar, (/sbar_dim, sbar_dim, nm_store, this%ntot/))
#else
      allocate (this%sbar(sbar_dim, sbar_dim, nm_store, this%ntot))
#endif
      call this%init_strux_storage(sbar_dim, nm_store)
      write (17, *) 'outmap', this%nmax, maxval(this%irec)
      call outmap(17, this%iz, this%nn, this%num, kk, nnmx, max(this%nmax, maxval(this%irec)))
      write (17, 10003) kk, nm
      flush(17)
      if (do_str) then
         if (rank==0) call g_logger%info('STRUCTB structure-constants, backend='//trim(this%strux_backend), __FILE__, __LINE__)
         if (this%control%lmax > 2 .and. trim(lower(this%strux_backend)) /= 'strux_lib') then
            call g_logger%fatal('lmax='//int2str(this%control%lmax)//' requires strux_lib backend for structure constants; legacy backend is only valid up to spd (lmax=2).', __FILE__, __LINE__)
         end if
         if (trim(lower(this%strux_backend)) == 'strux_lib') then
            call this%structb_strux()
         else
            do ii = 1, this%ntot
               ia = this%iu(ii)
               nr = this%nn(ia, 1)
               write (17, '(1x, a, i5, a, i5)') 'Sbar atom no:', ii, ' Ntot:', this%ntot
               call this%dbar1(ia, ncut*this%r2, this%wav, this%cr*this%alat, kk, kk, norb, nr, ii)
            end do
         end if
      end if
      call g_timer%stop('structb')
10000 format(i5)
10001 format(" LATTICE COORDINATES")
10002 format(2(i5, 3f8.4))
10003 format(3i5)
10004 format(3i5)
10005 format(7x, i7)
   end subroutine structb

   module subroutine init_strux_storage(this, sbar_dim, nm)
      class(lattice), intent(inout) :: this
      integer, intent(in) :: sbar_dim, nm

#ifdef USE_SAFE_ALLOC
      if (allocated(this%sdot)) call g_safe_alloc%deallocate('lattice.sdot', this%sdot)
      call g_safe_alloc%allocate('lattice.sdot', this%sdot, (/sbar_dim, sbar_dim, nm, this%ntot/))
      if (allocated(this%alpha)) call g_safe_alloc%deallocate('lattice.alpha', this%alpha)
      call g_safe_alloc%allocate('lattice.alpha', this%alpha, (/sbar_dim, this%kk, this%ntot/))
      if (allocated(this%alpha_dot)) call g_safe_alloc%deallocate('lattice.alpha_dot', this%alpha_dot)
      call g_safe_alloc%allocate('lattice.alpha_dot', this%alpha_dot, (/sbar_dim, this%kk, this%ntot/))
#else
      if (allocated(this%sdot)) deallocate (this%sdot)
      allocate (this%sdot(sbar_dim, sbar_dim, nm, this%ntot))
      if (allocated(this%alpha)) deallocate (this%alpha)
      allocate (this%alpha(sbar_dim, this%kk, this%ntot))
      if (allocated(this%alpha_dot)) deallocate (this%alpha_dot)
      allocate (this%alpha_dot(sbar_dim, this%kk, this%ntot))
#endif

      this%sbar(:, :, :, :) = czero
      this%sdot(:, :, :, :) = czero
      this%alpha(:, :, :) = 0.0_rp
      this%alpha_dot(:, :, :) = 0.0_rp
   end subroutine init_strux_storage

   module pure function default_screening_alpha(this, nl) result(alpha_default)
      class(lattice), intent(in) :: this
      integer, intent(in) :: nl
      real(rp) :: alpha_default(nl)
      real(rp), parameter :: default_values(4) = [0.3485_rp, 0.0530_rp, 0.0107_rp, 0.00674_rp]
      integer :: i

      alpha_default = 0.0_rp
      do i = 1, nl
         alpha_default(i) = default_values(min(i, size(default_values)))
      end do
   end function default_screening_alpha

   module integer function strux_mode(this)
      class(lattice), intent(in) :: this
      character(len=:), allocatable :: screening_mode

      screening_mode = trim(lower(this%screening))
      select case (screening_mode)
      case ('manual', 'default')
         strux_mode = STRUX_LMTO47_IALPHA_MANUAL
      case ('sigma')
         strux_mode = STRUX_LMTO47_IALPHA_SIGMA
      case ('fitted')
         strux_mode = STRUX_LMTO47_IALPHA_FITD
      case default
         call g_logger%fatal('Unsupported lattice.screening='//trim(this%screening), __FILE__, __LINE__)
      end select
   end function strux_mode

   module subroutine load_symbolic_atoms_if_needed(this)
      class(lattice), intent(inout) :: this
      
      ! 1. Declare the temporary array. 
      type(symbolic_atom), allocatable :: temp_atoms(:)

      if (allocated(this%symbolic_atoms)) then
         if (size(this%symbolic_atoms) == this%ntype) return
         deallocate(this%symbolic_atoms)
      end if
      
      if (.not. allocated(this%symbolic_atoms)) then
         ! 2. Assign the function result to the standard local variable
         temp_atoms = array_of_symbolic_atoms(this%control%fname, this%ntype)
         
         ! 3. Safely move the memory allocation into the derived type component
         call move_alloc(from=temp_atoms, to=this%symbolic_atoms)
      end if
   end subroutine load_symbolic_atoms_if_needed

   module subroutine build_rmt(this, nspec, species_labels, rmt)
      class(lattice), intent(inout) :: this
      integer, intent(in) :: nspec
      integer, intent(in) :: species_labels(nspec)
      real(rp), intent(out) :: rmt(nspec)
      integer :: is, label

      call this%load_symbolic_atoms_if_needed()

      do is = 1, nspec
         label = species_labels(is)
         if (label < 1 .or. label > size(this%symbolic_atoms)) then
            call g_logger%fatal('strux backend found invalid species label in lattice%no for WS radius lookup', __FILE__, __LINE__)
         end if
         if (this%symbolic_atoms(label)%potential%ws_r <= tiny(1.0_rp)) then
            call g_logger%fatal('strux backend requires positive potential%ws_r for sigma/fitted screening', __FILE__, __LINE__)
         end if
         ! NOTE: potential%ws_r is already stored in Bohr in this code path.
         ! Keep rmt in the same unit to remain consistent with avw_bohr.
         rmt(is) = this%symbolic_atoms(label)%potential%ws_r
      end do
   end subroutine build_rmt

   module subroutine build_hcr(this, nl, nspec, rmt, hcr)
      class(lattice), intent(in) :: this
      integer, intent(in) :: nl, nspec
      real(rp), intent(in) :: rmt(nspec)
      real(rp), intent(out) :: hcr(nl, nspec)
      integer :: is, l

      do is = 1, nspec
         do l = 1, nl
            hcr(l, is) = this%screening_sigma(l)*rmt(is)
         end do
      end do
   end subroutine build_hcr

   module subroutine build_strux_inputs(this, nspec, nl, species_labels, alpha_in, hcr, rmt)
      class(lattice), intent(inout) :: this
      integer, intent(in) :: nspec, nl
      integer, intent(in) :: species_labels(nspec)
      real(rp), intent(out) :: alpha_in(0:nl - 1, nspec)
      real(rp), intent(out) :: hcr(nl, nspec)
      real(rp), intent(out) :: rmt(nspec)
      integer :: is
      real(rp) :: alpha_global(nl)

      call this%build_rmt(nspec, species_labels, rmt)
      call this%build_hcr(nl, nspec, rmt, hcr)

      select case (trim(lower(this%screening)))
      case ('manual')
         alpha_global = this%screening_alpha(1:nl)
      case ('default')
         alpha_global = this%default_screening_alpha(nl)
      case default
         alpha_global = 0.0_rp
      end select

      alpha_in(:, :) = 0.0_rp
      do is = 1, nspec
         alpha_in(:, is) = alpha_global
      end do
   end subroutine build_strux_inputs

   module subroutine structb_strux(this)
      class(lattice), intent(inout) :: this

      integer, parameter :: max_orb = 16
      real(rp), parameter :: match_tol = 1.0e-5_rp
      integer :: ii, ia, ib, ja, jb, m, is, js, nspec, nl, nl2, sbar_dim, pair_idx, nt_store
      integer :: nbas, nttab
      integer :: label, species_idx
      integer, allocatable :: ips(:), lmxb(:), orb_map(:), species_labels(:)
      real(rp) :: pair_cutoff, solve_cutoff, pair_cutoff_bohr, solve_cutoff_bohr, alat_bohr, wav_bohr
      real(rp) :: t_total_start, t_total_end, t_map_start, t_map_end, t_inputs_start, t_inputs_end, t_store_start, t_store_end, t_remap_end
      real(rp) :: alpha_debug(4)
      character(len=16) :: effective_screening
      real(rp) :: vec_target(3)
      real(rp), allocatable :: pos(:,:), cralat(:,:), rmt(:), alpha_in(:,:), hcr(:,:), alpha_site(:,:), &
         adot_site(:,:), alpha_l_out(:,:), tral(:,:,:), trad(:,:,:)
      type(strux_options) :: opts
      type(strux_result) :: result

      nl = this%control%lmax + 1
      nl2 = nl*nl
      nbas = this%nbas
      if (nbas <= 0) nbas = max(1, maxval(this%no(1:min(this%kk, size(this%no)))))
      sbar_dim = size(this%sbar, 1)
      pair_cutoff = sqrt(this%r2)
      solve_cutoff = sqrt(max(this%strux_solve_scale, 1.0_rp)*this%r2)
      pair_cutoff_bohr = pair_cutoff*ang2au
      solve_cutoff_bohr = solve_cutoff*ang2au
      alat_bohr = this%alat*ang2au
      wav_bohr = this%wav*ang2au

      if (nl2 > sbar_dim) then
         call g_logger%fatal('strux basis size exceeds allocated sbar dimension', __FILE__, __LINE__)
      end if
      if (nbas <= 0) then
         call g_logger%fatal('strux backend requires lattice basis coordinates', __FILE__, __LINE__)
      end if
      if (.not. allocated(this%crd) .or. .not. allocated(this%izp)) then
         call g_logger%fatal('strux backend requires primitive-cell coordinates and species labels', __FILE__, __LINE__)
      end if
      if (size(this%crd, 2) < nbas .or. size(this%izp) < nbas) then
         call g_logger%fatal('strux inferred basis size exceeds primitive-cell storage', __FILE__, __LINE__)
      end if
      if (any(this%izp(1:nbas) <= 0)) then
         call g_logger%fatal('strux backend found non-positive primitive species labels in lattice%izp', __FILE__, __LINE__)
      end if

      allocate(orb_map(max_orb))
      call build_orbital_map(sbar_dim, orb_map)

      opts%method = STRUX_METHOD_LMTO47
      opts%auto_alpha = .false.
      opts%want_sdot = this%strux_want_sdot
      opts%lmaxw = -1
      opts%pair_cutoff = pair_cutoff_bohr
      opts%solve_cutoff = solve_cutoff_bohr
      effective_screening = trim(lower(this%screening))
      select case (effective_screening)
      case ('manual', 'default')
         opts%screening_mode = STRUX_LMTO47_IALPHA_MANUAL
      case ('sigma')
         opts%screening_mode = STRUX_LMTO47_IALPHA_SIGMA
      case ('fitted')
         opts%screening_mode = STRUX_LMTO47_IALPHA_FITD
      case default
         call g_logger%fatal('Unsupported effective screening mode='//trim(effective_screening), __FILE__, __LINE__)
      end select

      if (effective_screening == 'sigma' .or. effective_screening == 'fitted') then
         opts%auto_alpha = .true.
      end if

      allocate(species_labels(nbas))
      species_labels = 0
      nspec = 0
      allocate(pos(3, nbas), ips(nbas))
      allocate(cralat(3, this%kk))

      pos(:, :) = this%crd(:, 1:nbas)
      cralat(:, :) = this%cr(:, 1:this%kk)*this%alat
      do ib = 1, nbas
         label = this%izp(ib)
         if (label <= 0) then
            call g_logger%fatal('strux backend found non-positive primitive species labels in lattice%izp', __FILE__, __LINE__)
         end if
         species_idx = 0
         do is = 1, nspec
            if (species_labels(is) == label) then
               species_idx = is
               exit
            end if
         end do
         if (species_idx == 0) then
            nspec = nspec + 1
            species_labels(nspec) = label
            species_idx = nspec
         end if
         ips(ib) = species_idx
      end do
      allocate(lmxb(nspec), rmt(nspec))
      allocate(alpha_in(0:nl - 1, nspec), hcr(nl, nspec))
      lmxb(:) = this%control%lmax
      call this%build_strux_inputs(nspec, nl, species_labels(1:nspec), alpha_in, hcr, rmt)
      write (17, '(a, i6, a, i6, a, f10.4, a, f10.4)') 'STRUX periodic solve nbas=', nbas, ' kk=', this%kk, ' pair_cutoff=', pair_cutoff, ' solve_cutoff=', solve_cutoff

      ! DEBUG_STRUX_DIAG_BEGIN (temporary instrumentation; safe to remove as one block)
      write (17, '(a)') 'DEBUG_STRUX_DIAG: begin'
      write (17, '(a,a)') 'DEBUG_STRUX_DIAG: effective_screening=', trim(effective_screening)
      write (17, '(a,i0)') 'DEBUG_STRUX_DIAG: opts%screening_mode=', opts%screening_mode
      write (17, '(a,l1)') 'DEBUG_STRUX_DIAG: opts%auto_alpha=', opts%auto_alpha
      write (17, '(a,f12.6,a,f12.6)') 'DEBUG_STRUX_DIAG: hcr[min,max]=', minval(hcr), ',', maxval(hcr)
      write (17, '(a,f12.6,a,f12.6)') 'DEBUG_STRUX_DIAG: alpha_in[min,max]=', minval(alpha_in), ',', maxval(alpha_in)
      write (17, '(a,f12.6,a,f12.6)') 'DEBUG_STRUX_DIAG: rmt[min,max]=', minval(rmt), ',', maxval(rmt)
      write (17, '(a,f12.6,a,f12.6)') 'DEBUG_STRUX_DIAG: avw_bohr, solve_cutoff_bohr=', wav_bohr, ',', solve_cutoff_bohr
      write (17, '(a)') 'DEBUG_STRUX_DIAG: end'
      flush(17)
      if (rank == 0) then
         call g_logger%info('DEBUG_STRUX_DIAG effective_screening='//trim(effective_screening)// &
                            ' mode='//int2str(opts%screening_mode)// &
                            ' auto='//merge('T','F',opts%auto_alpha), __FILE__, __LINE__)
      end if
      ! DEBUG_STRUX_DIAG_END

      if (effective_screening == 'manual' .or. effective_screening == 'default') then
         if (this%strux_want_sdot) then
            call g_logger%fatal('strux manual/default screening with sdot is unsupported; use sigma or fitted', __FILE__, __LINE__)
         end if
      end if

      if (allocated(result%iax)) deallocate(result%iax)
      if (allocated(result%alpha)) deallocate(result%alpha)
      if (allocated(result%alpha_l)) deallocate(result%alpha_l)
      if (allocated(result%s)) deallocate(result%s)
      if (allocated(result%sdot)) deallocate(result%sdot)

      select case (effective_screening)
      case ('manual', 'default')
         call strux_compute(opts, nbas, nspec, nl, alat_bohr, this%a, pos, ips, lmxb, wav_bohr, rmt, result, alpha_in=alpha_in)
         call strux_lmto47_screening(nbas, nspec, nl, wav_bohr, ips, lmxb, rmt, alpha_in, alpha_site, adot_site, &
            tral, trad, screening_mode=opts%screening_mode)
      case ('sigma')
         call strux_compute(opts, nbas, nspec, nl, alat_bohr, this%a, pos, ips, lmxb, wav_bohr, rmt, result, hcr=hcr)
         call strux_lmto47_autoalpha_screening(nbas, nspec, nl, wav_bohr, ips, lmxb, rmt, hcr, &
            alpha_l_out=alpha_l_out, alpha_out=alpha_site, adot=adot_site, tral=tral, trad=trad, &
            screening_mode=opts%screening_mode)
      case ('fitted')
         call strux_compute(opts, nbas, nspec, nl, alat_bohr, this%a, pos, ips, lmxb, wav_bohr, rmt, result)
         call strux_lmto47_autoalpha_screening(nbas, nspec, nl, wav_bohr, ips, lmxb, rmt, hcr, &
            alpha_l_out=alpha_l_out, alpha_out=alpha_site, adot=adot_site, tral=tral, trad=trad, &
            screening_mode=opts%screening_mode)
      end select
      do is = 1, nspec
         alpha_debug = 0.0_rp
         select case (effective_screening)
         case ('manual', 'default')
            alpha_debug(1:min(nl, 4)) = alpha_in(0:min(nl - 1, 3), is)
         case ('sigma', 'fitted')
            if (allocated(alpha_l_out)) then
               alpha_debug(1:min(nl, 4)) = alpha_l_out(0:min(nl - 1, 3), is)
            else if (allocated(result%alpha_l)) then
               alpha_debug(1:min(nl, 4)) = result%alpha_l(0:min(nl - 1, 3), is)
            end if
         end select
         write (17, '(a, a, a, i4, a, i6, a, 4f12.6)') 'STRUX screening=', trim(effective_screening), &
            ' species=', is, ' label=', species_labels(is), ' alpha=', alpha_debug
      end do
      flush(17)

      nttab = result%nttab
      write (17, '(a, i8, a, f10.3)') 'STRUX periodic result nttab=', nttab, ' cpu_s=', t_total_end - t_total_start

      do ii = 1, this%ntot
         this%alpha(:, :, ii) = 0.0_rp
         this%alpha_dot(:, :, ii) = 0.0_rp
         this%alpha(1:nl2, 1:nbas, ii) = alpha_site(:, 1:nbas)
         this%alpha_dot(1:nl2, 1:nbas, ii) = adot_site(:, 1:nbas)
      end do

      do is = 1, nspec
         label = species_labels(is)
         if (label < 1 .or. label > size(this%symbolic_atoms)) then
            call g_logger%fatal('strux backend found invalid species label for screening alpha storage', __FILE__, __LINE__)
         end if
         if (allocated(this%symbolic_atoms(label)%potential%screening_alpha)) then
            if (lbound(this%symbolic_atoms(label)%potential%screening_alpha, 1) /= 0 .or. &
                ubound(this%symbolic_atoms(label)%potential%screening_alpha, 1) /= nl - 1) then
               deallocate(this%symbolic_atoms(label)%potential%screening_alpha)
               allocate(this%symbolic_atoms(label)%potential%screening_alpha(0:nl - 1))
            end if
         else
            allocate(this%symbolic_atoms(label)%potential%screening_alpha(0:nl - 1))
         end if
         select case (effective_screening)
         case ('manual', 'default')
            this%symbolic_atoms(label)%potential%screening_alpha(0:nl - 1) = alpha_in(0:nl - 1, is)
         case ('sigma', 'fitted')
            if (allocated(alpha_l_out)) then
               this%symbolic_atoms(label)%potential%screening_alpha(0:nl - 1) = alpha_l_out(0:nl - 1, is)
            else if (allocated(result%alpha_l)) then
               this%symbolic_atoms(label)%potential%screening_alpha(0:nl - 1) = result%alpha_l(0:nl - 1, is)
            else
               call g_logger%fatal('strux backend did not produce screening alpha data for potential transform', __FILE__, __LINE__)
            end if
         end select
      end do
      do ii = 1, this%ntot
         ia = representative_atom_index(this, ii)
         ib = primitive_basis_label(this, ia)
         nt_store = this%kk
         call this%clusba(this%r2, cralat, ia, this%kk, this%kk, nt_store)
         call cpu_time(t_map_start)
         write (17, '(a, i5, a, i5, a, i5, a, i6)') 'STRUX map   center ', ii, ' atom ', ia, ' basis=', ib, ' nt=', nt_store
         call write_neighbor_vector_dump(17, this%sbarvec, nt_store)

         do m = 1, nt_store
            vec_target(:) = this%sbarvec(:, m)
            if (m == 1) then
               ja = ia
               jb = ib
            else
               ja = this%nn(ia, m)
               if (ja == 0) then
                  ja = find_neighbor_atom_by_vector(this, ia, vec_target, cralat, match_tol)
                  if (ja == 0) cycle
               end if
               jb = primitive_basis_label(this, ja)
            end if
            pair_idx = find_pair_by_vector(nttab, result%iax, this%a, pos, this%alat, ib, jb, vec_target, match_tol)
            if (pair_idx == 0) then
               write (17, '(a,3f14.8,a,2i6,a,2i6)') 'STRUX pair miss vec=', vec_target, ' center/neigh=', ia, ja, ' basis=', ib, jb
               call g_logger%fatal('Failed to map strux periodic pair to nn ordering', __FILE__, __LINE__)
            end if

            do is = 1, nl2
               if (orb_map(is) <= 0) cycle
               do js = 1, nl2
                  if (orb_map(js) <= 0) cycle
                  this%sbar(is, js, m, ii) = cmplx(result%s(orb_map(is), orb_map(js), pair_idx), 0.0_rp, kind=rp)
                  if (this%strux_want_sdot) then
                     this%sdot(is, js, m, ii) = cmplx(result%sdot(orb_map(is), orb_map(js), pair_idx), 0.0_rp, kind=rp)
                  end if
               end do
            end do
            call write_strux_block(this, nl2, m, ia, ja)
         end do
         call cpu_time(t_map_end)
         write (17, '(a, i5, a, f10.3, a, i6)') 'STRUX done  center ', ii, ' cpu_s=', t_map_end - t_map_start, ' stored_nn=', nt_store
         flush(17)
      end do

      if (allocated(alpha_l_out)) deallocate(alpha_l_out)
      if (allocated(tral)) deallocate(tral)
      if (allocated(trad)) deallocate(trad)
      if (allocated(alpha_site)) deallocate(alpha_site)
      if (allocated(adot_site)) deallocate(adot_site)
      deallocate(pos, cralat, ips, lmxb, rmt, alpha_in, hcr, orb_map, species_labels)
   end subroutine structb_strux

   module subroutine write_strux_block(this, nl2, m, ii, jclus)
      class(lattice), intent(in) :: this
      integer, intent(in) :: nl2, m, ii, jclus
      integer :: is, js

      ! write (*, '(" SBAR neighbor center=",i5," slot=",i5," iclus=",i5," vec=",3f12.6)') ii, m, iclus, &
      !    this%sbarvec(1, m), this%sbarvec(2, m), this%sbarvec(3, m)
      write (16, '(" VECTOR=",3f12.6,"   ICLUS=",i5," JCLUS=",i5)') this%sbarvec(1, m), this%sbarvec(2, m), this%sbarvec(3, m), ii, jclus
      if (this%strux_want_sdot) then
         write (15, '(" VECTOR=",3f12.6,"   ICLUS=",i5," JCLUS=",i5)') this%sbarvec(1, m), this%sbarvec(2, m), this%sbarvec(3, m), ii, jclus
      end if

      do is = 1, nl2
         do js = 1, nl2
            write (13) this%sbar(is, js, m, ii)
            if (this%strux_want_sdot) write (14) this%sdot(is, js, m, ii)
         end do
         write (16, '(*(f10.4))') (real(this%sbar(is, js, m, ii), rp), js=1, nl2)
         if (this%strux_want_sdot) then
            write (15, '(*(f10.4))') (real(this%sdot(is, js, m, ii), rp), js=1, nl2)
         end if
      end do
      flush(16)
      flush(17)
      if (this%strux_want_sdot) flush(15)
   end subroutine write_strux_block

   module subroutine write_neighbor_vector_dump(iunit, sbarvec, nt)
      integer, intent(in) :: iunit, nt
      real(rp), intent(in) :: sbarvec(:, :)
      integer :: i

      write (iunit, '(i5)') nt
      do i = 1, nt
         write (iunit, '(3f8.4)') sbarvec(1, i), sbarvec(2, i), sbarvec(3, i)
      end do
   end subroutine write_neighbor_vector_dump

   module integer function representative_atom_index(this, ii)
      class(lattice), intent(in) :: this
      integer, intent(in) :: ii

      representative_atom_index = 0

      if (allocated(this%iu)) then
         if (ii <= size(this%iu)) then
            if (this%iu(ii) > 0) representative_atom_index = this%iu(ii)
         end if
      end if
      if (representative_atom_index == 0 .and. allocated(this%ib)) then
         if (ii <= size(this%ib)) then
            if (this%ib(ii) > 0) representative_atom_index = this%ib(ii)
         end if
      end if
      if (representative_atom_index == 0 .and. allocated(this%irec)) then
         if (ii <= size(this%irec)) then
            if (this%irec(ii) > 0) representative_atom_index = this%irec(ii)
         end if
      end if
      if (representative_atom_index == 0) then
         if (ii <= this%kk) representative_atom_index = ii
      end if
      if (representative_atom_index <= 0 .or. representative_atom_index > this%kk) then
         call g_logger%fatal('Failed to resolve representative atom index for strux mapping', __FILE__, __LINE__)
      end if
   end function representative_atom_index

   module integer function primitive_basis_label(this, ia)
      class(lattice), intent(in) :: this
      integer, intent(in) :: ia

      primitive_basis_label = 0
      ! For strux_lib mapping, primitive labels must be species labels.
      ! Use izp first so this remains consistent with legacy setups where
      ! lattice%no may encode prototype classes (not species IDs).
      if (allocated(this%izp)) then
         if (ia >= 1 .and. ia <= size(this%izp)) then
            primitive_basis_label = this%izp(ia)
         end if
      end if
      if (primitive_basis_label <= 0 .and. allocated(this%no)) then
         if (ia >= 1 .and. ia <= size(this%no)) then
            primitive_basis_label = this%no(ia)
         end if
      end if
      if (primitive_basis_label <= 0 .and. allocated(this%num)) then
         if (ia >= 1 .and. ia <= size(this%num)) then
            primitive_basis_label = this%num(ia)
         end if
      end if
      if (primitive_basis_label <= 0) then
         call g_logger%fatal('Failed to resolve primitive basis label for strux mapping', __FILE__, __LINE__)
      end if
   end function primitive_basis_label

   module integer function find_pair_by_vector(nttab, iax, plat, pos, alat, ib, jb, vec_target, tol)
      integer, intent(in) :: nttab, ib, jb
      integer, intent(in) :: iax(:,:)
      real(rp), intent(in) :: plat(3, 3), pos(:, :), alat, vec_target(3), tol
      integer :: i, n1, n2, n3
      real(rp) :: vec(3)
      real(rp) :: dmax, best_d
      integer :: best_i

      find_pair_by_vector = 0
      best_d = huge(1.0_rp)
      best_i = 0
      do i = 1, nttab
         if (iax(1, i) /= ib) cycle
         if (iax(2, i) /= jb) cycle
         n1 = iax(3, i)
         n2 = iax(4, i)
         n3 = iax(5, i)
         vec(:) = alat*(pos(:, jb) - pos(:, ib) + real(n1, rp)*plat(:, 1) + real(n2, rp)*plat(:, 2) + real(n3, rp)*plat(:, 3))
         dmax = maxval(abs(vec - vec_target))
         if (dmax < best_d) then
            best_d = dmax
            best_i = i
         end if
         if (dmax <= tol) then
            find_pair_by_vector = i
            return
         end if
      end do
      if (best_i /= 0) find_pair_by_vector = best_i
   end function find_pair_by_vector

   module integer function find_neighbor_atom_by_vector(this, ia, vec_target, cralat, tol)
      class(lattice), intent(in) :: this
      integer, intent(in) :: ia
      real(rp), intent(in) :: vec_target(:), cralat(:, :)
      real(rp), intent(in) :: tol
      integer :: ja
      real(rp) :: dv(3)

      find_neighbor_atom_by_vector = 0
      do ja = 1, this%kk
         if (ja == ia) cycle
         dv(:) = cralat(:, ja) - cralat(:, ia)
         if (maxval(abs(dv - vec_target)) <= tol) then
            find_neighbor_atom_by_vector = ja
            return
         end if
      end do
   end function find_neighbor_atom_by_vector

   module subroutine build_orbital_map(norb, orb_map)
      integer, intent(in) :: norb
      integer, intent(out) :: orb_map(16)
      integer :: i
      integer, parameter :: full_map(16) = [ &
         1, 4, 2, 3, 5, 6, 8, 9, 7, 13, 14, 12, 15, 11, 16, 10 ]

      orb_map = 0
      do i = 1, min(norb, size(full_map))
         orb_map(i) = full_map(i)
      end do
   end subroutine build_orbital_map

   module subroutine dbar1(this, ia, r2, wav, crd, nat, ndi, np, nr, ii)
      implicit none
      class(lattice), intent(inout) :: this
      ! Inputs
      integer, intent(in) :: ia, nat, ndi, np
      integer, intent(in) :: nr, ii
      real(rp), intent(in) :: r2, wav
      real(rp), dimension(3, ndi), intent(in) :: crd
      ! Local Scalars
      integer :: i, j, k, m, na, nrl, nt, jclus_dbg
      real(rp), dimension(:), allocatable :: bet, wk
      real(rp), dimension(:), allocatable :: a
      real(rp), dimension(:, :), allocatable :: cr
      real(rp), dimension(:, :), allocatable :: s
      real(rp), dimension(:, :, :), allocatable :: sbar
      real(rp), dimension(:, :), allocatable :: sbarvec
      !
      ! External Calls
      !external CLUSBA, MICHA

      nt = this%kk ! Use cluster size instead of fixed 500 for neighbours
      ! allocate (cr(3, nt))
      allocate (sbar(np, np, nt))
      allocate(sbarvec(3, nt))
      call this%clusba(r2, crd, ia, nat, ndi, nt, sbarvec)
      write (17, 10000) nt
      write (17, 10002) ((sbarvec(j, i), j=1, 3), i=1, nt)
      !call write_neighbor_vector_dump(17, this%sbarvec, nt)
      !write (17, 10001) ((this%sbarvec(j, i), j=1, 3), i=1, nt)
      nrl = np*nt
      na = (nrl*(nrl + 1))/2
      allocate (a(na))
      allocate (bet(nrl))
      allocate (wk(nrl))
      allocate (s(nrl, nrl))
      call micha(wav, sbarvec, nt, np, nrl, na, sbar, a, wk, bet, s, ia, r2)
      !call micha(wav, this%sbarvec, nt, np, nrl, na, sbar, a, wk, bet, s, ia, r2)

      ! Saving parameters to be used in the Hamiltonian build
      ! Call clusba again for proper number of TB neighbours
      call this%clusba((r2/9.0d0), crd, ia, nat, ndi, nt, sbarvec)
      ! Store the number of neighbours
      this%nn_max = nt

      do m = 1, nt
         do i = 1, np
            do j = 1, np
               this%sbar(i, j, m, ii) = sbar(i, j, m)
            end do
         end do
         jclus_dbg = ia
         if (m > 1 .and. allocated(this%nn)) then
            if (ia >= 1 .and. ia <= size(this%nn, 1) .and. m <= size(this%nn, 2)) then
               if (this%nn(ia, m) > 0) jclus_dbg = this%nn(ia, m)
            end if
         end if
         write (16, '(" VECTOR=",3f12.6,"   ICLUS=",i5," JCLUS=",i5)') this%sbarvec(1, m), this%sbarvec(2, m), this%sbarvec(3, m), ia, jclus_dbg
         do i = 1, np
            write (16, '(*(f10.4))') (real(this%sbar(i, j, m, ii), rp), j=1, np)
         end do
      end do
      flush(17)

      !do m=1, nt
      !  write(*, *) ia, m
      !  do i=1, 9
      !    write(*, ´(9f10.4)´)((real(this%sbar(i, j, m, ia))), j=1, 9)
      !  end do
      !end do
      deallocate (a, bet, wk, s, sbar, sbarvec)
      return

10000 format(i5)
10001 format(3f8.4)
10002 format(3f8.4)
   end subroutine dbar1

   module subroutine clusba(this, r2, crd, ia, nat, ndi, n, sbarvec_out)
      implicit none
      class(lattice), intent(inout) :: this
      ! Inputs
      integer, intent(in) :: ia, nat, ndi
      real(rp), intent(in) :: r2
      real(rp), dimension(3, ndi), intent(in) :: crd
      ! Output
      integer, intent(inout) :: n
      real(rp), dimension(:, :), intent(inout), optional :: sbarvec_out
      ! Local variables
      integer :: i, ii, k, nn
      real(rp) :: s1
      real(rp), dimension(3) :: dum

#ifdef USE_SAFE_ALLOC
      if (allocated(this%sbarvec)) call g_safe_alloc%deallocate('lattice.sbarvec', this%sbarvec)
      call g_safe_alloc%allocate('lattice.sbarvec', this%sbarvec, (/3, this%kk/))
#else
      if (allocated(this%sbarvec)) deallocate (this%sbarvec)
      allocate (this%sbarvec(3, this%kk))
#endif

      this%sbarvec(:, :) = 0.0d0
      if (present(sbarvec_out)) sbarvec_out(:, :) = 0.0d0

      ii = 1
      do k = 1, 3
         this%sbarvec(k, 1) = 0.0d0
      end do
      if (present(sbarvec_out)) sbarvec_out(:, 1) = 0.0d0

      do nn = 1, nat
         s1 = 0.0
         do i = 1, 3
            dum(i) = (crd(i, nn) - crd(i, ia))**2
            s1 = s1 + dum(i)
         end do
         if (s1 < r2 .and. s1 > 0.0001) then
            ii = ii + 1
            this%sbarvec(1, ii) = crd(1, nn) - crd(1, ia)
            this%sbarvec(2, ii) = crd(2, nn) - crd(2, ia)
            this%sbarvec(3, ii) = crd(3, nn) - crd(3, ia)
            if (present(sbarvec_out)) then
               if (ii <= size(sbarvec_out, 2)) then
                  sbarvec_out(:, ii) = this%sbarvec(:, ii)
               end if
            end if
         end if
      end do
      n = ii
   end subroutine clusba

   module subroutine micha(rws, r, nr, nlm, nrl, na, sbar, a, wk, bet, s, iclus, r2)
      implicit none
      ! Inputs              .
      integer, intent(in) :: iclus, na, nlm, nr, nrl
      real(rp), intent(in) :: rws, r2
      real(rp), dimension(3, nr), intent(in) :: r
      ! Outputs
      real(rp), dimension(na), intent(inout) :: a
      real(rp), dimension(nrl), intent(inout) :: bet, wk
      real(rp), dimension(nrl, nrl), intent(inout) :: s
      real(rp), dimension(nlm, nlm, nr), intent(inout) :: sbar
      ! Local variables
      real(rp) :: fak, pi
      real(rp), dimension(4) :: q
      ! External Calls
      !external SHLDCH, STREZE
      ! Intrinsic Functions
      intrinsic ATAN

      pi = 4.d0*ATAN(1.d0)
      ! ----------------------------------
      call STREZE(rws, r, nr, s, nrl, nlm)
      ! -------------------------------------------
      fak = 2.d0
      !Original faktors
      q(1) = 0.3485d0*fak
      q(2) = 0.05303d0*fak
      q(3) = 0.010714d0*fak
      q(4) = 0.00337d0*fak   ! f-channel screening parameter
      ! Factors from LMTO47
      !q(1) = 0.33727d0 * fak
      !q(2) = 0.05115d0 * fak
      !q(3) = 0.01397d0 * fak
  !! ! New from LMTO47 (fcc Pt/Es)    data QM/0.3500000d0, 0.0571667d0, 0.0168070d0/
      !q(1) = 0.3500000d0* fak
      !q(2) = 0.0571667d0* fak
      !q(3) = 0.0168070d0* fak
      !     NA=(NRL*(NRL+1))/2
      call SHLDCH(r, nr, nlm, nrl, s, a, na, q, bet, wk, sbar, iclus, r2)
   end subroutine micha

   module subroutine STREZE(w, r, nr, s, nrl, nlm)
      implicit none
      ! Input
      integer, intent(in) :: nlm, nr, nrl
      real(rp), intent(in) :: w
      real(rp), dimension(3, nr), intent(in) :: r
      ! Output
      real(rp), dimension(nrl, nrl), intent(inout) :: s
      ! Local variables
      integer :: ilm, ir, irl0, jlm, jr, jrl0
      real(rp) :: rr, w1
      real(rp), dimension(3) :: dr
      real(rp), dimension(16, 16) :: s0
      ! External calls
      !external CANSO
      ! Intrinsic Functions
      intrinsic SQRT
      do ir = 1, nr
         irl0 = (ir - 1)*nlm
         do jr = 1, nr
            jrl0 = (jr - 1)*nlm
            dr(1) = (r(1, jr) - r(1, ir))/w
            dr(2) = (r(2, jr) - r(2, ir))/w
            dr(3) = (r(3, jr) - r(3, ir))/w
            rr = SQRT(dr(1)**2 + dr(2)**2 + dr(3)**2)
            w1 = 1.d0
            call CANSO(w1, dr, s0)
            do jlm = 1, nlm
               do ilm = 1, nlm
                  s(ilm + irl0, jlm + jrl0) = s0(ilm, jlm)
               end do
            end do
         end do
      end do
   end subroutine streze

   module subroutine shldch(r, nr, nlm, nrl, s, a, na, q, bet, wk, sbar, iclus, r2)
      implicit none
      !parameter for cutoff of sbar construction
      integer, parameter :: ncut = 9
      ! Input
      integer, intent(in) :: iclus, na, nlm, nr, nrl
      real(rp), intent(in) ::r2
      real(rp), dimension(4), intent(in) :: q
      real(rp), dimension(3, nr), intent(in) :: r
      ! Output
      real(rp), dimension(na), intent(inout) :: a
      real(rp), dimension(nrl), intent(inout) :: bet, wk
      real(rp), dimension(nrl, nrl), intent(inout) :: s
      real(rp), dimension(nlm, nlm, nr), intent(inout) :: sbar
      ! Local variables
      integer :: i, ia, ilm, ir, irl, irl0, isb, j, jlm, jsb, l, l2, lmax, m, ndef, hitc, info
      real(rp), dimension(:, :), allocatable :: s_temp
      ! External Calls
      !external chlr2f, chlr2s
      ! External Functions  .
      !integer, external :: LL
      ndef = 0
      lmax = LL(nlm)
      allocate (s_temp(nrl, nrl))
      write (17, 10000) lmax, q
      irl = 0
      do ir = 1, nr
         do l = 0, lmax
            do m = 1, 2*l + 1
               irl = irl + 1
               bet(irl) = 1.d0/q(l + 1)
            end do
         end do
      end do
      s_temp = s
      do i = 1, nrl
         s_temp(i, i) = s_temp(i, i) + bet(i)
      end do
      call DPOTRF('U', nrl, s_temp, nrl, info)
      write (17, 10001) ndef
      call DPOTRS('U', nrl, nlm, s_temp, nrl, s, nrl, INFO)
      deallocate (s_temp)
      do ilm = 1, nlm
         do irl = 1, nrl
            s(irl, ilm) = -bet(irl)*s(irl, ilm)
         end do
      end do
      ! --------------------------------
      ir = 0
      hitc = 0
      !print *, ´ nr = ´, nr
      do ir = 1, nr
         if (abs(R(1, IR)**2 + R(2, IR)**2 + R(3, IR)**2 - &
                 R(1, 1)**2 + R(2, 1)**2 + R(3, 1)**2) <= (R2/ncut)) then
            ! Legacy view.sbar printing is emitted in dbar1 where neighbour
            ! atom indices are available from nn(ia,m).
            hitc = hitc + 1
            irl0 = (ir - 1)*nlm
            do ilm = 1, nlm
               do jlm = 1, nlm
                  sbar(ilm, jlm, hitc) = 2.*s(ilm + irl0, jlm)
               end do
            end do
          !!! Changed for ´symmetrizing´ Sbar
          !! NOT NEEDED SINCE WITH LARGER CUTOFF
          !! SBAR IS CALCULATED PROPERLY
            do isb = 1, nlm
               write (13) (sbar(isb, jsb, hitc), jsb=1, nlm)
            end do
            ! 208 FORMAT(5F12.6:/(12X, 4F12.6))
            do isb = 1, nlm
               ! Legacy view.sbar printing is emitted in dbar1.
            end do
         end if
      end do
      !print *, ´ hitc = ´, hitc
      return
      !
      ! ... Format Declarations ...
      !
10000 format(" LMAX=", i2, "   Q=", 4f10.6)
10001 format(" NDEF=", i10)
10002 format(" VECTOR=", 3f12.6, "   ICLUS=", i5)
10003 format(9f10.4)
   end subroutine shldch

   module subroutine canso(w, dr, sc)
      implicit none
      ! Input
      real(rp), intent(in) :: w
      real(rp), dimension(3), intent(in) :: dr
      ! Output
      real(rp), dimension(16, 16), intent(out) :: sc
      ! Local variables
      integer :: i, j, l, ll
      real(rp) :: el, el2, elem, elen, em, em2, emen, en, en2, r1, r2, r3, rr, s2, s3, s4, s5, s6, s7, &
                  sbyr, sq3, sq5, sq7
      integer, dimension(16) :: ip
      real(rp), dimension(16, 16) :: s
      ! Intrinsic Functions
      intrinsic SQRT
      !.. Data Declarations ..
      ! original and correct
      !     S=1, X=2, Y=3, Z=4, XY=5, YZ=6, ZX=7, X**2-Y**2=8, 3Z*Z-R*R=9
      !     f-orbitals (10-16): fz3, fxz2, fyz2, fz(x2-y2), fxyz, fx3, fy3
      data ip/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16/
      !     S=1, X=2, Y=3, Z=4, XY=5, YZ=6, ZX=7, X**2-Y**2=8, 3Z*Z-R*R=9
      !     (10-16 reserved for f-orbitals)
      r1 = dr(1)
      r2 = dr(2)
      r3 = dr(3)
      rr = SQRT(r1*r1 + r2*r2 + r3*r3)
      do i = 1, 16
         do j = 1, 16
            sc(i, j) = 0.0d0
         end do
      end do
      if (rr/w <= 0.30d0) return
      sbyr = w/rr
      s2 = sbyr*sbyr
      s3 = s2*sbyr
      s4 = s3*sbyr
      s5 = s4*sbyr
      s6 = s5*sbyr
      s7 = s6*sbyr
      sq3 = SQRT(3.d0)
      sq5 = SQRT(5.d0)
      sq7 = SQRT(7.d0)
      el = r1/rr
      em = r2/rr
      en = r3/rr
      el2 = el*el
      em2 = em*em
      en2 = en*en
      elem = el*em
      elen = el*en
      emen = em*en
      !
      sc(1, 1) = -2.0d0*sbyr
      sc(1, 2) = el*s2*2.d0*sq3
      sc(1, 3) = em*s2*2.d0*sq3
      sc(1, 4) = en*s2*2.d0*sq3
      sc(1, 5) = -2.d0*sq3*sq5*elem*s3
      sc(1, 6) = -2.d0*sq3*sq5*emen*s3
      sc(1, 7) = -2.d0*sq3*sq5*elen*s3
      sc(1, 8) = -sq3*sq5*s3*(el2 - em2)
      sc(1, 9) = sq5*s3*(1.d0 - 3.d0*en2)
      !-----------------------------------------------------------------------
      sc(2, 2) = (3.0d0*el2 - 1.0)*6.d0*s3
      sc(2, 3) = 18.0d0*s3*elem
      sc(2, 4) = 18.0d0*s3*elen
      sc(2, 5) = 6.d0*sq5*s4*em*(1.0d0 - 5.0d0*el2)
      sc(2, 6) = -30.0d0*sq5*s4*elem*en
      sc(2, 7) = 6.d0*sq5*s4*en*(1.0d0 - 5.0d0*el2)
      sc(2, 8) = 6.0d0*sq5*s4*el*(1.0d0 - 2.5d0*el2 + 2.5d0*em2)
      sc(2, 9) = 3.d0*sq3*sq5*s4*el*(1.0d0 - 5.0d0*en2)
      !-----------------------------------------------------------------------
      sc(3, 3) = 6.d0*s3*(3.0d0*em2 - 1.d0)
      sc(3, 4) = 18.d0*s3*emen
      sc(3, 5) = 6.0d0*sq5*s4*el*(1.d0 - 5.d0*em2)
      sc(3, 6) = 6.d0*sq5*s4*en*(1.d0 - 5.d0*em2)
      sc(3, 7) = sc(2, 6)
      sc(3, 8) = -6.d0*sq5*s4*em*(1.d0 - 2.5d0*em2 + 2.5d0*el2)
      sc(3, 9) = 3.d0*sq3*sq5*s4*em*(1.d0 - 5.d0*en2)
      !----------------------------------------------------------------------
      sc(4, 4) = 6.d0*s3*(3.d0*en2 - 1.d0)
      sc(4, 5) = sc(2, 6)
      sc(4, 6) = 6.d0*sq5*s4*em*(1.d0 - 5.d0*en2)
      sc(4, 7) = 6.d0*sq5*s4*el*(1.d0 - 5.d0*en2)
      sc(4, 8) = -15.d0*sq5*s4*en*(el*el - em2)
      sc(4, 9) = 3.d0*sq3*sq5*s4*en*(3.0d0 - 5.0d0*en2)
      !-----------------------------------------------------------------------
      sc(5, 5) = 10.d0*s5*(-35.d0*el2*em2 - 5.d0*en2 + 4.d0)
      sc(5, 6) = -50.d0*s5*elen*(7.d0*em2 - 1.d0)
      sc(5, 7) = -50.d0*s5*emen*(7.d0*el2 - 1.d0)
      sc(5, 8) = -175.d0*s5*elem*(el2 - em2)
      sc(5, 9) = -25.d0*sq3*s5*elem*(7.d0*en2 - 1.d0)
      !----------------------------------------------------------------------
      sc(6, 6) = 10.d0*s5*(-35.d0*em2*en2 - 5.d0*el2 + 4.d0)
      sc(6, 7) = -50.d0*s5*elem*(7.d0*en2 - 1.d0)
      sc(6, 8) = 50.d0*s5*emen*(3.5d0*em2 - 3.5d0*el2 - 1.d0)
      sc(6, 9) = -25.d0*sq3*s5*emen*(7.d0*en2 - 3.d0)
      !-----------------------------------------------------------------------
      sc(7, 7) = 10.d0*s5*(-35.d0*el2*en2 - 5.d0*em2 + 4.d0)
      sc(7, 8) = -50.d0*s5*elen*(3.5d0*el2 - 3.5d0*em2 - 1.0d0)
      sc(7, 9) = -25.d0*sq3*s5*elen*(7.d0*en2 - 3.d0)
      !-----------------------------------------------------------------------
      sc(8, 8) = 10.d0*s5*(-8.75d0*(el2 - em2)**2 - 5.d0*en2 + 4.d0)
      sc(8, 9) = -12.5d0*sq3*s5*(7.d0*en2 - 1.d0)*(el2 - em2)
      !-----------------------------------------------------------------------
      sc(9, 9) = -7.5d0*s5*(35.d0*en2*en2 - 30.d0*en2 + 3.d0)
      !-----------------------------------------------------------------------
      ! F-ORBITAL BLOCK (indices 10-16)
      ! From Methfessel, Mossner & Springborg, J. Phys. C 20, 1069 (1987), Table 1
      ! Orbital order: ilm 10=fz3, 11=fxz2, 12=fyz2, 13=fz(x2-y2), 14=fxyz, 15=fx3, 16=fy3
      !
      ! s-f interactions (row 1, cols 10-16)
      sc(1, 10) = sq7*s3*(5.d0*en2*en2 - 3.d0*en2)
      sc(1, 11) = -2.d0*sq7*sq5*elen*s4*(1.d0 - 5.d0*en2)
      sc(1, 12) = -2.d0*sq7*sq5*emen*s4*(1.d0 - 5.d0*en2)
      sc(1, 13) = -2.d0*sq7*sq5*(el2 - em2)*en*s4
      sc(1, 14) = -4.d0*sq7*sq5*elem*en*s4
      sc(1, 15) = sq7*s4*el*(5.d0*el2 - 3.d0)
      sc(1, 16) = sq7*s4*em*(5.d0*em2 - 3.d0)
      !
      ! p-f interactions (rows 2-4, cols 10-16)
      sc(2, 10) = 6.d0*sq7*sq3*s4*el*(5.d0*en2*en2 - 1.d0)
      sc(2, 11) = 3.d0*sq7*sq5*s5*em*(35.d0*en2*en2 - 30.d0*en2 + 3.d0)
      sc(2, 12) = -30.d0*sq7*sq5*s5*elem*en*(en2 - 1.d0)
      sc(2, 13) = 3.d0*sq7*sq5*s5*el*(7.d0*en2 - 1.d0)*(el2 - em2)
      sc(2, 14) = 6.d0*sq7*sq5*s5*elem*en*(7.d0*en2 - 3.d0)
      sc(2, 15) = 3.d0*sq7*s5*el*(35.d0*el2*em2 - 5.d0*el2 - 20.d0*em2 + 4.d0)
      sc(2, 16) = 3.d0*sq7*s5*el2*em*(35.d0*em2 - 15.d0)
      !
      sc(3, 10) = 6.d0*sq7*sq3*s4*em*(5.d0*en2*en2 - 1.d0)
      sc(3, 11) = -30.d0*sq7*sq5*s5*elem*en*(en2 - 1.d0)
      sc(3, 12) = 3.d0*sq7*sq5*s5*el*(35.d0*en2*en2 - 30.d0*en2 + 3.d0)
      sc(3, 13) = 3.d0*sq7*sq5*s5*em*(7.d0*en2 - 1.d0)*(el2 - em2)
      sc(3, 14) = 6.d0*sq7*sq5*s5*emen*en*(7.d0*en2 - 3.d0)
      sc(3, 15) = 3.d0*sq7*s5*el*em2*(35.d0*el2 - 15.d0)
      sc(3, 16) = 3.d0*sq7*s5*em*(35.d0*em2*el2 - 20.d0*el2 - 5.d0*em2 + 4.d0)
      !
      sc(4, 10) = 3.d0*sq7*sq5*s4*en*(35.d0*en2 - 28.d0)
      sc(4, 11) = 6.d0*sq7*sq5*s5*elen*en*(7.d0*en2 - 3.d0)
      sc(4, 12) = 6.d0*sq7*sq5*s5*emen*en*(7.d0*en2 - 3.d0)
      sc(4, 13) = -12.d0*sq7*sq5*s5*(el2 - em2)*en*(7.d0*en2 - 1.d0)
      sc(4, 14) = -12.d0*sq7*sq5*s5*elem*en*(7.d0*en2 - 1.d0)
      sc(4, 15) = 3.d0*sq7*s5*el*(5.d0*el2 - 1.d0)*(7.d0*en2 - 3.d0)
      sc(4, 16) = 3.d0*sq7*s5*em*(5.d0*em2 - 1.d0)*(7.d0*en2 - 3.d0)
      !
      ! d-f interactions (rows 5-9, cols 10-16)
      sc(5, 10) = -15.d0*sq7*sq5*s5*elem*(7.d0*en2*en2 - 6.d0*en2 + 1.d0)
      sc(5, 11) = -35.d0*sq7*s6*emen*elem*(1.d0 - 7.d0*en2)
      sc(5, 12) = -35.d0*sq7*s6*el2*en*elem*(1.d0 - 7.d0*en2)
      sc(5, 13) = -35.d0*sq7*s6*elem*(el2 - em2)*(7.d0*en2 - 1.d0)
      sc(5, 14) = -70.d0*sq7*s6*elem*(el2*en2 - em2*(1.d0 - en2))
      sc(5, 15) = -15.d0*sq7*s6*elem*(el2 - 3.d0*em2)*(7.d0*en2 - 1.d0)
      sc(5, 16) = -15.d0*sq7*s6*elem*(em2 - 3.d0*el2)*(7.d0*en2 - 1.d0)
      !
      sc(6, 10) = -15.d0*sq7*sq5*s5*emen*(7.d0*en2*en2 - 6.d0*en2 + 1.d0)
      sc(6, 11) = -35.d0*sq7*s6*el2*en*emen*(1.d0 - 7.d0*en2)
      sc(6, 12) = -35.d0*sq7*s6*em2*en*emen*(1.d0 - 7.d0*en2)
      sc(6, 13) = 35.d0*sq7*s6*emen*(el2 - em2)*(7.d0*en2 - 1.d0)
      sc(6, 14) = -70.d0*sq7*s6*emen*el2*(7.d0*en2 - 1.d0)
      sc(6, 15) = 15.d0*sq7*s6*emen*el*(em2 - 7.d0*el2)
      sc(6, 16) = -15.d0*sq7*s6*emen*em*(el2 - 7.d0*em2)
      !
      sc(7, 10) = -15.d0*sq7*sq5*s5*elen*(7.d0*en2*en2 - 6.d0*en2 + 1.d0)
      sc(7, 11) = -35.d0*sq7*s6*el2*en*elen*(1.d0 - 7.d0*en2)
      sc(7, 12) = 35.d0*sq7*s6*em2*elen*(1.d0 - 7.d0*en2)
      sc(7, 13) = 35.d0*sq7*s6*elen*(el2 - em2)*(7.d0*en2 - 1.d0)
      sc(7, 14) = -70.d0*sq7*s6*elen*em2*(7.d0*en2 - 1.d0)
      sc(7, 15) = 15.d0*sq7*s6*elen*el*(em2 - 7.d0*el2)
      sc(7, 16) = -15.d0*sq7*s6*elen*em*(el2 - 7.d0*em2)
      !
      sc(8, 10) = -15.d0*sq7*sq5*s5*(el2 - em2)*(7.d0*en2*en2 - 6.d0*en2 + 1.d0)
      sc(8, 11) = 70.d0*sq7*s6*en*(el2 - em2)*elen
      sc(8, 12) = 70.d0*sq7*s6*en*(el2 - em2)*emen
      sc(8, 13) = -70.d0*sq7*s6*(el2 - em2)*(el2*7.d0*en2 - el2 - 7.d0*em2*en2 + em2)
      sc(8, 14) = 140.d0*sq7*s6*elem*(el2 - em2)*(3.d0 - 7.d0*en2)
      sc(8, 15) = 35.d0*sq7*s6*en*el*(el2 - em2)*(7.d0*el2 - 9.d0)
      sc(8, 16) = -35.d0*sq7*s6*en*em*(el2 - em2)*(7.d0*em2 - 9.d0)
      !
      sc(9, 10) = 5.d0*sq7*sq5*s5*en*(21.d0*en2*en2 - 14.d0*en2 + 1.d0)
      sc(9, 11) = 35.d0*sq7*s6*elen*(21.d0*en2 - 5.d0)
      sc(9, 12) = 35.d0*sq7*s6*emen*(21.d0*en2 - 5.d0)
      sc(9, 13) = -70.d0*sq7*s6*en*(el2 - em2)*(21.d0*en2 - 5.d0)
      sc(9, 14) = -140.d0*sq7*s6*elem*en*(21.d0*en2 - 5.d0)
      sc(9, 15) = 35.d0*sq7*s6*el*(21.d0*en2 - 3.d0)*(7.d0*en2 - 1.d0)
      sc(9, 16) = 35.d0*sq7*s6*em*(21.d0*en2 - 3.d0)*(7.d0*en2 - 1.d0)
      !
      ! f-f interactions (rows 10-16, cols 10-16)
      ! Diagonal terms
      sc(10, 10) = -7.d0*s7*(99.d0*en2*en2*en2 - 135.d0*en2*en2 + 55.d0*en2 - 5.d0)
      sc(11, 11) = s7*(-385.d0*en2*en2*(el2 + em2) + 70.d0*en2*(el2 + em2) + 245.d0*el2*em2 + 5.d0)
      sc(12, 12) = sc(11, 11)  ! same for y
      sc(13, 13) = s7*(245.d0*en2*(el2 - em2)**2 - 70.d0*(el2 - em2)**2 - 385.d0*en2*en2*(el2 - em2)**2 + 5.d0)
      sc(14, 14) = s7*(-1225.d0*en2*en2*el2*em2 + 350.d0*en2*el2*em2 - 35.d0*el2*em2 + 5.d0)
      sc(15, 15) = s7*(-231.d0*el2*el2*em2 - 385.d0*el2*el2*en2*en2 + 70.d0*el2*el2*en2 - 55.d0*el2*el2 + &
                       245.d0*el2*em2*em2 + 5.d0)
      sc(16, 16) = s7*(-231.d0*em2*em2*el2 - 385.d0*em2*em2*en2*en2 + 70.d0*em2*em2*en2 - 55.d0*em2*em2 + &
                       245.d0*em2*el2*el2 + 5.d0)
      ! Off-diagonal f-f terms (selected important ones)
      sc(11, 12) = -35.d0*s7*elem*(35.d0*en2*en2 - 10.d0*en2 + 1.d0)
      sc(11, 13) = s7*el*(-245.d0*en2*en2*el2 + 245.d0*en2*en2*em2 + 70.d0*en2*el2 - 70.d0*en2*em2 - 5.d0*el2 + 5.d0*em2)
      sc(11, 14) = -70.d0*s7*elen*em*(35.d0*en2*en2 - 10.d0*en2 + 1.d0)
      sc(11, 15) = s7*el*em*(245.d0*en2*(el2 - em2) - 49.d0*(el2 - em2))
      sc(11, 16) = 70.d0*s7*em2*el*(35.d0*en2 - 5.d0)*(en2 - 1.d0)
      sc(12, 13) = s7*em*(-245.d0*en2*en2*el2 + 245.d0*en2*en2*em2 + 70.d0*en2*el2 - 70.d0*en2*em2 - 5.d0*el2 + 5.d0*em2)
      sc(12, 14) = -70.d0*s7*emen*el*(35.d0*en2*en2 - 10.d0*en2 + 1.d0)
      sc(12, 15) = 70.d0*s7*el2*em*(35.d0*en2 - 5.d0)*(en2 - 1.d0)
      sc(12, 16) = s7*el*em*(245.d0*en2*(el2 - em2) - 49.d0*(el2 - em2))
      sc(13, 14) = -140.d0*s7*elem*(el2 - em2)*(7.d0*en2*en2 - 5.d0*en2 + 1.d0)
      sc(13, 15) = -35.d0*s7*el*(el2 - em2)*(el2 - 7.d0*em2)*(7.d0*en2 - 1.d0)
      sc(13, 16) = -35.d0*s7*em*(el2 - em2)*(em2 - 7.d0*el2)*(7.d0*en2 - 1.d0)
      sc(14, 15) = -35.d0*s7*elem*em*(el2 - 3.d0*em2)*(7.d0*en2 - 1.d0)
      sc(14, 16) = -35.d0*s7*elem*el*(em2 - 3.d0*el2)*(7.d0*en2 - 1.d0)
      sc(15, 16) = 35.d0*s7*el*em*(el2 - em2)*(el2 + em2 - 5.d0*en2*(el2 + em2))
      !-----------------------------------------------------------------------
      do l = 2, 16
         ll = l - 1
         do j = 1, ll
            sc(l, j) = sc(j, l)
         end do
      end do
      do l = 1, 3
         sc(l + 1, 1) = -sc(l + 1, 1)
      end do
      do l = 5, 9
         do j = 2, 4
            sc(l, j) = -sc(l, j)
         end do
      end do
      do l = 10, 16
         do j = 2, 4
            sc(l, j) = -sc(l, j)
         end do
      end do
      ! ------ THIS PART CHANGES THE YLM ORDER AND MULTIPLIES BY -0.5 ----
      do i = 1, 16
         do j = 1, 16
            s(ip(j), ip(i)) = -0.5d0*sc(j, i)
         end do
      end do
      do i = 1, 16
         do j = 1, 16
            sc(j, i) = s(j, i)
         end do
      end do
   end subroutine canso

   module subroutine chlr2f(c, na, w, n, ndef)
      implicit none
      ! Input
      integer, intent(in) :: n, na
      ! Output
      integer, intent(out) :: ndef
      real(rp), dimension(n), intent(inout) :: w
      real(rp), dimension(na), intent(inout) :: c
      ! Local variables
      integer :: i, ic0, j, jc0, k
      real(rp) :: csum
      ! Intrinsic Functions
      intrinsic SQRT

      write (17, *) 'chol', na, n
      ic0 = 0
      do i = 1, n
         jc0 = 0
         do j = 1, i - 1
            csum = c(ic0 + j)
            do k = 1, j - 1
               csum = csum - w(k)*c(jc0 + k)
            end do
            jc0 = jc0 + j
            w(j) = csum/c(jc0)
         end do
         csum = c(ic0 + i)
         do k = 1, i - 1
            csum = csum - w(k)*w(k)
         end do
         if (csum <= 0.d0) then
            ndef = i - 1
            goto 1000
         else
            do j = 1, i - 1
               c(ic0 + j) = w(j)
            end do
            !
            ic0 = ic0 + i
            c(ic0) = SQRT(csum)
         end if
      end do
      ndef = n
1000  if (ndef < n) then
         write (17, 10000) n, ndef
      end if
      return
10000 format(" CHLR2F:    N=", i4, "    NDEF=", i4)
   end subroutine chlr2f

   module subroutine chlr2s(c, na, v, n, m)
      implicit none
      ! Input
      integer, intent(in) :: m, n, na
      real(rp), dimension(na), intent(in) :: c
      ! Output
      real(rp), dimension(n, m), intent(inout) :: v
      ! Local Variables
      integer :: i, ic0, ind1, ind2, k, mm
      real(rp) :: csum

      do mm = 1, m
         ic0 = 0
         do i = 1, n
            csum = v(i, mm)
            do k = 1, i - 1
               csum = csum - v(k, mm)*c(ic0 + k)
            end do
            ic0 = ic0 + i
            v(i, mm) = csum/c(ic0)
         end do
         !
         do i = n, 1, -1
            ind2 = (i*(i + 1))/2
            csum = v(i, mm)
            ind1 = ind2 + i
            do k = i + 1, n
               csum = csum - v(k, mm)*c(ind1)
               ind1 = ind1 + k
            end do
            v(i, mm) = csum/c(ind2)
         end do
      end do
   end subroutine chlr2s

   module function LL(ilm)
      implicit none
      ! Input
      integer, intent(in) :: ilm
      ! Function Declaration
      integer :: LL
      ! Local variables
      integer, dimension(100) :: lla
      ! Data Declarations
      data lla/0, 3*1, 5*2, 7*3, 9*4, 11*5, 13*6, 15*7, 17*8, 19*9/

      LL = lla(ilm)
   end function LL

   module subroutine outmap(IM, IZP, NN, NO, ND, NM, NTOT)
      implicit none
      ! Input
      integer, intent(in) :: IM, ND, NM, NTOT
      integer, dimension(ND), intent(in) :: NO
      integer, dimension(ND) :: IZP
      integer, dimension(ND, NM), intent(in) :: NN
      ! Local variables
      integer :: I, ID, IDM, J, K
      ! Intrinsic functions
      intrinsic MIN

      write (IM, 10000)
      do I = 1, NTOT
!!    K = NO(I)
         K = IZP(I)
         ID = NN(I, 1)
!!    IDM = MIN(ID, 21)
         IDM = MIN(ID, 16)
         write (IM, 10001) I, K, ID, (NN(I, J), J=2, IDM)
         if (IDM /= NN(I, 1)) then
            !     ID=NN(I, 1)
            write (IM, 10002) (NN(I, J), J=17, ID)
         end if
      end do
      return

10000 format( &
         " NEAREST NEIGHBOUR MAP"/, "      ATOM   TYPE  CONNECTIVITY", 5x, &
         "NEIGHBOURS")
10001 format(1x, 3i4, 2x, 20i5)
10002 format(29x, 16i5)
   end subroutine outmap

   module integer function mapa(I, J, R2, DD, CT)
      implicit none
      ! Input
      integer, intent(in) :: I, J
      real(rp), intent(in) :: R2
      real(rp) :: DD
      real(rp), dimension(50), intent(in) :: CT
      ! Local variables
      real(rp) :: CTM, CTSM

      CTM = (CT(1) + CT(1))/2.
      CTSM = CTM**2
      if (R2 >= CTSM) then
         MAPA = 0
      else
         MAPA = 1
      end if
   end function mapa

end submodule lattice_strux
