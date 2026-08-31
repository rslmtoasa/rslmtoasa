submodule (lattice_mod) lattice_lifecycle
   implicit none

contains

   !> @brief Construct a lattice object from an initialized control object.
   !> @details Wires the lattice to the run control state, restores default
   !>          geometry/structure-constant storage, and reads lattice input.
   !> @param[in] control_obj Control object that owns the input-file name and run options.
   !> @return Initialized lattice object.
   module function constructor(control_obj) result(obj)
      type(lattice) :: obj
      type(control), target, intent(in) :: control_obj

      obj%control => control_obj

      call obj%restore_to_default()
      call obj%build_from_file()
   end function constructor

   !> @brief Finalize a lattice object.
   !> @details Releases cluster, neighbor, surface, impurity, screening, and
   !>          structure-constant arrays owned by the object.
   !> @param[inout] this Lattice object being finalized.
   module subroutine destructor(this)
      type(lattice) :: this
#ifdef USE_SAFE_ALLOC
      if (allocated(this%ct)) call g_safe_alloc%deallocate('lattice.ct', this%ct)
      if (allocated(this%cr)) call g_safe_alloc%deallocate('lattice.cr', this%cr)
      if (allocated(this%crd)) call g_safe_alloc%deallocate('lattice.crd', this%crd)
      if (allocated(this%iu)) call g_safe_alloc%deallocate('lattice.iu', this%iu)
      if (allocated(this%ib)) call g_safe_alloc%deallocate('lattice.ib', this%ib)
      if (allocated(this%irec)) call g_safe_alloc%deallocate('lattice.irec', this%irec)
      if (allocated(this%izp)) call g_safe_alloc%deallocate('lattice.izp', this%izp)
      if (allocated(this%iz)) call g_safe_alloc%deallocate('lattice.iz', this%iz)
      if (allocated(this%num)) call g_safe_alloc%deallocate('lattice.num', this%num)
      if (allocated(this%no)) call g_safe_alloc%deallocate('lattice.no', this%no)
      if (allocated(this%z)) call g_safe_alloc%deallocate('lattice.z', this%z)
      if (allocated(this%izpsurf)) call g_safe_alloc%deallocate('lattice.izpsurf', this%izpsurf)
      if (allocated(this%izsurf)) call g_safe_alloc%deallocate('lattice.izsurf', this%izsurf)
      if (allocated(this%nosurf)) call g_safe_alloc%deallocate('lattice.nosurf', this%nosurf)
      if (allocated(this%crsurf)) call g_safe_alloc%deallocate('lattice.crsurf', this%crsurf)
      if (allocated(this%inclu)) call g_safe_alloc%deallocate('lattice.inclu', this%inclu)
      if (allocated(this%izpo)) call g_safe_alloc%deallocate('lattice.izpo', this%izpo)
      if (allocated(this%acr)) call g_safe_alloc%deallocate('lattice.acr', this%acr)
      if (allocated(this%atlist)) call g_safe_alloc%deallocate('lattice.atlist', this%atlist)
      if (allocated(this%nn)) call g_safe_alloc%deallocate('lattice.nn', this%nn)
      if (allocated(this%sbar)) call g_safe_alloc%deallocate('lattice.sbar', this%sbar)
      if (allocated(this%sdot)) call g_safe_alloc%deallocate('lattice.sdot', this%sdot)
      if (allocated(this%sbarvec)) call g_safe_alloc%deallocate('lattice.sbarvec', this%sbarvec)
      if (allocated(this%alpha)) call g_safe_alloc%deallocate('lattice.alpha', this%alpha)
      if (allocated(this%alpha_dot)) call g_safe_alloc%deallocate('lattice.alpha_dot', this%alpha_dot)
      if (allocated(this%screening_alpha)) call g_safe_alloc%deallocate('lattice.screening_alpha', this%screening_alpha)
      if (allocated(this%screening_sigma)) call g_safe_alloc%deallocate('lattice.screening_sigma', this%screening_sigma)
      if (allocated(this%ijpair)) call g_safe_alloc%deallocate('lattice.ijpair', this%ijpair)
      if (allocated(this%ijktrio)) call g_safe_alloc%deallocate('lattice.ijktrio', this%ijktrio)
      if (allocated(this%natoms_layer)) call g_safe_alloc%deallocate('lattice.natoms_layer', this%natoms_layer)
#else
      if (allocated(this%ct)) deallocate (this%ct)
      if (allocated(this%cr)) deallocate (this%cr)
      if (allocated(this%crd)) deallocate (this%crd)
      if (allocated(this%iu)) deallocate (this%iu)
      if (allocated(this%ib)) deallocate (this%ib)
      if (allocated(this%irec)) deallocate (this%irec)
      if (allocated(this%izp)) deallocate (this%izp)
      if (allocated(this%iz)) deallocate (this%iz)
      if (allocated(this%num)) deallocate (this%num)
      if (allocated(this%no)) deallocate (this%no)
      if (allocated(this%z)) deallocate (this%z)
      if (allocated(this%izpsurf)) deallocate (this%izpsurf)
      if (allocated(this%izsurf)) deallocate (this%izsurf)
      if (allocated(this%nosurf)) deallocate (this%nosurf)
      if (allocated(this%crsurf)) deallocate (this%crsurf)
      if (allocated(this%inclu)) deallocate (this%inclu)
      if (allocated(this%izpo)) deallocate (this%izpo)
      if (allocated(this%acr)) deallocate (this%acr)
      if (allocated(this%atlist)) deallocate (this%atlist)
      if (allocated(this%nn)) deallocate (this%nn)
      if (allocated(this%sbar)) deallocate (this%sbar)
      if (allocated(this%sdot)) deallocate (this%sdot)
      if (allocated(this%sbarvec)) deallocate (this%sbarvec)
      if (allocated(this%alpha)) deallocate (this%alpha)
      if (allocated(this%alpha_dot)) deallocate (this%alpha_dot)
      if (allocated(this%screening_alpha)) deallocate (this%screening_alpha)
      if (allocated(this%screening_sigma)) deallocate (this%screening_sigma)
      if (allocated(this%ijpair)) deallocate (this%ijpair)
      if (allocated(this%ijktrio)) deallocate (this%ijktrio)
      if (allocated(this%natoms_layer)) deallocate (this%natoms_layer)
#endif
   end subroutine destructor

   !> @brief Read the &lattice namelist and build derived lattice data.
   !> @details Parses bulk/surface/impurity geometry, cluster cutoffs, periodic
   !>          boundary flags, and strux-screening options, then prepares
   !>          primitive-cell and cluster state for later stages.
   !> @param[inout] this Lattice object to populate.
   !> @param[in] fname Optional input file; defaults to this%control%fname.
   !> @note This is an input boundary and may raise fatal diagnostics for invalid options.
   module subroutine build_from_file(this, fname)
      class(lattice), intent(inout) :: this
      character(len=*), intent(in), optional :: fname
      ! variables associated with the reading processes
      integer :: iostatus, funit, i
      character(len=sl) :: fname_
      real(rp), allocatable :: ct_first(:), screening_alpha_first(:), screening_sigma_first(:)

      include 'include_codes/namelists/lattice.f90'

      if (present(fname)) then
         fname_ = fname
         this%control%fname = fname
      else
         fname_ = this%control%fname
      end if

      ! Save previous values
      ! Bulk initialization
      ndim = this%ndim
      npe = this%npe
      rc = this%rc
      r2 = this%r2
      ntype = this%ntype
      crystal_sym = this%crystal_sym
      a = this%a
      wav = this%wav
      celldm = this%celldm
      strux_backend = this%strux_backend
      screening = this%screening
      strux_want_sdot = this%strux_want_sdot
      strux_solve_scale = this%strux_solve_scale
      b1 = this%b1
      b2 = this%b2
      b3 = this%b3
      n1 = this%n1 
      n2 = this%n2
      n3 = this%n3
      pbc = this%pbc
      morton_sfc = this%morton_sfc

      if (.not. allocated(this%izp) .or. size(this%izp) .ne. this%ndim) then
#ifdef USE_SAFE_ALLOC
         if (allocated(this%izp)) call g_safe_alloc%deallocate('lattice.izp', this%izp)
         call g_safe_alloc%allocate('lattice.izp', this%izp, (/this%ndim/))
#else
         if (allocated(this%izp)) deallocate (this%izp)
         allocate (this%izp(this%ndim))
#endif
      end if
      if (.not. allocated(this%no) .or. size(this%no) .ne. this%ndim) then
#ifdef USE_SAFE_ALLOC
         if (allocated(this%no)) call g_safe_alloc%deallocate('lattice.no', this%no)
         call g_safe_alloc%allocate('lattice.no', this%no, (/this%ndim/))
#else
         if (allocated(this%no)) deallocate (this%no)
         allocate (this%no(this%ndim))
#endif
      end if
      if (.not. allocated(this%crd) .or. size(this%crd) .ne. 3*this%ndim) then
#ifdef USE_SAFE_ALLOC
         if (allocated(this%crd)) call g_safe_alloc%deallocate('lattice.crd', this%crd)
         call g_safe_alloc%allocate('lattice.crd', this%crd, (/3, this%ndim/))
#else
         if (allocated(this%crd)) deallocate (this%crd)
         allocate (this%crd(3, this%ndim))
#endif
      end if

      call move_alloc(this%izp, izp)
      call move_alloc(this%no, no)
      call move_alloc(this%crd, crd)

      ! Impurity initialization
      nclu = this%nclu
      if (.not. allocated(this%inclu) .or. size(this%inclu) .ne. 3*this%nclu) then
#ifdef USE_SAFE_ALLOC
         if (allocated(this%inclu)) call g_safe_alloc%deallocate('lattice.inclu', this%inclu)
         call g_safe_alloc%allocate('lattice.inclu', this%inclu, (/this%nclu, 3/))
#else
         if (allocated(this%inclu)) deallocate (this%inclu)
         allocate (this%inclu(this%nclu, 3))
#endif
      end if
      call move_alloc(this%inclu, inclu)

      ! Surface initialization
      surftype = this%surftype
      nlay = this%nlay
      ! B7.5: frozen-boundary layer counts for the two-sided (calctype=’L’)
      ! interface geometry. Zero for every other calctype.
      nlay_a = this%nlay_a
      nlay_b = this%nlay_b
      region_b_kind = this%region_b_kind
      ntype = this%ntype
      call move_alloc(this%ct, ct)
      call move_alloc(this%screening_alpha, screening_alpha)
      call move_alloc(this%screening_sigma, screening_sigma)
      njij = this%njij
      call move_alloc(this%ijpair, ijpair)
      njijk = this%njijk
      call move_alloc(this%ijktrio, ijktrio)

      ! Pre-size ct before the first read. ntype is unknown at this point
      ! (local ntype = 0 from restore_to_default), so ct has size 0. 1000 is
      ! a safe upper bound for atom types; the resize check below will shrink
      ! it to the actual ntype read from the file.
      if (.not. allocated(ct)) then
         allocate (ct(size(izp)))
      else if (size(ct) == 0) then
         deallocate (ct)
         allocate (ct(size(izp)))
      end if
      if (.not. allocated(screening_alpha)) then
         allocate (screening_alpha(4))
         screening_alpha = this%default_screening_alpha(size(screening_alpha))
      else if (size(screening_alpha) < 4) then
         deallocate (screening_alpha)
         allocate (screening_alpha(4))
         screening_alpha = this%default_screening_alpha(size(screening_alpha))
      end if
      if (.not. allocated(screening_sigma)) then
         allocate (screening_sigma(4))
         screening_sigma = 0.7_rp
      else if (size(screening_sigma) < 4) then
         deallocate (screening_sigma)
         allocate (screening_sigma(4))
         screening_sigma = 0.7_rp
      end if
      if (.not. allocated(z)) then
         allocate (z(1))
         z = 0.0_rp
      end if
      if (.not. allocated(primcell)) then
         allocate (primcell(1, 1))
         primcell = 0.0_rp
      end if
      if (.not. allocated(crsurf)) then
         allocate (crsurf(1, 1))
         crsurf = 0.0_rp
      end if
      if (.not. allocated(cr)) then
         allocate (cr(1, 1))
         cr = 0.0_rp
      end if
      if (.not. allocated(acr)) then
         allocate (acr(1, 1))
         acr = 0.0_rp
      end if
      if (.not. allocated(reduced_acr)) then
         allocate (reduced_acr(1))
         reduced_acr = 0
      end if
      if (.not. allocated(num)) then
         allocate (num(1))
         num = 0
      end if
      if (.not. allocated(izpsurf)) then
         allocate (izpsurf(1))
         izpsurf = 0
      end if
      if (.not. allocated(izsurf)) then
         allocate (izsurf(1))
         izsurf = 0
      end if
      if (.not. allocated(nosurf)) then
         allocate (nosurf(1))
         nosurf = 0
      end if
      if (.not. allocated(izpo)) then
         allocate (izpo(1))
         izpo = 0
      end if
      if (.not. allocated(iz)) then
         allocate (iz(1))
         iz = 0
      end if
      if (.not. allocated(iu)) then
         allocate (iu(1))
         iu = 0
      end if
      if (.not. allocated(irec)) then
         allocate (irec(1))
         irec = 0
      end if
      if (.not. allocated(ib)) then
         allocate (ib(1))
         ib = 0
      end if
      if (.not. allocated(ijpair)) then
         allocate (ijpair(1, 2))
         ijpair = 0
      end if
      if (.not. allocated(ijktrio)) then
         allocate (ijktrio(1, 6))
         ijktrio = 0.0_rp
      end if

      open (newunit=funit, file=fname_, action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) then
         call g_logger%fatal('file '//fmt('A', trim(fname_))//'not found', __FILE__, __LINE__)
      end if

      read (funit, nml=lattice, iostat=iostatus)

      if (allocated(ct)) then
         allocate (ct_first(size(ct)))
         ct_first = ct
      end if
      if (allocated(screening_alpha)) then
         allocate (screening_alpha_first(size(screening_alpha)))
         screening_alpha_first = screening_alpha
      end if
      if (allocated(screening_sigma)) then
         allocate (screening_sigma_first(size(screening_sigma)))
         screening_sigma_first = screening_sigma
      end if

      if (size(izp) .ne. ndim) then
         deallocate (izp)
         allocate (izp(ndim))
      end if
      if (size(no) .ne. ndim) then
         deallocate (no)
         allocate (no(ndim))
      end if
      if (size(crd) .ne. 3*ndim) then
         deallocate (crd)
         allocate (crd(3, ndim))
      end if
      if (size(inclu) .ne. 3*nclu) then
         deallocate (inclu)
         allocate (inclu(nclu, 3))
      end if

      if (size(ijpair) .ne. 2*njij) then
         deallocate (ijpair)
         allocate (ijpair(njij, 2))
      end if

      if (size(ct) .ne. ntype) then
         deallocate (ct)
         allocate (ct(ntype))
      end if
      if (size(screening_alpha) < 4) then
         deallocate (screening_alpha)
         allocate (screening_alpha(4))
         screening_alpha = this%default_screening_alpha(size(screening_alpha))
      end if
      if (size(screening_sigma) < 4) then
         deallocate (screening_sigma)
         allocate (screening_sigma(4))
         screening_sigma = 0.7_rp
      end if

      if (size(ijktrio) .ne. 2*njijk) then
         deallocate (ijktrio)
         allocate (ijktrio(njijk, 6))
      end if

      rewind (funit)
      read (funit, nml=lattice, iostat=iostatus)
      if (iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
         call g_logger%error('Error while reading namelist', __FILE__, __LINE__)
         call g_logger%error('iostatus = '//fmt('I0', iostatus), __FILE__, __LINE__)
      end if
      close (funit)

      if (allocated(ct_first)) then
         if (allocated(ct)) then
            if (all(abs(ct) <= tiny(1.0_rp)) .and. any(abs(ct_first) > tiny(1.0_rp))) then
               ct(1:min(size(ct), size(ct_first))) = ct_first(1:min(size(ct), size(ct_first)))
            end if
         end if
         deallocate (ct_first)
      end if
      if (allocated(screening_alpha_first)) then
         if (allocated(screening_alpha)) then
            if (all(abs(screening_alpha) <= tiny(1.0_rp)) .and. any(abs(screening_alpha_first) > tiny(1.0_rp))) then
               screening_alpha(1:min(size(screening_alpha), size(screening_alpha_first))) = &
                  screening_alpha_first(1:min(size(screening_alpha), size(screening_alpha_first)))
            end if
         end if
         deallocate (screening_alpha_first)
      end if
      if (allocated(screening_sigma_first)) then
         if (allocated(screening_sigma)) then
            if (all(abs(screening_sigma) <= tiny(1.0_rp)) .and. any(abs(screening_sigma_first) > tiny(1.0_rp))) then
               screening_sigma(1:min(size(screening_sigma), size(screening_sigma_first))) = &
                  screening_sigma_first(1:min(size(screening_sigma), size(screening_sigma_first)))
            end if
         end if
         deallocate (screening_sigma_first)
      end if

      ! General intialization

      this%r2 = r2
      this%alat = alat
      this%celldm = celldm
      this%ntype = ntype
      this%strux_backend = trim(adjustl(strux_backend))
      this%screening = trim(adjustl(screening))
      this%strux_want_sdot = strux_want_sdot
      this%strux_solve_scale = strux_solve_scale
      if (len_trim(this%strux_backend) == 0) this%strux_backend = 'legacy'
      if (len_trim(this%screening) == 0) this%screening = 'default'
      if (this%strux_solve_scale <= 0.0_rp) this%strux_solve_scale = 9.0_rp
      this%pbc = pbc
      this%morton_sfc = morton_sfc
      this%b1 = b1
      this%b2 = b2
      this%b3 = b3
      this%n1 = n1
      this%n2 = n2
      this%n3 = n3

      call move_alloc(ct, this%ct)
      if (allocated(this%ct)) then
         if (size(this%ct) > 0) then
            if (all(abs(this%ct) <= tiny(1.0_rp)) .and. this%r2 > 0.0_rp) then
               this%ct = sqrt(this%r2)
            end if
         end if
      end if
      call move_alloc(screening_alpha, this%screening_alpha)
      call move_alloc(screening_sigma, this%screening_sigma)
      ! Reads Wigner-Seitz radius if available
      this%wav = wav
      ! Bulk initialization
      this%ndim = ndim
      this%npe = npe
      this%rc = rc
      this%crystal_sym = crystal_sym
      this%a = a
      call move_alloc(izp, this%izp)
      call move_alloc(no, this%no)
      call move_alloc(crd, this%crd)

      ! Impurity initialization
      this%nclu = nclu
      call move_alloc(inclu, this%inclu)

      ! Surface initialization
      this%surftype = surftype
      this%nlay = nlay
      this%nlay_a = nlay_a
      this%nlay_b = nlay_b

      ! B7.6: region B’s physical kind. A true input boundary, so the check
      ! belongs here (Phase-1 rule 3) rather than at the call sites, which then
      ! only ever see one of the two valid values. A silent fallback would be
      ! the worst option available: a misspelt ’vaccum’ would run the metallic
      ! path and report a plausible wrong barrier, which is precisely the
      ! failure mode B7 §1.3 exists to prevent.
      region_b_kind = lower(region_b_kind)
      select case (trim(region_b_kind))
      case ('metal', 'vacuum')
         this%region_b_kind = region_b_kind
      case default
         call g_logger%fatal('lattice: unknown &lattice region_b_kind='''// &
                             trim(region_b_kind)//'''. Valid values are ''metal'' '// &
                             '(region B is a second metallic reference, loaded from &atoms '// &
                             'label(:)) and ''vacuum'' (region B is semi-infinite vacuum, '// &
                             'parameters generated per run by vacuum_lead).', __FILE__, __LINE__)
      end select

      ! Exchange calculation initialization
      call move_alloc(ijpair, this%ijpair)
      this%njij = njij

      ! Spin-lattice calculation initialization
      call move_alloc(ijktrio, this%ijktrio)
      this%njijk = njijk

      ! Condition to store the trios in the original njij type

      if ((njijk .ne. 0) .and. (njij .eq. 0)) then
         this%njij = 3*this%njijk
#ifdef USE_SAFE_ALLOC
         if (allocated(this%ijpair)) call g_safe_alloc%deallocate('lattice.ijpair', this%ijpair)
         call g_safe_alloc%allocate('lattice.ijpair', this%ijpair, (/this%njij, 2/))
#else
         if (allocated(this%ijpair)) deallocate (this%ijpair)
         allocate (this%ijpair(this%njij, 2))
#endif
         do i = 1, this%njijk
            this%ijpair(3*(i - 1) + 1, 1) = int(this%ijktrio(i, 1))
            this%ijpair(3*(i - 1) + 1, 2) = int(this%ijktrio(i, 2))
            this%ijpair(3*(i - 1) + 2, 1) = int(this%ijktrio(i, 1))
            this%ijpair(3*(i - 1) + 2, 2) = int(this%ijktrio(i, 3))
            this%ijpair(3*(i - 1) + 3, 1) = int(this%ijktrio(i, 2))
            this%ijpair(3*(i - 1) + 3, 2) = int(this%ijktrio(i, 3))
         end do
      else if ((this%njijk .ne. 0) .and. (this%njij .ne. 0)) then
         call g_logger%fatal('not possible to calculate Jij and Jijk at the same time', __FILE__, __LINE__)
      end if

      ! checks
      call this%check_all()

   end subroutine build_from_file

   !> @brief Rebuild lattice state from already installed lattice members.
   !> @details Reuses current geometry and namelist-derived settings to run the
   !>          same cluster, primitive-cell, and structure setup normally reached
   !>          from build_from_file.
   !> @param[inout] this Lattice object whose current members provide the input state.
   module subroutine build_from_lattice(this)
      class(lattice), intent(inout) :: this
      integer :: nbulk_bulk, ntot, nbas, nrec, funit, iostatus, nsite_guess
      real(rp) :: r2, strux_solve_scale
      real(rp), dimension(3, 3) :: a
      real(rp), dimension(:), allocatable :: ct, screening_alpha, screening_sigma
      real(rp), dimension(:), allocatable :: ct_first, screening_alpha_first, screening_sigma_first
      integer, dimension(:), allocatable :: izp, no, iu, ib, irec
      real(rp), dimension(:, :), allocatable :: crd
      character(len=16) :: strux_backend, screening
      logical :: strux_want_sdot
      namelist /lattice/ r2, nbulk_bulk, ntot, nbas, nrec, &
         a, crd, &
         ct, izp, no, iu, ib, irec, strux_backend, screening, strux_want_sdot, &
         strux_solve_scale, screening_alpha, screening_sigma

      strux_backend = this%strux_backend
      screening = this%screening
      strux_want_sdot = this%strux_want_sdot
      strux_solve_scale = this%strux_solve_scale

      call move_alloc(this%crd, crd)
      call move_alloc(this%izp, izp)
      call move_alloc(this%no, no)
      call move_alloc(this%ct, ct)
      call move_alloc(this%screening_alpha, screening_alpha)
      call move_alloc(this%screening_sigma, screening_sigma)

      open (newunit=funit, file='lattice.nml', action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) then
         call g_logger%fatal('file lattice.nml not found', __FILE__, __LINE__)
      end if

      ! Pre-allocate ib, iu, irec before the sizing read to avoid undefined
      ! behaviour on compilers that write into unallocated storage when these
      ! variables appear in the namelist file. crd was moved from this%crd and
      ! has size (3, ntot), so its second dimension is a safe upper bound.
      nsite_guess = 1
      if (allocated(crd)) then
         nsite_guess = max(1, size(crd, 2))
         allocate (ib(size(crd, 2)), iu(size(crd, 2)), irec(size(crd, 2)))
      else if (allocated(izp)) then
         nsite_guess = max(1, size(izp))
         allocate (ib(nsite_guess), iu(nsite_guess), irec(nsite_guess))
      else
         allocate (ib(1), iu(1), irec(1))
      end if
      if (.not. allocated(ct)) then
         allocate (ct(nsite_guess))
         ct = 0.0_rp
      else if (size(ct) == 0) then
         deallocate (ct)
         allocate (ct(nsite_guess))
         ct = 0.0_rp
      end if
      if (.not. allocated(screening_alpha)) then
         allocate (screening_alpha(4))
         screening_alpha = this%default_screening_alpha(size(screening_alpha))
      else if (size(screening_alpha) < 4) then
         deallocate (screening_alpha)
         allocate (screening_alpha(4))
         screening_alpha = this%default_screening_alpha(size(screening_alpha))
      end if
      if (.not. allocated(screening_sigma)) then
         allocate (screening_sigma(4))
         screening_sigma = 0.7_rp
      else if (size(screening_sigma) < 4) then
         deallocate (screening_sigma)
         allocate (screening_sigma(4))
         screening_sigma = 0.7_rp
      end if
      read (funit, nml=lattice, iostat=iostatus)
      if (allocated(ct)) then
         allocate (ct_first(size(ct)))
         ct_first = ct
      end if
      if (allocated(screening_alpha)) then
         allocate (screening_alpha_first(size(screening_alpha)))
         screening_alpha_first = screening_alpha
      end if
      if (allocated(screening_sigma)) then
         allocate (screening_sigma_first(size(screening_sigma)))
         screening_sigma_first = screening_sigma
      end if
      deallocate (ib, iu, irec)
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('lattice.ib', ib, (/ntot/))
      call g_safe_alloc%allocate('lattice.iu', iu, (/ntot/))
      call g_safe_alloc%allocate('lattice.irec', irec, (/nrec/))
#else
      allocate (ib(ntot), iu(ntot), irec(nrec))
#endif
      if (allocated(izp)) nsite_guess = max(1, size(izp))
      if (size(ct) /= nsite_guess) then
         deallocate (ct)
         allocate (ct(nsite_guess))
         ct = 0.0_rp
      end if
      if (size(screening_alpha) < 4) then
         deallocate (screening_alpha)
         allocate (screening_alpha(4))
         screening_alpha = this%default_screening_alpha(size(screening_alpha))
      end if
      if (size(screening_sigma) < 4) then
         deallocate (screening_sigma)
         allocate (screening_sigma(4))
         screening_sigma = 0.7_rp
      end if
      rewind (funit)
      read (funit, nml=lattice, iostat=iostatus)
      if (iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
         call g_logger%error('Error while reading namelist', __FILE__, __LINE__)
         call g_logger%error('iostatus = '//fmt('I0', iostatus), __FILE__, __LINE__)
      end if
      close (funit)
      if (allocated(ct_first)) then
         if (all(abs(ct) <= tiny(1.0_rp)) .and. any(abs(ct_first) > tiny(1.0_rp))) then
            ct(1:min(size(ct), size(ct_first))) = ct_first(1:min(size(ct), size(ct_first)))
         end if
         deallocate (ct_first)
      end if
      if (allocated(screening_alpha_first)) then
         if (all(abs(screening_alpha) <= tiny(1.0_rp)) .and. any(abs(screening_alpha_first) > tiny(1.0_rp))) then
            screening_alpha(1:min(size(screening_alpha), size(screening_alpha_first))) = &
               screening_alpha_first(1:min(size(screening_alpha), size(screening_alpha_first)))
         end if
         deallocate (screening_alpha_first)
      end if
      if (allocated(screening_sigma_first)) then
         if (all(abs(screening_sigma) <= tiny(1.0_rp)) .and. any(abs(screening_sigma_first) > tiny(1.0_rp))) then
            screening_sigma(1:min(size(screening_sigma), size(screening_sigma_first))) = &
               screening_sigma_first(1:min(size(screening_sigma), size(screening_sigma_first)))
         end if
         deallocate (screening_sigma_first)
      end if

      this%nbulk_bulk = nbulk_bulk
      this%ntot = ntot
      !this%ntype = ntype
      this%nbas = nbas
      this%nrec = nrec
      !this%r2 = r2
      this%nbulk = 0
      this%strux_backend = trim(adjustl(strux_backend))
      this%screening = trim(adjustl(screening))
      this%strux_want_sdot = strux_want_sdot
      this%strux_solve_scale = strux_solve_scale
      if (len_trim(this%strux_backend) == 0) this%strux_backend = 'legacy'
      if (len_trim(this%screening) == 0) this%screening = 'default'
      if (this%strux_solve_scale <= 0.0_rp) this%strux_solve_scale = 9.0_rp

      call move_alloc(ct, this%ct)
      if (allocated(this%ct)) then
         if (size(this%ct) > 0) then
            if (all(abs(this%ct) <= tiny(1.0_rp)) .and. this%r2 > 0.0_rp) then
               this%ct = sqrt(this%r2)
            end if
         end if
      end if
      call move_alloc(screening_alpha, this%screening_alpha)
      call move_alloc(screening_sigma, this%screening_sigma)
      call move_alloc(izp, this%izp)
      call move_alloc(no, this%no)
      call move_alloc(ib, this%ib)
      call move_alloc(irec, this%irec)
      call move_alloc(iu, this%iu)
      call move_alloc(crd, this%crd)
      this%a = a
   end subroutine build_from_lattice

   !> @brief Compute primitive-cell derived quantities.
   !> @details Converts lattice vectors to internal units, computes cell volume
   !>          and Wigner-Seitz radius, and prepares inverse Cartesian lattice
   !>          data used by periodic wrapping and neighbor searches.
   !> @param[inout] this Lattice object whose a/alat/celldm state is updated.
   module subroutine build_data(this)
      class(lattice), intent(inout) :: this
      !> Local variables
      real(rp), dimension(3, 3) :: a
      integer :: i, j

      call g_timer%start('build_data')

      select case (this%crystal_sym)
      case ('bcc')
         a(:, 1) = [-0.50000000, 0.50000000, 0.50000000]
         a(:, 2) = [0.50000000, -0.50000000, 0.50000000]
         a(:, 3) = [0.50000000, 0.50000000, -0.50000000]
         this%crd(:, 1) = [0.0, 0.0, 0.0]
         this%izp(1) = 1
         this%no(1) = 1
         this%nbulk_bulk = 1
         this%ntot = 1
         if (this%control%calctype == 'B') then
            this%nrec = this%nbulk_bulk
            this%nbulk = 0
            this%ntype = 1
            this%nbas = this%ntot
#ifdef USE_SAFE_ALLOC
            call g_safe_alloc%allocate('lattice.ib', this%ib, (/this%ntot/))
            call g_safe_alloc%allocate('lattice.iu', this%iu, (/this%ntot/))
            call g_safe_alloc%allocate('lattice.irec', this%irec, (/this%nrec/))
#else
            allocate (this%ib(this%ntot), this%iu(this%ntot), this%irec(this%nrec))!, this%ct(this%ntype)) ! Now ct is defined at &lattice
#endif
            this%iu(1) = 1
            this%ib(1) = 1
            this%irec(1) = 1
         end if
         this%a = a
      case ('b2')
         a(:, 1) = [1.00000000, 0.00000000, 0.00000000]
         a(:, 2) = [0.00000000, 1.00000000, 0.00000000]
         a(:, 3) = [0.00000000, 0.00000000, 1.00000000]
         this%crd(:, 1) = [0.0, 0.0, 0.0]
         this%crd(:, 2) = [0.5, 0.5, 0.5]
         this%izp(1) = 1
         this%izp(2) = 2
         this%no(1) = 1
         this%no(2) = 2
         this%nbulk_bulk = 2
         this%ntot = 2
         if (this%control%calctype == 'B') then
            this%nrec = this%nbulk_bulk
            this%nbulk = 0
            this%ntype = 2
            this%nbas = this%ntot
#ifdef USE_SAFE_ALLOC
            call g_safe_alloc%allocate('lattice.ib', this%ib, (/this%ntot/))
            call g_safe_alloc%allocate('lattice.iu', this%iu, (/this%ntot/))
            call g_safe_alloc%allocate('lattice.irec', this%irec, (/this%nrec/))
#else
            allocate (this%ib(this%ntot), this%iu(this%ntot), this%irec(this%nrec)) !, this%ct(this%ntype)) Now ct is defined at &lattice
#endif
            this%iu(1) = 1
            this%iu(2) = 2
            this%ib(1) = 1
            this%ib(2) = 2
            this%irec(1) = 1
            this%irec(2) = 2
         end if
         this%a = a
      case ('fcc')
         a(:, 1) = [0.00000000, 0.50000000, 0.50000000]
         a(:, 2) = [0.50000000, 0.00000000, 0.50000000]
         a(:, 3) = [0.50000000, 0.50000000, 0.00000000]
         this%crd(:, 1) = [0.0, 0.0, 0.0]
         this%izp(1) = 1
         this%no(1) = 1
         this%nbulk_bulk = 1
         this%ntot = 1
         if (this%control%calctype == 'B') then
            this%nrec = this%nbulk_bulk
            this%nbulk = 0
            this%ntype = 1
            this%nbas = this%ntot
#ifdef USE_SAFE_ALLOC
            call g_safe_alloc%allocate('lattice.ib', this%ib, (/this%ntot/))
            call g_safe_alloc%allocate('lattice.iu', this%iu, (/this%ntot/))
            call g_safe_alloc%allocate('lattice.irec', this%irec, (/this%nrec/))
#else
            allocate (this%ib(this%ntot), this%iu(this%ntot), this%irec(this%nrec))!, this%ct(this%ntype)) Now ct is defined at &lattice
#endif
            this%iu(1) = 1
            this%ib(1) = 1
            this%irec(1) = 1
         end if
         this%a = a
      case ('fcc2')
         a(:, 1) = [0.00000000, 0.50000000, 0.50000000]
         a(:, 2) = [0.50000000, 0.00000000, 0.50000000]
         a(:, 3) = [0.50000000, 0.50000000, 0.00000000]
         this%crd(:, 1) = [0.00, 0.00, 0.00]
         this%crd(:, 2) = [0.50, 0.50, 0.50]
         this%izp(:) = [1, 2]
         this%no(:) = [1, 2]
         this%nbulk_bulk = 2
         this%ntot = 2
         if (this%control%calctype == 'B') then
            this%nrec = this%nbulk_bulk
            this%nbulk = 0
            this%ntype = 2
            this%nbas = this%ntot
#ifdef USE_SAFE_ALLOC
            call g_safe_alloc%allocate('lattice.ib', this%ib, (/this%ntot/))
            call g_safe_alloc%allocate('lattice.iu', this%iu, (/this%ntot/))
            call g_safe_alloc%allocate('lattice.irec', this%irec, (/this%nrec/))
#else
            allocate (this%ib(this%ntot), this%iu(this%ntot), this%irec(this%nrec))!, this%ct(this%ntype)) Now ct is defined at &lattice
#endif
            this%iu(:) = [1, 2]
            this%ib(:) = [1, 2]
            this%irec(:) = [1, 2]
         end if
         this%a = a
      case ('fcc3')
         a(:, 1) = [0.50000000, 0.50000000, 0.00000000]
         a(:, 2) = [0.50000000, 0.00000000, 0.50000000]
         a(:, 3) = [0.00000000, 0.50000000, 0.50000000]
         this%crd(:, 1) = [0.25000000, 0.25000000, 0.25000000]
         this%crd(:, 2) = [0.00000000, 0.00000000, 0.00000000]
         this%crd(:, 3) = [0.50000000, 0.50000000, 0.50000000]
         this%crd(:, 4) = [-0.25000000, -0.25000000, -0.25000000]
         this%izp(1:4) = [1, 2, 3, 4]
         this%no(1:4) = [1, 2, 3, 4]
         this%nbulk_bulk = 4
         this%ntot = 4
         if (this%control%calctype == 'B') then
            this%nrec = this%nbulk_bulk
            this%nbulk = 0
            this%ntype = 4
            this%nbas = this%ntot
#ifdef USE_SAFE_ALLOC
            call g_safe_alloc%allocate('lattice.ib', this%ib, (/this%ntot/))
            call g_safe_alloc%allocate('lattice.iu', this%iu, (/this%ntot/))
            call g_safe_alloc%allocate('lattice.irec', this%irec, (/this%nrec/))
#else
            allocate (this%ib(this%ntot), this%iu(this%ntot), this%irec(this%nrec)) !, this%ct(this%ntype)) Now ct is defined at &lattice
#endif
            this%iu(1:4) = [1, 2, 3, 4]
            this%ib(1:4) = [1, 2, 3, 4]
            this%irec(1:4) = [1, 2, 3, 4]
         end if
         this%a = a
      case ('hcp')
         if (this%celldm == 0.0d0) then
            print *, 'WARNING: hcp structure using as default c/a = 1.633'
            this%celldm = 1.633d0
         end if
         a(:, 1) = [1.00000000, 0.00000000, 0.00000000]
         a(:, 2) = [-.50000000, 0.86602500, 0.00000000]
         a(:, 3) = [0.00000000, 0.00000000, 0.00000000]
         a(3, 3) = this%celldm
         this%crd(:, 1) = [0.0, 0.0, 0.0]
         this%crd(:, 2) = [0.0, 0.577350000, 0.00000000]
         this%crd(3, 2) = (0.5d0)*this%celldm
         this%nbulk_bulk = 2
         this%ntot = 2
         this%izp(:) = [1, 2]
         this%no(:) = [1, 2]
         if (this%control%calctype == 'B') then
            this%nrec = this%nbulk_bulk
            this%nbulk = 0
            this%ntype = 2
            this%nbas = this%ntot
#ifdef USE_SAFE_ALLOC
            call g_safe_alloc%allocate('lattice.ib', this%ib, (/this%ntot/))
            call g_safe_alloc%allocate('lattice.iu', this%iu, (/this%ntot/))
            call g_safe_alloc%allocate('lattice.irec', this%irec, (/this%nrec/))
#else
            allocate (this%ib(this%ntot), this%iu(this%ntot), this%irec(this%nrec)) !, this%ct(this%ntype)) Now ct is defined at &lattice
#endif
            this%iu(:) = [1, 2]
            this%ib(:) = [1, 2]
            this%irec(:) = [1, 2]
         end if
         this%a = a
      case ('file')
         ! reads from input
         call this%build_from_lattice()
      end select

      ! Ensure neighbour-cutoff data is valid for built-in crystal presets.
      ! For legacy bulk/surface inputs, ct is often omitted and should default
      ! to include first + second shells (historical behaviour: alat + 0.1).
      if (this%crystal_sym /= 'file') then
         if (.not. allocated(this%ct) .or. size(this%ct) /= this%ntype) then
            if (allocated(this%ct)) deallocate(this%ct)
            allocate(this%ct(this%ntype))
            this%ct(:) = 0.0_rp
         end if
         if (all(abs(this%ct) <= tiny(1.0_rp))) then
            this%ct(:) = this%alat + 0.1_rp
         end if
         if (this%r2 <= 0.0_rp) then
            this%r2 = this%ct(1)**2
         end if
      end if

      ! Volume in cubic Angstroms
      this%vol = abs(dot_product(this%a(:, 3), cross_product(this%a(:, 1), this%a(:, 2))))*(this%alat**3)
      if (this%wav .eq. 0) then
         this%wav = (this%vol/((16.0d0/3.0d0)*atan(1.0d0)*this%ntot))**(1.0d0/3.0d0)
         write (*, *) 'wav', this%wav
      end if
      if (this%control%calctype == 'B' .or. this%control%calctype == 'S' &
          .or. this%control%calctype == 'L') this%nmax = 0
      call g_timer%stop('build_data')
   end subroutine build_data

   !> @brief Reset lattice members to their default values.
   !> @details Clears allocatable storage and restores scalar defaults; with the
   !>          optional full flag it also clears persistent input-facing state.
   !> @param[inout] this Lattice object to reset.
   !> @param[in] full Optional flag requesting a full reset.
   module subroutine restore_to_default(this, full)
      class(lattice) :: this
      logical, intent(in), optional :: full

      this%ndim = 9900000
      this%npe = 49
      this%nclu = 0
      this%surftype = 'none'
      this%nlay = 0
      this%nlay_a = 0
      this%nlay_b = 0
      this%region_b_kind = 'metal'
      this%ntype = 0
      this%nbas = 0
      this%wav = 0
      this%celldm = 0.0d0
      this%njij = 0
      this%njijk = 0
      this%strux_backend = 'legacy'
      this%screening = 'default'
      this%strux_want_sdot = .false.
      this%strux_solve_scale = 9.0_rp
      this%nn_max = 0
      this%a_cart_inv = 0.0_rp
      this%a_cart_inv_ready = .false.
      this%b1 = .false.
      this%b2 = .false.
      this%b3 = .false.
      this%pbc = .false.
      this%morton_sfc = .false.
      this%n1 = 0
      this%n2 = 0
      this%n3 = 0
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('lattice.izp', this%izp, (/this%ndim/))
      call g_safe_alloc%allocate('lattice.no', this%no, (/this%ndim/))
      call g_safe_alloc%allocate('lattice.crd', this%crd, (/3, this%ndim/))
      call g_safe_alloc%allocate('lattice.inclu', this%inclu, (/this%nclu, 3/))
      call g_safe_alloc%allocate('lattice.ijpair', this%ijpair, (/this%njij, 2/))
      call g_safe_alloc%allocate('lattice.ijktrio', this%ijktrio, (/this%njijk, 6/))
      call g_safe_alloc%allocate('lattice.chargetrf_type', this%chargetrf_type, (/this%nbas/))
      call g_safe_alloc%allocate('lattice.ct', this%ct, (/this%ntype/))
      call g_safe_alloc%allocate('lattice.screening_alpha', this%screening_alpha, (/max(1, this%control%lmax + 1)/))
      call g_safe_alloc%allocate('lattice.screening_sigma', this%screening_sigma, (/max(1, this%control%lmax + 1)/))
#else
      allocate (this%izp(this%ndim), this%no(this%ndim))
      allocate (this%crd(3, this%ndim))
      allocate (this%inclu(this%nclu, 3))
      allocate (this%ijpair(this%njij, 2))
      allocate (this%ijktrio(this%njijk, 6))
      allocate (this%chargetrf_type(this%nbas))
      allocate (this%ct(this%ntype))
      allocate (this%screening_alpha(max(1, this%control%lmax + 1)))
      allocate (this%screening_sigma(max(1, this%control%lmax + 1)))
#endif

      this%izp = 0.0d0
      this%no = 0
      this%crd = 0.d0
      this%inclu = 0.0d0
      this%ijpair = 0.0d0
      this%ijktrio = 0.0d0
      this%chargetrf_type = 0.0d0
      this%ct = 0.0d0
      this%screening_alpha = this%default_screening_alpha(size(this%screening_alpha))
      this%screening_sigma = 0.7_rp
      if (associated(this%control)) then
         if (present(full)) then
            if (full) then
               call this%control%restore_to_default()
            end if
         end if
      end if

   end subroutine restore_to_default

   !> @brief Build the bulk Bravais cluster.
   !> @details Expands primitive-cell sites through lattice translations, cuts
   !>          them by the configured radius, and fills bulk cluster coordinates,
   !>          atomic labels, and optional Morton ordering.
   !> @param[inout] this Lattice object receiving cr/iz/num/no cluster state.
   module subroutine bravais(this)
      class(lattice), intent(inout) :: this
      ! Local variables
      real(rp) :: rc, rs, lc, lcx, lcy, lcz
      integer, dimension(:), allocatable :: iz, num
      real(rp), dimension(:, :), allocatable :: cr, crbravais
   real(rp), allocatable :: tmp(:, :)
      integer :: npe, ndim, nx, ny, nz, npr, l, n, i, nl, k, kk
      logical :: isopen
      integer :: iostatus

      call g_timer%start('bravais')

      inquire (unit=10, opened=isopen)
      if (isopen) then
         call g_logger%fatal('lattice%bravais, file clust: Unit 10 is already open', __FILE__, __LINE__)
      else
         open (unit=10, file='clust')
      end if

      n = this%ntot
      ! Defining a size-like parameter for the clust before the cut
      npe = this%npe
      ndim = this%ndim
      rc = this%rc
      allocate (iz(ndim), num(ndim))
      allocate (cr(3, ndim), crbravais(3, ndim))
      ! Defining the cut radius
      rs = (0.8d0*int((npe/(1.0d0))/2.0d0))**2
      rs = dmin1(rs, rc)

      npr = int((ndim/(n*1.0d0))**(1.0d0/3.0d0))

      crbravais(:, :) = this%crd(:, :)
      if (this%pbc) then
         lcx = (this%n1 + 1)/2
         lcy = (this%n2 + 1)/2
         lcz = (this%n3 + 1)/2
      else
         lc = (npr + 1)/2
         lcx = lc
         lcy = lc
         lcz = lc
      end if
      l = n
      if (.not.this%pbc) then
         this%n1 = npr; this%n2 = npr; this%n3 = npr
      end if

      do i = 1, l
         do nx = 1, this%n1        
            do ny = 1, this%n2 
               do nz = 1, this%n3  
                  if (nx .eq. lcx .and. ny .eq. lcy .and. nz .eq. lcz) go to 13
                  n = n + 1
!  ...........Verifies dimension NDIM.........................
                  if (n .gt. ndim) then
                     write (6, 600) n
                     stop
                  end if
! ...............................................................
                  crbravais(1, n) = this%crd(1, i) + (nx - lcx)*this%a(1, 1) + (ny - lcy)*this%a(1, 2) + (nz - lcz)*this%a(1, 3)
                  crbravais(2, n) = this%crd(2, i) + (nx - lcx)*this%a(2, 1) + (ny - lcy)*this%a(2, 2) + (nz - lcz)*this%a(2, 3)
                  crbravais(3, n) = this%crd(3, i) + (nx - lcx)*this%a(3, 1) + (ny - lcy)*this%a(3, 2) + (nz - lcz)*this%a(3, 3)
                  this%no(n) = this%no(i)
                  this%izp(n) = this%izp(i)
13                continue
               end do
            end do
         end do
      end do
      nl = l*(npr**3)

      kk = 0
      if (rc == 0.0d0) rs = npr**3

      if (this%pbc) then
         ndim = this%n1*this%n2*this%n3*l 
         iz = this%izp
         num = this%no
         cr = crbravais
         kk = n
      else
          do i = 1, nl
             call cut(i, l, ndim, crbravais, cr, this%izp, iz, num, this%no, rs, kk)
          end do
      end if 

      if (int(kk/2) /= kk/2.d0) kk = kk - 1

      write (10, 300) kk
      do k = 1, kk, 2
         write (10, 200) (cr(1, k + i - 1), cr(2, k + i - 1), cr(3, k + i - 1), iz(k + i - 1), num(k + i - 1), i=1, 2)
      end do
      write (800, *) kk
      do k = 1, kk
         write (800, '(i5,3f14.8)') iz(k), cr(:, k)
      end do

200   format(3F14.8, 2I4, 3F14.8, 2I4)
300   format(3X, "II =", I7)
600   format(1X, 'K =', I7, 5X, 'redefine the dimension ndim')

      this%kk = kk
      call move_alloc(cr, this%cr)
      ! Shrink this%cr to the actual cluster size kk to avoid keeping the
      ! large temporary allocation of shape (3, ndim). We copy the first kk
      ! columns into a smaller array and move that allocation into this%cr.
      if (allocated(this%cr)) then
         if (size(this%cr, 2) > kk) then
            allocate(tmp(3, kk))
            tmp(:, :) = this%cr(:, 1:kk)
            call move_alloc(tmp, this%cr)
         end if
      end if
      call move_alloc(iz, this%iz)
      call move_alloc(num, this%num)
      close (10)

      ! Optional Morton (Z-order) reordering of the bulk cluster to improve
      ! data locality of the recursion/Chebyshev SpMV kernels. Atom 1 (the
      ! central reference site) is kept first; only atoms 2..kk are permuted.
      if (this%morton_sfc) then
         call this%morton_reorder_bulk()
      end if

      ! Test to set nmax to the whole cluster
      ! this%nmax = kk
      call g_timer%stop('bravais')
   end subroutine bravais

   !> @brief Reorder the bulk cluster by Morton space-filling-curve key.
   !> @details Sorts cluster atoms in normalized Cartesian space to improve
   !>          locality for real-space recursion kernels while preserving the
   !>          per-atom coordinate/type arrays as a consistent permutation.
   !> @param[inout] this Lattice object whose bulk cluster arrays are reordered.
   module subroutine morton_reorder_bulk(this)
      class(lattice), intent(inout) :: this
      integer :: kk, i, n
      integer, allocatable :: perm(:)
      integer(8), allocatable :: key(:)
      real(rp) :: cmin(3), cmax(3), span(3), inv_span(3)
      integer, parameter :: nbits = 21      ! 21 bits/axis -> 63-bit Morton key
      integer, parameter :: gmax = 2**nbits - 1
      real(rp), allocatable :: cr_new(:, :)
      integer, allocatable :: iz_new(:), num_new(:)
      integer :: q(3)

      kk = this%kk
      if (kk <= 2) return

      ! Bounding box of the full cluster (used for quantization).
      do i = 1, 3
         cmin(i) = minval(this%cr(i, 1:kk))
         cmax(i) = maxval(this%cr(i, 1:kk))
      end do
      span(:) = cmax(:) - cmin(:)
      do i = 1, 3
         if (span(i) > tiny(1.0_rp)) then
            inv_span(i) = real(gmax, rp)/span(i)
         else
            inv_span(i) = 0.0_rp
         end if
      end do

      ! Compute a Morton key for atoms 2..kk; atom 1 stays first via a key of 0.
      allocate (key(kk), perm(kk))
      key(1) = -1_8                          ! sorts before any non-negative key
      perm(1) = 1
      do n = 2, kk
         do i = 1, 3
            q(i) = nint((this%cr(i, n) - cmin(i))*inv_span(i))
            if (q(i) < 0) q(i) = 0
            if (q(i) > gmax) q(i) = gmax
         end do
         key(n) = morton_encode3(q(1), q(2), q(3), nbits)
         perm(n) = n
      end do

      ! Sort indices 1..kk by Morton key (stable not required; ties are
      ! spatially coincident bins).
      call sort_by_key(key, perm, kk)

      ! Apply the permutation to the cluster arrays.
      allocate (cr_new(3, kk), iz_new(kk), num_new(kk))
      do n = 1, kk
         cr_new(:, n) = this%cr(:, perm(n))
         iz_new(n) = this%iz(perm(n))
         num_new(n) = this%num(perm(n))
      end do
      this%cr(:, 1:kk) = cr_new(:, 1:kk)
      this%iz(1:kk) = iz_new(1:kk)
      this%num(1:kk) = num_new(1:kk)

      deallocate (cr_new, iz_new, num_new, key, perm)

      call g_logger%info('lattice%bravais: bulk cluster reordered along Morton (Z-order) curve', __FILE__, __LINE__)
   end subroutine morton_reorder_bulk

   !> @brief Encode a 3D integer grid coordinate as a Morton key.
   !> @param[in] x Grid coordinate on the first axis.
   !> @param[in] y Grid coordinate on the second axis.
   !> @param[in] z Grid coordinate on the third axis.
   !> @param[in] nbits Number of bits to interleave from each coordinate.
   !> @return Interleaved 64-bit Morton code.
   module pure function morton_encode3(x, y, z, nbits) result(code)
      integer, intent(in) :: x, y, z, nbits
      integer(8) :: code
      integer :: b
      code = 0_8
      do b = 0, nbits - 1
         code = ior(code, ishft(iand(int(x, 8), ishft(1_8, b)), 2*b))
         code = ior(code, ishft(iand(int(y, 8), ishft(1_8, b)), 2*b + 1))
         code = ior(code, ishft(iand(int(z, 8), ishft(1_8, b)), 2*b + 2))
      end do
   end function morton_encode3

   !> @brief Sort an integer permutation by associated Morton keys.
   !> @details Uses an in-place heapsort so key and permutation arrays are kept
   !>          synchronized without extra sorting dependencies.
   !> @param[inout] key Keys to sort in ascending order.
   !> @param[inout] perm Permutation entries carried with each key.
   !> @param[in] n Number of active entries in key and perm.
   module subroutine sort_by_key(key, perm, n)
      integer(8), intent(inout) :: key(:)
      integer, intent(inout) :: perm(:)
      integer, intent(in) :: n
      integer :: i, start, bottom
      integer(8) :: tkey
      integer :: tperm

      ! Build heap then sift down (standard heapsort, O(n log n)).
      do start = n/2, 1, -1
         call sift_down(key, perm, start, n)
      end do
      do bottom = n, 2, -1
         tkey = key(1); key(1) = key(bottom); key(bottom) = tkey
         tperm = perm(1); perm(1) = perm(bottom); perm(bottom) = tperm
         call sift_down(key, perm, 1, bottom - 1)
      end do
   contains
      subroutine sift_down(key, perm, start, end)
         integer(8), intent(inout) :: key(:)
         integer, intent(inout) :: perm(:)
         integer, intent(in) :: start, end
         integer :: root, child
         integer(8) :: sk
         integer :: sp
         root = start
         do while (2*root <= end)
            child = 2*root
            if (child < end) then
               if (key(child) < key(child + 1)) child = child + 1
            end if
            if (key(root) < key(child)) then
               sk = key(root); key(root) = key(child); key(child) = sk
               sp = perm(root); perm(root) = perm(child); perm(child) = sp
               root = child
            else
               return
            end if
         end do
      end subroutine sift_down
   end subroutine sort_by_key

   !> @brief Count reduced basis/species entries from the cluster atom labels.
   !> @details Computes nbas/reduced_nbas-style counts used by charge and
   !>          Hamiltonian setup after the cluster atom types are known.
   !> @param[inout] this Lattice object whose basis counters are updated.
   module subroutine calculate_nbas(this)
      implicit none
      class(lattice) :: this
      integer :: size_iz
      integer, dimension(:) :: atype(this%nbas), amount(this%nbas)
      integer :: i, j

      amount(:) = 0
      atype(:) = 0

      j = 1

      amount(j) = 1
      atype(j) = this%iz(j)

      do i = 2, this%nbas
         if (this%iz(i) .eq. this%iz(i - 1)) then
            amount(j) = amount(j) + 1
         else
            j = j + 1
            atype(j) = this%iz(i)
            amount(j) = 1
         end if
      end do

      this%reduced_nbas = 0

      do i = 1, size(amount)
         if (amount(i) .eq. 0) exit
         this%reduced_nbas = this%reduced_nbas + 1
      end do

#ifdef USE_SAFE_ALLOC
      if (allocated(this%reduced_acr)) call g_safe_alloc%deallocate('lattice.reduced_acr', this%reduced_acr)
      call g_safe_alloc%allocate('lattice.reduced_acr', this%reduced_acr, (/this%reduced_nbas/))
#else
      if (allocated(this%reduced_acr)) deallocate (this%reduced_acr)
      allocate (this%reduced_acr(this%reduced_nbas))
#endif

      do i = 1, size(this%reduced_acr)
         this%reduced_acr(i) = atype(i)
      end do

   end subroutine calculate_nbas

   !> @brief Check lattice option consistency after input parsing.
   !> @param[inout] this Lattice object whose options are validated.
   module subroutine check_all(this)
      implicit none
      class(lattice) :: this

      if (this%crystal_sym /= 'bcc' &
          .and. this%crystal_sym /= 'b2' &
          .and. this%crystal_sym /= 'fcc' &
          .and. this%crystal_sym /= 'hcp' &
          .and. this%crystal_sym /= 'fcc2' &
          .and. this%crystal_sym /= 'fcc3' &
          .and. this%crystal_sym /= 'file') then
         call g_logger%fatal('lattice%crystal_sym must be one of: ''bcc'', ''fcc'', ''hcp'' or  ''file''")', __FILE__, __LINE__)
      end if
   end subroutine check_all

end submodule lattice_lifecycle
