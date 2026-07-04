submodule (lattice_mod) lattice_print
   implicit none

contains

   module subroutine print_state_full(this, unit, file)
      class(lattice), intent(in) :: this

      integer, intent(in), optional :: unit
      character(len=*), intent(in), optional :: file
      integer :: newunit

      include 'include_codes/namelists/lattice.f90'

      ! scalar

      zmin = this%zmin
      zmax = this%zmax
      zstep = this%zstep
      wav = this%wav
      vol = this%vol
      rc = this%rc
      r2 = this%r2
      celldm = this%celldm
      alat = this%alat
      reduced_nbas = this%reduced_nbas
      ntype = this%ntype
      ntot = this%ntot
      nrec = this%nrec
      nmax = this%nmax
      nlay = this%nlay
      ndim = this%ndim
      npe = this%npe
      nclu = this%nclu
      nbulk_bulk = this%nbulk_bulk
      nbulk = this%nbulk
      nbas = this%nbas
      kk = this%kk
      dx = this%dx
      dy = this%dy
      dz = this%dz
      dw = this%dw
      crystal_sym = this%crystal_sym
      surftype = this%surftype
      strux_backend = this%strux_backend
      screening = this%screening
      strux_want_sdot = this%strux_want_sdot
      strux_solve_scale = this%strux_solve_scale
      a = this%a

      ! one dimensional allocatables

      if (allocated(this%z)) then
         allocate (z, mold=this%z)
         z = this%z
      else
         allocate (z(0))
      end if
      if (allocated(this%ct)) then
         allocate (ct, mold=this%ct)
         ct = this%ct
      else
         allocate (ct(0))
      end if
      if (allocated(this%screening_alpha)) then
         allocate (screening_alpha, mold=this%screening_alpha)
         screening_alpha = this%screening_alpha
      else
         allocate (screening_alpha(0))
      end if
      if (allocated(this%screening_sigma)) then
         allocate (screening_sigma, mold=this%screening_sigma)
         screening_sigma = this%screening_sigma
      else
         allocate (screening_sigma(0))
      end if
      if (allocated(this%reduced_acr)) then
         allocate (reduced_acr, mold=this%reduced_acr)
         reduced_acr = this%reduced_acr
      else
         allocate (reduced_acr(0))
      end if
      if (allocated(this%num)) then
         allocate (num, mold=this%num)
         num = this%num
      else
         allocate (num(0))
      end if
      if (allocated(this%no)) then
         allocate (no, mold=this%no)
         no = this%no
      else
         allocate (no(0))
      end if
      if (allocated(this%izpsurf)) then
         allocate (izpsurf, mold=this%izpsurf)
         izpsurf = this%izpsurf
      else
         allocate (izpsurf(0))
      end if
      if (allocated(this%izsurf)) then
         allocate (izsurf, mold=this%izsurf)
         izsurf = this%izsurf
      else
         allocate (izsurf(0))
      end if
      if (allocated(this%nosurf)) then
         allocate (nosurf, mold=this%nosurf)
         nosurf = this%nosurf
      else
         allocate (nosurf(0))
      end if
      if (allocated(this%izpo)) then
         allocate (izpo, mold=this%izpo)
         izpo = this%izpo
      else
         allocate (izpo(0))
      end if
      if (allocated(this%izp)) then
         allocate (izp, mold=this%izp)
         izp = this%izp
      else
         allocate (izp(0))
      end if
      if (allocated(this%iz)) then
         allocate (iz, mold=this%iz)
         iz = this%iz
      else
         allocate (iz(0))
      end if
      if (allocated(this%iu)) then
         allocate (iu, mold=this%iu)
         iu = this%iu
      else
         allocate (iu(0))
      end if
      if (allocated(this%irec)) then
         allocate (irec, mold=this%irec)
         irec = this%irec
      else
         allocate (irec(0))
      end if
      if (allocated(this%ib)) then
         allocate (ib, mold=this%ib)
         ib = this%ib
      else
         allocate (ib(0))
      end if

      ! two dimensional allocatables

      if (allocated(this%primcell)) then
         allocate (primcell, mold=this%primcell)
         primcell = this%primcell
      else
         allocate (primcell(0, 0))
      end if
      if (allocated(this%inclu)) then
         allocate (inclu, mold=this%inclu)
         inclu = this%inclu
      else
         allocate (inclu(0, 0))
      end if
      if (allocated(this%crsurf)) then
         allocate (crsurf, mold=this%crsurf)
         crsurf = this%crsurf
      else
         allocate (crsurf(0, 0))
      end if
      if (allocated(this%crd)) then
         allocate (crd, mold=this%crd)
         crd = this%crd
      else
         allocate (crd(0, 0))
      end if
      if (allocated(this%cr)) then
         allocate (cr, mold=this%cr)
         cr = this%cr
      else
         allocate (cr(0, 0))
      end if
      if (allocated(this%acr)) then
         allocate (acr, mold=this%acr)
         acr = this%acr
      else
         allocate (acr(0, 0))
      end if

      if (present(unit) .and. present(file)) then
         call g_logger%fatal('Argument error: both unit and file are present', __FILE__, __LINE__)
      else if (present(unit)) then
         write (unit, nml=lattice)
      else if (present(file)) then
         open (newunit=newunit, file=file)
         write (newunit, nml=lattice)
         close (newunit)
      else
         write (*, nml=lattice)
      end if

   end subroutine print_state_full

   module subroutine print_state(this, unit, file)
      implicit none
      class(lattice), intent(in) :: this

      integer, intent(in), optional :: unit
      character(len=*), intent(in), optional :: file
      integer :: newunit

      include 'include_codes/namelists/lattice.f90'

      ! scalar

      wav = this%wav
      rc = this%rc
      celldm = this%celldm
      alat = this%alat
      nlay = this%nlay
      ndim = this%ndim
      npe = this%npe
      nclu = this%nclu
      crystal_sym = this%crystal_sym
      surftype = this%surftype
      a = this%a

      ! ! one dimensional allocatables

      if (allocated(this%no)) then
         allocate (no, mold=this%no)
         no = this%no
      else
         allocate (no(0))
      end if
      if (allocated(this%izp)) then
         allocate (izp, mold=this%izp)
         izp = this%izp
      else
         allocate (izp(0))
      end if

      ! ! two dimensional allocatables

      if (allocated(this%inclu)) then
         allocate (inclu, mold=this%inclu)
         inclu = this%inclu
      else
         allocate (inclu(0, 0))
      end if
      if (allocated(this%crd)) then
         allocate (crd, mold=this%crd)
         crd = this%crd
      else
         allocate (crd(0, 0))
      end if

      if (present(unit) .and. present(file)) then
         call g_logger%fatal('Argument error: both unit and file are present', __FILE__, __LINE__)
      else if (present(unit)) then
         call this%print_state_formatted(unit=unit)
      else if (present(file)) then
         call this%print_state_formatted(file=file)
      else
         call this%print_state_formatted()
      end if

   end subroutine print_state

   module subroutine print_state_formatted(this, unit, file)
      class(lattice), intent(in) :: this

      integer, intent(in), optional :: unit
      character(len=*), intent(in), optional :: file

      type(namelist_generator) :: nml
      integer :: i

      nml = namelist_generator('lattice')

      ! scalar

      call nml%add('zmin', this%zmin)
      call nml%add('zmax', this%zmax)
      call nml%add('zstep', this%zstep)
      call nml%add('wav', this%wav)
      call nml%add('vol', this%vol)
      call nml%add('rc', this%rc)
      call nml%add('r2', this%r2)
      call nml%add('celldm', this%celldm)
      call nml%add('alat', this%alat)
      call nml%add('reduced_nbas', this%reduced_nbas)
      call nml%add('ntype', this%ntype)
      call nml%add('ntot', this%ntot)
      call nml%add('nrec', this%nrec)
      call nml%add('nmax', this%nmax)
      call nml%add('nlay', this%nlay)
      call nml%add('ndim', this%ndim)
      call nml%add('npe', this%npe)
      call nml%add('nclu', this%nclu)
      call nml%add('nbulk_bulk', this%nbulk_bulk)
      call nml%add('nbulk', this%nbulk)
      call nml%add('nbas', this%nbas)
      call nml%add('kk', this%kk)
      call nml%add('dx', this%dx)
      call nml%add('dy', this%dy)
      call nml%add('dz', this%dz)
      call nml%add('dw', this%dw)
      call nml%add('crystal_sym', this%crystal_sym)
      call nml%add('surftype', this%surftype)
      call nml%add('njij', this%njij)
      call nml%add('njijk', this%njijk)
      call nml%add('strux_backend', trim(this%strux_backend))
      call nml%add('screening', trim(this%screening))
      call nml%add('strux_want_sdot', this%strux_want_sdot)
      call nml%add('strux_solve_scale', this%strux_solve_scale)
      call nml%add('morton_sfc', this%morton_sfc)

      ! one dimensional allocatables
      ! TODO: implement test inside namelist_generator
      if (allocated(this%z)) call nml%add('z', this%z)
      if (allocated(this%ct)) call nml%add('ct', this%ct)
      if (allocated(this%screening_alpha)) call nml%add('screening_alpha', this%screening_alpha)
      if (allocated(this%screening_sigma)) call nml%add('screening_sigma', this%screening_sigma)
      if (allocated(this%reduced_acr)) call nml%add('reduced_acr', this%reduced_acr)
      if (allocated(this%num)) call nml%add('num', this%num)
      if (allocated(this%no)) call nml%add('no', this%no)
      if (allocated(this%izpsurf)) call nml%add('izpsurf', this%izpsurf)
      if (allocated(this%izsurf)) call nml%add('izsurf', this%izsurf)
      if (allocated(this%nosurf)) call nml%add('nosurf', this%nosurf)
      if (allocated(this%izpo)) call nml%add('izpo', this%izpo)
      if (allocated(this%izp)) call nml%add('izp', this%izp)
      if (allocated(this%iz)) call nml%add('iz', this%iz)
      if (allocated(this%iu)) call nml%add('iu', this%iu)
      if (allocated(this%irec)) call nml%add('irec', this%irec)
      if (allocated(this%ib)) call nml%add('ib', this%ib)

      ! ! two dimensional allocatables
      ! TODO: implement test inside namelist_generator
      call nml%add('a', this%a)
      if (allocated(this%primcell)) call nml%add('primcell', this%primcell)
      if (allocated(this%inclu)) call nml%add('inclu', this%inclu)
      if (allocated(this%crsurf)) call nml%add('crsurf', this%crsurf)
      if (allocated(this%crd)) call nml%add('crd', this%crd)
      if (allocated(this%cr)) call nml%add('cr', this%cr)
      if (allocated(this%acr)) call nml%add('acr', this%acr)

      if (present(unit) .and. present(file)) then
         call g_logger%fatal('Argument error: both unit and file are present', __FILE__, __LINE__)
      else if (present(unit)) then
         call nml%generate_namelist(unit=unit)
      else if (present(file)) then
         call nml%generate_namelist(file=file)
      else
         call nml%generate_namelist()
      end if
   end subroutine print_state_formatted

end submodule lattice_print
