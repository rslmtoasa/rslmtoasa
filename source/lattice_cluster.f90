submodule (lattice_mod) lattice_cluster
   implicit none

contains

   !> @brief Prepare surface-layer cluster metadata from a bulk cluster.
   !> @details Builds layer z positions and surface index arrays used by the
   !>          buildsurf path before selecting the active surface atoms.
   !> @param[inout] this Lattice object whose surface helper arrays are filled.
   module subroutine build_clusup(this)
      class(lattice), intent(inout) :: this
      character(len=10) :: surftype
      ! Local variables
      real(rp), dimension(:), allocatable :: z
      integer, dimension(:), allocatable :: izpsurf
      real(rp) :: ds, ds2, new, one
      integer :: i, j, n

      select case (this%crystal_sym)
      case ('hcp')
         read (this%surftype, *) this%dx, this%dy, this%dz, this%dw
         if ((this%dz) .ne. (-1*(this%dx + this%dy))) then
            this%dz = -1*(this%dx + this%dy)
            call g_logger%fatal('Hexagonal Miller indices not right. Does it should be ['// &
                                fmt('I0', this%dx)//fmt('I0', this%dy)//fmt('I0', this%dz)//fmt('I0', this%dw)//']?', &
                                __FILE__, __LINE__)
         end if
         this%dx = (2*this%dx + this%dy)
         this%dy = (this%dx + 2*this%dy)
         this%dz = this%dw
      case default
         read (this%surftype, *) this%dx, this%dy, this%dz
      end select

!  ............Find zstep, zmin, zmax....................
      ds = 1000.0d0; ds2 = 1000.0d0
      select case (this%crystal_sym)
      case ('hcp')
         do i = 1, this%kk
            do j = 1, this%kk
               new = this%dx*this%cr(1, i) + this%dy*this%cr(2, i) + this%dz*this%cr(3, i)
               one = this%dx*this%cr(2, j) + this%dy*this%cr(2, j) + this%dz*this%cr(3, j)
               if ((abs(new - one) .gt. 1.d-6) .and. (abs(new - one) .le. ds)) then
                  ds = abs(new - one)
               end if
               if ((abs(new) .le. ds2)) then
                  ds2 = abs(new)
               end if
            end do
         end do
      case default
         do i = 1, this%kk
            new = this%dx*this%cr(1, i) + this%dy*this%cr(2, i) + this%dz*this%cr(3, i)
            if ((abs(new) .gt. 1.d-6) .and. (abs(new) .le. ds)) then
               ds = abs(new)
            end if
            if ((abs(new) .le. ds2)) then
               ds2 = abs(new)
            end if
         end do
      end select

      this%zstep = ds
      this%zmin = ds2 - this%zstep
      this%zmax = ds2 + 15*this%zstep

!  ......................................................

      n = int((this%zmax - this%zmin)/this%zstep) + 1

      allocate (z(n), izpsurf(n))
      do i = 1, n
         z(i) = this%zmin + ((i - 1)*this%zstep)
      end do
      call move_alloc(z, this%z)

      do i = 1, n
         izpsurf(i) = mod(i, this%ntot) + 1
      end do

      j = this%ntot
      do i = 1, this%nlay
         j = j + 1
         izpsurf(i) = j
      end do
      call move_alloc(izpsurf, this%izpsurf)
      if (this%control%calctype == 'S') then
         this%ntype = this%nbulk_bulk + this%nlay
         this%nbulk = this%nbulk_bulk
         this%nrec = this%ntype - this%nbulk
         this%nbas = 51
#ifdef USE_SAFE_ALLOC
         call g_safe_alloc%allocate('lattice.ib', this%ib, (/this%nbulk/))
         call g_safe_alloc%allocate('lattice.irec', this%irec, (/this%nrec/))
         call g_safe_alloc%allocate('lattice.iu', this%iu, (/this%ntot/))
         call g_safe_alloc%allocate('lattice.ct', this%ct, (/this%ntype/))
#else
         allocate (this%ib(this%nbulk), this%irec(this%nrec), this%iu(this%ntot)) !, this%ct(this%ntype)) Now ct is defined at &lattice
#endif
!         this%ct(:) = this%alat + 0.1d0
!         this%r2 = this%ct(1)**2
      end if

      this%surftype = clean_str(this%surftype)
   end subroutine build_clusup

   !> @brief Build the full surface cluster representation.
   !> @details Selects atoms by surface orientation and layer bounds, identifies
   !>          unique layer/type representatives, and installs the full surface
   !>          coordinate/type arrays for downstream surface calculations.
   !> @param[inout] this Lattice object receiving the full surface cluster.
   module subroutine build_surf_full(this)
      class(lattice), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, natoms, nsurf, currentType, newType, newCrystalType
      integer :: nTypesTotal, nUnique, atomIdx, nTypesInLayer, ichoice, ichoicetype
      real(rp), dimension(:, :), allocatable :: crsurf
      real(rp), dimension(:), allocatable :: crh, crhd, z
      integer, dimension(:), allocatable :: atomType, crystalType, typesurf, crystalsurf, uniqueTypes, ichoicen, ichoicetypen
      integer, dimension(:), allocatable :: nTypesForCurrentLayer
      real(rp) :: dx, dy, dz, new, one, ds, ds2, disi, disi_min
      real(rp) :: zstep, zmin, zmax
      integer :: n, atomCount, maxType, nlay
      logical :: isUnique, isopen
      character(20) :: header
      ! Variables
      real(rp) :: rotated_cr(3, this%kk)
      real(rp) :: rotation_matrix(3, 3)
      real(rp) :: axis(3), theta
      real(rp) :: norm

      ! Initial definitions
      natoms = this%kk

      ! Open the file and read header
      inquire (unit=10, opened=isopen)
      if (isopen) then
         call g_logger%fatal('lattice%build_surf, file clust: Unit 10 is already open', __FILE__, __LINE__)
      else
         open (unit=10, file='clust')
      end if
      ! Open the output for surface information
      inquire (unit=20, opened=isopen)
      if (isopen) then
         call g_logger%fatal('lattice%build_surf, file clust: Unit 10 is already open', __FILE__, __LINE__)
      else
         open (unit=20, file='surfclu.out')
      end if

      select case (this%crystal_sym)
      case ('hcp')
         read (this%surftype, *) this%dx, this%dy, this%dz, this%dw
         if ((this%dz) .ne. (-1*(this%dx + this%dy))) then
            this%dz = -1*(this%dx + this%dy)
            call g_logger%fatal('Hexagonal Miller indices not right. Does it should be ['// &
                                fmt('I0', this%dx)//fmt('I0', this%dy)//fmt('I0', this%dz)//fmt('I0', this%dw)//']?', &
                                __FILE__, __LINE__)
         end if
         this%dx = (2*this%dx + this%dy)
         this%dy = (this%dx + 2*this%dy)
         this%dz = this%dw
      case default
         read (this%surftype, *) this%dx, this%dy, this%dz
      end select

      ! Allocate arrays to store atom details
      allocate (crsurf(3, this%kk), crh(this%kk), typesurf(this%kk), crystalsurf(this%kk))
      allocate (uniqueTypes(this%kk), ichoicen(this%kk), ichoicetypen(this%kk), crhd(this%kk))

      ds = 1000.0d0
      ds2 = 1000.0d0
      do i = 1, this%kk
         do j = 1, this%kk
            new = this%dx*this%cr(1, i) + this%dy*this%cr(2, i) + this%dz*this%cr(3, i)
            one = this%dx*this%cr(1, j) + this%dy*this%cr(2, j) + this%dz*this%cr(3, j)
            if ((abs(new - one) .gt. 1.d-6) .and. (abs(new - one) .le. ds)) then
               ds = abs(new - one)
            end if
            if ((abs(new) .le. ds2)) then
               ds2 = abs(new)
            end if
         end do
      end do

      this%zstep = ds
      this%zmin = ds2 - this%zstep
      this%zmax = ds2 + 20*this%zstep

      n = int((this%zmax - this%zmin)/this%zstep) + 1

      ! Allocate the number of layers
      allocate (z(n))
      ! Allocate the number of atoms per layer variable
      allocate (nTypesForCurrentLayer(n))

      do i = 1, n
         z(i) = this%zmin + ((i - 1)*this%zstep)
      end do
      call move_alloc(z, this%z)

      ! Determining max atom type from input
      maxType = maxval(this%iz(:))
      call move_alloc(this%iz, atomType)
      call move_alloc(this%num, crystalType)

      nsurf = 0
      do i = 1, n
         nTypesForCurrentLayer(i) = 0
         disi_min = sqrt(this%z(i)**2) + 1.0d0
         do k = 1, natoms
            crh(k) = 0.0d0
            crh(k) = dot_product([this%dx, this%dy, this%dz], this%cr(:, k))
            if (abs(crh(k) - this%z(i)) < 1.0d-6) then
               nsurf = nsurf + 1
               crsurf(:, nsurf) = this%cr(:, k)
               crhd(nsurf) = dot_product([this%dx, this%dy, this%dz], crsurf(:, nsurf))
               if (i <= this%nlay) then
                  ! Determine the unique atom types and update typesurf and crystalsurf.
                  isUnique = .true.
                  do j = 1, nTypesForCurrentLayer(i)
                     if (atomType(k) == uniqueTypes(j)) then
                        isUnique = .false.
                        exit
                     end if
                  end do

                  if (isUnique) then
                     nTypesForCurrentLayer(i) = nTypesForCurrentLayer(i) + 1
                     uniqueTypes(nTypesForCurrentLayer(i)) = atomType(k)
                     maxType = maxType + 1
                     typesurf(nsurf) = maxType
                     crystalsurf(nsurf) = crystalType(k) !maxType
                     disi = norm2(crsurf(:, nsurf))
                     if (disi < disi_min) then
                        ichoicen(typesurf(nsurf)) = nsurf
                        ichoicetypen(typesurf(nsurf)) = typesurf(nsurf)
                     end if
                  else
                     typesurf(nsurf) = maxType - nTypesForCurrentLayer(i) + findloc(uniqueTypes, atomType(k), dim=1)
                     crystalsurf(nsurf) = crystalType(k)
                     disi = norm2(crsurf(:, nsurf))
                     if (disi < disi_min) then
                        ichoicen(typesurf(nsurf)) = nsurf
                        ichoicetypen(typesurf(nsurf)) = typesurf(nsurf)
                     end if
                  end if
               else
                  typesurf(nsurf) = atomType(k)
                  crystalsurf(nsurf) = crystalType(k)
                  disi = norm2(crsurf(:, nsurf))
                  if (i .le. (this%nlay + this%nbulk_bulk)) then
                     if (disi < disi_min) then
                        ichoicen(typesurf(nsurf)) = nsurf
                        ichoicetypen(typesurf(nsurf)) = typesurf(nsurf)
                     end if
                  end if
               end if
            end if
         end do
      end do

      ! Passing onto the global variable
      call move_alloc(nTypesForCurrentLayer, this%natoms_layer)
      write (20, *) 'Maxtype:', maxType
      write (20, *) 'Type of atoms chosen:', ichoicetypen(1:maxType)
      write (20, *) 'Atoms chosen:', ichoicen(1:maxType)
      write (20, *) 'Unique atoms type per layer:', this%natoms_layer(:)
      write (20, *) 'Layers are:', this%z
      write (20, *) 'Step is:', this%zstep, this%zmax, this%zmin
      if (this%control%calctype == 'S') then
         this%ntype = maxType
         this%nbulk = this%nbulk_bulk
         this%nrec = this%ntype - this%nbulk
         this%nbas = 49

         if (allocated(this%chargetrf_type)) deallocate (this%chargetrf_type)
         if (allocated(this%ib)) deallocate (this%ib)
         if (allocated(this%irec)) deallocate (this%irec)
         if (allocated(this%iu)) deallocate (this%iu)
         allocate (this%ib(this%nbulk), this%irec(this%nrec), this%iu(this%ntot))!, this%ct(this%ntype)) Now ct is defined at &lattice
         allocate (this%chargetrf_type(this%nbas))

         do i = 1, this%nrec
            this%irec(i) = ichoicen(this%nbulk + i)
         end do
         do i = 1, this%ntot
            this%iu(ichoicetypen(i)) = ichoicen(i)
         end do
         do i = 1, this%nbulk
            this%ib(ichoicetypen(i)) = ichoicen(i)
         end do
      end if

      if (int(nsurf/2) /= nsurf/2.d0) nsurf = nsurf - 1

      ! Rotate the surface cluster so that z is perpendicular to the plane
      ! Calculate the unit vector normal to the crystallographic plane
    !! Setup the rotation matrix using Rodrigues´ formula
    !! Rotate the primitive lattice vectors

      write (10, 10004) nsurf
      do k = 1, nsurf - 1, 2
         write (10, 10002) &
            crsurf(1, k), crsurf(2, k), crsurf(3, k), typesurf(k), crystalsurf(k), crsurf(1, k + 1), crsurf(2, k + 1), &
            crsurf(3, k + 1), typesurf(k + 1), crystalsurf(k + 1)
      end do

      do i = 1, maxtype
         write (20, '(3f12.6, 2i5, f12.6)') crsurf(:, ichoicen(i)), ichoicetypen(i), crystalsurf(ichoicen(i)), &
                    & dot_product([this%dx, this%dy, this%dz], crsurf(:, ichoicen(i)))   
      end do

      ! Cleanup
      deallocate (crh, uniqueTypes)
      deallocate (this%cr)

      call move_alloc(crsurf, this%cr)
      call move_alloc(typesurf, this%iz)
      call move_alloc(crystalsurf, this%num)

      this%kk = nsurf
    !!!! FIX LATER
      if (this%control%calctype == 'I') this%nlay = maxtype - this%nbulk_bulk
    !!!!!
10002 format(3f14.6, 2i4, 3f14.6, 2i4)
10004 format(3x, "II =", i7)
      close (10)
      close (20)
   end subroutine build_surf_full

   !> @brief Build the legacy compact surface cluster.
   !> @details Chooses representative atoms from the bulk cluster for each
   !>          requested surface layer and writes the compact surface cluster
   !>          arrays consumed by newclu/buildsurf workflows.
   !> @param[inout] this Lattice object receiving compact surface state.
   module subroutine build_surf(this)
      class(lattice), intent(inout) :: this
      ! Local variables
      integer :: i, j, n, k, kk, ichoice, icont, ichoicetype
      real(rp), dimension(this%kk) :: crh, crhd
      real(rp) :: disf, disi_min, disi
      real(rp), dimension(:, :), allocatable :: crsurf
      integer, dimension(:), allocatable :: izp, no, ichoicen, ichoicetypen
      logical :: isopen

      n = int((this%zmax - this%zmin)/this%zstep) + 1
      kk = this%kk
      allocate (crsurf(3, kk), izp(kk), no(kk), ichoicen(kk), ichoicetypen(kk))
      icont = 0
      crh = 0.0d0
      crhd = 0.0d0
      inquire (unit=10, opened=isopen)
      if (isopen) then
         call g_logger%fatal('lattice%build_surf, file clust: Unit 10 is already open', __FILE__, __LINE__)
      else
         open (unit=10, file='clust')
      end if

      do i = 1, n
         disf = 1.0d0
         disi_min = sqrt(this%z(i)**2) + 0.5d0
         do k = 1, kk
            crh(k) = 0.0d0
            crh(k) = this%dx*this%cr(1, k) + this%dy*this%cr(2, k) + this%dz*this%cr(3, k)
            if (abs(crh(k) - this%z(i)) < 1.0d-6) then
               icont = icont + 1
               do j = 1, 3
                  crsurf(j, icont) = this%cr(j, k)
               end do
               crhd(icont) = this%dx*crsurf(1, icont) + this%dy*crsurf(2, icont) + this%dz*crsurf(3, icont)
               if (abs(crhd(icont) - this%z(i)) < 1.0d-6) then
                  izp(icont) = this%izpsurf(i)
                  disi = norm2(crsurf(:, icont))
                  if (disi < disi_min) then
                     disi_min = disi
                     ichoice = icont
                     ichoicetype = izp(icont)
                  end if
                  no(icont) = this%num(k)
               else
               end if
            end if
         end do
         ichoicetypen(i) = ichoicetype
         ichoicen(i) = ichoice
      end do

      if (this%control%calctype == 'S') then
         do i = 1, this%nrec
            this%irec(i) = ichoicen(i)
         end do
         do i = 1, this%ntot
            this%iu(ichoicetypen(i + this%nrec)) = ichoicen(i + this%nrec)
         end do
         do i = 1, this%nbulk
            this%ib(ichoicetypen(i + this%nrec)) = ichoicen(i + this%nrec)
         end do
      end if

      if (int(icont/2) /= icont/2.d0) icont = icont - 1

      do k = 1, icont
         write (805, *) no(k), crsurf(:, k)
      end do
      write (10, 10004) icont
      do k = 1, icont - 1, 2
         write (10, 10002) &
            crsurf(1, k), crsurf(2, k), crsurf(3, k), izp(k), no(k), crsurf(1, k + 1), crsurf(2, k + 1), &
            crsurf(3, k + 1), izp(k + 1), no(k + 1)
      end do

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%deallocate('lattice.cr', this%cr)
      call g_safe_alloc%deallocate('lattice.iz', this%iz)
      call g_safe_alloc%deallocate('lattice.num', this%num)
#else
      deallocate (this%cr, this%iz, this%num)
#endif

      call move_alloc(crsurf, this%cr)
      call move_alloc(izp, this%iz)
      call move_alloc(no, this%num)

      this%kk = icont

10002 format(3f14.6, 2i4, 3f14.6, 2i4)
10004 format(3x, "II =", i7)
      close (10)
   end subroutine build_surf

   !> @brief Build an impurity/local cluster from a bulk or surface host.
   !> @details Combines host cluster coordinates with impurity inclusions,
   !>          creates local atom/type maps, and builds neighbor tables for
   !>          impurity and defect calculations.
   !> @param[inout] this Lattice object receiving impurity cluster arrays.
   module subroutine newclu(this)
      class(lattice), intent(inout) :: this
      ! Local variables
      integer :: ndi, nnmx, ncnt, nmax
      logical :: isopen
      integer :: i, j, kk, k, ntypecount, ireccount, inclucheck
      integer, dimension(:), allocatable ::  ibulk
      integer, dimension(:), allocatable :: izpo, izp, no, nnmax, izimp, noimp
      integer, dimension(:, :), allocatable :: nn, nn2
      real(rp) :: nnscale
      real(rp), dimension(:, :), allocatable :: acr, crd, crimp
      real(rp), dimension(:), allocatable :: ctnew

      ! Set clust variables
      this%nbulk = this%nbulk_bulk + this%nlay
      this%ntype = this%nbulk_bulk + this%nlay + this%nclu
      this%nrec = this%ntype - this%nbulk
      write (*, *) this%nbulk, this%ntype, this%nrec
      ! Set the clust dimension
      kk = this%kk
      ndi = 150000
      nnmx = 200

      ! Allocating clust dimension
      if (allocated(this%ib)) deallocate (this%ib)
      if (allocated(this%irec)) deallocate (this%irec)
      if (allocated(this%iu)) deallocate (this%iu)
      !if (allocated(this%ct)) deallocate (this%ct)

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('lattice.ib', this%ib, (/this%nbulk/))
      call g_safe_alloc%allocate('lattice.iu', this%iu, (/this%ntot/))
      call g_safe_alloc%allocate('lattice.irec', this%irec, (/this%nrec/))
!      call g_safe_alloc%allocate('lattice.ct', this%ct, (/this%ntype/))
#else
      allocate (this%ib(this%nbulk), this%iu(this%ntot), this%irec(this%nrec))!, this%ct(this%ntype)) Now ct is defined at &lattice
#endif
      allocate (ctnew(this%ntype), ibulk(this%nbulk))
      allocate (nn(ndi, nnmx), nn2(ndi, nnmx))
      allocate (izpo(kk), izp(kk), no(kk), nnmax(kk), izimp(kk), noimp(kk))
      allocate (acr(kk, 7), crd(3, kk), crimp(3, kk))
      ! Setting ct values for impurity
      !this%ct(:) = this%alat+0.1d0
      ! Identify impurity atoms from ´inclu´
      do i = 1, kk
         izpo(i) = this%iz(i)
      end do
      ntypecount = this%nbulk
      inclucheck = 0
      do j = 1, this%nclu
         do i = 1, kk
            if (abs(this%cr(1, i) - this%inclu(j, 1)) < 1.0d-6 .and. &
                abs(this%cr(2, i) - this%inclu(j, 2)) < 1.0d-6 .and. &
                abs(this%cr(3, i) - this%inclu(j, 3)) < 1.0d-6) then
               inclucheck = inclucheck + 1
               ntypecount = ntypecount + 1
               this%iz(i) = ntypecount
            end if
         end do
      end do

      if (inclucheck /= this%nclu) then
         call g_logger%fatal('Atoms chosen for impurity were not found inside the clust. Please, check the inclu input', __FILE__, __LINE__)
      end if

      acr = 0.0d0
      inquire (unit=10, opened=isopen)
      if (isopen) then
         call g_logger%fatal('lattice%newclu, file clust: Unit 10 is already open', __FILE__, __LINE__)
      else
         open (unit=10, file='clust')
      end if

      open (unit=11, file='outnewclu')

      do k = 1, kk
         do i = 1, 3
            acr(k, i) = this%cr(i, k)
            acr(k, 6) = acr(k, 6) + (this%cr(i, k) - this%inclu(1, i))**2
         end do
         acr(k, 4) = (this%iz(k))
         acr(k, 5) = (this%num(k))
         acr(k, 7) = (izpo(k))
      end do
      call bubble(this%nclu, this%nclu, acr(1:this%nclu, 1:7), 7, 4)
      call bubble(kk - this%nclu, kk - this%nclu, acr(this%nclu + 1:kk, 1:7), 7, 6)

      write (10, 10004) kk
      do k = 1, kk - 1, 2
         write (10, 10002) &
            ( &
            acr(k + i - 1, 1), acr(k + i - 1, 2), acr(k + i - 1, 3), int(acr(k + i - 1, 4)) &
            , int(acr(k + i - 1, 5)), i=1, 2)
      end do

      rewind (10)

#ifdef USE_MPI
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif

      call leia(this%alat, kk, crd, izp, no, 10)

      close (10)

      nnscale = 0.95d0

      ctnew(:) = this%ct(:)*nnscale

      call this%nncal(ctnew, crd, 3, kk, izp, nn, ndi, nnmx, mapa, this%ntot)
      nnmx = 200
      call this%nncal(this%ct, crd, 3, kk, izp, nn2, ndi, nnmx, mapa, this%ntot)
      nnmx = 200

      do i = 1, kk
         if (acr(i, 4) > this%nbulk .and. acr(i, 4) <= this%ntype) then
            do j = 2, nn2(i, 1)
               if (acr(nn2(i, j), 4) <= this%nbulk) acr(nn2(i, j), 4) = 2000 + acr(i, 7)
            end do
         end if
      end do
      do i = 1, kk
         if (acr(i, 4) > this%nbulk .and. acr(I, 4) <= this%ntype) then
            do j = 2, nn(i, 1)
               if (acr(nn(i, j), 4) <= this%nbulk .or. acr(nn(i, j), 4) > 2000) acr(nn(i, j), 4) = 1000 + acr(i, 7)
            end do
         end if
      end do
      do i = 1, kk
         if (acr(i, 4) == 1) acr(i, 4) = 4000 + acr(I, 7)
      end do
      do i = 1, kk
         if (acr(i, 4) <= this%nbulk) acr(i, 4) = 3000 + acr(I, 7)
      end do
      call bubble(kk, kk, acr(1:kk, 1:7), 7, 4)
      ncnt = 0
      do I = 1, kk
         if (acr(i, 4) < 2000) ncnt = ncnt + 1
      end do
      do i = 1, kk
         if (acr(i, 4) > this%ntype) acr(i, 4) = acr(i, 7)
      end do
      call bubble(kk - ncnt, kk - ncnt, acr(ncnt + 1:(kk), 1:7), 7, 6)

      open (unit=10, file='clust')
      write (10, 10004) kk
      do k = 1, kk, 2
         write (10, 10002) &
            ( &
            acr(k + i - 1, 1), acr(k + i - 1, 2), acr(k + i - 1, 3), int(acr(k + i - 1, 4)) &
            , int(acr(k + i - 1, 5)), i=1, 2)
      end do
      rewind (10)
#ifdef USE_MPI
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
      call leia(this%alat, kk, crd, izp, no, 10)
      call this%nncal(this%ct, crd, 3, kk, izp, nn, ndi, nnmx, mapa, this%ntot)
      call outmap(11, izp, nn, no, ndi, nnmx, this%nrec)
      close (10)
      nmax = 0
      do i = 1, this%nrec
         do j = 2, nn(i, 1)
            if (nmax < nn(i, j)) nmax = nn(i, j)
         end do
      end do
      write (11, *) "--Info-for-bulcri------------------"
      !
      if (allocated(this%chargetrf_type)) deallocate (this%chargetrf_type)
      allocate (this%chargetrf_type(ncnt))
      !
      do k = 1, ncnt
         write (11, '(i5)') int(acr(k, 7))   !, int(ACR(K, 7))
         this%chargetrf_type(k) = int(acr(k, 7))
      end do
      this%nbas = ncnt
      this%nmax = nmax
      ! Temporary hack for full HALL
      !this%nmax = kk
      write (11, *) "--Info-for-control-----------------"
      write (11, '(a6, i4, a6, i6, a6, i6, a6, i6)') 'NTYPE=', this%ntype, 'NMAX=', this%nmax, 'NBAS=', this%nbas, 'NREC=', this%nrec

      nnmax = 0
      do i = nmax + 1, kk
         if (nnmax(izp(i)) < nn(i, 1)) then
            nnmax(izp(i)) = nn(i, 1)
            ibulk(izp(i)) = i
         end if
      end do
      write (11, '(50i7)') (nnmax(k), k=1, this%nbulk)
      write (11, '(50i7)') (ibulk(k), k=1, this%nbulk)
      do i = 1, this%nbulk
         do j = 2, nn(ibulk(i), 1)
            if (nn(ibulk(i), j) == 0) write (*, *) "Warning! Atom", ibulk(i), " is probably wrong."
         end do
      end do

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%deallocate('lattice.cr', this%cr)
      call g_safe_alloc%deallocate('lattice.cr', this%cr)
      call g_safe_alloc%deallocate('lattice.iz', this%iz)
      call g_safe_alloc%deallocate('lattice.num', this%num)
#else
      deallocate (this%cr, this%iz, this%num)
#endif

      do k = 1, kk
         crimp(1:3, k) = acr(k, 1:3)
      end do

      izimp(:) = int(acr(:, 4))
      noimp(:) = int(acr(:, 5))

      call move_alloc(izimp, this%iz)
      call move_alloc(noimp, this%num)
      call move_alloc(crimp, this%cr)

      do i = 1, this%nbulk
         this%ib(i) = ibulk(i)
      end do

      do i = 1, this%ntot
         this%iu(i) = ibulk(i)
      end do

      ireccount = 0
      do j = 1, this%nclu
         do i = 1, kk
            if (abs(this%cr(1, i) - this%inclu(j, 1)) < 1.0d-6 .and. &
                abs(this%cr(2, i) - this%inclu(j, 2)) < 1.0d-6 .and. &
                abs(this%cr(3, i) - this%inclu(j, 3)) < 1.0d-6) then
               ireccount = ireccount + 1
               this%irec(ireccount) = i
            end if
         end do
      end do

      close (11)
10002 format(3f14.6, 2i4, 3f14.6, 2i4)
10004 format(3x, "II =", i7)
   end subroutine newclu

   !> @brief Build atom-list metadata for self-consistent atom types.
   !> @details Fills atlist/chargetrf_type-style mappings used by impurity and
   !>          local-cluster charge transfer paths.
   !> @param[inout] this Lattice object whose atom-list arrays are updated.
   module subroutine atomlist(this)
      class(lattice), intent(inout) :: this
      ! Local variables
      integer :: i, j, itype
      real(rp) :: mom_tmp(3)

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('lattice.atlist', this%atlist, (/this%ntype/))
      call g_safe_alloc%allocate('lattice.ham_i', this%ham_i, (/this%kk/))
#else
      allocate (this%atlist(this%ntype))
      allocate (this%ham_i(this%kk))
#endif

      do i = 1, this%kk
         this%ham_i(i) = i
      end do

      j = 0
      do i = 1, this%nbulk
         j = j + 1
         this%atlist(j) = this%ib(i)
      end do
      do i = 1, this%nrec
         j = j + 1
         this%atlist(j) = this%irec(i)
      end do
      call this%load_symbolic_atoms_if_needed()

      
      write (805, *) this%kk
      write (805, *)
      do i = 1, this%kk
         itype = 1
         if (allocated(this%iz)) then
            if (i <= size(this%iz)) then
               if (this%iz(i) >= 1 .and. this%iz(i) <= size(this%symbolic_atoms)) itype = this%iz(i)
            end if
         end if
         mom_tmp(:) = [0.0_rp, 0.0_rp, 1.0_rp]
         if (allocated(this%symbolic_atoms(itype)%potential%mom)) then
            if (size(this%symbolic_atoms(itype)%potential%mom) >= 3) then
               mom_tmp(:) = this%symbolic_atoms(itype)%potential%mom(1:3)
            end if
         end if
         write (805, '(A6,6F16.6)') (elem_var(int(this%symbolic_atoms(itype)%element%atomic_number))), &
            this%cr(:, i), mom_tmp(:)
      end do
   end subroutine atomlist

   !> @brief Select cluster atoms inside a primitive-cell volume.
   !> @details Tests positions relative to a central atom against the three
   !>          primitive vectors and returns the atom indices that lie inside.
   !> @param[inout] this Lattice object providing volume-test helper methods.
   !> @param[in] cr Cluster coordinates.
   !> @param[in] num Cluster atom labels.
   !> @param[in] num_atoms Number of atoms in cr/num.
   !> @param[in] central_atom Index of the atom used as the volume origin.
   !> @param[in] a1 First primitive vector.
   !> @param[in] a2 Second primitive vector.
   !> @param[in] a3 Third primitive vector.
   !> @param[in] plane_constant Boundary tolerance/plane constant.
   !> @param[out] atoms_in_volume Atom indices found inside the volume.
   !> @param[out] atom_count Number of returned atoms.
   module subroutine check_atoms_in_volume(this, cr, num, num_atoms, central_atom, a1, a2, a3, plane_constant, atoms_in_volume, atom_count)
       class(lattice), intent(inout) :: this
       real(rp), intent(in) :: cr(3, num_atoms), a1(3), a2(3), a3(3)
       integer, intent(in) :: num(num_atoms), num_atoms, central_atom
       real(rp), intent(in) :: plane_constant
       integer, allocatable, intent(out) :: atoms_in_volume(:)
       integer, intent(out) :: atom_count
       real(rp) :: relative_pos(3)
       logical :: inside
       integer :: i
   
       ! Initialize array for atoms inside the primitive cell volume
       allocate(atoms_in_volume(num_atoms))
       atom_count = 0
       atoms_in_volume = 0
       ! Loop over all atoms to find those inside the primitive cell volume
       do i = 1, num_atoms
          ! Calculate the position relative to the central atom
          relative_pos = cr(:, i) - cr(:, central_atom)
          ! Check if the atom is within the parallelepiped defined by a1, a2, and a3
          call this%check_within_volume(relative_pos, a1, a2, a3, inside)
   
          if (inside) then
             atom_count = atom_count + 1
             atoms_in_volume(atom_count) = i
          end if
       end do
       !write(*,*) 'atoms_in_volume', atoms_in_volume
       ! Resize atoms_in_volume array to the actual number of atoms found
       !if (atom_count > 0) then
       !   deallocate(atoms_in_volume); allocate(atoms_in_volume(atom_count))
       !else
       !   deallocate(atoms_in_volume)
       !end if
   
   end subroutine check_atoms_in_volume

   !> @brief Test whether a relative position is inside a primitive parallelepiped.
   !> @param[inout] this Lattice object providing the type-bound helper context.
   !> @param[in] relative_pos Position measured from the cell origin.
   !> @param[in] a1 First primitive vector.
   !> @param[in] a2 Second primitive vector.
   !> @param[in] a3 Third primitive vector.
   !> @param[out] inside True if the point is inside the primitive volume.
   module subroutine check_within_volume(this, relative_pos, a1, a2, a3, inside)
       class(lattice), intent(inout) :: this
       real(rp), intent(in) :: relative_pos(3), a1(3), a2(3), a3(3)
       logical, intent(out) :: inside
       real(rp) :: dot11, dot12, dot13, dot22, dot23, dot33
       real(rp) :: dot1r, dot2r, dot3r, inv_denom, u, v, w
   
       ! Calculate dot products between the vectors
       dot11 = dot_product(a1, a1)
       dot12 = dot_product(a1, a2)
       dot13 = dot_product(a1, a3)
       dot22 = dot_product(a2, a2)
       dot23 = dot_product(a2, a3)
       dot33 = dot_product(a3, a3)
       dot1r = dot_product(a1, relative_pos)
       dot2r = dot_product(a2, relative_pos)
       dot3r = dot_product(a3, relative_pos)
   
       ! Calculate inverse of the denominator for the linear combination
       inv_denom = 1.0_rp / (dot11 * (dot22 * dot33 - dot23 * dot23) &
                           - dot12 * (dot12 * dot33 - dot23 * dot13) &
                           + dot13 * (dot12 * dot23 - dot22 * dot13))
   
       ! Calculate coordinates (u, v, w) in the basis defined by a1, a2, and a3
       u = ((dot22 * dot33 - dot23 * dot23) * dot1r + &
           (dot13 * dot23 - dot12 * dot33) * dot2r + &
           (dot12 * dot23 - dot13 * dot22) * dot3r) * inv_denom
   
       v = ((dot13 * dot23 - dot12 * dot33) * dot1r + &
           (dot11 * dot33 - dot13 * dot13) * dot2r + &
           (dot12 * dot13 - dot11 * dot23) * dot3r) * inv_denom
   
       w = ((dot12 * dot23 - dot13 * dot22) * dot1r + &
           (dot12 * dot13 - dot11 * dot23) * dot2r + &
           (dot11 * dot22 - dot12 * dot12) * dot3r) * inv_denom
   
       ! Check if the atom is inside the parallelepiped
       inside = (u >= 0.0_rp .and. u <= 1.0_rp .and. &
                 v >= 0.0_rp .and. v <= 1.0_rp .and. &
                 w >= 0.0_rp .and. w <= 1.0_rp)
   end subroutine check_within_volume

   !> @brief Return unique integer structure/type labels.
   !> @param[inout] this Lattice object providing the type-bound helper context.
   !> @param[in] num Input labels.
   !> @param[in] num_atoms Number of active labels.
   !> @param[out] unique_nums Unique labels in first-seen order.
   module subroutine find_unique_struct(this, num, num_atoms, unique_nums)
       class(lattice), intent(inout) :: this
       integer, intent(in) :: num(:), num_atoms
       integer, allocatable, intent(out) :: unique_nums(:)
       integer, allocatable :: temp_nums(:)
       integer :: i, j, num_unique
       logical :: found
   
       ! Allocate temporary array for unique numbers
       allocate(temp_nums(num_atoms))
       num_unique = 0
   
       ! Loop over all atoms to find unique structure types
       do i = 1, num_atoms
          found = .false.
          do j = 1, num_unique
             if (num(i) == temp_nums(j)) then
                found = .true.
                exit
             end if
          end do
          if (.not. found) then
             num_unique = num_unique + 1
             temp_nums(num_unique) = num(i)
          end if
       end do
   
       ! Resize array to the actual number of unique nums
       allocate(unique_nums(num_unique))
       unique_nums(:) = temp_nums(1:num_unique)
       ! Deallocate temporary array
       deallocate(temp_nums)
   end subroutine find_unique_struct

   !> @brief Identify symmetry-unique atoms inside a primitive volume.
   !> @details Removes atoms related by primitive translations from an input
   !>          volume selection and returns one representative per unique site.
   !> @param[inout] this Lattice object providing helper context.
   !> @param[in] cr Cluster coordinates.
   !> @param[in] num_atoms Number of atoms in cr.
   !> @param[in] atoms_in_volume Candidate atom indices.
   !> @param[in] atom_count Number of candidate atoms.
   !> @param[in] a1 First primitive vector.
   !> @param[in] a2 Second primitive vector.
   !> @param[in] a3 Third primitive vector.
   !> @param[out] unique_atoms Representative atom indices.
   !> @param[out] unique_atom_count Number of representatives.
   module subroutine identify_unique_atoms(this, cr, num_atoms, atoms_in_volume, atom_count, a1, a2, a3, unique_atoms, unique_atom_count)
       class(lattice), intent(inout) :: this
       real(rp), intent(in) :: cr(3, num_atoms), a1(3), a2(3), a3(3)
       integer, intent(in) :: num_atoms, atoms_in_volume(:), atom_count
       integer, allocatable, intent(out) :: unique_atoms(:)
       integer, intent(out) :: unique_atom_count
       logical :: found, is_transformed
       real(rp) :: trans_atom(3), delta(3)
       integer :: i, j, k, n, m, p
       integer, allocatable :: temp_unique_atoms(:)
   
       ! Initialize temporary array for unique atoms list
       allocate(temp_unique_atoms(atom_count))
       unique_atom_count = 0
       temp_unique_atoms = -1
       ! First pass to identify all unique atoms within the primitive cell volume
       do i = 1, atom_count
           found = .false.
   
           ! Check if this atom can be generated by translating another atom
           do j = 1, unique_atom_count
               ! Translate atom by all combinations of a1, a2, and a3 within the cell
               do k = -1, 1
                   do n = -1, 1
                       do p = -1, 1
                           trans_atom = cr(:, temp_unique_atoms(j)) + k * a1 + n * a2 + p * a3
                           delta = cr(:, atoms_in_volume(i)) - trans_atom
                           if (norm2(delta) < 1.0d-6) then
                               found = .true.
                               exit
                           end if
                       end do
                       if (found) exit
                   end do
                   if (found) exit
               end do
           end do
   
           ! If the atom is not redundant, add it to the temporary list of unique atoms
           if (.not. found) then
               unique_atom_count = unique_atom_count + 1
               temp_unique_atoms(unique_atom_count) = atoms_in_volume(i)
           end if
       end do
   
       ! Allocate the final unique_atoms array to the correct size
       allocate(unique_atoms(unique_atom_count))
       unique_atoms(:) = temp_unique_atoms(1:unique_atom_count)
       ! Deallocate temporary array
       deallocate(temp_unique_atoms)
   end subroutine identify_unique_atoms

   !> @brief Build representative neighbor vectors for each atom type.
   !> @details Compares cluster coordinates around representative atoms and
   !>          fills neighbor maps plus displacement sets used by structb.
   !> @param[inout] this Lattice object providing tolerance and helper context.
   !> @param[in] crd Cluster coordinates.
   !> @param[in] no Atom type labels for cluster atoms.
   !> @param[in] iu Representative atom indices.
   !> @param[inout] nn Neighbor table to fill.
   !> @param[in] nat Number of cluster atoms.
   !> @param[in] ntot Number of representative atoms.
   !> @param[in] nomx Number of representative/type slots.
   !> @param[in] ndi Leading dimension/capacity of crd.
   !> @param[in] nnmx Maximum neighbor slots.
   !> @param[inout] set Neighbor displacement vectors.
   !> @param[inout] idnn Neighbor identifier list.
   !> @param[out] ret Return/status vector used by legacy callers.
   module subroutine remd(this, crd, no, iu, nn, nat, ntot, nomx, ndi, nnmx, set, idnn, ret)
      implicit none
      ! Inputs
      class(lattice), intent(inout) :: this
      integer, intent(in) :: nat, ndi, nnmx, nomx, ntot
      integer, dimension(ndi), intent(in) :: no
      real(rp), dimension(3, ndi), intent(in) :: crd
      integer, dimension(nomx), intent(in) :: iu
      ! Output
      integer, dimension(nnmx), intent(inout) :: idnn
      integer, dimension(ndi, nnmx), intent(inout) :: nn
      real(rp), dimension(3), intent(out) :: ret
      real(rp), dimension(3, nomx, nnmx), intent(inout) :: set
      ! Local variables
      integer :: i, ii, iii, imax, inn, ino, j, jj, jnn, jsz, k, la, lk, lm, m, n
      real(rp) :: a1, a2, a3, aaa, eps
      !-BUILDS VECTORS SET(3, NOMX, NNMX) CONNECTING NEIGHBORS OF EACH TYPE NO-
      do i = 1, ntot
         la = iu(i)
         jsz = nn(la, 1)
         do j = 2, jsz
            jj = nn(la, j)
            if (this%pbc) then 
               call this%f_wrap_coord_diff(nat,crd,la,jj,set(:,i,j))
            else
               do m = 1, 3
                  set(m, i, j) = crd(m, la) - crd(m, jj)
               end do
            end if
         end do
      end do
      do i = 1, nat
         n = no(i)
         do lk = 1, ntot
            ino = iu(lk)
            lm = no(ino)
            if (n == lm) goto 1000
         end do
         write (17, 10003)
1000     imax = nn(ino, 1)
         do iii = 1, imax
            idnn(iii) = 0
         end do
         jsz = nn(i, 1)
         do j = 2, jsz
            jj = nn(i, j)
            if (this%pbc) then
               call this%f_wrap_coord_diff(nat,crd,i,jj,ret)
            else
               do m = 1, 3
                  ret(m) = crd(m, i) - crd(m, jj)
               end do
            end if
            !----------FINDS EQUIVALENT VECTOR------------------------
            eps = .0001
            do ii = 2, imax
               a1 = ret(1) - set(1, n, ii)
               a2 = ret(2) - set(2, n, ii)
               a3 = ret(3) - set(3, n, ii)
               aaa = a1**2 + a2**2 + a3**2
               if (aaa < eps) goto 1100
            end do
            goto 1200
1100        k = ii
            idnn(k) = jj
         end do
         nn(i, 1) = imax
         do j = 2, imax
            nn(i, j) = idnn(j)
         end do
      end do
      do inn = 1, nat
         write (12) (nn(inn, jnn), jnn=1, nn(inn, 1))
      end do
      return
1200  write (17, 10002)
      print *, 'atom:', i
      stop
      !  9  WRITE(6, 223)I, SET(1, I, J), SET(2, I, J), SET(3, I, J)
10000 format(i5, 3f9.4)
      !--REORDERS NEIGHBORS FOR I LARGER THAN NTOT ACCORDING TO TYPICAL NO--
10001 format(3i5)
10002 format(" VECTOR   NOT FOUND ")
10003 format(" TYPE NO  NOT FOUND ")
   end subroutine remd

   !> @brief Compute the minimum-image coordinate difference between two atoms.
   !> @param[inout] this Lattice object containing periodic-boundary settings.
   !> @param[in] Natom Number of atoms in coord.
   !> @param[in] coord Atomic coordinates.
   !> @param[in] i_atom First atom index.
   !> @param[in] j_atom Second atom index.
   !> @param[out] cdiff Minimum-image displacement from i_atom to j_atom.
   module subroutine f_wrap_coord_diff(this,Natom,coord,i_atom,j_atom,cdiff)
      implicit none
      class(lattice), intent(inout) :: this
      integer, intent(in) :: Natom
      real(rp), dimension(3,Natom), intent(in) :: coord
      integer, intent(in) :: i_atom
      integer, intent(in) :: j_atom
      real(rp), dimension(3), intent(out) :: cdiff
      !
      real(rp), dimension(3) :: odiff, oshift, mdiff
      integer :: x,y,z
      integer :: xmin,xmax,ymin,ymax,zmin,zmax
      !
      odiff=coord(:,j_atom) - coord(:,i_atom)
      !
      xmax=0;xmin=0;ymax=0;ymin=0;zmax=0;zmin=0
      if(this%b1)then
         xmax=1
         xmin=-1
      end if
      if(this%b2)then
         ymax=1
         ymin=-1
      end if
      if(this%b3)then
         zmax=1
         zmin=-1
      end if
      
      mdiff=odiff
      do z=zmin,zmax
         do y=ymin,ymax
            do x=xmin,xmax
               oshift = odiff + x*(this%n1)*this%a(:, 1)*this%alat& 
                              + y*(this%n2)*this%a(:, 2)*this%alat&
                              + z*(this%n3)*this%a(:, 3)*this%alat
               if(norm2(oshift)<norm2(mdiff))  mdiff = oshift
            end do
         end do
      end do
      cdiff=mdiff
      return
      !
   end subroutine f_wrap_coord_diff

   !> @brief Build the atom neighbor table using shell cutoffs.
   !> @details Uses spatial binning and optional periodic wrapping to find
   !>          neighbors within the configured cutoff shells and populate NN.
   !> @param[inout] this Lattice object containing cell/PBC state.
   !> @param[inout] ct Shell cutoff radii.
   !> @param[in] crd Atomic coordinates.
   !> @param[in] ndim Coordinate leading dimension.
   !> @param[in] nat Number of atoms.
   !> @param[in] izp Atomic numbers/labels.
   !> @param[inout] nn Neighbor table.
   !> @param[in] nd Number of atoms represented in nn.
   !> @param[inout] nm Maximum/actual neighbor slots.
   !> @param[in] ngbr Legacy neighbor-shell classification function.
   !> @param[in] ntot Number of representative atoms.
   module subroutine nncal(this,ct, crd, ndim, nat, izp, nn, nd, nm, ngbr, ntot)
      implicit none
      class(lattice), intent(inout) :: this
      ! Input
      integer, intent(in) :: NAT, ND, NDIM, NTOT
      integer, dimension(NAT), intent(in) :: IZP
      real(rp), dimension(NDIM, NAT), intent(in) :: CRD
      ! Output
      integer, intent(inout) :: NM
      integer, dimension(ND, NM), intent(inout) :: NN
      real(rp), dimension(50), intent(inout) :: CT
      ! External function
      integer, external :: NGBR
      ! Intrinsic function
      intrinsic ABS, MAX, MIN, FLOOR, NINT
      ! Local variables
      integer :: I, IADD, ID, II, IIP, ILJ, J, JJP, L, NNMAX
      integer :: NX, NY, NZ, NBIN, BX, BY, BZ, DBX, DBY, DBZ, BIN_ID
      integer :: CAP, CAND_COUNT, K, IX, IY, IZ, RX, RY, RZ
      real(rp) :: R2, RCUT, RCUT2, DETC
      real(rp), dimension(3) :: DDUM
      real(rp), dimension(3) :: MINC, MAXC, SPAN, BINW, DS
      real(rp), dimension(3, 3) :: CELL, CELL_INV
      real(rp), dimension(3, 2) :: CROSS_TMP
      real(rp), dimension(NM) :: DUM
      real(rp), allocatable :: FRAC(:, :)
      integer, allocatable :: HEAD(:), NEXT_ATOM(:), BIN_X(:), BIN_Y(:), BIN_Z(:), CANDIDATES(:)

      CAP = NM
      NN(:, :) = 0
      NN(:, 1) = 1

      NNMAX = 0
      IADD = 1
      if (IZP(1) < 0) then
         NN(1, 1) = 1
      end if
      do II = 1, NAT
         if (IZP(II) < 0) goto 1000
      end do
      NN(1, 1) = 1
      do I = 1, NAT
         do J = 2, NM
            NN(I, J) = 0
         end do
      end do
      IADD = 0
      II = 2
1000  if (II <= 1) then
         II = 2
      end if

      RCUT = CT(1)
      RCUT2 = RCUT*RCUT
      if (RCUT2 <= 0.0_rp) then
         NM = 0
         return
      end if

      if (this%pbc .and. this%b1 .and. this%b2 .and. this%b3) then
         CELL(:, 1) = real(this%n1, rp)*this%a(:, 1)*this%alat
         CELL(:, 2) = real(this%n2, rp)*this%a(:, 2)*this%alat
         CELL(:, 3) = real(this%n3, rp)*this%a(:, 3)*this%alat
         CELL_INV = inverse_3x3(CELL)
         DETC = abs(determinant(CELL))

         CROSS_TMP(:, 1) = [ &
            CELL(2, 2)*CELL(3, 3) - CELL(3, 2)*CELL(2, 3), &
            CELL(3, 2)*CELL(1, 3) - CELL(1, 2)*CELL(3, 3), &
            CELL(1, 2)*CELL(2, 3) - CELL(2, 2)*CELL(1, 3) ]
         CROSS_TMP(:, 2) = [ &
            CELL(2, 3)*CELL(3, 1) - CELL(3, 3)*CELL(2, 1), &
            CELL(3, 3)*CELL(1, 1) - CELL(1, 3)*CELL(3, 1), &
            CELL(1, 3)*CELL(2, 1) - CELL(2, 3)*CELL(1, 1) ]
         DDUM = [ &
            CELL(2, 1)*CELL(3, 2) - CELL(3, 1)*CELL(2, 2), &
            CELL(3, 1)*CELL(1, 2) - CELL(1, 1)*CELL(3, 2), &
            CELL(1, 1)*CELL(2, 2) - CELL(2, 1)*CELL(1, 2) ]

         SPAN(1) = DETC/max(sqrt(sum(CROSS_TMP(:, 1)**2)), tiny(1.0_rp))
         SPAN(2) = DETC/max(sqrt(sum(CROSS_TMP(:, 2)**2)), tiny(1.0_rp))
         SPAN(3) = DETC/max(sqrt(sum(DDUM(:)**2)), tiny(1.0_rp))
         NX = max(1, int(SPAN(1)/RCUT))
         NY = max(1, int(SPAN(2)/RCUT))
         NZ = max(1, int(SPAN(3)/RCUT))
         NBIN = NX*NY*NZ

         allocate(FRAC(3, NAT), HEAD(NBIN), NEXT_ATOM(NAT), BIN_X(NAT), BIN_Y(NAT), BIN_Z(NAT), CANDIDATES(NAT))
         HEAD = 0
         NEXT_ATOM = 0
         do I = 1, NAT
            FRAC(:, I) = matmul(CELL_INV, CRD(:, I))
            do L = 1, 3
               FRAC(L, I) = FRAC(L, I) - floor(FRAC(L, I))
            end do
            BIN_X(I) = 1 + min(NX - 1, int(FRAC(1, I)*NX))
            BIN_Y(I) = 1 + min(NY - 1, int(FRAC(2, I)*NY))
            BIN_Z(I) = 1 + min(NZ - 1, int(FRAC(3, I)*NZ))
            BIN_ID = ((BIN_Z(I) - 1)*NY + (BIN_Y(I) - 1))*NX + BIN_X(I)
            NEXT_ATOM(I) = HEAD(BIN_ID)
            HEAD(BIN_ID) = I
         end do

         do I = II, NAT
            if (IZP(I)*IADD > 0) cycle
            CAND_COUNT = 0
            BX = BIN_X(I); BY = BIN_Y(I); BZ = BIN_Z(I)
            RX = merge(1, 0, NX > 1)
            RY = merge(1, 0, NY > 1)
            RZ = merge(1, 0, NZ > 1)
            do DBZ = -RZ, RZ
               IZ = modulo(BZ - 1 + DBZ, NZ) + 1
               do DBY = -RY, RY
                  IY = modulo(BY - 1 + DBY, NY) + 1
                  do DBX = -RX, RX
                     IX = modulo(BX - 1 + DBX, NX) + 1
                     BIN_ID = ((IZ - 1)*NY + (IY - 1))*NX + IX
                     J = HEAD(BIN_ID)
                     do while (J > 0)
                        if (J < I) then
                           DS(:) = FRAC(:, J) - FRAC(:, I)
                           DS(:) = DS(:) - real(nint(DS(:)), rp)
                           DDUM(:) = matmul(CELL, DS)
                           R2 = dot_product(DDUM, DDUM)
                           if (R2 < RCUT2) then
                              CAND_COUNT = CAND_COUNT + 1
                              CANDIDATES(CAND_COUNT) = J
                           end if
                        end if
                        J = NEXT_ATOM(J)
                     end do
                  end do
               end do
            end do

            call sort_integer_list(CANDIDATES, CAND_COUNT)
            do K = 1, CAND_COUNT
               J = CANDIDATES(K)
               ID = NN(I, 1) + 1
               NN(I, 1) = ID
               NN(I, ID) = J
               NNMAX = MAX(NNMAX, ID)
               ID = NN(J, 1) + 1
               NN(J, 1) = ID
               NN(J, ID) = I
               NNMAX = MAX(NNMAX, ID)
               if (NNMAX > CAP) then
                  write (6, '(" TOO MANY NEIGHBOURS")')
                  write (6, '(" NEIGHBOUR MAP AS FAR AS", i6, "TH SITE")') I
                  write (6, *) NNMAX, ID, CAP
                  stop
               end if
            end do
         end do

         deallocate(FRAC, HEAD, NEXT_ATOM, BIN_X, BIN_Y, BIN_Z, CANDIDATES)
      else if (.not. this%pbc) then
         do L = 1, 3
            MINC(L) = minval(CRD(L, 1:NAT))
            MAXC(L) = maxval(CRD(L, 1:NAT))
            SPAN(L) = max(MAXC(L) - MINC(L), RCUT)
         end do
         NX = max(1, int(SPAN(1)/RCUT))
         NY = max(1, int(SPAN(2)/RCUT))
         NZ = max(1, int(SPAN(3)/RCUT))
         BINW(1) = SPAN(1)/real(NX, rp)
         BINW(2) = SPAN(2)/real(NY, rp)
         BINW(3) = SPAN(3)/real(NZ, rp)
         NBIN = NX*NY*NZ

         allocate(HEAD(NBIN), NEXT_ATOM(NAT), BIN_X(NAT), BIN_Y(NAT), BIN_Z(NAT), CANDIDATES(NAT))
         HEAD = 0
         NEXT_ATOM = 0
         do I = 1, NAT
            BIN_X(I) = 1 + min(NX - 1, int((CRD(1, I) - MINC(1))/BINW(1)))
            BIN_Y(I) = 1 + min(NY - 1, int((CRD(2, I) - MINC(2))/BINW(2)))
            BIN_Z(I) = 1 + min(NZ - 1, int((CRD(3, I) - MINC(3))/BINW(3)))
            BIN_ID = ((BIN_Z(I) - 1)*NY + (BIN_Y(I) - 1))*NX + BIN_X(I)
            NEXT_ATOM(I) = HEAD(BIN_ID)
            HEAD(BIN_ID) = I
         end do

         do I = II, NAT
            if (IZP(I)*IADD > 0) cycle
            CAND_COUNT = 0
            BX = BIN_X(I); BY = BIN_Y(I); BZ = BIN_Z(I)
            do DBZ = max(-1, 1 - BZ), min(1, NZ - BZ)
               IZ = BZ + DBZ
               do DBY = max(-1, 1 - BY), min(1, NY - BY)
                  IY = BY + DBY
                  do DBX = max(-1, 1 - BX), min(1, NX - BX)
                     IX = BX + DBX
                     BIN_ID = ((IZ - 1)*NY + (IY - 1))*NX + IX
                     J = HEAD(BIN_ID)
                     do while (J > 0)
                        if (J < I) then
                           DDUM(:) = CRD(:, J) - CRD(:, I)
                           R2 = dot_product(DDUM, DDUM)
                           if (R2 < RCUT2) then
                              CAND_COUNT = CAND_COUNT + 1
                              CANDIDATES(CAND_COUNT) = J
                           end if
                        end if
                        J = NEXT_ATOM(J)
                     end do
                  end do
               end do
            end do

            call sort_integer_list(CANDIDATES, CAND_COUNT)
            do K = 1, CAND_COUNT
               J = CANDIDATES(K)
               ID = NN(I, 1) + 1
               NN(I, 1) = ID
               NN(I, ID) = J
               NNMAX = MAX(NNMAX, ID)
               ID = NN(J, 1) + 1
               NN(J, 1) = ID
               NN(J, ID) = I
               NNMAX = MAX(NNMAX, ID)
               if (NNMAX > CAP) then
                  write (6, '(" TOO MANY NEIGHBOURS")')
                  write (6, '(" NEIGHBOUR MAP AS FAR AS", i6, "TH SITE")') I
                  write (6, *) NNMAX, ID, CAP
                  stop
               end if
            end do
         end do

         deallocate(HEAD, NEXT_ATOM, BIN_X, BIN_Y, BIN_Z, CANDIDATES)
      else
         do I = II, NAT
            if (IZP(I)*IADD <= 0) then
               NN(I, 1) = 1
               IIP = IZP(I)
               IIP = ABS(IIP)
               ILJ = I - 1
               do J = 1, ILJ
                  JJP = IZP(J)
                  JJP = ABS(JJP)
                  R2 = 0.0
                  if (this%pbc) then
                     call this%f_wrap_coord_diff(nat, crd, i, j, ddum)
                     r2 = sum(ddum(:)**2)
                  else        
                     do L = 1, 3
                        DDUM(L) = CRD(L, I) - CRD(L, J)
                        R2 = R2 + DDUM(L)*DDUM(L)
                     end do
                  end if
                  ID = NGBR(IIP, JJP, R2, DUM, CT)
                  if (ID /= 0) then
                     ID = NN(I, 1) + 1
                     NN(I, 1) = ID
                     NN(I, ID) = J
                     NNMAX = MAX(NNMAX, ID)
                     ID = NN(J, 1) + 1
                     NN(J, 1) = ID
                     NN(J, ID) = I
                     NNMAX = MAX(NNMAX, ID)
                     if (NNMAX > NM) then
                        write (6, '(" TOO MANY NEIGHBOURS")')
                        write (6, '(" NEIGHBOUR MAP AS FAR AS", i6, "TH SITE")') I
                        write (6, *) NNMAX, ID, NM
                        stop
                     end if
                  end if
               end do
            end if
         end do
      end if
      nm = nnmax
      return

contains

      subroutine sort_integer_list(list, n)
         integer, intent(inout) :: list(:)
         integer, intent(in) :: n
         integer :: a, b, value

         do a = 2, n
            value = list(a)
            b = a - 1
            do while (b >= 1)
               if (list(b) <= value) exit
               list(b + 1) = list(b)
               b = b - 1
            end do
            list(b + 1) = value
         end do
      end subroutine sort_integer_list
   end subroutine

   !> @brief Read legacy cluster coordinate records from a unit.
   !> @param[in] alat Lattice parameter used for coordinate scaling.
   !> @param[in] nndim Maximum number of cluster atoms.
   !> @param[inout] cr Coordinate array to fill.
   !> @param[inout] iz Atomic numbers/labels.
   !> @param[inout] n Atom-number labels.
   !> @param[in] ip Input unit.
   module subroutine leia(alat, nndim, cr, iz, n, ip)
      implicit none
      ! Input
      integer, intent(in) :: ip, nndim
      real(rp), intent(in) :: alat
      ! Output
      integer, dimension(nndim + 10), intent(inout) :: iz, n
      real(rp), dimension(3, nndim + 10), intent(inout) :: cr
      ! Local variables
      integer :: i, j

      read (IP, *)
      do I = 1, nndim, 2
         read (IP, *) &
            CR(1, I), CR(2, I), CR(3, I), IZ(I), N(I), CR(1, I + 1), CR(2, I + 1), &
            CR(3, I + 1), IZ(I + 1), N(I + 1)
      end do
      do I = 1, nndim
         do J = 1, 3
            CR(J, I) = ALAT*CR(J, I)
         end do
      end do
      return
   end subroutine

   !> @brief Sort a legacy matrix block by its first column.
   !> @param[in] nl First active row/starting index.
   !> @param[in] ndim Leading dimension of m.
   !> @param[inout] m Matrix block to sort.
   !> @param[in] nd Number of columns.
   !> @param[in] nt Number of active rows.
   module subroutine bubble(nl, ndim, m, nd, nt)
      implicit none
      ! Inputs
      integer, intent(in) :: ndim, nl, nd, nt
      ! Output
      real(rp), dimension(ndim, nd), intent(inout) :: m
      ! Local variables
      integer :: ind, inic, j, k
      real(rp) :: fim, z

      IND = 1
      INIC = 2
      FIM = NL
      do while ((IND == 1) .and. (INIC <= FIM))
         IND = 0
         do J = INT(FIM), INIC, -1
            if (M(J, NT) < M(J - 1, NT)) then
            do K = 1, ND
               Z = M(J, K)
               M(J, K) = M(J - 1, K)
               M(J - 1, K) = Z
            end do
            end if
         end do
         FIM = FIM - 1
         do J = INIC, INT(FIM)
            if (M(J + 1, NT) < M(J, NT)) then
            do K = 1, ND
               Z = M(J + 1, K)
               M(J + 1, K) = M(J, K)
               M(J, K) = Z
            end do
            IND = 1
            end if
         end do
         INIC = INIC + 1
      end do
   end subroutine bubble

   !> @brief Cut translated primitive-cell atoms to a spherical cluster.
   !> @param[in] i Center atom index.
   !> @param[in] l Number of candidate atoms.
   !> @param[in] ndim Array capacity.
   !> @param[in] crd Candidate coordinates.
   !> @param[out] cr Coordinates that survive the cut.
   !> @param[in] izp Candidate atomic labels.
   !> @param[out] iz Atomic labels that survive the cut.
   !> @param[out] num Atom-number labels that survive the cut.
   !> @param[in] no Candidate atom-number labels.
   !> @param[in] rs Cut radius.
   !> @param[out] ii Number of atoms that survived the cut.
   module subroutine cut(i, l, ndim, crd, cr, izp, iz, num, no, rs, ii)
      ! Inputs
      real(rp), intent(in) :: rs
      integer, intent(in) :: i, l, ndim
      real(rp), dimension(3, ndim), intent(in) :: crd
      integer, dimension(ndim), intent(in) :: izp, no
      ! Output
      integer, intent(out) :: ii
      real(rp), dimension(3, ndim), intent(out) :: cr
      integer, dimension(ndim), intent(out) :: iz, num
      ! Local variables
      real(rp) :: r2
      real(rp), dimension(3) :: dum
      integer :: na, j

      do na = 1, l
         r2 = 0.0d0
         do j = 1, 3
            dum(j) = (crd(j, na) - crd(j, i))**2
            r2 = r2 + dum(j)
         end do
         if (r2 .le. rs) then
            ii = ii + 1
            cr(1, ii) = crd(1, i)
            cr(2, ii) = crd(2, i)
            cr(3, ii) = crd(3, i)
            iz(ii) = izp(i)
            num(ii) = no(i)
            return
         end if
      end do
      return
   end subroutine cut

end submodule lattice_cluster
