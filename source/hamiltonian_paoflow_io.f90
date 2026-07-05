submodule(hamiltonian_mod) hamiltonian_paoflow_io

contains

   !> @brief Export the real-space Hamiltonian in legacy PAOFLOW layout.
   !> @details Writes hopping records with lattice-cell indices and global orbital
   !>          numbers so external PAOFLOW-style tooling can consume the RS-LMTO
   !>          real-space Hamiltonian.
   !> @param[inout] this Hamiltonian object containing built hopping blocks.
   !> @note This is the legacy export entry point; export_rs_tb_all is the newer dispatcher.
   module subroutine rs2pao(this)
      implicit none
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      real(rp), dimension(3) :: rij, rijtest
      integer :: i, j, k, l, idxi, idxj, idxk, itype, ino, ja, jo, ji, nr, ia, iia, jja, ipao, jpao
      integer :: jj, jt, max_orbital, n_atoms
      integer :: ntype, iostat1, iostat2, iostatus
      real(rp), dimension(3) :: vet, vetpao, idx
      real(rp), dimension(3, 3) :: a_inv
      complex(rp), dimension(nb, nb) :: dum
      n_atoms = this%charge%lattice%ntype
      max_orbital = norb

      open (unit=92, file='rs2paoham.dat', action='write', iostat=iostatus, status='replace')
      do ntype = 1, this%charge%lattice%ntype
         ia = this%charge%lattice%atlist(ntype) ! Atom number in clust
         ino = this%charge%lattice%num(ia) ! Atom bravais type of ia
         nr = this%charge%lattice%nn(ia, 1) ! Number of neighbours considered
         do k = 1, nr
            jj = this%charge%lattice%nn(ia, k)
            !write(123, *)´ia, ii´, ia, m, this%charge%lattice%nn(ia, m)
            if (k == 1) then
               jj = ia
            end if
            if (jj /= 0) then
               rij(:) = this%charge%lattice%cr(:, ia) - this%charge%lattice%cr(:, jj)

               rijtest(:) = 0.0d0
               do idxi = -5, 5
                  do idxj = -5, 5
                     do idxk = -5, 5
                        rijtest(:) = this%charge%lattice%cr(:, ia) - (this%charge%lattice%cr(:, this%charge%lattice%iz(jj)) &
                                                                      + idxi*this%charge%lattice%a(:, 1) + idxj*this%charge%lattice%a(:, 2) + idxk*this%charge%lattice%a(:, 3))
                        if (norm2(rij(:) - rijtest(:)) < 1.0d-3) then
                           idx(:) = [idxi, idxj, idxk]
                        end if
                     end do
                  end do
               end do

               if (k == 1) this%ee(:, :, k, ntype) = this%ee(:, :, k, ntype) + this%lsham(:, :, ntype) !+ this%enim(:,:,ntype)

               call hcpx(this%ee(1:norb,1:norb,k,ntype), 'sph2cart')
               call hcpx(this%ee(1:norb,norb+1:nb,k,ntype), 'sph2cart')
               call hcpx(this%ee(norb+1:nb,1:norb,k,ntype), 'sph2cart')
               call hcpx(this%ee(norb+1:nb,norb+1:nb,k,ntype), 'sph2cart')
               do i = 1, nb
                  do j = 1, nb
                     ipao = 0; jpao = 0
                     call site2orb(i, ia, ipao, n_atoms, max_orbital)
                     call site2orb(j, this%charge%lattice%iz(jj), jpao, n_atoms, max_orbital)
                     write (92, '(3I4,2I7,2F22.14)') int(idx(:)), ipao, jpao, real(this%ee(i, j, k, ntype))*ry2ev, aimag(this%ee(i, j, k, ntype))*ry2ev
                  end do
               end do
            end if
         end do
      end do
      close (92)
   end subroutine rs2pao

   !> @brief Import a PAOFLOW-format real-space Hamiltonian from paoham.dat.
   !> @details Reads legacy seven-column hopping records and fills the bulk
   !>          Hamiltonian blocks for PAOFLOW-to-real-space post-processing routes.
   !> @param[inout] this Hamiltonian object; replaces bulk hopping arrays from file data.
   !> @note PAOFLOW import is used by paoflow2rs, exchange_p2rs, and conductivity_p2rs.
   module subroutine build_from_paoflow_opt(this)
      implicit none
      type hamData
         integer :: idxi, idxj, idxk
         integer :: orbl, orbm
         real :: dumre, dumcmplx
      end type hamData
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m, n, itype, ino, ja, jo, ji, nr, ia, iia, jja
      integer :: jj, jt, orbl, orbm, idxi, idxj, idxk
      integer :: ntype, iostat1, iostat2, iostatus, n_atoms, max_orbital, numLines
      real(rp), dimension(3) :: vet, vetpao
      real(rp) :: dumre, dumcmplx
      integer, dimension(maxval(this%charge%lattice%nn(:, 1)) + 1, 3) :: idxup, idxdw, idx
      type(hamData) :: ham
      type(hamData), allocatable :: hamArray(:)

      n_atoms = this%charge%lattice%ntype
      max_orbital = norb
      numLines = countLines('paoham.dat')
      allocate (hamArray(numLines))

      ! Reading Hamiltonian data once
      open (unit=92, file='paoham.dat', action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) then
         call g_logger%fatal('file containing the paoflow Hamiltonian not found', __FILE__, __LINE__)
      end if
      do i = 1, numLines
         read (92, *, iostat=iostat1) hamArray(i)%idxi, hamArray(i)%idxj, hamArray(i)%idxk, &
            hamArray(i)%orbl, hamArray(i)%orbm, hamArray(i)%dumre, hamArray(i)%dumcmplx
      end do
      close (92)
      idx = 0
      write (*, *) 'PAOFLOW Hamiltonian has been read'
      call g_timer%start('Hamiltonian allocation')
      !$omp parallel do private(ntype, ia, ino, nr, jj, k, l, ham, vet, vetpao, i, j, iia, jja, idx) shared(this, numLines, hamArray, n_atoms, max_orbital)
      do ntype = 1, this%charge%lattice%ntype
         ia = this%charge%lattice%atlist(ntype)
         ino = this%charge%lattice%num(ia)
         nr = this%charge%lattice%nn(ia, 1)
         write (*, *) ia, nr, n_atoms, max_orbital
         do k = 1, nr
            jj = this%charge%lattice%nn(ia, k)
            if (k == 1) then
               jj = ia
            end if
            if (jj /= 0) then
               vet(:) = (this%charge%lattice%cr(:, ia) - this%charge%lattice%cr(:, jj))!*this%charge%lattice%alat
               do l = 1, numLines
                  ham = hamArray(l)
                  call orb2site(ham%orbl, i, iia, n_atoms, max_orbital)
                  call orb2site(ham%orbm, j, jja, n_atoms, max_orbital)
                  if (iia == ia) then
                     vetpao(:) = this%charge%lattice%cr(:, iia) - (this%charge%lattice%cr(:, jja) &
                                                                   + ham%idxi*this%charge%lattice%a(:, 1) + ham%idxj*this%charge%lattice%a(:, 2) + ham%idxk*this%charge%lattice%a(:, 3))!*this%charge%lattice%alat
                     if (norm2(vet(:) - vetpao(:)) < 1.0d-3) then
                        idx(k, :) = [ham%idxi, ham%idxj, ham%idxk]
                        this%ee(i, j, k, ntype) = cmplx(ham%dumre, ham%dumcmplx)/13.605703976
                        ! if(ntype==1.or.ntype==2)then
                        !   if(i==j.and.k==1.and.i>=5.and.i<=9)then
                        !     this%ee(i, j, k, ntype) = cmplx(ham%dumre-1.113 , ham%dumcmplx)/13.605703976
                        !   else if(i==j.and.k==1.and.i>=14.and.i<=18)then
                        !     this%ee(i, j, k, ntype) = cmplx(ham%dumre+1.113 , ham%dumcmplx)/13.605703976
                        !   end if
                        ! end if
                     end if
                  end if
               end do
            end if
            write (128, *) 'm=', k, 'Atom=', jj, 'Coordinates=', this%charge%lattice%cr(:, jj), 'Ntype=', ntype, 'Index=', idx(k, :)
            write (129, *) 'm=', k, 'Atom=', jj, 'Coordinates=', this%charge%lattice%cr(:, jj), 'Ntype=', ntype, 'Index=', idx(k, :)
            write (128, '(18f10.6)') real(this%EE(1:nb, 1:nb, k, ntype))*13.605703976
            write (129, '(18f10.6)') aimag(this%EE(1:nb, 1:nb, k, ntype))*13.605703976
            write (128, *) sum(real(this%ee(:, :, k, ntype)))
            write (129, *) sum(real(this%ee(:, :, k, ntype)))
         end do
      end do
      !$omp end parallel do
      deallocate (hamArray)
      call g_timer%stop('Hamiltonian allocation')
   end subroutine build_from_paoflow_opt

   !> @brief Import a PAOFLOW-format real-space Hamiltonian with the legacy reader.
   !> @details Maps PAOFLOW orbital/cell records back onto RS-LMTO neighbor blocks.
   !>          Kept for compatibility with older import paths.
   !> @param[inout] this Hamiltonian object; fills bulk hopping arrays from file data.
   module subroutine build_from_paoflow(this)
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m, n, itype, ino, ja, jo, ji, nr, ia, iia, jja
      integer :: jj, jt, orbl, orbm, idxi, idxj, idxk
      integer :: ntype, iostat1, iostat2, iostatus, n_atoms, max_orbital
      real(rp), dimension(3) :: vet, vetpao, cri_dir, crj_dir, cri_cart, crj_cart
      integer, dimension(maxval(this%charge%lattice%nn(:, 1)) + 1, 3) :: idxup, idxdw, idx
      real(rp) :: dumre, dumcmplx

      n_atoms = this%charge%lattice%ntype
      max_orbital = norb
      open (unit=90, file='paoup.dat', action='read', iostat=iostatus, status='old')
      open (unit=91, file='paodw.dat', action='read', iostat=iostatus, status='old')
      open (unit=92, file='paoham.dat', action='read', iostat=iostatus, status='old')

      if (iostatus /= 0) then
         call g_logger%fatal('file containing the paoflow Hamiltonian not found', __FILE__, __LINE__)
      end if
      do ntype = 1, this%charge%lattice%ntype
         ia = this%charge%lattice%atlist(ntype) ! Atom number in clust
         ino = this%charge%lattice%num(ia) ! Atom bravais type of ia
         nr = this%charge%lattice%nn(ia, 1) ! Number of neighbours considered
         write (*, *) ia, nr, n_atoms, max_orbital
         do k = 1, nr
            jj = this%charge%lattice%nn(ia, k)
            !write(123, *)´ia, ii´, ia, m, this%charge%lattice%nn(ia, m)
            if (k == 1) then
               jj = ia
            end if
            if (jj /= 0) then
               !cri_cart(:) = this%charge%lattice%cr(:, ia)
               !crj_cart(:) = this%charge%lattice%cr(:, jj)

               !cri_dir(:) = cartesian_to_direct(this%charge%lattice%a,cri_cart)
               !crj_dir(:) = cartesian_to_direct(this%charge%lattice%a,crj_cart)

               !vet(:) = cri_dir(:) - crj_dir(:)
               vet(:) = (this%charge%lattice%cr(:, ia) - this%charge%lattice%cr(:, jj))!*this%charge%lattice%alat
               do
                  read (92, *, iostat=iostat1) idxi, idxj, idxk, orbl, orbm, dumre, dumcmplx
                  if (iostat1 /= 0) then
                     exit
                  end if
                  if (orbl <= n_atoms*max_orbital) then
                     i = modulo(orbl - 1, max_orbital) + 1
                     iia = int((orbl - 1)/max_orbital) + 1
                  else
                     i = modulo(orbl - 1, max_orbital) + 10
                     iia = int((orbl - 1 - n_atoms*max_orbital)/max_orbital) + 1
                  end if

                  if (orbm <= n_atoms*max_orbital) then
                     j = modulo(orbm - 1, max_orbital) + 1
                     jja = int((orbm - 1)/max_orbital) + 1
                  else
                     j = modulo(orbm - 1, max_orbital) + 10
                     jja = int((orbm - 1 - n_atoms*max_orbital)/max_orbital) + 1
                  end if

                  !cri_cart(:) = this%charge%lattice%cr(:, iia)
                  !crj_cart(:) = this%charge%lattice%cr(:, jja)

                  !cri_dir(:) = cartesian_to_direct(this%charge%lattice%a,cri_cart)
                  !crj_dir(:) = cartesian_to_direct(this%charge%lattice%a,crj_cart)
                  if (iia == ia) then
                     !vetpao(:) = cri_dir(:) - (crj_dir(:) + [idxi,idxj,idxk])

                     vetpao(:) = this%charge%lattice%cr(:, iia) - (this%charge%lattice%cr(:, jja) &
                                                                   + idxi*this%charge%lattice%a(:, 1) + idxj*this%charge%lattice%a(:, 2) + idxk*this%charge%lattice%a(:, 3))!*this%charge%lattice%alat
                     if (norm2(vet(:) - vetpao(:)) < 1.0d-3) then
                        idx(k, :) = [idxi, idxj, idxk]
                        this%ee(i, j, k, ntype) = cmplx(dumre, dumcmplx)/13.605703976
                     end if
                  end if
               end do
            end if
            write (128, *) 'm=', k, 'Atom=', jj, 'Coordinates=', this%charge%lattice%cr(:, jj), 'Ntype=', ntype, 'Index=', idx(k, :)
            write (128, '(18f10.6)') real(this%EE(1:nb, 1:nb, k, ntype))!*13.605703976
            write (128, *) sum(real(this%ee(:, :, k, ntype)))
            rewind (90)
            rewind (91)
            rewind (92)
         end do
      end do
   end subroutine build_from_paoflow

   !> @brief Export real-space tight-binding data in all selected formats.
   !> @details Dispatcher for metadata and hopping-record writers used by the
   !>          hamiltonian export namelist option. Supports PAOFLOW legacy and
   !>          Python-friendly real-space records.
   !> @param[in] this Hamiltonian object containing built hopping blocks.
   !> @param[in] basename Optional output basename.
   !> @param[in] tol Optional threshold for skipping tiny records.
   !> @param[in] include_lsham Optional flag to include onsite SOC blocks.
   !> @param[in] transform_sph2cart Optional flag to transform orbital basis.
   module subroutine export_rs_tb_all(this, basename, tol, include_lsham, transform_sph2cart)
      implicit none
      class(hamiltonian), intent(in) :: this
      character(len=*), intent(in), optional :: basename
      real(rp), intent(in), optional :: tol
      logical, intent(in), optional :: include_lsham, transform_sph2cart
   
      character(len=512) :: base
      real(rp) :: eps
      logical :: add_lsham, do_sph2cart
   
      base = 'rs_tb'
      if (present(basename)) base = trim(basename)
   
      eps = 1.0e-3_rp
      if (present(tol)) eps = tol
   
      add_lsham = .true.
      if (present(include_lsham)) add_lsham = include_lsham
   
      do_sph2cart = .true.
      if (present(transform_sph2cart)) do_sph2cart = transform_sph2cart
   
      call export_rs_tb_metadata(this, trim(base)//'.meta')
      call export_rs_tb_hoppings(this, trim(base)//'.tb', eps, add_lsham, do_sph2cart)
      call export_rs_paoflow_legacy(this, trim(base)//'_paoham.dat', eps, add_lsham, do_sph2cart)
   end subroutine export_rs_tb_all

   !> @brief Write metadata for real-space tight-binding exports.
   !> @details Records atom/orbital layout information needed to interpret the
   !>          exported hopping records.
   !> @param[in] this Hamiltonian object containing lattice and basis metadata.
   !> @param[in] filename Metadata output path.
   module subroutine export_rs_tb_metadata(this, filename)
      implicit none
      class(hamiltonian), intent(in) :: this
      character(len=*), intent(in) :: filename
   
      integer :: u, i, ios, n_atoms, max_orbital
   
      n_atoms = this%charge%lattice%ntype
      max_orbital = norb
   
      open(newunit=u, file=filename, action='write', status='replace', iostat=ios)
      if (ios /= 0) then
         call g_logger%fatal('could not open TB metadata file for writing', __FILE__, __LINE__)
      end if
   
      write(u,'(A)') '# rs-tb-meta-v1'
      write(u,'(A)') 'energy_unit eV'
      write(u,'(A)') 'length_unit Angstrom'
      write(u,'(A,I0)') 'n_atoms ', n_atoms
      write(u,'(A,I0)') 'norb_scalar ', max_orbital
      write(u,'(A,I0)') 'block_size ', nb
      ! write(u,'(A,f12.6)') 'Fermi level', this
      write(u,'(A)') 'index_base 1'
      write(u,'(A)') 'lattice_vectors_cart'
      do i = 1, 3
         write(u,'(3ES26.16)') this%charge%lattice%a(:, i)
      end do
      write(u,'(A)') 'positions_cart'
      do i = 1, n_atoms
         ! atlist(i) maps type/basis index to the representative cluster atom.
         write(u,'(I8,3ES26.16)') i, this%charge%lattice%cr(:, this%charge%lattice%atlist(i))
      end do
      close(u)
   end subroutine export_rs_tb_metadata

   !> @brief Write PAOFLOW legacy-format hopping records.
   !> @details Emits the seven-column record layout read by build_from_paoflow_opt:
   !>          cell index, global orbital indices, and complex hopping value.
   !> @param[in] this Hamiltonian object containing built hopping blocks.
   !> @param[in] filename Output path.
   !> @param[in] tol Threshold for skipping tiny records.
   !> @param[in] include_lsham Include onsite SOC blocks in the output.
   !> @param[in] transform_sph2cart Transform orbital blocks before writing.
   module subroutine export_rs_paoflow_legacy(this, filename, tol, include_lsham, transform_sph2cart)
      implicit none
      class(hamiltonian), intent(in) :: this
      character(len=*), intent(in) :: filename
      real(rp), intent(in) :: tol
      logical, intent(in) :: include_lsham, transform_sph2cart
   
      integer :: u, ios
   
      open(newunit=u, file=filename, action='write', status='replace', iostat=ios)
      if (ios /= 0) then
         call g_logger%fatal('could not open legacy PAOFLOW-like Hamiltonian file', __FILE__, __LINE__)
      end if
   
      call write_rs_tb_records(this, u, 'legacy7', tol, include_lsham, transform_sph2cart)
      close(u)
   end subroutine export_rs_paoflow_legacy

   !> @brief Write Python-friendly real-space hopping records.
   !> @details Emits hopping records plus metadata conventions intended for direct
   !>          parsing by scripts and post-processing tools.
   !> @param[in] this Hamiltonian object containing built hopping blocks.
   !> @param[in] filename Output path.
   !> @param[in] tol Threshold for skipping tiny records.
   !> @param[in] include_lsham Include onsite SOC blocks in the output.
   !> @param[in] transform_sph2cart Transform orbital blocks before writing.
   module subroutine export_rs_tb_hoppings(this, filename, tol, include_lsham, transform_sph2cart)
      implicit none
      class(hamiltonian), intent(in) :: this
      character(len=*), intent(in) :: filename
      real(rp), intent(in) :: tol
      logical, intent(in) :: include_lsham, transform_sph2cart
   
      integer :: u, ios
   
      open(newunit=u, file=filename, action='write', status='replace', iostat=ios)
      if (ios /= 0) then
         call g_logger%fatal('could not open canonical TB hopping file', __FILE__, __LINE__)
      end if
   
      write(u,'(A)') '# rs-tb-hoppings-v1'
      write(u,'(A)') '# columns:'
      write(u,'(A)') '# R1 R2 R3 ia jbasis k_neigh i_local j_local ipao jpao real_eV imag_eV'
      call write_rs_tb_records(this, u, 'canonical13', tol, include_lsham, transform_sph2cart)
      close(u)
   end subroutine export_rs_tb_hoppings

   !> @brief Write real-space hopping records to an open unit.
   !> @details Shared implementation for legacy PAOFLOW and Python export modes,
   !>          walking local/bulk neighbor blocks and optionally transforming basis
   !>          or adding onsite spin-orbit terms.
   !> @param[in] this Hamiltonian object containing built hopping blocks.
   !> @param[in] u Open output unit.
   !> @param[in] mode Record format selector.
   !> @param[in] tol Threshold for skipping tiny records.
   !> @param[in] include_lsham Include onsite SOC blocks in the output.
   !> @param[in] transform_sph2cart Transform orbital blocks before writing.
   module subroutine write_rs_tb_records(this, u, mode, tol, include_lsham, transform_sph2cart)
      implicit none
      class(hamiltonian), intent(in) :: this
      integer, intent(in) :: u
      character(len=*), intent(in) :: mode
      real(rp), intent(in) :: tol
      logical, intent(in) :: include_lsham, transform_sph2cart
   
      integer :: ntype, ia, nr, k, jj, jbasis
      integer :: i, j, ipao, jpao, n_atoms, max_orbital
      integer :: idx(3)
      logical :: found
      complex(rp) :: hblock(nb, nb)
   
      n_atoms = this%charge%lattice%ntype
      max_orbital = norb
   
      do ntype = 1, this%charge%lattice%ntype
         ia = this%charge%lattice%atlist(ntype)
         nr = this%charge%lattice%nn(ia, 1)
      
         do k = 1, nr
            jj = this%charge%lattice%nn(ia, k)
            if (k == 1) jj = ia
            if (jj == 0) cycle
         
            call rs_neighbor_lattice_index(this, ia, jj, idx, found, tol)
            if (.not. found) then
               write(*,'(A,2I8,A)') 'WARNING: no lattice index found for neighbour pair ia,jj=', ia, jj, '; skipping'
               cycle
            end if
         
            jbasis = this%charge%lattice%iz(jj)
         
            hblock(:, :) = this%ee(:, :, k, ntype)
            if (include_lsham .and. k == 1) hblock(:, :) = hblock(:, :) + this%lsham(:, :, ntype)
         
            if (transform_sph2cart) then
               call hcpx(hblock(1:norb,       1:norb      ), 'sph2cart')
               call hcpx(hblock(1:norb,       norb+1:nb   ), 'sph2cart')
               call hcpx(hblock(norb+1:nb,    1:norb      ), 'sph2cart')
               call hcpx(hblock(norb+1:nb,    norb+1:nb   ), 'sph2cart')
            end if
         
            do i = 1, nb
               do j = 1, nb
                  if (abs(hblock(i,j)) <= 0.0_rp) cycle
               
                  ipao = 0
                  jpao = 0
                  call site2orb(i, ia,     ipao, n_atoms, max_orbital)
                  call site2orb(j, jbasis, jpao, n_atoms, max_orbital)
               
                  select case (trim(mode))
                  case ('legacy7')
                     ! PAOFLOW-like/local legacy format used by build_from_paoflow_opt:
                     ! R1 R2 R3 global_orb_i global_orb_j Re[eV] Im[eV]
                     write(u,'(3I6,2I10,2ES26.16)') idx(:), ipao, jpao, &
                        real(hblock(i,j))*ry2ev, aimag(hblock(i,j))*ry2ev
                  case ('canonical13')
                     ! Extended self-describing scalar record. ia and jbasis are 1-based
                     ! host-code atom/basis identifiers; i,j and ipao,jpao are 1-based.
                     write(u,'(3I6,7I10,2ES26.16)') idx(:), ia, jbasis, k, i, j, ipao, jpao, &
                        real(hblock(i,j))*ry2ev, aimag(hblock(i,j))*ry2ev
                  end select
               end do
            end do
         end do
      end do
   end subroutine write_rs_tb_records

   !> @brief Resolve a neighbor into integer lattice-cell indices.
   !> @details Matches a real-space neighbor displacement against lattice
   !>          translations so exported hopping records can carry PAOFLOW-style
   !>          integer cell offsets.
   !> @param[in] this Hamiltonian object containing lattice vectors and coordinates.
   !> @param[in] ia Source atom index.
   !> @param[in] jj Neighbor-list entry.
   !> @param[out] idx Integer lattice-cell offset.
   !> @param[out] found True when a matching offset was found within tolerance.
   !> @param[in] tol Matching tolerance.
   module subroutine rs_neighbor_lattice_index(this, ia, jj, idx, found, tol)
      implicit none
      class(hamiltonian), intent(in) :: this
      integer, intent(in) :: ia, jj
      integer, intent(out) :: idx(3)
      logical, intent(out) :: found
      real(rp), intent(in) :: tol
   
      integer :: idxi, idxj, idxk, jbasis
      real(rp) :: rij(3), rijtest(3), err, best_err
   
      idx = 0
      found = .false.
      best_err = huge(1.0_rp)
   
      ! For onsite blocks, do not search a supercell image.
      if (ia == jj) then
         idx = [0, 0, 0]
         found = .true.
         return
      end if
   
      rij(:) = this%charge%lattice%cr(:, ia) - this%charge%lattice%cr(:, jj)
      jbasis = this%charge%lattice%iz(jj)
   
      do idxi = -5, 5
         do idxj = -5, 5
            do idxk = -5, 5
               rijtest(:) = this%charge%lattice%cr(:, ia) - &
                  (this%charge%lattice%cr(:, jbasis) + &
                   idxi*this%charge%lattice%a(:, 1) + &
                   idxj*this%charge%lattice%a(:, 2) + &
                   idxk*this%charge%lattice%a(:, 3))
               err = norm2(rij(:) - rijtest(:))
               if (err < best_err) then
                  best_err = err
                  idx = [idxi, idxj, idxk]
               end if
            end do
         end do
      end do
   
      found = (best_err < tol)
   end subroutine rs_neighbor_lattice_index

end submodule hamiltonian_paoflow_io
