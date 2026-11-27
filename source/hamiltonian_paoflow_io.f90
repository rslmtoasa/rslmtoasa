submodule(hamiltonian_mod) hamiltonian_paoflow_io

contains

   module subroutine rs2pao(this)
      class(hamiltonian), intent(inout) :: this
      call this%export_rs_tb_all()
   end subroutine rs2pao

   !> @brief Import a PAOFLOW-format real-space Hamiltonian from paoham.dat.
   !> @details Reads legacy seven-column hopping records and fills the bulk
   !>          Hamiltonian blocks for PAOFLOW-to-real-space post-processing routes.
   !> @param[inout] this Hamiltonian object; replaces bulk hopping arrays from file data.
   !> @note PAOFLOW import is used by paoflow2rs, exchange_p2rs, and conductivity_p2rs.
   module subroutine build_from_paoflow_opt(this)
      implicit none
      ! The file is copied into ee unchanged, i.e. it must be in the spherical basis (tests/KNOWN_ISSUES.md, Orbital basis audit)
      type hamData
         integer :: idxi, idxj, idxk
         integer :: orbl, orbm
         real(rp) :: dumre, dumcmplx
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
                        this%ee(i, j, k, ntype) = cmplx(ham%dumre, ham%dumcmplx, rp)/ry2ev
                        ! if(ntype==1.or.ntype==2)then
                        !   if(i==j.and.k==1.and.i>=5.and.i<=9)then
                        !     this%ee(i, j, k, ntype) = cmplx(ham%dumre-1.113 , ham%dumcmplx)/ry2ev
                        !   else if(i==j.and.k==1.and.i>=14.and.i<=18)then
                        !     this%ee(i, j, k, ntype) = cmplx(ham%dumre+1.113 , ham%dumcmplx)/ry2ev
                        !   end if
                        ! end if
                     end if
                  end if
               end do
            end if
            write (128, *) 'm=', k, 'Atom=', jj, 'Coordinates=', this%charge%lattice%cr(:, jj), 'Ntype=', ntype, 'Index=', idx(k, :)
            write (129, *) 'm=', k, 'Atom=', jj, 'Coordinates=', this%charge%lattice%cr(:, jj), 'Ntype=', ntype, 'Index=', idx(k, :)
            write (128, '(18f10.6)') real(this%EE(1:nb, 1:nb, k, ntype))*ry2ev
            write (129, '(18f10.6)') aimag(this%EE(1:nb, 1:nb, k, ntype))*ry2ev
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
                        this%ee(i, j, k, ntype) = cmplx(dumre, dumcmplx, rp)/ry2ev
                     end if
                  end if
               end do
            end if
            write (128, *) 'm=', k, 'Atom=', jj, 'Coordinates=', this%charge%lattice%cr(:, jj), 'Ntype=', ntype, 'Index=', idx(k, :)
            write (128, '(18f10.6)') real(this%EE(1:nb, 1:nb, k, ntype))!*ry2ev
            write (128, *) sum(real(this%ee(:, :, k, ntype)))
            rewind (90)
            rewind (91)
            rewind (92)
         end do
      end do
   end subroutine build_from_paoflow

end submodule hamiltonian_paoflow_io
