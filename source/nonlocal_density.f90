!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: nonlocal_density_mod
!
!> @brief Non-local density matrix from the inter-site Chebyshev moments of the SCF recursion.
!> @details Holds the (site, neighbour) pairs requested from the one-shot recursion
!>          when control%nonlocal_dm is set, and their density matrices n_ij.
!------------------------------------------------------------------------------
module nonlocal_density_mod
   use lattice_mod, only: lattice
   use logger_mod, only: g_logger
   use math_mod, only: pi
   use mpi_mod, only: rank, numprocs
   use precision_mod, only: rp
#ifdef USE_MPI
   use mpi
#endif
   implicit none
   private

   type, public :: nonlocal_density
      !> Requested (i, j) pairs of the local recursion sites, self pair first.
      integer, allocatable :: pairs(:, :)
      !> Density matrix n_ij (nb, nb, npairs) of each pair, aligned with pairs.
      complex(rp), allocatable :: n(:, :, :)
   contains
      procedure :: neighbour_pairs
      procedure :: integrate
      procedure :: write
      procedure :: restore_to_default
      final :: destructor
   end type nonlocal_density

contains

   !> For each recursion site i = first..last: (irec(i), irec(i)), then (irec(i), j) for every
   !> neighbour j = nn(irec(i), 2:nn(irec(i),1)) in table order, zeros and repeated entries skipped.
   subroutine neighbour_pairs(this, lat, first, last)
      class(nonlocal_density), intent(inout) :: this
      type(lattice), intent(in) :: lat
      integer, intent(in) :: first, last
      integer :: i, m, k, n, site, j
      integer, allocatable :: tmp(:, :)

      call this%restore_to_default()
      allocate (tmp(max(0, last - first + 1)*size(lat%nn, 2), 2))
      n = 0
      do i = first, last
         site = lat%irec(i)
         n = n + 1
         tmp(n, :) = [site, site]
         k = n
         do m = 2, lat%nn(site, 1)
            j = lat%nn(site, m)
            if (j == 0) cycle
            if (any(tmp(k:n, 2) == j)) cycle
            n = n + 1
            tmp(n, :) = [site, j]
         end do
      end do
      this%pairs = tmp(1:n, :)
   end subroutine neighbour_pairs

   !> n_ij = -1/(2 pi i) int f(E) [G_ij - G_ji^H] dE, with G_ij = sum_m g_m M_ij(m), M_ij(m) = M_ji(m)^H and
   !> M_ji(m) = os_mu(:, :, m, entry) (target j bra, seed i ket, see the gather of chebyshev_oneshot_pairs).
   !> G_ij - G_ji^H = sum_m (g_m - conjg(g_m)) M_ji(m)^H = 2i sum_m Im(g_m) M_ji(m)^H, and c_m = int f(E) g_m(E) dE
   !> (cheb_linear_weights), so  n_ij = -(1/pi) sum_m Im(c_m) M_ji(m)^H.
   !> req must have been built from pairs(1:npairs), so that req%pair_entry(p) is the entry of pair p.
   subroutine integrate(this, req, os_mu, c)
      use intersite_request_mod, only: intersite_request
      class(nonlocal_density), intent(inout) :: this
      type(intersite_request), intent(in) :: req
      complex(rp), intent(in) :: os_mu(:, :, :, :), c(:)
      integer :: p, m

      if (req%npair /= size(this%pairs, 1)) call g_logger%fatal('nonlocal_density: request does not match pairs', &
                                                                __FILE__, __LINE__)
      if (allocated(this%n)) deallocate (this%n)
      allocate (this%n(size(os_mu, 1), size(os_mu, 2), req%npair), source=(0.0_rp, 0.0_rp))
      do p = 1, req%npair
         do m = 1, size(c)
            this%n(:, :, p) = this%n(:, :, p) - aimag(c(m))/pi*conjg(transpose(os_mu(:, :, m, req%pair_entry(p))))
         end do
      end do
   end subroutine integrate

   !> Write nonlocal_dm.out from rank 0 (pairs gathered in rank order = global site order): per pair a line
   !> "site_i atom_j Rx Ry Rz" (R = cr(:, j) - cr(:, i) in alat) and nb lines of nb "Re Im" values of n_ij.
   subroutine write(this, lat)
      class(nonlocal_density), intent(in) :: this
      type(lattice), intent(in) :: lat
      integer, allocatable :: counts(:), displs(:), pg(:, :), pt(:, :)
      complex(rp), allocatable :: ng(:, :, :)
      integer :: nb, np, p, b, r, u, ierr

      nb = size(this%n, 1)
      np = size(this%pairs, 1)
      allocate (counts(numprocs), displs(numprocs))
      counts = np
#ifdef USE_MPI
      call MPI_ALLGATHER(np, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#endif
      displs(1) = 0
      do r = 2, numprocs
         displs(r) = displs(r - 1) + counts(r - 1)
      end do
      allocate (pg(2, sum(counts)), ng(nb, nb, sum(counts)))
      pt = transpose(this%pairs)
#ifdef USE_MPI
      call MPI_GATHERV(pt, 2*np, MPI_INTEGER, pg, 2*counts, 2*displs, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_GATHERV(this%n, nb*nb*np, MPI_DOUBLE_COMPLEX, ng, nb*nb*counts, nb*nb*displs, MPI_DOUBLE_COMPLEX, 0, &
                       MPI_COMM_WORLD, ierr)
#else
      pg = pt
      ng = this%n
#endif
      if (rank /= 0) return
      open (newunit=u, file='nonlocal_dm.out', status='replace', action='write')
      write (u, '(a)') '# site_i atom_j Rx Ry Rz [alat]; then nb rows of nb (Re Im) of n_ij [electrons]'
      do p = 1, size(pg, 2)
         write (u, '(2i6,3es16.8)') pg(1, p), pg(2, p), lat%cr(:, pg(2, p)) - lat%cr(:, pg(1, p))
         do b = 1, nb
            write (u, '(*(2es16.8))') ng(b, :, p)
         end do
      end do
      close (u)
   end subroutine write

   subroutine restore_to_default(this)
      class(nonlocal_density), intent(inout) :: this
      if (allocated(this%pairs)) deallocate (this%pairs)
      if (allocated(this%n)) deallocate (this%n)
   end subroutine restore_to_default

   subroutine destructor(this)
      type(nonlocal_density), intent(inout) :: this
      call this%restore_to_default()
   end subroutine destructor

end module nonlocal_density_mod
