!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: self_energy_mod
!
!> @brief One-shot local Green function G_ii(z) = <i|(z - H - Sigma_t(z))^-1|i> for a site-local,
!> non-Hermitian self-energy Sigma_t(z) per atom type t, read from sigma.dat.
!> @details Chebyshev expansion of K = H + Sigma in Kt = (K - b)/a with complex centre b (Re b folded into the
!> Hamiltonian shift, i Im b folded into the on-site blocks). No moment doubling: order n costs n applications.
!> All energies in Ry.
!------------------------------------------------------------------------------
module self_energy_mod
   use precision_mod, only: rp
   use basis_mod, only: nb
   use logger_mod, only: g_logger
   use string_mod, only: int2str, real2str
   use recursion_mod, only: recursion
   use cheb_sigma_kernel_mod, only: cheb_sigma_window, cheb_sigma_w, cheb_sigma_sum
   use mpi_mod, only: rank, start_atom, end_atom
#ifdef USE_MPI
   use mpi
#endif
   implicit none
   private

   type, public :: self_energy
      !> Relative tail bound that marks a Green function as converged
      real(rp) :: tol
      !> Number of energies, atom types and block size of sigma.dat
      integer :: nz, ntype, nb
      !> Energy offsets x + i y relative to E_F (Ry)
      complex(rp), allocatable :: zoff(:)
      !> Self-energy blocks (nb,nb,ntype,nz) (Ry)
      complex(rp), allocatable :: sig(:, :, :, :)
      !> Local Green function blocks (nb,nb,nrec,nz) (1/Ry)
      complex(rp), allocatable :: gloc(:, :, :, :)
      !> Tail bound satisfied (nrec,nz)
      logical, allocatable :: converged(:, :)
      !> Number of Chebyshev moments used (nrec,nz); 0 if skipped
      integer, allocatable :: order(:, :)
   contains
      procedure :: read => read_sigma
      procedure :: compute_gloc
      procedure :: write => write_gloc
      procedure :: restore_to_default
   end type self_energy

   interface self_energy
      procedure :: constructor
   end interface self_energy

contains

   !> Constructor: default state.
   function constructor() result(obj)
      type(self_energy) :: obj
      call obj%restore_to_default()
   end function constructor

   !> Reset all members to default.
   subroutine restore_to_default(this)
      class(self_energy), intent(inout) :: this
      this%tol = 1.0e-6_rp
      this%nz = 0
      this%ntype = 0
      this%nb = 0
      if (allocated(this%zoff)) deallocate (this%zoff)
      if (allocated(this%sig)) deallocate (this%sig)
      if (allocated(this%gloc)) deallocate (this%gloc)
      if (allocated(this%converged)) deallocate (this%converged)
      if (allocated(this%order)) deallocate (this%order)
   end subroutine restore_to_default

   !> Read sigma.dat: header line, "nz ntype nb", then per z "x y" and ntype blocks of nb rows of nb (re im) pairs.
   !> The block size must match the run.
   subroutine read_sigma(this, fname)
      class(self_energy), intent(inout) :: this
      character(len=*), intent(in) :: fname
      real(rp) :: x, y, row(2*nb)
      integer :: unit, ios, iz, t, l

      open (newunit=unit, file=fname, action='read', status='old', iostat=ios)
      if (ios /= 0) call g_logger%fatal('self-energy file '//trim(fname)//' not found', __FILE__, __LINE__)
      read (unit, *)
      read (unit, *, iostat=ios) this%nz, this%ntype, this%nb
      if (ios /= 0) call g_logger%fatal('cannot read nz ntype nb from '//trim(fname), __FILE__, __LINE__)
      if (this%nb /= nb) &
         call g_logger%fatal(trim(fname)//': nb = '//int2str(this%nb)//' but the run has '//int2str(nb), __FILE__, __LINE__)
      allocate (this%zoff(this%nz), this%sig(nb, nb, this%ntype, this%nz))
      do iz = 1, this%nz
         read (unit, *, iostat=ios) x, y
         if (ios /= 0) call g_logger%fatal('cannot read energy '//int2str(iz)//' from '//trim(fname), __FILE__, __LINE__)
         this%zoff(iz) = cmplx(x, y, rp)
         do t = 1, this%ntype
            do l = 1, nb
               read (unit, *, iostat=ios) row
               if (ios /= 0) call g_logger%fatal('cannot read Sigma block of energy '//int2str(iz)//' from '//trim(fname), &
                                                 __FILE__, __LINE__)
               this%sig(l, :, t, iz) = cmplx(row(1::2), row(2::2), rp)
            end do
         end do
      end do
      close (unit)
   end subroutine read_sigma

   !> G_ii(z = fermi + zoff) of every recursion atom; moments by the Chebyshev recurrence in Kt, log-scaled.
   !> With gpu_plugin the moments of all z of one atom come from one batched device run.
   subroutine compute_gloc(this, rec, fermi)
      class(self_energy), intent(inout) :: this
      class(recursion), intent(inout) :: rec
      real(rp), intent(in) :: fermi
      real(rp) :: emin, emax, rho_s, rho_z, tail, m, lscale, x
      real(rp), allocatable :: logscale(:), az(:), rz(:), ag(:), bg(:), lsg(:, :)
      complex(rp) :: b, z
      complex(rp), allocatable :: bz(:), wz(:), dz(:, :, :, :), dg(:, :, :, :), mu(:, :, :), mug(:, :, :, :), &
                                  prev(:, :, :), cur(:, :, :), nxt(:, :, :), tmp(:, :, :), work(:, :, :, :)
      integer, allocatable :: nmz(:), idx(:), ng(:)
      logical :: on_gpu
#ifdef USE_MPI
      integer, allocatable :: nconv(:, :)
      integer :: ierr
#endif
      integer :: i, j, iz, t, k, n, g, kk, lld, nrec

      if (this%ntype /= rec%lattice%ntype) &
         call g_logger%fatal('sigma.dat has '//int2str(this%ntype)//' atom types but the run has '// &
                             int2str(rec%lattice%ntype), __FILE__, __LINE__)
      nrec = rec%lattice%nrec
      kk = rec%lattice%kk
      lld = rec%lattice%control%lld
      call rec%resolve_chebyshev_window(emin, emax)
      if (allocated(this%gloc)) deallocate (this%gloc, this%converged, this%order)
      allocate (this%gloc(nb, nb, nrec, this%nz), source=(0.0_rp, 0.0_rp))
      allocate (this%converged(nrec, this%nz), source=.false.)
      allocate (this%order(nrec, this%nz), source=0)
      allocate (az(this%nz), rz(this%nz), bz(this%nz), wz(this%nz), nmz(this%nz), dz(nb, nb, this%ntype, this%nz))

      on_gpu = rec%gpu_plugin_ready('self_energy%compute_gloc()', allow_hoh=.true.)
      if (on_gpu) then
         call rec%gpu_plugin_upload_hamiltonian()
         call g_logger%info('self-energy on GPU: '//real2str(real(merge(7, 4, rec%hamiltonian%hoh)*merge(16, 8, &
                            rec%gpu_precision_code() == 1), rp)*real(nb, rp)**2*real(kk, rp)/1.0e6_rp, '(f10.1)')// &
                            ' MB device memory per energy', __FILE__, __LINE__)
      else
         allocate (prev(nb, nb, kk), cur(nb, nb, kk), nxt(nb, nb, kk), tmp(nb, nb, kk), work(nb, nb, kk, 3), logscale(lld))
      end if

      do i = start_atom, end_atom
         j = rec%lattice%irec(i)
         nmz(:) = 0
         do iz = 1, this%nz
            z = cmplx(fermi, 0.0_rp, rp) + this%zoff(iz)
            call cheb_sigma_window(emin, emax, this%sig(:, :, :, iz), az(iz), b, rho_s)
            bz(iz) = b
            call cheb_sigma_w(z, az(iz), b, wz(iz), rho_z)
            rz(iz) = rho_s/rho_z
            if (rz(iz) >= 1.0_rp) then
               call g_logger%warning('self-energy: z offset '//real2str(real(this%zoff(iz), rp), '(f9.4)')//' '// &
                                     real2str(aimag(this%zoff(iz)), '(f9.4)')//' Ry is inside the ellipse of H + Sigma (atom '// &
                                     int2str(j)//'); not computed', __FILE__, __LINE__)
               cycle
            end if
            x = log(this%tol*(1.0_rp - rz(iz))/(2.0_rp*rz(iz)))/log(rz(iz))
            nmz(iz) = min(lld, max(2, ceiling(min(real(lld, rp), x)) + 1))

            ! On-site part of Kt: (Sigma_t - i Im b)/a; the Hamiltonian part is (H - Re b)/a.
            do t = 1, this%ntype
               dz(:, :, t, iz) = this%sig(:, :, t, iz)
               do k = 1, nb
                  dz(k, k, t, iz) = dz(k, k, t, iz) - cmplx(0.0_rp, aimag(b), rp)
               end do
            end do
            dz(:, :, :, iz) = dz(:, :, :, iz)/az(iz)
         end do

         if (on_gpu) then
            idx = pack([(iz, iz=1, this%nz)], nmz > 0)
            if (size(idx) == 0) cycle
            ag = az(idx)
            bg = real(bz(idx), rp)
            dg = dz(:, :, :, idx)
            ng = nmz(idx)
            allocate (mug(nb, nb, maxval(nmz), size(idx)), lsg(maxval(nmz), size(idx)))
            call rec%gpu_backend%chebyshev_sigma_batch(j, ng, ag, bg, dg, mug, lsg)
            do g = 1, size(idx)
               iz = idx(g)
               call cheb_sigma_sum(mug(:, :, 1:nmz(iz), g), lsg(1:nmz(iz), g), wz(iz), az(iz), rz(iz), &
                                   this%gloc(:, :, i, iz), tail)
               call record(iz, tail)
            end do
            deallocate (mug, lsg)
            cycle
         end if

         do iz = 1, this%nz
            if (nmz(iz) == 0) cycle
            allocate (mu(nb, nb, nmz(iz)))
            prev = (0.0_rp, 0.0_rp)
            do k = 1, nb
               prev(k, k, j) = (1.0_rp, 0.0_rp)
            end do
            rec%izero(:) = 0
            rec%izero(j) = 1
            mu(:, :, 1) = prev(:, :, j)
            logscale(1) = 0.0_rp
            lscale = 0.0_rp
            do n = 2, nmz(iz)
               call apply_kt(n == 2, prev, cur, nxt)
               mu(:, :, n) = cur(:, :, j)
               logscale(n) = lscale
            end do
            call cheb_sigma_sum(mu, logscale(1:nmz(iz)), wz(iz), az(iz), rz(iz), this%gloc(:, :, i, iz), tail)
            call record(iz, tail)
            deallocate (mu)
         end do
      end do

#ifdef USE_MPI
      allocate (nconv(nrec, this%nz))
      nconv = merge(1, 0, this%converged)
      call MPI_ALLREDUCE(MPI_IN_PLACE, this%gloc, size(this%gloc), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, nconv, size(nconv), MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, this%order, size(this%order), MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      this%converged = nconv > 0
#endif

   contains

      !> Store order and convergence of (atom i, energy iz) and log them.
      subroutine record(iz_, tail_)
         integer, intent(in) :: iz_
         real(rp), intent(in) :: tail_
         this%order(i, iz_) = nmz(iz_)
         this%converged(i, iz_) = tail_ <= this%tol
         call g_logger%info('self-energy atom '//int2str(j)//' z offset '//real2str(real(this%zoff(iz_), rp), '(f9.4)')// &
                            ' '//real2str(aimag(this%zoff(iz_)), '(f9.4)')//' Ry: moments '//int2str(nmz(iz_))// &
                            ', tail bound '//real2str(tail_, '(es10.2)'), __FILE__, __LINE__)
      end subroutine record

      !> nxt = Kt cur (first step: cur = Kt prev), later 2 Kt cur - prev, then shift prev <- cur, cur <- nxt.
      !> The common log scale of (prev, cur) grows when the new vector exceeds 1e50.
      subroutine apply_kt(first, p, c, xn)
         logical, intent(in) :: first
         complex(rp), intent(inout) :: p(:, :, :), c(:, :, :), xn(:, :, :)
         if (first) then
            call rec%onsite_vec_matmul(dz(:, :, :, iz), p, tmp)
            call rec%kubo_ham_vec_matmul(p, c, work, a=az(iz), b=real(bz(iz), rp))
            rec%izero(:) = rec%idum(:)
            c = c + tmp
            return
         end if
         call rec%onsite_vec_matmul(dz(:, :, :, iz), c, tmp)
         call rec%kubo_ham_vec_matmul(c, xn, work, a=az(iz), b=real(bz(iz), rp))
         rec%izero(:) = rec%idum(:)
         xn = 2.0_rp*(xn + tmp) - p
         m = maxval(abs(xn))
         if (m > 1.0e50_rp) then
            xn = xn/m
            c = c/m
            lscale = lscale + log(m)
         end if
         p = c
         c = xn
      end subroutine apply_kt

   end subroutine compute_gloc

   !> Write gloc_sigma.dat (rank 0): per z the offset and one nb x nb block per recursion atom; n/a if not converged.
   subroutine write_gloc(this, fname)
      class(self_energy), intent(in) :: this
      character(len=*), intent(in) :: fname
      integer :: unit, iz, i, l, kc

      if (rank /= 0) return
      open (newunit=unit, file=fname, action='write', status='replace')
      write (unit, '(a)') '# z = E_F + x + i y (Ry): x y, then G_ii (1/Ry) blocks per recursion atom '// &
         '(nb rows of nb re im pairs); n/a if not converged'
      write (unit, '(3(i0,1x))') this%nz, size(this%gloc, 3), nb
      do iz = 1, this%nz
         write (unit, '(2(es24.16e3,1x))') this%zoff(iz)
         do i = 1, size(this%gloc, 3)
            do l = 1, nb
               if (this%converged(i, iz)) then
                  write (unit, '(*(es24.16e3,1x))') (real(this%gloc(l, kc, i, iz), rp), aimag(this%gloc(l, kc, i, iz)), kc=1, nb)
               else
                  write (unit, '(a)') 'n/a'
               end if
            end do
         end do
      end do
      close (unit)
   end subroutine write_gloc

end module self_energy_mod
