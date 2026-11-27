!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: bdg
!
!> DESCRIPTION:
!> Bogoliubov-de Gennes kernels on plain arrays (k-space, on-site singlet pairing).
!> Conventions (oracles: tests/bdg/oracles/README.md):
!>
!>   H_BdG(k) = [[H(k) - mu, D], [D^dag, -H(-k)^* + mu]],  D(iL up, iL dn) = +Delta_iL, D(iL dn, iL up) = -Delta_iL
!>   chi_iL   = sum_k w_k sum_n u_{iL up,n} conj(v_{iL dn,n}) (1 - f(E_n)),  Delta_iL = sum_Lp K_L,Lp chi_iLp
!>
!> Matrix index a runs over nmat states; site_of, orb_of and spin_of (1 = up, 2 = dn) map it to (i, L, sigma).
!------------------------------------------------------------------------------
module bdg_mod
   use precision_mod, only: rp
   use math_mod, only: hermitian_eigensystem, ry2ev
   use basis_mod, only: lmax_basis, nb, norb, spin_off
   use mpi_mod, only: rank, numprocs, ierr, g_parallel_context
   use lattice_mod, only: lattice
   use reciprocal_mod, only: reciprocal, fermi_dirac_occupation
   use logger_mod, only: g_logger
   use string_mod, only: int2str
#ifdef USE_MPI
   use mpi
#endif
   implicit none
   private

   public :: bdg, bdg_kmat_from_lambda
   public :: bdg_reciprocal_source, bdg_hamiltonian_source, bdg_assemble, bdg_accumulate, bdg_gap, bdg_iterate
   public :: bdg_dos_accumulate

   !> @brief Supplies the normal-state Hamiltonian at k and -k for one k index.
   !> @details reduce (sum chi and nel over ranks, after each k sweep) and report (once per iteration after
   !>          mixing) and eigenpairs (after the eigensolve of each k point) do nothing unless overridden.
   type, abstract :: bdg_hamiltonian_source
   contains
      procedure(hk_interface), deferred :: hk
      procedure :: reduce => reduce_none
      procedure :: report => report_none
      procedure :: eigenpairs => eigenpairs_none
   end type bdg_hamiltonian_source

   !> @brief Production source: H(k) and H(-k) from the reciprocal object on the share of this rank of its k mesh.
   type, extends(bdg_hamiltonian_source) :: bdg_reciprocal_source
      class(reciprocal), pointer :: reciprocal => null()
      !> Fractional k points and weights owned by this rank.
      real(rp), allocatable :: k_local(:, :), wk_local(:)
      integer, allocatable :: type_of(:), site_of(:), orb_of(:), spin_of(:)
      real(rp) :: n_normal = 0.0_rp
      !> DOS pass: when dos_on, eigenpairs adds to dos(n_e, ntype) on the grid dos_energies (Ry, from mu) with width dos_sigma.
      logical :: dos_on = .false.
      real(rp) :: dos_sigma = 0.0_rp
      real(rp), allocatable :: dos_energies(:), dos(:, :)
   contains
      procedure :: init => source_init
      procedure :: hk => source_hk
      procedure :: reduce => source_reduce
      procedure :: report => source_report
      procedure :: eigenpairs => source_eigenpairs
      procedure :: dos_reduce => source_dos_reduce
   end type bdg_reciprocal_source

   !> @brief &bdg input: pairing kernel K (Ry) and iteration, DOS and seed settings.
   type :: bdg
      class(reciprocal), pointer :: reciprocal => null()
      class(lattice), pointer :: lattice => null()
      !> K(L, Lp, type) in Ry.
      real(rp), allocatable :: kmat(:, :, :)
      !> Ry
      real(rp) :: delta_seed, tol, dos_window, dos_sigma
      character(len=256) :: delta_file
      real(rp) :: mix
      integer :: max_iter, n_dos
   contains
      procedure :: build_from_file
      procedure :: restore_to_default
   end type bdg

   interface bdg
      procedure :: constructor
   end interface bdg

   abstract interface
      subroutine hk_interface(this, ik, hk, hmk)
         import :: bdg_hamiltonian_source, rp
         class(bdg_hamiltonian_source), intent(in) :: this
         integer, intent(in) :: ik
         complex(rp), intent(out) :: hk(:, :), hmk(:, :)
      end subroutine hk_interface
   end interface

contains

   subroutine reduce_none(this, chi, nel)
      class(bdg_hamiltonian_source), intent(in) :: this
      complex(rp), intent(inout) :: chi(:, :)
      real(rp), intent(inout) :: nel(2)
   end subroutine reduce_none

   subroutine report_none(this, iter, dmax, delta, chi, nel)
      class(bdg_hamiltonian_source), intent(in) :: this
      integer, intent(in) :: iter
      real(rp), intent(in) :: dmax, nel(2)
      complex(rp), intent(in) :: delta(:, :), chi(:, :)
   end subroutine report_none

   subroutine eigenpairs_none(this, ik, evals, evecs)
      class(bdg_hamiltonian_source), intent(inout) :: this
      integer, intent(in) :: ik
      real(rp), intent(in) :: evals(:)
      complex(rp), intent(in) :: evecs(:, :)
   end subroutine eigenpairs_none

   subroutine source_eigenpairs(this, ik, evals, evecs)
      class(bdg_reciprocal_source), intent(inout) :: this
      integer, intent(in) :: ik
      real(rp), intent(in) :: evals(:)
      complex(rp), intent(in) :: evecs(:, :)

      if (this%dos_on) call bdg_dos_accumulate(evals, evecs, this%wk_local(ik), this%dos_energies, this%dos_sigma, &
                                               this%site_of, this%type_of, this%dos)
   end subroutine source_eigenpairs

   !> @brief Sums the DOS over ranks.
   subroutine source_dos_reduce(this)
      class(bdg_reciprocal_source), intent(inout) :: this

#ifdef USE_MPI
      if (numprocs > 1) call MPI_ALLREDUCE(MPI_IN_PLACE, this%dos, size(this%dos), MPI_DOUBLE_PRECISION, MPI_SUM, &
                                           MPI_COMM_WORLD, ierr)
#endif
   end subroutine source_dos_reduce

   !> @brief Index maps a = (i-1)*nb + L + (sigma-1)*spin_off and the slice of this rank of the reciprocal k mesh.
   subroutine source_init(this, recip)
      class(bdg_reciprocal_source), intent(inout) :: this
      class(reciprocal), target, intent(inout) :: recip
      integer :: i, L, a, first, last, count

      this%reciprocal => recip
      allocate (this%site_of(nb*recip%lattice%nrec), this%orb_of(nb*recip%lattice%nrec), this%spin_of(nb*recip%lattice%nrec))
      allocate (this%type_of(recip%lattice%nrec))
      do i = 1, recip%lattice%nrec
         this%type_of(i) = recip%lattice%ib(i)
         do L = 1, norb
            a = (i - 1)*nb + L
            this%site_of(a) = i
            this%orb_of(a) = L
            this%spin_of(a) = 1
            this%site_of(a + spin_off) = i
            this%orb_of(a + spin_off) = L
            this%spin_of(a + spin_off) = 2
         end do
      end do
      call g_parallel_context%split_range(recip%k_workset%nk_global, first, last, count)
      this%k_local = recip%k_workset%points(:, first:last)
      this%wk_local = recip%k_workset%weights(first:last)
   end subroutine source_init

   subroutine source_hk(this, ik, hk, hmk)
      class(bdg_reciprocal_source), intent(in) :: this
      integer, intent(in) :: ik
      complex(rp), intent(out) :: hk(:, :), hmk(:, :)

      call this%reciprocal%build_hamiltonian_at_kpoint(this%k_local(:, ik), hk)
      call this%reciprocal%build_hamiltonian_at_kpoint(-this%k_local(:, ik), hmk)
   end subroutine source_hk

   subroutine source_reduce(this, chi, nel)
      class(bdg_reciprocal_source), intent(in) :: this
      complex(rp), intent(inout) :: chi(:, :)
      real(rp), intent(inout) :: nel(2)

#ifdef USE_MPI
      if (numprocs > 1) then
         call MPI_ALLREDUCE(MPI_IN_PLACE, chi, size(chi), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, nel, 2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      end if
#endif
   end subroutine source_reduce

   !> @brief Rank 0 appends: iteration, max|dDelta|, N - N_normal, then |Delta_l| and |chi_l| per type and l (orbital and site average).
   subroutine source_report(this, iter, dmax, delta, chi, nel)
      class(bdg_reciprocal_source), intent(in) :: this
      integer, intent(in) :: iter
      real(rp), intent(in) :: dmax, nel(2)
      complex(rp), intent(in) :: delta(:, :), chi(:, :)
      integer :: funit, t, l, i, L1, nt, nl, count
      real(rp) :: dsum, csum
      character(len=32) :: tag
      character(len=:), allocatable :: header

      if (rank /= 0) return
      nt = this%reciprocal%lattice%ntype
      nl = lmax_basis + 1
      if (iter == 1) then
         header = '# iter max_dDelta[Ry] N-N_normal[e]'
         do t = 1, nt
            do l = 0, nl - 1
               write (tag, '(a,i0,a,i0)') ' |Delta|[Ry]_t', t, '_l', l
               header = header//trim(tag)
               write (tag, '(a,i0,a,i0)') ' |chi|[Ry]_t', t, '_l', l
               header = header//trim(tag)
            end do
         end do
         open (newunit=funit, file='bdg_state.dat', action='write', status='replace')
         write (funit, '(a)') header
      else
         open (newunit=funit, file='bdg_state.dat', action='write', status='old', position='append')
      end if
      write (funit, '(i6,2es24.16)', advance='no') iter, dmax, sum(nel) - this%n_normal
      do t = 1, nt
         do l = 0, nl - 1
            dsum = 0.0_rp
            csum = 0.0_rp
            count = 0
            do i = 1, size(this%type_of)
               if (this%type_of(i) /= t) cycle
               do L1 = l*l + 1, (l + 1)**2
                  dsum = dsum + abs(delta(L1, i))
                  csum = csum + abs(chi(L1, i))
                  count = count + 1
               end do
            end do
            if (count > 0) then
               write (funit, '(2es24.16)', advance='no') dsum/count, csum/count
            else
               write (funit, '(2a24)', advance='no') 'n/a', 'n/a'
            end if
         end do
      end do
      write (funit, '(a)') ''
      close (funit)
   end subroutine source_report

   !> @brief K(L, Lp, t) = lambda(l(L)+1, l(Lp)+1, t)/ry2ev, with L in the order s, p (3), d (5), f (7).
   subroutine bdg_kmat_from_lambda(lambda, kmat)
      real(rp), intent(in) :: lambda(:, :, :)
      real(rp), intent(out) :: kmat(:, :, :)
      integer :: i, j, l_of(size(kmat, 1))

      do i = 1, size(kmat, 1)
         l_of(i) = int(sqrt(real(i - 1, rp)))
      end do
      do j = 1, size(kmat, 2)
         do i = 1, size(kmat, 1)
            kmat(i, j, :) = lambda(l_of(i) + 1, l_of(j) + 1, :)/ry2ev
         end do
      end do
   end subroutine bdg_kmat_from_lambda

   function constructor(reciprocal_obj) result(obj)
      type(bdg) :: obj
      type(reciprocal), target, intent(inout) :: reciprocal_obj

      obj%reciprocal => reciprocal_obj
      obj%lattice => reciprocal_obj%lattice
      call obj%restore_to_default()
      call obj%build_from_file()
   end function constructor

   subroutine restore_to_default(this)
      class(bdg), intent(inout) :: this

      this%delta_seed = 1.0e-3_rp/ry2ev
      this%delta_file = ''
      this%mix = 0.3_rp
      this%tol = 1.0e-6_rp/ry2ev
      this%max_iter = 50
      this%dos_window = 0.01_rp/ry2ev
      this%n_dos = 401
      this%dos_sigma = 1.0e-4_rp/ry2ev
   end subroutine restore_to_default

   subroutine build_from_file(this)
      class(bdg), intent(inout) :: this

      integer :: iostatus, funit

      include 'include_codes/namelists/bdg.f90'

      allocate (lambda(lmax_basis + 1, lmax_basis + 1, this%lattice%ntype))
      lambda = 0.0_rp
      delta_seed = this%delta_seed*ry2ev
      delta_file = this%delta_file
      mix = this%mix
      tol = this%tol*ry2ev
      max_iter = this%max_iter
      dos_window = this%dos_window*ry2ev
      n_dos = this%n_dos
      dos_sigma = this%dos_sigma*ry2ev

      open (newunit=funit, file=this%lattice%control%fname, action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) then
         call g_logger%fatal('file '//trim(this%lattice%control%fname)//' not found', __FILE__, __LINE__)
      end if
      read (funit, nml=bdg, iostat=iostatus)
      if (iostatus /= 0 .and. .not. is_iostat_end(iostatus)) then
         call g_logger%fatal('error while reading namelist bdg, iostat = '//int2str(iostatus), __FILE__, __LINE__)
      end if
      close (funit)

      if (mix <= 0.0_rp .or. tol <= 0.0_rp .or. max_iter < 0 .or. n_dos < 1 .or. dos_window <= 0.0_rp .or. dos_sigma <= 0.0_rp) then
         call g_logger%fatal('bdg needs mix > 0, tol > 0, max_iter >= 0, n_dos >= 1, dos_window > 0 and dos_sigma > 0', &
                             __FILE__, __LINE__)
      end if
      this%delta_seed = delta_seed/ry2ev
      this%delta_file = delta_file
      this%mix = mix
      this%tol = tol/ry2ev
      this%max_iter = max_iter
      this%dos_window = dos_window/ry2ev
      this%n_dos = n_dos
      this%dos_sigma = dos_sigma/ry2ev
      allocate (this%kmat((lmax_basis + 1)**2, (lmax_basis + 1)**2, size(lambda, 3)))
      call bdg_kmat_from_lambda(lambda, this%kmat)
   end subroutine build_from_file

   !> @brief Index of the spin-down partner of an up state with the same site and orbital.
   pure integer function partner(a, site_of, orb_of, spin_of)
      integer, intent(in) :: a, site_of(:), orb_of(:), spin_of(:)
      do partner = 1, size(site_of)
         if (spin_of(partner) == 2 .and. site_of(partner) == site_of(a) .and. orb_of(partner) == orb_of(a)) return
      end do
   end function partner

   !> @brief Builds H_BdG(k) (2 nmat x 2 nmat) from H(k), H(-k), mu and Delta(norb, nsite).
   subroutine bdg_assemble(hk, hmk, mu, delta, site_of, orb_of, spin_of, hbdg)
      complex(rp), intent(in) :: hk(:, :), hmk(:, :), delta(:, :)
      real(rp), intent(in) :: mu
      integer, intent(in) :: site_of(:), orb_of(:), spin_of(:)
      complex(rp), intent(out) :: hbdg(:, :)
      integer :: n, a, b

      n = size(hk, 1)
      hbdg = (0.0_rp, 0.0_rp)
      hbdg(1:n, 1:n) = hk
      hbdg(n + 1:2*n, n + 1:2*n) = -conjg(hmk)
      do a = 1, n
         hbdg(a, a) = hbdg(a, a) - mu
         hbdg(n + a, n + a) = hbdg(n + a, n + a) + mu
         if (spin_of(a) /= 1) cycle
         b = partner(a, site_of, orb_of, spin_of)
         hbdg(a, n + b) = delta(orb_of(a), site_of(a))
         hbdg(b, n + a) = -delta(orb_of(a), site_of(a))
         hbdg(n + b, a) = conjg(delta(orb_of(a), site_of(a)))
         hbdg(n + a, b) = -conjg(delta(orb_of(a), site_of(a)))
      end do
   end subroutine bdg_assemble

   !> @brief Adds one k point (weight wk) to chi(norb, nsite) and to the electron counts nel(2) = (up, dn).
   subroutine bdg_accumulate(evals, evecs, wk, kT, site_of, orb_of, spin_of, chi, nel)
      real(rp), intent(in) :: evals(:), wk, kT
      complex(rp), intent(in) :: evecs(:, :)
      integer, intent(in) :: site_of(:), orb_of(:), spin_of(:)
      complex(rp), intent(inout) :: chi(:, :)
      real(rp), intent(inout) :: nel(2)
      integer :: n, a, b, m
      real(rp) :: f

      n = size(site_of)
      do m = 1, size(evals)
         f = fermi_dirac_occupation(evals(m), 0.0_rp, kT)
         do a = 1, n
            nel(spin_of(a)) = nel(spin_of(a)) + wk*f*abs(evecs(a, m))**2
            if (spin_of(a) /= 1) cycle
            b = partner(a, site_of, orb_of, spin_of)
            chi(orb_of(a), site_of(a)) = chi(orb_of(a), site_of(a)) &
                                         + wk*evecs(a, m)*conjg(evecs(n + b, m))*(1.0_rp - f)
         end do
      end do
   end subroutine bdg_accumulate

   !> @brief Adds one k point (weight wk) to the DOS per type: dos(ie, t) += wk sum_n sum_{a in t} |u_{a,n}|^2 g(energies(ie) - E_n).
   !> @details u = evecs(1:nmat, :) is the electron part, g the normalized Gaussian of width sigma. evals and energies
   !>          are measured from mu, in the same units as sigma.
   subroutine bdg_dos_accumulate(evals, evecs, wk, energies, sigma, site_of, type_of, dos)
      real(rp), intent(in) :: evals(:), wk, energies(:), sigma
      complex(rp), intent(in) :: evecs(:, :)
      integer, intent(in) :: site_of(:), type_of(:)
      real(rp), intent(inout) :: dos(:, :)
      integer :: m, a, t
      real(rp) :: wt(size(dos, 2)), g(size(energies))
      real(rp), parameter :: sqrt_two_pi = 2.5066282746310002_rp

      do m = 1, size(evals)
         wt = 0.0_rp
         do a = 1, size(site_of)
            t = type_of(site_of(a))
            wt(t) = wt(t) + abs(evecs(a, m))**2
         end do
         g = exp(-0.5_rp*((energies - evals(m))/sigma)**2)/(sigma*sqrt_two_pi)
         do t = 1, size(dos, 2)
            dos(:, t) = dos(:, t) + wk*wt(t)*g
         end do
      end do
   end subroutine bdg_dos_accumulate

   !> @brief Gap equation: delta(:, i) = K(type_of(i)) chi(:, i), with kmat(norb, norb, ntype) in Ry.
   subroutine bdg_gap(kmat, type_of, chi, delta)
      real(rp), intent(in) :: kmat(:, :, :)
      integer, intent(in) :: type_of(:)
      complex(rp), intent(in) :: chi(:, :)
      complex(rp), intent(out) :: delta(:, :)
      integer :: i

      do i = 1, size(type_of)
         delta(:, i) = matmul(kmat(:, :, type_of(i)), chi(:, i))
      end do
   end subroutine bdg_gap

   !> @brief Self-consistent loop in Delta: pass over all k, gap equation, linear mixing.
   !> @details max_iter = 0 does one pass at the given delta and leaves delta unchanged. chi and nel are those of
   !>          the last pass; dmax is max|delta_new - delta_old| of the last iteration.
   subroutine bdg_iterate(src, wk, mu, kT, kmat, type_of, site_of, orb_of, spin_of, mix, tol, max_iter, &
                          delta, chi, nel, niter, dmax)
      class(bdg_hamiltonian_source), intent(inout) :: src
      real(rp), intent(in) :: wk(:), mu, kT, kmat(:, :, :), mix, tol
      integer, intent(in) :: type_of(:), site_of(:), orb_of(:), spin_of(:), max_iter
      complex(rp), intent(inout) :: delta(:, :)
      complex(rp), intent(out) :: chi(:, :)
      real(rp), intent(out) :: nel(2), dmax
      integer, intent(out) :: niter
      complex(rp), allocatable :: hk(:, :), hmk(:, :), hbdg(:, :), delta_new(:, :)
      real(rp), allocatable :: evals(:)
      integer :: n, ik

      n = size(site_of)
      allocate (hk(n, n), hmk(n, n), hbdg(2*n, 2*n), evals(2*n), delta_new(size(delta, 1), size(delta, 2)))
      niter = 0
      dmax = 0.0_rp
      do
         chi = (0.0_rp, 0.0_rp)
         nel = 0.0_rp
         do ik = 1, size(wk)
            call src%hk(ik, hk, hmk)
            call bdg_assemble(hk, hmk, mu, delta, site_of, orb_of, spin_of, hbdg)
            call hermitian_eigensystem(hbdg, evals)
            call src%eigenpairs(ik, evals, hbdg)
            call bdg_accumulate(evals, hbdg, wk(ik), kT, site_of, orb_of, spin_of, chi, nel)
         end do
         call src%reduce(chi, nel)
         if (max_iter == 0) exit
         call bdg_gap(kmat, type_of, chi, delta_new)
         dmax = maxval(abs(delta_new - delta))
         delta = delta + mix*(delta_new - delta)
         niter = niter + 1
         call src%report(niter, dmax, delta, chi, nel)
         if (dmax < tol .or. niter >= max_iter) exit
      end do
   end subroutine bdg_iterate

end module bdg_mod
