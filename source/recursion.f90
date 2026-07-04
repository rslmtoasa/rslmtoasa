 !------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Recursion
!
!> @author
!> Angela Klautau
!> Ramon Cardias
!> Lucas P. Campagna
!> S. Frota-Pessôa
!> Pascoal R. Peduto
!> Anders Bergman
!> S. B. Legoas
!> H. M. Petrilli
!> Ivan P. Miranda
!
! DESCRIPTION:
!> Module to handle recursion calculations
!------------------------------------------------------------------------------

module recursion_mod

   use hamiltonian_mod
   use chebyshev_fast_mod, only: cheb_moments_fast, cheb_moments_fast_batched, &
      cheb_moments_fast_mkl_batch, &
      cheb_moments_fast_mkl_sparse, cheb_fast_reset_cache, &
      cheb_moments_stochastic_fast, cheb_moments_orbital_fast
   use haydock_fast_mod, only: block_lanczos_fast
   use lattice_mod
   use control_mod
   use energy_mod
   use sparse_mod, only: sparse
   use precision_mod, only: rp
   use math_mod
   use string_mod
   use logger_mod, only: g_logger
   use rsrec_cuda_plugin_mod, only: rsrec_cuda_backend, get_gpu_context, &
      rsrec_cuda_plugin_compiled, decode_gpu_backend, gpu_backend_bsr, &
      gpu_backend_conv, gpu_backend_fft
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc
#endif
   use timer_mod, only: g_timer
   use basis_mod, only: nb, norb, spin_off, lmax_basis
   implicit none

   private

   !> Module´s main structure
   type, public :: recursion
      !> Hamiltonian
      class(hamiltonian), pointer :: hamiltonian
      !> Sparse algebra helper
      type(sparse), allocatable :: sparse_obj
      !> Lattice
      class(lattice), pointer :: lattice
      !> Energy
      class(energy), pointer :: en
      !> Control
      class(control), pointer :: control
      ! General variables

      !> Scalar recursion coefficients
      real(rp), dimension(:, :, :, :), allocatable :: a, b2
      real(rp), dimension(:), allocatable :: atemp, b2temp
      !> Block recursion coefficients
      complex(rp), dimension(:, :, :, :), allocatable :: a_b, b2_b
      complex(rp), dimension(:, :, :), allocatable :: atemp_b, b2temp_b
      !> Atom list in hopping region
      integer, dimension(:), allocatable :: izero, idum, irlist, nzero, mzero, ndum, mdum
      integer :: irnum !< Atoms within recursion shells
      !> Atom list in hopping region as a function of recursion step
      integer, dimension(:, :), allocatable :: izeroll
      !> Wave functions for recursion hopping (Lanczos)
      complex(rp), dimension(:, :), allocatable :: psi, pmn
      !> Wave functions for block recursion hopping (Lanczos)
      complex(rp), dimension(:, :, :), allocatable :: psi_b, pmn_b, hpsi, hohpsi, enupsi, socpsi
      !> Wave functions for recursion hopping (Chebyshev)
      complex(rp), dimension(:, :, :), allocatable :: psi0, psi1, psi2
      !> Cached full first-order operator blocks for no-hoh fast/GPU CCOR paths
      complex(rp), dimension(:, :, :, :), allocatable :: ee_ccor_work, hall_ccor_work
      !> Chebyshev moments
      complex(rp), dimension(:, :, :, :), allocatable :: mu_n, mu_ng
      complex(rp), dimension(:, :, :, :, :), allocatable :: mu_nm_stochastic
      complex(rp), allocatable :: cheb_mom_temp(:,:)
      !> Variable to save H|psi>
      complex(rp), dimension(:, :), allocatable :: v
      !> Variable to save T(H)
      complex(rp), dimension(:, :, :, :, :), allocatable :: t_h
      type(rsrec_cuda_backend), pointer :: gpu_backend => null()
   contains
      procedure :: hop
      procedure :: crecal
      procedure :: recur
      procedure :: hop_b
      procedure :: hop_b_hoh
      procedure :: crecal_b
      procedure :: recur_b
      procedure :: recur_b_ij
      procedure :: zsqr
      procedure :: bpopt
      procedure :: get_terminf
      procedure :: get_cinf
      procedure :: emami
      procedure :: create_ll_map
      procedure :: cheb_0th_mom
      procedure :: cheb_1st_mom
      procedure :: cheb_1st_mom_hoh
      procedure :: resolve_chebyshev_window
      procedure :: chebyshev_recur_ij
      procedure :: chebyshev_recur_ll
      procedure :: chebyshev_recur_ll_hoh
      procedure :: chebyshev_orbital_mod
      procedure :: chebyshev_recur
      procedure :: compute_moments_stochastic
      procedure :: ham_vec_matmul
      procedure :: ham_hoh_vec_matmul
      procedure :: velo_vec_matmul
      procedure :: velo_hoh_vec_matmul
      procedure :: restore_to_default
      final :: destructor
   end type recursion

   interface recursion
      procedure :: constructor
   end interface recursion


   interface
   module subroutine cheb_moments_cpu(this, psi0, lld, a, b, mu)
      class(recursion), intent(inout) :: this
      complex(rp), intent(in) :: psi0(:, :, :)
      integer, intent(in) :: lld
      real(rp), intent(in) :: a, b
      complex(rp), intent(inout) :: mu(:, :, :)
      character(len=:), allocatable :: backend
      character(len=256) :: msg

   end subroutine cheb_moments_cpu

   module subroutine cheb_moments_legacy(this, psi0, lld, a, b, mu)
      class(recursion), intent(inout) :: this
      complex(rp), intent(in) :: psi0(:, :, :)
      integer, intent(in) :: lld
      real(rp), intent(in) :: a, b
      complex(rp), intent(inout) :: mu(:, :, :)
      ! Local variables
      integer :: k, ll

      ! Rebuild the initial lightcone from the support of psi0. The leaf
      ! routines propagate it dynamically through this%izero/this%idum, so the
      ! izeroll/create_ll_map bookkeeping of the old inline drivers is not
      ! needed (those arrays are never read by the leaf routines).
   end subroutine cheb_moments_legacy

   module subroutine resolve_chebyshev_window(this, emin, emax)
      class(recursion), intent(inout) :: this
      real(rp), intent(out) :: emin, emax
      character(len=256) :: msg

   end subroutine resolve_chebyshev_window

   module function constructor(hamiltonian_obj, energy_obj, sparse_obj) result(obj)
      type(recursion) :: obj
      type(hamiltonian), target, intent(in) :: hamiltonian_obj
      type(energy), target, intent(in) :: energy_obj
      type(sparse), intent(in) :: sparse_obj

   end function constructor

   module logical function gpu_plugin_enabled(this)
      class(recursion), intent(in) :: this

   end function gpu_plugin_enabled

   module logical function gpu_periodic_backend(this)
      class(recursion), intent(in) :: this
      integer :: backend

   end function gpu_periodic_backend

   module logical function gpu_plugin_ready(this, feature, require_nsp1, allow_hoh)
      class(recursion), intent(in) :: this
      character(len=*), intent(in) :: feature
      logical, intent(in), optional :: require_nsp1
      logical, intent(in), optional :: allow_hoh
      integer :: backend
      logical :: nsp1_only, hoh_ok

   end function gpu_plugin_ready

   module subroutine gpu_plugin_upload_hamiltonian(this)
      class(recursion), intent(inout) :: this
      integer :: backend

   end subroutine gpu_plugin_upload_hamiltonian

   module subroutine ensure_ccor_operator_blocks(this)
      class(recursion), intent(inout) :: this

   end subroutine ensure_ccor_operator_blocks

   module subroutine block_lanczos_fast_nohoh(this, llmax, use_ccor)
      class(recursion), intent(inout) :: this
      integer, intent(in) :: llmax
      logical, intent(in) :: use_ccor

   end subroutine block_lanczos_fast_nohoh

   module subroutine destructor(this)
      type(recursion) :: this
   end subroutine destructor

   module subroutine velo_vec_matmul(this, c_or_n, v_op, psi_in, psi_out)
      class(recursion), intent(inout) :: this
      complex(rp), dimension(:, :, :, :), intent(in) :: v_op
      complex(rp), dimension(:, :, :), intent(in) :: psi_in
      character :: c_or_n 
      complex(rp), dimension(:, :, :), intent(out) :: psi_out
      ! Local variables
      integer :: ineigh, ih, j, k, nr, m, n, l, hblocksize, nat, nnmap, nlimplus1
      complex(rp), dimension(nb, nb) :: dum1, dum2, locham
   end subroutine velo_vec_matmul

   module subroutine velo_hoh_vec_matmul(this, v_op, vo_op, psi_in, psi_out)
      class(recursion), intent(inout) :: this
      complex(rp), dimension(:, :, :), intent(in) :: psi_in
      complex(rp), dimension(:, :, :), intent(out) :: psi_out
      complex(rp), dimension(:, :, :, :), intent(in) :: v_op
      complex(rp), dimension(:, :, :, :), intent(in) :: vo_op
      ! Local variables
      integer :: ineigh, ih, j, k, nr, m, n, l, hblocksize, nat, nnmap, nlimplus1
      complex(rp), dimension(nb, nb) :: dum1, dum2, locham

   end subroutine velo_hoh_vec_matmul

   module subroutine ham_hoh_vec_matmul(this, psi_in, psi_out, a, b)
      class(recursion), intent(inout) :: this
      complex(rp), dimension(:, :, :), intent(in) :: psi_in
      complex(rp), dimension(:, :, :), intent(out) :: psi_out
      real(rp), intent(in) :: a, b
      ! Local variables
      integer :: ineigh, ih, j, k, nr, m, n, l, hblocksize, nat, nnmap, nlimplus1
      complex(rp), dimension(nb, nb) :: dum1, dum2, locham

   end subroutine ham_hoh_vec_matmul

   module subroutine ham_vec_matmul(this, psi_in, psi_out, a, b)
      class(recursion), intent(inout) :: this
      complex(rp), dimension(:, :, :), intent(in) :: psi_in
      complex(rp), dimension(:, :, :), intent(out) :: psi_out
      real(rp), intent(in) :: a, b
      ! Local variables
      integer :: ineigh, ih, j, k, nr, m, n, l, hblocksize, nat, nnmap, nlimplus1
      complex(rp), dimension(nb, nb) :: dum1, dum2, locham
   end subroutine ham_vec_matmul

   module subroutine compute_moments_stochastic(this)
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: ineigh, ih, i, j, k, nr, ll, m, n, l, hblocksize, nat, nnmap, loop_over, ie, lmax, ntype
      complex(rp), dimension(nb, nb) :: dum, dum1, dum2
      complex(rp), dimension(:, :), allocatable :: S_op, L_op
      complex(rp), dimension(norb, norb) :: mLx, mLy, mLz
      complex(rp), dimension(:, :, :), allocatable :: psiref, w0, w1, w2, right_vec, v0, v1, v2, cn
      complex(rp), dimension(:, :, :, :), allocatable :: left_vec
      real(rp), dimension(this%en%channels_ldos + 10) :: w, wscale
      real(rp), dimension(this%control%cond_ll) :: kernel
      complex(rp), dimension(nb, nb, this%en%channels_ldos + 10) :: g0
      real(rp) :: a, b, rng, emin_win, emax_win
      complex(rp) :: exp_factor

   end subroutine compute_moments_stochastic

   module subroutine hop_b_hoh(this, ll)
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m, n, nr, nnmap, nlimplus1
      integer :: ll ! Recursion step
      integer :: ino ! Atom type
      integer :: ih ! Atom number in the clust
      integer, dimension(0:this%lattice%kk) :: idum
      complex(rp), dimension(nb, nb) :: summ
      complex(rp), dimension(nb, nb) :: locham

   end subroutine hop_b_hoh

   module subroutine hop_b(this, ll)
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m, n, nr, nnmap, nlimplus1
      integer :: ll ! Recursion step
      integer :: ino ! Atom type
      integer :: ih ! Atom number in the clust
      integer, dimension(0:this%lattice%kk) :: idum
      complex(rp), dimension(nb, nb) :: summ
      complex(rp), dimension(nb, nb) :: locham

   end subroutine hop_b

   module subroutine recur_b_ij(this)
      use mpi_mod
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, ij, j, l, ll, kk, m, reci
      integer :: llmax ! Recursion steps
      complex(rp) :: asign, bsign

      integer :: ij_loc

   end subroutine recur_b_ij

   module subroutine recur_b(this)
      use mpi_mod
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, l, ll, kk, m
      integer :: llmax ! Recursion steps
      real(rp) :: rng

      ! Determine how many atoms each process should handle
   end subroutine recur_b

   module subroutine crecal_b(this)
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m, info, lwork
      integer :: ll, llmax ! Recursion step
      integer :: nm1 ! LL-1
      integer :: nat ! Clust size
      integer :: hblocksize ! Hamiltonian size (18)
      complex(rp), dimension(nb, nb) :: sum_b, sum_a
      complex(rp), dimension(nb, nb) :: dum, b, b_i, u, lam, lam_i, u_ls, u_rs
!    complex(rp), dimension(nb, nb,this%lattice%kk) :: psi_t
      complex(rp), dimension(:, :, :), allocatable :: psi_t
      real(rp), dimension(nb) :: ev
      complex(rp), dimension(nb) :: zev
      real(rp), dimension(3*nb - 2) ::rwork
      complex(rp), dimension(3*nb - 2) ::zwork

   end subroutine crecal_b

   module subroutine zsqr(this)
      use mpi_mod
      implicit none
      class(recursion), intent(inout) :: this
      !
      integer :: na, LDIM, I, J, K, L, N, M, info, n_glob
      real(rp), dimension(nb) :: ev
      real(rp), dimension(3*nb - 2) ::rwork
      complex(rp), dimension(nb*nb) ::zwork
      complex(rp), dimension(nb, nb) :: B2t, U, DUM, lam

      !
   end subroutine zsqr

   module subroutine get_cinf(this, alpha, beta, ll, np, nw, alpha_inf, beta_inf)
      implicit none
      class(recursion), intent(inout) :: this
      !.. Formal Arguments ..
      integer, intent(in) :: np
      integer, intent(inout) :: ll, nw
      real(rp), dimension(np, ll), intent(inout) :: alpha
      real(rp), dimension(np, ll), intent(in) :: beta
      real(rp), dimension(np), intent(out) :: alpha_inf
      real(rp), dimension(np), intent(out) :: beta_inf
      !
      !.. Local Scalars ..
      integer :: icode, iii, k, l, ll1, ineigh, nbp1, nl, nt, eidx, ne, nq, ifail
      real(rp) :: a1, a2, emax, emin, eps, err, e_shift, pi
      complex(rp) :: g_e
      !
      !.. Local Arrays ..
      integer, dimension(np) :: nb2
      integer, dimension(nw) :: jc
      integer, dimension(nw) :: iwk
      real(rp), dimension(ll) :: aa, am, bb, bm, sqbb
      real(rp), dimension(10) :: edge, weight, width
      real(rp), dimension(200, 10) :: bwk
      ! real(rp), dimension(np,ll) :: am2,bm2
      real(rp), dimension(np, 10) :: edge2, width2, weight2
      real(rp), dimension(nw, 2, 5) :: work
      real(rp), dimension(2) :: bandedges

      !
      !.. External Calls ..
      !
      !.. External Functions ..
      !
      !.. Intrinsic Functions ..
      intrinsic MAX, MIN
      !
      ! ... Executable Statements ...
      !
      !**************************************************************************
   end subroutine get_cinf

   module subroutine get_terminf(this, Acoef_b, B2coef_b, na, ll, ldim, nw, a_inf, b_inf, a_inf0, b_inf0)
      implicit none
      class(recursion), intent(inout) :: this
      integer, intent(in) :: na
      integer, intent(in) :: ll
      integer, intent(inout) :: nw
      integer, intent(in) :: ldim
      complex(rp), dimension(ldim, ldim, ll, na), intent(in) :: Acoef_b, B2coef_b
      real(rp), dimension(ldim, ldim, na), intent(out) :: a_inf, b_inf
      real(rp), dimension(na), intent(out) ::  a_inf0, b_inf0
      !
      integer :: n, i, j, ll_t, k
      complex(rp), dimension(ldim, ldim) :: MatIn, MatOut
      real(rp), dimension(ldim, ldim, ll) :: Acoef_r, B2coef_r
      real(rp) :: maxA, maxB, maxAinf, maxBinf, tmpval
      logical :: foundNaN_in, foundNaN_out
      !
   end subroutine get_terminf

   module subroutine cheb_0th_mom(this, psiref)
      class(recursion), intent(inout) :: this
      complex(rp), dimension(nb, nb, this%lattice%kk), intent(in) :: psiref
      complex(rp), dimension(nb, nb) :: dum
      integer :: nat, hblocksize, nlimplus1, k

   end subroutine cheb_0th_mom

   module subroutine cheb_1st_mom(this, psiref, a, b)
      class(recursion), intent(inout) :: this
      real(rp), intent(in) :: a, b ! scale and shift
      complex(rp), dimension(nb, nb, this%lattice%kk), intent(in) :: psiref
      integer :: nat, hblocksize, nlimplus1, k, ih, nr, n, ineigh, nnmap, i

   end subroutine cheb_1st_mom

   module subroutine cheb_1st_mom_hoh(this, psiref, a, b)
      class(recursion), intent(inout) :: this
      real(rp), intent(in) :: a, b ! scale and shift
      complex(rp), dimension(nb, nb, this%lattice%kk), intent(in) :: psiref
      integer :: nat, hblocksize, nlimplus1, k, ih, nr, n, ineigh, nnmap, i

   end subroutine cheb_1st_mom_hoh

   module subroutine chebyshev_recur_ij(this)
      use mpi_mod
      class(recursion), intent(inout) :: this
      ! Local variables
      real(rp) :: a, b, emin_win, emax_win
      integer :: i, ij, j, l, ll, kk, m, reci
      integer :: llmax ! Recursion steps
      complex(rp) :: asign, bsign

      integer :: ij_loc

   end subroutine chebyshev_recur_ij

   module subroutine chebyshev_recur_ll(this, mu, ll, a, b)
      class(recursion), intent(inout) :: this
      complex(rp), intent(inout) :: mu(:, :, :)
      integer, intent(in) :: ll
      real(rp), intent(in) :: a, b
      ! Local variables
      integer :: ineigh, ih, j, k, nr, m, n, l, hblocksize, nat, nnmap, nlimplus1
      complex(rp), dimension(nb, nb) :: dum1, dum2, locham

   end subroutine chebyshev_recur_ll

   module subroutine chebyshev_recur_ll_hoh(this, mu, ll, a, b)
      class(recursion), intent(inout) :: this
      complex(rp), intent(inout) :: mu(:, :, :)
      integer, intent(in) :: ll
      real(rp), intent(in) :: a, b
      ! Local variables
      integer :: ineigh, ih, j, k, nr, m, n, l, hblocksize, nat, nnmap, nlimplus1
      complex(rp), dimension(nb, nb) :: dum1, dum2, locham

   end subroutine chebyshev_recur_ll_hoh

   module subroutine chebyshev_orbital_mod(this)
      use mpi_mod
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: atom_neighbor, ih, i, j, k, nr, ll, m, n, l, hblocksize, ntype, nnmap, nlimplus1, llcheb, lmax, nv, ie, random, iz, ia
      complex(rp), dimension(nb, nb) :: dum, dum1, dum2
      complex(rp), dimension(nb, this%lattice%kk) :: v
      complex(rp), dimension(:, :, :, :), allocatable :: hcheb
      complex(rp), dimension(:, :, :), allocatable :: psiref, w0, w1, w2, right_vec, v0, v1, v2, left_vec, left_vec1, left_vec2
      complex(rp), dimension(:, :, :), allocatable :: mu_n_orb, mu_tmp
      complex(rp), dimension(nb, nb, this%en%channels_ldos + 10) :: g0
      complex(rp), dimension(norb, norb) :: mLx, mLy, mLz
      complex(rp), dimension(nb, nb) :: L_op
      complex(rp) :: exp_factor
      real(rp) :: lz
      real(rp), dimension(this%control%lld) :: kernel
      real(rp), dimension(this%en%channels_ldos + 10) :: w, wscale, lzi
      real(rp), dimension(3) :: rij, crossrij
      
      real(rp) :: a, b, start, finish, rng, emin_win, emax_win
      ! External functions
      complex(rp), external :: zdotc

      integer :: i_loc
      real(rp) :: maxg_tmp
      integer :: ie_tmp

   end subroutine chebyshev_orbital_mod

   module subroutine chebyshev_recur(this)
      use mpi_mod
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: ineigh, ih, i, j, k, nr, ll, m, n, l, hblocksize, nat, nnmap, nlimplus1, llcheb
      integer, dimension(0:this%lattice%kk) :: idumll
      !complex(rp) :: cone = (1.0d0, 0.0d0)
      complex(rp), dimension(nb, nb) :: dum, dum1, dum2
      complex(rp), dimension(nb, this%lattice%kk) :: v
      complex(rp), dimension(:, :, :, :), allocatable :: hcheb
      real(rp) :: a, b, start, finish, emin_win, emax_win
      ! External functions
      complex(rp), external :: zdotc

      integer :: i_loc

   end subroutine chebyshev_recur

   module subroutine create_ll_map(this)
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, nr, nnmap, ll
      integer, dimension(0:this%lattice%kk) :: idumll

   end subroutine create_ll_map

   module subroutine hop(this, ll)
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m, n, nr, nnmap, nlimplus1
      integer :: ll ! Recursion step
      integer :: ino ! Atom type
      integer :: ih ! Atom number in the clust
      integer, dimension(0:this%lattice%kk) :: idum
      complex(rp), dimension(nb) :: dum
      real(rp) :: summ, start, finish

   end subroutine hop

   module subroutine crecal(this)
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m
      integer :: ll, llmax ! Recursion step
      integer :: nm1 ! LL-1
      integer :: nat ! Clust size
      integer :: hblocksize ! Hamiltonian size (18)
      real(rp) :: s, summ, start, finish
      complex(rp) :: dum, ajc
      complex(rp), dimension(nb, this%lattice%kk) :: thpsi
      character(len=1) :: transa
      character, dimension(6) :: matdescra
      ! External functions
      complex(rp), external :: zdotc

   end subroutine crecal

   module subroutine recur(this)
      use mpi_mod
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, l, ll
      integer :: llmax ! Recursion steps
      integer :: i_loc

   end subroutine recur

   module subroutine bpopt(this, ll, a, rb, n, ainf, rbinf, ifail)
      class(recursion), intent(inout) :: this
      ! Input
      integer, intent(in) :: n, ll
      real(rp), dimension(ll), intent(in) :: A, RB
      ! Output
      integer, intent(out) :: ifail
      real(rp), intent(out) :: ainf, rbinf
      ! Local variables
      integer :: I, JITER, NDIME
      real(rp) :: BM, BMAX, BMIN, EPS
      real(rp), dimension(ll) :: AZ, RBZ

   end subroutine bpopt

   module subroutine emami(this, nl, as, bs, n, emax, emin)
      class(recursion), intent(inout) :: this
      ! Input
      integer, intent(in) :: N, NL
      real(rp), dimension(NL), intent(in) :: AS, BS
      ! Output
      real(rp), intent(out) :: EMAX, EMIN
      ! Local Variables
      integer :: I, ISTOP, NUM
      real(rp) :: DELE, E, E1, E2, EMAX0, EMIN0, EPS, P, RELFEH, X1, X2
      real(rp), dimension(NL) :: A, B

   end subroutine emami

   module subroutine restore_to_default(this, full)
      class(recursion) :: this
      logical, intent(in), optional :: full
      integer :: lmax

   end subroutine restore_to_default

   end interface

end module recursion_mod
