!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: recursion_gpu_mod
!
! DESCRIPTION:
!> GPU (or threaded-CPU) backend for recursion_mod via librsrec.
!>
!> The library exposes one context per MPI rank holding the device-resident
!> Hamiltonian, and four drivers replacing the inner work of:
!>
!>     recur()                     -> rsgpu_scalar_lanczos
!>     recur_b(), recur_b_ij()     -> rsgpu_block_lanczos
!>     chebyshev_recur(),
!>     chebyshev_recur_ij()        -> rsgpu_chebyshev_moments
!>     compute_moments_stochastic()-> rsgpu_stochastic_moments
!>
!> Integration pattern inside recursion_mod (e.g. in chebyshev_recur):
!>
!>     use recursion_gpu_mod
!>     type(rsgpu) :: gpu
!>     ...
!>     call gpu%init(this%lattice%kk, nb, size(this%lattice%nn, 2), &
!>                   this%lattice%ntype, this%lattice%nmax, device=0)
!>     call gpu%set_hamiltonian(this%hamiltonian%ee, this%hamiltonian%hall, &
!>                              this%hamiltonian%lsham, this%lattice%nn,    &
!>                              this%lattice%iz)
!>     do i = start_atom, end_atom          ! MPI decomposition unchanged
!>        j = this%lattice%irec(i)
!>        psi0 = (0.0_rp, 0.0_rp)
!>        do l = 1, nb
!>           psi0(l, l, j) = (1.0_rp, 0.0_rp)
!>        end do
!>        call gpu%chebyshev_moments(psi0, this%lattice%control%lld, a, b, &
!>                                   this%mu_n(:, :, :, g2l_map(i)))
!>     end do
!>     call gpu%finalize()
!>
!> The izero/idum lightcone bookkeeping of the CPU path is not needed:
!> atoms outside the lightcone carry exactly zero wavefunction, so the
!> dense GPU sweep returns identical results.
!>
!> The hoh path (ham_hoh_vec_matmul / hop_b_hoh) is NOT yet implemented in
!> the library; drivers correspond to the hoh = .false. branches. Guard the
!> GPU dispatch with `if (.not. this%hamiltonian%hoh)` until then.
!------------------------------------------------------------------------------
module recursion_gpu_mod

   use, intrinsic :: iso_c_binding
   implicit none

   private
   public :: rsgpu

   integer, parameter :: rp = selected_real_kind(15, 307)

   type :: rsgpu
      type(c_ptr) :: ctx = c_null_ptr
      integer :: kk = 0, nb = 0, nnmax = 0, ntype = 0, nmax = 0
   contains
      procedure :: init => rsgpu_init
      procedure :: set_hamiltonian => rsgpu_set_hamiltonian
      procedure :: set_velocity => rsgpu_set_velocity
      procedure :: set_grid => rsgpu_set_grid
      procedure :: set_precision => rsgpu_set_precision
      procedure :: chebyshev_dos => rsgpu_chebyshev_dos
      procedure :: ham_apply => rsgpu_ham_apply
      procedure :: chebyshev_moments => rsgpu_chebyshev_moments
      procedure :: block_lanczos => rsgpu_block_lanczos
      procedure :: scalar_lanczos => rsgpu_scalar_lanczos
      procedure :: stochastic_moments => rsgpu_stochastic_moments
      procedure :: finalize => rsgpu_finalize
   end type rsgpu

   interface
      function c_rsrec_create(kk, nb, nnmax, ntype, nmax, device) &
         bind(C, name="rsrec_create") result(ctx)
         import :: c_ptr, c_int
         integer(c_int), value :: kk, nb, nnmax, ntype, nmax, device
         type(c_ptr) :: ctx
      end function

      subroutine c_rsrec_destroy(ctx) bind(C, name="rsrec_destroy")
         import :: c_ptr
         type(c_ptr), value :: ctx
      end subroutine

      function c_rsrec_last_error() bind(C, name="rsrec_last_error") &
         result(msg)
         import :: c_ptr
         type(c_ptr) :: msg
      end function

      function c_rsrec_set_hamiltonian(ctx, ee, hall, lsham, nn, iz, &
                                       eeo, hallo, enim) &
         bind(C, name="rsrec_set_hamiltonian") result(ierr)
         import :: c_ptr, c_int
         type(c_ptr), value :: ctx, ee, hall, lsham
         type(c_ptr), value :: nn, iz
         type(c_ptr), value :: eeo, hallo, enim
         integer(c_int) :: ierr
      end function

      function c_rsrec_set_velocity(ctx, va, vb) &
         bind(C, name="rsrec_set_velocity") result(ierr)
         import :: c_ptr, c_int
         type(c_ptr), value :: ctx, va, vb
         integer(c_int) :: ierr
      end function

      function c_rsrec_set_grid(ctx, coords, use_structured) &
         bind(C, name="rsrec_set_grid") result(ierr)
         import :: c_ptr, c_int
         type(c_ptr), value :: ctx, coords
         integer(c_int), value :: use_structured
         integer(c_int) :: ierr
      end function

      function c_rsrec_set_precision(ctx, prec) &
         bind(C, name="rsrec_set_precision") result(ierr)
         import :: c_ptr, c_int
         type(c_ptr), value :: ctx
         integer(c_int), value :: prec
         integer(c_int) :: ierr
      end function

      function c_rsrec_chebyshev_dos(ctx, mu, n_mom, natoms, ene, nv, &
                                     a, b, g0) &
         bind(C, name="rsrec_chebyshev_dos") result(ierr)
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: ctx, mu, ene, g0
         integer(c_int), value :: n_mom, natoms, nv
         real(c_double), value :: a, b
         integer(c_int) :: ierr
      end function

      function c_rsrec_op_apply(ctx, which, x, y, nrhs, a, b) &
         bind(C, name="rsrec_op_apply") result(ierr)
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: ctx, x, y
         integer(c_int), value :: which, nrhs
         real(c_double), value :: a, b
         integer(c_int) :: ierr
      end function

      function c_rsrec_chebyshev_moments(ctx, psi0, lld, a, b, mu) &
         bind(C, name="rsrec_chebyshev_moments") result(ierr)
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: ctx, psi0, mu
         integer(c_int), value :: lld
         real(c_double), value :: a, b
         integer(c_int) :: ierr
      end function

      function c_rsrec_block_lanczos(ctx, psi0, lld, a_b, b2_b) &
         bind(C, name="rsrec_block_lanczos") result(ierr)
         import :: c_ptr, c_int
         type(c_ptr), value :: ctx, psi0, a_b, b2_b
         integer(c_int), value :: lld
         integer(c_int) :: ierr
      end function

      function c_rsrec_scalar_lanczos(ctx, site_j, lld, a, b2) &
         bind(C, name="rsrec_scalar_lanczos") result(ierr)
         import :: c_ptr, c_int
         type(c_ptr), value :: ctx, a, b2
         integer(c_int), value :: site_j, lld
         integer(c_int) :: ierr
      end function

      function c_rsrec_stochastic_moments(ctx, psiref, lld, a, b, mu_nm) &
         bind(C, name="rsrec_stochastic_moments") result(ierr)
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: ctx, psiref, mu_nm
         integer(c_int), value :: lld
         real(c_double), value :: a, b
         integer(c_int) :: ierr
      end function
   end interface

contains

   subroutine die(where)
      use, intrinsic :: iso_c_binding
      character(len=*), intent(in) :: where
      character(len=512) :: msg
      character(kind=c_char), pointer :: cmsg(:)
      type(c_ptr) :: p
      integer :: i
      p = c_rsrec_last_error()
      msg = ""
      if (c_associated(p)) then
         call c_f_pointer(p, cmsg, [512])
         do i = 1, 512
            if (cmsg(i) == c_null_char) exit
            msg(i:i) = cmsg(i)
         end do
      end if
      write (*, "(a)") "recursion_gpu_mod: "//where//": "//trim(msg)
      error stop 1
   end subroutine die

   subroutine rsgpu_init(this, kk, nb, nnmax, ntype, nmax, device)
      class(rsgpu), intent(inout) :: this
      integer, intent(in) :: kk, nb, nnmax, ntype, nmax
      integer, intent(in), optional :: device
      integer :: dev
      dev = 0
      if (present(device)) dev = device
      this%kk = kk; this%nb = nb; this%nnmax = nnmax
      this%ntype = ntype; this%nmax = nmax
      this%ctx = c_rsrec_create(int(kk, c_int), int(nb, c_int), &
                                int(nnmax, c_int), int(ntype, c_int), &
                                int(nmax, c_int), int(dev, c_int))
      if (.not. c_associated(this%ctx)) call die("init")
   end subroutine rsgpu_init

   subroutine rsgpu_set_hamiltonian(this, ee, hall, lsham, nn, iz, nmax, eeo, hallo, enim)
      class(rsgpu), intent(inout) :: this
      complex(rp), dimension(:, :, :, :), intent(in), target :: ee
      complex(rp), dimension(:, :, :, :), intent(in), target, optional :: hall
      complex(rp), dimension(:, :, :), intent(in), target, optional :: lsham
      integer(c_int), dimension(:, :), intent(in), target :: nn
      integer(c_int), dimension(:), intent(in), target :: iz
      integer, intent(in), optional :: nmax   ! struct member; accepted for caller compat
      complex(rp), dimension(:, :, :, :), intent(in), target, optional :: eeo
      complex(rp), dimension(:, :, :, :), intent(in), target, optional :: hallo
      complex(rp), dimension(:, :, :), intent(in), target, optional :: enim
      type(c_ptr) :: p_hall, p_ls, p_eeo, p_hallo, p_enim
      p_hall = c_null_ptr
      p_ls = c_null_ptr
      p_eeo = c_null_ptr
      p_hallo = c_null_ptr
      p_enim = c_null_ptr
      if (present(hall)) p_hall = c_loc(hall)
      if (present(lsham)) p_ls = c_loc(lsham)
      if (present(eeo)) p_eeo = c_loc(eeo)
      if (present(hallo)) p_hallo = c_loc(hallo)
      if (present(enim)) p_enim = c_loc(enim)
      if (c_rsrec_set_hamiltonian(this%ctx, c_loc(ee), p_hall, p_ls, &
                                  c_loc(nn), c_loc(iz), p_eeo, p_hallo, p_enim) /= 0) &
         call die("set_hamiltonian")
   end subroutine rsgpu_set_hamiltonian

   subroutine rsgpu_set_velocity(this, v_a, v_b)
      class(rsgpu), intent(inout) :: this
      complex(rp), dimension(:, :, :, :), intent(in), target :: v_a, v_b
      if (c_rsrec_set_velocity(this%ctx, c_loc(v_a), c_loc(v_b)) /= 0) &
         call die("set_velocity")
   end subroutine rsgpu_set_velocity


   !> Enable the structured FFT-stencil + correction backend. coords are
   !> (3, kk) integer lattice coordinates of the clust atoms; pass
   !> use_structured = .false. to revert to the block-ELL backend.
   !> Call after set_hamiltonian (and set_velocity, if used).
   subroutine rsgpu_set_grid(this, coords, use_structured)
      class(rsgpu), intent(inout) :: this
      integer(c_int), dimension(:, :), intent(in), target :: coords
      logical, intent(in), optional :: use_structured
      integer(c_int) :: flag
      flag = 1_c_int
      if (present(use_structured)) then
         if (.not. use_structured) flag = 0_c_int
      end if
      if (c_rsrec_set_grid(this%ctx, c_loc(coords), flag) /= 0) &
         call die("set_grid")
   end subroutine rsgpu_set_grid

   !> Chebyshev arithmetic precision: 0 = fp32 engine (default, fast on
   !> all GPUs, ~1e-6 moment accuracy), 1 = fp64 (bit-comparable reference)
   subroutine rsgpu_set_precision(this, prec)
      class(rsgpu), intent(inout) :: this
      integer, intent(in) :: prec
      if (c_rsrec_set_precision(this%ctx, int(prec, c_int)) /= 0) &
         call die("set_precision")
   end subroutine rsgpu_set_precision

   !> Green function / DOS reconstruction, drop-in for the per-atom body
   !> of chebyshev_green: pass the LOCAL moment and g0 slices, e.g.
   !>   call gpu%chebyshev_dos(rec%mu_n(:, :, :, 1:n_loc), &
   !>                          en%ene(1:nv), a, b, green%g0(:, :, 1:nv, 1:n_loc))
   !> with a = (emax_win - emin_win)/(2 - 0.3_rp), b = (emax_win + emin_win)/2
   !> from resolve_chebyshev_window, replacing the whole n/ie/i loop nest.
   subroutine rsgpu_chebyshev_dos(this, mu_n, ene, a, b, g0)
      class(rsgpu), intent(inout) :: this
      complex(rp), dimension(:, :, :, :), intent(in), target :: mu_n
      real(rp), dimension(:), intent(in), target :: ene
      real(rp), intent(in) :: a, b
      complex(rp), dimension(:, :, :, :), intent(inout), target :: g0
      integer(c_int) :: n_mom, natoms, nv
      n_mom = int(size(mu_n, 3), c_int)
      natoms = int(size(mu_n, 4), c_int)
      nv = int(size(ene), c_int)
      if (size(g0, 3) < nv .or. size(g0, 4) < natoms) &
         call die("chebyshev_dos: g0 too small")
      if (c_rsrec_chebyshev_dos(this%ctx, c_loc(mu_n), n_mom, natoms, &
                                c_loc(ene), nv, real(a, c_double), &
                                real(b, c_double), c_loc(g0)) /= 0) &
         call die("chebyshev_dos")
   end subroutine rsgpu_chebyshev_dos

   !> psi_out = (H psi_in - b psi_in)/a  -- drop-in for ham_vec_matmul
   subroutine rsgpu_ham_apply(this, psi_in, psi_out, a, b)
      class(rsgpu), intent(inout) :: this
      complex(rp), dimension(:, :, :), intent(in), target :: psi_in
      complex(rp), dimension(:, :, :), intent(out), target :: psi_out
      real(rp), intent(in) :: a, b
      if (c_rsrec_op_apply(this%ctx, 0_c_int, c_loc(psi_in), &
                           c_loc(psi_out), int(size(psi_in, 2), c_int), &
                           real(a, c_double), real(b, c_double)) /= 0) &
         call die("ham_apply")
   end subroutine rsgpu_ham_apply

   !> mu_n(:,:,1:2*lld+2) for one starting block state (chebyshev_recur[_ij])
   subroutine rsgpu_chebyshev_moments(this, psi0, lld, a, b, mu_n)
      class(rsgpu), intent(inout) :: this
      complex(rp), dimension(:, :, :), intent(in), target :: psi0
      integer, intent(in) :: lld
      real(rp), intent(in) :: a, b
      complex(rp), dimension(:, :, :), intent(inout), target :: mu_n
      if (size(mu_n, 3) < 2 * lld + 2) call die("chebyshev_moments: mu_n too small")
      if (c_rsrec_chebyshev_moments(this%ctx, c_loc(psi0), int(lld, c_int), &
                                    real(a, c_double), real(b, c_double), &
                                    c_loc(mu_n)) /= 0) &
         call die("chebyshev_moments")
   end subroutine rsgpu_chebyshev_moments

   !> atemp_b/b2temp_b for one starting block state (recur_b / recur_b_ij)
   subroutine rsgpu_block_lanczos(this, psi0, lld, atemp_b, b2temp_b)
      class(rsgpu), intent(inout) :: this
      complex(rp), dimension(:, :, :), intent(in), target :: psi0
      integer, intent(in) :: lld
      complex(rp), dimension(:, :, :), intent(inout), target :: atemp_b, b2temp_b
      if (c_rsrec_block_lanczos(this%ctx, c_loc(psi0), int(lld, c_int), &
                                c_loc(atemp_b), c_loc(b2temp_b)) /= 0) &
         call die("block_lanczos")
   end subroutine rsgpu_block_lanczos

   !> a(ll,l), b2(ll,l) for all nb orbital chains of site j (recur)
   subroutine rsgpu_scalar_lanczos(this, site_j, lld, a, b2)
      class(rsgpu), intent(inout) :: this
      integer, intent(in) :: site_j, lld
      real(rp), dimension(:, :), intent(inout), target :: a, b2
      if (c_rsrec_scalar_lanczos(this%ctx, int(site_j, c_int), &
                                 int(lld, c_int), c_loc(a), c_loc(b2)) /= 0) &
         call die("scalar_lanczos")
   end subroutine rsgpu_scalar_lanczos

   !> mu_nm(:,:,n,m) for one reference state (compute_moments_stochastic)
   subroutine rsgpu_stochastic_moments(this, psiref, lld, a, b, mu_nm)
      class(rsgpu), intent(inout) :: this
      complex(rp), dimension(:, :, :), intent(in), target :: psiref
      integer, intent(in) :: lld
      real(rp), intent(in) :: a, b
      complex(rp), dimension(:, :, :, :), intent(inout), target :: mu_nm
      if (c_rsrec_stochastic_moments(this%ctx, c_loc(psiref), &
                                     int(lld, c_int), real(a, c_double), &
                                     real(b, c_double), c_loc(mu_nm)) /= 0) &
         call die("stochastic_moments")
   end subroutine rsgpu_stochastic_moments

   subroutine rsgpu_finalize(this)
      class(rsgpu), intent(inout) :: this
      if (c_associated(this%ctx)) call c_rsrec_destroy(this%ctx)
      this%ctx = c_null_ptr
   end subroutine rsgpu_finalize

end module recursion_gpu_mod
