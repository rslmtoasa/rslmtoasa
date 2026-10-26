! Spinor U insertion into the nsp=2 Hamiltonian (hubbard_soc_experimental).
!
! Two bcc-Fe nsp=2 states are built from one scratch deck: input_base.nml (no U) and input_u.nml
! (U=2 eV, J=0.5 eV, flag on). A frozen Hermitian spinor P with complex spin coherences is placed
! on the site, and the difference of the two built Hamiltonians (RS ee, KS H(k) at two k points)
! must equal the S1 functional's V on the d rows and columns (sigma-1)*spin_off+4+i and be zero
! elsewhere. V is evaluated here directly with hubbard_spinor_functional. The insertion count is
! one when the difference equals V and not 2V. Run from the scratch directory.
program test_hubbard_soc_insertion
   use precision_mod, only: rp
   use basis_mod, only: nb, spin_off
   use control_mod, only: control
   use lattice_mod, only: lattice
   use charge_mod, only: charge
   use energy_mod, only: energy
   use hamiltonian_mod, only: hamiltonian
   use reciprocal_mod, only: reciprocal
   use hubbard_spinor_mod, only: hubbard_d_tensor, hubbard_spinor_functional
   use mpi_mod, only: g_parallel_context, parallel_context
   use timer_mod, only: g_timer, timer
   use spin_response_fe_state, only: build_fe_state
   implicit none

   type(control), target :: ctl_b, ctl_u
   type(lattice), target :: lat_b, lat_u
   type(charge), target :: chg_b, chg_u
   type(energy), target :: en_b, en_u
   type(hamiltonian), target :: ham_b, ham_u
   type(reciprocal), target :: rec_b, rec_u
   complex(rp) :: p(10, 10), v10(10, 10), t(5, 5, 5, 5)
   complex(rp), allocatable :: vfull(:, :), vtrans(:, :), vdiag(:, :), hk_b(:, :), hk_u(:, :), diff(:, :), ref(:, :), coll_pot(:, :)
   real(rp) :: u, jbar, e_int, e_dc, e_hub, kpts(3, 2)
   integer(8) :: seed = 20261027_8
   logical :: failed = .false.
   integer :: i, j, s, r, ik

   g_parallel_context = parallel_context()
   g_timer = timer()
   call build_fe_state(.false., ctl_b, lat_b, chg_b, en_b, ham_b, 'input_base.nml')
   call build_fe_state(.false., ctl_u, lat_u, chg_u, en_u, ham_u, 'input_u.nml')

   allocate (hk_b(nb, nb), hk_u(nb, nb), diff(nb, nb), ref(nb, nb), coll_pot(nb, nb))
   u = lat_u%symbolic_atoms(1)%potential%hubbard_u(3)
   jbar = lat_u%symbolic_atoms(1)%potential%hubbard_j(3)
   p = rand_physical()
   allocate (lat_u%symbolic_atoms(1)%potential%ldm_spinor(10, 10))
   lat_u%symbolic_atoms(1)%potential%ldm_spinor = p
   call hubbard_d_tensor(u, jbar, t)
   call hubbard_spinor_functional(t, u, jbar, p, e_int, e_dc, e_hub, v10)
   vfull = embed(v10)
   vtrans = embed(transpose(v10))
   v10(1:5, 6:10) = (0.0_rp, 0.0_rp)
   v10(6:10, 1:5) = (0.0_rp, 0.0_rp)
   vdiag = embed(v10)
   call check('model P has spin coherences: max |V_up,dn|', &
      maxval(abs(vfull(5:9, spin_off + 5:spin_off + 9))), 1.0e-3_rp, .true.)

   ! RS: bulk on-site block ee(:,:,1,1)
   call ham_u%build_bulkham()
   diff = ham_u%ee(:, :, 1, 1) - ham_b%ee(:, :, 1, 1)
   call check('RS ee: inserted d block vs S1 V (all other entries zero)', maxval(abs(diff - vfull)), 1.0e-14_rp)
   call check('RS ee: control, count of insertions (difference vs 2V)', maxval(abs(diff - 2.0_rp*vfull)), 1.0e-14_rp, .true.)
   call check('RS ee: control, transposed V', maxval(abs(diff - vtrans)), 1.0e-14_rp, .true.)
   call check('RS ee: control, spin-diagonal blocks only', maxval(abs(diff - vdiag)), 1.0e-14_rp, .true.)

   ! KS: H(k) at Gamma and a generic point; the on-site V is k independent
   kpts(:, 1) = [0.0_rp, 0.0_rp, 0.0_rp]
   kpts(:, 2) = [0.13_rp, -0.21_rp, 0.07_rp]
   rec_b = reciprocal(ham_b)
   rec_u = reciprocal(ham_u)
   call rec_b%generate_mp_mesh()
   call rec_u%generate_mp_mesh()
   do ik = 1, 2
      call rec_b%build_hamiltonian_at_kpoint(kpts(:, ik), hk_b)
      call rec_u%build_hamiltonian_at_kpoint(kpts(:, ik), hk_u)
      diff = hk_u - hk_b
      call check('KS H(k): inserted d block vs S1 V (all other entries zero)', maxval(abs(diff - vfull)), 1.0e-14_rp)
      call check('KS H(k): control, transposed V', maxval(abs(diff - vtrans)), 1.0e-14_rp, .true.)
      call check('KS H(k): control, spin-diagonal blocks only', maxval(abs(diff - vdiag)), 1.0e-14_rp, .true.)
   end do

   ! (b) collinear limit: zero spin coherences reproduce the collinear potential
   p(1:5, 6:10) = (0.0_rp, 0.0_rp)
   p(6:10, 1:5) = (0.0_rp, 0.0_rp)
   lat_u%symbolic_atoms(1)%potential%ldm_spinor = p
   do s = 1, 2
      lat_u%symbolic_atoms(1)%potential%ldm(3, s, 1:5, 1:5) = p((s - 1)*5 + 1:s*5, (s - 1)*5 + 1:s*5)
   end do
   call ham_u%calculate_hubbard_u_potential_general()
   ref = ham_u%hubbard_u_pot(:, :, 1)
   ham_u%hubbard_soc_experimental = .false.
   call ham_u%calculate_hubbard_u_potential_general()
   coll_pot = ham_u%hubbard_u_pot(:, :, 1)
   ham_u%hubbard_soc_experimental = .true.
   call check('(b) zero spin coherences: spinor fill vs the collinear path', maxval(abs(ref - coll_pot)), 1.0e-14_rp)
   p = rand_physical()
   lat_u%symbolic_atoms(1)%potential%ldm_spinor = p
   call ham_u%calculate_hubbard_u_potential_general()
   call check('(b) control, spin coherences present: spinor fill vs the collinear path', &
      maxval(abs(ham_u%hubbard_u_pot(:, :, 1) - coll_pot)), 1.0e-14_rp, .true.)

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   function uniform() result(x)
      real(rp) :: x

      seed = mod(seed*1103515245_8 + 12345_8, 2147483648_8)
      x = real(seed, rp)/2147483648.0_rp
   end function uniform

   function rand_physical() result(pm)
      complex(rp) :: pm(10, 10), q(10, 10)
      integer :: a, b, c

      do a = 1, 10
         do b = 1, 10
            q(b, a) = cmplx(uniform() - 0.5_rp, uniform() - 0.5_rp, rp)
         end do
         do c = 1, a - 1
            q(:, a) = q(:, a) - dot_product(q(:, c), q(:, a))*q(:, c)
         end do
         q(:, a) = q(:, a)/sqrt(real(dot_product(q(:, a), q(:, a)), rp))
      end do
      pm = (0.0_rp, 0.0_rp)
      do a = 1, 10
         pm = pm + (0.1_rp + 0.8_rp*uniform())*spread(q(:, a), 2, 10)*spread(conjg(q(:, a)), 1, 10)
      end do
   end function rand_physical

   !> Place a 10x10 spinor block on the d rows and columns (sigma-1)*spin_off + 4 + i.
   function embed(m10) result(m)
      complex(rp), intent(in) :: m10(10, 10)
      complex(rp) :: m(nb, nb)
      integer :: sr, sc

      m = (0.0_rp, 0.0_rp)
      do sr = 1, 2
         do sc = 1, 2
            m((sr - 1)*spin_off + 5:(sr - 1)*spin_off + 9, (sc - 1)*spin_off + 5:(sc - 1)*spin_off + 9) = &
               m10((sr - 1)*5 + 1:sr*5, (sc - 1)*5 + 1:sc*5)
         end do
      end do
   end function embed

   subroutine check(name, measured, tol, control_check)
      character(len=*), intent(in) :: name
      real(rp), intent(in) :: measured, tol
      logical, intent(in), optional :: control_check
      logical :: ok, is_control

      is_control = .false.
      if (present(control_check)) is_control = control_check
      ok = merge(measured > tol, measured <= tol, is_control)
      if (.not. ok) failed = .true.
      if (is_control) then
         write (*, '(a,1x,a,2(2x,es10.3))') merge('PASS', 'FAIL', ok), trim(name)//' [control: must exceed]', measured, tol
      else
         write (*, '(a,1x,a,2(2x,es10.3))') merge('PASS', 'FAIL', ok), trim(name), measured, tol
      end if
   end subroutine check

end program test_hubbard_soc_insertion
