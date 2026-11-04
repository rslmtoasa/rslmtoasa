! Spinor U insertion into the nsp=2 surface (calctype S) operator, SOC on.
!
! Two Ni(001) surface states are built from one scratch deck: input_base.nml (no U) and input_u.nml
! (U=3 eV, J=0.7 eV on the frozen host type 1 and the surface types 3 and 4; type 2 is the empty sphere).
! Frozen Hermitian spinor P matrices with complex spin coherences are placed on the U types. The
! Chebyshev recursion applies, per cluster site of type t, ee(:,:,m,t) + lsham(:,:,t) (recursion_chebyshev.f90),
! so the type-based ee and lsham of the two builds are compared: ee(:,:,1,t) of the U build minus the base build
! must equal the S1 functional's V on the d rows and columns (sigma-1)*spin_off+4+i (inserted once, not twice),
! lsham must be identical, and the off-site blocks ee(:,:,m>1,t) must be unchanged. V is evaluated here
! directly with hubbard_spinor_functional. Run from the scratch directory.
program test_hubbard_soc_surface_insertion
   use precision_mod, only: rp
   use basis_mod, only: nb, spin_off
   use control_mod, only: control
   use lattice_mod, only: lattice
   use charge_mod, only: charge
   use hamiltonian_mod, only: hamiltonian
   use hubbard_spinor_mod, only: hubbard_d_tensor, hubbard_spinor_functional
   use mpi_mod, only: g_parallel_context, parallel_context
   use timer_mod, only: g_timer, timer
   implicit none

   type(control), target :: ctl_b, ctl_u
   type(lattice), target :: lat_b, lat_u
   type(charge), target :: chg_b, chg_u
   type(hamiltonian), target :: ham_b, ham_u
   complex(rp) :: p(10, 10), v10(10, 10), t(5, 5, 5, 5)
   complex(rp), allocatable :: vfull(:, :), vtrans(:, :), vdiag(:, :), diff(:, :)
   real(rp) :: u, jbar, e_int, e_dc, e_hub
   integer(8) :: seed = 20261104_8
   logical :: failed = .false.
   integer :: it, m
   character(len=48) :: tag

   g_parallel_context = parallel_context()
   g_timer = timer()
   call build_surface_state('input_base.nml', ctl_b, lat_b, chg_b, ham_b)
   call build_surface_state('input_u.nml', ctl_u, lat_u, chg_u, ham_u)
   allocate (diff(nb, nb))

   ! Frozen Hermitian spinor P with spin coherences on every U-carrying type (host and surface layers).
   do it = 1, lat_u%ntype
      if (lat_u%symbolic_atoms(it)%potential%hubbard_u(3) == 0.0_rp .and. &
          lat_u%symbolic_atoms(it)%potential%hubbard_j(3) == 0.0_rp) cycle
      allocate (lat_u%symbolic_atoms(it)%potential%ldm_spinor(10, 10))
      lat_u%symbolic_atoms(it)%potential%ldm_spinor = rand_physical()
   end do

   call ham_u%build_bulkham()
   do it = 1, lat_u%ntype
      write (tag, '(a,i0)') 'type ', it
      u = lat_u%symbolic_atoms(it)%potential%hubbard_u(3)
      jbar = lat_u%symbolic_atoms(it)%potential%hubbard_j(3)
      diff = ham_u%ee(:, :, 1, it) - ham_b%ee(:, :, 1, it)
      if (u == 0.0_rp .and. jbar == 0.0_rp) then
         call check(trim(tag)//' (no U, ES): on-site ee difference is zero', maxval(abs(diff)), 1.0e-14_rp)
      else
         p = lat_u%symbolic_atoms(it)%potential%ldm_spinor
         call hubbard_d_tensor(u, jbar, t)
         call hubbard_spinor_functional(t, u, jbar, p, e_int, e_dc, e_hub, v10)
         vfull = embed(v10)
         vtrans = embed(transpose(v10))
         v10(1:5, 6:10) = (0.0_rp, 0.0_rp)
         v10(6:10, 1:5) = (0.0_rp, 0.0_rp)
         vdiag = embed(v10)
         call check(trim(tag)//': model P has spin coherences: max |V_up,dn|', &
            maxval(abs(vfull(5:9, spin_off + 5:spin_off + 9))), 1.0e-3_rp, .true.)
         call check(trim(tag)//': on-site ee(U) - ee(U=0) vs S1 V', maxval(abs(diff - vfull)), 1.0e-14_rp)
         call check(trim(tag)//': control, inserted twice (vs 2V)', maxval(abs(diff - 2.0_rp*vfull)), 1.0e-14_rp, .true.)
         call check(trim(tag)//': control, transposed V', maxval(abs(diff - vtrans)), 1.0e-14_rp, .true.)
         call check(trim(tag)//': control, spin-diagonal blocks only', maxval(abs(diff - vdiag)), 1.0e-14_rp, .true.)
      end if
      call check(trim(tag)//': l.s block lsham(U) - lsham(U=0)', maxval(abs(ham_u%lsham(:, :, it) - ham_b%lsham(:, :, it))), &
         1.0e-14_rp)
      call check(trim(tag)//': lsham is nonzero (SOC on)', maxval(abs(ham_b%lsham(:, :, it))), 1.0e-6_rp, .true.)
      diff(1, 1) = (0.0_rp, 0.0_rp)
      do m = 2, size(ham_u%ee, 3)
         diff(1, 1) = diff(1, 1) + sum(abs(ham_u%ee(:, :, m, it) - ham_b%ee(:, :, m, it)))
      end do
      call check(trim(tag)//': off-site ee blocks unchanged', abs(diff(1, 1)), 1.0e-14_rp)
   end do

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine build_surface_state(fname, control_obj, lattice_obj, charge_obj, hamiltonian_obj)
      character(*), intent(in) :: fname
      type(control), target, intent(inout) :: control_obj
      type(lattice), target, intent(inout) :: lattice_obj
      type(charge), target, intent(inout) :: charge_obj
      type(hamiltonian), target, intent(inout) :: hamiltonian_obj
      integer :: i

      control_obj = control(fname)
      lattice_obj = lattice(control_obj)
      call lattice_obj%build_data()
      call lattice_obj%bravais()
      call lattice_obj%build_surf_full()
      call lattice_obj%structb(.true.)
      call lattice_obj%atomlist()
      charge_obj = charge(lattice_obj)
      call charge_obj%build_alelay
      call charge_obj%surfmat
      hamiltonian_obj = hamiltonian(charge_obj)
      do i = 1, lattice_obj%ntype
         call lattice_obj%symbolic_atoms(i)%build_pot()
      end do
      call hamiltonian_obj%build_lsham()
      call hamiltonian_obj%build_bulkham()
   end subroutine build_surface_state

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

end program test_hubbard_soc_surface_insertion
