!> ACBN0 potential of a generic Hermitian d-shell density matrix (code basis) under rotations.
!> The renormalisation of ACBN0 uses the diagonal of n, so only rotations that map the shell basis onto itself up to phases
!> and m -> -m keep it covariant: Rz(alpha), Ry(pi) and their products. For those, V[D n D^H] = D V[n] D^H must hold.
!> Expected values: D from ROTMAT (checked against the L operators in test_rotmat_orbital.f90); n chosen below.
!> Catches: the real-harmonic Coulomb tensor contracted with the spherical-basis n (not covariant, complex n dropped to Re n).
program test_acbn0_covariance
   use precision_mod, only: rp
   use basis_mod, only: basis_init, nb, spin_off
   use control_mod, only: control
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use logger_mod, only: g_logger
   use math_mod, only: init_math_operators, ROTMAT
   implicit none

   integer, parameter :: nd = 5
   real(rp), parameter :: tol = 1.0e-12_rp
   type(hamiltonian), target :: ham
   type(lattice), target :: lat
   type(control), target :: ctl
   complex(rp) :: n(nd, nd, 2), nr(nd, nd, 2), v(nd, nd, 2), vr(nd, nd, 2), d(nd, nd)
   real(rp) :: worst, res
   logical :: complex_n
   integer :: icase
   real(rp), parameter :: angs(2, 4) = reshape([0.0_rp, 3.141592653589793_rp, 0.8_rp, 3.141592653589793_rp, &
                                               1.1_rp, 0.0_rp, 2.3_rp, 3.141592653589793_rp], [2, 4])

   call g_logger%init()
   call basis_init(2)
   call init_math_operators()
   call ctl%restore_to_default()
   ctl%lmax = 2
   ctl%nsp = 2
   lat%ntype = 1
   allocate (lat%symbolic_atoms(1))
   call lat%symbolic_atoms(1)%restore_to_default()
   ham%lattice => lat
   ham%control => ctl
   allocate (ham%hubbard_u_pot(nb, nb, 1))
   ham%hubbard_u_pot = 0.0_rp
   ham%hubbard_u_potential_form = 'acbn0'
   lat%symbolic_atoms(1)%potential%hubbard_u(3) = 0.45_rp
   lat%symbolic_atoms(1)%potential%hubbard_j(3) = 0.06_rp

   worst = 0.0_rp
   do icase = 1, 4
      complex_n = icase > 1
      call generic_n(complex_n, n)
      call potential(n, v)
      call orbital_d(angs(1, icase), angs(2, icase), d)
      nr(:, :, 1) = matmul(d, matmul(n(:, :, 1), conjg(transpose(d))))
      nr(:, :, 2) = matmul(d, matmul(n(:, :, 2), conjg(transpose(d))))
      call potential(nr, vr)
      res = max(maxval(abs(vr(:, :, 1) - matmul(d, matmul(v(:, :, 1), conjg(transpose(d)))))), &
                maxval(abs(vr(:, :, 2) - matmul(d, matmul(v(:, :, 2), conjg(transpose(d)))))))
      print '(a,i2,a,2f8.4,a,es10.2,a,es10.2)', 'case', icase, ' alpha beta', angs(1, icase), angs(2, icase), &
         ' |V|max', maxval(abs(v)), ' covariance residual', res
      worst = max(worst, res)
   end do
   if (worst > tol) then
      print '(a,es10.2)', 'FAIL: residual', worst
      error stop 1
   end if
   print '(a)', 'test_acbn0_covariance: all checks passed'

contains

   !> Generic Hermitian n: unequal diagonal, not m <-> -m symmetric, not cubic; real symmetric when cplx is false.
   subroutine generic_n(cplx, nn)
      logical, intent(in) :: cplx
      complex(rp), intent(out) :: nn(nd, nd, 2)
      integer :: i, j, s

      nn = 0.0_rp
      do s = 1, 2
         do i = 1, nd
            nn(i, i, s) = 0.85_rp - 0.13_rp*i + 0.17_rp*(s - 1) + 0.02_rp*i*i
            do j = i + 1, nd
               nn(i, j, s) = cmplx(0.04_rp*(1 + mod(2*i + 3*j + s, 4)) - 0.05_rp, 0.0_rp, rp)
               if (cplx) nn(i, j, s) = nn(i, j, s) + cmplx(0.0_rp, 0.03_rp*(mod(i + 2*j + s, 5) - 2), rp)
               nn(j, i, s) = conjg(nn(i, j, s))
            end do
         end do
      end do
   end subroutine generic_n

   !> ACBN0 potential of the d shell for density matrices nn (spin up, spin down) through the production routine.
   subroutine potential(nn, vv)
      complex(rp), intent(in) :: nn(nd, nd, 2)
      complex(rp), intent(out) :: vv(nd, nd, 2)

      lat%symbolic_atoms(1)%potential%ldm(3, 1, 1:nd, 1:nd) = nn(:, :, 1)
      lat%symbolic_atoms(1)%potential%ldm(3, 2, 1:nd, 1:nd) = nn(:, :, 2)
      call ham%calculate_hubbard_u_potential_general()
      vv(:, :, 1) = ham%hubbard_u_pot(5:9, 5:9, 1)
      vv(:, :, 2) = ham%hubbard_u_pot(5 + spin_off:9 + spin_off, 5 + spin_off:9 + spin_off, 1)
   end subroutine potential

   !> d block of the orbital part of ROTMAT(alpha, beta, 0); the spin block is divided out through a large element.
   subroutine orbital_d(alpha, beta, dd)
      real(rp), intent(in) :: alpha, beta
      complex(rp), intent(out) :: dd(nd, nd)
      complex(rp) :: dm(2*9, 2*9)
      integer :: o

      call ROTMAT(dm, alpha, beta, 0.0_rp)
      o = merge(1, 10, abs(dm(1, 1)) > 0.5_rp)
      dd = dm(5:9, o + 4:o + 8)/dm(1, o)
   end subroutine orbital_d
end program test_acbn0_covariance
