!------------------------------------------------------------------------------
! DRESP-03R -- live LMTO exchange-vertex mapping audit fixture
!------------------------------------------------------------------------------
!
! This test is deliberately a boundary test.  It proves the live pieces that
! can be proved without inventing a path-operator -> physical-GF adapter:
!   * d_matrix(E) is dele_up*dele_down*(P_up(E)-P_down(E)), with source ordering;
!   * an independent local rotation differentiates the complete live magnetic
!     bond builder, including spin-dependent hopping and the onsite term.
!
! The missing native g=(P-S)^(-1), lambda/mu physical-GF construction is
! documented by DRESP_03R_LMTO_EXCHANGE_VERTEX_MAPPING.md and is intentionally
! not approximated here.
!------------------------------------------------------------------------------
program test_dresp03r_lmto_mapping
   use basis_mod, only: basis_init
   use lmto_magnetic_tangent_mod, only: lmto_bond_value, lmto_hhmag_to_spinor
   use math_mod, only: i_unit
   use precision_mod, only: rp
   use symbolic_atom_mod, only: symbolic_atom
   implicit none

   integer, parameter :: lmax = 1
   integer, parameter :: norb = (lmax + 1)**2
   integer, parameter :: nsite = 2
   integer, parameter :: site_spinor = 2*norb
   integer, parameter :: nmat = nsite*site_spinor
   integer, parameter :: ne = 5
   real(rp), parameter :: derivative_step = 1.0e-5_rp
   real(rp), parameter :: derivative_tolerance = 3.0e-9_rp
   real(rp), parameter :: d_formula_tolerance = 5.0e-14_rp
   logical :: failed

   failed = .false.
   call basis_init(lmax)
   call test_d_matrix_is_delta_p()
   call test_full_local_rotation_derivative()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains

   subroutine test_d_matrix_is_delta_p()
      type(symbolic_atom) :: atom
      complex(rp) :: dmat(norb, norb), pmat(2*norb, 2*norb, ne)
      real(rp) :: energies(ne), max_error, raw_error
      integer :: ie, i, l, m, idx

      call atom%restore_to_default()
      atom%potential%lmax = lmax
      atom%potential%vmad = 0.19_rp
      do l = 0, lmax
         atom%potential%c(l, 1) = -0.71_rp + 0.13_rp*real(l, rp)
         atom%potential%c(l, 2) = -0.43_rp - 0.07_rp*real(l, rp)
         atom%potential%dele(l, 1) = 0.66_rp + 0.05_rp*real(l, rp)
         atom%potential%dele(l, 2) = 0.81_rp - 0.03_rp*real(l, rp)
      end do
      energies = [-0.87_rp, -0.21_rp, 0.17_rp, 0.74_rp, 1.31_rp]
      pmat = cmplx(0.0_rp, 0.0_rp, rp)
      call atom%p_matrix(pmat, lmax, energies)

      max_error = 0.0_rp
      raw_error = 0.0_rp
      do ie = 1, ne
         call atom%d_matrix(dmat, energies(ie))
         do l = 0, lmax
            do m = 1, 2*l + 1
               idx = l*l + m
               max_error = max(max_error, abs(dmat(idx, idx) - &
                  atom%potential%dele(l, 1)*atom%potential%dele(l, 2)* &
                  (pmat(idx, idx, ie) - pmat(norb + idx, norb + idx, ie))))
               raw_error = max(raw_error, abs(dmat(idx, idx) - &
                  (pmat(idx, idx, ie) - pmat(norb + idx, norb + idx, ie))))
            end do
         end do
      end do
      write (*, '(a,es12.4)') 'd_matrix(E) - dele_up*dele_down DeltaP max_err    = ', max_error
      write (*, '(a,es12.4)') 'd_matrix(E) - DeltaP raw difference max_err        = ', raw_error
      call check_true('native d_matrix is endpoint-scaled DeltaP', max_error <= d_formula_tolerance)
      call check_true('native d_matrix is not raw DeltaP for unequal widths', raw_error > 1.0e-3_rp)
   end subroutine test_d_matrix_is_delta_p

   subroutine test_full_local_rotation_derivative()
      complex(rp) :: hhh_ii(norb, norb), hhh_jj(norb, norb), hhh_ij(norb, norb)
      complex(rp) :: h0(nmat, nmat), hplus(nmat, nmat), hminus(nmat, nmat), dh_analytic(nmat, nmat)
      complex(rp) :: hmag(norb, norb, 4), dhmag(norb, norb, 4), block(site_spinor, site_spinor)
      complex(rp) :: wx0_i(norb), wx1_i(norb), wx0_j(norb), wx1_j(norb)
      complex(rp) :: c0_i(norb), c1_i(norb), c0_j(norb), c1_j(norb)
      real(rp) :: mom_i(3), mom_j(3), mom_plus(3), mom_minus(3), dm(3), axis(3)
      real(rp) :: error, offsite_derivative, commutator_norm
      integer :: first_i, last_i, first_j, last_j

      ! A non-diagonal structure block and orbital-dependent widths prevent a
      ! commuting/scalar reduction.  All inputs have the live lmto_bond_value
      ! type and units; this is not an H_up-H_down surrogate.
      hhh_ii = cmplx(0.0_rp, 0.0_rp, rp)
      hhh_jj = cmplx(0.0_rp, 0.0_rp, rp)
      hhh_ij = cmplx(0.0_rp, 0.0_rp, rp)
      hhh_ii(1, 1) = 0.41_rp; hhh_ii(2, 2) = -0.23_rp; hhh_ii(3, 3) = 0.17_rp; hhh_ii(4, 4) = -0.08_rp
      hhh_jj(1, 1) = -0.37_rp; hhh_jj(2, 2) = 0.29_rp; hhh_jj(3, 3) = -0.12_rp; hhh_jj(4, 4) = 0.33_rp
      hhh_ij(1, 1) = cmplx(0.21_rp, 0.03_rp, rp)
      hhh_ij(2, 2) = cmplx(-0.16_rp, -0.02_rp, rp)
      hhh_ij(3, 3) = cmplx(0.11_rp, 0.04_rp, rp)
      hhh_ij(4, 4) = cmplx(-0.09_rp, 0.01_rp, rp)
      hhh_ij(1, 2) = cmplx(0.07_rp, 0.02_rp, rp)
      hhh_ij(2, 1) = cmplx(-0.05_rp, 0.01_rp, rp)
      hhh_ij(3, 4) = cmplx(0.06_rp, -0.03_rp, rp)
      hhh_ij(4, 3) = cmplx(-0.04_rp, 0.02_rp, rp)

      wx0_i = [1.00_rp, 0.83_rp, 1.17_rp, 0.91_rp]
      wx1_i = [0.18_rp, -0.12_rp, 0.09_rp, -0.07_rp]
      wx0_j = [0.94_rp, 1.11_rp, 0.88_rp, 1.06_rp]
      wx1_j = [-0.11_rp, 0.08_rp, -0.05_rp, 0.13_rp]
      c0_i = [-0.62_rp, -0.31_rp, 0.12_rp, 0.28_rp]
      c1_i = [0.16_rp, -0.09_rp, 0.22_rp, -0.13_rp]
      c0_j = [-0.49_rp, -0.27_rp, 0.19_rp, 0.34_rp]
      c1_j = [0.11_rp, -0.14_rp, 0.18_rp, -0.08_rp]
      mom_i = [0.0_rp, 0.0_rp, 1.0_rp]
      mom_j = [sin(0.37_rp), 0.0_rp, cos(0.37_rp)]
      axis = [0.0_rp, 1.0_rp, 0.0_rp]
      dm = cross3(axis, mom_i)
      mom_plus = rotate_moment(mom_i, axis, derivative_step)
      mom_minus = rotate_moment(mom_i, axis, -derivative_step)

      call build_two_site_hamiltonian(mom_i, mom_j, h0, hhh_ii, hhh_jj, hhh_ij, wx0_i, wx1_i, wx0_j, wx1_j, &
                                      c0_i, c1_i, c0_j, c1_j)
      call build_two_site_hamiltonian(mom_plus, mom_j, hplus, hhh_ii, hhh_jj, hhh_ij, wx0_i, wx1_i, wx0_j, wx1_j, &
                                      c0_i, c1_i, c0_j, c1_j)
      call build_two_site_hamiltonian(mom_minus, mom_j, hminus, hhh_ii, hhh_jj, hhh_ij, wx0_i, wx1_i, wx0_j, wx1_j, &
                                      c0_i, c1_i, c0_j, c1_j)
      call analytic_two_site_derivative(mom_i, mom_j, dm, dh_analytic, hhh_ii, hhh_ij, wx0_i, wx1_i, wx0_j, wx1_j, c1_i)
      error = maxval(abs((hplus - hminus)/(2.0_rp*derivative_step) - dh_analytic))

      first_i = 1; last_i = site_spinor
      first_j = last_i + 1; last_j = nmat
      offsite_derivative = max(maxval(abs(dh_analytic(first_i:last_i, first_j:last_j))), &
                               maxval(abs(dh_analytic(first_j:last_j, first_i:last_i))))
      commutator_norm = sqrt(sum(abs(matmul(diagonal_matrix(wx0_i), hhh_ij) - &
                                  matmul(hhh_ij, diagonal_matrix(wx0_j)))**2))

      write (*, '(a,es12.4)') 'full local-rotation derivative central error = ', error
      write (*, '(a,es12.4)') 'offsite derivative norm (not onsite-only)   = ', offsite_derivative
      write (*, '(a,es12.4)') 'width/structure noncommutator norm          = ', commutator_norm
      call check_true('live bond builder local derivative', error <= derivative_tolerance)
      call check_true('local rotation produces offsite vertex blocks', offsite_derivative > 1.0e-5_rp)
      call check_true('fixture has noncommuting structure/width factors', commutator_norm > 1.0e-5_rp)
   end subroutine test_full_local_rotation_derivative

   subroutine build_two_site_hamiltonian(mom_i, mom_j, hamiltonian, hhh_ii, hhh_jj, hhh_ij, &
                                         wx0_i, wx1_i, wx0_j, wx1_j, c0_i, c1_i, c0_j, c1_j)
      real(rp), intent(in) :: mom_i(3), mom_j(3)
      complex(rp), intent(out) :: hamiltonian(:, :)
      complex(rp), intent(in) :: hhh_ii(:, :), hhh_jj(:, :), hhh_ij(:, :)
      complex(rp), intent(in) :: wx0_i(:), wx1_i(:), wx0_j(:), wx1_j(:), c0_i(:), c1_i(:), c0_j(:), c1_j(:)
      complex(rp) :: hmag(norb, norb, 4)
      complex(rp) :: bond(site_spinor, site_spinor)
      complex(rp) :: zero_potential(norb)

      hamiltonian = cmplx(0.0_rp, 0.0_rp, rp)
      zero_potential = cmplx(0.0_rp, 0.0_rp, rp)
      call lmto_bond_value(hhh_ii, wx0_i, wx1_i, wx0_i, wx1_i, c0_i, c1_i, mom_i, mom_i, .true., hmag)
      call lmto_hhmag_to_spinor(hmag, bond)
      hamiltonian(1:site_spinor, 1:site_spinor) = bond
      call lmto_bond_value(hhh_jj, wx0_j, wx1_j, wx0_j, wx1_j, c0_j, c1_j, mom_j, mom_j, .true., hmag)
      call lmto_hhmag_to_spinor(hmag, bond)
      hamiltonian(site_spinor + 1:nmat, site_spinor + 1:nmat) = bond
      call lmto_bond_value(hhh_ij, wx0_i, wx1_i, wx0_j, wx1_j, zero_potential, zero_potential, mom_i, mom_j, .false., hmag)
      call lmto_hhmag_to_spinor(hmag, bond)
      hamiltonian(1:site_spinor, site_spinor + 1:nmat) = bond
      hamiltonian(site_spinor + 1:nmat, 1:site_spinor) = transpose(conjg(bond))
   end subroutine build_two_site_hamiltonian

   subroutine analytic_two_site_derivative(mom_i, mom_j, dm, derivative, hhh_ii, hhh_ij, &
                                           wx0_i, wx1_i, wx0_j, wx1_j, c1_i)
      real(rp), intent(in) :: mom_i(3), mom_j(3), dm(3)
      complex(rp), intent(out) :: derivative(:, :)
      complex(rp), intent(in) :: hhh_ii(:, :), hhh_ij(:, :), wx0_i(:), wx1_i(:), wx0_j(:), wx1_j(:), c1_i(:)
      complex(rp) :: dhmag(norb, norb, 4), bond(site_spinor, site_spinor)
      complex(rp) :: zero_potential(norb)

      derivative = cmplx(0.0_rp, 0.0_rp, rp)
      zero_potential = cmplx(0.0_rp, 0.0_rp, rp)
      call bond_derivative(hhh_ii, wx0_i, wx1_i, wx0_i, wx1_i, c1_i, mom_i, mom_i, dm, dm, .true., dhmag)
      call lmto_hhmag_to_spinor(dhmag, bond)
      derivative(1:site_spinor, 1:site_spinor) = bond
      call bond_derivative(hhh_ij, wx0_i, wx1_i, wx0_j, wx1_j, zero_potential, mom_i, mom_j, dm, &
                           [0.0_rp, 0.0_rp, 0.0_rp], .false., dhmag)
      call lmto_hhmag_to_spinor(dhmag, bond)
      derivative(1:site_spinor, site_spinor + 1:nmat) = bond
      derivative(site_spinor + 1:nmat, 1:site_spinor) = transpose(conjg(bond))
   end subroutine analytic_two_site_derivative

   pure subroutine bond_derivative(hhh, wx0_i, wx1_i, wx0_j, wx1_j, c1_i, mom_i, mom_j, dm_i, dm_j, onsite, dhmag)
      complex(rp), intent(in) :: hhh(:, :), wx0_i(:), wx1_i(:), wx0_j(:), wx1_j(:), c1_i(:)
      real(rp), intent(in) :: mom_i(3), mom_j(3), dm_i(3), dm_j(3)
      logical, intent(in) :: onsite
      complex(rp), intent(out) :: dhmag(:, :, :)
      real(rp) :: dcross(3), cross_i(3), cross_j(3), ddot
      integer :: i, j, idir

      dhmag = cmplx(0.0_rp, 0.0_rp, rp)
      ddot = dot_product(dm_i, mom_j) + dot_product(mom_i, dm_j)
      cross_i = cross3(dm_i, mom_j)
      cross_j = cross3(mom_i, dm_j)
      dcross = cross_i + cross_j
      do j = 1, size(hhh, 2)
         do i = 1, size(hhh, 1)
            dhmag(i, j, 4) = wx1_i(i)*hhh(i, j)*wx1_j(j)*ddot
            do idir = 1, 3
               dhmag(i, j, idir) = wx1_i(i)*hhh(i, j)*wx0_j(j)*dm_i(idir) + &
                  wx0_i(i)*hhh(i, j)*wx1_j(j)*dm_j(idir) + &
                  i_unit*wx1_i(i)*hhh(i, j)*wx1_j(j)*dcross(idir)
            end do
         end do
      end do
      if (.not. onsite) return
      do i = 1, size(hhh, 1)
         do idir = 1, 3
            dhmag(i, i, idir) = dhmag(i, i, idir) + c1_i(i)*cmplx(dm_i(idir), 0.0_rp, rp)
         end do
      end do
   end subroutine bond_derivative

   pure function rotate_moment(moment, axis, angle) result(rotated)
      real(rp), intent(in) :: moment(3), axis(3), angle
      real(rp) :: rotated(3)
      rotated = moment*cos(angle) + cross3(axis, moment)*sin(angle) + axis*dot_product(axis, moment)*(1.0_rp-cos(angle))
   end function rotate_moment

   pure function cross3(a, b) result(c)
      real(rp), intent(in) :: a(3), b(3)
      real(rp) :: c(3)
      c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
   end function cross3

   pure function diagonal_matrix(values) result(matrix)
      complex(rp), intent(in) :: values(:)
      complex(rp) :: matrix(size(values), size(values))
      integer :: i
      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, size(values)
         matrix(i, i) = values(i)
      end do
   end function diagonal_matrix

   subroutine check_true(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      if (.not. condition) then
         write (*, '(a)') 'FAIL: '//trim(label)
         failed = .true.
      end if
   end subroutine check_true

end program test_dresp03r_lmto_mapping
