!------------------------------------------------------------------------------
! LR-REP-00 representation audit
!------------------------------------------------------------------------------
!
!> @brief Test the live LMTO Green-function representation conversions.
!> @details The test deliberately keeps the representation labels explicit:
!>          the reciprocal and native-RS bridge is first checked in the
!>          orthogonal coefficient space, then both blocks are sent through
!>          the same endpoint and screening conversions before comparison.
!------------------------------------------------------------------------------
program test_lr_rs_gf_representation
   use basis_mod, only: basis_init, nb
   use dyson_kernel_mod, only: dyson_kspace_inverse
   use energy_mod, only: energy
   use green_mod, only: green
   use lattice_mod, only: lattice
   use lehmann_kernel_mod, only: lehmann_pair_block
   use math_mod, only: i_unit, two_pi
   use precision_mod, only: rp
   use symbolic_atom_mod, only: symbolic_atom
   implicit none

   integer, parameter :: lmax = 1
   integer, parameter :: norb_site = (lmax + 1)**2
   integer, parameter :: nsite = 2
   integer, parameter :: nmat = nsite*2*norb_site
   integer, parameter :: nk = 32
   integer, parameter :: ne = 5
   real(rp), parameter :: tol = 2.0e-11_rp
   real(rp), parameter :: bridge_tol = 3.0e-10_rp
   real(rp), parameter :: negative_floor = 1.0e-5_rp
   logical :: failed
   type(green) :: green_obj
   type(lattice), target :: lattice_obj
   type(energy), target :: energy_obj
   type(symbolic_atom), target :: atoms(2)
   real(rp) :: alpha(0:lmax, 2), beta(0:lmax, 2)
   complex(rp) :: z_grid(ne)
   complex(rp), allocatable :: p_alpha_i(:, :, :), p_beta_i(:, :, :), p_alpha_j(:, :, :), p_beta_j(:, :, :)
   complex(rp), allocatable :: g_alpha(:, :, :), g_beta(:, :, :), g_roundtrip(:, :, :), g_wrong(:, :, :)

   failed = .false.
   call basis_init(lmax)
   call setup_green(green_obj, lattice_obj, energy_obj, atoms, z_grid)
   call setup_screening(alpha, beta)
   allocate (p_alpha_i(nb, nb, ne), p_beta_i(nb, nb, ne), p_alpha_j(nb, nb, ne), p_beta_j(nb, nb, ne), &
             g_alpha(nb, nb, ne), g_beta(nb, nb, ne), g_roundtrip(nb, nb, ne), g_wrong(nb, nb, ne))
   call build_potential_matrices(z_grid, alpha, beta, 0.11_rp, p_alpha_i, p_beta_i)
   call build_potential_matrices(z_grid, alpha, beta, 0.37_rp, p_alpha_j, p_beta_j)
   call build_deterministic_gf(g_alpha)

   call test_auxiliary_endpoint(green_obj, g_alpha)
   call test_round_trip(green_obj, alpha, beta, p_alpha_i, p_beta_i, p_alpha_j, p_beta_j, &
                        g_alpha, g_beta, g_roundtrip, g_wrong)
   call test_bridge(green_obj, alpha, beta, p_alpha_i, p_beta_i, p_alpha_j, p_beta_j, z_grid)
   call test_high_energy(green_obj, z_grid, alpha, beta)

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains

   subroutine setup_green(gobj, lat, en, symbolic_atoms, energies)
      type(green), intent(out) :: gobj
      type(lattice), target, intent(out) :: lat
      type(energy), target, intent(out) :: en
      type(symbolic_atom), target, intent(out) :: symbolic_atoms(:)
      complex(rp), intent(out) :: energies(:)
      integer :: ia, ie, l, s

      do ia = 1, size(symbolic_atoms)
         call symbolic_atoms(ia)%restore_to_default()
         symbolic_atoms(ia)%potential%lmax = lmax
         do s = 1, 2
            do l = 0, lmax
               symbolic_atoms(ia)%potential%dele(l, s) = 0.70_rp + 0.06_rp*real(l, rp) + 0.025_rp*real(s, rp)
            end do
         end do
      end do

      allocate (lat%iz(nsite))
      lat%iz = 1
      allocate (en%ene(ne))
      do ie = 1, ne
         energies(ie) = cmplx(-0.83_rp + 0.29_rp*real(ie - 1, rp), 0.13_rp + 0.01_rp*real(ie, rp), rp)
         en%ene(ie) = real(energies(ie), rp)
      end do

      gobj%lattice => lat
      gobj%symbolic_atom => symbolic_atoms
      gobj%en => en
   end subroutine setup_green

   subroutine setup_screening(screen_in, screen_out)
      real(rp), intent(out) :: screen_in(0:, :), screen_out(0:, :)
      integer :: l, s

      do s = 1, 2
         do l = 0, lmax
            screen_in(l, s) = 0.12_rp + 0.031_rp*real(l, rp) + 0.017_rp*real(s, rp)
            screen_out(l, s) = -0.08_rp + 0.019_rp*real(l, rp) - 0.013_rp*real(s, rp)
         end do
      end do
   end subroutine setup_screening

   subroutine build_potential_matrices(energies, screen_in, screen_out, cshift, p_in, p_out)
      complex(rp), intent(in) :: energies(:)
      real(rp), intent(in) :: screen_in(0:, :), screen_out(0:, :)
      real(rp), intent(in) :: cshift
      complex(rp), intent(out) :: p_in(:, :, :), p_out(:, :, :)
      integer :: ie, l, m, s, idx
      real(rp) :: delta, center

      p_in = cmplx(0.0_rp, 0.0_rp, rp)
      p_out = cmplx(0.0_rp, 0.0_rp, rp)
      do ie = 1, size(energies)
         do s = 1, 2
            do l = 0, lmax
               delta = 0.73_rp + 0.04_rp*real(l, rp) + 0.02_rp*real(s, rp)
               center = cshift + 0.09_rp*real(l, rp) - 0.03_rp*real(s, rp)
               do m = 1, 2*l + 1
                  idx = l*l + m + norb_site*(s - 1)
                  p_in(idx, idx, ie) = (energies(ie) - cmplx(center, 0.0_rp, rp))/(delta*delta)
                  p_out(idx, idx, ie) = p_in(idx, idx, ie)/(1.0_rp + &
                     (screen_in(l, s) - screen_out(l, s))*p_in(idx, idx, ie))
               end do
            end do
         end do
      end do
   end subroutine build_potential_matrices

   subroutine build_deterministic_gf(gf)
      complex(rp), intent(out) :: gf(:, :, :)
      integer :: ie, i, j

      do ie = 1, size(gf, 3)
         do j = 1, size(gf, 2)
            do i = 1, size(gf, 1)
               gf(i, j, ie) = cmplx(0.017_rp*real(i + 2*j + ie, rp), &
                                    -0.013_rp*real(2*i - j + ie, rp), rp)
               if (i == j) gf(i, j, ie) = gf(i, j, ie) + cmplx(0.41_rp, 0.07_rp, rp)
            end do
         end do
      end do
   end subroutine build_deterministic_gf

   subroutine test_auxiliary_endpoint(gobj, physical_gf)
      type(green), intent(inout) :: gobj
      complex(rp), intent(in) :: physical_gf(:, :, :)
      complex(rp) :: auxiliary(size(physical_gf, 1), size(physical_gf, 2), size(physical_gf, 3))
      complex(rp) :: expected(size(physical_gf, 1), size(physical_gf, 2), size(physical_gf, 3))
      complex(rp) :: left(size(physical_gf, 1), size(physical_gf, 2))
      complex(rp) :: right(size(physical_gf, 1), size(physical_gf, 2))
      integer :: ie, l, m, s, idx

      left = cmplx(0.0_rp, 0.0_rp, rp)
      right = cmplx(0.0_rp, 0.0_rp, rp)
      do s = 1, 2
         do l = 0, lmax
            do m = 1, 2*l + 1
               idx = l*l + m + norb_site*(s - 1)
               left(idx, idx) = cmplx(gobj%symbolic_atom(1)%potential%dele(l, s), 0.0_rp, rp)
               right(idx, idx) = cmplx(gobj%symbolic_atom(2)%potential%dele(l, s), 0.0_rp, rp)
            end do
         end do
      end do
      do ie = 1, size(physical_gf, 3)
         expected(:, :, ie) = matmul(left, matmul(physical_gf(:, :, ie), right))
      end do
      call gobj%auxiliary_gij(physical_gf, auxiliary, 1, 2)
      call report('sqrt(Delta) endpoint scaling', maxval(abs(auxiliary - expected)), tol)
   end subroutine test_auxiliary_endpoint

   subroutine test_round_trip(gobj, screen_in, screen_out, p_in_i, p_out_i, p_in_j, p_out_j, &
                              gf_in, gf_out, gf_back, gf_wrong)
      type(green), intent(inout) :: gobj
      real(rp), intent(in) :: screen_in(0:, :), screen_out(0:, :)
      complex(rp), intent(in) :: p_in_i(:, :, :), p_out_i(:, :, :), p_in_j(:, :, :), p_out_j(:, :, :)
      complex(rp), intent(in) :: gf_in(:, :, :)
      complex(rp), intent(out) :: gf_out(:, :, :), gf_back(:, :, :), gf_wrong(:, :, :)
      complex(rp) :: ratio(size(gf_in, 1), size(gf_in, 2), size(gf_in, 3))
      complex(rp) :: missing_term(size(gf_in, 1), size(gf_in, 2), size(gf_in, 3))
      integer :: ie, i

      call gobj%transform_auxiliary_gij(p_in_i, p_out_i, p_in_j, p_out_j, gf_in, gf_out, &
                                        screen_in, screen_out, 1, 2)
      call gobj%transform_auxiliary_gij(p_out_i, p_in_i, p_out_j, p_in_j, gf_out, gf_back, &
                                        screen_out, screen_in, 1, 2)
      call report('screened offsite round trip', maxval(abs(gf_back - gf_in)), tol)

      do ie = 1, size(gf_in, 3)
         ratio(:, :, ie) = cmplx(0.0_rp, 0.0_rp, rp)
         do i = 1, size(gf_in, 1)
            ratio(i, i, ie) = p_in_i(i, i, ie)/p_out_i(i, i, ie)
         end do
         gf_wrong(:, :, ie) = matmul(ratio(:, :, ie), matmul(gf_in(:, :, ie), ratio(:, :, ie)))
      end do

      call gobj%transform_auxiliary_gij(p_in_i, p_out_i, p_in_i, p_out_i, gf_in, gf_out, &
                                        screen_in, screen_out, 1, 1)
      missing_term = gf_out - gf_wrong
      call report('onsite additive term is nonzero', maxval(abs(missing_term)), negative_floor, require_greater=.true.)
      call gobj%transform_auxiliary_gij(p_out_i, p_in_i, p_out_i, p_in_i, gf_out, gf_back, &
                                        screen_out, screen_in, 1, 1)
      call report('screened onsite round trip', maxval(abs(gf_back - gf_in)), tol)
   end subroutine test_round_trip

   subroutine test_bridge(gobj, screen_in, screen_out, p_in_i, p_out_i, p_in_j, p_out_j, energies)
      type(green), intent(inout) :: gobj
      real(rp), intent(in) :: screen_in(0:, :), screen_out(0:, :)
      complex(rp), intent(in) :: p_in_i(:, :, :), p_out_i(:, :, :), p_in_j(:, :, :), p_out_j(:, :, :)
      complex(rp), intent(in) :: energies(:)
      complex(rp) :: hk(nmat, nmat), gk(nmat, nmat), sigma(nmat, nmat)
      complex(rp) :: eigenvectors(nmat, nmat, nk), reciprocal_gf(nb, nb, ne), rs_gf(nb, nb, ne)
      complex(rp) :: reciprocal_aux(nb, nb, ne), rs_aux(nb, nb, ne)
      complex(rp) :: reciprocal_screened(nb, nb, ne), rs_screened(nb, nb, ne)
      complex(rp) :: d0(nmat, nmat), m0(nmat, nmat), seed1(nb, nb), seed2(nb, nb), seed3(nb, nb), seed4(nb, nb)
      real(rp) :: eigenvalues(nmat, nk), kfrac(3, nk), dr(3), theta, kdotr
      complex(rp) :: phase
      integer :: ik, ie

      d0 = cmplx(0.0_rp, 0.0_rp, rp)
      m0 = cmplx(0.0_rp, 0.0_rp, rp)
      do ie = 1, nb
         d0(ie, ie) = cmplx(-0.43_rp + 0.014_rp*real(ie, rp), 0.0_rp, rp)
         d0(nb + ie, nb + ie) = cmplx(0.27_rp - 0.011_rp*real(ie, rp), 0.0_rp, rp)
         m0(ie, nb + ie) = cmplx(0.08_rp + 0.003_rp*real(ie, rp), 0.012_rp*real(mod(ie, 3), rp), rp)
      end do
      do ik = 1, nk
         kfrac(:, ik) = 0.0_rp
         kfrac(1, ik) = real(ik - 1, rp)/real(nk, rp)
         theta = two_pi*kfrac(1, ik)
         hk = d0 + m0*exp(i_unit*theta) + conjg(transpose(m0))*exp(-i_unit*theta)
         call hermitian_eig(hk, eigenvalues(:, ik), eigenvectors(:, :, ik))
      end do
      dr = [ -1.0_rp, 0.0_rp, 0.0_rp ]
      call lehmann_pair_block(eigenvalues, eigenvectors, kfrac, energies, dr, 0, nb, nb, reciprocal_gf)

      sigma = cmplx(0.0_rp, 0.0_rp, rp)
      rs_gf = cmplx(0.0_rp, 0.0_rp, rp)
      do ie = 1, ne
         do ik = 1, nk
            theta = two_pi*kfrac(1, ik)
            hk = d0 + m0*exp(i_unit*theta) + conjg(transpose(m0))*exp(-i_unit*theta)
            call dyson_kspace_inverse(hk, energies(ie), sigma, gk)
            seed1 = 0.5_rp*(gk(1:nb, 1:nb) + gk(1:nb, nb + 1:nmat) + &
                            gk(nb + 1:nmat, 1:nb) + gk(nb + 1:nmat, nb + 1:nmat))
            seed2 = 0.5_rp*(gk(1:nb, 1:nb) - gk(1:nb, nb + 1:nmat) - &
                            gk(nb + 1:nmat, 1:nb) + gk(nb + 1:nmat, nb + 1:nmat))
            seed3 = 0.5_rp*(gk(1:nb, 1:nb) + i_unit*gk(1:nb, nb + 1:nmat) - &
                            i_unit*gk(nb + 1:nmat, 1:nb) + gk(nb + 1:nmat, nb + 1:nmat))
            seed4 = 0.5_rp*(gk(1:nb, 1:nb) - i_unit*gk(1:nb, nb + 1:nmat) + &
                            i_unit*gk(nb + 1:nmat, 1:nb) + gk(nb + 1:nmat, nb + 1:nmat))
            kdotr = two_pi*dot_product(kfrac(:, ik), dr)
            phase = cmplx(cos(kdotr), sin(kdotr), rp)
            rs_gf(:, :, ie) = rs_gf(:, :, ie) + phase*(0.5_rp*(seed1 - seed2 + (seed3 - seed4)/i_unit))
         end do
         rs_gf(:, :, ie) = rs_gf(:, :, ie)/real(nk, rp)
      end do
      call report('reciprocal Lehmann vs native-RS phase bridge', maxval(abs(reciprocal_gf - rs_gf)), bridge_tol)

      call gobj%auxiliary_gij(reciprocal_gf, reciprocal_aux, 1, 2)
      call gobj%auxiliary_gij(rs_gf, rs_aux, 1, 2)
      call gobj%transform_auxiliary_gij(p_in_i, p_out_i, p_in_j, p_out_j, reciprocal_aux, reciprocal_screened, &
                                        screen_in, screen_out, 1, 2)
      call gobj%transform_auxiliary_gij(p_in_i, p_out_i, p_in_j, p_out_j, rs_aux, rs_screened, &
                                        screen_in, screen_out, 1, 2)
      call report('endpoint and screening transformed bridge', maxval(abs(reciprocal_screened - rs_screened)), bridge_tol)
   end subroutine test_bridge

   subroutine test_high_energy(gobj, energies, screen_in, screen_out)
      type(green), intent(inout) :: gobj
      complex(rp), intent(in) :: energies(:)
      real(rp), intent(in) :: screen_in(0:, :), screen_out(0:, :)
      complex(rp) :: h(nmat, nmat), gk(nmat, nmat), sigma(nmat, nmat)
      complex(rp) :: on_mid(nb, nb), on_high(nb, nb), off_high(nb, nb)
      complex(rp) :: aux_mid(nb, nb, ne), aux_high(nb, nb, ne), aux_off(nb, nb, ne)
      complex(rp) :: p_in(nb, nb, 1), p_out(nb, nb, 1)
      complex(rp) :: path_in(nb, nb, 1), path_out(nb, nb, 1), path_expected(nb, nb, 1)
      complex(rp) :: expected(nb, nb), identity(nb, nb)
      complex(rp) :: z_mid, z_high, phase
      real(rp) :: err_raw_mid, err_raw_high, err_aux_mid, err_aux_high, err_off
      integer :: ie, i, local_orbital, orbital_l, spin

      sigma = cmplx(0.0_rp, 0.0_rp, rp)
      h = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, nb
         h(i, i) = cmplx(-0.20_rp + 0.01_rp*real(i, rp), 0.0_rp, rp)
         h(nb + i, nb + i) = cmplx(0.31_rp - 0.01_rp*real(i, rp), 0.0_rp, rp)
         h(i, nb + i) = cmplx(0.17_rp, -0.02_rp, rp)
         h(nb + i, i) = conjg(h(i, nb + i))
      end do
      z_mid = cmplx(80.0_rp, 0.4_rp, rp)
      z_high = cmplx(1000.0_rp, 0.4_rp, rp)
      call dyson_kspace_inverse(h, z_mid, sigma, gk)
      on_mid = gk(1:nb, 1:nb)
      call dyson_kspace_inverse(h, z_high, sigma, gk)
      on_high = gk(1:nb, 1:nb)
      off_high = gk(1:nb, nb + 1:nmat)

      identity = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, nb
         identity(i, i) = cmplx(1.0_rp, 0.0_rp, rp)
      end do
      expected = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, nb
         spin = 1 + (i - 1)/norb_site
         local_orbital = 1 + mod(i - 1, norb_site)
         orbital_l = merge(0, 1, local_orbital == 1)
         expected(i, i) = cmplx(gobj%symbolic_atom(1)%potential%dele(orbital_l, spin), 0.0_rp, rp)**2
      end do
      call gobj%auxiliary_gij(reshape(on_mid, [nb, nb, 1]), aux_mid, 1, 1)
      call gobj%auxiliary_gij(reshape(on_high, [nb, nb, 1]), aux_high, 1, 1)
      call gobj%auxiliary_gij(reshape(off_high, [nb, nb, 1]), aux_off, 1, 2)
      err_raw_mid = maxval(abs(z_mid*on_mid - identity))
      err_raw_high = maxval(abs(z_high*on_high - identity))
      err_aux_mid = maxval(abs(z_mid*aux_mid(:, :, 1) - expected))
      err_aux_high = maxval(abs(z_high*aux_high(:, :, 1) - expected))
      err_off = maxval(abs(z_high*aux_off(:, :, 1)))
      call report('raw onsite high-energy normalization improves', err_raw_high, err_raw_mid, require_less=.true.)
      call report('sqrt(Delta) onsite high-energy normalization improves', err_aux_high, err_aux_mid, require_less=.true.)
      call report('sqrt(Delta) offsite high-energy decay', err_off, 5.0e-3_rp)

      call build_potential_matrices([z_high], screen_in, screen_out, 0.11_rp, p_in, p_out)
      path_in = cmplx(0.0_rp, 0.0_rp, rp)
      path_expected = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, nb
         path_in(i, i, 1) = 1.0_rp/p_in(i, i, 1)
         path_expected(i, i, 1) = 1.0_rp/p_out(i, i, 1)
      end do
      call gobj%transform_auxiliary_gij(p_in, p_out, p_in, p_out, path_in, path_out, &
                                        screen_in, screen_out, 1, 1)
      call report('screened high-energy onsite identity', maxval(abs(path_out - path_expected)), tol)
   end subroutine test_high_energy

   subroutine hermitian_eig(matrix_in, eigenvalues, eigenvectors)
      complex(rp), intent(in) :: matrix_in(:, :)
      real(rp), intent(out) :: eigenvalues(:)
      complex(rp), intent(out) :: eigenvectors(:, :)
      complex(rp) :: work_query(1)
      complex(rp), allocatable :: work(:)
      real(rp) :: rwork(max(1, 3*size(eigenvalues) - 2))
      integer :: info, lwork

      eigenvectors = matrix_in
      call zheev('V', 'U', size(eigenvalues), eigenvectors, size(eigenvalues), eigenvalues, work_query, -1, rwork, info)
      lwork = max(1, nint(real(work_query(1), rp)))
      allocate (work(lwork))
      call zheev('V', 'U', size(eigenvalues), eigenvectors, size(eigenvalues), eigenvalues, work, lwork, rwork, info)
      if (info /= 0) error stop 'LR-REP-00 Hermitian eigensolve failed'
      deallocate (work)
   end subroutine hermitian_eig

   subroutine report(label, value, limit, require_greater, require_less)
      character(len=*), intent(in) :: label
      real(rp), intent(in) :: value, limit
      logical, intent(in), optional :: require_greater, require_less
      logical :: greater, less

      greater = .false.
      less = .false.
      if (present(require_greater)) greater = require_greater
      if (present(require_less)) less = require_less
      if (greater) then
         write (*, '(a,es16.8,a,es16.8)') trim(label)//' = ', value, ' > ', limit
         if (value <= limit) failed = .true.
      else if (less) then
         write (*, '(a,es16.8,a,es16.8)') trim(label)//' = ', value, ' < ', limit
         if (value >= limit) failed = .true.
      else
         write (*, '(a,es16.8,a,es16.8)') trim(label)//' = ', value, ' <= ', limit
         if (value > limit) failed = .true.
      end if
   end subroutine report

end program test_lr_rs_gf_representation
