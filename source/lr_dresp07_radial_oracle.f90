!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Numerical oracles for the DRESP-07R local-field variational audit.
!
!> The local spherical operator projection and the radial inverse are kept
!> separate deliberately.  The first is an exact Frobenius projection onto
!> site/l projectors.  The second is a real DGELSS minimum-norm solve of the
!> complete radial map; no shell-by-shell division is permitted.
!------------------------------------------------------------------------------
module lr_dresp07_radial_oracle_mod

   use precision_mod, only: rp
   use radial_ground_state_mod, only: radial_simpson_weight
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, lmto_orbital_l
   implicit none
   private

   public :: dresp07_build_radial_map
   public :: dresp07_project_local_spherical_operator
   public :: dresp07_solve_real_least_squares

   interface
      subroutine dgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info)
         import :: rp
         integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
         real(rp), intent(inout) :: a(lda, *), b(ldb, *), work(*)
         real(rp), intent(out) :: s(*)
         real(rp), intent(in) :: rcond
         integer, intent(out) :: rank, info
      end subroutine dgelss
   end interface

contains

   !> Build the complete site/l by site/r radial map used by the accepted
   !> large-component augmentation convention.
   subroutine dresp07_build_radial_map(radial_bases, radial_map)
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      real(rp), intent(out) :: radial_map(:, :)
      integer :: site, l, ir, row, col, nsite, lmax, npoint

      nsite = size(radial_bases)
      if (nsite < 1) error stop 'DRESP-07R radial map: empty site set'
      lmax = radial_bases(1)%lmax
      npoint = radial_bases(1)%npoint
      if (any(shape(radial_map) /= [nsite*(lmax + 1), nsite*npoint])) then
         error stop 'DRESP-07R radial map: shape mismatch'
      end if
      if (any([(radial_bases(site)%lmax, site=1,nsite)] /= lmax) .or. &
          any([(radial_bases(site)%npoint, site=1,nsite)] /= npoint)) then
         error stop 'DRESP-07R radial map: radial bases are not conforming'
      end if

      radial_map = 0.0_rp
      do site = 1, nsite
         do l = 0, lmax
            row = (site - 1)*(lmax + 1) + l + 1
            do ir = 2, npoint
               col = (site - 1)*npoint + ir
               radial_map(row, col) = radial_simpson_weight(ir, npoint)*radial_bases(site)%mesh_a* &
                  (radial_bases(site)%rofi(ir) + radial_bases(site)%mesh_b)* &
                  radial_bases(site)%phi_large(ir, l + 1, 1)**2
            end do
         end do
      end do
   end subroutine dresp07_build_radial_map

   !> Project a coefficient-space field onto the exact site/l spherical
   !> operator subspace.  The returned coefficients are the Frobenius-optimal
   !> projector coefficients and are m-degenerate by construction.
   subroutine dresp07_project_local_spherical_operator(target, norb_site, nsite, lmax, best, coefficients, residual_norm, &
                                                        relative_residual, offdiag_norm, within_l_anisotropy_norm, &
                                                        cross_l_norm, intersite_norm)
      complex(rp), intent(in) :: target(:, :)
      integer, intent(in) :: norb_site, nsite, lmax
      complex(rp), intent(out) :: best(:, :), coefficients(:, :)
      real(rp), intent(out) :: residual_norm, relative_residual, offdiag_norm
      real(rp), intent(out) :: within_l_anisotropy_norm, cross_l_norm, intersite_norm
      integer :: site, other_site, iorb, jorb, i, j, l, nmode
      complex(rp) :: coefficient
      real(rp) :: target_norm, offdiag_sum, anisotropy_sum, cross_l_sum, intersite_sum, mean_value

      if (any(shape(target) /= [norb_site*nsite, norb_site*nsite]) .or. &
          any(shape(best) /= shape(target)) .or. any(shape(coefficients) /= [nsite, lmax + 1])) then
         error stop 'DRESP-07R local projection: shape mismatch'
      end if
      best = cmplx(0.0_rp, 0.0_rp, rp)
      coefficients = cmplx(0.0_rp, 0.0_rp, rp)
      offdiag_sum = 0.0_rp
      anisotropy_sum = 0.0_rp
      cross_l_sum = 0.0_rp
      intersite_sum = 0.0_rp

      do site = 1, nsite
         do l = 0, lmax
            coefficient = cmplx(0.0_rp, 0.0_rp, rp)
            nmode = 0
            do iorb = 1, norb_site
               if (lmto_orbital_l(iorb) == l) then
                  i = (site - 1)*norb_site + iorb
                  coefficient = coefficient + target(i, i)
                  nmode = nmode + 1
               end if
            end do
            coefficient = coefficient/real(max(1, nmode), rp)
            coefficients(site, l + 1) = coefficient
            do iorb = 1, norb_site
               if (lmto_orbital_l(iorb) == l) then
                  i = (site - 1)*norb_site + iorb
                  best(i, i) = coefficient
               end if
            end do
         end do
      end do

      do site = 1, nsite
         do iorb = 1, norb_site
            i = (site - 1)*norb_site + iorb
            do other_site = 1, nsite
               do jorb = 1, norb_site
                  j = (other_site - 1)*norb_site + jorb
                  if (site /= other_site) then
                     intersite_sum = intersite_sum + abs(target(i, j))**2
                  else if (iorb == jorb) then
                     l = lmto_orbital_l(iorb)
                     mean_value = real(coefficients(site, l + 1), rp)
                     anisotropy_sum = anisotropy_sum + (real(target(i, j), rp) - mean_value)**2 + aimag(target(i, j))**2
                  else if (lmto_orbital_l(iorb) == lmto_orbital_l(jorb)) then
                     offdiag_sum = offdiag_sum + abs(target(i, j))**2
                  else
                     cross_l_sum = cross_l_sum + abs(target(i, j))**2
                  end if
               end do
            end do
         end do
      end do
      residual_norm = sqrt(sum(abs(target - best)**2))
      target_norm = sqrt(sum(abs(target)**2))
      relative_residual = residual_norm/max(target_norm, tiny(1.0_rp))
      offdiag_norm = sqrt(offdiag_sum)
      within_l_anisotropy_norm = sqrt(anisotropy_sum)
      cross_l_norm = sqrt(cross_l_sum)
      intersite_norm = sqrt(intersite_sum)
   end subroutine dresp07_project_local_spherical_operator

   !> Solve A x = b by DGELSS.  The caller owns the design matrix and target;
   !> this routine reports the actual singular spectrum, rank, condition, and
   !> minimum-norm residual of the complete coupled solve.
   subroutine dresp07_solve_real_least_squares(design, target, solution, rank_out, singular_values, condition, residual_norm, &
                                               relative_residual, solution_norm)
      real(rp), intent(in) :: design(:, :), target(:)
      real(rp), intent(out) :: solution(:), singular_values(:)
      integer, intent(out) :: rank_out
      real(rp), intent(out) :: condition, residual_norm, relative_residual, solution_norm
      real(rp), allocatable :: a_work(:, :), b_work(:, :), work(:)
      real(rp) :: work_query(1), target_norm, smax, smin
      integer :: m, n, ldb, lwork, info, i
      real(rp), parameter :: rcond = -1.0_rp

      m = size(design, 1)
      n = size(design, 2)
      ldb = max(m, n)
      if (size(target) /= m .or. size(solution) /= n .or. size(singular_values) /= min(m, n)) then
         error stop 'DRESP-07R DGELSS: shape mismatch'
      end if
      allocate(a_work(m, n), b_work(ldb, 1))
      a_work = design
      b_work = 0.0_rp
      b_work(1:m, 1) = target
      call dgelss(m, n, 1, a_work, m, b_work, ldb, singular_values, rcond, rank_out, work_query, -1, info)
      if (info /= 0) error stop 'DRESP-07R DGELSS: workspace query failed'
      lwork = max(1, nint(work_query(1)))
      allocate(work(lwork))
      a_work = design
      b_work = 0.0_rp
      b_work(1:m, 1) = target
      call dgelss(m, n, 1, a_work, m, b_work, ldb, singular_values, rcond, rank_out, work, lwork, info)
      if (info /= 0) error stop 'DRESP-07R DGELSS: solve failed'
      solution = b_work(1:n, 1)
      residual_norm = sqrt(sum((matmul(design, solution) - target)**2))
      target_norm = sqrt(sum(target**2))
      relative_residual = residual_norm/max(target_norm, tiny(1.0_rp))
      solution_norm = sqrt(sum(solution**2))
      if (rank_out > 0) then
         smax = maxval(singular_values(1:rank_out))
         smin = minval(singular_values(1:rank_out))
         condition = smax/max(smin, tiny(1.0_rp))
      else
         condition = huge(1.0_rp)
      end if
      deallocate(a_work, b_work, work)
   end subroutine dresp07_solve_real_least_squares

end module lr_dresp07_radial_oracle_mod
