!> Pair-list kernels on a bcc cluster built from integer combinations of the primitive vectors.
!> Expected counts are closed-form: shells of bcc (|r|^2 in a^2) 0 (1), 3/4 (8), 1 (6), 2 (12), 11/4 (24), 3 (8).
program test_pair_list
   use precision_mod, only: rp
   use pair_list_mod, only: pairs_within_radius, site_cart_rotations, pair_stars
   implicit none

   integer, parameter :: nmax = 4
   real(rp), parameter :: a(3, 3) = reshape([-0.5_rp, 0.5_rp, 0.5_rp, 0.5_rp, -0.5_rp, 0.5_rp, 0.5_rp, 0.5_rp, -0.5_rp], [3, 3])
   integer, parameter :: shell_count(6) = [1, 8, 6, 12, 24, 8]
   real(rp), parameter :: shell_r(6) = [0.0_rp, 0.8660254037844386_rp, 1.0_rp, 1.4142135623730951_rp, 1.6583123951777_rp, &
                                       1.7320508075688772_rp]
   real(rp), allocatable :: cr(:, :), rot48(:, :, :), rot8(:, :, :), rc(:, :, :), tr(:, :)
   real(rp) :: v(3)
   integer, allocatable :: jl(:), rep(:), types(:), rf(:, :, :)
   integer :: n1, n2, n3, n, seed, i, k, s, nbad, nrep, perm(3, 6), sg, p, q, n48, n8
   logical :: failed = .false.
   real(rp) :: d, ainv(3, 3), r(3, 3)


   ! cluster
   allocate (cr(3, (2*nmax + 1)**3))
   n = 0
   do n3 = -nmax, nmax
      do n2 = -nmax, nmax
         do n1 = -nmax, nmax
            n = n + 1
            cr(:, n) = n1*a(:, 1) + n2*a(:, 2) + n3*a(:, 3)
         end do
      end do
   end do
   seed = (n + 1)/2
   call require(maxval(abs(cr(:, seed))) == 0.0_rp, 'seed is the origin')
   allocate (types(n))
   types = 1

   ! catches: wrong radius comparison (< instead of <=), missing tolerance at the shell edge, wrong shell content
   call pairs_within_radius(cr, seed, 1.75_rp, jl)
   call require(size(jl) == sum(shell_count), 'radius 1.75: 59 pairs')
   call pairs_within_radius(cr, seed, sqrt(3.0_rp), jl)
   call require(size(jl) == sum(shell_count), 'radius sqrt(3) a on the shell: 59 pairs')
   call pairs_within_radius(cr, seed, 1.73_rp, jl)
   call require(size(jl) == sum(shell_count(1:5)), 'radius 1.73: 51 pairs')
   call pairs_within_radius(cr, seed, sqrt(11.0_rp)/2.0_rp, jl)
   call require(size(jl) == sum(shell_count(1:5)), 'radius sqrt(11)/2 a on the shell: 51 pairs')

   ! catches: self pair not first, order not by distance, ties not by atom index, wrong shell sizes
   call pairs_within_radius(cr, seed, 1.75_rp, jl)
   call require(jl(1) == seed, 'self pair first')
   k = 0
   do s = 1, 6
      do i = 1, shell_count(s)
         k = k + 1
         d = norm2(cr(:, jl(k)))
         if (abs(d - shell_r(s)) > 1.0e-9_rp) call require(.false., 'shell order and size')
         if (i > 1) call require(jl(k) > jl(k - 1), 'ties by atom index')
      end do
   end do

   ! 48 cubic operations as signed permutation matrices
   perm = reshape([1, 2, 3, 1, 3, 2, 2, 1, 3, 2, 3, 1, 3, 1, 2, 3, 2, 1], [3, 6])
   allocate (rot48(3, 3, 48), rf(3, 3, 48), tr(3, 48))
   tr = 0.0_rp
   ainv = inv3(a)
   n48 = 0
   do p = 1, 6
      do sg = 0, 7
         n48 = n48 + 1
         rot48(:, :, n48) = 0.0_rp
         do q = 1, 3
            rot48(q, perm(q, p), n48) = merge(-1.0_rp, 1.0_rp, btest(sg, q - 1))
         end do
         rf(:, :, n48) = nint(matmul(ainv, matmul(rot48(:, :, n48), a)))
      end do
   end do

   ! catches: operations not mapped to the site-invariant set, wrong Cartesian conversion
   call site_cart_rotations(a, rf, tr, [0.0_rp, 0.0_rp, 0.0_rp], rc)
   call require(size(rc, 3) == 48, 'site symmetry at the origin: 48 operations')
   r = 0.0_rp
   do k = 1, 48
      r = max(r, abs(rc(:, :, k) - rot48(:, :, k)))
   end do
   call require(maxval(r) < 1.0e-12_rp, 'Cartesian rotations equal the input')
   v = matmul(ainv, [0.25_rp, 0.0_rp, 0.0_rp])
   call site_cart_rotations(a, rf, tr, v, rc)
   call require(size(rc, 3) == 8, 'site symmetry on the x axis: 8 operations')

   ! catches: stars that do not cover a shell, wrong multiplicity, representative not the first member
   call pair_stars(cr, types, seed, jl, rot48, rep, nbad)
   call require(nbad == 0, 'stars: all images found')
   nrep = 0
   do k = 1, size(jl)
      if (rep(k) == k) nrep = nrep + 1
   end do
   call require(nrep == 6, 'cubic group: 6 representatives')
   k = 0
   do s = 1, 6
      call require(count(rep == k + 1) == shell_count(s), 'star multiplicity of a shell')
      k = k + shell_count(s)
   end do

   ! C4v about z: 8-shell splits 4 + 4, 6-shell 4 + 1 + 1
   n8 = 0
   do k = 1, 48
      if (rot48(3, 3, k) == 1.0_rp) n8 = n8 + 1
   end do
   allocate (rot8(3, 3, n8))
   n8 = 0
   do k = 1, 48
      if (rot48(3, 3, k) == 1.0_rp) then
         n8 = n8 + 1
         rot8(:, :, n8) = rot48(:, :, k)
      end if
   end do
   call require(n8 == 8, 'C4v: 8 operations')
   call pair_stars(cr, types, seed, jl, rot8, rep, nbad)
   call require(nbad == 0, 'C4v stars: all images found')
   call require(nstars(rep, 2, 9) == 2 .and. maxstar(rep, 2, 9) == 4 .and. minstar(rep, 2, 9) == 4, &
                '8-shell under C4v: 2 stars of 4')
   call require(nstars(rep, 10, 15) == 3 .and. maxstar(rep, 10, 15) == 4 .and. minstar(rep, 10, 15) == 1, &
                '6-shell under C4v: stars 4 + 1 + 1')

   ! catches: a type that is ignored (C4v maps the type-2 neighbour onto type-1 atoms, so images must be reported missing)
   types = 1
   types(jl(2)) = 2
   call pair_stars(cr, types, seed, jl, rot8, rep, nbad)
   call require(nbad > 0, 'a rotation that changes the type is not a symmetry')

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains
   integer function nstars(rp_, lo, hi)
      integer, intent(in) :: rp_(:), lo, hi
      integer :: kk
      nstars = 0
      do kk = lo, hi
         if (rp_(kk) == kk) nstars = nstars + 1
      end do
   end function nstars

   integer function maxstar(rp_, lo, hi)
      integer, intent(in) :: rp_(:), lo, hi
      integer :: kk
      maxstar = 0
      do kk = lo, hi
         if (rp_(kk) == kk) maxstar = max(maxstar, count(rp_(lo:hi) == kk))
      end do
   end function maxstar

   integer function minstar(rp_, lo, hi)
      integer, intent(in) :: rp_(:), lo, hi
      integer :: kk
      minstar = huge(1)
      do kk = lo, hi
         if (rp_(kk) == kk) minstar = min(minstar, count(rp_(lo:hi) == kk))
      end do
   end function minstar

   subroutine require(ok, label)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: label
      if (.not. ok) then
         failed = .true.
         write (*, '(2a)') 'FAIL: ', label
      end if
   end subroutine require

   function inv3(m) result(w)
      real(rp), intent(in) :: m(3, 3)
      real(rp) :: w(3, 3), dt
      dt = m(1, 1)*(m(2, 2)*m(3, 3) - m(2, 3)*m(3, 2)) - m(1, 2)*(m(2, 1)*m(3, 3) - m(2, 3)*m(3, 1)) &
           + m(1, 3)*(m(2, 1)*m(3, 2) - m(2, 2)*m(3, 1))
      w(1, 1) = m(2, 2)*m(3, 3) - m(2, 3)*m(3, 2); w(1, 2) = m(1, 3)*m(3, 2) - m(1, 2)*m(3, 3)
      w(1, 3) = m(1, 2)*m(2, 3) - m(1, 3)*m(2, 2); w(2, 1) = m(2, 3)*m(3, 1) - m(2, 1)*m(3, 3)
      w(2, 2) = m(1, 1)*m(3, 3) - m(1, 3)*m(3, 1); w(2, 3) = m(1, 3)*m(2, 1) - m(1, 1)*m(2, 3)
      w(3, 1) = m(2, 1)*m(3, 2) - m(2, 2)*m(3, 1); w(3, 2) = m(1, 2)*m(3, 1) - m(1, 1)*m(3, 2)
      w(3, 3) = m(1, 1)*m(2, 2) - m(1, 2)*m(2, 1)
      w = w/dt
   end function inv3
end program test_pair_list
