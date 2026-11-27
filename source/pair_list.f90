!> Pair lists of the exchange post-processing generated from a radius around the recursion atoms.
!> @details The kernels take plain arrays: pairs_within_radius selects and orders the neighbours of one seed,
!>          site_cart_rotations turns space-group operations into the Cartesian rotations that leave a seed site
!>          invariant, pair_stars groups the pairs of one seed into symmetry stars. generate_exchange_pairs fills
!>          lattice%ijpair (and pair_full, pair_rep with pair_unique) from the lattice keys pair_radius, pair_unique.
module pair_list_mod
   use precision_mod, only: rp
   use logger_mod, only: g_logger
   use string_mod, only: int2str
   use math_mod, only: inverse_3x3
   use lattice_mod, only: lattice
   use symmetry_mod, only: symmetry
   use mpi_mod, only: rank
   implicit none
   private

   !> Distance tolerance (units of alat) for the radius, the distance ties and the equivalence of positions.
   real(rp), parameter :: pair_tol = 1.0e-6_rp

   public :: pairs_within_radius, site_cart_rotations, pair_stars, generate_exchange_pairs

contains

   !> Cluster atoms j with |r_j - r_seed| <= radius (+ pair_tol), the seed first, then by distance (ties within pair_tol by
   !> atom index).
   subroutine pairs_within_radius(cr, seed, radius, jlist)
      real(rp), intent(in) :: cr(:, :), radius
      integer, intent(in) :: seed
      integer, allocatable, intent(out) :: jlist(:)
      real(rp), allocatable :: d(:)
      real(rp) :: dk
      integer :: j, k, n, jk

      d = norm2(cr - spread(cr(:, seed), 2, size(cr, 2)), dim=1)
      jlist = pack([(j, j=1, size(cr, 2))], d <= radius + pair_tol)
      n = size(jlist)
      do j = 2, n
         jk = jlist(j)
         dk = d(jk)
         k = j - 1
         do while (k >= 1)
            if (.not. (d(jlist(k)) > dk + pair_tol .or. (abs(d(jlist(k)) - dk) <= pair_tol .and. jlist(k) > jk))) exit
            jlist(k + 1) = jlist(k)
            k = k - 1
         end do
         jlist(k + 1) = jk
      end do
   end subroutine pairs_within_radius

   !> Cartesian rotations of the space-group operations (f' = rot_frac f + trans, fractional coordinates of the primitive
   !> vectors a(:, i)) that map the site at fractional position frac_seed onto itself modulo a lattice vector.
   subroutine site_cart_rotations(a, rot_frac, trans, frac_seed, rot_cart)
      real(rp), intent(in) :: a(3, 3), trans(:, :), frac_seed(3)
      integer, intent(in) :: rot_frac(:, :, :)
      real(rp), allocatable, intent(out) :: rot_cart(:, :, :)
      real(rp) :: d(3), ainv(3, 3)
      integer :: k, n
      logical, allocatable :: keep(:)

      ainv = inverse_3x3(a)
      allocate (keep(size(rot_frac, 3)))
      do k = 1, size(rot_frac, 3)
         d = matmul(real(rot_frac(:, :, k), rp), frac_seed) + trans(:, k) - frac_seed
         keep(k) = maxval(abs(d - nint(d))) < 1.0e-4_rp
      end do
      allocate (rot_cart(3, 3, count(keep)))
      n = 0
      do k = 1, size(rot_frac, 3)
         if (.not. keep(k)) cycle
         n = n + 1
         rot_cart(:, :, n) = matmul(a, matmul(real(rot_frac(:, :, k), rp), ainv))
      end do
   end subroutine site_cart_rotations

   !> Symmetry stars of the pairs (seed, jlist(k)): rep(k) is the first position in jlist of the star of k. Two pairs are in
   !> one star when a rotation of rot_cart maps r_j - r_seed onto r_j' - r_seed and types(j) == types(j'). irot(k) is the
   !> rotation that maps the vector of rep(k) onto that of k (0 for a representative). nbad counts the (rotation, pair)
   !> images that fall on no atom of jlist (zero for a complete list and a true symmetry).
   subroutine pair_stars(cr, types, seed, jlist, rot_cart, rep, irot, nbad)
      real(rp), intent(in) :: cr(:, :), rot_cart(:, :, :)
      integer, intent(in) :: types(:), seed, jlist(:)
      integer, allocatable, intent(out) :: rep(:), irot(:)
      integer, intent(out) :: nbad
      real(rp) :: v(3), w(3)
      integer :: k, l, m, r

      allocate (rep(size(jlist)), irot(size(jlist)))
      rep = 0
      irot = 0
      nbad = 0
      do k = 1, size(jlist)
         if (rep(k) /= 0) cycle
         rep(k) = k
         v = cr(:, jlist(k)) - cr(:, seed)
         do r = 1, size(rot_cart, 3)
            w = matmul(rot_cart(:, :, r), v)
            m = 0
            do l = 1, size(jlist)
               if (maxval(abs(cr(:, jlist(l)) - cr(:, seed) - w)) <= pair_tol .and. types(jlist(l)) == types(jlist(k))) then
                  m = l
                  exit
               end if
            end do
            if (m == 0) then
               nbad = nbad + 1
            else if (rep(m) == 0) then
               rep(m) = k
               irot(m) = r
            end if
         end do
      end do
   end subroutine pair_stars

   !> Fill lattice%ijpair and njij from pair_radius for every recursion atom of the calculation. With pair_unique the
   !> ijpair rows are the star representatives; the full list is kept in pair_full and pair_rep.
   subroutine generate_exchange_pairs(lat)
      type(lattice), target, intent(inout) :: lat
      type(symmetry) :: sym
      integer, allocatable :: jl(:), rep(:), irot(:), seeds(:), full(:, :), fullrep(:), rot_frac(:, :, :), flips(:), newidx(:)
      logical, allocatable :: isrep(:)
      real(rp), allocatable :: trans(:, :), rot_cart(:, :, :), frac(:, :), moments(:, :), rots(:, :, :), rs(:, :, :)
      real(rp) :: r_cluster
      real(rp), parameter :: identity(3, 3) = reshape([1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp], [3, 3])
      integer :: s, k, n, nfull, nops, nbad, nrep, nsp

      if (lat%pair_radius <= 0.0_rp) then
         if (lat%pair_unique) call g_logger%fatal('pair_unique needs pair_radius > 0.', __FILE__, __LINE__)
         return
      end if
      if (lat%njij /= 0 .or. lat%njijk /= 0) &
         call g_logger%fatal('pair_radius generates the pair list: do not set njij, ijpair or njijk together with it.', &
                             __FILE__, __LINE__)
      nsp = lat%control%nsp
      if (lat%pair_unique .and. nsp /= 1 .and. nsp /= 2 .and. nsp /= 4) &
         call g_logger%fatal('pair_unique needs nsp = 1, 2 or 4.', __FILE__, __LINE__)
      if (lat%pair_unique .and. trim(lat%control%calctype) /= 'B') &
         call g_logger%fatal('pair_unique needs calctype = ''B'' (the crystal space group acts on the bulk cluster).', &
                             __FILE__, __LINE__)

      seeds = lat%irec(1:lat%nrec)
      ! Cluster radius: the largest distance of a cluster atom from atom 1 (r2 is the Hamiltonian cutoff, not the cluster)
      r_cluster = maxval(norm2(lat%cr(:, 1:lat%kk) - spread(lat%cr(:, 1), 2, lat%kk), dim=1))
      do s = 1, size(seeds)
         if (norm2(lat%cr(:, seeds(s)) - lat%cr(:, 1)) + lat%pair_radius > r_cluster + pair_tol) &
            call g_logger%fatal('pair_radius = ' // trim(adjustl(real_str(lat%pair_radius))) // ' of seed atom ' // &
                                trim(int2str(seeds(s))) // ' reaches the cluster edge (cluster radius ' // &
                                trim(adjustl(real_str(r_cluster))) // ' from atom 1); increase rc or reduce pair_radius.', &
                                __FILE__, __LINE__)
      end do

      ! Full list, seed by seed
      allocate (full(0, 2), fullrep(0), rots(3, 3, 0))
      nfull = 0
      allocate (isrep(0))
      if (lat%pair_unique) then
#ifdef USE_SPGLIB
         n = lat%nrec
         call sym%initialize(lat)
         allocate (frac(3, n), moments(3, n))
         do k = 1, n
            frac(:, k) = matmul(inverse_3x3(lat%a), lat%cr(:, k))
            moments(:, k) = site_axis(lat%symbolic_atoms(lat%nbulk + k))
         end do
         ! nsp = 2, 4: magnetic group of the moments as axial vectors (time-reversal-combined operations included);
         ! nsp = 1: operations that map every site onto a site with the same type and moment
         nops = sym%spglib%get_magnetic_symmetry(frac, lat%iz(1:n), moments, nsp == 2 .or. nsp == 4, rot_frac, trans, flips)
         if (nops == 0) call g_logger%fatal('pair_unique: spglib returned no symmetry operations.', __FILE__, __LINE__)
#else
         call g_logger%fatal('pair_unique needs spglib, which is not compiled in.', __FILE__, __LINE__)
#endif
      end if
      do s = 1, size(seeds)
         call pairs_within_radius(lat%cr, seeds(s), lat%pair_radius, jl)
         n = size(jl)
         if (lat%pair_unique) then
            call site_cart_rotations(lat%a, rot_frac, trans, matmul(inverse_3x3(lat%a), lat%cr(:, seeds(s))), rot_cart)
            call pair_stars(lat%cr, lat%iz, seeds(s), jl, rot_cart, rep, irot, nbad)
            if (nbad /= 0) call g_logger%fatal('pair_unique: ' // trim(int2str(nbad)) // ' symmetry images of the pairs of atom ' // &
                                               trim(int2str(seeds(s))) // ' are not in the cluster.', __FILE__, __LINE__)
            allocate (rs(3, 3, n))
            do k = 1, n
               rs(:, :, k) = identity
               if (irot(k) /= 0) rs(:, :, k) = rot_cart(:, :, irot(k))
            end do
            rots = reshape([rots, rs], [3, 3, nfull + n])
            deallocate (rs)
         else
            rep = [(k, k=1, n)]
         end if
         full = reshape([full(:, 1), spread(seeds(s), 1, n), full(:, 2), jl], [nfull + n, 2])
         fullrep = [fullrep, nfull + rep]
         isrep = [isrep, rep == [(k, k=1, n)]]
         nfull = nfull + n
      end do

      nrep = count(isrep)
      if (lat%pair_unique) then
         ! Row of ijpair of every full pair
         allocate (newidx(nfull))
         newidx = 0
         n = 0
         do k = 1, nfull
            if (isrep(k)) then
               n = n + 1
               newidx(k) = n
            end if
         end do
         if (allocated(lat%ijpair)) deallocate (lat%ijpair)
         allocate (lat%ijpair(nrep, 2), lat%pair_rep(nfull))
         lat%ijpair = full(pack([(k, k=1, nfull)], isrep), :)
         lat%pair_rep = newidx(fullrep)
         lat%pair_full = full
         lat%pair_rot = rots
         lat%njij = nrep
         if (rank == 0) then
            do k = 1, nfull
               if (.not. isrep(k)) cycle
               call g_logger%info('pair_unique: pair ' // trim(int2str(full(k, 1))) // ' ' // trim(int2str(full(k, 2))) // &
                                  ' multiplicity ' // trim(int2str(count(lat%pair_rep == newidx(k)))), __FILE__, __LINE__)
            end do
            call g_logger%info('pair_radius: ' // trim(int2str(nfull)) // ' pairs, ' // trim(int2str(nrep)) // &
                               ' symmetry-unique (nops = ' // trim(int2str(nops)) // ')', __FILE__, __LINE__)
         end if
      else
         if (allocated(lat%ijpair)) deallocate (lat%ijpair)
         lat%ijpair = full
         lat%njij = nfull
         if (rank == 0) call g_logger%info('pair_radius: ' // trim(int2str(nfull)) // ' pairs', __FILE__, __LINE__)
      end if
   end subroutine generate_exchange_pairs

   !> Unit moment direction of a site (constraint target, else the potential moment), zero for a nonmagnetic site.
   function site_axis(atom) result(axis)
      use symbolic_atom_mod, only: symbolic_atom
      type(symbolic_atom), intent(in) :: atom
      real(rp) :: axis(3), nrm

      axis = atom%potential%mom
      if (norm2(atom%constraint_target) > 1.0e-12_rp) axis = atom%constraint_target
      nrm = norm2(axis)
      axis = merge(axis/nrm, [0.0_rp, 0.0_rp, 1.0_rp], nrm > tiny(1.0_rp))
      if (atom%force_nonmagnetic) axis = 0.0_rp
   end function site_axis

   function real_str(x) result(str)
      real(rp), intent(in) :: x
      character(len=24) :: str
      write (str, '(f0.6)') x
   end function real_str

end module pair_list_mod
