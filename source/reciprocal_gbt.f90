!> k-space assembly of the GBT (gbt_single_q) Hamiltonian, relation R1 (see gbt_kspace_mod). It reads the structure
!> constants, neighbour vectors and potential parameters and never the real-space blocks `ee`.
submodule (reciprocal_mod) reciprocal_gbt
   use gbt_kspace_mod, only: gbt_kspace_h1, gbt_kspace_h2
   use gbt_structure_mod, only: gbt_reference_frame
   implicit none

contains

   !> Fill workspace%h(:, :, 1:nk) for direct k_points. k - q/2 and k + q/2 come from the k given: nothing is folded
   !> separately. A site with potential%mom(3) < 0 takes theta = pi - theta_ss and phi = pi with W and C as for a
   !> majority site, so the frame alone carries the orientation (R1; this avoids finding F1 of the GBT restart spec).
   module subroutine reciprocal_assembler_assemble_gbt(this, k_points, workspace)
      class(reciprocal_assembler), intent(in) :: this
      real(rp), intent(in) :: k_points(:, :)
      class(reciprocal_workspace), intent(inout) :: workspace

      integer :: nrec, nbond, isite, jsite, ntype_i, ia, ja, nr, ineigh, it, i, j, ik, s, t, i0, j0
      integer, allocatable :: bond_a(:), bond_b(:)
      real(rp), allocatable :: bond_d(:, :)
      complex(rp), allocatable :: sblock(:, :, :), frame(:, :, :), w(:, :, :), c(:, :, :), h(:, :), h2(:, :)
      complex(rp), allocatable :: obar(:, :, :), enu(:, :, :)
      real(rp) :: q_dir(3), theta, phi
      logical :: second

      if (lower(trim(this%lattice%strux_backend)) /= 'strux_lib') &
         call g_logger%fatal("gbt_single_q k-space requires strux_backend='strux_lib'.", __FILE__, __LINE__)
      if (this%control%nsp == 2 .or. this%control%nsp == 4) &
         call g_logger%fatal('gbt_single_q k-space does not support SOC (nsp = 2 or 4).', __FILE__, __LINE__)
      if (this%hamiltonian%ccor_2c) call g_logger%fatal('gbt_single_q k-space does not support ccor_2c.', __FILE__, __LINE__)
      if (this%hamiltonian%hubbard_u_general_check .or. this%hamiltonian%hubbard_test_operator_check) &
         call g_logger%fatal('gbt_single_q k-space does not support onsite U.', __FILE__, __LINE__)
      if (this%hamiltonian%hubbard_v_check) &
         call g_logger%fatal('gbt_single_q k-space does not support intersite V.', __FILE__, __LINE__)
      if (this%control%constraints_enable) &
         call g_logger%fatal('gbt_single_q k-space does not support constraining fields.', __FILE__, __LINE__)
      if (this%hamiltonian%local_axis) call g_logger%fatal('gbt_single_q k-space does not support local_axis.', __FILE__, __LINE__)
      if (allocated(this%hamiltonian%theta_ss_sublattice)) &
         call g_logger%fatal('gbt_single_q k-space supports one cone angle for all sites only.', __FILE__, __LINE__)

      second = this%use_second_order .and. this%hamiltonian%hoh
      if (second .and. .not. (allocated(this%hamiltonian%obarm) .and. allocated(this%hamiltonian%enim))) &
         call g_logger%fatal('gbt_single_q k-space second order needs obarm and enim from the Hamiltonian build.', __FILE__, __LINE__)
      nrec = this%lattice%nrec
      allocate (bond_a(nrec*this%lattice%nn_max), bond_b(nrec*this%lattice%nn_max), bond_d(3, nrec*this%lattice%nn_max), &
                sblock(norb, norb, nrec*this%lattice%nn_max), frame(2, 2, nrec), w(norb, 2, nrec), c(norb, 2, nrec), &
                h(nb*nrec, nb*nrec), h2(nb*nrec, nb*nrec), obar(nb, nb, nrec), enu(nb, nb, nrec))
      nbond = 0
      do isite = 1, nrec
         ntype_i = this%lattice%ib(isite)
         ia = this%lattice%atlist(ntype_i)
         nr = this%lattice%nn(ia, 1)
         it = this%lattice%iz(ia)
         associate (pot => this%lattice%symbolic_atoms(it)%potential)
            if (abs(pot%mom(1)) + abs(pot%mom(2)) > 1.0e-12_rp) &
               call g_logger%fatal('gbt_single_q k-space needs potential%mom along +-z.', __FILE__, __LINE__)
            if (size(pot%wx0) < norb .or. size(pot%cx0) < norb) &
               call g_logger%fatal('gbt_single_q k-space: potential factors smaller than the orbital basis.', __FILE__, __LINE__)
            theta = this%hamiltonian%theta_ss
            phi = 0.0_rp
            if (pot%mom(3) < 0.0_rp) then
               theta = pi - theta
               phi = pi
            end if
            call gbt_reference_frame(theta, phi, frame(:, :, isite))
            w(:, 1, isite) = pot%wx0(1:norb) + pot%wx1(1:norb)
            w(:, 2, isite) = pot%wx0(1:norb) - pot%wx1(1:norb)
            c(:, 1, isite) = pot%cx0(1:norb) + pot%cx1(1:norb)
            c(:, 2, isite) = pot%cx0(1:norb) - pot%cx1(1:norb)
            if (second) then
               ! As the collinear route at second order: onsite cex in h, E_nu = cx - cex added by build_enim. obarm and
               ! enim hold the spin blocks for the sign of mom(3); swapping them gives the sign +1 blocks.
               c(:, 1, isite) = pot%cex0(1:norb) + pot%cex1(1:norb)
               c(:, 2, isite) = pot%cex0(1:norb) - pot%cex1(1:norb)
               obar(:, :, isite) = this%hamiltonian%obarm(:, :, it)
               enu(:, :, isite) = this%hamiltonian%enim(:, :, it)
               if (pot%mom(3) < 0.0_rp) then
                  obar(1:norb, 1:norb, isite) = this%hamiltonian%obarm(norb + 1:nb, norb + 1:nb, it)
                  obar(norb + 1:nb, norb + 1:nb, isite) = this%hamiltonian%obarm(1:norb, 1:norb, it)
                  enu(1:norb, 1:norb, isite) = this%hamiltonian%enim(norb + 1:nb, norb + 1:nb, it)
                  enu(norb + 1:nb, norb + 1:nb, isite) = this%hamiltonian%enim(1:norb, 1:norb, it)
               end if
            end if
         end associate
         do ineigh = 1, nr
            if (ineigh == 1) then
               jsite = isite
            else
               ja = this%lattice%nn(ia, ineigh)
               if (ja < 1 .or. ja > this%lattice%kk) cycle
               jsite = this%lattice%iz(ja)
               if (jsite < 1 .or. jsite > nrec) cycle
            end if
            nbond = nbond + 1
            bond_a(nbond) = isite
            bond_b(nbond) = jsite
            bond_d(:, nbond) = 0.0_rp
            if (ineigh > 1) bond_d(:, nbond) = this%ham_vec_type_direct(:, ineigh, ntype_i)
            do j = 1, norb
               do i = 1, norb
                  sblock(i, j, nbond) = this%lattice%sbar(j, i, ineigh, this%lattice%num(ia))   ! as raw_orb in build_gbt_bulkham
               end do
            end do
         end do
      end do
      do i = 1, 3
         q_dir(i) = dot_product(this%lattice%a(:, i), this%hamiltonian%q_ss)   ! q_ss is Cartesian in 2 pi/alat
      end do

      do ik = 1, size(k_points, 2)
         call gbt_kspace_h1(norb, nrec, nbond, bond_a(1:nbond), bond_b(1:nbond), bond_d(:, 1:nbond), sblock(:, :, 1:nbond), &
                            frame, w, c, k_points(:, ik), q_dir, h)
         do isite = 1, nrec
            do jsite = 1, nrec
               do s = 0, 1
                  do t = 0, 1
                     i0 = (isite - 1)*nb + s*norb
                     j0 = (jsite - 1)*nb + t*norb
                     call hcpx(h(i0 + 1:i0 + norb, j0 + 1:j0 + norb), 'cart2sph')   ! as build_gbt_bulkham
                  end do
               end do
            end do
         end do
         if (second) then   ! h, obar and enu are all in the spherical-harmonic basis here
            call gbt_kspace_h2(norb, nrec, h, obar, enu, h2)
            h = h2
         end if
         workspace%h(:, :, ik) = h
      end do
   end subroutine reciprocal_assembler_assemble_gbt

end submodule reciprocal_gbt
