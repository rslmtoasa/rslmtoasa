!> k-space assembly of the GBT (gbt_single_q) Hamiltonian, relation R1 (see gbt_kspace_mod). It reads the structure
!> constants, neighbour vectors and potential parameters and never the real-space blocks `ee`.
submodule (reciprocal_mod) reciprocal_gbt
   use gbt_kspace_mod, only: gbt_kspace_h1, gbt_kspace_h2, gbt_kspace_hcc, gbt_kspace_horth
   use gbt_structure_mod, only: gbt_reference_frame
   use magnetic_representation_mod, only: gbt_single_q
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
      complex(rp), allocatable :: sdotblock(:, :, :), hcc(:, :)
      real(rp), allocatable :: ccd(:, :, :)
      real(rp) :: ccd_full(norb, 0:2), lam
      complex(rp), allocatable :: obar(:, :, :), enu(:, :, :), wunit(:, :, :), czero3(:, :, :)
      complex(rp) :: blk(nb, nb)
      complex(rp), allocatable :: cmat(:, :)
      real(rp), allocatable :: cvec(:), dvec(:), gvec(:)
      real(rp) :: wow
      integer :: lmax_o, l, o, sp, iv
      real(rp) :: q_dir(3), theta, phi, mom_ref(3)
      logical :: second, gbt_rep

      gbt_rep = trim(this%hamiltonian%magnetic_representation) == gbt_single_q
      if (.not. gbt_rep .and. .not. this%use_orthogonal) &
         call g_logger%fatal("k-space GBT/orthogonal assembler called for neither gbt_single_q nor kspace_ham_order='orthogonal'.", __FILE__, __LINE__)
      if (lower(trim(this%lattice%strux_backend)) /= 'strux_lib') &
         call g_logger%fatal("k-space GBT/orthogonal requires strux_backend='strux_lib'.", __FILE__, __LINE__)
      if (this%control%nsp == 2 .or. this%control%nsp == 4) &
         call g_logger%fatal('k-space GBT/orthogonal does not support SOC (nsp = 2 or 4).', __FILE__, __LINE__)
      if (this%hamiltonian%ccor_2c) then   ! Hcc is added after the order step; the orthogonal H^gamma has no CCOR formulation here
         if (this%use_orthogonal) &
            call g_logger%fatal("k-space GBT/orthogonal: kspace_ham_order='orthogonal' does not support ccor_2c.", __FILE__, __LINE__)
         if (trim(this%hamiltonian%ccor_vmt_mode) == 'pair_surface') &
            call g_logger%fatal("k-space GBT/orthogonal: the channel-wise pair average of ccor_vmt_mode='pair_surface' is not defined for rotated local frames; use 'surface_scalar'.", __FILE__, __LINE__)
         call this%hamiltonian%validate_ccor_inputs()
         lam = this%hamiltonian%ccor_vmt_scalar() - this%hamiltonian%ccor_elin
      end if
      if (this%hamiltonian%hubbard_u_general_check .or. this%hamiltonian%hubbard_test_operator_check) &
         call g_logger%fatal('k-space GBT/orthogonal does not support onsite U.', __FILE__, __LINE__)
      if (this%hamiltonian%hubbard_v_check) &
         call g_logger%fatal('k-space GBT/orthogonal does not support intersite V.', __FILE__, __LINE__)
      if (this%hamiltonian%local_axis) call g_logger%fatal('k-space GBT/orthogonal does not support local_axis.', __FILE__, __LINE__)
      if (this%control%constraints_enable) then   ! the field is stored per site and applied per type: they must coincide
         if (this%use_orthogonal) &
            call g_logger%fatal("k-space GBT/orthogonal: kspace_ham_order='orthogonal' does not support constraining fields.", __FILE__, __LINE__)
         if (this%lattice%nrec /= this%lattice%ntype .or. this%lattice%nbulk /= 0) &
            call g_logger%fatal('k-space GBT/orthogonal constraints need one type per site and nbulk = 0.', __FILE__, __LINE__)
      end if

      second = this%use_second_order .and. this%hamiltonian%hoh
      if (second .and. .not. (allocated(this%hamiltonian%obarm) .and. allocated(this%hamiltonian%enim))) &
         call g_logger%fatal('k-space GBT/orthogonal second order needs obarm and enim from the Hamiltonian build.', __FILE__, __LINE__)
      if (this%use_orthogonal .and. this%hamiltonian%hoh) &
         call g_logger%fatal("kspace_ham_order='orthogonal' replaces h - h o h; hoh must be .false.", __FILE__, __LINE__)
      nrec = this%lattice%nrec
      allocate (bond_a(nrec*this%lattice%nn_max), bond_b(nrec*this%lattice%nn_max), bond_d(3, nrec*this%lattice%nn_max), &
                sblock(norb, norb, nrec*this%lattice%nn_max), frame(2, 2, nrec), w(norb, 2, nrec), c(norb, 2, nrec), &
                h(nb*nrec, nb*nrec), h2(nb*nrec, nb*nrec), obar(nb, nb, nrec), enu(nb, nb, nrec), &
                cvec(nb*nrec), dvec(nb*nrec), gvec(nb*nrec), cmat(nb*nrec, nb*nrec), wunit(norb, 2, nrec), czero3(norb, 2, nrec), &
                sdotblock(norb, norb, nrec*this%lattice%nn_max), hcc(nb*nrec, nb*nrec), ccd(norb, 0:1, nrec))
      wunit = cmplx(1.0_rp, 0.0_rp, rp)
      czero3 = cmplx(0.0_rp, 0.0_rp, rp)
      nbond = 0
      do isite = 1, nrec
         ntype_i = this%lattice%ib(isite)
         ia = this%lattice%atlist(ntype_i)
         nr = this%lattice%nn(ia, 1)
         it = this%lattice%iz(ia)
         associate (pot => this%lattice%symbolic_atoms(it)%potential)
            if (size(pot%wx0) < norb .or. size(pot%cx0) < norb) &
               call g_logger%fatal('k-space GBT/orthogonal: potential factors smaller than the orbital basis.', __FILE__, __LINE__)
            mom_ref = pot%mom
            if (allocated(this%hamiltonian%gbt_ref_mom)) mom_ref = this%hamiltonian%gbt_ref_mom(:, it)
            if (gbt_rep) then
               call gbt_site_frame(mom_ref, this%hamiltonian%theta_ss, theta, phi)
            else   ! 'orthogonal' outside GBT: nsp = 1 keeps the lab spin channels, nsp = 3 follows the current moment direction
               call mom_frame(pot%mom, this%control%nsp == 1, theta, phi)
            end if
            call gbt_reference_frame(theta, phi, frame(:, :, isite))
            if (this%control%constraints_enable) then   ! the target is the local axis, whatever the site orientation
               if (maxval(abs(this%lattice%symbolic_atoms(it)%constraint_target - [0.0_rp, 0.0_rp, 1.0_rp])) > 1.0e-10_rp) &
                  call g_logger%fatal('k-space GBT/orthogonal needs constraints_mom_ref = 0, 0, 1 (local frame) for every site.', &
                                      __FILE__, __LINE__)
            end if
            w(:, 1, isite) = pot%wx0(1:norb) + pot%wx1(1:norb)
            w(:, 2, isite) = pot%wx0(1:norb) - pot%wx1(1:norb)
            c(:, 1, isite) = pot%cx0(1:norb) + pot%cx1(1:norb)
            c(:, 2, isite) = pot%cx0(1:norb) - pot%cx1(1:norb)
            if (this%hamiltonian%ccor_2c) then   ! sign +1 for every site: the frame carries the orientation, as for w and c
               call this%hamiltonian%build_ccor_coefficients(ia, it, ccd_full)
               ccd(:, :, isite) = ccd_full(:, 0:1)
            end if
            if (this%use_orthogonal) then   ! exact ASA in the gamma representation: local-frame diagonals of C, sqrt(Delta), gamma - alpha
               lmax_o = int(sqrt(real(norb - 1, rp)))
               if (.not. allocated(pot%screening_alpha)) &
                  call g_logger%fatal("kspace_ham_order='orthogonal' needs the strux screening alpha.", __FILE__, __LINE__)
               if (size(pot%screening_alpha) /= lmax_o + 1 .or. size(pot%srdel, 1) < lmax_o + 1 .or. &
                   size(pot%qpar, 1) < lmax_o + 1 .or. pot%ws_r <= 0.0_rp) &
                  call g_logger%fatal("kspace_ham_order='orthogonal': potential parameters do not match the orbital basis.", __FILE__, __LINE__)
               wow = this%lattice%wav*ang2au/pot%ws_r   ! as predls (symbolic_atom.f90:297-306), index l + 1 in the exponents
               do sp = 1, 2
                  do o = 1, norb
                     l = int(sqrt(real(o - 1, rp)))
                     iv = (isite - 1)*nb + (sp - 1)*norb + o
                     cvec(iv) = pot%c(l, sp) + pot%vmad
                     dvec(iv) = pot%srdel(l, sp)*wow**(0.5_rp - real(l + 1, rp))
                     gvec(iv) = pot%qpar(l, sp)*wow**(1 - 2*(l + 1)) - pot%screening_alpha(l)
                  end do
               end do
            end if
            if (second) then
               ! As the collinear route at second order: onsite cex in h, E_nu = cx - cex added by build_enim. obarm and
               ! enim hold the spin blocks for the sign of mom(3); swapping them gives the sign +1 blocks.
               c(:, 1, isite) = pot%cex0(1:norb) + pot%cex1(1:norb)
               c(:, 2, isite) = pot%cex0(1:norb) - pot%cex1(1:norb)
               obar(:, :, isite) = this%hamiltonian%obarm(:, :, it)
               enu(:, :, isite) = this%hamiltonian%enim(:, :, it)
               if (mom_ref(3) < 0.0_rp) then
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
            if (this%hamiltonian%ccor_2c) then   ! as raw_sdot in build_gbt_bulkham
               do j = 1, norb
                  do i = 1, norb
                     sdotblock(i, j, nbond) = this%hamiltonian%normalize_ccor_sdot( &
                        this%lattice%sdot(j, i, ineigh, this%lattice%num(ia)))
                  end do
               end do
            end if
         end do
      end do
      cmat = cmplx(0.0_rp, 0.0_rp, rp)
      if (this%use_orthogonal) then
         do iv = 1, nb*nrec
            cmat(iv, iv) = cvec(iv)
         end do
      end if
      q_dir = 0.0_rp
      if (gbt_rep) then
         do i = 1, 3
            q_dir(i) = dot_product(this%lattice%a(:, i), this%hamiltonian%q_ss)   ! q_ss is Cartesian in 2 pi/alat
         end do
      end if

      do ik = 1, size(k_points, 2)
         if (this%use_orthogonal) then   ! h1 with w = 1, c = 0 is S^alpha_rot
            call gbt_kspace_h1(norb, nrec, nbond, bond_a(1:nbond), bond_b(1:nbond), bond_d(:, 1:nbond), sblock(:, :, 1:nbond), &
                               frame, wunit, czero3, k_points(:, ik), q_dir, h2)
            call gbt_kspace_horth(nb*nrec, h2, cmat, dvec, gvec, h)
         else
            call gbt_kspace_h1(norb, nrec, nbond, bond_a(1:nbond), bond_b(1:nbond), bond_d(:, 1:nbond), sblock(:, :, 1:nbond), &
                               frame, w, c, k_points(:, ik), q_dir, h)
         end if
         if (this%control%constraints_enable) then   ! -> +B.sigma onsite in the local spin basis, as the real-space route
            do isite = 1, nrec
               blk = cmplx(0.0_rp, 0.0_rp, rp)
               call this%hamiltonian%add_constraining_field_spin_block(this%lattice%iz(this%lattice%atlist(this%lattice%ib(isite))), blk)
               i0 = (isite - 1)*nb
               h(i0 + 1:i0 + nb, i0 + 1:i0 + nb) = h(i0 + 1:i0 + nb, i0 + 1:i0 + nb) + blk
            end do
         end if
         call to_spherical(h)
         if (.not. gbt_rep .and. this%control%nsp /= 1) call to_lab(h)   ! outside GBT the eigenvectors are in the lab spin basis
         if (second) then   ! h, obar and enu are all in the spherical-harmonic basis here
            call gbt_kspace_h2(norb, nrec, h, obar, enu, h2)
            h = h2
         end if
         if (this%hamiltonian%ccor_2c) then   ! after the order step: Hcc never enters h - h o h
            call gbt_kspace_hcc(norb, nrec, nbond, bond_a(1:nbond), bond_b(1:nbond), bond_d(:, 1:nbond), sblock(:, :, 1:nbond), &
                                sdotblock(:, :, 1:nbond), frame, w, ccd, lam, k_points(:, ik), q_dir, hcc)
            call to_spherical(hcc)
            h = h + hcc
         end if
         workspace%h(:, :, ik) = h
      end do
   contains

      !> Cartesian-to-spherical orbital basis of every (site, spin) block, as build_gbt_bulkham.
      subroutine to_spherical(m)
         complex(rp), intent(inout) :: m(:, :)
         integer :: isite, jsite, s, t, i0, j0

         do isite = 1, nrec
            do jsite = 1, nrec
               do s = 0, 1
                  do t = 0, 1
                     i0 = (isite - 1)*nb + s*norb
                     j0 = (jsite - 1)*nb + t*norb
                     call hcpx(m(i0 + 1:i0 + norb, j0 + 1:j0 + norb), 'cart2sph')
                  end do
               end do
            end do
         end do
      end subroutine to_spherical

      !> Local-frame spin basis to lab spin basis, H_lab(sigma, sigma') = sum_{s,t} F_a(sigma, s) H(s, t) conj(F_b(sigma', t)). At q = 0 the
      !> frames are the rotations that take z to the moment direction, so this is the same H as C0 + C1 sigma.m with a rotated S^gamma.
      subroutine to_lab(m)
         complex(rp), intent(inout) :: m(:, :)
         complex(rp) :: blk(norb, norb, 2, 2), rot(norb, norb, 2, 2)
         integer :: isite, jsite, s, t, sg, sp, i0, j0

         do isite = 1, nrec
            do jsite = 1, nrec
               i0 = (isite - 1)*nb
               j0 = (jsite - 1)*nb
               do s = 1, 2
                  do t = 1, 2
                     blk(:, :, s, t) = m(i0 + (s - 1)*norb + 1:i0 + s*norb, j0 + (t - 1)*norb + 1:j0 + t*norb)
                  end do
               end do
               rot = cmplx(0.0_rp, 0.0_rp, rp)
               do sg = 1, 2
                  do sp = 1, 2
                     do s = 1, 2
                        do t = 1, 2
                           rot(:, :, sg, sp) = rot(:, :, sg, sp) + frame(sg, s, isite)*blk(:, :, s, t)*conjg(frame(sp, t, jsite))
                        end do
                     end do
                  end do
               end do
               do sg = 1, 2
                  do sp = 1, 2
                     m(i0 + (sg - 1)*norb + 1:i0 + sg*norb, j0 + (sp - 1)*norb + 1:j0 + sp*norb) = rot(:, :, sg, sp)
                  end do
               end do
            end do
         end do
      end subroutine to_lab

   end subroutine reciprocal_assembler_assemble_gbt

   !> Frame angles of one site for 'orthogonal' outside GBT: the direction of the current moment (theta from z, phi azimuth), or the
   !> identity frame when the lab spin channels are kept (nsp = 1). The moment direction is the majority axis, as for ham0m_nc.
   subroutine mom_frame(mom, keep_lab, theta, phi)
      real(rp), intent(in) :: mom(3)
      logical, intent(in) :: keep_lab
      real(rp), intent(out) :: theta, phi
      real(rp) :: norm

      theta = 0.0_rp
      phi = 0.0_rp
      norm = sqrt(sum(mom**2))
      if (keep_lab .or. norm <= 1.0e-12_rp) return
      theta = acos(max(-1.0_rp, min(1.0_rp, mom(3)/norm)))
      if (abs(mom(1)) + abs(mom(2)) > 1.0e-12_rp) phi = atan2(mom(2), mom(1))
   end subroutine mom_frame

   !> Frame angles of one site from its reference moment (D2): +z gives (theta_ss, 0), -z gives (pi - theta_ss, pi), so the
   !> antiparallel site stays at m_b = -m_a for every cone angle. A moment off the z axis is not supported yet.
   subroutine gbt_site_frame(mom, theta_ss, theta, phi)
      real(rp), intent(in) :: mom(3), theta_ss
      real(rp), intent(out) :: theta, phi

      if (abs(mom(1)) + abs(mom(2)) > 1.0e-12_rp) &
         call g_logger%fatal('k-space GBT/orthogonal needs potential%mom along +-z.', __FILE__, __LINE__)
      theta = theta_ss
      phi = 0.0_rp
      if (mom(3) < 0.0_rp) then
         theta = pi - theta
         phi = pi
      end if
   end subroutine gbt_site_frame

end submodule reciprocal_gbt
