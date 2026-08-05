submodule(hamiltonian_mod) hamiltonian_ccor
   use magnetic_representation_mod, only: explicit_texture, gbt_single_q, select_endpoint_moments
   use gbt_structure_mod, only: gbt_lift_orbital_block, gbt_contract_collinear

contains

   !> @brief Add two-centre combined-correction blocks to the bulk Hamiltonian.
   !> @details Builds CCOR contributions for each bulk atom type and neighbor
   !>          shell, using scalar or noncollinear pair-block helpers as required
   !>          by the magnetic state.
   !> @param[inout] this Hamiltonian object; fills eecc.
   !> @note validate_ccor_inputs must accept the current CCOR mode before use.
   module subroutine build_ccor_bulk(this)
      class(hamiltonian), intent(inout) :: this
      integer :: ntype, ia, ino, nr, m, jj, it, jt
      complex(rp), dimension(nb, nb) :: hcc

         this%eecc(:, :, :, :) = czero
         if (.not. this%ccor_2c) return
         call validate_ccor_inputs(this)

      do ntype = 1, this%charge%lattice%ntype
         ia = this%charge%lattice%atlist(ntype)
         ino = this%charge%lattice%num(ia)
         it = this%charge%lattice%iz(ia)
         nr = this%charge%lattice%nn(ia, 1)
         do m = 1, nr
            if (m == 1) then
               jj = ia
            else
               jj = this%charge%lattice%nn(ia, m)
            end if
            if (jj <= 0) cycle
            jt = this%charge%lattice%iz(jj)
            call build_ccor_pair_block_noncollinear(this, ia, jj, it, jt, ino, m, hcc)
            this%eecc(:, :, m, ntype) = hcc(:, :)
         end do
         end do

         if (maxval(abs(this%eecc)) <= tiny(1.0_rp)) then
            if (this%ccor_strict) then
               call g_logger%fatal('ccor_2c built zero bulk Hcc; check lattice%sdot, VMTZ/vrmax, and ccor_elin.', __FILE__, __LINE__)
            else if (rank == 0) then
               call g_logger%warning('ccor_2c built zero bulk Hcc; k-space/real-space results will be unchanged.', __FILE__, __LINE__)
            end if
         end if
         if (this%ccor_debug) call log_ccor_debug(this, this%eecc, 'bulk')
      end subroutine build_ccor_bulk

   !> @brief Add two-centre combined-correction blocks to local Hamiltonian regions.
   !> @details Builds CCOR contributions for local/surface/impurity neighbor
   !>          blocks in the same layout as hall, preserving local atom indexing.
   !> @param[inout] this Hamiltonian object; fills hallcc.
   module subroutine build_ccor_local(this)
      class(hamiltonian), intent(inout) :: this
      integer :: nlim, ino, nr, m, jj, it, jt
      complex(rp), dimension(nb, nb) :: hcc

      this%hallcc(:, :, :, :) = czero
      if (.not. this%ccor_2c) return
      call validate_ccor_inputs(this)

      do nlim = 1, this%charge%lattice%nmax
         ino = this%charge%lattice%num(nlim)
         it = this%charge%lattice%iz(nlim)
         nr = this%charge%lattice%nn(nlim, 1)
         do m = 1, nr
            if (m == 1) then
               jj = nlim
            else
               jj = this%charge%lattice%nn(nlim, m)
            end if
            if (jj <= 0) cycle
            jt = this%charge%lattice%iz(jj)
            call build_ccor_pair_block_noncollinear(this, nlim, jj, it, jt, ino, m, hcc)
            this%hallcc(:, :, m, nlim) = hcc(:, :)
         end do
         end do

         if (maxval(abs(this%hallcc)) <= tiny(1.0_rp)) then
            if (this%ccor_strict) then
               call g_logger%fatal('ccor_2c built zero local Hcc; check lattice%sdot, VMTZ/vrmax, and ccor_elin.', __FILE__, __LINE__)
            else if (rank == 0) then
               call g_logger%warning('ccor_2c built zero local Hcc; local results will be unchanged.', __FILE__, __LINE__)
            end if
         end if
         if (this%ccor_debug) call log_ccor_debug(this, this%hallcc, 'local')
      end subroutine build_ccor_local

   !> @brief Validate the current CCOR namelist and lattice state.
   !> @details Checks that two-centre combined-correction parameters are meaningful
   !>          for the active calculation and emits warnings or fatal diagnostics
   !>          according to ccor_strict.
   !> @param[in] this Hamiltonian object with CCOR flags and lattice state.
   module subroutine validate_ccor_inputs(this)
      class(hamiltonian), intent(in) :: this

      if (trim(this%magnetic_representation) == gbt_single_q .and. &
          (.not. allocated(this%charge%lattice%sdot) .or. &
           .not. this%charge%lattice%strux_want_sdot)) then
         call g_logger%fatal('gbt_single_q with CCOR requires generated Sdot and strux_want_sdot=.true.', &
                             __FILE__, __LINE__)
      end if
      if (.not. allocated(this%charge%lattice%sdot)) then
         if (this%ccor_strict) then
            call g_logger%fatal('ccor_2c requires lattice%sdot; set strux_want_sdot=.true. with strux_lib.', __FILE__, __LINE__)
         else
            call g_logger%warning('ccor_2c requested but lattice%sdot is not allocated; Hcc will be zero.', __FILE__, __LINE__)
         end if
      end if
      if (.not. this%charge%lattice%strux_want_sdot) then
         if (this%ccor_strict) then
            call g_logger%fatal('ccor_2c requires strux_want_sdot=.true.; refusing to continue in strict mode.', __FILE__, __LINE__)
         else
            call g_logger%warning('ccor_2c requires strux_want_sdot=.true.; check that sdot was generated.', __FILE__, __LINE__)
         end if
      end if
      if (trim(lower(this%charge%lattice%screening)) /= 'sigma' .and. &
          trim(lower(this%charge%lattice%screening)) /= 'fitted') then
         if (this%ccor_strict) then
            call g_logger%fatal('ccor_2c requires screened alpha/adot from lattice screening=sigma or fitted.', __FILE__, __LINE__)
         else
            call g_logger%warning('ccor_2c with lattice screening other than sigma/fitted may lack valid adot.', __FILE__, __LINE__)
         end if
      end if
   end subroutine validate_ccor_inputs

   !> Select endpoint moments for ordinary noncollinear/texture CCOR.
   !> GBT is intentionally excluded because S and Sdot must be linked before
   !> either primitive factor is contracted with endpoint potential factors.
   module subroutine ccor_select_endpoint_moments(this, ia, ja, it, jt, mom_i, mom_j)
      class(hamiltonian), intent(in) :: this
      integer, intent(in) :: ia, ja, it, jt
      real(rp), dimension(3), intent(out) :: mom_i, mom_j
      logical :: valid

      if (trim(this%magnetic_representation) == gbt_single_q) then
         call g_logger%fatal('gbt_single_q must not enter the contracted CCOR pair builder.', __FILE__, __LINE__)
      end if
      if (trim(this%magnetic_representation) == explicit_texture) then
         call select_endpoint_moments(this%magnetic_representation, ia, ja, &
            this%charge%lattice%symbolic_atoms(it)%potential%mom(:), &
            this%charge%lattice%symbolic_atoms(jt)%potential%mom(:), &
            this%texture_moments, mom_i, mom_j, valid)
         if (.not. valid) then
            call g_logger%fatal('explicit_texture CCOR endpoint moment lookup failed.', __FILE__, __LINE__)
         end if
      else
         mom_i = this%charge%lattice%symbolic_atoms(it)%potential%mom(:)
         mom_j = this%charge%lattice%symbolic_atoms(jt)%potential%mom(:)
      end if
   end subroutine ccor_select_endpoint_moments

   !> @brief Build one scalar-relativistic CCOR pair block.
   !> @details Computes the two-centre combined-correction hopping contribution
   !>          between atom sites ia/ja and types it/jt for neighbor index m.
   !> @param[in] this Hamiltonian object containing potentials and CCOR settings.
   !> @param[in] ia Site index of the source atom.
   !> @param[in] ja Site index of the neighbor atom.
   !> @param[in] it Atom type of ia.
   !> @param[in] jt Atom type of ja.
   !> @param[in] ino Bravais/type index used by the target Hamiltonian layout.
   !> @param[in] m Neighbor-shell index.
   !> @param[out] hcc CCOR block in spin-orbital basis.
   module subroutine build_ccor_pair_block_scalar(this, ia, ja, it, jt, ino, m, hcc)
      class(hamiltonian), intent(in) :: this
      integer, intent(in) :: ia, ja, it, jt, ino, m
      complex(rp), dimension(nb, nb), intent(out) :: hcc

      call build_ccor_pair_block_noncollinear(this, ia, ja, it, jt, ino, m, hcc)
   end subroutine build_ccor_pair_block_scalar

   !> Build one GBT CCOR block from primitive directed S and Sdot factors.
   !> Both factors receive the already-computed common endpoint link before
   !> contraction.  The completed D/Ddot/Hcc objects are never rephased.
   module subroutine build_ccor_pair_block_gbt(this, ia, ja, it, jt, m, link, &
                                                sign_i, sign_j, s_block, sdot_block, hcc)
      class(hamiltonian), intent(in) :: this
      integer, intent(in) :: ia, ja, it, jt, m
      complex(rp), dimension(2, 2), intent(in) :: link
      real(rp), intent(in) :: sign_i, sign_j
      complex(rp), dimension(norb, norb), intent(in) :: s_block, sdot_block
      complex(rp), dimension(nb, nb), intent(out) :: hcc

      complex(rp), dimension(nb, nb) :: linked_s, linked_sdot, dmat, ddotmat
      real(rp), dimension(norb, 0:2) :: ccd_i, ccd_j
      real(rp), dimension(2) :: lambda_pair, lambda_i, lambda_j
      real(rp) :: lambda0, lambda1
      complex(rp) :: wi
      integer :: ilm, jlm, ispin, jspin, irow, jcol

      call build_ccor_coefficients(this, ia, it, ccd_i)
      call build_ccor_coefficients(this, ja, jt, ccd_j)

      ! Primitive directed factors: S_ij G_ij and Sdot_ij G_ij.
      call gbt_lift_orbital_block(s_block, link, linked_s)
      call gbt_lift_orbital_block(sdot_block, link, linked_sdot)
      call gbt_contract_collinear(linked_s, &
         this%lattice%symbolic_atoms(it)%potential%wx0(1:norb), &
         this%lattice%symbolic_atoms(it)%potential%wx1(1:norb), &
         this%lattice%symbolic_atoms(jt)%potential%wx0(1:norb), &
         this%lattice%symbolic_atoms(jt)%potential%wx1(1:norb), sign_i, sign_j, dmat)
      call gbt_contract_collinear(linked_sdot, &
         this%lattice%symbolic_atoms(it)%potential%wx0(1:norb), &
         this%lattice%symbolic_atoms(it)%potential%wx1(1:norb), &
         this%lattice%symbolic_atoms(jt)%potential%wx0(1:norb), &
         this%lattice%symbolic_atoms(jt)%potential%wx1(1:norb), sign_i, sign_j, ddotmat)

      hcc = czero
      if (trim(this%ccor_vmt_mode) == 'pair_surface') then
         lambda_pair = ccor_vmt_pair_surface(this, it, jt) - this%ccor_elin
         lambda0 = 0.5_rp*(lambda_pair(1) + lambda_pair(2))
         lambda1 = 0.5_rp*(lambda_pair(1) - lambda_pair(2))
         lambda_i = [lambda0 + sign_i*lambda1, lambda0 - sign_i*lambda1]
         lambda_j = [lambda0 + sign_j*lambda1, lambda0 - sign_j*lambda1]
         do ispin = 1, 2
            do jspin = 1, 2
               do ilm = 1, norb
                  irow = ilm + (ispin - 1)*spin_off
                  do jlm = 1, norb
                     jcol = jlm + (jspin - 1)*spin_off
                     hcc(irow, jcol) = 0.5_rp*(lambda_i(ispin) + lambda_j(jspin))*ddotmat(irow, jcol) + &
                        (ccd_i(ilm, 1)*lambda_i(ispin) + ccd_j(jlm, 1)*lambda_j(jspin))*dmat(irow, jcol)
                  end do
               end do
            end do
         end do
      else
         lambda0 = ccor_vmt_scalar(this) - this%ccor_elin
         do ispin = 1, 2
            do jspin = 1, 2
               do ilm = 1, norb
                  irow = ilm + (ispin - 1)*spin_off
                  do jlm = 1, norb
                     jcol = jlm + (jspin - 1)*spin_off
                     hcc(irow, jcol) = lambda0*(ddotmat(irow, jcol) + &
                        (ccd_i(ilm, 1) + ccd_j(jlm, 1))*dmat(irow, jcol))
                  end do
               end do
            end do
         end do
      end if

      ! Onsite CCOR is evaluated in the shared collinear rotating frame.
      ! d=0 gives alpha=0 and a common endpoint frame gives G_ii=I.
      if (m == 1) then
         do ispin = 1, 2
            do ilm = 1, norb
               irow = ilm + (ispin - 1)*spin_off
               wi = this%lattice%symbolic_atoms(it)%potential%wx0(ilm) + &
                    merge(sign_i, -sign_i, ispin == 1)* &
                    this%lattice%symbolic_atoms(it)%potential%wx1(ilm)
               if (trim(this%ccor_vmt_mode) == 'pair_surface') then
                  hcc(irow, irow) = hcc(irow, irow) + lambda_i(ispin)*ccd_i(ilm, 0)*wi*wi
               else
                  hcc(irow, irow) = hcc(irow, irow) + lambda0*ccd_i(ilm, 0)*wi*wi
               end if
            end do
         end do
      end if

      call hcpx(hcc(1:norb, 1:norb), 'cart2sph')
      call hcpx(hcc(1:norb, norb + 1:nb), 'cart2sph')
      call hcpx(hcc(norb + 1:nb, 1:norb), 'cart2sph')
      call hcpx(hcc(norb + 1:nb, norb + 1:nb), 'cart2sph')
   end subroutine build_ccor_pair_block_gbt

   !> @brief Build one noncollinear CCOR pair block.
   !> @details Computes spin-dependent two-centre combined-correction terms using
   !>          local moments, D components, and VMT coefficients for the atom pair.
   !> @param[in] this Hamiltonian object containing potentials and magnetic state.
   !> @param[in] ia Site index of the source atom.
   !> @param[in] ja Site index of the neighbor atom.
   !> @param[in] it Atom type of ia.
   !> @param[in] jt Atom type of ja.
   !> @param[in] ino Bravais/type index used by the target Hamiltonian layout.
   !> @param[in] m Neighbor-shell index.
   !> @param[out] hcc CCOR block in spin-orbital basis.
   module subroutine build_ccor_pair_block_noncollinear(this, ia, ja, it, jt, ino, m, hcc)
      class(hamiltonian), intent(in) :: this
      integer, intent(in) :: ia, ja, it, jt, ino, m
      complex(rp), dimension(nb, nb), intent(out) :: hcc

      complex(rp), dimension(norb, norb) :: s_block, sdot_block
      complex(rp), dimension(norb, norb, 4) :: dcomp, ddotcomp, kcomp
      real(rp), dimension(norb, 0:2) :: ccd_i, ccd_j
      real(rp), dimension(3) :: mom_i, mom_j
      real(rp) :: lambda
      integer :: ilm, jlm, idir

      hcc(:, :) = czero
      if (.not. allocated(this%charge%lattice%sdot)) return

      do ilm = 1, norb
         do jlm = 1, norb
            s_block(ilm, jlm) = this%charge%lattice%sbar(jlm, ilm, m, ino)
            sdot_block(ilm, jlm) = normalize_ccor_sdot(this, this%charge%lattice%sdot(jlm, ilm, m, ino))
         end do
      end do

      call build_ccor_coefficients(this, ia, it, ccd_i)
      call build_ccor_coefficients(this, ja, jt, ccd_j)
         call build_ccor_d_components(this, ia, ja, it, jt, s_block, dcomp)
         call build_ccor_d_components(this, ia, ja, it, jt, sdot_block, ddotcomp)

         if (trim(this%ccor_vmt_mode) == 'pair_surface') then
            call build_ccor_pair_surface_block(this, ia, ja, it, jt, m, dcomp, ddotcomp, ccd_i, ccd_j, hcc)
            return
         end if

         kcomp(:, :, :) = czero
      do ilm = 1, norb
         do jlm = 1, norb
            do idir = 1, 4
               kcomp(ilm, jlm, idir) = ddotcomp(ilm, jlm, idir) + &
                  (ccd_i(ilm, 1) + ccd_j(jlm, 1))*dcomp(ilm, jlm, idir)
            end do
         end do
      end do

      if (m == 1) then
         call ccor_select_endpoint_moments(this, ia, ja, it, jt, mom_i, mom_j)
         do ilm = 1, norb
            kcomp(ilm, ilm, 4) = kcomp(ilm, ilm, 4) + ccd_i(ilm, 0)* &
               (this%charge%lattice%symbolic_atoms(it)%potential%wx0(ilm)**2 + &
                this%charge%lattice%symbolic_atoms(it)%potential%wx1(ilm)**2)
            do idir = 1, 3
               kcomp(ilm, ilm, idir) = kcomp(ilm, ilm, idir) + ccd_i(ilm, 0)* &
                  2.0_rp*this%charge%lattice%symbolic_atoms(it)%potential%wx0(ilm)* &
                  this%charge%lattice%symbolic_atoms(it)%potential%wx1(ilm)*cmplx(mom_i(idir), 0.0_rp, rp)
            end do
         end do
      end if

      do idir = 1, 4
         call hcpx(kcomp(:, :, idir), 'cart2sph')
      end do

         lambda = ccor_vmt_scalar(this) - this%ccor_elin
      hcc(1:norb, 1:norb) = lambda*(kcomp(:, :, 4) + kcomp(:, :, 3))
      hcc(spin_off + 1:spin_off + norb, spin_off + 1:spin_off + norb) = lambda*(kcomp(:, :, 4) - kcomp(:, :, 3))
      hcc(1:norb, spin_off + 1:spin_off + norb) = lambda*(kcomp(:, :, 1) - i_unit*kcomp(:, :, 2))
         hcc(spin_off + 1:spin_off + norb, 1:norb) = lambda*(kcomp(:, :, 1) + i_unit*kcomp(:, :, 2))
      end subroutine build_ccor_pair_block_noncollinear

      !> @brief Build the surface-mode noncollinear CCOR pair block.
      !> @details Combines endpoint and pair VMT terms with D and D-dot components
      !>          for surfaces or local regions where global bulk assumptions do
      !>          not apply.
      !> @param[in] this Hamiltonian object containing CCOR mode and moments.
      !> @param[in] ia Site index of the source atom.
      !> @param[in] ja Site index of the neighbor atom.
      !> @param[in] it Atom type of ia.
      !> @param[in] jt Atom type of ja.
      !> @param[in] m Neighbor-shell index.
      !> @param[in] dcomp Spin decomposed D matrix components.
      !> @param[in] ddotcomp Spin decomposed D-dot matrix components.
      !> @param[in] ccd_i Combined-correction coefficients for type it.
      !> @param[in] ccd_j Combined-correction coefficients for type jt.
      !> @param[out] hcc CCOR block in spin-orbital basis.
      module subroutine build_ccor_pair_surface_block(this, ia, ja, it, jt, m, dcomp, ddotcomp, ccd_i, ccd_j, hcc)
         class(hamiltonian), intent(in) :: this
         integer, intent(in) :: ia, ja, it, jt, m
         complex(rp), dimension(norb, norb, 4), intent(in) :: dcomp, ddotcomp
         real(rp), dimension(norb, 0:2), intent(in) :: ccd_i, ccd_j
         complex(rp), dimension(nb, nb), intent(out) :: hcc

         complex(rp), dimension(norb, norb, 4) :: hcomp
         complex(rp), dimension(4) :: lambda_i, lambda_j, dloc, ddotloc, term_l, term_r, onsite
         real(rp), dimension(2) :: lambda_pair
         real(rp), dimension(3) :: mom_i, mom_j
         integer :: ilm, jlm, idir

         hcc(:, :) = czero
         hcomp(:, :, :) = czero
         lambda_pair(:) = (ccor_vmt_pair_surface(this, it, jt) - this%ccor_elin)

         call ccor_select_endpoint_moments(this, ia, ja, it, jt, mom_i, mom_j)
         call ccor_lambda_components(lambda_pair, mom_i, lambda_i)
         call ccor_lambda_components(lambda_pair, mom_j, lambda_j)

         do ilm = 1, norb
            do jlm = 1, norb
               dloc(:) = dcomp(ilm, jlm, :)
               ddotloc(:) = ddotcomp(ilm, jlm, :)

               call ccor_spin_product(lambda_i, ddotloc, term_l)
               call ccor_spin_product(ddotloc, lambda_j, term_r)
               hcomp(ilm, jlm, :) = 0.5_rp*(term_l(:) + term_r(:))

               call ccor_spin_product(lambda_i, dloc, term_l)
               call ccor_spin_product(dloc, lambda_j, term_r)
               hcomp(ilm, jlm, :) = hcomp(ilm, jlm, :) + ccd_i(ilm, 1)*term_l(:) + ccd_j(jlm, 1)*term_r(:)
            end do
         end do

         if (m == 1) then
            do ilm = 1, norb
               onsite(:) = czero
               onsite(4) = ccd_i(ilm, 0)*(this%charge%lattice%symbolic_atoms(it)%potential%wx0(ilm)**2 + &
                          this%charge%lattice%symbolic_atoms(it)%potential%wx1(ilm)**2)
               do idir = 1, 3
                  onsite(idir) = ccd_i(ilm, 0)*2.0_rp*this%charge%lattice%symbolic_atoms(it)%potential%wx0(ilm)* &
                                 this%charge%lattice%symbolic_atoms(it)%potential%wx1(ilm)*cmplx(mom_i(idir), 0.0_rp, rp)
               end do
               call ccor_spin_product(lambda_i, onsite, term_l)
               hcomp(ilm, ilm, :) = hcomp(ilm, ilm, :) + term_l(:)
            end do
         end if

         do idir = 1, 4
            call hcpx(hcomp(:, :, idir), 'cart2sph')
         end do

         hcc(1:norb, 1:norb) = hcomp(:, :, 4) + hcomp(:, :, 3)
         hcc(spin_off + 1:spin_off + norb, spin_off + 1:spin_off + norb) = hcomp(:, :, 4) - hcomp(:, :, 3)
         hcc(1:norb, spin_off + 1:spin_off + norb) = hcomp(:, :, 1) - i_unit*hcomp(:, :, 2)
         hcc(spin_off + 1:spin_off + norb, 1:norb) = hcomp(:, :, 1) + i_unit*hcomp(:, :, 2)
      end subroutine build_ccor_pair_surface_block

   !> @brief Decompose a structure-constant block into CCOR D components.
   !> @details Projects the orbital structure block onto scalar and spin-vector
   !>          components used by noncollinear CCOR pair-block construction.
   !> @param[in] this Hamiltonian object containing moments and spin-spiral state.
   !> @param[in] ia Site index of the source atom.
   !> @param[in] ja Site index of the neighbor atom.
   !> @param[in] it Atom type of ia.
   !> @param[in] jt Atom type of ja.
   !> @param[in] s_block Structure-constant block in orbital basis.
   !> @param[out] dcomp Decomposed D components.
   module subroutine build_ccor_d_components(this, ia, ja, it, jt, s_block, dcomp)
      class(hamiltonian), intent(in) :: this
      integer, intent(in) :: ia, ja, it, jt
      complex(rp), dimension(norb, norb), intent(in) :: s_block
      complex(rp), dimension(norb, norb, 4), intent(out) :: dcomp

      real(rp), dimension(3) :: mom_i, mom_j
      complex(rp), dimension(3) :: cross
      complex(rp) :: dotp
      integer :: ilm, jlm, idir

      call ccor_select_endpoint_moments(this, ia, ja, it, jt, mom_i, mom_j)
      dotp = cmplx(dot_product(mom_i, mom_j), 0.0_rp, rp)
      cross = cmplx(cross_product(mom_i, mom_j), 0.0_rp, rp)

      dcomp(:, :, :) = czero
      do ilm = 1, norb
         do jlm = 1, norb
            dcomp(ilm, jlm, 4) = &
               this%charge%lattice%symbolic_atoms(it)%potential%wx0(ilm)*s_block(ilm, jlm)* &
               this%charge%lattice%symbolic_atoms(jt)%potential%wx0(jlm) + &
               this%charge%lattice%symbolic_atoms(it)%potential%wx1(ilm)*s_block(ilm, jlm)* &
               this%charge%lattice%symbolic_atoms(jt)%potential%wx1(jlm)*dotp
            do idir = 1, 3
               dcomp(ilm, jlm, idir) = &
                  this%charge%lattice%symbolic_atoms(it)%potential%wx1(ilm)*s_block(ilm, jlm)* &
                  this%charge%lattice%symbolic_atoms(jt)%potential%wx0(jlm)*cmplx(mom_i(idir), 0.0_rp, rp) + &
                  this%charge%lattice%symbolic_atoms(it)%potential%wx0(ilm)*s_block(ilm, jlm)* &
                  this%charge%lattice%symbolic_atoms(jt)%potential%wx1(jlm)*cmplx(mom_j(idir), 0.0_rp, rp) + &
                  i_unit*this%charge%lattice%symbolic_atoms(it)%potential%wx1(ilm)*s_block(ilm, jlm)* &
                  this%charge%lattice%symbolic_atoms(jt)%potential%wx1(jlm)*cross(idir)
            end do
         end do
      end do
   end subroutine build_ccor_d_components

   !> @brief Compute orbital-resolved CCOR coefficients for one atom type.
   !> @details Evaluates the energy-linearized combined-correction coefficients
   !>          from potential parameters for the specified site/type.
   !> @param[in] this Hamiltonian object containing symbolic-atom potentials.
   !> @param[in] ia Site index used for local moment/context.
   !> @param[in] itype Atom type whose coefficients are requested.
   !> @param[out] ccd Coefficients indexed by orbital and angular channel.
   module subroutine build_ccor_coefficients(this, ia, itype, ccd)
      class(hamiltonian), intent(in) :: this
      integer, intent(in) :: ia, itype
      real(rp), dimension(norb, 0:2), intent(out) :: ccd

      integer :: ilm, l, alpha_idx
      real(rp) :: sw, dl, a, adot, l2p3, l2p1, l2m1, t1, t2, t3, avw

      ccd(:, :) = 0.0_rp
      avw = this%charge%lattice%wav*ang2au
      if (avw <= tiny(1.0_rp)) avw = this%charge%lattice%wav
      if (avw <= tiny(1.0_rp)) call g_logger%fatal('ccor_2c requires positive lattice%wav.', __FILE__, __LINE__)
      sw = this%charge%lattice%symbolic_atoms(itype)%potential%ws_r/avw
      if (sw <= tiny(1.0_rp)) call g_logger%fatal('ccor_2c requires positive ws_r/avw.', __FILE__, __LINE__)

      do ilm = 1, norb
         l = orbital_l_from_index(ilm)
         l2p3 = real(2*l + 3, rp)
         l2p1 = real(2*l + 1, rp)
         l2m1 = real(2*l - 1, rp)
         dl = sw**(-(2*l - 1))/(0.5_rp*l2m1)
         a = 0.0_rp
         if (allocated(this%charge%lattice%symbolic_atoms(itype)%potential%screening_alpha)) then
            if (lbound(this%charge%lattice%symbolic_atoms(itype)%potential%screening_alpha, 1) <= l .and. &
                ubound(this%charge%lattice%symbolic_atoms(itype)%potential%screening_alpha, 1) >= l) then
               a = this%charge%lattice%symbolic_atoms(itype)%potential%screening_alpha(l)
            end if
         else if (allocated(this%charge%lattice%alpha)) then
            alpha_idx = min(max(1, this%charge%lattice%num(max(1, min(ia, size(this%charge%lattice%num))))), size(this%charge%lattice%alpha, 2))
            a = this%charge%lattice%alpha(ilm, alpha_idx, max(1, min(ia, size(this%charge%lattice%alpha, 3))))
         end if
         adot = 0.0_rp
         if (allocated(this%charge%lattice%alpha_dot)) then
            alpha_idx = min(max(1, this%charge%lattice%num(max(1, min(ia, size(this%charge%lattice%num))))), size(this%charge%lattice%alpha_dot, 2))
            adot = this%charge%lattice%alpha_dot(ilm, alpha_idx, max(1, min(ia, size(this%charge%lattice%alpha_dot, 3))))
         end if
         t1 = -sw**(2*l + 3)/(2.0_rp*l2p1*l2p1*l2p3)
         t2 = a*sw*sw/l2p1
         t3 = a*a*2.0_rp*sw**(-(2*l - 1))/l2m1
         ccd(ilm, 0) = avw*avw*dl
         ccd(ilm, 1) = avw*avw*(sw*sw/(2.0_rp*l2p1) + a*dl)
         ccd(ilm, 2) = avw*avw*(adot + t1 + t2 + t3)
      end do
   end subroutine build_ccor_coefficients

   !> @brief Normalize a CCOR spin-dot product.
   !> @details Applies the selected CCOR normalization convention to a raw
   !>          spin-product scalar before it enters pair-block construction.
   !> @param[in] this Hamiltonian object with CCOR normalization settings.
   !> @param[in] sdot_raw Raw complex spin-dot value.
   !> @return Normalized spin-dot value.
   module function normalize_ccor_sdot(this, sdot_raw) result(sdot_cc)
      class(hamiltonian), intent(in) :: this
      complex(rp), intent(in) :: sdot_raw
      complex(rp) :: sdot_cc
      real(rp) :: avw

      avw = this%charge%lattice%wav*ang2au
      if (avw <= tiny(1.0_rp)) avw = this%charge%lattice%wav
      sdot_cc = -avw*avw*sdot_raw
   end function normalize_ccor_sdot

      !> @brief Return the scalar VMT value for CCOR.
      !> @details Selects the configured scalar VMT strategy and reduces any
      !>          spin-resolved surface estimate to the scalar value needed by
      !>          scalar CCOR blocks.
      !> @param[in] this Hamiltonian object with ccor_vmt_mode.
      !> @return Scalar VMT value in internal units.
      module function ccor_vmt_scalar(this) result(vmt)
         class(hamiltonian), intent(in) :: this
         real(rp) :: vmt
         real(rp), dimension(2) :: vmt_spin

         select case (trim(this%ccor_vmt_mode))
         case ('surface_scalar')
            vmt_spin = ccor_vmt_global_surface(this)
            vmt = 0.5_rp*(vmt_spin(1) + vmt_spin(2))
         case ('vmad_scalar')
            vmt = ccor_vmt_scalar_from_vmad(this)
         case default
            call g_logger%fatal('ccor_2c scalar VMT requested with incompatible ccor_vmt_mode.', __FILE__, __LINE__)
         end select
      end function ccor_vmt_scalar

      !> @brief Compute a global spin-resolved surface VMT estimate.
      !> @details Averages endpoint surface VMT values over symbolic atom types and
      !>          multiplicities to provide one spin-resolved correction scale.
      !> @param[in] this Hamiltonian object containing symbolic atoms.
      !> @return Spin-resolved VMT pair.
      module function ccor_vmt_global_surface(this) result(vmt_spin)
         class(hamiltonian), intent(in) :: this
         real(rp), dimension(2) :: vmt_spin
         real(rp), dimension(2) :: endpoint
         real(rp) :: weight, weight_sum
         integer :: itype, nsite, multiplicity

         vmt_spin(:) = 0.0_rp
         weight_sum = 0.0_rp
         nsite = min(this%charge%lattice%nbas, size(this%charge%lattice%iz))
         do itype = 1, this%charge%lattice%ntype
            multiplicity = count(this%charge%lattice%iz(1:nsite) == itype)
            if (multiplicity <= 0) multiplicity = 1
            weight = real(multiplicity, rp)*this%charge%lattice%symbolic_atoms(itype)%potential%ws_r**2
            endpoint = ccor_vmt_endpoint_surface(this, itype)
            vmt_spin(:) = vmt_spin(:) + weight*endpoint(:)
            weight_sum = weight_sum + weight
         end do
         if (weight_sum <= tiny(1.0_rp)) call g_logger%fatal('ccor_2c surface_scalar VMTZ has zero WS-radius weight.', __FILE__, __LINE__)
         vmt_spin(:) = vmt_spin(:)/weight_sum
      end function ccor_vmt_global_surface

      !> @brief Compute a pair-specific spin-resolved surface VMT estimate.
      !> @details Blends endpoint VMT values for two atom types with pair weights
      !>          for pair_surface CCOR mode.
      !> @param[in] this Hamiltonian object containing symbolic atoms.
      !> @param[in] itype First atom type.
      !> @param[in] jtype Second atom type.
      !> @return Spin-resolved VMT pair.
      module function ccor_vmt_pair_surface(this, itype, jtype) result(vmt_spin)
         class(hamiltonian), intent(in) :: this
         integer, intent(in) :: itype, jtype
         real(rp), dimension(2) :: vmt_spin, vi, vj
         real(rp) :: wi, wj

         vi = ccor_vmt_endpoint_surface(this, itype)
         vj = ccor_vmt_endpoint_surface(this, jtype)
         wi = this%charge%lattice%symbolic_atoms(itype)%potential%ws_r**2
         wj = this%charge%lattice%symbolic_atoms(jtype)%potential%ws_r**2
         if (wi + wj <= tiny(1.0_rp)) call g_logger%fatal('ccor_2c pair_surface VMTZ has zero endpoint WS-radius weight.', __FILE__, __LINE__)
         vmt_spin(:) = (wi*vi(:) + wj*vj(:))/(wi + wj)
      end function ccor_vmt_pair_surface

      !> @brief Compute an endpoint spin-resolved surface VMT estimate.
      !> @details Extracts the local surface VMT proxy for one symbolic atom type
      !>          from the available potential data.
      !> @param[in] this Hamiltonian object containing symbolic atoms.
      !> @param[in] itype Atom type whose endpoint value is requested.
      !> @return Spin-resolved VMT pair.
      module function ccor_vmt_endpoint_surface(this, itype) result(vmt_spin)
         class(hamiltonian), intent(in) :: this
         integer, intent(in) :: itype
         real(rp), dimension(2) :: vmt_spin
         real(rp) :: avg, diff

         avg = this%charge%lattice%symbolic_atoms(itype)%potential%vrmax(1)
         diff = this%charge%lattice%symbolic_atoms(itype)%potential%vrmax(2)
         if (abs(avg) <= tiny(1.0_rp) .and. abs(diff) <= tiny(1.0_rp) .and. this%ccor_2c) then
            if (this%ccor_strict) then
               call g_logger%fatal('ccor_2c surface VMTZ requested, but potential%vrmax is zero/unset.', __FILE__, __LINE__)
            else if (rank == 0) then
               call g_logger%warning('ccor_2c surface VMTZ has zero/unset potential%vrmax; falling back to potential%vmad for this endpoint. Regenerate *_out.nml to persist vrmax for production CCOR.', __FILE__, __LINE__)
            end if
            avg = this%charge%lattice%symbolic_atoms(itype)%potential%vmad
            diff = 0.0_rp
         end if
         vmt_spin(1) = avg + 0.5_rp*diff
         vmt_spin(2) = avg - 0.5_rp*diff
      end function ccor_vmt_endpoint_surface

      !> @brief Compute a scalar CCOR VMT estimate from Madelung potentials.
      !> @details Averages available vmad-like potential information over atom
      !>          types for the vmad_scalar CCOR mode.
      !> @param[in] this Hamiltonian object containing symbolic atoms.
      !> @return Scalar VMT value in internal units.
      module function ccor_vmt_scalar_from_vmad(this) result(vmt)
         class(hamiltonian), intent(in) :: this
         real(rp) :: vmt
         real(rp) :: weight, weight_sum
         integer :: itype

         vmt = 0.0_rp
         weight_sum = 0.0_rp
      do itype = 1, this%charge%lattice%ntype
         weight = this%charge%lattice%symbolic_atoms(itype)%potential%ws_r**2
         vmt = vmt + weight*this%charge%lattice%symbolic_atoms(itype)%potential%vmad
         weight_sum = weight_sum + weight
      end do
         if (weight_sum <= tiny(1.0_rp)) call g_logger%fatal('ccor_2c vmad_scalar VMT has zero WS-radius weight.', __FILE__, __LINE__)
         vmt = vmt/weight_sum
      end function ccor_vmt_scalar_from_vmad

      !> @brief Expand spin-resolved lambda values into Pauli components.
      !> @details Converts the up/down lambda pair and local moment direction into
      !>          scalar/vector components used by noncollinear CCOR algebra.
      !> @param[in] lambda_pair Spin-resolved lambda values.
      !> @param[in] mom Local magnetic moment direction.
      !> @param[out] lambda_comp Scalar and vector lambda components.
      module subroutine ccor_lambda_components(lambda_pair, mom, lambda_comp)
         real(rp), dimension(2), intent(in) :: lambda_pair
         real(rp), dimension(3), intent(in) :: mom
         complex(rp), dimension(4), intent(out) :: lambda_comp
         real(rp) :: lambda0, lambda1

         lambda0 = 0.5_rp*(lambda_pair(1) + lambda_pair(2))
         lambda1 = 0.5_rp*(lambda_pair(1) - lambda_pair(2))
         lambda_comp(1) = cmplx(lambda1*mom(1), 0.0_rp, rp)
         lambda_comp(2) = cmplx(lambda1*mom(2), 0.0_rp, rp)
         lambda_comp(3) = cmplx(lambda1*mom(3), 0.0_rp, rp)
         lambda_comp(4) = cmplx(lambda0, 0.0_rp, rp)
      end subroutine ccor_lambda_components

      !> @brief Multiply two scalar/vector spin-component tuples.
      !> @details Applies Pauli-matrix product algebra to compact four-component
      !>          spin decompositions used by noncollinear CCOR construction.
      !> @param[in] a Left scalar/vector spin tuple.
      !> @param[in] b Right scalar/vector spin tuple.
      !> @param[out] c Product tuple.
      module subroutine ccor_spin_product(a, b, c)
         complex(rp), dimension(4), intent(in) :: a, b
         complex(rp), dimension(4), intent(out) :: c

         c(4) = a(4)*b(4) + a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
         c(1) = a(4)*b(1) + b(4)*a(1) + i_unit*(a(2)*b(3) - a(3)*b(2))
         c(2) = a(4)*b(2) + b(4)*a(2) + i_unit*(a(3)*b(1) - a(1)*b(3))
         c(3) = a(4)*b(3) + b(4)*a(3) + i_unit*(a(1)*b(2) - a(2)*b(1))
      end subroutine ccor_spin_product

   !> @brief Return angular momentum l from a packed orbital index.
   !> @details Maps the LMTO orbital index ilm to its shell quantum number for
   !>          coefficient lookup and export helpers.
   !> @param[in] ilm Packed orbital index.
   !> @return Angular momentum shell l.
   module integer pure function orbital_l_from_index(ilm) result(l)
      integer, intent(in) :: ilm
      integer :: lp1

      l = 0
      do lp1 = 1, lmax_basis + 1
         if (ilm <= lp1*lp1) then
            l = lp1 - 1
            return
         end if
      end do
      l = lmax_basis
   end function orbital_l_from_index

   !> @brief Emit diagnostics for generated CCOR blocks.
   !> @details Summarizes CCOR coefficient and block magnitudes for debugging the
   !>          two-centre correction without changing Hamiltonian data.
   !> @param[in] this Hamiltonian object containing CCOR settings.
   !> @param[in] hcc CCOR block array to summarize.
   !> @param[in] label Human-readable label for the diagnostic message.
   module subroutine log_ccor_debug(this, hcc, label)
      class(hamiltonian), intent(in) :: this
      complex(rp), dimension(:, :, :, :), intent(in) :: hcc
      character(len=*), intent(in) :: label
      real(rp), dimension(norb, 0:2) :: ccd
         real(rp) :: vmt, lambda, rms_ccd2, max_ccd2, avw
         real(rp), dimension(2) :: vmt_spin
      integer :: itype, count_ccd
      character(len=256) :: msg

            if (this%ccor_vmt_mode == 'surface_scalar' .or. this%ccor_vmt_mode == 'pair_surface') then
               vmt_spin = ccor_vmt_global_surface(this)
               vmt = 0.5_rp*(vmt_spin(1) + vmt_spin(2))
            else
               vmt = ccor_vmt_scalar(this)
               vmt_spin(:) = vmt
            end if
         lambda = vmt - this%ccor_elin
      avw = this%charge%lattice%wav*ang2au
      if (avw <= tiny(1.0_rp)) avw = this%charge%lattice%wav
      rms_ccd2 = 0.0_rp
      max_ccd2 = 0.0_rp
      count_ccd = 0
      do itype = 1, this%charge%lattice%ntype
         call build_ccor_coefficients(this, this%charge%lattice%atlist(itype), itype, ccd)
         rms_ccd2 = rms_ccd2 + sum(ccd(:, 2)**2)
         max_ccd2 = max(max_ccd2, maxval(abs(ccd(:, 2))))
         count_ccd = count_ccd + size(ccd(:, 2))
      end do
      if (count_ccd > 0) rms_ccd2 = sqrt(rms_ccd2/real(count_ccd, rp))
         write(msg, '(a,a,a,f14.8,a,a,a,f14.8,a,f14.8,a,f14.8)') 'CCOR2C ', trim(label), &
            ' E_lin=', this%ccor_elin, ' VMT mode=', trim(this%ccor_vmt_mode), ' VMT=', vmt, ' lambda=', lambda, ' avw=', avw
         call g_logger%info(trim(msg), __FILE__, __LINE__)
         write(msg, '(a,a,a,f14.8,a,f14.8)') 'CCOR2C ', trim(label), &
            ' VMT_up=', vmt_spin(1), ' VMT_down=', vmt_spin(2)
         call g_logger%info(trim(msg), __FILE__, __LINE__)
      write(msg, '(a,a,a,es12.4,a,es12.4,a,es12.4,a,es12.4)') 'CCOR2C ', trim(label), &
         ' maxabs(S)=', maxval(abs(this%charge%lattice%sbar)), &
         ' maxabs(Sdot_raw)=', maxval(abs(this%charge%lattice%sdot)), &
         ' maxabs(Hcc)=', maxval(abs(hcc)), ' rms_ccd2=', rms_ccd2
      call g_logger%info(trim(msg), __FILE__, __LINE__)
      write(msg, '(a,a,a,es12.4)') 'CCOR2C ', trim(label), ' maxabs_ccd2=', max_ccd2
      call g_logger%info(trim(msg), __FILE__, __LINE__)
      if (max_ccd2 > 10.0_rp) then
         if (this%ccor_strict) then
            call g_logger%fatal('CCOR2C ccd2 diagnostic is unexpectedly large.', __FILE__, __LINE__)
         else
            call g_logger%warning('CCOR2C ccd2 diagnostic is unexpectedly large; two-centre approximation may be poor.', __FILE__, __LINE__)
         end if
      end if
   end subroutine log_ccor_debug

end submodule hamiltonian_ccor
