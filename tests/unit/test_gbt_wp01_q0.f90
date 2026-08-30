program test_gbt_wp01_q0
   ! WP01 is deliberately a fixed-potential operator test.  The ordinary side
   ! uses the ordinary LMTO value path; the GBT side uses only the primitive
   ! structure-constant link and contraction helpers.  No GBT helper is used
   ! to construct the ordinary reference.
   use precision_mod, only: rp
   use lmto_magnetic_tangent_mod, only: lmto_bond_value, lmto_hhmag_to_spinor
   use gbt_structure_mod, only: gbt_endpoint_link, gbt_lift_orbital_block, gbt_contract_collinear
   implicit none

   integer, parameter :: norb_fixture = 2
   integer, parameter :: nb_fixture = 2*norb_fixture
   integer, parameter :: nsite_fixture = 2
   integer, parameter :: nbond_fixture = 4
   integer, parameter :: nmat_fixture = nsite_fixture*nb_fixture
   integer, parameter :: nk_fixture = 3
   real(rp), parameter :: tol = 1.0e-12_rp
   real(rp), parameter :: negative_tol = 1.0e-3_rp
   real(rp), parameter :: alat = 2.8612_rp
   real(rp), parameter :: pi = acos(-1.0_rp)
   real(rp), parameter :: q_zero(3) = [0.0_rp, 0.0_rp, 0.0_rp]

   integer, parameter :: source_site(nbond_fixture) = [1, 2, 1, 2]
   integer, parameter :: target_site(nbond_fixture) = [1, 2, 2, 1]
   real(rp), parameter :: translation(3, nbond_fixture) = reshape([ &
      0.0_rp, 0.0_rp, 0.0_rp, &
      0.0_rp, 0.0_rp, 0.0_rp, &
      1.0_rp, 0.0_rp, 0.0_rp, &
     -1.0_rp, 0.0_rp, 0.0_rp], [3, nbond_fixture])
   real(rp), parameter :: k_points(3, nk_fixture) = reshape([ &
       0.137_rp, -0.219_rp, 0.083_rp, &
      -0.271_rp,  0.314_rp, 0.119_rp, &
       0.223_rp,  0.071_rp, -0.307_rp], [3, nk_fixture])

   complex(rp) :: raw(norb_fixture, norb_fixture, nbond_fixture)
   complex(rp) :: w0(norb_fixture, nsite_fixture), w1(norb_fixture, nsite_fixture)
   complex(rp) :: c0(norb_fixture, nsite_fixture), c1(norb_fixture, nsite_fixture)
   real(rp) :: theta(nsite_fixture), phi(nsite_fixture)
   complex(rp) :: ordinary(nmat_fixture, nmat_fixture, nbond_fixture)
   complex(rp) :: gbt(nmat_fixture, nmat_fixture, nbond_fixture)
   complex(rp) :: ordinary_collinear(nmat_fixture, nmat_fixture, nbond_fixture)
   complex(rp) :: gbt_global(nmat_fixture, nmat_fixture, nbond_fixture)
   complex(rp) :: ordinary_multi(nmat_fixture, nmat_fixture, nbond_fixture)
   complex(rp) :: gbt_multi(nmat_fixture, nmat_fixture, nbond_fixture)
   complex(rp) :: bad_gbt_multi(nmat_fixture, nmat_fixture, nbond_fixture)
   complex(rp) :: frame_global(nmat_fixture, nmat_fixture)
   complex(rp) :: frame_multi(nmat_fixture, nmat_fixture)
   logical :: failed
   real(rp) :: largest_abs, largest_rel, largest_eigen

   failed = .false.
   largest_abs = 0.0_rp
   largest_rel = 0.0_rp
   largest_eigen = 0.0_rp
   call initialize_fixture()
   call test_q0_collinear()
   call test_q0_global_rotation()
   call test_q0_multisublattice()

   write (*, '(a,3(es12.4,1x))') 'WP01 overall maxima (abs, rel, eigen): ', largest_abs, largest_rel, largest_eigen
   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine initialize_fixture()
      ! The orbital matrices are a deterministic finite structure-constant
      ! table.  The last bond is the exact reverse of the third bond.
      raw = cmplx(0.0_rp, 0.0_rp, rp)
      raw(1,1,1) = cmplx(0.41_rp, 0.0_rp, rp)
      raw(1,2,1) = cmplx(0.07_rp, 0.03_rp, rp)
      raw(2,1,1) = conjg(raw(1,2,1))
      raw(2,2,1) = cmplx(0.53_rp, 0.0_rp, rp)
      raw(:,:,2) = raw(:,:,1)
      raw(1,1,3) = cmplx(0.70_rp, 0.20_rp, rp)
      raw(1,2,3) = cmplx(-0.30_rp, 0.40_rp, rp)
      raw(2,1,3) = cmplx(0.10_rp, -0.50_rp, rp)
      raw(2,2,3) = cmplx(0.80_rp, -0.10_rp, rp)
      raw(:,:,4) = transpose(conjg(raw(:,:,3)))

      w0(:,1) = [cmplx(1.08_rp, 0.0_rp, rp), cmplx(0.74_rp, 0.0_rp, rp)]
      w1(:,1) = [cmplx(0.23_rp, 0.0_rp, rp), cmplx(-0.17_rp, 0.0_rp, rp)]
      w0(:,2) = [cmplx(0.86_rp, 0.0_rp, rp), cmplx(1.16_rp, 0.0_rp, rp)]
      w1(:,2) = [cmplx(-0.31_rp, 0.0_rp, rp), cmplx(0.28_rp, 0.0_rp, rp)]
      c0(:,1) = [cmplx(-0.19_rp, 0.0_rp, rp), cmplx(0.27_rp, 0.0_rp, rp)]
      c1(:,1) = [cmplx(0.39_rp, 0.0_rp, rp), cmplx(-0.13_rp, 0.0_rp, rp)]
      c0(:,2) = [cmplx(-0.11_rp, 0.0_rp, rp), cmplx(0.31_rp, 0.0_rp, rp)]
      c1(:,2) = [cmplx(0.22_rp, 0.0_rp, rp), cmplx(-0.26_rp, 0.0_rp, rp)]
   end subroutine initialize_fixture

   subroutine test_q0_collinear()
      real(rp) :: abs_blocks, rel_blocks, abs_onsite, rel_onsite
      real(rp) :: abs_k, rel_k, eig_k

      theta = 0.0_rp
      phi = 0.0_rp
      call build_operators(theta, phi, ordinary_collinear, gbt)
      call compare_block_sets(gbt, ordinary_collinear, abs_blocks, rel_blocks)
      call compare_onsite_blocks(gbt, ordinary_collinear, abs_onsite, rel_onsite)
      call compare_kspace_sets(gbt, ordinary_collinear, frame_identity(), abs_k, rel_k, eig_k)
      call report('q=0 collinear directed blocks', abs_blocks, rel_blocks, 0.0_rp)
      call report('q=0 collinear onsite blocks', abs_onsite, rel_onsite, 0.0_rp)
      call report('q=0 collinear reciprocal H(k)', abs_k, rel_k, eig_k)
   end subroutine test_q0_collinear

   subroutine test_q0_global_rotation()
      real(rp) :: abs_transform, rel_transform, abs_blocks, rel_blocks
      real(rp) :: abs_k, rel_k, eig_k
      real(rp) :: global_theta, global_phi
      complex(rp) :: rotated(nmat_fixture, nmat_fixture, nbond_fixture)

      global_theta = 37.0_rp*pi/180.0_rp
      global_phi = 0.41_rp
      theta = [global_theta, global_theta]
      phi = [global_phi, global_phi]
      call build_operators([0.0_rp, 0.0_rp], [0.0_rp, 0.0_rp], ordinary_collinear, gbt_global)
      call build_operators(theta, phi, rotated, gbt_global)
      ordinary = rotated
      call make_packed_frame(theta, phi, frame_global)
      call compare_block_sets_transformed(gbt_global, rotated, frame_global, abs_blocks, rel_blocks)
      call compare_block_sets_rotated(rotated, ordinary_collinear, frame_global, abs_transform, rel_transform)
      call compare_kspace_sets(gbt_global, rotated, frame_global, abs_k, rel_k, eig_k)
      call report('q=0 global SU(2) transformed blocks', abs_blocks, rel_blocks, 0.0_rp)
      call report('q=0 ordinary global rotation', abs_transform, rel_transform, 0.0_rp)
      call report('q=0 global SU(2) reciprocal H(k)', abs_k, rel_k, eig_k)
   end subroutine test_q0_global_rotation

   subroutine test_q0_multisublattice()
      real(rp) :: abs_blocks, rel_blocks, abs_onsite, rel_onsite
      real(rp) :: abs_k, rel_k, eig_k, negative_error
      complex(rp) :: expected_block(nb_fixture, nb_fixture)
      integer :: i0, j0

      theta = [0.63_rp, 1.11_rp]
      phi = [0.27_rp, -0.58_rp]
      call build_operators(theta, phi, ordinary_multi, gbt_multi)
      call make_packed_frame(theta, phi, frame_multi)
      call compare_block_sets_transformed(gbt_multi, ordinary_multi, frame_multi, abs_blocks, rel_blocks)
      call compare_onsite_blocks_transformed(gbt_multi, ordinary_multi, frame_multi, abs_onsite, rel_onsite)
      call compare_kspace_sets(gbt_multi, ordinary_multi, frame_multi, abs_k, rel_k, eig_k)
      call report('q=0 multi-sublattice transformed blocks', abs_blocks, rel_blocks, 0.0_rp)
      call report('q=0 multi-sublattice onsite blocks', abs_onsite, rel_onsite, 0.0_rp)
      call report('q=0 multi-sublattice reciprocal H(k)', abs_k, rel_k, eig_k)

      call build_operators(theta, phi, ordinary, bad_gbt_multi, 3)
      i0 = site_offset(source_site(3)); j0 = site_offset(target_site(3))
      expected_block = matmul(transpose(conjg(site_frame(theta(1), phi(1)))), &
         matmul(ordinary_multi(i0:i0+nb_fixture-1,j0:j0+nb_fixture-1,3), site_frame(theta(2), phi(2))))
      negative_error = maxval(abs(bad_gbt_multi(i0:i0+nb_fixture-1,j0:j0+nb_fixture-1,3)-expected_block))
      write (*, '(a,es12.4)') 'q=0 multi-sublattice negative-control residual: ', negative_error
      if (negative_error <= negative_tol) then
         write (*, '(a)') 'FAIL q=0 negative control did not detect a perturbed endpoint frame'
         failed = .true.
      end if
   end subroutine test_q0_multisublattice

   subroutine build_operators(theta_in, phi_in, ordinary_out, gbt_out, perturb_bond)
      real(rp), intent(in) :: theta_in(:), phi_in(:)
      complex(rp), intent(out) :: ordinary_out(:,:,:), gbt_out(:,:,:)
      integer, intent(in), optional :: perturb_bond
      complex(rp) :: hhmag(norb_fixture, norb_fixture, 4)
      complex(rp) :: linked(nb_fixture, nb_fixture), link(2,2), pair(nb_fixture, nb_fixture)
      real(rp) :: mom_i(3), mom_j(3), theta_i, theta_j
      integer :: bond, i

      ordinary_out = cmplx(0.0_rp, 0.0_rp, rp)
      gbt_out = cmplx(0.0_rp, 0.0_rp, rp)
      do bond = 1, nbond_fixture
         call independent_moment(theta_in(source_site(bond)), phi_in(source_site(bond)), mom_i)
         call independent_moment(theta_in(target_site(bond)), phi_in(target_site(bond)), mom_j)
         call lmto_bond_value(raw(:,:,bond), w0(:,source_site(bond)), w1(:,source_site(bond)), &
            w0(:,target_site(bond)), w1(:,target_site(bond)), c0(:,source_site(bond)), &
            c1(:,source_site(bond)), mom_i, mom_j, bond <= 2, hhmag)
         call lmto_hhmag_to_spinor(hhmag, ordinary_out(site_offset(source_site(bond)):site_offset(source_site(bond))+nb_fixture-1, &
            site_offset(target_site(bond)):site_offset(target_site(bond))+nb_fixture-1, bond))

         theta_i = theta_in(source_site(bond))
         theta_j = theta_in(target_site(bond))
         if (present(perturb_bond)) then
            if (bond == perturb_bond) theta_i = theta_i + 0.17_rp
         end if
         call gbt_endpoint_link(q_zero, translation(:,bond)*alat, alat, theta_i, phi_in(source_site(bond)), &
            theta_j, phi_in(target_site(bond)), link)
         call gbt_lift_orbital_block(raw(:,:,bond), link, linked)
         call gbt_contract_collinear(linked, w0(:,source_site(bond)), w1(:,source_site(bond)), &
            w0(:,target_site(bond)), w1(:,target_site(bond)), 1.0_rp, 1.0_rp, pair)
         if (bond <= 2) then
            do i = 1, norb_fixture
               pair(i,i) = pair(i,i) + c0(i,source_site(bond)) + c1(i,source_site(bond))
               pair(i+norb_fixture,i+norb_fixture) = pair(i+norb_fixture,i+norb_fixture) + &
                  c0(i,source_site(bond)) - c1(i,source_site(bond))
            end do
         end if
         gbt_out(site_offset(source_site(bond)):site_offset(source_site(bond))+nb_fixture-1, &
            site_offset(target_site(bond)):site_offset(target_site(bond))+nb_fixture-1, bond) = pair
      end do
   end subroutine build_operators

   subroutine compare_block_sets(candidate, reference, abs_error, rel_error)
      complex(rp), intent(in) :: candidate(:,:,:), reference(:,:,:)
      real(rp), intent(out) :: abs_error, rel_error
      abs_error = maxval(abs(candidate-reference))
      rel_error = norm_error_3d(candidate, reference)
   end subroutine compare_block_sets

   subroutine compare_onsite_blocks(candidate, reference, abs_error, rel_error)
      complex(rp), intent(in) :: candidate(:,:,:), reference(:,:,:)
      real(rp), intent(out) :: abs_error, rel_error
      complex(rp) :: left(nb_fixture,nb_fixture), right(nb_fixture,nb_fixture)
      integer :: site
      abs_error = 0.0_rp
      rel_error = 0.0_rp
      do site = 1, nsite_fixture
         left = candidate(site_offset(site):site_offset(site)+nb_fixture-1, &
                          site_offset(site):site_offset(site)+nb_fixture-1, site)
         right = reference(site_offset(site):site_offset(site)+nb_fixture-1, &
                           site_offset(site):site_offset(site)+nb_fixture-1, site)
         abs_error = max(abs_error, maxval(abs(left-right)))
         rel_error = max(rel_error, norm_error(left, right))
      end do
   end subroutine compare_onsite_blocks

   subroutine compare_block_sets_transformed(candidate, reference, frame, abs_error, rel_error)
      complex(rp), intent(in) :: candidate(:,:,:), reference(:,:,:), frame(:,:)
      real(rp), intent(out) :: abs_error, rel_error
      complex(rp) :: expected(nb_fixture,nb_fixture), ref_block(nb_fixture,nb_fixture)
      complex(rp) :: candidate_block(nb_fixture,nb_fixture), fs(nb_fixture,nb_fixture), ft(nb_fixture,nb_fixture)
      integer :: bond, i0, j0
      abs_error = 0.0_rp
      rel_error = 0.0_rp
      do bond = 1, nbond_fixture
         fs = frame(site_offset(source_site(bond)):site_offset(source_site(bond))+nb_fixture-1, &
                    site_offset(source_site(bond)):site_offset(source_site(bond))+nb_fixture-1)
         ft = frame(site_offset(target_site(bond)):site_offset(target_site(bond))+nb_fixture-1, &
                    site_offset(target_site(bond)):site_offset(target_site(bond))+nb_fixture-1)
         i0 = site_offset(source_site(bond)); j0 = site_offset(target_site(bond))
         ref_block = reference(i0:i0+nb_fixture-1,j0:j0+nb_fixture-1,bond)
         candidate_block = candidate(i0:i0+nb_fixture-1,j0:j0+nb_fixture-1,bond)
         expected = matmul(transpose(conjg(fs)), matmul(ref_block, ft))
         abs_error = max(abs_error, maxval(abs(candidate_block-expected)))
         rel_error = max(rel_error, norm_error(candidate_block, expected))
      end do
   end subroutine compare_block_sets_transformed

   subroutine compare_block_sets_rotated(candidate, collinear, frame, abs_error, rel_error)
      complex(rp), intent(in) :: candidate(:,:,:), collinear(:,:,:), frame(:,:)
      real(rp), intent(out) :: abs_error, rel_error
      complex(rp) :: expected(nb_fixture,nb_fixture), ref_block(nb_fixture,nb_fixture)
      complex(rp) :: candidate_block(nb_fixture,nb_fixture), fs(nb_fixture,nb_fixture), ft(nb_fixture,nb_fixture)
      integer :: bond, i0, j0
      abs_error = 0.0_rp
      rel_error = 0.0_rp
      do bond = 1, nbond_fixture
         fs = frame(site_offset(source_site(bond)):site_offset(source_site(bond))+nb_fixture-1, &
                    site_offset(source_site(bond)):site_offset(source_site(bond))+nb_fixture-1)
         ft = frame(site_offset(target_site(bond)):site_offset(target_site(bond))+nb_fixture-1, &
                    site_offset(target_site(bond)):site_offset(target_site(bond))+nb_fixture-1)
         i0 = site_offset(source_site(bond)); j0 = site_offset(target_site(bond))
         ref_block = collinear(i0:i0+nb_fixture-1,j0:j0+nb_fixture-1,bond)
         candidate_block = candidate(i0:i0+nb_fixture-1,j0:j0+nb_fixture-1,bond)
         expected = matmul(fs, matmul(ref_block, transpose(conjg(ft))))
         abs_error = max(abs_error, maxval(abs(candidate_block-expected)))
         rel_error = max(rel_error, norm_error(candidate_block, expected))
      end do
   end subroutine compare_block_sets_rotated

   subroutine compare_onsite_blocks_transformed(candidate, reference, frame, abs_error, rel_error)
      complex(rp), intent(in) :: candidate(:,:,:), reference(:,:,:), frame(:,:)
      real(rp), intent(out) :: abs_error, rel_error
      complex(rp) :: expected(nb_fixture,nb_fixture), ref_block(nb_fixture,nb_fixture)
      complex(rp) :: candidate_block(nb_fixture,nb_fixture), fs(nb_fixture,nb_fixture)
      integer :: site, i0
      abs_error = 0.0_rp
      rel_error = 0.0_rp
      do site = 1, nsite_fixture
         i0 = site_offset(site)
         fs = frame(i0:i0+nb_fixture-1,i0:i0+nb_fixture-1)
         ref_block = reference(i0:i0+nb_fixture-1,i0:i0+nb_fixture-1,site)
         candidate_block = candidate(i0:i0+nb_fixture-1,i0:i0+nb_fixture-1,site)
         expected = matmul(transpose(conjg(fs)), matmul(ref_block, fs))
         abs_error = max(abs_error, maxval(abs(candidate_block-expected)))
         rel_error = max(rel_error, norm_error(candidate_block, expected))
      end do
   end subroutine compare_onsite_blocks_transformed

   subroutine compare_kspace_sets(candidate, reference, frame, abs_error, rel_error, eigen_error)
      complex(rp), intent(in) :: candidate(:,:,:), reference(:,:,:), frame(:,:)
      real(rp), intent(out) :: abs_error, rel_error, eigen_error
      complex(rp) :: hc(nmat_fixture,nmat_fixture), hr(nmat_fixture,nmat_fixture), expected(nmat_fixture,nmat_fixture)
      integer :: ik
      abs_error = 0.0_rp
      rel_error = 0.0_rp
      eigen_error = 0.0_rp
      do ik = 1, nk_fixture
         call assemble_kspace(candidate, k_points(:,ik), hc)
         call assemble_kspace(reference, k_points(:,ik), hr)
         expected = matmul(transpose(conjg(frame)), matmul(hr, frame))
         abs_error = max(abs_error, maxval(abs(hc-expected)))
         rel_error = max(rel_error, norm_error(hc, expected))
         eigen_error = max(eigen_error, eigenvalue_error(hc, hr))
      end do
   end subroutine compare_kspace_sets

   subroutine assemble_kspace(blocks, k, h)
      complex(rp), intent(in) :: blocks(:,:,:)
      real(rp), intent(in) :: k(3)
      complex(rp), intent(out) :: h(:,:)
      complex(rp) :: phase
      integer :: bond, i0, j0
      h = cmplx(0.0_rp, 0.0_rp, rp)
      do bond = 1, nbond_fixture
         phase = cmplx(cos(2.0_rp*pi*dot_product(k, translation(:,bond))), &
                       sin(2.0_rp*pi*dot_product(k, translation(:,bond))), rp)
         i0 = site_offset(source_site(bond)); j0 = site_offset(target_site(bond))
         h(i0:i0+nb_fixture-1,j0:j0+nb_fixture-1) = &
            h(i0:i0+nb_fixture-1,j0:j0+nb_fixture-1) + &
            blocks(i0:i0+nb_fixture-1,j0:j0+nb_fixture-1,bond)*phase
      end do
   end subroutine assemble_kspace

   function frame_identity() result(frame)
      complex(rp) :: frame(nmat_fixture,nmat_fixture)
      frame = cmplx(0.0_rp, 0.0_rp, rp)
      frame = frame + identity_matrix(nmat_fixture)
   end function frame_identity

   subroutine make_packed_frame(theta_in, phi_in, frame)
      real(rp), intent(in) :: theta_in(:), phi_in(:)
      complex(rp), intent(out) :: frame(:,:)
      integer :: site
      frame = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, nsite_fixture
         frame(site_offset(site):site_offset(site)+nb_fixture-1, &
              site_offset(site):site_offset(site)+nb_fixture-1) = site_frame(theta_in(site), phi_in(site))
      end do
   end subroutine make_packed_frame

   function site_frame(theta_in, phi_in) result(frame)
      real(rp), intent(in) :: theta_in, phi_in
      complex(rp) :: frame(nb_fixture,nb_fixture), spin_frame(2,2)
      complex(rp) :: eiphi
      real(rp) :: c, s
      integer :: orbital, spin_a, spin_b
      c = cos(0.5_rp*theta_in); s = sin(0.5_rp*theta_in)
      eiphi = cmplx(cos(phi_in), sin(phi_in), rp)
      spin_frame(1,1) = cmplx(c, 0.0_rp, rp)
      spin_frame(1,2) = -conjg(eiphi)*s
      spin_frame(2,1) = eiphi*s
      spin_frame(2,2) = cmplx(c, 0.0_rp, rp)
      frame = cmplx(0.0_rp, 0.0_rp, rp)
      do spin_a = 1, 2
         do spin_b = 1, 2
            do orbital = 1, norb_fixture
               frame(orbital+(spin_a-1)*norb_fixture, orbital+(spin_b-1)*norb_fixture) = spin_frame(spin_a,spin_b)
            end do
         end do
      end do
   end function site_frame

   subroutine independent_moment(theta_in, phi_in, moment)
      real(rp), intent(in) :: theta_in, phi_in
      real(rp), intent(out) :: moment(3)
      moment = [sin(theta_in)*cos(phi_in), sin(theta_in)*sin(phi_in), cos(theta_in)]
   end subroutine independent_moment

   pure integer function site_offset(site) result(offset)
      integer, intent(in) :: site
      offset = (site-1)*nb_fixture + 1
   end function site_offset

   pure function identity_matrix(n) result(identity)
      integer, intent(in) :: n
      complex(rp) :: identity(n,n)
      integer :: i
      identity = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, n
         identity(i,i) = cmplx(1.0_rp, 0.0_rp, rp)
      end do
   end function identity_matrix

   pure function norm_error(actual, reference) result(error)
      complex(rp), intent(in) :: actual(:,:), reference(:,:)
      real(rp) :: error, denominator
      denominator = max(1.0_rp, sqrt(sum(abs(reference)**2)))
      error = sqrt(sum(abs(actual-reference)**2))/denominator
   end function norm_error

   pure function norm_error_3d(actual, reference) result(error)
      complex(rp), intent(in) :: actual(:,:,:), reference(:,:,:)
      real(rp) :: error, denominator
      denominator = max(1.0_rp, sqrt(sum(abs(reference)**2)))
      error = sqrt(sum(abs(actual-reference)**2))/denominator
   end function norm_error_3d

   function eigenvalue_error(first, second) result(error)
      complex(rp), intent(in) :: first(:,:), second(:,:)
      complex(rp) :: work(4*nmat_fixture), first_copy(nmat_fixture,nmat_fixture), second_copy(nmat_fixture,nmat_fixture)
      real(rp) :: first_eigen(nmat_fixture), second_eigen(nmat_fixture), rwork(3*nmat_fixture-2)
      real(rp) :: error
      integer :: info, lwork
      external zheev
      first_copy = first; second_copy = second
      lwork = size(work)
      call zheev('N', 'U', nmat_fixture, first_copy, nmat_fixture, first_eigen, work, lwork, rwork, info)
      if (info /= 0) then
         failed = .true.; error = huge(1.0_rp); return
      end if
      call zheev('N', 'U', nmat_fixture, second_copy, nmat_fixture, second_eigen, work, lwork, rwork, info)
      if (info /= 0) then
         failed = .true.; error = huge(1.0_rp); return
      end if
      error = maxval(abs(first_eigen-second_eigen))
   end function eigenvalue_error

   subroutine report(label, abs_error, rel_error, eigen_error)
      character(len=*), intent(in) :: label
      real(rp), intent(in) :: abs_error, rel_error, eigen_error
      write (*, '(a,3(es12.4,1x))') trim(label)//' (abs, rel, eigen): ', abs_error, rel_error, eigen_error
      largest_abs = max(largest_abs, abs_error)
      largest_rel = max(largest_rel, rel_error)
      largest_eigen = max(largest_eigen, eigen_error)
      if (rel_error > tol .or. eigen_error > tol) then
         write (*, '(a,a)') 'FAIL ', trim(label)
         failed = .true.
      end if
   end subroutine report

end program test_gbt_wp01_q0
