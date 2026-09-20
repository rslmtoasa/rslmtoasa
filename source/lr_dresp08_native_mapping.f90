!------------------------------------------------------------------------------
! DRESP-08 native second-order LMTO field mapping.
!
! The routines in this module are deliberately algebraic.  They do not fit a
! coefficient-space field and they do not enter Kxc, BES/Halle, or a Dyson
! solve.  The material seam is the same one used by hamiltonian%build_bulkham:
!
!   structure block -> wx/cex bond -> ee, eeo, enim
!   Fourier transform -> h(k), h(k)o, E_nu
!   H^(2)(k) = E_nu + h(k) - h(k)o h(k)
!------------------------------------------------------------------------------
module lr_dresp08_native_mapping_mod

   use precision_mod, only: rp
   use basis_mod, only: nb, norb, spin_off
   use math_mod, only: hcpx
   use hamiltonian_mod, only: hamiltonian
   use lmto_magnetic_tangent_mod, only: lmto_bond_derivative, lmto_hhmag_to_spinor
   implicit none
   private

   public :: dresp08_build_h2
   public :: dresp08_build_product_tangent
   public :: dresp08_commutator_tangent
   public :: dresp08_build_spinor_tangent
   public :: dresp08_build_spinor_operator
   public :: dresp08_build_native_realspace_tangent
   public :: dresp08_build_native_onsite_tangents
   public :: dresp08_block_diagonal
   public :: dresp08_extract_spin_blocks
   public :: dresp08_rotate_collinear_operator
   public :: dresp08_relative_residual

contains

   !> Exact live second-order production algebra for an already Fourier
   !> transformed h(k), overlap-bar o, and onsite E_nu.
   subroutine dresp08_build_h2(h, o, enu, h2, hoh)
      complex(rp), intent(in) :: h(:, :), o(:, :), enu(:, :)
      complex(rp), intent(out) :: h2(:, :)
      complex(rp), intent(out), optional :: hoh(:, :)
      complex(rp), allocatable :: work(:, :)
      integer :: n

      n = size(h, 1)
      if (size(h, 2) /= n .or. any(shape(o) /= [n, n]) .or. any(shape(enu) /= [n, n]) .or. &
          any(shape(h2) /= [n, n])) error stop 'DRESP-08 H2: inconsistent matrix shapes'
      allocate(work(n, n))
      work = matmul(h, o)
      work = matmul(work, h)
      h2 = enu + h - work
      if (present(hoh)) then
         if (any(shape(hoh) /= [n, n])) error stop 'DRESP-08 H2: invalid HOH output shape'
         hoh = work
      end if
      deallocate(work)
   end subroutine dresp08_build_h2

   !> Differentiate H^(2)=E+h-hoh term by term.  The three HOH outputs retain
   !> their signs-free products so the audit can report each native source.
   subroutine dresp08_build_product_tangent(d_enu, d_h, d_o, h, o, d_h2, term_enu, term_h, &
                                             term_dh_o_h, term_h_do_h, term_h_o_dh)
      complex(rp), intent(in) :: d_enu(:, :), d_h(:, :), d_o(:, :), h(:, :), o(:, :)
      complex(rp), intent(out) :: d_h2(:, :), term_enu(:, :), term_h(:, :), term_dh_o_h(:, :), &
         term_h_do_h(:, :), term_h_o_dh(:, :)
      integer :: n

      n = size(h, 1)
      if (size(h, 2) /= n .or. any(shape(o) /= [n, n]) .or. any(shape(d_h) /= [n, n]) .or. &
          any(shape(d_o) /= [n, n]) .or. any(shape(d_enu) /= [n, n]) .or. any(shape(d_h2) /= [n, n])) then
         error stop 'DRESP-08 product tangent: inconsistent matrix shapes'
      end if
      if (any(shape(term_enu) /= [n, n]) .or. any(shape(term_h) /= [n, n]) .or. &
          any(shape(term_dh_o_h) /= [n, n]) .or. any(shape(term_h_do_h) /= [n, n]) .or. &
          any(shape(term_h_o_dh) /= [n, n])) error stop 'DRESP-08 product tangent: invalid output shapes'

      term_enu = d_enu
      term_h = d_h
      term_dh_o_h = matmul(matmul(d_h, o), h)
      term_h_do_h = matmul(matmul(h, d_o), h)
      term_h_o_dh = matmul(matmul(h, o), d_h)
      d_h2 = term_enu + term_h - term_dh_o_h - term_h_do_h - term_h_o_dh
   end subroutine dresp08_build_product_tangent

   pure subroutine dresp08_commutator_tangent(generator, matrix, delta)
      complex(rp), intent(in) :: generator(:, :), matrix(:, :)
      complex(rp), intent(out) :: delta(:, :)
      integer :: n

      n = size(matrix, 1)
      if (size(matrix, 2) /= n .or. any(shape(generator) /= [n, n]) .or. any(shape(delta) /= [n, n])) then
         error stop 'DRESP-08 commutator: inconsistent matrix shapes'
      end if
      delta = cmplx(0.0_rp, -1.0_rp, rp)*(matmul(generator, matrix) - matmul(matrix, generator))
   end subroutine dresp08_commutator_tangent

   !> Build the transverse variation of a diagonal-in-spin local magnetic
   !> coefficient.  The input is the spin-vector coefficient X_1.  The
   !> cartesian spin blocks are transformed with the same hcpx seam used by
   !> build_obarm/build_enim and chbar_nc.
   subroutine dresp08_build_spinor_tangent(vector_coefficient, delta_matrix)
      complex(rp), intent(in) :: vector_coefficient(:)
      complex(rp), intent(out) :: delta_matrix(:, :)
      integer :: i

      if (size(vector_coefficient) /= norb .or. any(shape(delta_matrix) /= [nb, nb])) then
         error stop 'DRESP-08 spinor tangent: inconsistent dimensions'
      end if
      delta_matrix = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, norb
         delta_matrix(i, i + spin_off) = vector_coefficient(i)
         delta_matrix(i + spin_off, i) = vector_coefficient(i)
      end do
      call hcpx(delta_matrix(1:norb, 1:norb), 'cart2sph')
      call hcpx(delta_matrix(norb+1:nb, norb+1:nb), 'cart2sph')
      call hcpx(delta_matrix(1:norb, norb+1:nb), 'cart2sph')
      call hcpx(delta_matrix(norb+1:nb, 1:norb), 'cart2sph')
   end subroutine dresp08_build_spinor_tangent

   !> Reproduce the live collinear onsite construction from scalar/vector
   !> transformed potential parameters and a local moment direction.
   subroutine dresp08_build_spinor_operator(scalar_coefficient, vector_coefficient, moment, matrix)
      complex(rp), intent(in) :: scalar_coefficient(:), vector_coefficient(:)
      real(rp), intent(in) :: moment(3)
      complex(rp), intent(out) :: matrix(:, :)
      complex(rp), parameter :: i_unit = cmplx(0.0_rp, 1.0_rp, rp)
      integer :: i

      if (size(scalar_coefficient) /= norb .or. size(vector_coefficient) /= norb .or. &
          any(shape(matrix) /= [nb, nb])) error stop 'DRESP-08 spinor operator: inconsistent dimensions'
      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, norb
         matrix(i, i) = scalar_coefficient(i) + vector_coefficient(i)*moment(3)
         matrix(i + spin_off, i + spin_off) = scalar_coefficient(i) - vector_coefficient(i)*moment(3)
         matrix(i, i + spin_off) = vector_coefficient(i)*(moment(1) - i_unit*moment(2))
         matrix(i + spin_off, i) = vector_coefficient(i)*(moment(1) + i_unit*moment(2))
      end do
      call hcpx(matrix(1:norb, 1:norb), 'cart2sph')
      call hcpx(matrix(norb+1:nb, norb+1:nb), 'cart2sph')
      call hcpx(matrix(1:norb, norb+1:nb), 'cart2sph')
      call hcpx(matrix(norb+1:nb, 1:norb), 'cart2sph')
   end subroutine dresp08_build_spinor_operator

   !> Rebuild the complete real-space first-order transverse tangent from the
   !> stored screened structure constants and the live endpoint parameters.
   !> This is the production chbar_nc/ham0m_nc algebra, with dm= y x m for a
   !> positive global rotation about y.
   subroutine dresp08_build_native_realspace_tangent(hamiltonian_obj, d_ee)
      type(hamiltonian), intent(in) :: hamiltonian_obj
      complex(rp), intent(out) :: d_ee(:, :, :, :)
      complex(rp) :: hhhc(norb, norb), dhmag(norb, norb, 4)
      real(rp) :: mom_i(3), mom_j(3), dm_i(3), dm_j(3)
      integer :: ntype, ia, nr, m, ja, it, jt, ino, ilm, jlm, idir
      logical :: onsite

      if (.not. allocated(hamiltonian_obj%ee) .or. any(shape(d_ee) /= shape(hamiltonian_obj%ee))) then
         error stop 'DRESP-08 native tangent: ee storage is incomplete'
      end if
      if (.not. associated(hamiltonian_obj%charge) .or. .not. associated(hamiltonian_obj%lattice)) then
         error stop 'DRESP-08 native tangent: Hamiltonian ownership is incomplete'
      end if
      if (size(d_ee, 1) /= nb .or. size(d_ee, 2) /= nb) error stop 'DRESP-08 native tangent: invalid block size'
      d_ee = cmplx(0.0_rp, 0.0_rp, rp)

      do ntype = 1, hamiltonian_obj%charge%lattice%ntype
         ia = hamiltonian_obj%charge%lattice%atlist(ntype)
         nr = hamiltonian_obj%charge%lattice%nn(ia, 1)
         ino = hamiltonian_obj%charge%lattice%num(ia)
         it = hamiltonian_obj%charge%lattice%iz(ia)
         mom_i = real(hamiltonian_obj%charge%lattice%symbolic_atoms(it)%potential%mom(:), rp)
         dm_i = [mom_i(3), 0.0_rp, -mom_i(1)]
         do m = 1, nr
            if (m == 1) then
               ja = ia
            else
               ja = hamiltonian_obj%charge%lattice%nn(ia, m)
            end if
            if (ja < 1 .or. ja > hamiltonian_obj%charge%lattice%kk) cycle
            jt = hamiltonian_obj%charge%lattice%iz(ja)
            mom_j = real(hamiltonian_obj%charge%lattice%symbolic_atoms(jt)%potential%mom(:), rp)
            dm_j = [mom_j(3), 0.0_rp, -mom_j(1)]
            hhhc = cmplx(0.0_rp, 0.0_rp, rp)
            do jlm = 1, norb
               do ilm = 1, norb
                  ! hmfind uses hhh(ilm,jlm)=real(sbar(jlm,ilm,...)).
                  hhhc(ilm, jlm) = cmplx(real(hamiltonian_obj%charge%lattice%sbar(jlm, ilm, m, ino), rp), 0.0_rp, rp)
               end do
            end do
            onsite = m == 1
            if (hamiltonian_obj%hoh) then
               call lmto_bond_derivative(hhhc, &
                  hamiltonian_obj%charge%lattice%symbolic_atoms(it)%potential%wx0(1:norb), &
                  hamiltonian_obj%charge%lattice%symbolic_atoms(it)%potential%wx1(1:norb), &
                  hamiltonian_obj%charge%lattice%symbolic_atoms(jt)%potential%wx0(1:norb), &
                  hamiltonian_obj%charge%lattice%symbolic_atoms(jt)%potential%wx1(1:norb), &
                  hamiltonian_obj%charge%lattice%symbolic_atoms(it)%potential%cex1(1:norb), &
                  mom_i, mom_j, dm_i, dm_j, onsite, dhmag)
            else
               call lmto_bond_derivative(hhhc, &
                  hamiltonian_obj%charge%lattice%symbolic_atoms(it)%potential%wx0(1:norb), &
                  hamiltonian_obj%charge%lattice%symbolic_atoms(it)%potential%wx1(1:norb), &
                  hamiltonian_obj%charge%lattice%symbolic_atoms(jt)%potential%wx0(1:norb), &
                  hamiltonian_obj%charge%lattice%symbolic_atoms(jt)%potential%wx1(1:norb), &
                  hamiltonian_obj%charge%lattice%symbolic_atoms(it)%potential%cx1(1:norb), &
                  mom_i, mom_j, dm_i, dm_j, onsite, dhmag)
            end if
            do idir = 1, 4
               call hcpx(dhmag(:, :, idir), 'cart2sph')
            end do
            call lmto_hhmag_to_spinor(dhmag, d_ee(:, :, m, ntype))
         end do
      end do
   end subroutine dresp08_build_native_realspace_tangent

   !> Build direct onsite transverse tangents from the live obx1 and
   !> (cx1-cex1) magnetic components.  The returned arrays retain the live
   !> atom-type indexing used by hamiltonian%obarm/enim.
   subroutine dresp08_build_native_onsite_tangents(hamiltonian_obj, d_obarm, d_enim)
      type(hamiltonian), intent(in) :: hamiltonian_obj
      complex(rp), intent(out) :: d_obarm(:, :, :), d_enim(:, :, :)
      complex(rp) :: coefficient(norb)
      integer :: ntype

      if (.not. allocated(hamiltonian_obj%obarm) .or. .not. allocated(hamiltonian_obj%enim)) then
         error stop 'DRESP-08 onsite tangent: live obarm/enim arrays are incomplete'
      end if
      if (any(shape(d_obarm) /= shape(hamiltonian_obj%obarm)) .or. any(shape(d_enim) /= shape(hamiltonian_obj%enim))) then
         error stop 'DRESP-08 onsite tangent: output shapes differ from live matrices'
      end if
      d_obarm = cmplx(0.0_rp, 0.0_rp, rp)
      d_enim = cmplx(0.0_rp, 0.0_rp, rp)
      do ntype = 1, hamiltonian_obj%charge%lattice%ntype
         coefficient = hamiltonian_obj%charge%lattice%symbolic_atoms(ntype)%potential%obx1(1:norb)
         call dresp08_build_spinor_tangent(coefficient, d_obarm(:, :, ntype))
         coefficient = hamiltonian_obj%charge%lattice%symbolic_atoms(ntype)%potential%cx1(1:norb) - &
            hamiltonian_obj%charge%lattice%symbolic_atoms(ntype)%potential%cex1(1:norb)
         call dresp08_build_spinor_tangent(coefficient, d_enim(:, :, ntype))
      end do
   end subroutine dresp08_build_native_onsite_tangents

   subroutine dresp08_block_diagonal(type_blocks, site_types, matrix)
      complex(rp), intent(in) :: type_blocks(:, :, :)
      integer, intent(in) :: site_types(:)
      complex(rp), intent(out) :: matrix(:, :)
      integer :: nblock, isite, first, last, itype

      nblock = size(type_blocks, 1)
      if (size(type_blocks, 2) /= nblock .or. any(shape(matrix) /= [nblock*size(site_types), nblock*size(site_types)])) then
         error stop 'DRESP-08 block diagonal: inconsistent shapes'
      end if
      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      do isite = 1, size(site_types)
         itype = site_types(isite)
         if (itype < 1 .or. itype > size(type_blocks, 3)) error stop 'DRESP-08 block diagonal: invalid site type'
         first = (isite - 1)*nblock + 1
         last = isite*nblock
         matrix(first:last, first:last) = type_blocks(:, :, itype)
      end do
   end subroutine dresp08_block_diagonal

   subroutine dresp08_extract_spin_blocks(matrix, norb_site, nsite, up, down)
      complex(rp), intent(in) :: matrix(:, :)
      integer, intent(in) :: norb_site, nsite
      complex(rp), intent(out) :: up(:, :), down(:, :)
      integer :: isite, jsite, first_i, first_j, ui, di, uj, dj

      if (any(shape(matrix) /= [2*norb_site*nsite, 2*norb_site*nsite]) .or. &
          any(shape(up) /= [norb_site*nsite, norb_site*nsite]) .or. &
          any(shape(down) /= [norb_site*nsite, norb_site*nsite])) error stop 'DRESP-08 spin blocks: inconsistent shapes'
      up = cmplx(0.0_rp, 0.0_rp, rp)
      down = cmplx(0.0_rp, 0.0_rp, rp)
      do isite = 1, nsite
         first_i = (isite - 1)*norb_site + 1
         do jsite = 1, nsite
            first_j = (jsite - 1)*norb_site + 1
            ui = (isite - 1)*2*norb_site + 1
            di = ui + norb_site
            uj = (jsite - 1)*2*norb_site + 1
            dj = uj + norb_site
            up(first_i:first_i+norb_site-1, first_j:first_j+norb_site-1) = matrix(ui:ui+norb_site-1, uj:uj+norb_site-1)
            down(first_i:first_i+norb_site-1, first_j:first_j+norb_site-1) = matrix(di:di+norb_site-1, dj:dj+norb_site-1)
         end do
      end do
   end subroutine dresp08_extract_spin_blocks

   !> Rotate a collinear spin operator by theta about y.  This independent
   !> finite-rotation helper is used only by the DRESP-08 fixture oracle.
   subroutine dresp08_rotate_collinear_operator(up, down, theta, rotated)
      complex(rp), intent(in) :: up(:, :), down(:, :)
      real(rp), intent(in) :: theta
      complex(rp), intent(out) :: rotated(:, :)
      integer :: n
      real(rp) :: c, s

      n = size(up, 1)
      if (size(up, 2) /= n .or. any(shape(down) /= [n, n]) .or. any(shape(rotated) /= [2*n, 2*n])) then
         error stop 'DRESP-08 finite rotation: inconsistent matrix shapes'
      end if
      c = cos(0.5_rp*theta)
      s = sin(0.5_rp*theta)
      rotated = cmplx(0.0_rp, 0.0_rp, rp)
      rotated(1:n, 1:n) = c*c*up + s*s*down
      rotated(n+1:2*n, n+1:2*n) = s*s*up + c*c*down
      rotated(1:n, n+1:2*n) = c*s*(up - down)
      rotated(n+1:2*n, 1:n) = c*s*(up - down)
   end subroutine dresp08_rotate_collinear_operator

   real(rp) function dresp08_relative_residual(left, right) result(value)
      complex(rp), intent(in) :: left(:, :), right(:, :)
      real(rp) :: denominator
      if (any(shape(left) /= shape(right))) error stop 'DRESP-08 residual: shape mismatch'
      denominator = max(sqrt(sum(abs(right)**2)), tiny(1.0_rp))
      value = sqrt(sum(abs(left - right)**2))/denominator
   end function dresp08_relative_residual

end module lr_dresp08_native_mapping_mod
