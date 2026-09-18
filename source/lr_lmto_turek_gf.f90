!------------------------------------------------------------------------------
! DRESP-03TG -- native screened LMTO/Turek path-operator backend.
!
! This module owns only the representation-specific algebra.  In particular,
! it never obtains a native path operator by scaling a finite-H resolvent.
! The latter remains an independent oracle for the DRESP-03TG tests.
!------------------------------------------------------------------------------
module lr_lmto_turek_gf_mod
   use lattice_mod, only: lattice
   use symbolic_atom_mod, only: symbolic_atom
   use precision_mod, only: rp
   use math_mod, only: i_unit, pi, inverse_3x3
   implicit none
   private

   real(rp), parameter :: legacy_alpha(0:3) = [0.3485_rp, 0.05303_rp, 0.010714_rp, 0.00337_rp]

   public :: native_complex_p_matrix
   public :: native_screening_alpha
   public :: native_screened_p_matrix
   public :: native_delta_p
   public :: native_structure_constants
   public :: native_path_operator
   public :: native_path_operator_from_structure
   public :: native_inverse
   public :: native_exchange_trace
   public :: native_exchange_integrand
   public :: native_finite_h_integrand
   public :: native_collinear_pauli_integrand

contains

   !> Invert a complex matrix through a LAPACK solve against the identity.
   subroutine native_inverse(matrix, inverse)
      complex(rp), intent(in) :: matrix(:, :)
      complex(rp), intent(out) :: inverse(:, :)
      complex(rp), allocatable :: a(:,:), rhs(:,:)
      integer, allocatable :: ipiv(:)
      integer :: n, info, i
      external :: zgesv

      n = size(matrix,1)
      if (size(matrix,2) /= n .or. size(inverse,1) /= n .or. size(inverse,2) /= n) then
         error stop 'native_inverse: invalid matrix shape'
      end if
      allocate(a(n,n),rhs(n,n),ipiv(n))
      a = matrix
      rhs = cmplx(0.0_rp,0.0_rp,rp)
      do i=1,n
         rhs(i,i)=cmplx(1.0_rp,0.0_rp,rp)
      end do
      call zgesv(n,n,a,n,ipiv,rhs,n,info)
      if (info /= 0) error stop 'native_inverse: complex solve failed'
      inverse = rhs
      deallocate(a,rhs,ipiv)
   end subroutine native_inverse

   !> Complex continuation of symbolic_atom%p_matrix.
   !> The active state is read without modifying the symbolic atom.  The
   !> post-predls `dele` width makes this the normalized P representation;
   !> its matching screening gamma is `potential%qi`, not raw `qpar`.
   !> The ordering is the live one: all up orbitals followed by all down orbitals.
   subroutine native_complex_p_matrix(atom, z, pmat)
      type(symbolic_atom), intent(in) :: atom
      complex(rp), intent(in) :: z
      complex(rp), intent(out) :: pmat(:, :)
      integer :: l, m, mls, norb, lmax, spin
      complex(rp) :: center, width

      lmax = atom%potential%lmax
      norb = (lmax + 1)**2
      if (size(pmat, 1) /= 2*norb .or. size(pmat, 2) /= 2*norb) then
         error stop 'native_complex_p_matrix: invalid matrix shape'
      end if
      pmat = cmplx(0.0_rp, 0.0_rp, rp)
      do l = 0, lmax
         do m = 1, 2*l + 1
            mls = l*l + m
            do spin = 1, 2
               center = cmplx(atom%potential%c(l,spin) + atom%potential%vmad, 0.0_rp, rp)
               width = cmplx(atom%potential%dele(l,spin), 0.0_rp, rp)
               pmat(mls + (spin-1)*norb, mls + (spin-1)*norb) = (z-center)/(width*width)
            end do
         end do
      end do
   end subroutine native_complex_p_matrix

   !> Return the spin-independent target screening used by production sbar.
   !> A stored alpha is authoritative.  The legacy fallback is the explicit
   !> LMTO47 value identified by the representation audit.
   subroutine native_screening_alpha(atom, alpha)
      type(symbolic_atom), intent(in) :: atom
      real(rp), intent(out) :: alpha(0:)
      integer :: lmax, ncopy

      lmax = atom%potential%lmax
      if (size(alpha) < lmax + 1) error stop 'native_screening_alpha: alpha array too small'
      alpha = 0.0_rp
      if (allocated(atom%potential%screening_alpha)) then
         ncopy = min(size(alpha), size(atom%potential%screening_alpha))
         alpha(0:ncopy-1) = atom%potential%screening_alpha(lbound(atom%potential%screening_alpha,1): &
            lbound(atom%potential%screening_alpha,1)+ncopy-1)
      else
         alpha(0:min(lmax,ubound(legacy_alpha,1))) = legacy_alpha(0:min(lmax,ubound(legacy_alpha,1)))
         if (lmax > ubound(legacy_alpha,1)) alpha(ubound(legacy_alpha,1)+1:lmax) = legacy_alpha(ubound(legacy_alpha,1))
      end if
   end subroutine native_screening_alpha

   !> Apply the live Turek screening transform to normalized P^gamma.
   !> `native_complex_p_matrix` uses `potential%dele`, so the input gamma is
   !> the normalized `potential%qi`.  Raw `potential%qpar` belongs to the
   !> pre-predls `srdel` representation and is not interchangeable here.
   subroutine native_screened_p_matrix(atom, z, alpha_out, pmat)
      type(symbolic_atom), intent(in) :: atom
      complex(rp), intent(in) :: z
      real(rp), intent(in) :: alpha_out(0:)
      complex(rp), intent(out) :: pmat(:, :)
      complex(rp) :: raw(size(pmat,1),size(pmat,2))
      complex(rp) :: pin, ain, aout
      integer :: l, m, mls, norb, lmax, spin

      lmax = atom%potential%lmax
      norb = (lmax + 1)**2
      if (size(alpha_out) < lmax + 1) error stop 'native_screened_p_matrix: alpha array too small'
      call native_complex_p_matrix(atom, z, raw)
      pmat = cmplx(0.0_rp, 0.0_rp, rp)
      do spin = 1, 2
         do l = 0, lmax
            ain = cmplx(atom%potential%qi(l,spin), 0.0_rp, rp)
            aout = cmplx(alpha_out(l), 0.0_rp, rp)
            do m = 1, 2*l + 1
               mls = l*l + m + (spin-1)*norb
               pin = raw(mls,mls)
               pmat(mls,mls) = pin/(1.0_rp + (ain-aout)*pin)
            end do
         end do
      end do
   end subroutine native_screened_p_matrix

   !> Extract the full-orbital spin vertex Delta P = P_up - P_down.
   subroutine native_delta_p(pmat, delta)
      complex(rp), intent(in) :: pmat(:, :)
      complex(rp), intent(out) :: delta(:, :)
      integer :: norb
      norb = size(pmat,1)/2
      if (size(pmat,2) /= 2*norb .or. size(delta,1) /= norb .or. size(delta,2) /= norb) then
         error stop 'native_delta_p: invalid matrix shape'
      end if
      delta = pmat(1:norb,1:norb) - pmat(norb+1:2*norb,norb+1:2*norb)
   end subroutine native_delta_p

   !> Fourier transform the production screened structure constants.
   !> The transpose is intentional: hamiltonian_build consumes sbar(j,i) as
   !> the row-i/column-j LMTO block.  k is in the same direct fractional
   !> coordinates used by reciprocal%calculate_structure_factors.
   subroutine native_structure_constants(lat, k, smat)
      type(lattice), intent(inout) :: lat
      real(rp), intent(in) :: k(3)
      complex(rp), intent(out) :: smat(:, :)
      real(rp), allocatable :: cralat(:, :), ham_vec(:, :), kfrac(:)
      complex(rp) :: phase
      integer :: nsite, norb, ntype, ia, ja, nr, m, isite, jsite, ino
      integer :: i0, j0, nn_local

      nsite = lat%nrec
      if (nsite < 1) error stop 'native_structure_constants: lattice has no reciprocal sites'
      norb = (lat%symbolic_atoms(lat%iz(lat%atlist(1)))%potential%lmax + 1)**2
      if (size(smat,1) /= norb*nsite .or. size(smat,2) /= norb*nsite) then
         error stop 'native_structure_constants: invalid matrix shape'
      end if
      smat = cmplx(0.0_rp, 0.0_rp, rp)
      allocate(cralat(3,lat%kk), ham_vec(3,max(1,lat%nn_max)), kfrac(3))
      cralat = lat%cr(:,1:lat%kk)*lat%alat
      kfrac = k-floor(k+0.5_rp)
      do isite = 1, nsite
         ntype = lat%ib(isite)
         ia = lat%atlist(ntype)
         nr = lat%nn(ia,1)
         nn_local = nr
         call lat%clusba(lat%r2, cralat, ia, lat%kk, lat%kk, nn_local, ham_vec(:,1:max(1,nr)))
         ino = lat%num(ia)
         i0 = (isite-1)*norb
         do m = 1, nn_local
            if (m == 1) then
               jsite = isite
               phase = cmplx(1.0_rp,0.0_rp,rp)
            else
               ja = lat%nn(ia,m)
               if (ja < 1 .or. ja > lat%kk) cycle
               jsite = lat%iz(ja)
               if (jsite < 1 .or. jsite > nsite) cycle
               if (lat%a_cart_inv_ready) then
                  phase = exp(i_unit*2.0_rp*pi*dot_product(kfrac, &
                     matmul(lat%a_cart_inv, ham_vec(:,m)/lat%alat)))
               else
                  phase = exp(i_unit*2.0_rp*pi*dot_product(kfrac, &
                     matmul(inverse_3x3(lat%a), ham_vec(:,m)/lat%alat)))
               end if
            end if
            j0 = (jsite-1)*norb
            smat(i0+1:i0+norb,j0+1:j0+norb) = smat(i0+1:i0+norb,j0+1:j0+norb) + &
               transpose(lat%sbar(1:norb,1:norb,m,ino))*phase
         end do
      end do
      deallocate(cralat,ham_vec,kfrac)
   end subroutine native_structure_constants

   !> Construct the production native spin path operator directly.
   subroutine native_path_operator(lat, z, k, pmat, smat, gmat)
      type(lattice), intent(inout) :: lat
      complex(rp), intent(in) :: z
      real(rp), intent(in) :: k(3)
      complex(rp), intent(out) :: pmat(:, :), smat(:, :), gmat(:, :)
      complex(rp), allocatable :: s_orb(:,:)
      integer :: nsite, norb, nspin, isite, it, p0, q0, ioff

      nsite = lat%nrec
      it = lat%iz(lat%atlist(1))
      norb = (lat%symbolic_atoms(it)%potential%lmax + 1)**2
      nspin = 2*norb*nsite
      if (size(pmat,1) /= nspin .or. size(smat,1) /= nspin .or. size(gmat,1) /= nspin) then
         error stop 'native_path_operator: invalid matrix shape'
      end if
      allocate(s_orb(norb*nsite,norb*nsite))
      call native_structure_constants(lat,k,s_orb)
      smat = cmplx(0.0_rp,0.0_rp,rp)
      do isite = 1,nsite
         do q0 = 1,nsite
            p0 = (isite-1)*norb; ioff = (q0-1)*norb
            smat((isite-1)*2*norb+1:(isite-1)*2*norb+norb, (q0-1)*2*norb+1:(q0-1)*2*norb+norb) = &
               s_orb(p0+1:p0+norb,ioff+1:ioff+norb)
            smat((isite-1)*2*norb+norb+1:(isite-1)*2*norb+2*norb, (q0-1)*2*norb+norb+1:(q0-1)*2*norb+2*norb) = &
               s_orb(p0+1:p0+norb,ioff+1:ioff+norb)
         end do
      end do
      call native_path_operator_from_structure(lat,z,smat,pmat,gmat)
      deallocate(s_orb)
   end subroutine native_path_operator

   !> Solve the native P(z)-S path operator with a prebuilt spin structure
   !> matrix.  Contour integrations call this repeatedly for one fixed k or
   !> k+q; separating the structure build from the energy solve avoids
   !> rebuilding the real-space Fourier sum at every contour node.
   subroutine native_path_operator_from_structure(lat, z, smat, pmat, gmat)
      type(lattice), intent(inout) :: lat
      complex(rp), intent(in) :: z, smat(:, :)
      complex(rp), intent(out) :: pmat(:, :), gmat(:, :)
      complex(rp), allocatable :: psite(:,:), work(:,:), rhs(:,:)
      real(rp), allocatable :: alpha(:)
      integer :: nsite, norb, nspin, isite, ntype, ia, it, p0, info
      integer, allocatable :: ipiv(:)
      external :: zgesv

      nsite = lat%nrec
      it = lat%iz(lat%atlist(1))
      norb = (lat%symbolic_atoms(it)%potential%lmax + 1)**2
      nspin = 2*norb*nsite
      if (any(shape(smat) /= [nspin,nspin]) .or. any(shape(pmat) /= [nspin,nspin]) .or. &
          any(shape(gmat) /= [nspin,nspin])) then
         error stop 'native_path_operator_from_structure: invalid matrix shape'
      end if
      allocate(psite(2*norb,2*norb), alpha(0:lat%symbolic_atoms(it)%potential%lmax), &
         work(nspin,nspin),rhs(nspin,nspin),ipiv(nspin))
      pmat = cmplx(0.0_rp,0.0_rp,rp)
      do isite = 1,nsite
         ntype = lat%ib(isite); ia = lat%atlist(ntype); it = lat%iz(ia)
         if ((lat%symbolic_atoms(it)%potential%lmax+1)**2 /= norb) then
            error stop 'native_path_operator_from_structure: mixed lmax unsupported'
         end if
         call native_screening_alpha(lat%symbolic_atoms(it),alpha)
         call native_screened_p_matrix(lat%symbolic_atoms(it),z,alpha,psite)
         p0 = (isite-1)*2*norb
         pmat(p0+1:p0+2*norb,p0+1:p0+2*norb) = psite
      end do
      work = pmat-smat
      rhs = cmplx(0.0_rp,0.0_rp,rp)
      do info = 1,nspin
         rhs(info,info) = cmplx(1.0_rp,0.0_rp,rp)
      end do
      call zgesv(nspin,nspin,work,nspin,ipiv,rhs,nspin,info)
      if (info /= 0) error stop 'native_path_operator_from_structure: complex solve failed'
      gmat = rhs
      deallocate(psite,alpha,work,rhs,ipiv)
   end subroutine native_path_operator_from_structure

   pure function native_exchange_trace(delta_i, delta_j, gup_ij, gdown_ji) result(value)
      complex(rp), intent(in) :: delta_i(:,:), delta_j(:,:), gup_ij(:,:), gdown_ji(:,:)
      complex(rp) :: product(size(delta_i,1),size(delta_i,2)), value
      integer :: i
      product = matmul(delta_i,matmul(gup_ij,matmul(delta_j,gdown_ji)))
      value = cmplx(0.0_rp,0.0_rp,rp)
      do i = 1, size(product,1)
         value = value + product(i,i)
      end do
   end function native_exchange_trace

   pure function native_exchange_integrand(delta_i, delta_j, gup_ij, gdown_ji) result(value)
      complex(rp), intent(in) :: delta_i(:,:), delta_j(:,:), gup_ij(:,:), gdown_ji(:,:)
      real(rp) :: value
      value = aimag(native_exchange_trace(delta_i,delta_j,gup_ij,gdown_ji))
   end function native_exchange_integrand

   pure function native_finite_h_integrand(d_i,d_j,gup_ij,gdown_ji) result(value)
      complex(rp), intent(in) :: d_i(:,:), d_j(:,:), gup_ij(:,:), gdown_ji(:,:)
      complex(rp) :: product(size(d_i,1),size(d_i,2))
      real(rp) :: value
      integer :: i
      product = matmul(d_i,matmul(gup_ij,matmul(d_j,gdown_ji)))
      value = sum([(aimag(product(i,i)),i=1,size(product,1))])
   end function native_finite_h_integrand

   ! Historical dGdG_Jnc reduction in the live Pauli convention.  The
   ! separate routine is intentionally explicit so tests can compare it with
   ! the native up/down contraction without importing exchange.f90.
   pure function native_collinear_pauli_integrand(delta_i,delta_j,gij,gji) result(value)
      complex(rp), intent(in) :: delta_i(:,:),delta_j(:,:),gij(:,:),gji(:,:)
      complex(rp) :: g0i(size(delta_i,1),size(delta_i,2)), gzi(size(delta_i,1),size(delta_i,2))
      complex(rp) :: g0j(size(delta_i,1),size(delta_i,2)), gzj(size(delta_i,1),size(delta_i,2)), product(size(delta_i,1),size(delta_i,2))
      real(rp) :: value
      integer :: n, i
      n = size(delta_i,1)
      g0i = 0.5_rp*(gij(1:n,1:n)+gij(n+1:2*n,n+1:2*n))
      gzi = 0.5_rp*(gij(1:n,1:n)-gij(n+1:2*n,n+1:2*n))
      g0j = 0.5_rp*(gji(1:n,1:n)+gji(n+1:2*n,n+1:2*n))
      gzj = 0.5_rp*(gji(1:n,1:n)-gji(n+1:2*n,n+1:2*n))
      product = matmul(delta_i,matmul(g0i,matmul(delta_j,g0j))) - &
                matmul(delta_i,matmul(gzi,matmul(delta_j,gzj)))
      value = 0.0_rp
      do i=1,n
         value = value + aimag(product(i,i))
      end do
   end function native_collinear_pauli_integrand

end module lr_lmto_turek_gf_mod
