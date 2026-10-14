! Independent two-level oracle for the fixed-Liechtenstein KS and RS
! occupation producers. It checks conjugation order and complex coherence.
! Production's site slice is [site][spin block][orbital], then each l-channel
! is sliced at basis_offset=l*l with dimensions 2*l+1; the test's n=2 block
! is one such site/spin/l slice. Both production paths call the exact scalar
! kernels tested below; their site/spin indexing and MPI packing are reviewed
! at the call sites and exercised by the VAL-04 runtime campaign.
program test_ldau_occupation_producers
   use precision_mod, only: rp
   use math_mod, only: hubbard_projector_element, hubbard_spectral_density_element, hubbard_density_diagnostics
   implicit none
   integer, parameter :: n = 2, ne = 20000
   real(rp), parameter :: eta = 0.01_rp, emin = -10.0_rp, efermi = 0.0_rp
   real(rp), parameter :: pi = 3.1415926535897932384626433832795_rp
   complex(rp) :: psi(n, n), p_ks(n,n), p_rs(n,n), p_exact(n,n), g(n,n), a(n,n)
   complex(rp) :: diagnostic_matrix(n,n), diagnostic_before(n,n)
   real(rp) :: eps(n), h, e, weight, err, broadening_error, exact_weight
   real(rp) :: diag_trace, eig_min, eig_max, antiherm
   integer :: diagnostic_info
   integer :: k, i, j, band, failures

   failures = 0
   ! Orthonormal eigenvectors of a Hermitian Hamiltonian. The occupied state
   ! has purely imaginary off-diagonal coherence in the selected basis.
   psi(:,1) = [cmplx(1.0_rp,0.0_rp,rp), cmplx(0.0_rp,1.0_rp,rp)]/sqrt(2.0_rp)
   psi(:,2) = [cmplx(1.0_rp,0.0_rp,rp), cmplx(0.0_rp,-1.0_rp,rp)]/sqrt(2.0_rp)
   eps = [-1.0_rp, 1.0_rp]

   ! KS projector definition P_ab = sum_occ psi_a psi_b*.
   p_ks = cmplx(0.0_rp,0.0_rp,rp)
   do i = 1, n
      do j = 1, n
         p_ks(i,j) = hubbard_projector_element(psi(i,1), psi(j,1))
      end do
   end do
   call check(abs(p_ks(1,2)-cmplx(0.0_rp,-0.5_rp,rp)) < 1.0e-14_rp, 'KS imaginary coherence sign')
   call check(maxval(abs(p_ks-transpose(conjg(p_ks)))) < 1.0e-14_rp, 'KS Hermiticity')
   call check(abs(real(sum([(p_ks(i,i),i=1,n)]),rp)-1.0_rp) < 1.0e-14_rp, 'KS trace')

   ! RS spectral projector P = integral i(G-G^dagger)/(2pi) dE, with the
   ! same retarded Green function convention as the production engine.
   p_rs = cmplx(0.0_rp,0.0_rp,rp)
   p_exact = cmplx(0.0_rp,0.0_rp,rp)
   do k = 1, n
      exact_weight = (atan((efermi-eps(k))/eta)-atan((emin-eps(k))/eta))/pi
      do i = 1, n
         do j = 1, n
            p_exact(i,j) = p_exact(i,j) + exact_weight * &
               hubbard_projector_element(psi(i,k), psi(j,k))
         end do
      end do
   end do
   h = (efermi-emin)/real(ne,rp)
   do k = 0, ne
      e = emin + real(k,rp)*h
      g = cmplx(0.0_rp,0.0_rp,rp)
      do i = 1, n
         do j = 1, n
            g(i,j) = sum([(hubbard_projector_element(psi(i,band),psi(j,band))/cmplx(e-eps(band),eta,rp), band=1,n)])
         end do
      end do
      do i = 1, n
         do j = 1, n
            a(i,j) = hubbard_spectral_density_element(g(i,j),g(j,i))
         end do
      end do
      if (k == 0 .or. k == ne) then
         weight = 0.5_rp
      else
         weight = 1.0_rp
      end if
      p_rs = p_rs + weight*h*a
   end do

   err = maxval(abs(p_rs-p_exact))
   broadening_error = maxval(abs(p_exact-p_ks))
   call check(err < 1.0e-8_rp, 'RS quadrature agrees with analytic finite-window Lorentz integral')
   call check(broadening_error > 0.0_rp .and. broadening_error < 1.0e-2_rp, &
      'finite eta/window broadening error is separately bounded')
   call check(abs(aimag(p_rs(1,2))) > 0.49_rp, 'RS retains imaginary coherence')
   call check(maxval(abs(p_rs-transpose(conjg(p_rs)))) < 1.0e-12_rp, 'RS Hermiticity')

   ! Production diagnostics must accept a physical Hermitian PSD matrix and
   ! expose invalid Hermitian occupation bounds without repairing the input.
   diagnostic_matrix = reshape([cmplx(0.5_rp,0.0_rp,rp), cmplx(0.0_rp,0.3_rp,rp), &
                                cmplx(0.0_rp,-0.3_rp,rp), cmplx(0.5_rp,0.0_rp,rp)], [n,n])
   call hubbard_density_diagnostics(diagnostic_matrix, diag_trace, eig_min, eig_max, antiherm, diagnostic_info)
   call check(diagnostic_info == 0 .and. abs(diag_trace-1.0_rp) < 1.0e-14_rp, 'physical diagnostic trace')
   call check(abs(eig_min-0.2_rp) < 1.0e-14_rp .and. abs(eig_max-0.8_rp) < 1.0e-14_rp, &
      'physical PSD diagnostic eigenvalue bounds')
   call check(antiherm < 1.0e-14_rp, 'physical diagnostic Hermiticity')
   diagnostic_matrix = reshape([cmplx(1.01_rp,0.0_rp,rp), cmplx(0.0_rp,0.0_rp,rp), &
                                cmplx(0.0_rp,0.0_rp,rp), cmplx(-0.01_rp,0.0_rp,rp)], [n,n])
   diagnostic_before = diagnostic_matrix
   call hubbard_density_diagnostics(diagnostic_matrix, diag_trace, eig_min, eig_max, antiherm, diagnostic_info)
   call check(diagnostic_info == 0 .and. (eig_min < 0.0_rp .or. eig_max > 1.0_rp), &
      'invalid Hermitian occupation bounds are reported')
   call check(maxval(abs(diagnostic_matrix-diagnostic_before)) < 1.0e-14_rp, &
      'diagnostics leave input matrix unchanged')
   if (failures == 0) then
      write(*,'(a,es12.4,a,es12.4)') 'UnitLdaUOccupationProducers: ALL TESTS PASSED; RS quadrature error=', &
         err, ', finite-eta/window difference=', broadening_error
      stop 0
   else
      write(*,'(a,i0,a)') 'UnitLdaUOccupationProducers: ', failures, ' TEST(S) FAILED'
      stop 1
   end if

contains
   subroutine check(ok, label)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: label
      if (.not. ok) then
         write(*,'(a,a)') 'FAIL ', label
         failures = failures + 1
      end if
   end subroutine check
end program test_ldau_occupation_producers
