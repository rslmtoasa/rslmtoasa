! DRESP-10 exact static Frechet/divided-difference regression.
program test_lr_static_frechet_product_response
   use precision_mod, only: rp
   use lr_exact_ks_ward_mod, only: exact_ks_fermi_occupation, exact_ks_fermi_derivative
   use lr_static_frechet_product_response_mod, only: lr_static_divided_difference, lr_static_frechet_density
   implicit none

   real(rp), parameter :: fermi = 0.0_rp, temperature = 300.0_rp, tol = 2.0e-12_rp
   real(rp) :: e(3), coefficient, expected
   complex(rp) :: vectors(3,3), delta_h(3,3), delta_rho(3,3)
   integer :: i

   e = [-0.15_rp, 0.35_rp, 0.35_rp + 5.0e-13_rp]
   vectors = cmplx(0.0_rp, 0.0_rp, rp)
   do i = 1, 3
      vectors(i,i) = cmplx(1.0_rp, 0.0_rp, rp)
   end do
   delta_h = cmplx(0.0_rp, 0.0_rp, rp)
   delta_h(1,2) = cmplx(0.7_rp, -0.2_rp, rp)
   delta_h(2,1) = conjg(delta_h(1,2))
   delta_h(2,3) = cmplx(-0.4_rp, 0.1_rp, rp)
   delta_h(3,2) = conjg(delta_h(2,3))

   coefficient = lr_static_divided_difference(e(2), e(3), fermi, temperature)
   expected = exact_ks_fermi_derivative(e(2), fermi, temperature)
   if (abs(coefficient - expected) > tol*max(1.0_rp, abs(expected))) error stop 'DRESP-10 Frechet degeneracy limit failed'

   call lr_static_frechet_density(e, vectors, fermi, temperature, delta_h, delta_rho)
   if (abs(delta_rho(1,2) - lr_static_divided_difference(e(1), e(2), fermi, temperature)*delta_h(1,2)) > tol) then
      error stop 'DRESP-10 Frechet off-diagonal action failed'
   end if
   if (abs(delta_rho(2,3) - expected*delta_h(2,3)) > tol*max(1.0_rp, abs(expected))) then
      error stop 'DRESP-10 Frechet near-degenerate action failed'
   end if
   if (maxval(abs(delta_rho - transpose(conjg(delta_rho)))) > tol) error stop 'DRESP-10 Frechet Hermiticity failed'
   write (*, '(a)') 'UnitLrStaticFrechetProductResponse: PASS'
end program test_lr_static_frechet_product_response
