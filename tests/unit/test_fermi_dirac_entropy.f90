program test_fermi_dirac_entropy
   use precision_mod, only: rp
   use reciprocal_mod, only: fermi_dirac_occupation, fermi_dirac_entropy, fermi_dirac_entropy_term
   implicit none

   real(rp), parameter :: kB = 6.3336814e-6_rp
   real(rp), parameter :: temperature = 300.0_rp
   real(rp), parameter :: mu = 0.0_rp
   real(rp), parameter :: kT = 0.00190010442_rp
   real(rp), parameter :: expected_s = 2.106254006907724_rp
   real(rp), parameter :: expected_ts = 0.004002102548168077_rp
   real(rp), parameter :: bound_s = 1.198e-13_rp
   real(rp), parameter :: bound_ts = 2.277e-16_rp
   real(rp), parameter :: eigenvalues(4, 2) = reshape([ &
      -0.003_rp, 0.001_rp, -0.0015_rp, 0.0025_rp, &
      -0.002_rp, 0.003_rp, -0.004_rp, 0.001_rp], [4, 2])
   real(rp), parameter :: weights(2) = [1.0_rp/3.0_rp, 2.0_rp/3.0_rp]
   real(rp) :: entropy, ts, entropy_reference, occupation
   integer :: ik, ib

   entropy = fermi_dirac_entropy(eigenvalues, weights, mu, kT)
   ts = kB*temperature*entropy
   entropy_reference = 0.0_rp
   do ik = 1, size(eigenvalues, 2)
      do ib = 1, size(eigenvalues, 1)
         occupation = 1.0_rp/(1.0_rp + exp((eigenvalues(ib, ik) - mu)/kT))
         entropy_reference = entropy_reference - weights(ik)*( &
            occupation*log(occupation) + (1.0_rp - occupation)*log(1.0_rp - occupation))
      end do
   end do

   if (abs(kB*temperature - kT) > 2.0_rp*epsilon(kT)) error stop 'temperature/kT fixture mismatch'
   if (abs(entropy - expected_s) > bound_s) error stop 'production entropy differs from locked analytic value'
   if (abs(entropy - entropy_reference) > bound_s) error stop 'production entropy differs from direct FD sum'
   if (abs(ts - expected_ts) > bound_ts) error stop 'production TS differs from locked analytic value'
   if (fermi_dirac_entropy_term(0.0_rp) /= 0.0_rp) error stop 'empty-state entropy is not zero'
   if (fermi_dirac_entropy_term(1.0_rp) /= 0.0_rp) error stop 'full-state entropy is not zero'
   if (fermi_dirac_entropy_term(fermi_dirac_occupation(100.0_rp, mu, kT)) /= 0.0_rp) &
      error stop 'cutoff empty-state entropy is not zero'
   if (fermi_dirac_entropy_term(fermi_dirac_occupation(-100.0_rp, mu, kT)) /= 0.0_rp) &
      error stop 'cutoff full-state entropy is not zero'

   write (*, '(a,es24.16)') 'S=', entropy
   write (*, '(a,es24.16)') 'TS_Ry=', ts
end program test_fermi_dirac_entropy
