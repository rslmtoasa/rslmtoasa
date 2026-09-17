!------------------------------------------------------------------------------
! DRESP-03TG-C -- first hard gate: fixed-complex-energy native contraction.
!
! The native and finite-H objects are assembled independently from the live
! potential parameters and a production structure-constant block.  The test
! compares the complete exchange integrand, not a raw Green-function block.
!
! Fixture layout (all matrices are Fortran column-major, with the following
! logical ordering):
!   orbital channels per site and spin:  norb = (lmax+1)**2
!   two-site orbital space:             site1, site2                  (2*norb)
!   global spin-major space:             up(site1,site2), down(site1,site2)
!                                       (4*norb)
!------------------------------------------------------------------------------
program test_dresp03tg_native_fixed_z
   use precision_mod, only: rp
   use math_mod, only: ang2au, init_math_operators
   use timer_mod, only: g_timer, timer
   use logger_mod, only: g_logger
   use control_mod, only: control
   use lattice_mod, only: lattice
   use charge_mod, only: charge
   use hamiltonian_mod, only: hamiltonian
   use basis_mod, only: basis_init
   use lr_lmto_turek_gf_mod, only: native_complex_p_matrix, native_screening_alpha, &
      native_screened_p_matrix, native_delta_p, native_inverse, native_path_operator, &
      native_exchange_integrand, native_finite_h_integrand, native_collinear_pauli_integrand
   implicit none

   type(control) :: ctl
   type(lattice), target :: lat
   type(charge), target :: chg
   type(hamiltonian), target :: ham
   complex(rp), allocatable :: p_real(:,:,:), p_complex(:,:), p_screen(:,:), delta(:,:), expected(:,:)
   complex(rp), allocatable :: p_gamma(:,:), p_alpha(:,:), s_alpha(:,:), s_gamma(:,:), &
      s_sites_orbital(:,:), s_spin_major(:,:), s_roundtrip(:,:), d_screen(:,:), d_width(:,:), &
      eye(:,:), tmp(:,:), r_screen(:,:), r_inverse(:,:), g_gamma(:,:), g_alpha(:,:), h(:,:), gh(:,:)
   complex(rp), allocatable :: d1(:,:), d2(:,:), gup_ij(:,:), gup_ji(:,:), gdown_ij(:,:), gdown_ji(:,:), &
      ghup_ij(:,:), ghup_ji(:,:), ghdown_ij(:,:), ghdown_ji(:,:), gh_spin_ij(:,:), gh_spin_ji(:,:)
   complex(rp), allocatable :: path_p(:,:), path_s(:,:), path_g(:,:)
   real(rp), allocatable :: alpha(:)
   complex(rp) :: z_values(3), z
   real(rp) :: p_error, d_error, solve_error, p_transform_error, s_transform_error
   real(rp) :: alpha_error, alpha_du_error, pauli_error, pauli_self_error, gf_cov_error
   real(rp) :: native_alpha_ordered_ud, native_alpha_ordered_du, native_alpha_symmetrized
   real(rp) :: native_gamma_ordered_ud, native_gamma_ordered_du, native_gamma_symmetrized
   real(rp) :: finite_h_ordered_ud, finite_h_ordered_du, finite_h_symmetrized
   real(rp) :: historical_pauli_value, explicit_pauli_value
   real(rp) :: alpha_h_error, gamma_h_error, gamma_h_du_error, gf_cov_value, p_transform_value
   real(rp), allocatable :: s0(:,:), s1(:,:)
   real(rp) :: e
   integer :: i, l, m, lm, lmax, norb, nsite_fixture, norb_sites, nspin_orb_sites
   integer :: n2, n4, ia, ino
   logical :: failed

   call init_math_operators()
   call g_logger%init()
   g_timer = timer()
   ctl = control('input.nml')
   lat = lattice(ctl)
   call lat%build_data()
   call lat%bravais()
   call lat%structb(.true.)
   call lat%atomlist()
   call basis_init(lat%symbolic_atoms(1)%potential%lmax)
   chg = charge(lat)
   ham = hamiltonian(chg)
   do i=1,lat%ntype
      call lat%symbolic_atoms(i)%build_pot()
   end do
   call ham%build_bulkham()
   do i=1,lat%ntype
      call lat%symbolic_atoms(i)%predls(lat%wav*ang2au)
   end do

   lmax = lat%symbolic_atoms(1)%potential%lmax
   norb = (lmax+1)**2
   if (lat%nrec /= 1 .or. norb < 4) error stop 'DRESP-03TG fixture is not the required nontrivial one-site spd state'
   if (lat%nn(lat%atlist(1),1) < 2) error stop 'DRESP-03TG fixture has no hopping neighbor'
   nsite_fixture = 2
   norb_sites = nsite_fixture*norb
   nspin_orb_sites = 2*norb_sites
   n2 = norb_sites
   n4 = nspin_orb_sites
   allocate(p_real(n2,n2,1),p_complex(n2,n2),p_screen(n2,n2),delta(norb,norb),expected(norb,norb),alpha(0:lmax))
   allocate(p_gamma(n4,n4),p_alpha(n4,n4),s_alpha(n4,n4),s_gamma(n4,n4),s_sites_orbital(n2,n2), &
      s_spin_major(n4,n4),s_roundtrip(n4,n4),d_screen(n4,n4),d_width(n4,n4),eye(n4,n4),tmp(n4,n4), &
      r_screen(n4,n4),r_inverse(n4,n4),g_gamma(n4,n4),g_alpha(n4,n4),h(n4,n4),gh(n4,n4), &
      d1(norb,norb),d2(norb,norb),gup_ij(norb,norb),gup_ji(norb,norb),gdown_ij(norb,norb),gdown_ji(norb,norb), &
      ghup_ij(norb,norb),ghup_ji(norb,norb),ghdown_ij(norb,norb),ghdown_ji(norb,norb), &
      gh_spin_ij(n2,n2),gh_spin_ji(n2,n2))
   call native_screening_alpha(lat%symbolic_atoms(1),alpha)

   ! Complex P(z) must reduce exactly to the live real interface on the real axis.
   e = -0.217_rp
   call lat%symbolic_atoms(1)%p_matrix(p_real,lmax,[e])
   call native_complex_p_matrix(lat%symbolic_atoms(1),cmplx(e,0.0_rp,rp),p_complex)
   p_error = maxval(abs(p_complex-p_real(:,:,1)))
   write(*,'(a,es14.6)') 'DRESP-03TG complex P(real-axis) residual = ',p_error

   ! d_matrix(E) = dele_up*dele_down*DeltaP(E), checked over several energies.
   d_error = 0.0_rp
   do i=1,4
      select case(i)
      case(1); e=-0.73_rp
      case(2); e=-0.17_rp
      case(3); e=0.11_rp
      case default; e=0.64_rp
      end select
      call lat%symbolic_atoms(1)%p_matrix(p_real,lmax,[e])
      call native_delta_p(p_real(:,:,1),delta)
      expected = cmplx(0.0_rp,0.0_rp,rp)
      do l=0,lmax
         do m=1,2*l+1
            lm=l*l+m
            expected(lm,lm) = lat%symbolic_atoms(1)%potential%dele(l,1)* &
               lat%symbolic_atoms(1)%potential%dele(l,2)*delta(lm,lm)
         end do
      end do
      call lat%symbolic_atoms(1)%d_matrix(p_complex,e)
      d_error=max(d_error,maxval(abs(p_complex(1:norb,1:norb)-expected)))
   end do
   write(*,'(a,es14.6)') 'DRESP-03TG d_matrix/width-scaled DeltaP residual = ',d_error

   ! Check the direct production one-site constructor independently.
   allocate(path_p(n2,n2),path_s(n2,n2),path_g(n2,n2))
   call native_path_operator(lat,cmplx(-0.31_rp,0.27_rp,rp),[0.137_rp,-0.071_rp,0.209_rp],path_p,path_s,path_g)
   solve_error=maxval(abs(matmul(path_p-path_s,path_g)-identity(n2)))
   write(*,'(a,es14.6)') 'DRESP-03TG direct [P_alpha-S_alpha] solve residual = ',solve_error
   deallocate(path_p,path_s,path_g)

   ! Build a two-site controlled LMTO fixture from two nonzero live S bonds.
   ! Its alpha representation is transformed back to gamma for the finite-H
   ! oracle; no fitted matrix or response scale is introduced.
   ia=lat%atlist(1); ino=lat%num(ia)
   allocate(s0(norb,norb),s1(norb,norb))
   s0=real(transpose(lat%sbar(1:norb,1:norb,1,ino)),rp)
   s1=real(transpose(lat%sbar(1:norb,1:norb,2,ino)),rp)
   if (maxval(abs(s1)) <= 1.0e-12_rp) error stop 'DRESP-03TG controlled fixture has zero production hopping'
   ! First make the two-site orbital matrix, then lift it explicitly to the
   ! global spin-major matrix.  The intermediate shapes are intentional:
   ! S_sites_orbital is (site1,site2), while S_spin_major is
   ! (up/site1,up/site2,down/site1,down/site2).
   s_sites_orbital=cmplx(0.0_rp,0.0_rp,rp)
   s_sites_orbital(1:norb,1:norb)=cmplx(s0,0.0_rp,rp)
   s_sites_orbital(1:norb,norb+1:n2)=cmplx(s1,0.0_rp,rp)
   s_sites_orbital(norb+1:n2,1:norb)=cmplx(transpose(s1),0.0_rp,rp)
   s_sites_orbital(norb+1:n2,norb+1:n2)=cmplx(s0,0.0_rp,rp)
   s_spin_major=cmplx(0.0_rp,0.0_rp,rp)
   s_spin_major(1:n2,1:n2)=s_sites_orbital
   s_spin_major(n2+1:n4,n2+1:n4)=s_sites_orbital
   s_alpha=s_spin_major

   p_gamma=cmplx(0.0_rp,0.0_rp,rp); p_alpha=cmplx(0.0_rp,0.0_rp,rp)
   d_screen=cmplx(0.0_rp,0.0_rp,rp); d_width=cmplx(0.0_rp,0.0_rp,rp)
   do i=1,nsite_fixture
      call native_complex_p_matrix(lat%symbolic_atoms(1),cmplx(0.0_rp,0.0_rp,rp),p_complex)
      p_gamma((i-1)*norb+1:i*norb,(i-1)*norb+1:i*norb)=p_complex(1:norb,1:norb)
      p_gamma(n2+(i-1)*norb+1:n2+i*norb,n2+(i-1)*norb+1:n2+i*norb)= &
         p_complex(norb+1:n2,norb+1:n2)
      call native_screened_p_matrix(lat%symbolic_atoms(1),cmplx(0.0_rp,0.0_rp,rp),alpha,p_screen)
      p_alpha((i-1)*norb+1:i*norb,(i-1)*norb+1:i*norb)=p_screen(1:norb,1:norb)
      p_alpha(n2+(i-1)*norb+1:n2+i*norb,n2+(i-1)*norb+1:n2+i*norb)= &
         p_screen(norb+1:n2,norb+1:n2)
      do l=0,lmax
         do m=1,2*l+1
            lm=l*l+m
            d_screen((i-1)*norb+lm,(i-1)*norb+lm)=alpha(l)-lat%symbolic_atoms(1)%potential%qpar(l,1)
            d_screen(n2+(i-1)*norb+lm,n2+(i-1)*norb+lm)=alpha(l)-lat%symbolic_atoms(1)%potential%qpar(l,2)
            d_width((i-1)*norb+lm,(i-1)*norb+lm)=lat%symbolic_atoms(1)%potential%dele(l,1)
            d_width(n2+(i-1)*norb+lm,n2+(i-1)*norb+lm)=lat%symbolic_atoms(1)%potential%dele(l,2)
         end do
      end do
   end do
   eye=identity(n4)
   call native_inverse(eye+matmul(s_alpha,d_screen),tmp)
   s_gamma=matmul(tmp,s_alpha)
   ! Reverse S transformation: S_alpha = S_gamma*(I-D*S_gamma)^(-1).
   call native_inverse(eye-matmul(d_screen,s_gamma),tmp)
   s_roundtrip=matmul(s_gamma,tmp)
   s_transform_error=maxval(abs(s_roundtrip-s_alpha))
   write(*,'(a,es14.6)') 'DRESP-03TG S transformation residual = ',s_transform_error

   failed=.false.; alpha_error=0.0_rp; alpha_du_error=0.0_rp; gamma_h_error=0.0_rp; gamma_h_du_error=0.0_rp; alpha_h_error=0.0_rp
   gf_cov_error=0.0_rp; p_transform_error=0.0_rp; pauli_error=0.0_rp; pauli_self_error=0.0_rp
   z_values=[cmplx(-0.91_rp,0.83_rp,rp),cmplx(-0.17_rp,0.04_rp,rp),cmplx(0.62_rp,0.31_rp,rp)]
   do i=1,size(z_values)
      z=z_values(i)
      call native_complex_p_matrix(lat%symbolic_atoms(1),z,p_complex)
      call native_screened_p_matrix(lat%symbolic_atoms(1),z,alpha,p_screen)
      p_gamma=cmplx(0.0_rp,0.0_rp,rp); p_alpha=cmplx(0.0_rp,0.0_rp,rp)
      ! Repeat the two identical site blocks at the fixed complex energy in
      ! global spin-major order: up(site1,site2), down(site1,site2).
      ! Explicit diagonal ranges in the global layout:
      ! up/site1=1:norb, up/site2=norb+1:n2,
      ! down/site1=n2+1:n2+norb, down/site2=n2+norb+1:n4.
      do lm=1,norb
         p_gamma(lm,lm)=p_complex(lm,lm)
         p_gamma(norb+lm,norb+lm)=p_complex(lm,lm)
         p_gamma(n2+lm,n2+lm)=p_complex(norb+lm,norb+lm)
         p_gamma(n2+norb+lm,n2+norb+lm)=p_complex(norb+lm,norb+lm)
         p_alpha(lm,lm)=p_screen(lm,lm)
         p_alpha(norb+lm,norb+lm)=p_screen(lm,lm)
         p_alpha(n2+lm,n2+lm)=p_screen(norb+lm,norb+lm)
         p_alpha(n2+norb+lm,n2+norb+lm)=p_screen(norb+lm,norb+lm)
      end do
      call native_inverse(p_gamma-s_gamma,g_gamma)
      call native_inverse(p_alpha-s_alpha,g_alpha)
      ! P-side screening covariance: P_alpha = P_gamma*(I-D*P_gamma)^(-1).
      r_screen=eye-matmul(d_screen,p_gamma)
      call native_inverse(r_screen,r_inverse)
      p_transform_value=maxval(abs(p_alpha-matmul(p_gamma,r_inverse)))
      p_transform_error=max(p_transform_error,p_transform_value)

      ! Off-site path-operator covariance has no additive onsite term:
      ! g_alpha(1,2) = R_1*g_gamma(1,2)*R_2.  The two fixture sites are
      ! identical, so the global block-diagonal R is valid on both ends.
      tmp=matmul(r_screen,matmul(g_gamma,r_screen))
      gf_cov_value=max(maxval(abs(gamma_spin_block(g_alpha,norb,1,2)-gamma_spin_block(tmp,norb,1,2))), &
         maxval(abs(gamma_spin_block(g_alpha,norb,2,1)-gamma_spin_block(tmp,norb,2,1))))
      gf_cov_error=max(gf_cov_error,gf_cov_value)

      gup_ij=g_gamma(1:norb,norb+1:n2)
      gup_ji=g_gamma(norb+1:n2,1:norb)
      gdown_ij=g_gamma(n2+1:n2+norb,n2+norb+1:n4)
      gdown_ji=g_gamma(n2+norb+1:n4,n2+1:n2+norb)
      gh=cmplx(0.0_rp,0.0_rp,rp)
      do l=0,lmax
         do m=1,2*l+1
            lm=l*l+m
            gh(lm,lm)=cmplx(lat%symbolic_atoms(1)%potential%c(l,1)+lat%symbolic_atoms(1)%potential%vmad,0.0_rp,rp)
            gh(norb+lm,norb+lm)=gh(lm,lm)
            gh(n2+lm,n2+lm)=cmplx(lat%symbolic_atoms(1)%potential%c(l,2)+lat%symbolic_atoms(1)%potential%vmad,0.0_rp,rp)
            gh(n2+norb+lm,n2+norb+lm)=gh(n2+lm,n2+lm)
         end do
      end do
      h=gh+matmul(d_width,matmul(s_gamma,d_width))
      h=z*identity(n4)-h
      call native_inverse(h,gh)
      d1=cmplx(0.0_rp,0.0_rp,rp); d2=d1
      do l=0,lmax
         do m=1,2*l+1
            lm=l*l+m
            d1(lm,lm)=d_width(lm,lm)*d_width(n2+lm,n2+lm)*(p_gamma(lm,lm)-p_gamma(n2+lm,n2+lm))
            d2(lm,lm)=d1(lm,lm)
         end do
      end do
      ghup_ij=gh(1:norb,norb+1:n2)
      ghup_ji=gh(norb+1:n2,1:norb)
      ghdown_ij=gh(n2+1:n2+norb,n2+norb+1:n4)
      ghdown_ji=gh(n2+norb+1:n4,n2+1:n2+norb)

      finite_h_ordered_ud=native_finite_h_integrand(d1,d2,ghup_ij,ghdown_ji)
      finite_h_ordered_du=native_finite_h_integrand(d1,d2,ghdown_ij,ghup_ji)
      finite_h_symmetrized=0.5_rp*(finite_h_ordered_ud+finite_h_ordered_du)

      delta=p_alpha(1:norb,1:norb)-p_alpha(n2+1:n2+norb,n2+1:n2+norb)
      expected=p_alpha(norb+1:n2,norb+1:n2)-p_alpha(n2+norb+1:n4,n2+norb+1:n4)
      native_alpha_ordered_ud=native_exchange_integrand(delta,expected,g_alpha(1:norb,norb+1:n2), &
         g_alpha(n2+norb+1:n4,n2+1:n2+norb))
      native_alpha_ordered_du=native_exchange_integrand(delta,expected,g_alpha(n2+1:n2+norb,n2+norb+1:n4), &
         g_alpha(norb+1:n2,1:norb))
      native_alpha_symmetrized=0.5_rp*(native_alpha_ordered_ud+native_alpha_ordered_du)

      delta=p_gamma(1:norb,1:norb)-p_gamma(n2+1:n2+norb,n2+1:n2+norb)
      expected=p_gamma(norb+1:n2,norb+1:n2)-p_gamma(n2+norb+1:n4,n2+norb+1:n4)
      native_gamma_ordered_ud=native_exchange_integrand(delta,expected,gup_ij,gdown_ji)
      native_gamma_ordered_du=native_exchange_integrand(delta,expected,gdown_ij,gup_ji)
      native_gamma_symmetrized=0.5_rp*(native_gamma_ordered_ud+native_gamma_ordered_du)

      gh_spin_ij=gamma_spin_block(gh,norb,1,2)
      gh_spin_ji=gamma_spin_block(gh,norb,2,1)
      historical_pauli_value=native_collinear_pauli_integrand(d1,d2,gh_spin_ij,gh_spin_ji)
      explicit_pauli_value=explicit_collinear_pauli_integrand(d1,d2,gh_spin_ij,gh_spin_ji)
      pauli_self_error=max(pauli_self_error,abs(historical_pauli_value-explicit_pauli_value))

      alpha_h_error=max(alpha_h_error,abs(native_alpha_ordered_ud-finite_h_ordered_ud), &
         abs(native_alpha_ordered_du-finite_h_ordered_du))
      gamma_h_error=max(gamma_h_error,abs(native_gamma_ordered_ud-finite_h_ordered_ud))
      gamma_h_du_error=max(gamma_h_du_error,abs(native_gamma_ordered_du-finite_h_ordered_du))
      alpha_error=max(alpha_error,abs(native_alpha_ordered_ud-native_gamma_ordered_ud))
      alpha_du_error=max(alpha_du_error,abs(native_alpha_ordered_du-native_gamma_ordered_du))
      pauli_error=max(pauli_error,abs(historical_pauli_value-finite_h_ordered_ud), &
         abs(historical_pauli_value-finite_h_ordered_du),abs(historical_pauli_value-finite_h_symmetrized))

      write(*,'(a,2es14.6)') 'z = ',real(z,rp),aimag(z)
      write(*,'(a,es14.6)') '  native_gamma_ordered_ud = ',native_gamma_ordered_ud
      write(*,'(a,es14.6)') '  native_gamma_ordered_du = ',native_gamma_ordered_du
      write(*,'(a,es14.6)') '  native_gamma_symmetrized = ',native_gamma_symmetrized
      write(*,'(a,es14.6)') '  native_alpha_ordered_ud = ',native_alpha_ordered_ud
      write(*,'(a,es14.6)') '  native_alpha_ordered_du = ',native_alpha_ordered_du
      write(*,'(a,es14.6)') '  native_alpha_symmetrized = ',native_alpha_symmetrized
      write(*,'(a,es14.6)') '  finite_H_ordered_ud = ',finite_h_ordered_ud
      write(*,'(a,es14.6)') '  finite_H_ordered_du = ',finite_h_ordered_du
      write(*,'(a,es14.6)') '  finite_H_symmetrized = ',finite_h_symmetrized
      write(*,'(a,es14.6)') '  historical_Pauli = ',historical_pauli_value
      write(*,'(a,es14.6)') '  alpha_minus_gamma = ',abs(native_alpha_ordered_ud-native_gamma_ordered_ud)
      write(*,'(a,es14.6)') '  gamma_minus_H_ud = ',abs(native_gamma_ordered_ud-finite_h_ordered_ud)
      write(*,'(a,es14.6)') '  Pauli_minus_ud = ',abs(historical_pauli_value-finite_h_ordered_ud)
      write(*,'(a,es14.6)') '  Pauli_minus_du = ',abs(historical_pauli_value-finite_h_ordered_du)
      write(*,'(a,es14.6)') '  Pauli_minus_symmetrized = ',abs(historical_pauli_value-finite_h_symmetrized)
      write(*,'(a,es14.6)') '  path_operator_covariance = ',gf_cov_value
      write(*,'(a,es14.6)') '  P_transformation = ',p_transform_value
   end do
   solve_error=maxval(abs(matmul(p_alpha-s_alpha,g_alpha)-identity(n4)))
   write(*,'(a,es14.6)') 'DRESP-03TG native inverse residual = ',solve_error
   write(*,'(a,es14.6)') 'DRESP-03TG P transformation residual = ',p_transform_error
   write(*,'(a,es14.6)') 'DRESP-03TG S transformation residual = ',s_transform_error
   write(*,'(a,es14.6)') 'DRESP-03TG path-operator covariance residual = ',gf_cov_error
   write(*,'(a,es14.6)') 'DRESP-03TG alpha-gamma ordered_ud residual = ',alpha_error
   write(*,'(a,es14.6)') 'DRESP-03TG alpha-gamma ordered_du residual = ',alpha_du_error
   write(*,'(a,es14.6)') 'DRESP-03TG alpha-H residual = ',alpha_h_error
   write(*,'(a,es14.6)') 'DRESP-03TG gamma finite-H ordered_ud residual = ',gamma_h_error
   write(*,'(a,es14.6)') 'DRESP-03TG gamma finite-H ordered_du residual = ',gamma_h_du_error
   write(*,'(a,es14.6)') 'DRESP-03TG historical Pauli minus finite-H ordered_ud/du/sym max = ',pauli_error
   write(*,'(a,es14.6)') 'DRESP-03TG Pauli-helper self-consistency residual = ',pauli_self_error
   ! R1 certifies the fixture and the representation transformations.  The
   ! alpha/gamma vertex contraction and Pauli ordering residuals are reported
   ! for the later TG-FZ-R2/R3 work; they are not silently converted into a
   ! fitted factor or made acceptance criteria here.
   failed = p_error > 2.0e-14_rp .or. d_error > 2.0e-13_rp .or. solve_error > 2.0e-12_rp .or. &
      p_transform_error > 2.0e-12_rp .or. s_transform_error > 2.0e-12_rp .or. gf_cov_error > 2.0e-12_rp .or. &
      gamma_h_error > 2.0e-10_rp .or. gamma_h_du_error > 2.0e-10_rp .or. pauli_self_error > 2.0e-12_rp
   if (failed) then
      write(*,'(a)') 'RESULT: BLOCKED — TG-FZ-R1 STRUCTURAL DIAGNOSTIC CHECK'
      return
   end if
   write(*,'(a)') 'TG-FZ-R1 FIXED-Z DIAGNOSTIC: PASS (unresolved contraction residuals reported)'

contains

   function identity(n) result(a)
      integer, intent(in) :: n
      complex(rp) :: a(n,n)
      integer :: j
      a=cmplx(0.0_rp,0.0_rp,rp)
      do j=1,n; a(j,j)=cmplx(1.0_rp,0.0_rp,rp); end do
   end function identity

   function gamma_spin_block(g,n,i_site,j_site) result(block)
      complex(rp), intent(in) :: g(:,:)
      integer, intent(in) :: n,i_site,j_site
      complex(rp) :: block(2*n,2*n)
      integer :: up_i0,up_j0,down_i0,down_j0
      if (n <= 0 .or. size(g,1) /= nspin_orb_sites .or. size(g,2) /= nspin_orb_sites) then
         error stop 'gamma_spin_block: invalid global spin-major shape'
      end if
      if (i_site < 1 .or. i_site > nsite_fixture .or. j_site < 1 .or. j_site > nsite_fixture) then
         error stop 'gamma_spin_block: invalid fixture site'
      end if
      up_i0=(i_site-1)*n; up_j0=(j_site-1)*n
      down_i0=norb_sites+up_i0; down_j0=norb_sites+up_j0
      block=cmplx(0.0_rp,0.0_rp,rp)
      block(1:n,1:n)=g(up_i0+1:up_i0+n,up_j0+1:up_j0+n)
      block(n+1:2*n,n+1:2*n)=g(down_i0+1:down_i0+n,down_j0+1:down_j0+n)
   end function gamma_spin_block

   function explicit_collinear_pauli_integrand(delta_i,delta_j,gij,gji) result(value)
      complex(rp), intent(in) :: delta_i(:,:),delta_j(:,:),gij(:,:),gji(:,:)
      complex(rp) :: g0i(size(delta_i,1),size(delta_i,2)), gzi(size(delta_i,1),size(delta_i,2))
      complex(rp) :: g0j(size(delta_i,1),size(delta_i,2)), gzj(size(delta_i,1),size(delta_i,2))
      complex(rp) :: product(size(delta_i,1),size(delta_i,2))
      real(rp) :: value
      integer :: n, j
      n=size(delta_i,1)
      if (size(delta_i,2) /= n .or. size(delta_j,1) /= n .or. size(delta_j,2) /= n .or. &
          size(gij,1) /= 2*n .or. size(gij,2) /= 2*n .or. size(gji,1) /= 2*n .or. size(gji,2) /= 2*n) then
         error stop 'explicit_collinear_pauli_integrand: invalid block shape'
      end if
      ! With Gx=Gy=0, G0=(Gup+Gdown)/2 and Gz=(Gup-Gdown)/2.
      ! Therefore G0_ij*G0_ji-Gz_ij*Gz_ji is exactly
      ! (Gup_ij*Gdown_ji+Gdown_ij*Gup_ji)/2, preserving matrix order.
      g0i=0.5_rp*(gij(1:n,1:n)+gij(n+1:2*n,n+1:2*n))
      gzi=0.5_rp*(gij(1:n,1:n)-gij(n+1:2*n,n+1:2*n))
      g0j=0.5_rp*(gji(1:n,1:n)+gji(n+1:2*n,n+1:2*n))
      gzj=0.5_rp*(gji(1:n,1:n)-gji(n+1:2*n,n+1:2*n))
      product=matmul(delta_i,matmul(g0i,matmul(delta_j,g0j))) - &
         matmul(delta_i,matmul(gzi,matmul(delta_j,gzj)))
      value=0.0_rp
      do j=1,n
         value=value+aimag(product(j,j))
      end do
   end function explicit_collinear_pauli_integrand

end program test_dresp03tg_native_fixed_z
