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
   use lr_kl_hessian_mod, only: lmto_live_hamiltonian_fixture, lmto_fixture_init, &
      assemble_lmto_hamiltonian, assemble_lmto_rotation_terms, assemble_lmto_mixed_derivative, &
      force_theorem_integrand, force_theorem_pi
   use lr_lmto_turek_gf_mod, only: native_complex_p_matrix, native_screening_alpha, &
      native_screened_p_matrix, native_delta_p, native_inverse, native_path_operator, &
      native_exchange_integrand, native_finite_h_integrand, native_collinear_pauli_integrand
   implicit none

   type(control) :: ctl
   type(lattice), target :: lat
   type(charge), target :: chg
   type(hamiltonian), target :: ham
   type(lmto_live_hamiltonian_fixture) :: r4_fixture
   complex(rp), allocatable :: p_real(:,:,:), p_complex(:,:), p_screen(:,:), delta(:,:), expected(:,:)
   complex(rp), allocatable :: p_gamma(:,:), p_alpha(:,:), s_alpha(:,:), s_gamma(:,:), &
      s_sites_orbital(:,:), s_spin_major(:,:), s_roundtrip(:,:), d_screen(:,:), d_width(:,:), &
      eye(:,:), tmp(:,:), r_screen(:,:), r_inverse(:,:), g_gamma(:,:), g_alpha(:,:), h(:,:), gh(:,:)
   complex(rp), allocatable :: d1(:,:), d2(:,:), gup_ij(:,:), gup_ji(:,:), gdown_ij(:,:), gdown_ji(:,:), &
      ghup_ij(:,:), ghup_ji(:,:), ghdown_ij(:,:), ghdown_ji(:,:), gh_spin_ij(:,:), gh_spin_ji(:,:)
   complex(rp), allocatable :: delta_alpha_i(:,:), delta_alpha_j(:,:), delta_gamma_i(:,:), delta_gamma_j(:,:), &
      correction_i(:,:), correction_j(:,:), candidate_i(:,:), candidate_j(:,:), &
      r_up_i(:,:), r_up_j(:,:), r_down_i(:,:), r_down_j(:,:), &
      vertex_i_ud(:,:), vertex_j_ud(:,:), vertex_i_du(:,:), vertex_j_du(:,:), &
      expanded_i_ud(:,:), expanded_j_ud(:,:), expanded_i_du(:,:), expanded_j_du(:,:), &
      d_tilde_i(:,:), d_tilde_j(:,:), delta_d_i(:,:), delta_d_j(:,:)
   complex(rp), allocatable :: path_p(:,:), path_s(:,:), path_g(:,:)
   real(rp), allocatable :: alpha(:)
   complex(rp) :: z_values(3), z
   real(rp) :: p_error, d_error, solve_error, p_transform_error, s_transform_error
   real(rp) :: alpha_error, alpha_du_error, pauli_error, pauli_self_error, gf_cov_error
   real(rp) :: vertex_transform_error, expanded_vertex_error, alpha_transformed_error, raw_gamma_transformed_error
   real(rp) :: delta_d_error, common_gamma_error, common_gamma_vertex_error
   real(rp) :: native_alpha_ordered_ud, native_alpha_ordered_du, native_alpha_symmetrized
   real(rp) :: native_gamma_ordered_ud, native_gamma_ordered_du, native_gamma_symmetrized
   real(rp) :: transformed_gamma_ordered_ud, transformed_gamma_ordered_du
   real(rp) :: finite_h_ordered_ud, finite_h_ordered_du, finite_h_symmetrized
   real(rp) :: historical_pauli_value, explicit_pauli_value
   real(rp) :: alpha_h_error, gamma_h_error, gamma_h_du_error, gf_cov_value, p_transform_value
   real(rp) :: raw_gamma_vertex_magnitude, screening_correction_magnitude, transformed_vertex_magnitude
   real(rp) :: common_gamma_alpha_ud, common_gamma_alpha_du, common_gamma_raw_ud, common_gamma_raw_du
   real(rp) :: vertex_transform_value, expanded_vertex_value, alpha_transformed_value, raw_gamma_transformed_value
   real(rp) :: delta_d_value, common_gamma_vertex_value, common_gamma_correction_magnitude
   real(rp) :: common_gamma_correction_error
   real(rp), allocatable :: s0(:,:), s1(:,:)
   real(rp) :: e
   integer :: i, l, m, lm, lmax, norb, nsite_fixture, norb_sites, nspin_orb_sites
   integer :: n2, n4, ia, ino, isite
   logical :: failed
   real(rp) :: r4_max_fd_error, r4_max_contact_norm, r4_max_vertex_raw_error, r4_max_vertex_transformed_error

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
      gh_spin_ij(n2,n2),gh_spin_ji(n2,n2),delta_alpha_i(norb,norb),delta_alpha_j(norb,norb), &
      delta_gamma_i(norb,norb),delta_gamma_j(norb,norb),correction_i(norb,norb),correction_j(norb,norb), &
      candidate_i(norb,norb),candidate_j(norb,norb),r_up_i(norb,norb),r_up_j(norb,norb), &
      r_down_i(norb,norb),r_down_j(norb,norb),vertex_i_ud(norb,norb),vertex_j_ud(norb,norb), &
      vertex_i_du(norb,norb),vertex_j_du(norb,norb),expanded_i_ud(norb,norb),expanded_j_ud(norb,norb), &
      expanded_i_du(norb,norb),expanded_j_du(norb,norb),d_tilde_i(norb,norb),d_tilde_j(norb,norb), &
      delta_d_i(norb,norb),delta_d_j(norb,norb))
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
   vertex_transform_error=0.0_rp; expanded_vertex_error=0.0_rp; alpha_transformed_error=0.0_rp
   raw_gamma_transformed_error=0.0_rp; delta_d_error=0.0_rp; common_gamma_error=0.0_rp
   common_gamma_vertex_error=0.0_rp; common_gamma_correction_error=0.0_rp
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
      ! Restore the physical spin-dependent D matrix after the isolated
      ! common-gamma control from the preceding z, then rebuild S_gamma.
      d_screen=cmplx(0.0_rp,0.0_rp,rp)
      do isite=1,nsite_fixture
         do l=0,lmax
            do m=1,2*l+1
               lm=l*l+m
               d_screen((isite-1)*norb+lm,(isite-1)*norb+lm)= &
                  alpha(l)-lat%symbolic_atoms(1)%potential%qpar(l,1)
               d_screen(n2+(isite-1)*norb+lm,n2+(isite-1)*norb+lm)= &
                  alpha(l)-lat%symbolic_atoms(1)%potential%qpar(l,2)
            end do
         end do
      end do
      call native_inverse(eye+matmul(s_alpha,d_screen),tmp)
      s_gamma=matmul(tmp,s_alpha)
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

      ! Construct the local raw vertices and endpoint transformations.  For
      ! ordered ud, cyclicity of the trace gives
      !   DeltaP_i^a g^a_up,ij DeltaP_j^a g^a_down,ji
      ! = (R_i,down DeltaP_i^a R_i,up) g^g_up,ij
      !   (R_j,up DeltaP_j^a R_j,down) g^g_down,ji.
      delta_alpha_i=p_alpha(1:norb,1:norb)-p_alpha(n2+1:n2+norb,n2+1:n2+norb)
      delta_alpha_j=p_alpha(norb+1:n2,norb+1:n2)-p_alpha(n2+norb+1:n4,n2+norb+1:n4)
      delta_gamma_i=p_gamma(1:norb,1:norb)-p_gamma(n2+1:n2+norb,n2+1:n2+norb)
      delta_gamma_j=p_gamma(norb+1:n2,norb+1:n2)-p_gamma(n2+norb+1:n4,n2+norb+1:n4)
      r_up_i=r_screen(1:norb,1:norb)
      r_up_j=r_screen(norb+1:n2,norb+1:n2)
      r_down_i=r_screen(n2+1:n2+norb,n2+1:n2+norb)
      r_down_j=r_screen(n2+norb+1:n4,n2+norb+1:n4)
      vertex_i_ud=matmul(r_down_i,matmul(delta_alpha_i,r_up_i))
      vertex_j_ud=matmul(r_up_j,matmul(delta_alpha_j,r_down_j))
      vertex_i_du=matmul(r_up_i,matmul(delta_alpha_i,r_down_i))
      vertex_j_du=matmul(r_down_j,matmul(delta_alpha_j,r_up_j))

      ! The independently expanded form keeps the two transformed P terms
      ! separate before comparing with the raw-gamma-plus-correction form.
      expanded_i_ud=matmul(r_down_i,matmul(p_alpha(1:norb,1:norb),r_up_i))- &
         matmul(r_down_i,matmul(p_alpha(n2+1:n2+norb,n2+1:n2+norb),r_up_i))
      expanded_j_ud=matmul(r_up_j,matmul(p_alpha(norb+1:n2,norb+1:n2),r_down_j))- &
         matmul(r_up_j,matmul(p_alpha(n2+norb+1:n4,n2+norb+1:n4),r_down_j))
      expanded_i_du=matmul(r_up_i,matmul(p_alpha(1:norb,1:norb),r_down_i))- &
         matmul(r_up_i,matmul(p_alpha(n2+1:n2+norb,n2+1:n2+norb),r_down_i))
      expanded_j_du=matmul(r_down_j,matmul(p_alpha(norb+1:n2,norb+1:n2),r_up_j))- &
         matmul(r_down_j,matmul(p_alpha(n2+norb+1:n4,n2+norb+1:n4),r_up_j))

      correction_i=cmplx(0.0_rp,0.0_rp,rp)
      correction_j=cmplx(0.0_rp,0.0_rp,rp)
      ! D_sigma = alpha-gamma_sigma and R_sigma = I-D_sigma*P_gamma_sigma.
      ! For commuting diagonal local blocks,
      ! R_down*DeltaP_alpha*R_up
      !   = DeltaP_gamma + (D_up-D_down)*P_gamma_up*P_gamma_down
      !   = DeltaP_gamma + (gamma_down-gamma_up)*P_gamma_up*P_gamma_down.
      do l=0,lmax
         do m=1,2*l+1
            lm=l*l+m
            correction_i(lm,lm)=(lat%symbolic_atoms(1)%potential%qpar(l,2)- &
               lat%symbolic_atoms(1)%potential%qpar(l,1))*p_gamma(lm,lm)*p_gamma(n2+lm,n2+lm)
            correction_j(lm,lm)=(lat%symbolic_atoms(1)%potential%qpar(l,2)- &
               lat%symbolic_atoms(1)%potential%qpar(l,1))*p_gamma(norb+lm,norb+lm)*p_gamma(n2+norb+lm,n2+norb+lm)
         end do
      end do
      candidate_i=delta_gamma_i+correction_i
      candidate_j=delta_gamma_j+correction_j

      vertex_transform_value=max(maxval(abs(vertex_i_ud-expanded_i_ud)), &
         maxval(abs(vertex_j_ud-expanded_j_ud)),maxval(abs(vertex_i_du-expanded_i_du)), &
         maxval(abs(vertex_j_du-expanded_j_du)))
      expanded_vertex_value=max(maxval(abs(expanded_i_ud-candidate_i)), &
         maxval(abs(expanded_j_ud-candidate_j)),maxval(abs(expanded_i_du-candidate_i)), &
         maxval(abs(expanded_j_du-candidate_j)))
      vertex_transform_error=max(vertex_transform_error,vertex_transform_value)
      expanded_vertex_error=max(expanded_vertex_error,expanded_vertex_value)
      raw_gamma_vertex_magnitude=max(maxval(abs(delta_gamma_i)),maxval(abs(delta_gamma_j)))
      screening_correction_magnitude=max(maxval(abs(correction_i)),maxval(abs(correction_j)))
      transformed_vertex_magnitude=max(maxval(abs(vertex_i_ud)),maxval(abs(vertex_j_ud)), &
         maxval(abs(vertex_i_du)),maxval(abs(vertex_j_du)))

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

      native_alpha_ordered_ud=native_exchange_integrand(delta_alpha_i,delta_alpha_j,g_alpha(1:norb,norb+1:n2), &
         g_alpha(n2+norb+1:n4,n2+1:n2+norb))
      native_alpha_ordered_du=native_exchange_integrand(delta_alpha_i,delta_alpha_j, &
         g_alpha(n2+1:n2+norb,n2+norb+1:n4),g_alpha(norb+1:n2,1:norb))
      native_alpha_symmetrized=0.5_rp*(native_alpha_ordered_ud+native_alpha_ordered_du)

      native_gamma_ordered_ud=native_exchange_integrand(delta_gamma_i,delta_gamma_j,gup_ij,gdown_ji)
      native_gamma_ordered_du=native_exchange_integrand(delta_gamma_i,delta_gamma_j,gdown_ij,gup_ji)
      native_gamma_symmetrized=0.5_rp*(native_gamma_ordered_ud+native_gamma_ordered_du)

      transformed_gamma_ordered_ud=native_exchange_integrand(vertex_i_ud,vertex_j_ud,gup_ij,gdown_ji)
      transformed_gamma_ordered_du=native_exchange_integrand(vertex_i_du,vertex_j_du,gdown_ij,gup_ji)
      alpha_transformed_value=max(abs(native_alpha_ordered_ud-transformed_gamma_ordered_ud), &
         abs(native_alpha_ordered_du-transformed_gamma_ordered_du))
      raw_gamma_transformed_value=max(abs(native_gamma_ordered_ud-transformed_gamma_ordered_ud), &
         abs(native_gamma_ordered_du-transformed_gamma_ordered_du))
      alpha_transformed_error=max(alpha_transformed_error,alpha_transformed_value)
      raw_gamma_transformed_error=max(raw_gamma_transformed_error,raw_gamma_transformed_value)

      ! Map the transformed vertex into the coefficient-space units used by
      ! the historical finite-H d_matrix, without changing that oracle.
      d_tilde_i=matmul(d_width(1:norb,1:norb),matmul(vertex_i_ud,d_width(n2+1:n2+norb,n2+1:n2+norb)))
      d_tilde_j=matmul(d_width(norb+1:n2,norb+1:n2),matmul(vertex_j_ud, &
         d_width(n2+norb+1:n4,n2+norb+1:n4)))
      delta_d_i=d_tilde_i-d1
      delta_d_j=d_tilde_j-d2
      delta_d_value=max(maxval(abs(delta_d_i)),maxval(abs(delta_d_j)))
      delta_d_error=max(delta_d_error,delta_d_value)

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
      write(*,'(a,es14.6)') '  vertex_endpoint_expansion = ',vertex_transform_value
      write(*,'(a,es14.6)') '  vertex_screening_identity = ',expanded_vertex_value
      write(*,'(a,es14.6)') '  raw_DeltaP_gamma_magnitude = ',raw_gamma_vertex_magnitude
      write(*,'(a,es14.6)') '  screening_correction_magnitude = ',screening_correction_magnitude
      write(*,'(a,es14.6)') '  transformed_vertex_magnitude = ',transformed_vertex_magnitude
      write(*,'(a,es14.6)') '  alpha_minus_transformed_gamma_ud = ',abs(native_alpha_ordered_ud-transformed_gamma_ordered_ud)
      write(*,'(a,es14.6)') '  alpha_minus_transformed_gamma_du = ',abs(native_alpha_ordered_du-transformed_gamma_ordered_du)
      write(*,'(a,es14.6)') '  raw_gamma_minus_transformed_gamma_ud = ',abs(native_gamma_ordered_ud-transformed_gamma_ordered_ud)
      write(*,'(a,es14.6)') '  raw_gamma_minus_transformed_gamma_du = ',abs(native_gamma_ordered_du-transformed_gamma_ordered_du)
      write(*,'(a,es14.6)') '  delta_d_screen = d_tilde - d (max abs) = ',delta_d_value
      write(*,'(a,es14.6)') '  raw_gamma_minus_finite_H_du = ',abs(native_gamma_ordered_du-finite_h_ordered_du)

      ! Isolated common-gamma control: retain the physical, spin-dependent
      ! P_gamma blocks (including their c/dele dependence), but make the two
      ! local screening matrices equal.  This changes only local diagnostic
      ! arrays; no live potential input is modified.
      d_screen=cmplx(0.0_rp,0.0_rp,rp)
      do isite=1,nsite_fixture
         do l=0,lmax
            do m=1,2*l+1
               lm=l*l+m
               d_screen((isite-1)*norb+lm,(isite-1)*norb+lm)=alpha(l)- &
                  lat%symbolic_atoms(1)%potential%qpar(l,1)
               d_screen(n2+(isite-1)*norb+lm,n2+(isite-1)*norb+lm)=alpha(l)- &
                  lat%symbolic_atoms(1)%potential%qpar(l,1)
            end do
         end do
      end do
      r_screen=eye-matmul(d_screen,p_gamma)
      call native_inverse(r_screen,r_inverse)
      p_alpha=matmul(p_gamma,r_inverse)
      call native_inverse(eye+matmul(s_alpha,d_screen),tmp)
      s_gamma=matmul(tmp,s_alpha)
      call native_inverse(p_gamma-s_gamma,g_gamma)
      call native_inverse(p_alpha-s_alpha,g_alpha)
      gup_ij=g_gamma(1:norb,norb+1:n2)
      gup_ji=g_gamma(norb+1:n2,1:norb)
      gdown_ij=g_gamma(n2+1:n2+norb,n2+norb+1:n4)
      gdown_ji=g_gamma(n2+norb+1:n4,n2+1:n2+norb)
      delta_alpha_i=p_alpha(1:norb,1:norb)-p_alpha(n2+1:n2+norb,n2+1:n2+norb)
      delta_alpha_j=p_alpha(norb+1:n2,norb+1:n2)-p_alpha(n2+norb+1:n4,n2+norb+1:n4)
      delta_gamma_i=p_gamma(1:norb,1:norb)-p_gamma(n2+1:n2+norb,n2+1:n2+norb)
      delta_gamma_j=p_gamma(norb+1:n2,norb+1:n2)-p_gamma(n2+norb+1:n4,n2+norb+1:n4)
      r_up_i=r_screen(1:norb,1:norb)
      r_up_j=r_screen(norb+1:n2,norb+1:n2)
      r_down_i=r_screen(n2+1:n2+norb,n2+1:n2+norb)
      r_down_j=r_screen(n2+norb+1:n4,n2+norb+1:n4)
      vertex_i_ud=matmul(r_down_i,matmul(delta_alpha_i,r_up_i))
      vertex_j_ud=matmul(r_up_j,matmul(delta_alpha_j,r_down_j))
      vertex_i_du=matmul(r_up_i,matmul(delta_alpha_i,r_down_i))
      vertex_j_du=matmul(r_down_j,matmul(delta_alpha_j,r_up_j))
      common_gamma_vertex_value=max(maxval(abs(vertex_i_ud-delta_gamma_i)), &
         maxval(abs(vertex_j_ud-delta_gamma_j)),maxval(abs(vertex_i_du-delta_gamma_i)), &
         maxval(abs(vertex_j_du-delta_gamma_j)))
      common_gamma_vertex_error=max(common_gamma_vertex_error,common_gamma_vertex_value)
      correction_i=cmplx(0.0_rp,0.0_rp,rp)
      correction_j=cmplx(0.0_rp,0.0_rp,rp)
      do l=0,lmax
         do m=1,2*l+1
            lm=l*l+m
            correction_i(lm,lm)=(d_screen(lm,lm)-d_screen(n2+lm,n2+lm))* &
               p_gamma(lm,lm)*p_gamma(n2+lm,n2+lm)
            correction_j(lm,lm)=(d_screen(norb+lm,norb+lm)-d_screen(n2+norb+lm,n2+norb+lm))* &
               p_gamma(norb+lm,norb+lm)*p_gamma(n2+norb+lm,n2+norb+lm)
         end do
      end do
      common_gamma_correction_magnitude=max(maxval(abs(correction_i)),maxval(abs(correction_j)))
      common_gamma_correction_error=max(common_gamma_correction_error,common_gamma_correction_magnitude)
      common_gamma_alpha_ud=native_exchange_integrand(delta_alpha_i,delta_alpha_j, &
         g_alpha(1:norb,norb+1:n2),g_alpha(n2+norb+1:n4,n2+1:n2+norb))
      common_gamma_alpha_du=native_exchange_integrand(delta_alpha_i,delta_alpha_j, &
         g_alpha(n2+1:n2+norb,n2+norb+1:n4),g_alpha(norb+1:n2,1:norb))
      common_gamma_raw_ud=native_exchange_integrand(delta_gamma_i,delta_gamma_j,gup_ij,gdown_ji)
      common_gamma_raw_du=native_exchange_integrand(delta_gamma_i,delta_gamma_j,gdown_ij,gup_ji)
      common_gamma_error=max(common_gamma_error,abs(common_gamma_alpha_ud-common_gamma_raw_ud), &
         abs(common_gamma_alpha_du-common_gamma_raw_du))
      write(*,'(a,es14.6)') '  common_gamma_screening_correction = ',common_gamma_correction_magnitude
      write(*,'(a,es14.6)') '  common_gamma_vertex_residual = ',common_gamma_vertex_value
      write(*,'(a,es14.6)') '  common_gamma_alpha_ordered_ud = ',common_gamma_alpha_ud
      write(*,'(a,es14.6)') '  common_gamma_raw_ordered_ud = ',common_gamma_raw_ud
      write(*,'(a,es14.6)') '  common_gamma_alpha_minus_raw_ud = ',abs(common_gamma_alpha_ud-common_gamma_raw_ud)
      write(*,'(a,es14.6)') '  common_gamma_alpha_ordered_du = ',common_gamma_alpha_du
      write(*,'(a,es14.6)') '  common_gamma_raw_ordered_du = ',common_gamma_raw_du
      write(*,'(a,es14.6)') '  common_gamma_alpha_minus_raw_du = ',abs(common_gamma_alpha_du-common_gamma_raw_du)
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
   write(*,'(a,es14.6)') 'DRESP-03TG vertex endpoint-expansion residual = ',vertex_transform_error
   write(*,'(a,es14.6)') 'DRESP-03TG vertex screening-identity residual = ',expanded_vertex_error
   write(*,'(a,es14.6)') 'DRESP-03TG alpha direct/transformed-gamma residual = ',alpha_transformed_error
   write(*,'(a,es14.6)') 'DRESP-03TG raw-gamma/transformed-gamma residual = ',raw_gamma_transformed_error
   write(*,'(a,es14.6)') 'DRESP-03TG coefficient delta_d_screen max = ',delta_d_error
   write(*,'(a,es14.6)') 'DRESP-03TG common-gamma vertex residual = ',common_gamma_vertex_error
   write(*,'(a,es14.6)') 'DRESP-03TG common-gamma screening correction max = ',common_gamma_correction_error
   write(*,'(a,es14.6)') 'DRESP-03TG common-gamma contraction residual = ',common_gamma_error
   call run_r4_finite_h_bridge(lat, s_alpha, r4_fixture, r4_max_fd_error, r4_max_contact_norm, &
      r4_max_vertex_raw_error, r4_max_vertex_transformed_error)
   failed = p_error > 2.0e-14_rp .or. d_error > 2.0e-13_rp .or. solve_error > 2.0e-12_rp .or. &
      p_transform_error > 2.0e-12_rp .or. s_transform_error > 2.0e-12_rp .or. gf_cov_error > 2.0e-12_rp .or. &
      gamma_h_error > 2.0e-10_rp .or. gamma_h_du_error > 2.0e-10_rp .or. pauli_self_error > 2.0e-12_rp .or. &
      vertex_transform_error > 2.0e-12_rp .or. expanded_vertex_error > 2.0e-12_rp .or. &
      alpha_transformed_error > 2.0e-10_rp .or. common_gamma_vertex_error > 2.0e-12_rp .or. &
      common_gamma_correction_error > 2.0e-12_rp .or. common_gamma_error > 2.0e-10_rp .or. &
      raw_gamma_transformed_error < 1.0e-12_rp
   if (failed) then
      write(*,'(a)') 'RESULT: BLOCKED — TG-FZ-R2 SPIN-SCREENING VERTEX CHECK'
      return
   end if
   write(*,'(a)') 'TG-FZ-R2 SPIN-SCREENING VERTEX COVARIANCE: PASS-A'

contains

   ! TG-FZ-R4 is deliberately a diagnostic extension of the R2 program.  The
   ! native P/S objects and the finite-H objects are kept on their own paths:
   ! the former use the two-site screened S matrices constructed above, while
   ! the latter use the certified live B-QB+E_nu fixture below.
   subroutine run_r4_finite_h_bridge(lat, s_alpha_in, fix, max_fd_error, max_contact_norm, &
                                     max_raw_error, max_transformed_error)
      type(lattice), intent(in) :: lat
      complex(rp), intent(in) :: s_alpha_in(:,:)
      type(lmto_live_hamiltonian_fixture), intent(out) :: fix
      real(rp), intent(out) :: max_fd_error, max_contact_norm, max_raw_error, max_transformed_error
      real(rp), parameter :: kpoint(3) = [0.0_rp, 0.0_rp, 0.0_rp]
      real(rp), parameter :: axis(3) = [0.0_rp, 1.0_rp, 0.0_rp]
      real(rp), parameter :: eps_values(6) = [1.0e-2_rp, 1.0e-3_rp, 1.0e-4_rp, 1.0e-5_rp, 1.0e-6_rp, 1.0e-7_rp]
      integer :: nmat, n2loc, iz, site, ie, l, m, lm
      complex(rp), allocatable :: b(:,:), q(:,:), enu(:,:), bi(:,:), qi(:,:), enui(:,:), torque(:,:)
      complex(rp), allocatable :: tij(:,:), tj(:,:), cij(:,:), h(:,:), gh(:,:), eye_r4(:,:)
      complex(rp), allocatable :: pga(:,:), paa(:,:), sga(:,:), rmat(:,:), gng(:,:), gna(:,:)
      complex(rp), allocatable :: bi_i(:,:), qi_i(:,:), enui_i(:,:), bi_j(:,:), qi_j(:,:), enui_j(:,:)
      complex(rp), allocatable :: dg_i(:,:), dg_j(:,:), da_i(:,:), da_j(:,:), raw_i(:,:), raw_j(:,:), &
         tilde_i(:,:), tilde_j(:,:), qcorr_i(:,:), qcorr_j(:,:), tud_i(:,:), tdu_j(:,:)
      complex(rp), allocatable :: gup_ij(:,:), gup_ji(:,:), gdown_ij(:,:), gdown_ji(:,:)
      complex(rp), allocatable :: hplus(:,:), hminus(:,:)
      real(rp), allocatable :: saved_mom(:,:)
      real(rp) :: raw_norm, screen_norm, tilde_norm, t_norm, bi_norm, q_norm, qib_norm, qbi_norm, enui_norm
      real(rp) :: raw_residual, tilde_residual, q_screen_residual, fd_error, tt, contact, complete
      real(rp) :: raw_contract, alpha_contract, transformed_contract, actual_contract, full_tt, full_contact, full_complete
      real(rp) :: native_h_error, split_h_error, decomp_error, min_candidate_residual

      n2loc = nsite_fixture*norb
      nmat = 2*n2loc
      if (size(s_alpha_in,1) /= nmat .or. size(s_alpha_in,2) /= nmat) then
         error stop 'TG-FZ-R4: screening fixture shape mismatch'
      end if

      ! This is the same two-site sbar-derived fixture as R1/R2, now retained
      ! in its common-screened representation and fed to the certified live
      ! second-order H assembler.  No candidate vertex is used to build H.
      call lmto_fixture_init(fix, nsite_fixture, norb, 4, .true.)
      fix%include_enu = .true.
      fix%cartesian_to_spherical = .false.
      fix%moments(:,1) = [0.0_rp,0.0_rp,1.0_rp]
      fix%moments(:,2) = [0.0_rp,0.0_rp,1.0_rp]
      fix%bond_source = [1,2,1,2]
      fix%bond_target = [1,2,2,1]
      fix%onsite = [.true.,.true.,.false.,.false.]
      fix%bond_vector = 0.0_rp
      fix%hhh(:,:,1) = s_alpha_in(1:norb,1:norb)
      fix%hhh(:,:,2) = s_alpha_in(norb+1:n2loc,norb+1:n2loc)
      fix%hhh(:,:,3) = s_alpha_in(1:norb,norb+1:n2loc)
      fix%hhh(:,:,4) = s_alpha_in(norb+1:n2loc,1:norb)
      do site = 1, nsite_fixture
         ! predls writes the transformed l-channel arrays.  Expand those
         ! arrays into the orbital channels exactly as build_pot/build_obarm
         ! do; the pre-predls obx0/obx1 cache is intentionally not reused.
         do l = 0, lmax
            do m = 1, 2*l + 1
               lm = l*l + m
               fix%wx0(lm,site) = 0.5_rp*(lat%symbolic_atoms(1)%potential%width_band(l+1,1) + &
                  lat%symbolic_atoms(1)%potential%width_band(l+1,2))
               fix%wx1(lm,site) = 0.5_rp*(lat%symbolic_atoms(1)%potential%width_band(l+1,1) - &
                  lat%symbolic_atoms(1)%potential%width_band(l+1,2))
               fix%c0(lm,site) = 0.5_rp*(lat%symbolic_atoms(1)%potential%shifted_band(l+1,1) + &
                  lat%symbolic_atoms(1)%potential%shifted_band(l+1,2))
               fix%c1(lm,site) = 0.5_rp*(lat%symbolic_atoms(1)%potential%shifted_band(l+1,1) - &
                  lat%symbolic_atoms(1)%potential%shifted_band(l+1,2))
               fix%obar0(lm,site) = 0.5_rp*(lat%symbolic_atoms(1)%potential%obar(l+1,1) + &
                  lat%symbolic_atoms(1)%potential%obar(l+1,2))
               fix%obar1(lm,site) = 0.5_rp*(lat%symbolic_atoms(1)%potential%obar(l+1,1) - &
                  lat%symbolic_atoms(1)%potential%obar(l+1,2))
               fix%enu0(lm,site) = 0.5_rp*((lat%symbolic_atoms(1)%potential%center_band(l+1,1) - &
                  lat%symbolic_atoms(1)%potential%shifted_band(l+1,1)) + &
                  (lat%symbolic_atoms(1)%potential%center_band(l+1,2) - &
                  lat%symbolic_atoms(1)%potential%shifted_band(l+1,2)))
               fix%enu1(lm,site) = 0.5_rp*((lat%symbolic_atoms(1)%potential%center_band(l+1,1) - &
                  lat%symbolic_atoms(1)%potential%shifted_band(l+1,1)) - &
                  (lat%symbolic_atoms(1)%potential%center_band(l+1,2) - &
                  lat%symbolic_atoms(1)%potential%shifted_band(l+1,2)))
            end do
         end do
      end do

      allocate(b(nmat,nmat),q(nmat,nmat),enu(nmat,nmat),bi(nmat,nmat),qi(nmat,nmat),enui(nmat,nmat),torque(nmat,nmat), &
         tij(nmat,nmat),tj(nmat,nmat),cij(nmat,nmat),h(nmat,nmat),gh(nmat,nmat),eye_r4(nmat,nmat), &
         pga(nmat,nmat),paa(nmat,nmat),sga(nmat,nmat),rmat(nmat,nmat),gng(nmat,nmat),gna(nmat,nmat), &
         bi_i(nmat,nmat),qi_i(nmat,nmat),enui_i(nmat,nmat),bi_j(nmat,nmat),qi_j(nmat,nmat),enui_j(nmat,nmat), &
         dg_i(norb,norb),dg_j(norb,norb),da_i(norb,norb),da_j(norb,norb),raw_i(norb,norb),raw_j(norb,norb), &
         tilde_i(norb,norb),tilde_j(norb,norb),qcorr_i(norb,norb),qcorr_j(norb,norb), &
         tud_i(norb,norb),tdu_j(norb,norb),gup_ij(norb,norb),gup_ji(norb,norb),gdown_ij(norb,norb),gdown_ji(norb,norb), &
         hplus(nmat,nmat),hminus(nmat,nmat),saved_mom(3,nsite_fixture))
      eye_r4 = identity(nmat)

      ! Independent finite-difference oracle for the actual assembled H.  The
      ! analytic path and the FD path share only the fixture inputs, not a
      ! candidate Turek vertex.
      saved_mom = fix%moments
      max_fd_error = 0.0_rp
      do ie = 1, size(eps_values)
         do site = 1, nsite_fixture
            fix%moments = saved_mom
            call assemble_lmto_rotation_terms(fix,kpoint,site,axis,b,q,enu,bi,qi,enui,torque)
            fix%moments(:,site) = r4_rotate_moment(saved_mom(:,site),axis,eps_values(ie))
            call assemble_lmto_hamiltonian(fix,kpoint,hplus)
            fix%moments(:,site) = r4_rotate_moment(saved_mom(:,site),axis,-eps_values(ie))
            call assemble_lmto_hamiltonian(fix,kpoint,hminus)
            fd_error = maxval(abs((hplus-hminus)/(2.0_rp*eps_values(ie))-torque))
            ! The two coarsest steps intentionally expose the expected O(eps^2)
            ! central-difference truncation.  The oracle gate uses the
            ! converged eps <= 1e-4 subset while retaining every value in the
            ! printed audit table.
            if (eps_values(ie) <= 1.0e-4_rp) max_fd_error = max(max_fd_error,fd_error)
            write(*,'(a,es10.3,a,i0,a,es14.6)') 'R4 FD epsilon=',eps_values(ie),' site=',site,' residual=',fd_error
         end do
      end do
      fix%moments = saved_mom
      if (max_fd_error > 5.0e-8_rp) error stop 'TG-FZ-R4: finite-H rotation derivative oracle failed'
      write(*,'(a,es14.6)') 'TG-FZ-R4 FD derivative max residual = ',max_fd_error

      call assemble_lmto_mixed_derivative(fix,kpoint,1,axis,2,axis,cij)
      max_contact_norm = sqrt(sum(abs(cij)**2))
      write(*,'(a,es14.6)') 'TG-FZ-R4 ||C_12||_F = ',max_contact_norm
      if (max_contact_norm <= 1.0e-12_rp) error stop 'TG-FZ-R4: expected nonzero offsite contact derivative'

      max_raw_error = 0.0_rp
      max_transformed_error = 0.0_rp
      native_h_error = 0.0_rp
      split_h_error = 0.0_rp
      decomp_error = 0.0_rp
      min_candidate_residual = huge(1.0_rp)
      z_values = [cmplx(-0.91_rp,0.83_rp,rp),cmplx(-0.17_rp,0.04_rp,rp),cmplx(0.62_rp,0.31_rp,rp)]
      do iz = 1, size(z_values)
         z = z_values(iz)
         call r4_build_native_state(lat,z,s_alpha_in,pga,paa,sga,rmat,gng,gna)
         call assemble_lmto_hamiltonian(fix,kpoint,h)
         call native_inverse(z*eye_r4-h,gh)
         call assemble_lmto_rotation_terms(fix,kpoint,1,axis,b,q,enu,bi_i,qi_i,enui_i,tij)
         call assemble_lmto_rotation_terms(fix,kpoint,2,axis,b,q,enu,bi_j,qi_j,enui_j,tj)
         call extract_local_spin_flip(tij,1,1,tud_i)
         call extract_local_spin_flip(tj,2,2,tdu_j)
         call extract_local_spin_flip(bi_i,1,1,raw_i)
         call extract_local_spin_flip(bi_j,2,2,raw_j)
         qcorr_i = -extract_r4_local_product(qi_i,b,1,1) - extract_r4_local_product(q,bi_i,1,1)
         qcorr_j = -extract_r4_local_product(qi_j,b,2,2) - extract_r4_local_product(q,bi_j,2,2)

         ! Keep the local spin-flip blocks of every live product-rule term
         ! before raw_i/raw_j are reused for the candidate vertices below.
         call extract_local_spin_flip(bi_i,1,1,raw_i)
         call extract_local_spin_flip(bi_j,2,2,raw_j)
         call extract_local_spin_flip(enui_i,1,1,tilde_i)
         call extract_local_spin_flip(enui_j,2,2,tilde_j)
         bi_norm = max(sqrt(sum(abs(raw_i)**2)),sqrt(sum(abs(raw_j)**2)))
         q_norm = sqrt(sum(abs(q)**2))
         qib_norm = max(sqrt(sum(abs(extract_r4_local_product(qi_i,b,1,1))**2)), &
            sqrt(sum(abs(extract_r4_local_product(qi_j,b,2,2))**2)))
         qbi_norm = max(sqrt(sum(abs(extract_r4_local_product(q,bi_i,1,1))**2)), &
            sqrt(sum(abs(extract_r4_local_product(q,bi_j,2,2))**2)))
         enui_norm = max(sqrt(sum(abs(tilde_i)**2)),sqrt(sum(abs(tilde_j)**2)))
         t_norm = max(sqrt(sum(abs(tud_i)**2)),sqrt(sum(abs(tdu_j)**2)))
         decomp_error = max(decomp_error, max(maxval(abs(tud_i-(raw_i+qcorr_i+tilde_i))), &
            maxval(abs(tdu_j-(raw_j+qcorr_j+tilde_j)))))
         write(*,'(a,2es14.6)') 'z = ',real(z,rp),aimag(z)
         write(*,'(a,es14.6)') '  ||B_i^flip||_F = ',bi_norm
         write(*,'(a,es14.6)') '  ||Q||_F = ',q_norm
         write(*,'(a,es14.6)') '  ||(-Q_i B)^flip||_F = ',qib_norm
         write(*,'(a,es14.6)') '  ||(-Q B_i)^flip||_F = ',qbi_norm
         write(*,'(a,es14.6)') '  ||E_nu,i^flip||_F = ',enui_norm
         write(*,'(a,es14.6)') '  ||T_i^flip||_F = ',t_norm
         write(*,'(a,es14.6)') '  product-rule decomposition residual = ', &
            max(maxval(abs(tud_i-(raw_i+qcorr_i+tilde_i))),maxval(abs(tdu_j-(raw_j+qcorr_j+tilde_j))))

         call r4_native_vertices(lat,z,pga,paa,rmat,raw_i,raw_j,tilde_i,tilde_j,da_i,da_j,dg_i,dg_j)
         ! At theta_y around +z the live hhmag convention gives Hx-iHy
         ! without an extra i or 1/2.  The opposite block uses Hx+iHy.
         raw_residual = max(maxval(abs(tud_i-raw_i)),maxval(abs(tdu_j-raw_j)))
         tilde_residual = max(maxval(abs(tud_i-tilde_i)),maxval(abs(tdu_j-tilde_j)))
         q_screen_residual = max(maxval(abs(qcorr_i - (tilde_i-raw_i))), &
            maxval(abs(qcorr_j - (tilde_j-raw_j))))
         max_raw_error = max(max_raw_error,raw_residual)
         max_transformed_error = max(max_transformed_error,tilde_residual)
         min_candidate_residual = min(min_candidate_residual,raw_residual,tilde_residual)

         raw_norm = max(sqrt(sum(abs(raw_i)**2)),sqrt(sum(abs(raw_j)**2)))
         screen_norm = max(sqrt(sum(abs(tilde_i-raw_i)**2)),sqrt(sum(abs(tilde_j-raw_j)**2)))
         tilde_norm = max(sqrt(sum(abs(tilde_i)**2)),sqrt(sum(abs(tilde_j)**2)))
         write(*,'(a,es14.6)') '  ||d_raw||_F = ',raw_norm
         write(*,'(a,es14.6)') '  ||delta_d_screen||_F = ',screen_norm
         write(*,'(a,es14.6)') '  ||d_tilde_alpha||_F = ',tilde_norm
         write(*,'(a,es14.6)') '  ||T_i^ud-d_raw||_max = ',raw_residual
         write(*,'(a,es14.6)') '  ||T_i^ud-d_tilde_alpha||_max = ',tilde_residual
         write(*,'(a,es14.6)') '  Q-term vs screening-correction residual = ',q_screen_residual

         gup_ij = gng(1:norb,norb+1:n2loc)
         gup_ji = gng(norb+1:n2loc,1:norb)
         gdown_ij = gng(n2loc+1:n2loc+norb,n2loc+norb+1:nmat)
         gdown_ji = gng(n2loc+norb+1:nmat,n2loc+1:n2loc+norb)
         ! pga is replaced by the actual finite-H Green below; retain the
         ! native blocks first for the four representation-specific values.
         raw_contract = native_exchange_integrand(dg_i,dg_j,gup_ij,gdown_ji)
         alpha_contract = native_exchange_integrand(da_i,da_j, &
            gna(1:norb,norb+1:n2loc),gna(n2loc+norb+1:nmat,n2loc+1:n2loc+norb))
         transformed_contract = native_exchange_integrand( &
            matmul(rmat(n2loc+1:n2loc+norb,n2loc+1:n2loc+norb),matmul(da_i, &
               rmat(1:norb,1:norb))),matmul(rmat(norb+1:n2loc,norb+1:n2loc),matmul(da_j, &
               rmat(n2loc+norb+1:nmat,n2loc+norb+1:nmat))),gup_ij,gdown_ji)

         call extract_site_spin_block(gh,norb,1,1,2,gup_ij)
         call extract_site_spin_block(gh,norb,2,1,2,gdown_ij)
         call extract_site_spin_block(gh,norb,1,2,1,gup_ji)
         call extract_site_spin_block(gh,norb,2,2,1,gdown_ji)
         actual_contract = aimag(trace4(tud_i,gdown_ij,tdu_j,gup_ji))
         call force_theorem_integrand(tij,gh,tj,cij,tt,contact,complete)
         full_tt = -aimag(trace4(tij,gh,tj,gh))/force_theorem_pi
         full_contact = -aimag(trace2(cij,gh))/force_theorem_pi
         full_complete = full_tt+full_contact
         write(*,'(a,es14.6)') '  raw-gamma contraction = ',raw_contract
         write(*,'(a,es14.6)') '  alpha-direct contraction = ',alpha_contract
         write(*,'(a,es14.6)') '  transformed-gamma contraction = ',transformed_contract
         write(*,'(a,es14.6)') '  actual-T ud-block contraction = ',actual_contract
         write(*,'(a,es14.6)') '  TT contribution (-Im/pi, full T) = ',full_tt
         write(*,'(a,es14.6)') '  contact contribution (-Im/pi) = ',full_contact
         write(*,'(a,es14.6)') '  complete contribution (-Im/pi) = ',full_complete
         native_h_error = max(native_h_error,abs(alpha_contract-transformed_contract))
         split_h_error = max(split_h_error,abs(full_complete-(full_tt+full_contact)))
      end do
      write(*,'(a,es14.6)') 'TG-FZ-R4 max T-vs-raw-d residual = ',max_raw_error
      write(*,'(a,es14.6)') 'TG-FZ-R4 max T-vs-transformed-alpha-d residual = ',max_transformed_error
      write(*,'(a,es14.6)') 'TG-FZ-R4 max native alpha/transformed residual = ',native_h_error
      write(*,'(a,es14.6)') 'TG-FZ-R4 max product-rule decomposition residual = ',decomp_error
      write(*,'(a,es14.6)') 'TG-FZ-R4 max TT/contact split residual = ',split_h_error
      if (max_transformed_error <= 2.0e-10_rp) then
         write(*,'(a)') 'TG-FZ-R4 verdict: PASS-A'
      else if (max_raw_error <= 2.0e-10_rp) then
         write(*,'(a)') 'TG-FZ-R4 verdict: PASS-B'
      else
         write(*,'(a)') 'TG-FZ-R4 verdict: BLOCKED — no unfitted finite-H/Turek vertex identity'
      end if
      write(*,'(a,es14.6)') 'TG-FZ-R4 minimum tested vertex residual = ',min_candidate_residual
      fix%moments = saved_mom
      deallocate(b,q,enu,bi,qi,enui,torque,tij,tj,cij,h,gh,eye_r4,pga,paa,sga,rmat,gng,gna, &
         bi_i,qi_i,enui_i,bi_j,qi_j,enui_j, &
         dg_i,dg_j,da_i,da_j,raw_i,raw_j,tilde_i,tilde_j,qcorr_i,qcorr_j,tud_i,tdu_j, &
         gup_ij,gup_ji,gdown_ij,gdown_ji,hplus,hminus,saved_mom)
   end subroutine run_r4_finite_h_bridge

   subroutine r4_build_native_state(atom_lattice,zloc,s_alpha_loc,p_gamma_loc,p_alpha_loc,s_gamma_loc,r_loc,g_gamma_loc,g_alpha_loc)
      type(lattice), intent(in) :: atom_lattice
      complex(rp), intent(in) :: zloc,s_alpha_loc(:,:)
      complex(rp), intent(out) :: p_gamma_loc(:,:),p_alpha_loc(:,:),s_gamma_loc(:,:),r_loc(:,:),g_gamma_loc(:,:),g_alpha_loc(:,:)
      complex(rp), allocatable :: psite(:,:), asite(:,:), dloc(:,:), eye_loc(:,:), inv_loc(:,:)
      real(rp) :: alpha_loc(0:lmax)
      integer :: site, lm, l, m, nloc
      nloc = size(p_gamma_loc,1)
      allocate(psite(2*norb,2*norb),asite(2*norb,2*norb),dloc(nloc,nloc),eye_loc(nloc,nloc),inv_loc(nloc,nloc))
      call native_screening_alpha(atom_lattice%symbolic_atoms(1),alpha_loc)
      p_gamma_loc=cmplx(0.0_rp,0.0_rp,rp); p_alpha_loc=p_gamma_loc; dloc=p_gamma_loc
      do site=1,nsite_fixture
         call native_complex_p_matrix(atom_lattice%symbolic_atoms(1),zloc,psite)
         call native_screened_p_matrix(atom_lattice%symbolic_atoms(1),zloc,alpha_loc,asite)
         p_gamma_loc((site-1)*norb+1:site*norb,(site-1)*norb+1:site*norb)=psite(1:norb,1:norb)
         p_gamma_loc(nsite_fixture*norb+(site-1)*norb+1:nsite_fixture*norb+site*norb, &
            nsite_fixture*norb+(site-1)*norb+1:nsite_fixture*norb+site*norb)=psite(norb+1:2*norb,norb+1:2*norb)
         p_alpha_loc((site-1)*norb+1:site*norb,(site-1)*norb+1:site*norb)=asite(1:norb,1:norb)
         p_alpha_loc(nsite_fixture*norb+(site-1)*norb+1:nsite_fixture*norb+site*norb, &
            nsite_fixture*norb+(site-1)*norb+1:nsite_fixture*norb+site*norb)=asite(norb+1:2*norb,norb+1:2*norb)
         do l=0,lmax
            do m=1,2*l+1
               lm=l*l+m
               dloc((site-1)*norb+lm,(site-1)*norb+lm)=alpha_loc(l)-atom_lattice%symbolic_atoms(1)%potential%qpar(l,1)
               dloc(nsite_fixture*norb+(site-1)*norb+lm,nsite_fixture*norb+(site-1)*norb+lm)= &
                  alpha_loc(l)-atom_lattice%symbolic_atoms(1)%potential%qpar(l,2)
            end do
         end do
      end do
      eye_loc=identity(nloc)
      call native_inverse(eye_loc+matmul(s_alpha_loc,dloc),inv_loc)
      s_gamma_loc=matmul(inv_loc,s_alpha_loc)
      r_loc=eye_loc-matmul(dloc,p_gamma_loc)
      call native_inverse(r_loc,inv_loc)
      p_alpha_loc=matmul(p_gamma_loc,inv_loc)
      call native_inverse(p_gamma_loc-s_gamma_loc,g_gamma_loc)
      call native_inverse(p_alpha_loc-s_alpha_loc,g_alpha_loc)
      deallocate(psite,asite,dloc,eye_loc,inv_loc)
   end subroutine r4_build_native_state

   subroutine r4_native_vertices(atom_lattice,zloc,p_gamma_loc,p_alpha_loc,r_loc,raw_i,raw_j,tilde_i,tilde_j, &
                                 da_i,da_j,dg_i,dg_j)
      type(lattice), intent(in) :: atom_lattice
      complex(rp), intent(in) :: zloc,p_gamma_loc(:,:),p_alpha_loc(:,:),r_loc(:,:)
      complex(rp), intent(out) :: raw_i(:,:),raw_j(:,:),tilde_i(:,:),tilde_j(:,:),da_i(:,:),da_j(:,:),dg_i(:,:),dg_j(:,:)
      complex(rp) :: rup_i(norb,norb),rdn_i(norb,norb),rup_j(norb,norb),rdn_j(norb,norb)
      complex(rp) :: wup(norb,norb),wdn(norb,norb)
      real(rp) :: alpha_loc(0:lmax)
      integer :: l,m,lm,n2local,nmatlocal
      n2local=nsite_fixture*norb; nmatlocal=2*n2local
      call native_screening_alpha(atom_lattice%symbolic_atoms(1),alpha_loc)
      dg_i=p_gamma_loc(1:norb,1:norb)-p_gamma_loc(n2local+1:n2local+norb,n2local+1:n2local+norb)
      dg_j=p_gamma_loc(norb+1:n2local,norb+1:n2local)-p_gamma_loc(n2local+norb+1:nmatlocal,n2local+norb+1:nmatlocal)
      da_i=p_alpha_loc(1:norb,1:norb)-p_alpha_loc(n2local+1:n2local+norb,n2local+1:n2local+norb)
      da_j=p_alpha_loc(norb+1:n2local,norb+1:n2local)-p_alpha_loc(n2local+norb+1:nmatlocal,n2local+norb+1:nmatlocal)
      rup_i=r_loc(1:norb,1:norb); rup_j=r_loc(norb+1:n2local,norb+1:n2local)
      rdn_i=r_loc(n2local+1:n2local+norb,n2local+1:n2local+norb)
      rdn_j=r_loc(n2local+norb+1:nmatlocal,n2local+norb+1:nmatlocal)
      tilde_i=matmul(rdn_i,matmul(da_i,rup_i)); tilde_j=matmul(rup_j,matmul(da_j,rdn_j))
      wup=cmplx(0.0_rp,0.0_rp,rp); wdn=wup
      do l=0,lmax
         do m=1,2*l+1
            lm=l*l+m
            wup(lm,lm)=atom_lattice%symbolic_atoms(1)%potential%dele(l,1)
            wdn(lm,lm)=atom_lattice%symbolic_atoms(1)%potential%dele(l,2)
         end do
      end do
      raw_i=matmul(wup,matmul(dg_i,wdn)); raw_j=matmul(wdn,matmul(dg_j,wup))
      tilde_i=matmul(wup,matmul(tilde_i,wdn)); tilde_j=matmul(wdn,matmul(tilde_j,wup))
   end subroutine r4_native_vertices

   subroutine extract_local_spin_flip(matrix,site,spin_block,block)
      complex(rp), intent(in) :: matrix(:,:)
      integer, intent(in) :: site,spin_block
      complex(rp), intent(out) :: block(:,:)
      integer :: base
      base=(site-1)*2*norb
      if (spin_block == 1) then
         block=matrix(base+1:base+norb,base+norb+1:base+2*norb)
      else
         block=matrix(base+norb+1:base+2*norb,base+1:base+norb)
      end if
   end subroutine extract_local_spin_flip

   subroutine extract_site_spin_block(matrix,norb_in,spin,site_i,site_j,block)
      complex(rp), intent(in) :: matrix(:,:)
      integer, intent(in) :: norb_in,spin,site_i,site_j
      complex(rp), intent(out) :: block(:,:)
      integer :: row0,col0
      row0=(site_i-1)*2*norb_in+(spin-1)*norb_in
      col0=(site_j-1)*2*norb_in+(spin-1)*norb_in
      block=matrix(row0+1:row0+norb_in,col0+1:col0+norb_in)
   end subroutine extract_site_spin_block

   function extract_r4_local_product(left,right,site,spin_block) result(block)
      complex(rp), intent(in) :: left(:,:),right(:,:)
      integer, intent(in) :: site,spin_block
      complex(rp) :: block(norb,norb), product(size(left,1),size(left,2))
      integer :: base
      product=matmul(left,right)
      base=(site-1)*2*norb
      if (spin_block == 1) then
         block=product(base+1:base+norb,base+norb+1:base+2*norb)
      else
         block=product(base+norb+1:base+2*norb,base+1:base+norb)
      end if
   end function extract_r4_local_product

   pure function trace2(a,b) result(value)
      complex(rp), intent(in) :: a(:,:),b(:,:)
      complex(rp) :: value,product(size(a,1),size(a,2))
      integer :: i
      product=matmul(a,b)
      value=sum([(product(i,i),i=1,size(a,1))])
   end function trace2

   pure function trace4(a,b,c,d) result(value)
      complex(rp), intent(in) :: a(:,:),b(:,:),c(:,:),d(:,:)
      complex(rp) :: product(size(a,1),size(a,2)),value
      integer :: i
      product=matmul(a,matmul(b,matmul(c,d)))
      value=sum([(product(i,i),i=1,size(a,1))])
   end function trace4

   function identity(n) result(a)
      integer, intent(in) :: n
      complex(rp) :: a(n,n)
      integer :: j
      a=cmplx(0.0_rp,0.0_rp,rp)
      do j=1,n; a(j,j)=cmplx(1.0_rp,0.0_rp,rp); end do
   end function identity

   pure function r4_rotate_moment(moment,axis,angle) result(rotated)
      real(rp), intent(in) :: moment(3),axis(3),angle
      real(rp) :: rotated(3)
      real(rp) :: cross_axis_moment(3)
      cross_axis_moment=[axis(2)*moment(3)-axis(3)*moment(2), &
         axis(3)*moment(1)-axis(1)*moment(3), axis(1)*moment(2)-axis(2)*moment(1)]
      rotated=moment*cos(angle)+cross_axis_moment*sin(angle)+ &
         axis*dot_product(axis,moment)*(1.0_rp-cos(angle))
   end function r4_rotate_moment

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
