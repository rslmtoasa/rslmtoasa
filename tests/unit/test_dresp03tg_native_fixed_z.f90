!------------------------------------------------------------------------------
! DRESP-03TG-C -- first hard gate: fixed-complex-energy native contraction.
!
! The native and finite-H objects are assembled independently from the live
! potential parameters and a production structure-constant block.  The test
! compares the complete exchange integrand, not a raw Green-function block.
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
      d_screen(:,:), eye(:,:), tmp(:,:), g_gamma(:,:), g_alpha(:,:), h(:,:), gh(:,:)
   complex(rp), allocatable :: d1(:,:), d2(:,:), gup(:,:), gdown(:,:), ghup(:,:), ghdown(:,:)
   complex(rp), allocatable :: path_p(:,:), path_s(:,:), path_g(:,:)
   real(rp), allocatable :: alpha(:)
   complex(rp) :: z_values(3), z
   real(rp) :: p_error, d_error, solve_error, alpha_error, fixed_error, pauli_error
   real(rp) :: native_alpha_value, native_gamma_value, finite_h_value, gamma_h_error, alpha_h_error, gf_cov_error
   real(rp), allocatable :: s0(:,:), s1(:,:)
   real(rp) :: e
   integer :: i, l, m, lm, lmax, norb, nsite, n2, n4, ia, ino, ncopy
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
   n2 = 2*norb
   n4 = 4*norb
   allocate(p_real(n2,n2,1),p_complex(n2,n2),p_screen(n2,n2),delta(norb,norb),expected(norb,norb),alpha(0:lmax))
   allocate(p_gamma(n4,n4),p_alpha(n4,n4),s_alpha(n4,n4),s_gamma(n4,n4),d_screen(n4,n4),eye(n4,n4),tmp(n4,n4), &
      g_gamma(n4,n4),g_alpha(n4,n4),h(n4,n4),gh(n4,n4),d1(norb,norb),d2(norb,norb),gup(norb,norb),gdown(norb,norb), &
      ghup(norb,norb),ghdown(norb,norb))
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
   s_alpha=cmplx(0.0_rp,0.0_rp,rp)
   s_alpha(1:norb,1:norb)=s0
   s_alpha(1:norb,norb+1:2*norb)=s1
   s_alpha(norb+1:2*norb,1:norb)=transpose(s1)
   s_alpha(norb+1:2*norb,norb+1:2*norb)=s0
   ! Lift orbital S to spin-major order: (site1,site2)_up,(site1,site2)_down.
   tmp=cmplx(0.0_rp,0.0_rp,rp)
   tmp(1:2*norb,1:2*norb)=s_alpha
   tmp(2*norb+1:n4,2*norb+1:n4)=s_alpha
   s_alpha=tmp

   p_gamma=cmplx(0.0_rp,0.0_rp,rp); p_alpha=cmplx(0.0_rp,0.0_rp,rp); d_screen=cmplx(0.0_rp,0.0_rp,rp)
   do i=1,2
      call native_complex_p_matrix(lat%symbolic_atoms(1),cmplx(0.0_rp,0.0_rp,rp),p_complex)
      p_gamma((i-1)*norb+1:i*norb,(i-1)*norb+1:i*norb)=p_complex(1:norb,1:norb)
      p_gamma(2*norb+(i-1)*norb+1:2*norb+i*norb,2*norb+(i-1)*norb+1:2*norb+i*norb)= &
         p_complex(norb+1:n2,norb+1:n2)
      call native_screened_p_matrix(lat%symbolic_atoms(1),cmplx(0.0_rp,0.0_rp,rp),alpha,p_screen)
      p_alpha((i-1)*norb+1:i*norb,(i-1)*norb+1:i*norb)=p_screen(1:norb,1:norb)
      p_alpha(2*norb+(i-1)*norb+1:2*norb+i*norb,2*norb+(i-1)*norb+1:2*norb+i*norb)= &
         p_screen(norb+1:n2,norb+1:n2)
      do l=0,lmax
         do m=1,2*l+1
            lm=l*l+m
            d_screen((i-1)*norb+lm,(i-1)*norb+lm)=alpha(l)-lat%symbolic_atoms(1)%potential%qpar(l,1)
            d_screen(2*norb+(i-1)*norb+lm,2*norb+(i-1)*norb+lm)=alpha(l)-lat%symbolic_atoms(1)%potential%qpar(l,2)
         end do
      end do
   end do
   eye=identity(n4)
   call native_inverse(eye+matmul(s_alpha,d_screen),tmp)
   s_gamma=matmul(tmp,s_alpha)

   failed=.false.; fixed_error=0.0_rp; alpha_error=0.0_rp; gamma_h_error=0.0_rp; alpha_h_error=0.0_rp; gf_cov_error=0.0_rp; pauli_error=0.0_rp
   z_values=[cmplx(-0.91_rp,0.83_rp,rp),cmplx(-0.17_rp,0.04_rp,rp),cmplx(0.62_rp,0.31_rp,rp)]
   do i=1,size(z_values)
      z=z_values(i)
      call native_complex_p_matrix(lat%symbolic_atoms(1),z,p_complex)
      call native_screened_p_matrix(lat%symbolic_atoms(1),z,alpha,p_screen)
      p_gamma=cmplx(0.0_rp,0.0_rp,rp); p_alpha=cmplx(0.0_rp,0.0_rp,rp)
      ! Repeat the two identical site blocks at the fixed complex energy in
      ! global spin-major order: up(site1,site2), down(site1,site2).
      do lm=1,norb
         p_gamma(lm,lm)=p_complex(lm,lm)
         p_gamma(norb+lm,norb+lm)=p_complex(lm,lm)
         p_gamma(2*norb+lm,2*norb+lm)=p_complex(norb+lm,norb+lm)
         p_gamma(3*norb+lm,3*norb+lm)=p_complex(norb+lm,norb+lm)
         p_alpha(lm,lm)=p_screen(lm,lm)
         p_alpha(norb+lm,norb+lm)=p_screen(lm,lm)
         p_alpha(2*norb+lm,2*norb+lm)=p_screen(norb+lm,norb+lm)
         p_alpha(3*norb+lm,3*norb+lm)=p_screen(norb+lm,norb+lm)
      end do
      call native_inverse(p_gamma-s_gamma,g_gamma)
      call native_inverse(p_alpha-s_alpha,g_alpha)
      gup=g_gamma(1:norb,norb+1:2*norb)
      gdown=g_gamma(3*norb+1:4*norb,2*norb+1:3*norb)
      gh=cmplx(0.0_rp,0.0_rp,rp)
      do l=0,lmax
         do m=1,2*l+1
            lm=l*l+m
            gh(lm,lm)=cmplx(lat%symbolic_atoms(1)%potential%c(l,1)+lat%symbolic_atoms(1)%potential%vmad,0.0_rp,rp)
            gh(norb+lm,norb+lm)=gh(lm,lm)
            gh(2*norb+lm,2*norb+lm)=cmplx(lat%symbolic_atoms(1)%potential%c(l,2)+lat%symbolic_atoms(1)%potential%vmad,0.0_rp,rp)
            gh(3*norb+lm,3*norb+lm)=gh(2*norb+lm,2*norb+lm)
            d_screen(lm,lm)=cmplx(lat%symbolic_atoms(1)%potential%dele(l,1),0.0_rp,rp)
            d_screen(norb+lm,norb+lm)=d_screen(lm,lm)
            d_screen(2*norb+lm,2*norb+lm)=cmplx(lat%symbolic_atoms(1)%potential%dele(l,2),0.0_rp,rp)
            d_screen(3*norb+lm,3*norb+lm)=d_screen(2*norb+lm,2*norb+lm)
         end do
      end do
      h=gh+matmul(d_screen,matmul(s_gamma,d_screen))
      h=z*identity(n4)-h
      call native_inverse(h,gh)
      d1=cmplx(0.0_rp,0.0_rp,rp); d2=d1
      do l=0,lmax
         do m=1,2*l+1
            lm=l*l+m
            d1(lm,lm)=d_screen(lm,lm)*d_screen(2*norb+lm,2*norb+lm)*(p_gamma(lm,lm)-p_gamma(2*norb+lm,2*norb+lm))
            d2(lm,lm)=d1(lm,lm)
         end do
      end do
      ghup=gh(1:norb,norb+1:2*norb)
      ghdown=gh(3*norb+1:4*norb,2*norb+1:3*norb)
      delta=p_alpha(1:norb,1:norb)-p_alpha(2*norb+1:3*norb,2*norb+1:3*norb)
      expected=p_alpha(norb+1:2*norb,norb+1:2*norb)-p_alpha(3*norb+1:4*norb,3*norb+1:4*norb)
      native_alpha_value=native_exchange_integrand(delta,expected,g_alpha(1:norb,norb+1:2*norb), &
         g_alpha(3*norb+1:4*norb,2*norb+1:3*norb))
      delta=p_gamma(1:norb,1:norb)-p_gamma(2*norb+1:3*norb,2*norb+1:3*norb)
      expected=p_gamma(norb+1:2*norb,norb+1:2*norb)-p_gamma(3*norb+1:4*norb,3*norb+1:4*norb)
      native_gamma_value=native_exchange_integrand(delta,expected,g_gamma(1:norb,norb+1:2*norb), &
         g_gamma(3*norb+1:4*norb,2*norb+1:3*norb))
      finite_h_value=native_finite_h_integrand(d1,d2,ghup,ghdown)
      fixed_error=abs(native_alpha_value-finite_h_value)
      gamma_h_error=max(gamma_h_error,abs(native_gamma_value-finite_h_value))
      alpha_error=max(alpha_error,abs(native_alpha_value-native_gamma_value))
      pauli_error=max(pauli_error,abs(native_collinear_pauli_integrand(d1,d2, &
         gamma_spin_block(gh,norb,1,2),gamma_spin_block(gh,norb,2,1))-finite_h_value))
      write(*,'(a,2es14.6,a,3es14.6)') 'z=',real(z,rp),aimag(z),' alpha-H/gamma-H/alpha-gamma=', &
         fixed_error,abs(native_gamma_value-finite_h_value),abs(native_alpha_value-native_gamma_value)
   end do
   solve_error=maxval(abs(matmul(p_alpha-s_alpha,g_alpha)-identity(n4)))
   write(*,'(a,es14.6)') 'DRESP-03TG native inverse residual = ',solve_error
   write(*,'(a,es14.6)') 'DRESP-03TG screening-invariance residual = ',alpha_error
   write(*,'(a,es14.6)') 'DRESP-03TG gamma finite-H residual = ',gamma_h_error
   failed = fixed_error > 2.0e-10_rp .or. alpha_error > 2.0e-10_rp

   write(*,'(a,es14.6)') 'DRESP-03TG historical Pauli/native residual = ',pauli_error
   if (p_error > 2.0e-14_rp .or. d_error > 2.0e-13_rp .or. solve_error > 2.0e-12_rp .or. failed .or. &
       gamma_h_error > 2.0e-10_rp .or. pauli_error > 2.0e-10_rp) then
      write(*,'(a)') 'RESULT: BLOCKED — FIXED-Z REPRESENTATION CONTRACTION'
      return
   end if
   write(*,'(a)') 'DRESP-03TG FIXED-Z GATE: PASS'

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
      complex(rp) :: block(n,n)
      integer :: a0,b0
      a0=(i_site-1)*n; b0=(j_site-1)*n
      block=cmplx(0.0_rp,0.0_rp,rp)
      block(1:n,1:n)=g(a0+1:a0+n,b0+1:b0+n)
      block(n+1:2*n,n+1:2*n)=g(2*n+a0+1:2*n+a0+n,2*n+b0+1:2*n+b0+n)
   end function gamma_spin_block

end program test_dresp03tg_native_fixed_z
