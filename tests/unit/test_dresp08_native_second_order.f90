! DRESP-08 focused native second-order LMTO mapping regression.
!
! The fixture is deliberately nontrivial: orbital matrices are noncommuting,
! complex, and spin dependent.  It checks the real-space bond derivative,
! overlap/E_nu tangents, the complete H-h-o-h product rule, and an independent
! finite-rotation central-difference oracle.
program test_dresp08_native_second_order
   use precision_mod, only: rp
   use basis_mod, only: basis_init
   use math_mod, only: init_math_operators, hcpx
   use logger_mod, only: g_logger
   use lmto_magnetic_tangent_mod, only: lmto_bond_value, lmto_bond_derivative, lmto_hhmag_to_spinor
   use lr_dresp08_native_mapping_mod, only: dresp08_build_h2, dresp08_build_product_tangent, &
      dresp08_commutator_tangent, dresp08_build_spinor_tangent, dresp08_rotate_collinear_operator, &
      dresp08_relative_residual
   implicit none

   integer, parameter :: lmax = 2, norb = (lmax + 1)**2, n = 2*norb
   real(rp), parameter :: tol = 2.0e-11_rp
   complex(rp) :: h_up(norb,norb), h_down(norb,norb), o_up(norb,norb), o_down(norb,norb)
   complex(rp) :: e_up(norb,norb), e_down(norb,norb), h(n,n), o(n,n), enu(n,n), h2(n,n)
   complex(rp) :: dh(n,n), do_(n,n), de(n,n), d2h(n,n), comm(n,n)
   complex(rp) :: t_enu(n,n), t_h(n,n), t_left(n,n), t_mid(n,n), t_right(n,n)
   complex(rp) :: hp(n,n), hm(n,n), op(n,n), om(n,n), ep(n,n), em(n,n), h2p(n,n), h2m(n,n)
   complex(rp) :: generator(n,n), hhh(norb,norb), bond(n,n), dbond(n,n)
   complex(rp) :: hhmag(norb,norb,4), dhmag(norb,norb,4), delta_o(n,n), delta_e(n,n)
   complex(rp) :: coeff_o(norb), coeff_e(norb)
   real(rp) :: theta, err, previous, slope, max_bond_error, max_overlap_error, max_enu_error
   real(rp) :: max_product_error, angular_offdiag, angular_anisotropy
   integer :: i, j, idir, shell_start, shell_end
   logical :: failed

   call basis_init(lmax)
   call init_math_operators()
   call g_logger%init()

   call build_nontrivial_hermitian(h_up, 0.17_rp, 0.021_rp)
   call build_nontrivial_hermitian(h_down, -0.11_rp, -0.017_rp)
   call build_nontrivial_hermitian(o_up, 0.31_rp, 0.013_rp)
   call build_nontrivial_hermitian(o_down, 0.27_rp, -0.009_rp)
   call build_nontrivial_hermitian(e_up, -0.42_rp, 0.008_rp)
   call build_nontrivial_hermitian(e_down, -0.36_rp, -0.006_rp)
   call embed(h_up, h_down, h)
   call embed(o_up, o_down, o)
   call embed(e_up, e_down, enu)

   call make_generator(generator)
   call dresp08_build_h2(h, o, enu, h2)
   call dresp08_commutator_tangent(generator, h, dh)
   call dresp08_commutator_tangent(generator, o, do_)
   call dresp08_commutator_tangent(generator, enu, de)
   call dresp08_build_product_tangent(de, dh, do_, h, o, d2h, t_enu, t_h, t_left, t_mid, t_right)
   call dresp08_commutator_tangent(generator, h2, comm)
   max_product_error = dresp08_relative_residual(d2h, comm)

   coeff_o = cmplx(0.0_rp, 0.0_rp, rp)
   coeff_e = cmplx(0.0_rp, 0.0_rp, rp)
   do i = 1, norb
      coeff_o(i) = cmplx(0.041_rp + 0.007_rp*real(i, rp), 0.0_rp, rp)
      coeff_e(i) = cmplx(-0.083_rp + 0.011_rp*real(i, rp), 0.0_rp, rp)
   end do
   call dresp08_build_spinor_tangent(coeff_o, delta_o)
   call dresp08_build_spinor_tangent(coeff_e, delta_e)
   call dresp08_commutator_tangent(generator, make_collinear(coeff_o, 0.19_rp), comm)
   max_overlap_error = dresp08_relative_residual(delta_o, comm)
   call dresp08_commutator_tangent(generator, make_collinear(coeff_e, -0.23_rp), comm)
   max_enu_error = dresp08_relative_residual(delta_e, comm)

   hhh = cmplx(0.0_rp, 0.0_rp, rp)
   do i = 1, norb
      hhh(i,i) = cmplx(0.12_rp + 0.015_rp*real(i, rp), 0.0_rp, rp)
      do j = i + 1, norb
         hhh(i,j) = cmplx(0.008_rp*real(i+j, rp), 0.002_rp*real(j-i, rp), rp)
         hhh(j,i) = cmplx(-0.003_rp*real(i+j, rp), 0.001_rp*real(j-i, rp), rp)
      end do
   end do
   call lmto_bond_value(hhh, orbital_vector(0.93_rp), orbital_vector(0.17_rp), orbital_vector(1.07_rp), &
      orbital_vector(-0.13_rp), orbital_vector(0.04_rp), orbital_vector(-0.02_rp), [0.0_rp,0.0_rp,1.0_rp], &
      [0.0_rp,0.0_rp,1.0_rp], .false., hhmag)
   do idir = 1, 4
      call hcpx(hhmag(:,:,idir), 'cart2sph')
   end do
   call lmto_hhmag_to_spinor(hhmag, bond)
   call lmto_bond_derivative(hhh, orbital_vector(0.93_rp), orbital_vector(0.17_rp), orbital_vector(1.07_rp), &
      orbital_vector(-0.13_rp), orbital_vector(0.04_rp), [0.0_rp,0.0_rp,1.0_rp], [0.0_rp,0.0_rp,1.0_rp], &
      [1.0_rp,0.0_rp,0.0_rp], [1.0_rp,0.0_rp,0.0_rp], .false., dhmag)
   do idir = 1, 4
      call hcpx(dhmag(:,:,idir), 'cart2sph')
   end do
   call lmto_hhmag_to_spinor(dhmag, dbond)
   call dresp08_commutator_tangent(generator, bond, comm)
   max_bond_error = dresp08_relative_residual(dbond, comm)

   ! L=0 selection fixture: equal coefficients within each shell remain
   ! diagonal and m-degenerate after the cartesian/spherical transform.
   coeff_o = cmplx(0.0_rp, 0.0_rp, rp)
   do shell_start = 0, lmax
      shell_end = shell_start*shell_start + 2*shell_start + 1
      coeff_o(shell_start*shell_start+1:shell_end) = cmplx(0.2_rp + 0.1_rp*real(shell_start, rp), 0.0_rp, rp)
   end do
   call dresp08_build_spinor_tangent(coeff_o, delta_o)
   angular_offdiag = sqrt(sum(abs(delta_o(1:norb,1:norb) - diag_spin_part(delta_o(1:norb,1:norb)))**2) + &
      sum(abs(delta_o(norb+1:2*norb,norb+1:2*norb) - diag_spin_part(delta_o(norb+1:2*norb,norb+1:2*norb)))**2) + &
      sum(abs(delta_o(1:norb,norb+1:2*norb) - diag_spin_part(delta_o(1:norb,norb+1:2*norb)))**2) + &
      sum(abs(delta_o(norb+1:2*norb,1:norb) - diag_spin_part(delta_o(norb+1:2*norb,1:norb)))**2))
   angular_anisotropy = shell_anisotropy(delta_o)

   previous = huge(1.0_rp)
   slope = 0.0_rp
   do i = 1, 5
      theta = 0.20_rp/(2.0_rp**real(i-1, rp))
      call dresp08_rotate_collinear_operator(h_up, h_down, theta, hp)
      call dresp08_rotate_collinear_operator(h_up, h_down, -theta, hm)
      call dresp08_rotate_collinear_operator(o_up, o_down, theta, op)
      call dresp08_rotate_collinear_operator(o_up, o_down, -theta, om)
      call dresp08_rotate_collinear_operator(e_up, e_down, theta, ep)
      call dresp08_rotate_collinear_operator(e_up, e_down, -theta, em)
      call dresp08_build_h2(hp, op, ep, h2p)
      call dresp08_build_h2(hm, om, em, h2m)
      err = dresp08_relative_residual((h2p-h2m)/(2.0_rp*theta), d2h)
      if (i > 1) slope = max(slope, err/previous)
      previous = err
   end do

   failed = max_product_error > tol .or. max_bond_error > tol .or. max_overlap_error > tol .or. max_enu_error > tol .or. &
      angular_offdiag > tol .or. angular_anisotropy > tol .or. slope > 0.30_rp
   if (failed) then
      write (*, '(a)') 'UnitDresp08NativeSecondOrder: FAIL'
      write (*, '(a,es12.4)') '  bond_error=', max_bond_error
      write (*, '(a,es12.4)') '  overlap_error=', max_overlap_error
      write (*, '(a,es12.4)') '  enu_error=', max_enu_error
      write (*, '(a,es12.4)') '  product_error=', max_product_error
      write (*, '(a,es12.4)') '  angular_offdiag=', angular_offdiag
      write (*, '(a,es12.4)') '  angular_anisotropy=', angular_anisotropy
      write (*, '(a,es12.4)') '  central_difference_ratio=', slope
      error stop 1
   end if
   write (*, '(a,es12.4)') '  lmto_bond_derivative_vs_commutator=', max_bond_error
   write (*, '(a,es12.4)') '  overlap_tangent_vs_commutator=', max_overlap_error
   write (*, '(a,es12.4)') '  enu_tangent_vs_commutator=', max_enu_error
   write (*, '(a,es12.4)') '  product_rule_vs_commutator=', max_product_error
   write (*, '(a,es12.4)') '  finite_difference_final=', previous
   write (*, '(a,es12.4)') '  finite_difference_ratio_max=', slope
   write (*, '(a)') 'UnitDresp08NativeSecondOrder: PASS (native bond, overlap, E_nu, HOH product, finite rotation, L=0)'

contains

   subroutine build_nontrivial_hermitian(a, diagonal, scale)
      complex(rp), intent(out) :: a(:,:)
      real(rp), intent(in) :: diagonal, scale
      integer :: i, j
      a = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, size(a,1)
         a(i,i) = cmplx(diagonal + 0.03_rp*real(i, rp), 0.0_rp, rp)
         do j = i + 1, size(a,1)
            a(i,j) = cmplx(scale*real(i+j, rp), 0.004_rp*real(j-i, rp), rp)
            a(j,i) = conjg(a(i,j))
         end do
      end do
   end subroutine build_nontrivial_hermitian

   subroutine embed(up, down, full)
      complex(rp), intent(in) :: up(:,:), down(:,:)
      complex(rp), intent(out) :: full(:,:)
      integer :: nloc
      nloc = size(up,1)
      full = cmplx(0.0_rp, 0.0_rp, rp)
      full(1:nloc,1:nloc) = up
      full(nloc+1:2*nloc,nloc+1:2*nloc) = down
   end subroutine embed

   subroutine make_generator(g)
      complex(rp), intent(out) :: g(:,:)
      integer :: i, nloc
      nloc = size(g,1)/2
      g = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, nloc
         g(i,nloc+i) = cmplx(0.0_rp,-0.5_rp,rp)
         g(nloc+i,i) = cmplx(0.0_rp,0.5_rp,rp)
      end do
   end subroutine make_generator

   function make_collinear(vector, scalar) result(full)
      complex(rp), intent(in) :: vector(:)
      real(rp), intent(in) :: scalar
      complex(rp) :: full(2*norb,2*norb)
      integer :: i
      full = cmplx(0.0_rp,0.0_rp,rp)
      do i = 1, norb
         full(i,i) = cmplx(scalar,0.0_rp,rp) + vector(i)
         full(norb+i,norb+i) = cmplx(scalar,0.0_rp,rp) - vector(i)
      end do
      call hcpx(full(1:norb,1:norb), 'cart2sph')
      call hcpx(full(norb+1:2*norb,norb+1:2*norb), 'cart2sph')
      call hcpx(full(1:norb,norb+1:2*norb), 'cart2sph')
      call hcpx(full(norb+1:2*norb,1:norb), 'cart2sph')
   end function make_collinear

   function orbital_vector(base) result(v)
      real(rp), intent(in) :: base
      complex(rp) :: v(norb)
      integer :: i
      do i = 1, norb
         v(i) = cmplx(base + 0.009_rp*real(i, rp), 0.001_rp*real(mod(i,3),rp), rp)
      end do
   end function orbital_vector

   function diag_spin_part(a) result(d)
      complex(rp), intent(in) :: a(:,:)
      complex(rp) :: d(size(a,1),size(a,2))
      integer :: i
      d = cmplx(0.0_rp,0.0_rp,rp)
      do i = 1, size(a,1)
         d(i,i) = a(i,i)
      end do
   end function diag_spin_part

   function shell_anisotropy(a) result(value)
      complex(rp), intent(in) :: a(:,:)
      real(rp) :: value
      integer :: l, first, last, i
      value = 0.0_rp
      do l = 0, lmax
         first = l*l + 1
         last = (l+1)*(l+1)
         do i = first + 1, last
            value = max(value, abs(a(i,norb+i)-a(first,norb+first)))
         end do
      end do
   end function shell_anisotropy

end program test_dresp08_native_second_order
