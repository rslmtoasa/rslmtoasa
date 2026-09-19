! DRESP-06A-FINAL compact covariance regression.
! The q pair is a non-self-inverse one-step displacement on a 4x4x4 grid.
! Every response stage is constructed and checked independently of the
! material driver: bare response, Kxc transport, Dyson, interacting response,
! compact loss, site projection, and projected loss.
program test_dresp06a_commensurate_covariance
   use precision_mod, only: rp
   use basis_mod, only: basis_init
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use lr_response_space_mod, only: response_space_layout
   use lr_projected_site_spin_mod, only: projected_site_spin_contract
   use lr_lmto_product_response_basis_mod, only: lmto_product_response_basis, lmto_product_channel_plus, &
      lmto_product_channel_minus
   use lr_dresp06a_bridge_mod, only: build_compact_covariance_transport, transport_compact_matrix, &
      project_compact_response_to_sites, project_compact_loss_to_sites, compact_loss_matrix, &
      dresp06a_relative_matrix_difference
   use tddft_dyson_mod, only: solve_tddft_dyson_frequency
   implicit none

   integer, parameter :: nsite = 2, npoint = 51, orbital_lmax = 2, response_lmax = 4
   integer :: n, i, j, ik
   real(rp), parameter :: mesh_a = 0.035_rp, mesh_b = 0.11_rp, q = 0.25_rp, tolerance = 2.0e-11_rp
   real(rp) :: radius(npoint), qplus, qminus, mapped_q
   type(lmto_radial_basis), target :: radial(nsite)
   type(response_space_layout), target :: space
   type(projected_site_spin_contract), target :: contract
   type(lmto_product_response_basis), target :: product_plus, product_minus
   complex(rp), allocatable :: transport(:, :), minus_bare(:, :), plus_bare(:, :), minus_kernel(:, :), plus_kernel(:, :)
   complex(rp), allocatable :: minus_d(:, :), plus_d(:, :), minus_chi(:, :), plus_chi(:, :)
   complex(rp), allocatable :: minus_loss(:, :), plus_loss(:, :), mapped(:, :)
   complex(rp), allocatable :: plus_site(:, :), transported_site(:, :), plus_site_loss(:, :), transported_site_loss(:, :)
   real(rp) :: bare_residual, kxc_residual, denominator_residual, interacting_residual
   real(rp) :: compact_loss_residual, site_residual, site_loss_residual, diagonal_residual
   integer :: info

   call basis_init(2)
   call build_mesh(radius)
   call setup_radial_bases(radial, radius)
   call space%initialize(nsite, response_lmax, radius, mesh_a, mesh_b, 1)
   call contract%initialize(space, radial, 'spd')
   call product_plus%initialize(space, radial, lmto_product_channel_plus, .true.)
   call product_minus%initialize(space, radial, lmto_product_channel_minus, .true.)
   n = product_plus%product_dimension
   if (n /= product_minus%product_dimension .or. n < 1) error stop 'covariance fixture: product dimension mismatch'

   qplus = q
   qminus = -q
   do ik = 0, 3
      mapped_q = modulo(real(ik, rp)/4.0_rp + qplus, 1.0_rp)
      if (abs(mapped_q - modulo(real(ik + 1, rp)/4.0_rp, 1.0_rp)) > 1.0e-13_rp) then
         error stop 'covariance fixture: q does not map the finite k mesh'
      end if
   end do
   if (abs(qplus - qminus) < 1.0e-13_rp) error stop 'covariance fixture: q pair is self-inverse'

   allocate(transport(n, n), minus_bare(n, n), plus_bare(n, n), minus_kernel(n, n), plus_kernel(n, n), &
      minus_d(n, n), plus_d(n, n), minus_chi(n, n), plus_chi(n, n), minus_loss(n, n), plus_loss(n, n), &
      mapped(n, n), plus_site(nsite, nsite), transported_site(nsite, nsite), plus_site_loss(nsite, nsite), &
      transported_site_loss(nsite, nsite))
   call build_compact_covariance_transport(product_plus, product_minus, transport)
   do i = 1, n
      do j = 1, n
         minus_bare(i, j) = cmplx(0.0004_rp*real(1 + mod(5*i + 3*j, 31), rp), &
            0.0002_rp*real(i - j, rp), rp)
         minus_kernel(i, j) = cmplx(0.0_rp, 0.0_rp, rp)
      end do
      minus_kernel(i, i) = cmplx(0.03_rp + 0.0001_rp*real(i, rp), 0.0_rp, rp)
      if (i < n) minus_kernel(i, i + 1) = cmplx(0.0003_rp, 0.0_rp, rp)
   end do
   call transport_compact_matrix(product_plus, product_minus, minus_bare, plus_bare)
   call transport_compact_matrix(product_plus, product_minus, minus_kernel, plus_kernel, .false.)
   call transport_compact_matrix(product_plus, product_minus, minus_bare, mapped)
   bare_residual = dresp06a_relative_matrix_difference(plus_bare, mapped)
   call transport_compact_matrix(product_plus, product_minus, minus_kernel, mapped, .false.)
   kxc_residual = dresp06a_relative_matrix_difference(plus_kernel, mapped)

   call solve_tddft_dyson_frequency(minus_bare, minus_kernel, minus_chi, minus_d, info)
   if (info /= 0) error stop 'covariance fixture: minus Dyson solve failed'
   call solve_tddft_dyson_frequency(plus_bare, plus_kernel, plus_chi, plus_d, info)
   if (info /= 0) error stop 'covariance fixture: plus Dyson solve failed'
   call transport_compact_matrix(product_plus, product_minus, minus_d, mapped)
   denominator_residual = maxval(abs(plus_d - mapped))
   call transport_compact_matrix(product_plus, product_minus, minus_chi, mapped)
   interacting_residual = maxval(abs(plus_chi - mapped))
   call compact_loss_matrix(minus_chi, minus_loss)
   call compact_loss_matrix(plus_chi, plus_loss)
   call transport_compact_matrix(product_plus, product_minus, -minus_loss, mapped)
   compact_loss_residual = maxval(abs(plus_loss - mapped))

   call transport_compact_matrix(product_plus, product_minus, minus_chi, mapped)
   call project_compact_response_to_sites(contract, product_plus, plus_chi, plus_site)
   call project_compact_response_to_sites(contract, product_plus, mapped, transported_site)
   site_residual = maxval(abs(plus_site - transported_site))
   call transport_compact_matrix(product_plus, product_minus, -minus_loss, mapped)
   call project_compact_loss_to_sites(contract, product_plus, plus_chi, plus_site_loss)
   call project_compact_response_to_sites(contract, product_plus, mapped, transported_site_loss)
   site_loss_residual = maxval(abs(plus_site_loss - transported_site_loss))
   diagonal_residual = maxval(abs(matmul(conjg(transpose(transport)), transport) - &
      identity_matrix(n)))

   write (*, '(a,es12.4)') '  q-grid mapped covariance residual = ', maxval(abs([qplus - 0.25_rp, qminus + 0.25_rp]))
   write (*, '(a,7(es12.4,1x))') '  covariance stages bare/Kxc/D/interacting/loss/site/site-loss = ', bare_residual, &
      kxc_residual, denominator_residual, interacting_residual, compact_loss_residual, site_residual, site_loss_residual
   if (max(bare_residual, kxc_residual, denominator_residual, interacting_residual, compact_loss_residual, &
       site_residual, site_loss_residual, diagonal_residual) > tolerance) then
      error stop 'UnitDresp06aCommensurateCovariance: FAIL'
   end if
   write (*, '(a)') 'UnitDresp06aCommensurateCovariance: PASS (commensurate q covariance stages)'

contains

   function identity_matrix(size_in) result(identity)
      integer, intent(in) :: size_in
      complex(rp) :: identity(size_in, size_in)
      integer :: index
      identity = cmplx(0.0_rp, 0.0_rp, rp)
      do index = 1, size_in
         identity(index, index) = cmplx(1.0_rp, 0.0_rp, rp)
      end do
   end function identity_matrix

   subroutine build_mesh(mesh)
      real(rp), intent(out) :: mesh(:)
      integer :: ir
      mesh(1) = 0.0_rp
      do ir = 2, size(mesh)
         mesh(ir) = mesh_b*(exp(mesh_a*real(ir - 1, rp)) - 1.0_rp)
      end do
   end subroutine build_mesh

   subroutine setup_radial_bases(bases, mesh)
      type(lmto_radial_basis), intent(out) :: bases(:)
      real(rp), intent(in) :: mesh(:)
      integer :: site, spin, l, ir
      real(rp) :: scale
      do site = 1, size(bases)
         call bases(site)%initialize(size(mesh), orbital_lmax, 2)
         bases(site)%rofi = mesh
         bases(site)%mesh_a = mesh_a
         bases(site)%mesh_b = mesh_b
         bases(site)%channel_present = .true.
         bases(site)%gfac = 1.0_rp
         do spin = 1, 2
            do l = 0, orbital_lmax
               bases(site)%enu_work(l + 1, spin) = -0.28_rp + 0.035_rp*real(l, rp)
               scale = 0.75_rp + 0.025_rp*real(l, rp)
               do ir = 1, size(mesh)
                  bases(site)%phi_large(ir, l + 1, spin) = (0.65_rp + 0.04_rp*real(l, rp))* &
                     mesh(ir)**l*exp(-scale*mesh(ir))
                  bases(site)%phidot_large(ir, l + 1, spin) = (0.10_rp + 0.02_rp*real(l, rp))* &
                     mesh(ir)**l*(1.0_rp + 0.13_rp*mesh(ir))*exp(-0.6_rp*scale*mesh(ir))
                  bases(site)%phi_small(ir, l + 1, spin) = 0.0_rp
                  bases(site)%phidot_small(ir, l + 1, spin) = 0.0_rp
                  bases(site)%phiddot_large(ir, l + 1, spin) = 0.0_rp
                  bases(site)%phiddot_small(ir, l + 1, spin) = 0.0_rp
               end do
            end do
         end do
      end do
   end subroutine setup_radial_bases

end program test_dresp06a_commensurate_covariance
