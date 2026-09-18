! DRESP-06A representation gates: complete product space, local ALSDA action,
! exact DRESP-01 site projection, and loss commutation.
program test_dresp06a_projection
   use precision_mod, only: rp
   use basis_mod, only: basis_init
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use lr_response_space_mod, only: response_space_layout
   use lr_projected_site_spin_mod, only: projected_site_spin_contract
   use lr_lmto_product_response_basis_mod, only: lmto_product_response_basis, lmto_product_channel_plus
   use lr_compact_static_interaction_mod, only: compact_project_local_operator, compact_apply_local_operator, &
      compact_project_magnetization, compact_reconstruct_point_vector, compact_weighted_projection_diagnostics
   use lr_dresp06a_bridge_mod, only: project_compact_response_to_sites, project_compact_loss_to_sites, &
      compact_loss_matrix, compact_apply_local_operator_independent, dresp06a_relative_matrix_difference
   implicit none

   integer, parameter :: nsite = 2, npoint = 51, orbital_lmax = 2, response_lmax = 4
   real(rp), parameter :: mesh_a = 0.035_rp, mesh_b = 0.11_rp, tolerance = 3.0e-11_rp
   real(rp) :: radius(npoint), values(nsite, npoint), magnetization(nsite, npoint)
   type(lmto_radial_basis), target :: radial(nsite)
   type(response_space_layout), target :: space
   type(projected_site_spin_contract), target :: contract
   type(lmto_product_response_basis), target :: product
   complex(rp), allocatable :: compact(:, :), projected(:, :), expected(:, :), compact_loss(:, :)
   complex(rp), allocatable :: site_loss_a(:, :), site_loss_b(:, :), vector(:), action_a(:), action_b(:)
   complex(rp), allocatable :: magnetization_compact(:), point(:), reconstructed(:)
   real(rp) :: d_projection, d_loss, d_action, magnetization_norm, magnetization_residual, magnetization_relative
   integer :: i, j, ir, site, mode, flat
   logical :: failed

   call basis_init(2)
   call build_mesh(radius)
   call setup_radial_bases(radial, radius)
   call space%initialize(nsite, response_lmax, radius, mesh_a, mesh_b, 1)
   call contract%initialize(space, radial, 'spd')
   call product%initialize(space, radial, lmto_product_channel_plus, .true.)
   if (product%product_dimension < 1) error stop 'DRESP-06A fixture: complete product basis is empty'

   allocate(compact(product%product_dimension, product%product_dimension), projected(nsite, nsite), &
      expected(nsite, nsite), compact_loss(product%product_dimension, product%product_dimension), &
      site_loss_a(nsite, nsite), site_loss_b(nsite, nsite), vector(product%product_dimension), &
      action_a(product%product_dimension), action_b(product%product_dimension), &
      magnetization_compact(product%product_dimension), point(space%ndim), reconstructed(space%ndim))
   do i = 1, product%product_dimension
      vector(i) = cmplx(0.003_rp*real(1 + mod(7*i, 29), rp), -0.002_rp*real(mod(11*i, 17), rp), rp)
      do j = 1, product%product_dimension
         compact(i, j) = cmplx(2.0e-4_rp*real(1 + mod(i + 3*j, 23), rp), &
            7.0e-5_rp*real(i - j, rp), rp)
      end do
      compact(i, i) = compact(i, i) + cmplx(0.04_rp + 1.0e-5_rp*real(i, rp), -0.002_rp, rp)
   end do
   do site = 1, nsite
      do ir = 1, npoint
         values(site, ir) = 0.03_rp + 0.004_rp*real(site, rp) + 0.001_rp*real(ir, rp)
         magnetization(site, ir) = (0.7_rp + 0.1_rp*real(site, rp))*exp(-0.15_rp*radius(ir))
      end do
      magnetization(site, 1) = 0.0_rp
   end do

   call project_compact_response_to_sites(contract, product, compact, projected)
   call explicit_site_projection(contract, product, compact, expected)
   d_projection = maxval(abs(projected - expected))

   call compact_loss_matrix(compact, compact_loss)
   call project_compact_loss_to_sites(contract, product, compact, site_loss_a)
   call project_compact_response_to_sites(contract, product, compact, site_loss_b)
   call site_loss_matrix(site_loss_b, expected)
   site_loss_b = expected
   d_loss = maxval(abs(site_loss_a - site_loss_b))

   call compact_apply_local_operator(space, product, cmplx(values, 0.0_rp, rp), vector, action_a)
   call compact_apply_local_operator_independent(space, product, values, vector, action_b)
   d_action = maxval(abs(action_a - action_b))

   call compact_project_magnetization(space, product, magnetization, magnetization_compact, point)
   call compact_reconstruct_point_vector(space, product, magnetization_compact, reconstructed)
   call compact_weighted_projection_diagnostics(space, point, reconstructed, magnetization_norm, &
      magnetization_residual, magnetization_relative)

   failed = d_projection > tolerance .or. d_loss > tolerance .or. d_action > tolerance
   failed = failed .or. magnetization_norm <= 0.0_rp .or. magnetization_residual < 0.0_rp
   write (*, '(a,es14.6)') 'DRESP-06A bare compact-to-site max residual = ', d_projection
   write (*, '(a,es14.6)') 'DRESP-06A loss projection commutation residual = ', d_loss
   write (*, '(a,es14.6)') 'DRESP-06A independent Kxc action residual = ', d_action
   write (*, '(a,3(es14.6,1x))') 'DRESP-06A magnetization norm/residual/relative = ', magnetization_norm, &
      magnetization_residual, magnetization_relative
   if (failed) error stop 'UnitDresp06aProjection: FAIL'
   write (*, '(a)') 'UnitDresp06aProjection: PASS (F^H chi F, loss commutation, independent local action, magnetization projection)'

contains

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

   subroutine explicit_site_projection(contract_in, product_in, matrix, result)
      type(projected_site_spin_contract), intent(in) :: contract_in
      type(lmto_product_response_basis), intent(in) :: product_in
      complex(rp), intent(in) :: matrix(:, :)
      complex(rp), intent(out) :: result(:, :)
      complex(rp), allocatable :: f(:, :), image(:)
      integer :: i_site, j_site
      allocate(f(product_in%product_dimension, contract_in%nsite), image(product_in%product_dimension))
      call contract_in%site_integration_functional(product_in, f)
      do j_site = 1, contract_in%nsite
         image = matmul(matrix, f(:, j_site))
         do i_site = 1, contract_in%nsite
            result(i_site, j_site) = dot_product(f(:, i_site), image)
         end do
      end do
      deallocate(f, image)
   end subroutine explicit_site_projection

   subroutine site_loss_matrix(response, loss)
      complex(rp), intent(in) :: response(:, :)
      complex(rp), intent(inout) :: loss(:, :)
      loss = -(response - conjg(transpose(response)))/cmplx(0.0_rp, 2.0_rp*acos(-1.0_rp), rp)
   end subroutine site_loss_matrix

end program test_dresp06a_projection
