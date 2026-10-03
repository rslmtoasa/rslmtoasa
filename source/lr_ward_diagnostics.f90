! Read-only LR-METHOD-05R1 material export. Invoked only by an explicit
! RSLMTO_LR_WARD_DUMP environment variable, never used by the interaction.
module lr_ward_diagnostics_mod
   use precision_mod, only: rp
   use linear_response_mod
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, lmto_orbital_l
   use radial_ground_state_mod, only: radial_ground_state
   use reciprocal_mod, only: reciprocal
   use lattice_mod, only: lattice
   use hamiltonian_mod, only: hamiltonian
   use pauli_ground_state_projection_mod, only: compute_accepted_pauli_magnetization
   implicit none
   private
   public :: dump_lr_ward_inputs
contains
   subroutine dump_lr_ward_inputs(path, space, radial, states, electronic, product, recip, lat, ham)
      character(len=*), intent(in) :: path
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial(:)
      type(radial_ground_state), intent(in) :: states(:)
      type(lr_electronic_state), intent(in) :: electronic
      type(lmto_product_response_basis), intent(in) :: product
      type(reciprocal), intent(in) :: recip
      type(lattice), intent(in) :: lat
      type(hamiltonian), intent(in) :: ham
      type(lmto_live_hamiltonian_fixture) :: fixture, without_hoh
      type(lmto_radial_basis) :: pauli_radial
      type(pauli_endpoint_state) :: left, right
      real(rp), allocatable :: total(:,:), valence(:,:), core(:,:), bxc(:,:), branches(:,:,:,:)
      complex(rp), allocatable :: h(:,:,:), torque(:,:,:), torque_base(:,:,:), hf(:,:), t(:,:), tf(:,:), tc(:), projected(:,:), modes(:,:)
      integer, allocatable :: indices(:,:)
      integer :: nk, nb, nmat, norb, nr, nc, npair, ik, ib, jb, l, branch, ir, unit, site, m, mode, i
      real(rp) :: adapter_error, up_norm, down_norm
      if (space%nsite /= 1 .or. product%circular_channel /= lmto_product_channel_plus) &
         error stop '05R1 export requires the accepted one-site chi_plus Fe gate'
      nk=electronic%nk; nb=electronic%nbands; nmat=electronic%nbasis; norb=nmat/2
      nr=space%npoint; nc=product%product_dimension
      call compute_accepted_pauli_magnetization(recip, lat%symbolic_atoms, lat%nbulk, total, valence, core)
      allocate(bxc(1,nr), branches(nr,states(1)%pauli_lmax+1,6,2), modes(nr,nc), indices(4,nc), projected(nc,4))
      bxc(1,:)=0.5_rp*(states(1)%vxc_up-states(1)%vxc_down)
      call compact_project_magnetization(space,product,bxc,projected(:,1))
      call compact_project_magnetization(space,product,total,projected(:,2))
      call compact_project_magnetization(space,product,valence,projected(:,3))
      call compact_project_magnetization(space,product,core,projected(:,4))
      ! Reproduce the existing radial_points adapter exactly: large phi and
      ! Pauli phidot; its initialized phiddot/small stay zero and GFAC stays one.
      call pauli_radial%initialize(nr,states(1)%pauli_lmax,2)
      pauli_radial%rofi=states(1)%r; pauli_radial%mesh_a=states(1)%a; pauli_radial%mesh_b=states(1)%b
      pauli_radial%phi_large=states(1)%pauli_large
      pauli_radial%phidot_large=states(1)%pauli_large_dot
      pauli_radial%enu_work=states(1)%pauli_enu; pauli_radial%enu_radial=states(1)%pauli_enu
      pauli_radial%channel_present=.true.
      do branch=1,6
         do l=0,states(1)%pauli_lmax
            do ir=1,nr
               call lmto_product_second_order_radial_branch(radial(1),ir,l,l,1,2,branch,branches(ir,l+1,branch,1))
               call lmto_product_second_order_radial_branch(pauli_radial,ir,l,l,1,2,branch,branches(ir,l+1,branch,2))
            end do
         end do
      end do
      do i=1,nc
         call product%unflatten_index(i,site,l,m,mode)
         indices(:,i)=[site,l,m,mode]
         modes(:,i)=product%blocks(site,l)%weighted_modes(:,mode)
      end do
      call lmto_fixture_from_hamiltonian(ham,fixture)
      without_hoh=fixture; without_hoh%hoh=.false.
      allocate(h(nmat,nmat,nk), torque(nmat,nmat,nk), torque_base(nmat,nmat,nk), &
         hf(nmat,nmat),t(nmat,nmat),tf(nmat,nmat),tc(nc))
      adapter_error=0.0_rp
      do ik=1,nk
         h(:,:,ik)=recip%hk_bulk(:,:,ik)
         call assemble_lmto_hamiltonian(fixture,electronic%k_points(:,ik),hf)
         adapter_error=max(adapter_error,maxval(abs(hf-h(:,:,ik))))
         call assemble_lmto_finite_q_torque(fixture,electronic%k_points(:,ik),[0.0_rp,0.0_rp,0.0_rp], &
            1,[0.0_rp,1.0_rp,0.0_rp],t)
         call assemble_lmto_finite_q_torque(without_hoh,electronic%k_points(:,ik),[0.0_rp,0.0_rp,0.0_rp], &
            1,[0.0_rp,1.0_rp,0.0_rp],tf)
         torque(:,:,ik)=t; torque_base(:,:,ik)=tf
      end do
      npair=0
      do ik=1,nk
         do ib=1,nb
            up_norm=sum(abs(electronic%eigenvectors(:norb,ib,ik))**2)
            if (up_norm < 1.0e-12_rp) cycle
            do jb=1,nb
               down_norm=sum(abs(electronic%eigenvectors(norb+1:,jb,ik))**2)
               if (down_norm < 1.0e-12_rp) cycle
               npair=npair+1
            end do
         end do
      end do
      open(newunit=unit,file=path,status='replace',access='stream',form='unformatted',action='write')
      write(unit) [501,1,nk,nb,nmat,nr,nc,states(1)%pauli_lmax,npair]
      write(unit) electronic%fermi_level,electronic%temperature,states(1)%a,states(1)%b, &
         lat%symbolic_atoms(lat%nbulk+1)%potential%mtot,adapter_error
      write(unit) electronic%eigenvalues,electronic%eigenvectors,electronic%occupations, &
         electronic%k_weights,electronic%k_points,h,torque,torque_base
      write(unit) space%radius,space%radial_weights,states(1)%n_up-states(1)%n_down, &
         bxc,total,valence,core,projected,branches,indices,modes
      do ik=1,nk
         do ib=1,nb
            if (sum(abs(electronic%eigenvectors(:norb,ib,ik))**2) < 1.0e-12_rp) cycle
            call left%initialize(electronic%eigenvalues(ib,ik),electronic%eigenvectors(:,ib,ik))
            do jb=1,nb
               if (sum(abs(electronic%eigenvectors(norb+1:,jb,ik))**2) < 1.0e-12_rp) cycle
               call right%initialize(electronic%eigenvalues(jb,ik),electronic%eigenvectors(:,jb,ik))
               call product%transition_coordinates(left,right,tc)
               write(unit) [ik,ib,jb],tc
            end do
         end do
      end do
      close(unit)
   end subroutine dump_lr_ward_inputs
end module lr_ward_diagnostics_mod
