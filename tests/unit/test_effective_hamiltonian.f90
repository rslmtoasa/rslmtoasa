program test_effective_hamiltonian
   use precision_mod, only: rp
   use basis_mod, only: nb
   use hamiltonian_effective_mod, only: effective_row, resolve_translation
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use control_mod, only: control
   use math_mod, only: inverse_3x3
   use recursion_mod, only: recursion
   implicit none

   integer, parameter :: n = 3
   type(hamiltonian), target :: ham
   type(lattice), target :: lat
   type(control), target :: ctl
   type(recursion) :: rec
   type(effective_row) :: row
   integer, allocatable :: sites(:), keys(:,:,:)
   integer :: i, j, k, d, p, idx(3)
   real(rp) :: cell(3,3), inv(3,3), offset(3), err
   logical :: found
   complex(rp), allocatable :: h(:,:), o(:,:), expected(:,:), actual(:,:), cc(:,:), onsite(:,:), &
                               psi(:,:,:), out(:,:,:), reference(:,:)

   allocate(h(n*nb,n*nb), o(n*nb,n*nb), expected(n*nb,n*nb), actual(n*nb,n*nb), &
            cc(n*nb,n*nb), onsite(n*nb,n*nb), psi(nb,2,n), out(nb,2,n), reference(n*nb,2))
   lat%control => ctl
   ham%lattice => lat
   rec%lattice => lat
   rec%hamiltonian => ham
   lat%kk = n
   lat%pbc = .false.
   lat%nmax = 1
   lat%ntype = 2
   lat%alat = 1.0_rp
   lat%a = 0.0_rp
   allocate(lat%nn(n,3), lat%iz(n), lat%cr(3,n), rec%izero(0:n), rec%idum(0:n), &
            rec%psi2(nb,2,n), rec%hohpsi(nb,2,n))
   lat%iz = [1,2,1]
   lat%nn(1,:) = [2,2,0]
   lat%nn(2,:) = [3,1,3]
   lat%nn(3,:) = [2,2,0]
   lat%cr = 0.0_rp
   lat%cr(1,:) = [0.0_rp,1.0_rp,2.0_rp]
   allocate(ham%ee(nb,nb,3,2), ham%hall(nb,nb,3,1), ham%eecc(nb,nb,3,2), ham%hallcc(nb,nb,3,1), &
            ham%lsham(nb,nb,2), ham%enim(nb,nb,2), ham%obarm(nb,nb,2))
   ham%ee = 0
   ham%hall = 0
   ham%eecc = 0
   ham%hallcc = 0
   ham%lsham = 0
   ham%enim = 0
   ham%obarm = 0
   ham%local_axis = .false.
   ham%ccor_2c = .true.
   ham%lsham(1,2,:) = cmplx(0.0_rp,0.07_rp,rp)
   ham%lsham(2,1,:) = conjg(ham%lsham(1,2,:))
   ham%hall(1,1,1,1) = 0.77_rp
   ham%hall(1,2,2,1) = cmplx(0.17_rp,0.38_rp,rp)
   ham%ee(2,1,2,2) = conjg(ham%hall(1,2,2,1))
   ham%ee(1,2,3,2) = cmplx(0.28_rp,-0.47_rp,rp)
   ham%ee(2,1,2,1) = conjg(ham%ee(1,2,3,2))
   ham%ee(2,2,1,:) = [0.11_rp,0.26_rp]
   ! This entry gives the path 1 -> 2 -> 3 a complex two-hop-only block.
   ham%obarm(2,1,2) = cmplx(0.13_rp,-0.21_rp,rp)
   ham%obarm(1,2,2) = conjg(ham%obarm(2,1,2))
   ham%obarm(1,1,:) = [0.44_rp,0.32_rp]
   ham%obarm(2,2,:) = [0.54_rp,0.22_rp]
   ham%enim(1,1,:) = [0.12_rp,0.51_rp]
   ham%hallcc(1,1,2,1) = 0.03_rp
   ham%eecc(1,1,2,2) = 0.03_rp
   ham%eecc(2,2,1,:) = 0.09_rp

   call ham%effective_graph(.false., sites, keys)
   do p = 1, 2
      ham%hoh = p == 1
      if (.not. ham%hoh) ham%hall(1,1,1,1) = 0.59_rp
      h = 0
      o = 0
      onsite = 0
      cc = 0
      do i = 1, n
         d = lat%iz(i)
         o(ix(i):ix(i)+nb-1,ix(i):ix(i)+nb-1) = ham%obarm(:,:,d)
         onsite(ix(i):ix(i)+nb-1,ix(i):ix(i)+nb-1) = ham%lsham(:,:,d)
         if (ham%hoh) onsite(ix(i):ix(i)+nb-1,ix(i):ix(i)+nb-1) = ham%lsham(:,:,d) + ham%enim(:,:,d)
         do k = 1, lat%nn(i,1)
            j = i
            if (k > 1) j = lat%nn(i,k)
            if (i == 1) then
               h(ix(i):ix(i)+nb-1,ix(j):ix(j)+nb-1) = ham%hall(:,:,k,1)
               cc(ix(i):ix(i)+nb-1,ix(j):ix(j)+nb-1) = ham%hallcc(:,:,k,1)
            else
               h(ix(i):ix(i)+nb-1,ix(j):ix(j)+nb-1) = ham%ee(:,:,k,d)
               cc(ix(i):ix(i)+nb-1,ix(j):ix(j)+nb-1) = ham%eecc(:,:,k,d)
            end if
         end do
      end do
      expected = h + cc + onsite
      if (ham%hoh) expected = expected - matmul(h,matmul(o,h))
      actual = 0
      do i = 1, n
         call ham%assemble_effective_row(sites, keys, i, row)
         do k = 1, row%count
            j = row%key(1,k)
            actual(ix(i):ix(i)+nb-1,ix(j):ix(j)+nb-1) = row%block(:,:,k)
         end do
      end do
      err = maxval(abs(actual - expected))
      print '(a,l1,a,es12.4)', 'HOH=',ham%hoh,' dense matrix max error (Ry): ',err
      if (err > 1.e-13_rp) error stop 'Dense effective Hamiltonian oracle'
      if (ham%hoh .and. abs(actual(1,2*nb+2)) < 1.e-4_rp) error stop 'Missing two-hop-only path'
      rec%izero = 1
      rec%izero(0) = 0
      do i = 1, n
         do j = 1, 2
            do k = 1, nb
               psi(k,j,i) = cmplx(sin(real(k+i*j,rp)),cos(real(i*k-j,rp)),rp)
            end do
            reference(ix(i):ix(i)+nb-1,j) = psi(:,j,i)
         end do
      end do
      reference = matmul(expected,reference)
      if (ham%hoh) then
         call rec%ham_hoh_vec_matmul(psi,out,1.0_rp,0.0_rp)
      else
         call rec%ham_vec_matmul(psi,out,1.0_rp,0.0_rp)
      end if
      do i = 1, n
         do j = 1, 2
            err = max(err,maxval(abs(reference(ix(i):ix(i)+nb-1,j)-out(:,j,i))))
         end do
      end do
      print '(a,es12.4)', 'Native action max error (Ry): ',err
      if (err > 1.e-13_rp) error stop 'Native effective Hamiltonian action'
   end do

   cell = reshape([1.0_rp,0.2_rp,-0.1_rp,0.3_rp,1.2_rp,0.2_rp,0.1_rp,-0.2_rp,1.45_rp],[3,3])
   inv = inverse_3x3(cell)
   offset = matmul(cell, real([-13,7,3],rp))
   call resolve_translation(cell,inv,offset,1.e-9_rp,idx,found)
   if (.not. found .or. any(idx /= [-13,7,3])) error stop 'Skew cell translation'
   call resolve_translation(cell,inv,offset+[0.1_rp,0.0_rp,0.0_rp],1.e-9_rp,idx,found)
   if (found) error stop 'Invalid translation was accepted'
   ! Mixed PBC on a three-site ring. The wrapped direct link and
   ! the open two-hop link end at two distinct images of column 3.
   lat%pbc = .true.
   lat%b1 = .true.
   lat%b2 = .false.
   lat%b3 = .false.
   lat%n1 = 3
   lat%n2 = 1
   lat%n3 = 1
   lat%a = 0.0_rp
   do i = 1, 3
      lat%a(i,i) = 1.0_rp
   end do
   lat%nn(1,:) = [3,2,3]
   lat%nn(3,:) = [3,2,1]
   ham%hoh = .true.
   ham%hall(1,2,3,1) = cmplx(0.12_rp,0.18_rp,rp)
   ham%ee(2,1,3,1) = conjg(ham%hall(1,2,3,1))
   call ham%effective_graph(.false.,sites,keys)
   if (any(keys(:,3,1) /= [3,-1,0,0])) error stop 'Wrong wrapped direct bond key'
   call ham%assemble_effective_row(sites,keys,1,row)
   k = row%lookup([3,0,0,0])
   j = row%lookup([3,-1,0,0])
   if (k == 0 .or. j == 0) error stop 'Two images of the same destination were collapsed'
   if (abs(row%block(1,2,k)) < 1.e-4_rp) error stop 'Open two-hop image missing'
   if (abs(row%block(1,2,j)) < 1.e-4_rp) error stop 'Wrapped direct image missing'
   ! Types along 1 -> 2 -> 3 have weights [1,2,1].
   ! Per-leg weighting makes the derivative of the +2 path twice as large.
   err = abs(cmplx(0.0_rp,4.0_rp,rp)*row%block(1,2,k)-cmplx(0.0_rp,2.0_rp,rp)*row%block(1,2,k))
   if (err < 1.e-4_rp) error stop 'Nonuniform HOH velocity weighting example is degenerate'
   if (any(abs(ham%hall(1,2,2:3,1)-[cmplx(0.17_rp,0.38_rp,rp),cmplx(0.12_rp,0.18_rp,rp)]) > 0)) &
      error stop 'Assembly mutated local hopping blocks'
contains
   integer function ix(site)
      integer, intent(in) :: site
      ix = (site-1)*nb+1
   end function ix
end program test_effective_hamiltonian
