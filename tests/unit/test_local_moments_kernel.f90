program test_local_moments_kernel
   use precision_mod, only: rp
   use math_mod, only: correct_d_shell_moments, extract_hubbard_d_shell_potential, fixed_impurity_atom_carries_uj
   implicit none
   integer, parameter :: ncase = 4
   character(len=1024) :: path, magic
   integer :: unit, ios, nc, k, site, spin, i, j, status
   real(rp) :: q0, c, q2, c_direct, q2_direct, c_bound, q2_bound
   real(rp) :: cp, q2p, t, wrong_cp, wrong_q2, scalar
   real(rp) :: scalar_cp, scalar_q2, vout_cp, vout_q2, conj_cp, conj_q2
   complex(rp) :: p(5,5), m1(5,5), vin(5,5), vout(5,5), mapped(18,18,2), got(5,5)
   complex(rp) :: m1_scalar(5,5), vscalar(5,5)
   real(rp) :: centers(ncase), q2s(ncase), cbounds(ncase), qbounds(ncase)
   integer :: sites(ncase), spins(ncase)
   complex(rp) :: ps(5,5,ncase), m1s(5,5,ncase), vins(5,5,ncase)
   real(rp) :: q0s(ncase), cs(ncase), rawq2s(ncase)
   logical :: failed

   failed = .false.
   call get_command_argument(1, path)
   if (len_trim(path) == 0) error stop 'fixture path required'
   open(newunit=unit, file=trim(path), status='old', action='read', iostat=ios)
   if (ios /= 0) error stop 'cannot open fixture'
   read(unit, *, iostat=ios) magic, nc
   if (ios /= 0 .or. trim(magic) /= 'A3M14' .or. nc /= 1) error stop 'bad fixture version'
   read(unit, *, iostat=ios) nc
   if (ios /= 0 .or. nc /= ncase) error stop 'bad fixture case count'
   do k = 1, ncase
      read(unit, *, iostat=ios) site, spin, q0, c, q2, c_direct, q2_direct, c_bound, q2_bound
      if (ios /= 0) error stop 'bad fixture case header'
      call read_matrix(unit, p)
      call read_matrix(unit, m1)
      call read_matrix(unit, vin)
      call read_matrix(unit, vout)
      call correct_d_shell_moments(p, m1, vin, q0, c, q2, cp, q2p, t, status)
      if (status /= 0) then
         write(*,*) 'FAIL valid fixture status', k, status
         failed = .true.
      else
         call check(abs(cp-c_direct) <= c_bound, 'direct c prime', k)
         call check(abs(q2p-q2_direct) <= q2_bound, 'direct q2 prime', k)
      end if
      sites(k)=site; spins(k)=spin; centers(k)=cp; q2s(k)=q2p
      cbounds(k)=c_bound; qbounds(k)=q2_bound
      ps(:,:,k)=p; m1s(:,:,k)=m1; vins(:,:,k)=vin
      q0s(k)=q0; cs(k)=c; rawq2s(k)=q2

      ! A scalar operator has the preregistered exact action c'=c-v and q2'=q2.
      scalar = 0.137_rp
      vscalar = cmplx(0.0_rp,0.0_rp,rp)
      do i=1,5
         vscalar(i,i)=cmplx(scalar,0.0_rp,rp)
      end do
      call correct_d_shell_moments(p, m1, vscalar, q0, c, q2, wrong_cp, wrong_q2, t, status)
      call check(status == 0, 'scalar-shift status', k)
      call check(abs(wrong_cp-(c-scalar)) <= c_bound, 'scalar-operator center identity', k)
      call check(abs(wrong_q2-q2) <= q2_bound, 'scalar-operator q2 identity', k)
      scalar_cp=wrong_cp; scalar_q2=wrong_q2
      ! Wrong operator provenance controls must remain distinguishable using
      ! the case-specific preregistered bounds.
      call correct_d_shell_moments(p, m1, vout, q0, c, q2, wrong_cp, wrong_q2, t, status)
      call check(status == 0, 'Vout control status', k)
      call check(abs(wrong_cp-cp) > 2.0_rp*c_bound .or. abs(wrong_q2-q2p) > 2.0_rp*q2_bound, &
         'Vout control separation', k)
      vout_cp=wrong_cp; vout_q2=wrong_q2
      call correct_d_shell_moments(conjg(p), conjg(m1), vin, q0, c, q2, wrong_cp, wrong_q2, t, status)
      call check(status == 0, 'conjugated-eigenvector control status', k)
      call check(abs(wrong_cp-cp) > 2.0_rp*c_bound .or. abs(wrong_q2-q2p) > 2.0_rp*q2_bound, &
         'conjugated-eigenvector control separation', k)
      conj_cp=wrong_cp; conj_q2=wrong_q2
      write(*,*) 'M1M14_CASE', site, spin, 'TRUE', cp, q2p, 'ERROR', cp-c_direct, q2p-q2_direct, &
         'SCALAR', scalar_cp, scalar_q2, 'VOUT_DELTA', vout_cp-cp, vout_q2-q2p, &
         'CONJ_DELTA', conj_cp-cp, conj_q2-q2p, 'BOUNDS', c_bound, q2_bound
   end do
   close(unit)

   ! Exercise the production type/spin block extractor with tagged slots.
   mapped = cmplx(0.0_rp,0.0_rp,rp)
   do k=1,ncase
      open(newunit=unit, file=trim(path), status='old', action='read')
      read(unit,*) magic, nc
      read(unit,*) nc
      do i=1,k-1
         read(unit,*)
         do j=1,4
            call skip_matrix(unit)
         end do
      end do
      read(unit,*) site, spin, q0, c, q2, c_direct, q2_direct, c_bound, q2_bound
      call read_matrix(unit,p); call read_matrix(unit,m1); call read_matrix(unit,vin); call read_matrix(unit,vout)
      do i=1,5
         do j=1,5
            mapped(4+i+9*(spin-1),4+j+9*(spin-1),site)=vin(i,j)
         end do
      end do
      close(unit)
   end do

   do k=1,ncase
      do i=1,ncase
         if (sites(i) == sites(k) .or. spins(i) /= spins(k)) cycle
         call correct_d_shell_moments(ps(:,:,k),m1s(:,:,k),vins(:,:,i),q0s(k),cs(k),rawq2s(k), &
            wrong_cp,wrong_q2,t,status)
         call check(status == 0, 'site/type swap status', k)
         call check(abs(wrong_cp-centers(k)) > cbounds(k)+cbounds(i) .or. &
            abs(wrong_q2-q2s(k)) > qbounds(k)+qbounds(i), 'site/type swap separation', k)
         write(*,*) 'M1M14_SITE_SWAP', sites(k), spins(k), 'DELTA', wrong_cp-centers(k), &
            wrong_q2-q2s(k), 'SUM_BOUNDS', cbounds(k)+cbounds(i), qbounds(k)+qbounds(i)
      end do
   end do
   do k=1,ncase
      call extract_hubbard_d_shell_potential(mapped, sites(k), spins(k), 9, got, status)
      call check(status == 0, 'mapped type/spin status', k)
      ! Re-read the selected matrix from the fixture and compare the mapping.
      open(newunit=unit, file=trim(path), status='old', action='read')
      read(unit,*) magic, nc
      read(unit,*) nc
      do i=1,ncase
         read(unit,*) site, spin, q0, c, q2, c_direct, q2_direct, c_bound, q2_bound
         call read_matrix(unit,p); call read_matrix(unit,m1); call read_matrix(unit,vin); call read_matrix(unit,vout)
         if (site == sites(k) .and. spin == spins(k)) exit
      end do
      close(unit)
      call check(maxval(abs(got-vin)) <= 1.0e-14_rp, 'mapped operator block', k)
   end do

   ! Existing physicality policy rejects an appreciably negative eigenvalue;
   ! invalid q2 is also returned to the caller, never silently clamped.
   p = cmplx(0.0_rp,0.0_rp,rp); p(1,1)=cmplx(-0.01_rp,0.0_rp,rp)
   m1 = cmplx(0.0_rp,0.0_rp,rp); vin = m1
   call correct_d_shell_moments(p,m1,vin,1.0_rp,0.0_rp,1.0_rp,cp,q2p,t,status)
   call check(status == 1, 'non-PSD density rejection', 0)
   p = cmplx(0.0_rp,0.0_rp,rp); p(1,1)=cmplx(0.5_rp,0.0_rp,rp)
   p(1,2)=cmplx(0.0_rp,2.0e-3_rp,rp)
   call correct_d_shell_moments(p,m1,vin,1.0_rp,0.0_rp,1.0_rp,cp,q2p,t,status)
   call check(status == 1, 'anti-Hermitian density rejection', 0)
   p = cmplx(0.0_rp,0.0_rp,rp); p(1,1)=cmplx(0.5_rp,0.0_rp,rp)
   call correct_d_shell_moments(p,m1,vin,1.0_rp,0.0_rp,-0.01_rp,cp,q2p,t,status)
   call check(status == 4, 'negative q2 rejection', 0)
   call correct_d_shell_moments(p,m1,vin,1.0_rp,0.123_rp,0.456_rp,cp,q2p,t,status)
   call check(status == 0 .and. cp == 0.123_rp .and. q2p == 0.456_rp .and. t == 0.0_rp, &
      'zero-potential exact no-op', 0)
   call check(.not. fixed_impurity_atom_carries_uj([0.0_rp, 0.0_rp], [0.0_rp, 0.0_rp]), &
      'zero-U/J impurity is not a correction-bearing atom', 0)
   call check(fixed_impurity_atom_carries_uj([0.0_rp, 0.4_rp], [0.0_rp, 0.0_rp]), &
      'nonzero impurity U retains local correction', 0)
   call check(fixed_impurity_atom_carries_uj([0.0_rp, 0.0_rp], [0.0_rp, 0.2_rp]), &
      'nonzero impurity J retains local correction', 0)
   if (failed) error stop 'local moment kernel failed'
   write(*,'(a)') 'RESULT: PASS'
contains
   subroutine read_matrix(u,a)
      integer, intent(in) :: u
      complex(rp), intent(out) :: a(5,5)
      integer :: r,c
      real(rp) :: re,im
      do c=1,5
         do r=1,5
            read(u,*,iostat=ios) re,im
            if (ios /= 0) error stop 'truncated fixture matrix'
            a(r,c)=cmplx(re,im,rp)
         end do
      end do
   end subroutine
   subroutine skip_matrix(u)
      integer,intent(in)::u
      integer::r
      real(rp)::re,im
      do r=1,25
         read(u,*,iostat=ios) re,im
         if (ios /= 0) error stop 'truncated fixture while skipping'
      end do
   end subroutine
   subroutine check(ok,label,which)
      logical,intent(in)::ok
      character(len=*),intent(in)::label
      integer,intent(in)::which
      if (.not.ok) then
         write(*,'(a,1x,a,1x,i0)') 'FAIL',label,which
         failed=.true.
      end if
   end subroutine
end program test_local_moments_kernel
