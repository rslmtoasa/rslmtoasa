program test_local_moments_site_swap
   use precision_mod, only: rp
   use math_mod, only: correct_d_shell_moments
   implicit none
   integer, parameter :: ncase=4
   integer :: u, ios, nc, k, site, spin, status
   character(len=64) :: magic
   character(len=1024) :: fixture
   real(rp) :: q0,c,q2,cp_direct,q2_direct,cb,q2b,cp_wrong,q2_wrong,cwb,q2wb
   real(rp) :: cp,q2p,t,cw,qw
   complex(rp) :: p(5,5),m1(5,5),vin(5,5),vwrong(5,5)
   logical :: failed
   failed=.false.
   call get_command_argument(1,fixture)
   if (len_trim(fixture)==0) error stop 'fixture path required'
   open(newunit=u,file=trim(fixture),status='old',action='read',iostat=ios)
   if (ios /= 0) error stop 'cannot open fresh held-P site-swap fixture'
   read(u,*,iostat=ios) magic,nc
   if (ios /= 0 .or. trim(magic)/='A3M14SWAP' .or. nc/=1) error stop 'bad fixture version'
   read(u,*,iostat=ios) nc
   if (ios /= 0 .or. nc/=ncase) error stop 'bad fixture case count'
   do k=1,ncase
      read(u,*,iostat=ios) site,spin,q0,c,q2,cp_direct,q2_direct,cb,q2b,cp_wrong,q2_wrong,cwb,q2wb
      if (ios /= 0) error stop 'bad site-swap case header'
      call read_matrix(u,p); call read_matrix(u,m1); call read_matrix(u,vin); call read_matrix(u,vwrong)
      call correct_d_shell_moments(p,m1,vin,q0,c,q2,cp,q2p,t,status)
      if (status/=0) then
         write(*,*) 'FAIL correct case status',site,spin,status
         failed=.true.
      end if
      if (abs(cp-cp_direct)>cb .or. abs(q2p-q2_direct)>q2b) then
         write(*,*) 'FAIL direct correct oracle',site,spin,cp-cp_direct,q2p-q2_direct,cb,q2b
         failed=.true.
      end if
      call correct_d_shell_moments(p,m1,vwrong,q0,c,q2,cw,qw,t,status)
      if (status/=0) then
         write(*,*) 'FAIL alternate-V status',site,spin,status
         failed=.true.
      end if
      if (abs(cw-cp_wrong)>cwb .or. abs(qw-q2_wrong)>q2wb) then
         write(*,*) 'FAIL alternate-V independent prediction',site,spin,cw-cp_wrong,qw-q2_wrong,cwb,q2wb
         failed=.true.
      end if
      if (abs(cw-cp)<=cb+cwb .and. abs(qw-q2p)<=q2b+q2wb) then
         write(*,*) 'FAIL alternate-V detectability',site,spin,abs(cw-cp),cb+cwb,abs(qw-q2p),q2b+q2wb
         failed=.true.
      end if
      write(*,*) 'M1M14_FRESH_SWAP',site,spin,'TRUE',cp,q2p,'TRUE_ERR',cp-cp_direct,q2p-q2_direct, &
         'WRONG',cw,qw,'WRONG_PRED_ERR',cw-cp_wrong,qw-q2_wrong,'DELTA',cw-cp,qw-q2p, &
         'BOUNDS',cb,q2b,cwb,q2wb
   end do
   close(u)
   if (failed) error stop 'fresh held-P alternate-V site-swap test failed'
   write(*,'(a)') 'RESULT: PASS'
contains
   subroutine read_matrix(unit_id,a)
      integer,intent(in)::unit_id
      complex(rp),intent(out)::a(5,5)
      integer::i,j,read_status
      real(rp)::re,im
      do j=1,5
         do i=1,5
            read(unit_id,*,iostat=read_status) re,im
            if (read_status/=0) error stop 'truncated fresh fixture matrix'
            a(i,j)=cmplx(re,im,rp)
         end do
      end do
   end subroutine
end program test_local_moments_site_swap
