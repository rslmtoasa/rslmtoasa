! Private frozen-state M0.3 clone replay. Input is a strictly ordered text
! stream produced by the hash-validating Python adapter; all physics kernels
! run in the unchanged self_mod fixture seams.
program accepted_m03_replay_driver
   use precision_mod, only: rp
   use self_mod, only: legacy_accepted_m03_newrho_fixture, legacy_accepted_m03_components_fixture
   implicit none
   integer :: lmax, nr, ncore, nval, ncase, fcore, txc, i, s, ell, m, icase, uin, uout, ios
   real(rp) :: z, a, b, tol, vhrmax, b_fsm
   real(rp), allocatable :: rofi(:), v(:, :), ec_seed(:), ev_seed(:), gc(:, :), c(:, :), qpar(:, :), srdel(:, :)
   real(rp), allocatable :: pl_seed(:, :), pl_live(:, :), ql(:, :, :, :), shifts(:, :)
   real(rp), allocatable :: rho_in(:, :), rho(:, :, :), d_rho(:, :, :), pl(:, :, :)
   real(rp), allocatable :: sumec(:, :), sumev(:, :), ec(:, :), ev(:, :), cin(:, :), cout(:, :), rvin(:, :), rvout(:, :)
   logical, allocatable :: do_calc(:)
   character(len=4096) :: input_path, output_path
   call get_command_argument(1, input_path)
   call get_command_argument(2, output_path)
   if (len_trim(input_path) == 0 .or. len_trim(output_path) == 0) error stop 'usage: accepted_m03_replay_driver input output'
   open(newunit=uin, file=trim(input_path), status='old', action='read', iostat=ios)
   if (ios /= 0) error stop 'cannot open M0.3 input'
   read(uin, *, iostat=ios) lmax, nr, ncore, nval, ncase, fcore, txc
   if (ios /= 0 .or. lmax /= 2 .or. nr /= 495 .or. ncore < 1 .or. nval /= 6 .or. ncase < 1) error stop 'bad M0.3 dimensions'
   read(uin, *, iostat=ios) z, a, b, tol, vhrmax, b_fsm
   if (ios /= 0) error stop 'bad M0.3 scalar controls'
   allocate(rofi(nr), v(nr,2), ec_seed(ncore), ev_seed(nval), gc(3,2), c(3,2), qpar(3,2), srdel(3,2), &
      pl_seed(3,2), pl_live(3,2), ql(3,3,2,ncase), shifts(2,ncase), do_calc(ncase), rho_in(nr,2), &
      rho(nr,2,ncase), d_rho(nr,2,ncase), pl(3,2,ncase), sumec(2,ncase), sumev(2,ncase), &
      ec(ncore,ncase), ev(nval,ncase), cin(10,ncase), cout(10,ncase), rvin(2,ncase), rvout(2,ncase))
   read(uin, *, iostat=ios) rofi
   if (ios /= 0) error stop 'bad mesh'
   do i=1,nr
      read(uin, *, iostat=ios) v(i,1),v(i,2),rho_in(i,1),rho_in(i,2)
      if (ios /= 0) error stop 'bad radial input'
   end do
   read(uin, *, iostat=ios) ec_seed
   if (ios /= 0) error stop 'bad EC seeds'
   read(uin, *, iostat=ios) ev_seed
   if (ios /= 0) error stop 'bad EV seeds'
   do s=1,2
      do ell=1,3
         read(uin, *, iostat=ios) gc(ell,s),c(ell,s),qpar(ell,s),srdel(ell,s),pl_seed(ell,s),pl_live(ell,s)
         if (ios /= 0) error stop 'bad PL operands'
      end do
   end do
   do icase=1,ncase
      read(uin, *, iostat=ios) shifts(:,icase),do_calc(icase)
      if (ios /= 0) error stop 'bad case header'
      do s=1,2
         do ell=1,3
            do m=1,3
               read(uin, *, iostat=ios) ql(m,ell,s,icase)
               if (ios /= 0) error stop 'bad QL case'
            end do
         end do
      end do
   end do
   close(uin)
   call legacy_accepted_m03_newrho_fixture(lmax,z,a,b,tol,fcore,rofi,v,ec_seed,ev_seed,gc,c,qpar,srdel, &
      pl_seed,pl_live,ql,shifts,do_calc,pl,rho,d_rho,sumec,sumev,ec,ev)
   call legacy_accepted_m03_components_fixture(z,a,b,nr,rofi,v,vhrmax,b_fsm,txc,rho_in,rho,sumec,sumev, &
      cin,cout,rvin,rvout)
   open(newunit=uout,file=trim(output_path),status='replace',action='write',iostat=ios)
   if (ios /= 0) error stop 'cannot open M0.3 output'
   do icase=1,ncase
      do s=1,2
         do ell=1,3
            write(uout,'(a,",",i0,",",i0,",",i0,",",es25.16e3)') 'PL',icase,ell-1,s,pl(ell,s,icase)
         end do
         write(uout,'(a,",",i0,",",i0,",",es25.16e3)') 'SUMEC',icase,s,sumec(s,icase)
         write(uout,'(a,",",i0,",",i0,",",es25.16e3)') 'SUMEV',icase,s,sumev(s,icase)
         write(uout,'(a,",",i0,",",i0,",",es25.16e3)') 'RVH_IN',icase,s,rvin(s,icase)
         write(uout,'(a,",",i0,",",i0,",",es25.16e3)') 'RVH_OUT',icase,s,rvout(s,icase)
         do i=1,nr
            write(uout,'(a,",",i0,",",i0,",",i0,",",es25.16e3,",",es25.16e3)') 'RHO',icase,i-1,s,rho(i,s,icase),d_rho(i,s,icase)
         end do
      end do
      do i=1,10
         write(uout,'(a,",",i0,",",i0,",",es25.16e3,",",es25.16e3)') 'COMPONENT',icase,i,cin(i,icase),cout(i,icase)
      end do
      do i=1,nval
         write(uout,'(a,",",i0,",",i0,",",es25.16e3)') 'EV',icase,i,ev(i,icase)
      end do
   end do
   close(uout)
end program accepted_m03_replay_driver
