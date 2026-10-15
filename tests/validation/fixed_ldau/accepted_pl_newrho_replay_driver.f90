! Private driver for exact accepted-state clone replay. The only numerical
! routines invoked are bands%calculate_pl and the unchanged self%newrho.
program accepted_pl_newrho_replay_driver
   use precision_mod, only: rp
   use self_mod, only: legacy_accepted_pl_newrho_fixture
   implicit none
   integer :: lmax, nr, ncore, nval, nshifts, fcore
   integer :: i, s, ell, m, iv, unit_in, unit_out, ios
   real(rp) :: z, a, b, tol
   real(rp), allocatable :: shifts(:), rofi(:), potential(:, :), ec_seed(:), ev_seed(:)
   real(rp), allocatable :: gravity_center(:, :), c_in(:, :), qpar(:, :), srdel(:, :), pl_seed(:, :)
   real(rp), allocatable :: ql(:, :, :), pl_out(:, :, :), rho_out(:, :, :)
   real(rp), allocatable :: sumec(:, :), sumev(:, :), ec_out(:, :), ev_out(:, :)
   character(len=4096) :: input_path, output_path

   call get_command_argument(1, input_path)
   call get_command_argument(2, output_path)
   if (len_trim(input_path) == 0 .or. len_trim(output_path) == 0) &
      error stop 'usage: accepted_pl_newrho_replay_driver input.txt output.csv'
   open(newunit=unit_in, file=trim(input_path), status='old', action='read', iostat=ios)
   if (ios /= 0) error stop 'cannot open accepted replay input'
   read(unit_in, *, iostat=ios) lmax, nr, ncore, nval, nshifts, fcore
   if (ios /= 0 .or. lmax /= 2 .or. nr /= 495 .or. ncore /= 10 .or. nval /= 6 .or. &
       (nshifts /= 1 .and. nshifts /= 9 .and. nshifts /= 6)) &
      error stop 'accepted replay dimensions do not match reviewed schema'
   read(unit_in, *, iostat=ios) z, a, b, tol
   if (ios /= 0) error stop 'invalid accepted replay scalar header'
   allocate(shifts(nshifts), rofi(nr), potential(nr, 2), ec_seed(ncore), ev_seed(nval))
   allocate(gravity_center(lmax + 1, 2), c_in(lmax + 1, 2), qpar(lmax + 1, 2), &
            srdel(lmax + 1, 2), pl_seed(lmax + 1, 2), ql(3, lmax + 1, 2))
   allocate(pl_out(lmax + 1, 2, nshifts), rho_out(nr, 2, nshifts), &
            sumec(2, nshifts), sumev(2, nshifts), ec_out(ncore, nshifts), ev_out(nval, nshifts))
   read(unit_in, *, iostat=ios) shifts
   if (ios /= 0) error stop 'invalid shift vector'
   read(unit_in, *, iostat=ios) rofi
   if (ios /= 0) error stop 'invalid radial mesh'
   do i = 1, nr
      read(unit_in, *, iostat=ios) potential(i, 1), potential(i, 2)
      if (ios /= 0) error stop 'invalid radial potential'
   end do
   read(unit_in, *, iostat=ios) ec_seed
   if (ios /= 0) error stop 'invalid core energy seeds'
   read(unit_in, *, iostat=ios) ev_seed
   if (ios /= 0) error stop 'invalid valence energy seeds'
   do s = 1, 2
      do ell = 0, lmax
         read(unit_in, *, iostat=ios) gravity_center(ell + 1, s), c_in(ell + 1, s), &
              qpar(ell + 1, s), srdel(ell + 1, s), pl_seed(ell + 1, s)
         if (ios /= 0) error stop 'invalid calculate_pl operand tuple'
      end do
   end do
   do s = 1, 2
      do ell = 0, lmax
         do m = 1, 3
            read(unit_in, *, iostat=ios) ql(m, ell + 1, s)
            if (ios /= 0) error stop 'invalid QL moments'
         end do
      end do
   end do
   close(unit_in)

   call legacy_accepted_pl_newrho_fixture(lmax, z, a, b, tol, fcore, rofi, potential, ec_seed, ev_seed, &
        gravity_center, c_in, qpar, srdel, pl_seed, ql, shifts, pl_out, rho_out, sumec, sumev, ec_out, ev_out)

   open(newunit=unit_out, file=trim(output_path), status='replace', action='write', iostat=ios)
   if (ios /= 0) error stop 'cannot open accepted replay output'
   do iv = 1, nshifts
      do s = 1, 2
         do ell = 0, lmax
            write(unit_out, '(a,",",i0,",",i0,",",i0,",",es25.16e3)') 'PL', iv, ell, s, pl_out(ell + 1, s, iv)
         end do
         write(unit_out, '(a,",",i0,",",i0,",",es25.16e3)') 'SUMEC', iv, s, sumec(s, iv)
         write(unit_out, '(a,",",i0,",",i0,",",es25.16e3)') 'SUMEV', iv, s, sumev(s, iv)
      end do
      do i = 1, ncore
         write(unit_out, '(a,",",i0,",",i0,",",es25.16e3)') 'EC', iv, i, ec_out(i, iv)
      end do
      do i = 1, nval
         write(unit_out, '(a,",",i0,",",i0,",",es25.16e3)') 'EV', iv, i, ev_out(i, iv)
      end do
      do s = 1, 2
         do i = 1, nr
            write(unit_out, '(a,",",i0,",",i0,",",i0,",",es25.16e3)') 'RHO', iv, i, s, rho_out(i, s, iv)
         end do
      end do
   end do
   close(unit_out)
end program accepted_pl_newrho_replay_driver
