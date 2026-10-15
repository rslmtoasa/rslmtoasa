! Private A2.4 harness: feed one frozen radial channel to the unchanged
! self_mod audit seam and emit EVAL plus normalized G and PHDFSR GP.
program a2_radial_mapping_driver
   use precision_mod, only: rp
   use self_mod, only: legacy_radial_fixture
   implicit none
   integer :: nr, l, ir
   real(rp) :: z, a, b, slope, energy
   real(rp), allocatable :: rofi(:), potential(:), g(:, :), gp(:, :), gpp(:, :)
   character(len=4096) :: input_path

   call get_command_argument(1, input_path)
   if (len_trim(input_path) == 0) error stop 'usage: driver input.txt'
   open(unit=10, file=trim(input_path), status='old', action='read')
   read(10, *) nr, l, z, a, b, slope
   allocate(rofi(nr), potential(nr))
   do ir = 1, nr
      read(10, *) rofi(ir), potential(ir)
   end do
   close(10)
   call legacy_radial_fixture(z, l, a, b, rofi, potential, energy, g, gp, gpp, slope)
   write(*, '(a,1x,es25.16e3)') 'EVAL', energy
   do ir = 1, nr
      write(*, '(4(es25.16e3,1x))') g(ir, 1), g(ir, 2), gp(ir, 1), gp(ir, 2)
   end do
end program a2_radial_mapping_driver
