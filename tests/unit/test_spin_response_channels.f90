! spin_response channels s, p, d on bcc Fe, from the scratch deck of fe_scratch.py (which also ran the spd driver files).
! Usage: <oracle directory>. Per-channel occupations vs band_moments, E1 summed over channels, spd with U on d only vs d only.
program test_spin_response_channels
   use precision_mod, only: rp
   use math_mod, only: pi
   use control_mod, only: control
   use lattice_mod, only: lattice
   use charge_mod, only: charge
   use energy_mod, only: energy
   use hamiltonian_mod, only: hamiltonian
   use reciprocal_mod, only: reciprocal
   use spin_response_mod, only: spin_response
   use spin_response_fe_state
   implicit none
   type(control), target :: control_obj
   type(lattice), target :: lattice_obj
   type(charge), target :: charge_obj
   type(energy), target :: energy_obj
   type(hamiltonian), target :: ham
   type(reciprocal), target :: rec, rec3, rec_ref
   type(spin_response) :: sr, sr3
   character(len=512) :: dir
   character(len=22), parameter :: dfile(4) = [character(len=22) :: 'spin_response_q001.dat', 'spin_response_q002.dat', &
                                               'spin_response_q003.dat', 'static_q002.dat']
   character(len=22), parameter :: sfile(4) = [character(len=22) :: 'spd_q001.dat', 'spd_q002.dat', 'spd_q003.dat', 'spd_static_q002.dat']
   real(rp), allocatable :: bm(:, :, :, :), om(:), xd(:), xs(:)
   complex(rp), allocatable :: c3(:, :, :)
   real(rp) :: x, q0, cen, n_ref, m_ref, dev(4)
   integer :: l, s, nw, iw, i

   call get_command_argument(1, dir)
   if (len_trim(dir) == 0) error stop 'usage: test_spin_response_channels <oracle directory>'
   call read_tolerances(trim(dir))
   call require('channel_d_rel', channel_d_rel)
   call build_fe_state(.false., control_obj, lattice_obj, charge_obj, energy_obj, ham, 'input_driver.nml')
   rec = reciprocal(ham)
   rec%fermi_level = energy_obj%fermi
   sr = spin_response(rec)
   rec3 = reciprocal(ham)
   rec3%fermi_level = energy_obj%fermi
   sr3 = spin_response(rec3)
   sr3%channel_l = [0, 1, 2]
   call sr%prepare()
   call sr3%prepare()
   call reference_band_moments(ham, rec_ref, bm)

   ! Mills and Juelich occupation of each channel vs band_moments(site, l+1, spin, 1:3); the Juelich one is rebuilt from the
   ! moments about E_nu as in C2.3. Catches a shifted orbital range or channel index, and E_nu or p of another channel.
   m_ref = 0.0_rp
   do l = 0, 2
      do s = 1, 2
         associate (pot => lattice_obj%symbolic_atoms(lattice_obj%nbulk + 1)%potential)
            q0 = bm(1, l + 1, s, 1)
            cen = bm(1, l + 1, s, 2) - pot%enu(l, s)
            x = rec_ref%fermi_level - pot%enu(l, s)
            n_ref = (q0 + 2.0_rp*x*q0*cen/pot%ppar(l, s)**2 + x*x*q0*(bm(1, l + 1, s, 3)**2 + cen*cen)/pot%ppar(l, s)**4) &
                    /(1.0_rp + x*x/pot%ppar(l, s)**2)
         end associate
         call report('channel '//achar(48 + l)//' spin '//achar(48 + s)//' Mills |n - q0|/q0', abs(sr3%occ_mills(l + 1, s) - q0)/q0, moment_rel)
         call report('channel '//achar(48 + l)//' spin '//achar(48 + s)//' Juelich |n - ref|/ref', abs(sr3%occ_juelich(l + 1, s) - n_ref)/n_ref, moment_rel)
      end do
      m_ref = m_ref + bm(1, l + 1, 1, 1) - bm(1, l + 1, 2, 1)
   end do

   ! E1: integral of tr L over |omega| <= 5 Ry, eta = 0.02, Mills amplitudes, summed over channels, vs sum of the channel moments.
   nw = 2001
   allocate (om(nw), c3(3, 3, nw))
   om = [(-5.0_rp + (iw - 1)*10.0_rp/(nw - 1), iw=1, nw)]
   call sr3%bare_response([-1.0_rp, 1.0_rp, 1.0_rp]/12.0_rp, om, 0.02_rp, .false., c3)
   x = sum(0.5_rp*(-aimag(c3(1, 1, 1:nw - 1) + c3(2, 2, 1:nw - 1) + c3(3, 3, 1:nw - 1)) &
                   - aimag(c3(1, 1, 2:) + c3(2, 2, 2:) + c3(3, 3, 2:)))*(om(2:) - om(1:nw - 1)))/pi
   call report('E1 |int tr L - sum_l M_l|/sum_l M_l, spd, eta 0.02', abs(x - m_ref)/abs(m_ref), e1_fe_rel)

   ! spd with U on d only vs d only: U, d moment, and the smallest |eig(I + chi0 K)| of the driver files (q files and the static
   ! files at eta = 0), where the s and p eigenvalues are 1, so the expected value is min(1, d only).
   call sr%interaction()
   call sr3%interaction()
   call report('U_Juelich spd vs d, relative', abs(sr3%uj(1) - sr%uj(1))/sr%uj(1), channel_d_rel)
   call report('M_Juelich spd vs d, relative', abs(sr3%moment(1) - sr%moment(1))/sr%moment(1), channel_d_rel)
   do i = 1, 4
      xd = min(1.0_rp, eig_column(trim(dfile(i))))
      xs = eig_column(trim(sfile(i)))
      dev(i) = maxval(abs(xs - xd)/xd)
      call report('min|eig(I+chi0 K)| spd vs min(1, d), '//trim(sfile(i))//', max relative', dev(i), channel_d_rel)
   end do
   call finish()

contains

   !> 7th column of a q file.
   function eig_column(file) result(x)
      character(*), intent(in) :: file
      real(rp), allocatable :: x(:)
      real(rp) :: row(7), buf(400)
      integer :: u, n, ios

      n = 0
      open (newunit=u, file=file, status='old', action='read')
      read (u, *)
      do
         read (u, *, iostat=ios) row
         if (ios /= 0) exit
         n = n + 1
         buf(n) = row(7)
      end do
      close (u)
      x = buf(1:n)
   end function eig_column

end program test_spin_response_channels
