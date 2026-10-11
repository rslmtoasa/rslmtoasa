! Synthetic channel-kernel test, usage: <oracle directory>. accumulate_chi0 on two rows kappa, one pair with a nonzero vertex
! (band 1 at k, band 2 at k+q; the others have an exactly zero spin factor), against chi0_ij = w df T_i conj(T_j)/(w + de + i eta)
! written here: ordinary values, a vertex of 1e-100 and df = 1e-30 (a threshold on either would drop the pair), 41 omega
! (blocks of 32 and 9). Then the offsets of channels p and f on two sites at nb = 32.
program test_spin_response_channels_kernel
   use precision_mod, only: rp
   use basis_mod, only: basis_init
   use spin_response_mod, only: accumulate_chi0, channel_amplitudes_coefficient
   use spin_response_tolerances_mod
   implicit none
   character(len=512) :: dir
   logical :: ok
   real(rp) :: err, scale(3) = [1.0_rp, 1.0e-100_rp, 1.0_rp], f1(3) = [1.0_rp, 1.0_rp, 1.0e-30_rp], om(41), e_k(1, 2), f_k(1, 2)
   complex(rp) :: a_k(1, 2, 3, 2, 2), a_kq(1, 2, 3, 2, 2), chi0(2, 2, 41), ref(2, 2, 41), t(2), v(64, 1, 1)
   complex(rp), allocatable :: a(:, :, :, :, :)
   integer :: c, i, j, mu, iw, ic, s, l

   call get_command_argument(1, dir)
   call read_tolerance_file(trim(dir))
   call require('a1_rel', a1_rel)
   om = [(0.4_rp*(iw - 1)/40, iw=1, 41)]
   do c = 1, 3
      a_k = (0.0_rp, 0.0_rp)
      a_kq = (0.0_rp, 0.0_rp)
      do i = 1, 2
         do mu = 1, 3
            a_k(1, 1, mu, i, 1) = scale(c)*cmplx(mu + i, mu - 0.5_rp*i, rp)
            a_kq(1, 2, mu, i, 2) = cmplx(1.0_rp + 0.3_rp*mu, 0.2_rp*i, rp)
         end do
         t(i) = sum(conjg(a_k(1, 1, :, i, 1))*a_kq(1, 2, :, i, 2))
      end do
      do iw = 1, 41
         ref(:, :, iw) = 0.5_rp*f1(c)/cmplx(om(iw) - 0.3_rp, 0.02_rp, rp)*spread(t, 2, 2)*spread(conjg(t), 1, 2)
      end do
      chi0 = (0.0_rp, 0.0_rp)
      f_k(1, :) = [f1(c), 0.7_rp]
      call accumulate_chi0(0.01_rp, [0.5_rp], reshape([0.1_rp, 0.2_rp], [1, 2]), f_k, a_k, reshape([0.5_rp, 0.4_rp], [1, 2]), &
                           reshape([0.9_rp, 0.0_rp], [1, 2]), a_kq, om, 0.02_rp, chi0)
      err = maxval(abs(chi0 - ref))/maxval(abs(ref))
      write (*, '(a,i1,t45,es10.3,2x,es10.3)') 'pair case ', c, err, a1_rel
      if (.not. err <= a1_rel) error stop 1
   end do
   call basis_init(3)
   v(:, 1, 1) = [(cmplx(real(i, rp), 0.0_rp, rp), i=1, 64)]
   call channel_amplitudes_coefficient(v, 2, [1, 3], a)
   ok = all(shape(a) == [1, 1, 7, 4, 2])
   do i = 1, 2
      do ic = 1, 2
         l = merge(1, 3, ic == 1)
         do s = 1, 2
            do mu = 1, 7
               ok = ok .and. real(a(1, 1, mu, (i - 1)*2 + ic, s), rp) == merge(real((i - 1)*32 + (s - 1)*16 + l*l + mu, rp), 0.0_rp, mu <= 2*l + 1)
            end do
         end do
      end do
   end do
   write (*, '(a,l2)') 'offsets, channels p and f, two sites, nb = 32: ', ok
   if (.not. ok) error stop 1
   write (*, '(a)') 'RESULT: PASS'
end program test_spin_response_channels_kernel
