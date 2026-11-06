program test_nonmagnetic_projection
   use precision_mod, only: rp
   use symbolic_atom_mod, only: symbolic_atom
   use potential_mod, only: average_spin_channels
   implicit none
   type(symbolic_atom) :: atom, magnetic
   real(rp) :: qsum(3, 0:2), psum(0:2), density(5, 2), total(5)
   integer :: i

   call atom%potential%restore_to_default()
   atom%potential%ql = reshape([(real(i, rp)/7, i=1,18)], [3,3,2])
   atom%potential%pl(:,1) = [4.4_rp, 4.5_rp, 3.6_rp]
   atom%potential%pl(:,2) = [4.6_rp, 4.7_rp, 3.8_rp]
   atom%potential%center_band(:,1) = 1.0_rp
   atom%potential%center_band(:,2) = 3.0_rp
   atom%potential%width_band = 2*atom%potential%center_band
   atom%potential%shifted_band = 3*atom%potential%center_band
   atom%potential%obar = 4*atom%potential%center_band
   atom%potential%mom = [0.6_rp, 0.0_rp, 0.8_rp]
   atom%potential%mom0 = 2*atom%potential%mom
   atom%potential%mom1 = atom%potential%mom
   atom%potential%mtot = 2.0_rp
   atom%mag_cfield = 5.0_rp
   magnetic = atom
   qsum = sum(atom%potential%ql, dim=3)
   psum = sum(atom%potential%pl, dim=2)
   atom%force_nonmagnetic = .true.
   call atom%build_pot()
   call magnetic%build_pot()
   if (maxval(abs(sum(atom%potential%ql,dim=3)-qsum)) > 1e-14_rp) error stop 'QL charge lost'
   if (maxval(abs(sum(atom%potential%pl,dim=2)-psum)) > 1e-14_rp) error stop 'PL sum lost'
   if (any(atom%potential%ql(:,:,1) /= atom%potential%ql(:,:,2))) error stop 'QL split'
   if (any(atom%potential%pl(:,1) /= atom%potential%pl(:,2))) error stop 'PL split'
   if (any(abs(atom%potential%cx1) > 0)) error stop 'cx1 split'
   if (any(abs(atom%potential%wx1) > 0)) error stop 'wx1 split'
   if (any(abs(atom%potential%cex1) > 0)) error stop 'cex1 split'
   if (any(abs(atom%potential%obx1) > 0)) error stop 'obx1 split'
   if (any(atom%potential%center_band /= 2.0_rp)) error stop 'Wrong parameter average'
   if (any(atom%potential%width_band /= 4.0_rp)) error stop 'Wrong width average'
   if (any(atom%potential%shifted_band /= 6.0_rp)) error stop 'Wrong shifted average'
   if (any(atom%potential%obar /= 8.0_rp)) error stop 'Wrong obar average'
   if (any(atom%mag_cfield /= 0)) error stop 'Constraining field survives'
   if (norm2(atom%potential%mom) /= 1) error stop 'Invalid reference axis'
   if (norm2(atom%potential%mom0) /= 0 .or. atom%potential%mtot /= 0) error stop 'Nonzero moment'
   if (any(magnetic%potential%ql /= reshape([(real(i,rp)/7,i=1,18)], [3,3,2]))) &
      error stop 'Magnetic occupations changed'
   if (any(abs(magnetic%potential%cx1) /= 1) .or. magnetic%potential%mtot /= 2) &
      error stop 'Magnetic atom changed'
   density = reshape([(real(i,rp)/3,i=1,10)], [5,2])
   total = sum(density,dim=2)
   call average_spin_channels(density)
   if (any(density(:,1) /= density(:,2))) error stop 'Density split'
   if (maxval(abs(sum(density,dim=2)-total)) > 1e-14_rp) error stop 'Density charge lost'
   print *, 'PASS: QL/PL and radial charge preserved; potentials degenerate; magnetic atom unchanged'
end program
