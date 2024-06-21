!
module abSpinlib
   !
   use Depondt
   !
   real(selected_real_kind(15)), dimension(:, :), allocatable :: asd_cfield, asd_cfield_diff
   real(selected_real_kind(15)) :: asd_cf_tol
   real(selected_real_kind(15)), parameter :: ry2Tesla = 1.0d0 !235051.946397547d0
   real(selected_real_kind(15)), dimension(:, :), allocatable :: beff, b2eff, btorque, she_btorque
   real(selected_real_kind(15)), dimension(:, :), allocatable :: thermal_field
   real(selected_real_kind(15)), dimension(:), allocatable :: bxc
   real(selected_real_kind(15)), dimension(:, :), allocatable :: emom, emom2, emomM
   real(selected_real_kind(15)), dimension(:), allocatable :: mmom
   real(selected_real_kind(15)), dimension(:), allocatable :: asd_damp_array
   real(selected_real_kind(15)), dimension(:), allocatable :: temp_array
   character(len=1) :: do_stt = 'N'
   character(len=1) :: do_she = 'N'
   logical :: do_asd_pred
   logical :: do_asd_corr
   logical :: asd_cfd_conv
   logical :: scf_conv
   logical :: xc_int = .false. !.true.
   real(selected_real_kind(15))  :: lambda1 = 0.1d0
   real(selected_real_kind(15))  :: delta_t = 1.0d-16
   real(selected_real_kind(15))  :: temprescale = 1.0d0
   real(selected_real_kind(15))  :: asd_temp, asd_dt, asd_damp, asd_alpha, asd_alpha_org
   real(selected_real_kind(15))  :: asd_tote, asd_bande, asd_bande2
   integer :: asd_step
   integer :: asd_atom
   integer :: tseed = 1
   !
   logical, dimension(:), allocatable :: is_induced
   real(selected_real_kind(15)) :: induced_thresh = 0.5d0
   !
   private
   !
   public :: allocate_asd, asd_corr, asd_pred
   public :: asd_temp, asd_dt, read_asd_input
   !
contains

   subroutine allocate_asd(na, flag)
      !
      use RandomNumbers, only: mt_ran_init_c, setup_rng_hb, allocate_randomwork
      use Depondt, only: allocate_depondtfields
      !
      implicit none
      !
      integer, intent(in) :: na
      integer, intent(in) :: flag
      !
      integer :: i_stat, i_all
      if (flag > 0) then
         !print *,´Allocating for ASD´, na
         allocate (bxc(na), stat=i_stat)
         allocate (mmom(na), stat=i_stat)
         allocate (temp_array(na), stat=i_stat)
         allocate (beff(3, na), stat=i_stat)
         allocate (b2eff(3, na), stat=i_stat)
         allocate (btorque(3, na), stat=i_stat)
         allocate (she_btorque(3, na), stat=i_stat)
         allocate (thermal_field(3, na), stat=i_stat)
         allocate (asd_damp_array(na), stat=i_stat)
         allocate (emom(3, na), stat=i_stat)
         allocate (emom2(3, na), stat=i_stat)
         allocate (emomM(3, na), stat=i_stat)
         allocate (is_induced(na), stat=i_stat)
         !
         ! Initialize random number generators
         call setup_rng_hb(tseed, 'Y', 'N')
         call allocate_randomwork(na, 1, 1, 'N')
         call allocate_depondtfields(na, 1, 1)
         call mt_ran_init_c(tseed) ! Monte-Carlo

         open (2600, file='moments.asd.dat', status='replace')
         open (2605, file='moments.pred.dat', status='replace')
         open (2610, file='moments.corr.dat', status='replace')
         !
         asd_step = 0
         !
      else
         deallocate (bxc, stat=i_stat)
         deallocate (mmom, stat=i_stat)
         deallocate (temp_array, stat=i_stat)
         deallocate (beff, stat=i_stat)
         deallocate (b2eff, stat=i_stat)
         deallocate (btorque, stat=i_stat)
         deallocate (emom, stat=i_stat)
         deallocate (emom2, stat=i_stat)
         deallocate (emomM, stat=i_stat)
         deallocate (asd_damp_array, stat=i_stat)
         deallocate (is_induced, stat=i_stat)
         call allocate_randomwork(Na, 1, -1, 'N')
         call allocate_depondtfields(na, 1, -1)
      end if
      !
   end subroutine allocate_asd
   !
   subroutine asd_pred(inmom, infield, intemp, indt, na)
      !
      implicit none
      !
      integer, intent(in) :: na
      real(selected_real_kind(15)), dimension(3, na), intent(inout) :: inmom
      real(selected_real_kind(15)), dimension(3, na), intent(in) :: infield
      real(selected_real_kind(15)), intent(in) :: intemp
      real(selected_real_kind(15)), intent(in) :: indt
      !
      integer :: i, j, natom
      real(selected_real_kind(15)) :: mnorm
      !
      print *, ' asd_pred', asd_damp, na
      asd_step = asd_step + 1
      !mmom=1.0d0
      !emom=0.0d0
      emom2 = 0.0d0
      natom = na
      do i = 1, natom
         mnorm = sqrt(sum(inmom(:, i)*inmom(:, i)))
         mmom(i) = mnorm
         do j = 1, 3
            emom(j, i) = inmom(j, i)/mnorm
            emomM(j, i) = inmom(j, i)!*mmom(i)
            ! Careful about the conversion from input field
            beff(j, i) = -infield(j, i)*ry2Tesla/2.0d0!*mmom(i)  ! *0.5 for lande-g
            !beff(j,i)=-infield(j,i)/2.0d0!*ry2Tesla
         end do
         write (2600, '(2i5,3f12.6,2x,3f18.6,2x,3f18.6)') asd_step, i, emom(:, i), beff(:, i) &
            , asd_tote, asd_bande, asd_bande2
      end do
      !      emomm=emom
      b2eff = 0.0d0
      btorque = 0.0d0
      do_stt = 'N'
      !mmom=1.0d0
      !
      temp_array = intemp
      delta_t = indt
      asd_damp_array = asd_damp
      !
      !!! print *,´-> beff´,beff
      !!! print *,´-> b2eff´,b2eff
      !!! print *,´-> btorque´,btorque
      !!! print *,´-> emom´,emom
      !!! print *,´-> emom2´,emom2
      !!! print *,´-> emomM´,emomM
      !!! print *,´-> mmom´,mmom
      call depondt_evolve_first(Natom, 1, asd_damp_array, beff, b2eff, btorque, emom, emom2, emomM, mmom, delta_t, Temp_array, &
                                temprescale, do_stt, thermal_field, do_she, she_btorque)
      !!! print *,´---> beff´,beff
      !!! print *,´---> b2eff´,b2eff
      !!! print *,´---> btorque´,btorque
      !!! print *,´---> emom´,emom
      !!! print *,´---> emom2´,emom2
      !!! print *,´---> emomM´,emomM
      !!! print *,´---> mmom´,mmom
      !
      do_asd_pred = .false.
      do_asd_corr = .true.
      do i = 1, natom
         if (.not. is_induced(i)) then
            mnorm = sqrt(sum(inmom(:, i)*inmom(:, i)))
            do j = 1, 3
               inmom(j, i) = emom(j, i)*mnorm
            end do
         else
            do j = 1, 3
               emom(j, i) = inmom(j, i)/mnorm
            end do
         end if
         write (2605, '(2i5,3f12.6,2x,3f18.6)') asd_step, i, emom(:, i), beff(:, i)
      end do
      rewind (7)
      do i = 1, na
         write (7, 10001) inmom(1:3, i)
         !         write(*,10001) inmom(1:3,i)
      end do
      rewind (7)
      !
10001 format(1x, 3f14.9)
   end subroutine asd_pred

   subroutine asd_corr(inmom, infield, intemp, indt, na)
      !
      implicit none
      !
      integer, intent(in) :: na
      real(selected_real_kind(15)), dimension(3, na), intent(inout) :: inmom
      real(selected_real_kind(15)), dimension(3, na), intent(in) :: infield
      real(selected_real_kind(15)), intent(in) :: intemp
      real(selected_real_kind(15)), intent(in) :: indt
      !
      integer :: i, j, natom
      real(selected_real_kind(15)) :: mnorm
      !
      !      emom=0.0d0
      !      emom2=0.0d0
      !mmom=1.0d0
      print *, ' asd_corr', asd_damp
      natom = na
      do i = 1, natom
         mnorm = sqrt(sum(inmom(:, i)*inmom(:, i)))
         do j = 1, 3
            mmom(i) = mnorm
            emom(j, i) = inmom(j, i)/mnorm
            emomm(j, i) = inmom(j, i)!*mmom(i)
            !emomm(j,i)=inmom(j,i)*1.0d0!*mmom(i)
            ! Careful with unit of input field!
            beff(j, i) = -infield(j, i)*ry2Tesla/2.0d0!*mmom(i)   ! *0.5 for lande-g
            !beff(j,i)=-infield(j,i)/2.0d0
         end do
      end do
!     emomm=emom*mnor
      !      b2eff=0.0d0
      do_stt = 'N'
      btorque = 0.0d0
      do_she = 'N'
      she_btorque = 0.0d0
      !mmom=1.0d0
      !
      temp_array = intemp
      delta_t = indt
      asd_damp_array = asd_damp
      !
      !!! print *,´=> beff´,beff
      !!! print *,´=> b2eff´,b2eff
      !!! print *,´=> btorque´,btorque
      !!! print *,´=> emom´,emom
      !!! print *,´=> emom2´,emom2
      call depondt_evolve_second(Natom, 1, asd_damp_array, beff, b2eff, btorque, emom, emom2, delta_t, do_stt, do_she, she_btorque)
      !!! print *,´===> beff´,beff
      !!! print *,´===> b2eff´,b2eff
      !!! print *,´===> btorque´,btorque
      !!! print *,´===> emom´,emom
      !!! print *,´===> emom2´,emom2
      !!! print *,´asd_corr´,natom
      do i = 1, natom
         mnorm = sqrt(sum(inmom(:, i)*inmom(:, i)))
         if (.not. is_induced(i)) then
            do j = 1, 3
               inmom(j, i) = emom2(j, i)*mnorm
            end do
         else
            do j = 1, 3
               emom2(j, i) = inmom(j, i)/mnorm
            end do
         end if
         write (2610, '(2i5,3f12.6,2x,3f18.6)') asd_step, i, emom2(:, i), beff(:, i)
      end do
      rewind (7)
      do i = 1, na
         write (7, 10001) inmom(1:3, i)
         !         write(*,10001) inmom(1:3,i)
      end do
      rewind (7)
      !
      do_asd_pred = .true.
      do_asd_corr = .false.
      !
      !
10001 format(1x, 3f14.9)
   end subroutine asd_corr

   subroutine asd_check_induced(natom, mmom)
      !
      implicit none
      !
      integer, intent(in) :: natom
      real, dimension(natom) :: mmom
      !
      integer :: iatom
      !
      is_induced = .false.
      do iatom = 1, natom
         if (mmom(iatom) < induced_thresh) is_induced(iatom) = .true.
      end do
      !
      return
   end subroutine asd_check_induced

   subroutine read_asd_input()
      !
      implicit none
      !
      open (20, file='asd.in')
      read (20, *, END=100) asd_temp, asd_dt, asd_damp, asd_alpha!, xc_int
100   continue
      asd_alpha_org = asd_alpha
      close (20)
      !open(30,file=´asd.out´,status=´replace´)
      !write(30,´(a)´) ´---------------------------------------------------------------------------´
      !open(35,file=´asd.log´,position=´append´)
      !write(35,´(a)´) ´#  Iter.    Magnetization    Total_{energy}  Band_{energy}   RMS_{moments}   RMS_{energy}    RMS_{eband} [Ry]´
   end subroutine read_asd_input

end module abSpinlib
!

