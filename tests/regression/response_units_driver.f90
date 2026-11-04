! Production conductivity normalization, integration and output fixture.
program response_units_driver
 use precision_mod, only: rp
 use conductivity_mod, only: conductivity, conductivity_integral, response_operator_units
 use control_mod, only: control
 use lattice_mod, only: lattice
 use energy_mod, only: energy
 use recursion_mod, only: recursion
 use hamiltonian_mod, only: hamiltonian
 use logger_mod, only: g_logger
 implicit none
 type(conductivity) :: c
 type(control), target :: ctl
 type(lattice), target :: lat
 type(energy), target :: en
 type(recursion), target :: rec
 type(hamiltonian), target :: ham
 real(rp) :: x(12), y(12), weights(2), conversion, height
 complex(rp) :: at(18,18,12,2), total(18,18,12)
 character(len=160) :: units
 character(len=32) :: name, argument
 integer :: i, j, pe, pl, rank, u
 call g_logger%init()
 c%control=>ctl; c%lattice=>lat; c%en=>en; c%recursion=>rec
 lat%control=>ctl
 rec%hamiltonian=>ham; ham%pol_alpha="z"; ham%pol_beta="z"
 rec%transport_dir_a=[0.0_rp,1.0_rp,0.0_rp]; rec%transport_dir_b=[1.0_rp,0.0_rp,0.0_rp]
 ctl%fname='input.nml'; ctl%linear_in='charge'; ctl%linear_out='charge'
 call get_command_argument(1,argument)
 if (len_trim(argument)>0) ctl%linear_out=trim(argument)
 call get_command_argument(2,argument)
 if (len_trim(argument)>0) ctl%linear_in=trim(argument)
 height=2
 call get_command_argument(3,argument)
 if (len_trim(argument)>0) read(argument,*) height
 lat%kk=5; lat%ntot=5; lat%ntype=2; lat%n1=1; lat%n2=1; lat%n3=1
 lat%pbc=.true.; lat%b1=.true.; lat%b2=.true.; lat%b3=.true.
 lat%a=0; lat%a(1,1)=1; lat%a(2,2)=1; lat%a(3,3)=height
 lat%alat=2; lat%iz=[1,1,1,2,2]
 call get_command_argument(4,argument)
 if (trim(argument)=="open") lat%pbc=.false.
 ctl%cond_calctype='per_type'; ctl%random_vec_num=2; ctl%cond_ll=4
 call c%restore_to_default()
 call c%build_from_file()
 call c%response_normalization(weights,conversion,units)
 open(newunit=u,file='units.txt',status='replace')
 write(u,'(a)') trim(units)
 write(u,'(3es25.16)') weights,conversion
 call response_operator_units(trim(ctl%linear_out),pe,pl,rank)
 write(u,*) pe,pl,rank
 close(u)
 do i=1,12
  x(i)=-.7_rp+real(i-1,rp)*1.4_rp/11
 end do
 ! Integrate several polynomials on an explicitly nonuniform grid; every
 ! prefix includes the formerly out-of-bounds final and partial panels.
 open(newunit=u,file='quadrature.txt',status='replace')
 do j=0,2
  y=(x+1)**2
  y=y**j
  do i=1,12
   write(u,'(i3,i3,es25.16)') j,i,conductivity_integral(y,(x+1)**2,i)
  end do
 end do
 close(u)
 en%channels_ldos=2; en%nv1=3; en%energy_min=-1.5_rp; en%energy_max=1.9_rp; en%fermi=.1_rp
 en%ene=2*x+.2_rp
 call c%calculate_gamma_nm()
 open(newunit=u,file='gamma.bin',access='stream',form='unformatted',status='replace')
 write(u) c%gamma_nm
 close(u)
 ! A polynomial spectral integrand gives an analytic complete pipeline check.
 at=0
 at(1,1,:,1)=1+2*x+3*x*x
 at(2,2,:,2)=3*(1+2*x+3*x*x)
 total=at(:,:,:,1)+at(:,:,:,2)
 call c%write_physical_response(total,at,x,2.0_rp,.2_rp)
 ! Exercise the raw pipeline too, without requiring atomic potential files.
 ctl%cond_calctype='random_vec'
 allocate(rec%mu_nm_stochastic(18,18,4,4,2))
 rec%mu_nm_stochastic=0
 rec%mu_nm_stochastic(1,1,1,1,1)=1
 rec%mu_nm_stochastic(2,2,1,1,2)=3
 c%gamma_nm=0
 c%gamma_nm(:,1,1)=1+2*x+3*x*x
 call c%calculate_conductivity_tensor()
end program
