! Binary fixture runner for the actual production matrix-free transport action.
program transport_driver
 use precision_mod, only: rp
 use recursion_mod, only: recursion
 use hamiltonian_mod, only: hamiltonian
 use lattice_mod, only: lattice
 implicit none
 integer, parameter :: b=18,n=3,d=b*n
 integer,target :: nn(n,4),iz(n)
 real(rp),target :: cr(3,n),scale(n)
 complex(rp),target :: h(b,b,4,n),hl(b,b,4,n),o(b,b,n),e(b,b,n),s(b,b,n)
 complex(rp) :: psi(b,d,n),out(b,d,n),q(b,b)
 real(rp) :: direction(3)
 integer :: ih,iw,im,nq,k,j,u,ncases
 character(32) :: keyword
 character(10) :: pol
 character(16) :: kinds(6)=[character(16)::'charge','local','current','spin_torque','orbital_torque','soc_torque']
 type(recursion) :: t
 type(hamiltonian), target :: ham
 type(lattice), target :: lat
 complex(rp) :: work(b,d,n,6)
 open(newunit=u,file='fixture.bin',access='stream',form='unformatted',status='old')
 read(u) nn,iz,cr,scale,h,hl,o,e,s,psi,direction,nq
 t%lattice => lat
 t%hamiltonian => ham
 lat%nn=nn; lat%iz=iz; lat%cr=cr; lat%alat=1.3_rp; lat%kk=n
 ham%ee=h; ham%hall=hl; ham%obarm=o; ham%enim=e; ham%lsham=s
 ham%velocity_scale=scale
 allocate(t%idum(n))
 open(unit=20,file='result.bin',access='stream',form='unformatted',status='replace')
 do k=1,nq
  read(u) q
  do im=0,1
   lat%nmax=2*im
   do ih=0,1
    ham%hoh=ih==1
    call t%kubo_ham_vec_matmul(psi,out,work)
    write(20) out
    do iw=0,1
     do j=1,size(kinds)
      call t%kubo_operator_vec_matmul(trim(kinds(j)),q,direction,iw==1,psi,out,work)
      write(20) out
     end do
    end do
   end do
  end do
 end do
 ! Exercise both public response slots with independently normalized directions.
 read(u) ncases
 ham%v_alpha=2.0_rp*direction
 ham%v_beta=-3.0_rp*direction
 lat%nmax=2
 do k=1,ncases
  read(u) keyword,pol
  ham%pol_alpha=pol
  ham%pol_beta=pol
  do ih=0,1
   ham%hoh=ih==1
   call t%setup_kubo_operators(trim(keyword),'charge')
   call t%apply_kubo_operator('a',psi,out,work)
   write(20) out
   call t%setup_kubo_operators('charge',trim(keyword))
   call t%apply_kubo_operator('b',psi,out,work)
   write(20) out
  end do
 end do
 ! Reuse the same actions on zero input to detect stale scratch accumulation.
 psi=(0.0_rp,0.0_rp)
 do j=1,size(kinds)
  call t%kubo_operator_vec_matmul(trim(kinds(j)),q,direction,.true.,psi,out,work)
  if (maxval(abs(out)) > 1.e-14_rp) error stop 'Nonzero action on zero input'
 end do
 close(u)
 close(20)
end program
