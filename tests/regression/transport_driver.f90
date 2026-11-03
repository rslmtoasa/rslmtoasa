! Binary fixture runner for the actual production matrix-free transport action.
program transport_driver
 use precision_mod, only: rp
 use transport_operators_mod
 implicit none
 integer, parameter :: b=18,n=3,d=b*n
 integer,target :: nn(n,4),iz(n)
 real(rp),target :: cr(3,n),scale(n)
 complex(rp),target :: h(b,b,4,n),hl(b,b,4,n),o(b,b,n),e(b,b,n),s(b,b,n)
 complex(rp) :: psi(b,d,n),out(b,d,n),q(b,b)
 real(rp) :: direction(3)
 integer :: ih,iw,im,nq,k,j,u
 character(16) :: kinds(6)=[character(16)::'charge','local','current','spin_torque','orbital_torque','soc_torque']
 type(transport_operators) :: t
 open(newunit=u,file='fixture.bin',access='stream',form='unformatted',status='old')
 read(u) nn,iz,cr,scale,h,hl,o,e,s,psi,direction,nq
 t%nn=>nn; t%iz=>iz; t%cr=>cr; t%scale=>scale
 t%h=>h; t%h_local=>hl; t%overlap=>o; t%enu=>e; t%soc=>s
 t%alat=1.3_rp
 open(unit=20,file='result.bin',access='stream',form='unformatted',status='replace')
 do k=1,nq
  read(u) q
  do im=0,1
   t%nmax=2*im
   do ih=0,1
    t%hoh=ih==1
    call t%hamiltonian(psi,out)
    write(20) out
    do iw=0,1
     do j=1,size(kinds)
      call t%response(trim(kinds(j)),q,direction,iw==1,psi,out)
      write(20) out
     end do
    end do
   end do
  end do
 end do
 ! Reuse the same actions on zero input to detect stale scratch accumulation.
 psi=(0.0_rp,0.0_rp)
 do j=1,size(kinds)
  call t%response(trim(kinds(j)),q,direction,.true.,psi,out)
  if (maxval(abs(out)) > 1.e-14_rp) error stop 'Nonzero action on zero input'
 end do
 close(u)
 close(20)
end program
