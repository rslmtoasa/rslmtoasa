! Binary fixture runner for the actual production matrix-free transport action.
program transport_driver
 use precision_mod, only: rp
 use math_mod, only: cone, czero
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
 call check_legacy_products()
 t%izero=1
 call check_masked_products()
 t%izero=1
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
contains

 ! Compare sparse support with the all-site action for every operator kind,
 ! both slots, both Hamiltonians, impurity rows and zero neighbor placeholders.
 subroutine check_masked_products()
  complex(rp) :: input(b,d,n),masked(b,d,n),full(b,d,n),scratch(b,d,n,6)
  complex(rp) :: qtest(b,b),prev(b,d,n),state(b,d,n),next(b,d,n)
  complex(rp) :: refprev(b,d,n),refstate(b,d,n),refnext(b,d,n)
  integer :: im,ih,iw,ik,step,i,mask(n)
  qtest=0
  do i=1,b
   qtest(i,i)=real(i,rp)/b
  end do
  qtest(1,2)=cmplx(0.2_rp,0.3_rp,rp)
  qtest(2,1)=conjg(qtest(1,2))
  input=psi
  input(:,:,2:n)=0
  do im=0,1
   lat%nmax=2*im
   do ih=0,1
    ham%hoh=ih==1
    do iw=0,1
     do ik=1,size(kinds)
      t%izero=1
      call t%kubo_operator_vec_matmul(trim(kinds(ik)),qtest,direction,iw==1,input,full,scratch)
      t%izero=[1,0,0]
      call t%kubo_operator_vec_matmul(trim(kinds(ik)),qtest,direction,iw==1,input,masked,scratch)
      if(maxval(abs(full-masked))>5.e-13_rp) error stop 'Masked observable lost a path'
      if(any(t%izero/=[1,0,0])) error stop 'Observable changed input map'
      mask=[1,1,0]
      if(ham%hoh) mask=[1,1,1]
      if(kinds(ik)=='local'.or.kinds(ik)=='soc_torque') mask=[1,0,0]
      if(any(t%idum/=mask)) error stop 'Incorrect observable output map'
     end do
    end do
    ! Repeated shifted/scaled Chebyshev steps; observable a must not change
    ! the recurrence support. This detects stale scratch and map side effects.
    state=input; prev=0; refstate=input; refprev=0; mask=[1,0,0]
    do step=1,6
     t%izero=1
     call t%kubo_ham_vec_matmul(refstate,refnext,scratch,a=1.7_rp,b=0.13_rp)
     t%izero=mask
     call t%kubo_ham_vec_matmul(state,next,scratch,a=1.7_rp,b=0.13_rp)
     if(any(t%izero/=mask)) error stop 'Hamiltonian changed input map'
     mask=t%idum
     if(step>1) then
      next=2*next-prev
      refnext=2*refnext-refprev
     end if
     if(maxval(abs(next-refnext))>5.e-10_rp) error stop 'Masked recurrence differs'
     prev=state; state=next; refprev=refstate; refstate=refnext
     t%izero=mask
     call t%kubo_operator_vec_matmul('current',qtest,direction,.false.,state,masked,scratch)
     if(any(t%izero/=mask)) error stop 'Output current changed recurrence map'
    end do
    t%izero=0
    input=0
    call t%kubo_ham_vec_matmul(input,masked,scratch,a=1.7_rp,b=0.13_rp)
    if(any(t%idum/=0).or.maxval(abs(masked))>0) error stop 'Empty H support'
    call t%kubo_operator_vec_matmul('current',qtest,direction,.true.,input,masked,scratch)
    if(any(t%idum/=0).or.maxval(abs(masked))>0) error stop 'Empty current support'
    input=psi; input(:,:,2:n)=0
   end do
  end do
 end subroutine check_masked_products
   subroutine legacy_ham_hoh_vec_matmul(this, psi_in, psi_out, a, b)
      class(recursion), intent(inout) :: this
      complex(rp), dimension(:, :, :), intent(in) :: psi_in
      complex(rp), dimension(:, :, :), intent(out) :: psi_out
      real(rp), intent(in) :: a, b
      ! Local variables
      integer :: nb, ih, j, k, nr, m, n, l, hblocksize, nat, nnmap, nlimplus1
      complex(rp), dimension(18, 18) :: dum1, dum2, locham

      hblocksize = 18
      nat = this%lattice%kk

      this%hohpsi(:, :, :) = (0.0d0, 0.0d0)
      this%enupsi(:, :, :) = (0.0d0, 0.0d0)
      this%socpsi(:, :, :) = (0.0d0, 0.0d0)

      ! Write H*|phi_1> for the local Hamiltonian
      nlimplus1 = this%lattice%nmax + 1
      if (this%lattice%nmax /= 0) then ! In case of impurities using the local hamiltonian
         !$omp parallel do default(shared) private(k, ih, nr, nb, nnmap,locham)
         do k = 1, this%lattice%nmax ! Loop in the neighbouring
            this%idum(k) = this%izero(k)
            ih = this%lattice%iz(k)
            nr = this%lattice%nn(k, 1) ! Number of neighbours of atom i
            if (this%izero(k) /= 0) then
               locham = this%hamiltonian%hall(1:18, 1:18, 1, k)
               call zgemm('n', 'n', 18, 18, 18, cone, locham, 18, psi_in(:, :, k), 18, cone, this%psi2(:, :, k), 18)
               call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%lsham(:, :, ih), 18, psi_in(:, :, k), 18, cone, this%socpsi(:, :, k), 18)
               call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%enim(:, :, ih), 18, psi_in(:, :, k), 18, cone, this%enupsi(:, :, k), 18)
            end if
            if (nr >= 2) then
               do nb = 2, nr ! Loop on the neighbouring
                  nnmap = this%lattice%nn(k, nb)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%hall(:, :, nb, k), 18, psi_in(:, :, nnmap), 18, cone, this%psi2(:, :, k), 18)
                        this%idum(k) = 1
                     end if
                  end if
               end do ! End of loop in the neighbouring
            end if
         end do ! End of loop in the neighbouring
         !$omp end parallel do
      end if ! End of local Hamiltonian loop

      ! Write H*|phi_1> for the bulk Hamiltonian
      !$omp parallel do default(shared) private(k, ih, nr, nb, nnmap,locham)
      do k = nlimplus1, this%lattice%kk ! Loop in the clust
         this%idum(k) = this%izero(k)
         ih = this%lattice%iz(k)
         nr = this%lattice%nn(k, 1)
         if (this%izero(k) /= 0) then
            locham = this%hamiltonian%ee(1:18, 1:18, 1, ih)
            call zgemm('n', 'n', 18, 18, 18, cone, locham, 18, psi_in(:, :, k), 18, cone, this%psi2(:, :, k), 18)
            call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%enim(:, :, ih), 18, psi_in(:, :, k), 18, cone, this%enupsi(:, :, k), 18)
            call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%lsham(:, :, ih), 18, psi_in(:, :, k), 18, cone, this%socpsi(:, :, k), 18)
         end if
         if (nr >= 2) then
            do nb = 2, nr ! Loop in the neighbouring
               nnmap = this%lattice%nn(k, nb)
               if (nnmap /= 0 .and. this%izero(nnmap) /= 0) then
                  call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%ee(:, :, nb, ih), 18, psi_in(:, :, nnmap), 18, cone, this%psi2(:, :, k), 18)
                  this%idum(k) = 1
               end if
            end do ! End of the loop in the neighbouring
         end if
      end do ! End loop in the clust
      !$omp end parallel do

      ! Mapping update for hoh calculation
      this%izero(:) = this%idum(:)

      if (this%lattice%nmax /= 0) then ! In case of impurities using the local hamiltonian
         !$omp parallel do default(shared) private(k, ih, nr, nb, nnmap,locham)
         do k = 1, this%lattice%nmax ! Loop in the neighbouring
            this%idum(k) = this%izero(k)
            ih = this%lattice%iz(k)
            nr = this%lattice%nn(k, 1) ! Number of neighbours of atom i
            if (this%izero(k) /= 0) then
               call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%hallo(:, :, 1, k), 18, this%psi2(:, :, k), 18, cone, this%hohpsi(:, :, k), 18)
            end if
            if (nr >= 2) then
               do nb = 2, nr ! Loop on the neighbouring
                  nnmap = this%lattice%nn(k, nb)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%hallo(:, :, nb, k), 18, this%psi2(:, :, nnmap), 18, cone, this%hohpsi(:, :, k), 18)
                        this%idum(k) = 1
                     end if
                  end if
               end do ! End of loop in the neighbouring
            end if
         end do ! End of loop in the neighbouring
         !$omp end parallel do
      end if ! End of local Hamiltonian loop

      !$omp parallel do default(shared) private(k, ih, nr, nb, nnmap,locham)
      do k = nlimplus1, this%lattice%kk ! Loop to find the bulk atoms using the bulk Hamiltonian
         this%idum(k) = this%izero(k)
         ih = this%lattice%iz(k) ! Atom type
         nr = this%lattice%nn(k, 1) ! Number of neighbours
         if (this%izero(k) /= 0) then
            call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%eeo(:, :, 1, ih), 18, this%psi2(:, :, k), 18, cone, this%hohpsi(:, :, k), 18)
         end if
         if (nr >= 2) then
            do nb = 2, nr ! Loop on the neighbouring
               nnmap = this%lattice%nn(k, nb)
               if (nnmap /= 0) then
                  if (this%izero(nnmap) /= 0) then
                     call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%eeo(:, :, nb, ih), 18, this%psi2(:, :, nnmap), 18, cone, this%hohpsi(:, :, k), 18)
                     this%idum(k) = 1
                  end if
               end if
            end do ! End of loop in the neighbouring
         end if
      end do
      !$omp end parallel do

      ! H = h - hoh + e_nu + l.s
      psi_out(:, :, :) = this%psi2(:, :, :) - this%hohpsi(:, :, :) + this%enupsi(:, :, :) + this%socpsi(:, :, :)

      ! Do the scaling and shifting
      psi_out(:, :, :) = psi_out(:, :, :) - b*psi_in(:, :, :)
      psi_out(:, :, :) = psi_out(:, :, :)/a

      this%psi2(:, :, :) = (0.0_rp, 0.0_rp)
   end subroutine legacy_ham_hoh_vec_matmul
   subroutine legacy_ham_vec_matmul(this, psi_in, psi_out, a, b)
      class(recursion), intent(inout) :: this
      complex(rp), dimension(:, :, :), intent(in) :: psi_in
      complex(rp), dimension(:, :, :), intent(out) :: psi_out
      real(rp), intent(in) :: a, b
      ! Local variables
      integer :: nb, ih, j, k, nr, m, n, l, hblocksize, nat, nnmap, nlimplus1
      complex(rp), dimension(18, 18) :: dum1, dum2, locham
      hblocksize = 18
      nat = this%lattice%kk

      psi_out(:, :, :) = (0.0d0, 0.0d0)

      ! Write H*|phi_1> for the local Hamiltonian
      nlimplus1 = this%lattice%nmax + 1
      if (this%lattice%nmax /= 0) then ! In case of impurities using the local hamiltonian
         !$omp parallel do default(shared) private(k, ih, nr, nb, nnmap,locham)
         do k = 1, this%lattice%nmax ! Loop in the neighbouring
            this%idum(k) = this%izero(k)
            ih = this%lattice%iz(k)
            nr = this%lattice%nn(k, 1) ! Number of neighbours of atom i
            if (this%izero(k) /= 0) then
               locham = this%hamiltonian%hall(1:18, 1:18, 1, k) + this%hamiltonian%lsham(1:18, 1:18, ih)
               call zgemm('n', 'n', 18, 18, 18, cone, locham, 18, psi_in(:, :, k), 18, cone, psi_out(:, :, k), 18)
            end if
            if (nr >= 2) then
               do nb = 2, nr ! Loop on the neighbouring
                  nnmap = this%lattice%nn(k, nb)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%hall(:, :, nb, k), 18, psi_in(:, :, nnmap), 18, cone, psi_out(:, :, k), 18)
                        this%idum(k) = 1
                     end if
                  end if
               end do ! End of loop in the neighbouring
            end if
         end do ! End of loop in the neighbouring
         !$omp end parallel do
      end if ! End of local Hamiltonian loop

      ! Write H*|phi_1> for the bulk Hamiltonian
      !$omp parallel do default(shared) private(k, ih, nr, nb, nnmap,locham)
      do k = nlimplus1, this%lattice%kk ! Loop in the clust
         this%idum(k) = this%izero(k)
         ih = this%lattice%iz(k)
         nr = this%lattice%nn(k, 1)
         if (this%izero(k) /= 0) then
            locham = this%hamiltonian%ee(1:18, 1:18, 1, ih) + this%hamiltonian%lsham(1:18, 1:18, ih)
            call zgemm('n', 'n', 18, 18, 18, cone, locham, 18, psi_in(:, :, k), 18, cone, psi_out(:, :, k), 18)
         end if
         if (nr >= 2) then
            do nb = 2, nr ! Loop in the neighbouring
               nnmap = this%lattice%nn(k, nb)
               if (nnmap /= 0 .and. this%izero(nnmap) /= 0) then
                  call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%ee(:, :, nb, ih), 18, psi_in(:, :, nnmap), 18, cone, psi_out(:, :, k), 18)
                  this%idum(k) = 1
               end if
            end do ! End of the loop in the neighbouring
         end if
      end do ! End loop in the clust
      !$omp end parallel do
      ! Do the scaling and shifting
      psi_out(:, :, :) = psi_out(:, :, :) - b*psi_in(:, :, :)
      psi_out(:, :, :) = psi_out(:, :, :)/a
   end subroutine legacy_ham_vec_matmul

 ! Compare the production H kernels directly with the archived legacy routines.
 ! Compare V with the archived velocity kernels below; the HOH difference must
 ! equal exactly the missing -h o v term, not zero.
 subroutine check_legacy_products()
  complex(rp) :: input(b,b,n), old(b,b,n), current(b,b,n), scratch(b,b,n,6)
  complex(rp) :: vp(b,b,4,n), vop(b,b,4,n), missing(b,b,n), tmp(b,b,n), tmp2(b,b,n)
  real(rp), parameter :: a_scale=1.7_rp, b_shift=0.13_rp, tolerance=5.e-13_rp
  real(rp) :: error, correction_size, distance
  integer :: i,j,m,it,im,ih,localized,repeat,iw,report
  allocate(t%izero(n),t%psi1(b,b,n),t%psi2(b,b,n))
  allocate(t%hohpsi(b,b,n),t%enupsi(b,b,n),t%socpsi(b,b,n))
  allocate(ham%eeo(b,b,4,n),ham%hallo(b,b,4,n))
  ! The old routines evaluate izero(nnmap) in a compound .and. expression.
  ! Use compact neighbor lists for their valid-domain comparison. New kernels
  ! are separately checked with zero neighbor entries by the dense tests.
  lat%nn(:,1)=[2,3,2]
  ham%eeo=0; ham%hallo=0
  do i=1,n
   it=lat%iz(i)
   do m=1,lat%nn(i,1)
    j=i
    if(m>1) j=lat%nn(i,m)
    ham%eeo(:,:,m,it)=matmul(h(:,:,m,it),o(:,:,lat%iz(j)))
    ham%hallo(:,:,m,i)=matmul(hl(:,:,m,i),o(:,:,lat%iz(j)))
   end do
  end do
  open(newunit=report,file='legacy_comparison.csv',status='replace')
  write(report,'(a)') 'kind,nmax,localized,repeat,hoh,weighted,max_abs_error,missing_term_size'
  do im=0,1
   lat%nmax=2*im
   do ih=0,1
    ham%hoh=ih==1
    do localized=0,1
     do repeat=1,2
      do i=1,n
       do j=1,b
        do m=1,b
         input(m,j,i)=cmplx(sin(real(m+3*j+7*i+repeat,rp)),cos(real(2*m+j+i+repeat,rp)),rp)
        end do
       end do
      end do
      if(localized==1) input(:,:,2:n)=0
      t%izero=1
      if(localized==1) t%izero(2:n)=0
      t%idum=0; t%psi2=0
      if(ham%hoh) then
       call legacy_ham_hoh_vec_matmul(t,input,old,a_scale,b_shift)
      else
       call legacy_ham_vec_matmul(t,input,old,a_scale,b_shift)
      end if
      t%izero=1
      if(localized==1) t%izero(2:n)=0
      call t%kubo_ham_vec_matmul(input,current,scratch,a=a_scale,b=b_shift)
      error=maxval(abs(current-old))
      write(report,'(a,5(",",i0),2(",",es24.16))') 'H',lat%nmax,localized,repeat,ih,0, error,0.0_rp
      if(error>tolerance) error stop 'Legacy H and new H differ after matching scale/shift'
      if(ham%hoh) then
       ! Original orbital_modern call signature uses existing object scratch.
       t%izero=1
       if(localized==1) t%izero(2:n)=0
       call t%ham_hoh_vec_matmul(input,current,a_scale,b_shift)
       if(maxval(abs(current-old))>tolerance) error stop 'Object scratch HOH differs'
      end if
      ! The old velocity has no impurity rows: compare bulk only.
      if(im/=0) cycle
      do iw=0,1
       vp=0; vop=0
       do i=1,n
        it=lat%iz(i)
        do m=2,lat%nn(i,1)
         j=lat%nn(i,m)
         distance=dot_product(direction,cr(:,i)-cr(:,j))*lat%alat
         if(iw==1) distance=distance*max(scale(it),scale(lat%iz(j)))
         vp(:,:,m,it)=cmplx(0.0_rp,-distance,rp)*h(:,:,m,it)
         vop(:,:,m,it)=matmul(vp(:,:,m,it),o(:,:,lat%iz(j)))
        end do
       end do
       t%izero=1
       if(localized==1) t%izero(2:n)=0
       t%idum=0; t%psi1=0; t%psi2=0
       if(ham%hoh) then
        call legacy_velo_hoh_vec_matmul(t,vp,vop,input,old)
       else
        call legacy_velo_vec_matmul(t,'n',vp,input,old)
       end if
       t%izero=1
       if(localized==1) t%izero(2:n)=0
       if(ham%hoh) then
        call t%velo_hoh_vec_matmul(input,current,scratch,direction,iw==1)
       else
        call t%velo_vec_matmul(input,current,direction,iw==1)
       end if
       missing=0
       if(ham%hoh) then
        call t%velo_vec_matmul(input,tmp,direction,iw==1)
        t%izero=t%idum
        call t%onsite_vec_matmul(ham%obarm,tmp,tmp2)
        call t%ham_vec_matmul(tmp2,missing,1.0_rp,0.0_rp,.false.)
       end if
       error=maxval(abs(current-old+missing))
       correction_size=maxval(abs(missing))
       write(report,'(a,5(",",i0),2(",",es24.16))') 'V',lat%nmax,localized,repeat,ih,iw, error,correction_size
       if(error>tolerance) error stop 'New velocity does not equal old velocity minus h o v'
      end do
     end do
    end do
   end do
  end do
  close(report)
  lat%nn=nn
 end subroutine check_legacy_products

! Test-only archived velocity routines from e604f4029e59d35b6e53551dd2d95919a9baf8e0.
! Active statements and OpenMP directives are unchanged, except procedure names.
! Comment-only lines are omitted. These routines are never linked into the solver.
   subroutine legacy_velo_vec_matmul(this, c_or_n, v_op, psi_in, psi_out)
      class(recursion), intent(inout) :: this
      complex(rp), dimension(:, :, :, :), intent(in) :: v_op
      complex(rp), dimension(:, :, :), intent(in) :: psi_in
      character :: c_or_n
      complex(rp), dimension(:, :, :), intent(out) :: psi_out
      integer :: nb, ih, j, k, nr, m, n, l, hblocksize, nat, nnmap, nlimplus1
      complex(rp), dimension(18, 18) :: dum1, dum2, locham
      hblocksize = 18
      nat = this%lattice%kk

      psi_out(:, :, :) = (0.0d0, 0.0d0)

      nlimplus1 = this%lattice%nmax + 1

      !$omp parallel do default(shared) private(k, ih, nr, nb, nnmap)
      do k = nlimplus1, this%lattice%kk ! Loop in the clust
         this%idum(k) = this%izero(k)
         ih = this%lattice%iz(k)
         nr = this%lattice%nn(k, 1)
         if (this%izero(k) /= 0) then
            call zgemm('n', 'n', 18, 18, 18, cone, v_op(:, :, 1, ih), 18, psi_in(:, :, k), 18, cone, psi_out(:, :, k), 18)
         end if
         if (nr >= 2) then
            do nb = 2, nr ! Loop in the neighbouring
               nnmap = this%lattice%nn(k, nb)
               if (nnmap /= 0 .and. this%izero(nnmap) /= 0) then
                  call zgemm(c_or_n, 'n', 18, 18, 18, cone, v_op(:, :, nb, ih), 18, psi_in(:, :, nnmap), 18, cone, psi_out(:, :, k), 18)
                  this%idum(k) = 1
               end if
            end do ! End of the loop in the neighbouring
         end if
      end do ! End loop in the clust
      !$omp end parallel do
   end subroutine legacy_velo_vec_matmul

   subroutine legacy_velo_hoh_vec_matmul(this, v_op, vo_op, psi_in, psi_out)
      class(recursion), intent(inout) :: this
      complex(rp), dimension(:, :, :), intent(in) :: psi_in
      complex(rp), dimension(:, :, :), intent(out) :: psi_out
      complex(rp), dimension(:, :, :, :), intent(in) :: v_op
      complex(rp), dimension(:, :, :, :), intent(in) :: vo_op
      integer :: nb, ih, j, k, nr, m, n, l, hblocksize, nat, nnmap, nlimplus1
      complex(rp), dimension(18, 18) :: dum1, dum2, locham

      hblocksize = 18
      nat = this%lattice%kk

      this%hohpsi(:, :, :) = (0.0d0, 0.0d0)
      this%enupsi(:, :, :) = (0.0d0, 0.0d0)
      this%socpsi(:, :, :) = (0.0d0, 0.0d0)

      nlimplus1 = this%lattice%nmax + 1

      !$omp parallel do default(shared) private(k, ih, nr, nb, nnmap,locham)
      do k = nlimplus1, this%lattice%kk ! Loop in the clust
         this%idum(k) = this%izero(k)
         ih = this%lattice%iz(k)
         nr = this%lattice%nn(k, 1)
         if (this%izero(k) /= 0) then
            call zgemm('n', 'n', 18, 18, 18, cone, v_op(:, :, 1, ih), 18, psi_in(:, :, k), 18, cone, this%psi2(:, :, k), 18)
            call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%ee(:, :, 1, ih), 18, psi_in(:, :, k), 18, cone, this%psi1(:, :, k), 18)
         end if
         if (nr >= 2) then
            do nb = 2, nr ! Loop in the neighbouring
               nnmap = this%lattice%nn(k, nb)
               if (nnmap /= 0 .and. this%izero(nnmap) /= 0) then
                  call zgemm('n', 'n', 18, 18, 18, cone, v_op(:, :, nb, ih), 18, psi_in(:, :, nnmap), 18, cone, this%psi2(:, :, k), 18)
                  call zgemm('n', 'n', 18, 18, 18, cone, this%hamiltonian%ee(:, :, nb, ih), 18, psi_in(:, :, nnmap), 18, cone, this%psi1(:, :, k), 18)
                  this%idum(k) = 1
               end if
            end do ! End of the loop in the neighbouring
         end if
      end do ! End loop in the clust
      !$omp end parallel do

      this%izero(:) = this%idum(:)


      !$omp parallel do default(shared) private(k, ih, nr, nb, nnmap,locham)
      do k = nlimplus1, this%lattice%kk ! Loop to find the bulk atoms using the bulk Hamiltonian
         this%idum(k) = this%izero(k)
         ih = this%lattice%iz(k) ! Atom type
         nr = this%lattice%nn(k, 1) ! Number of neighbours
         if (this%izero(k) /= 0) then
         end if
         if (nr >= 2) then
            do nb = 2, nr ! Loop on the neighbouring
               nnmap = this%lattice%nn(k, nb)
               if (nnmap /= 0) then
                  if (this%izero(nnmap) /= 0) then
                     call zgemm('n', 'n', 18, 18, 18, cone, vo_op(:, :, nb, ih), 18, this%psi1(:, :, nnmap), 18, cone, this%hohpsi(:, :, k), 18)
                     this%idum(k) = 1
                  end if
               end if
            end do ! End of loop in the neighbouring
         end if
      end do
      !$omp end parallel do

      psi_out(:, :, :) = this%psi2(:, :, :) - this%hohpsi(:, :, :) + this%enupsi(:, :, :) + this%socpsi(:, :, :)

      this%psi2(:, :, :) = (0.0_rp, 0.0_rp)
      this%psi1(:, :, :) = (0.0_rp, 0.0_rp)

   end subroutine legacy_velo_hoh_vec_matmul
end program
