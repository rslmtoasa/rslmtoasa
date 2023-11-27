!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Green
!
!> @author
!> Angela Klautau
!> Ramon Cardias
!> Lucas P. Campagna
!> S. Frota-Pessôa
!> Pascoal R. Peduto
!> Anders Bergman
!> S. B. Legoas
!> H. M. Petrilli
!> Ivan P. Miranda
!
! DESCRIPTION:
!> Module to handle the calculation of the density of states
!------------------------------------------------------------------------------

module density_of_states_mod

   use control_mod
   use lattice_mod
   use energy_mod
   use symbolic_atom_mod
   use recursion_mod
   use precision_mod, only: rp
   use math_mod
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc
#endif
   implicit none

   private

   !> Module's main structure
   type, public :: dos
      ! General variables
      real(rp), dimension(:, :, :), allocatable :: doscheb

      !> Recursion
      class(recursion), pointer :: recursion
      !> Symbolic atom
      class(symbolic_atom), dimension(:), pointer :: symbolic_atom
      !> Lattice
      class(lattice), pointer :: lattice
      !> Control
      class(control), pointer :: control
      !> Energy
      class(energy), pointer :: en
   contains
      procedure :: density
      procedure :: chebyshev_dos
      procedure :: chebyshev_dos_full
      procedure :: restore_to_default
      final :: destructor
   end type

   interface dos
      procedure :: constructor
   end interface dos

contains

   function constructor(recursion_obj, energy_obj) result(obj)
      type(dos) :: obj
      class(recursion), target, intent(in) :: recursion_obj
      class(energy), target, intent(in) :: energy_obj

      obj%recursion => recursion_obj
      obj%en => energy_obj
      obj%symbolic_atom => recursion_obj%hamiltonian%charge%lattice%symbolic_atoms
      obj%lattice => recursion_obj%lattice
      obj%control => recursion_obj%lattice%control

      call obj%restore_to_default()
   end function constructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Destructor
   !---------------------------------------------------------------------------
   subroutine destructor(this)
      type(dos) :: this
#ifdef USE_SAFE_ALLOC
      if (allocated(this%doscheb)) call g_safe_alloc%deallocate('density_of_states.doscheb', this%doscheb)
#else
      if (allocated(this%doscheb)) deallocate (this%doscheb)
#endif
   end subroutine destructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Uses the moments form the Chebyshev recursions to calculate the DOS
   !---------------------------------------------------------------------------
   subroutine chebyshev_dos(this)
      class(dos), intent(inout) :: this
      ! Local variables
      real(rp), dimension(this%control%lld*2 + 2) :: kernel
      real(rp), dimension(this%en%channels_ldos + 10, 0:this%control%lld*2 + 2) :: polycheb
      real(rp), dimension(this%en%channels_ldos + 10) :: w, wscale
      real(rp) :: eps, a, b
      integer :: i, l, n
      eps = 0.0001D0
      ! Defining rescaling coeficients
      a = (this%en%energy_max - this%en%energy_min) / (2 - 0.3)
      b = (this%en%energy_max + this%en%energy_min) / 2

      wscale(:) = (this%en%ene(:) - b) / a

      ! Calculating the Jackson Kernel
      call jackson_kernel((this%control%lld) * 2 + 2, kernel)

      ! Calculating the Lorentz Kernel
!    call lorentz_kernel(this%control%lld, kernel, 4.0d0)

      do n = 1, this%lattice%nrec ! Loop on self-consistent atoms
         ! Multiply the moments with the kernel
         do i = 1, 18
            this%recursion%mu_ng(i, i, :, n) = this%recursion%mu_n(i, i, :, n) * kernel(:)
         end do

         this%recursion%mu_ng(:, :, 2:size(kernel), n) = this%recursion%mu_ng(:, :, 2:size(kernel), n) * 2.0_RP

         do i = 1, size(kernel)
            write (400 + n, *) i, sum(this%recursion%mu_n(:, :, i, n))
         end do

         ! Calculate the Chebyshev polynomials
         call t_polynomial(size(w), size(kernel), wscale(:), polycheb)

         ! Calculate the density of states
         do i = 1, size(kernel)
            do l = 1, 18
               this%doscheb(l, :, n) = this%doscheb(l, :, n) + this%recursion%mu_ng(l, l, i, n) * polycheb(:, i - 1)
               !do m=1, 18
               !  green(n, l, m, :) = green(n, l, m, :) + (-i_unit*(this%recursion%mu_ng(n, i, l, m)*exp(-i_unit*(i-1)*acos(wscale(:)))))
               !end do
            end do
         end do
         do l = 1, 18
            this%doscheb(l, :, n) = this%doscheb(l, :, n) / ((sqrt((a**2) - ((this%en%ene(:) - b)**2))) * pi)
            !do m=1, 18
            !green(n, l, m, :) = green(n, l, m, :)/((sqrt((a**2)-((this%self%ene(:)-b)**2))))
            !end do
         end do
         do l = 1, 18
            do i = 1, this%en%channels_ldos + 10
               if (isnan(this%doscheb(l, i, n))) this%doscheb(l, i, n) = 0.0D0
            end do
         end do

         do i = 1, this%en%channels_ldos + 10
            write (200 + n, '(8f16.4)') this%en%ene(i), (this%doscheb(1, i, n)), sum(this%doscheb(2:4, i, n)), sum(this%doscheb(5:9, i, n)), &
               (this%doscheb(10, i, n)), sum(this%doscheb(11:13, i, n)), sum(this%doscheb(14:18, i, n)), &
               sum(this%doscheb(1:18, i, n))
            !write(125+n, '(10f10.6)') this%self%ene(i), real(green(n, 1, 1, i)), real(green(n, 2, 2, i)), real(green(n, 3, 3, i)), real(green(n, 4, 4, i)), real(green(n, 5, 5, i)), &
            !                                           &real(green(n, 6, 6, i)), real(green(n, 7, 7, i)), real(green(n, 8, 8, i)), real(green(n, 9, 9, i))
         end do
      end do  ! End loop on self-consistent atoms
   end subroutine chebyshev_dos

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Uses the moments form the Chebyshev recursions to calculate the DOS
   !---------------------------------------------------------------------------
   subroutine chebyshev_dos_full(this)
      class(dos), intent(inout) :: this
      ! Local variables
      real(rp), dimension(this%control%lld) :: kernel
      real(rp), dimension(this%en%channels_ldos + 10, 0:this%control%lld) :: polycheb
      real(rp), dimension(this%en%channels_ldos + 10) :: w, wscale
      real(rp) :: eps, a, b
      
      complex(rp), dimension(this%lattice%nrec, 18, 18, this%en%channels_ldos + 10) :: green
      integer :: i, l, m, n
      eps = 0.0001D0

      ! Defining rescaling coeficients
      a = (this%en%energy_max - this%en%energy_min) / (2 - 0.3)
      b = (this%en%energy_max + this%en%energy_min) / 2

      wscale(:) = (this%en%ene(:) - b) / a

      ! Calculating the Jackson Kernel
      call jackson_kernel(this%control%lld, kernel)

      ! Calculating the Lorentz Kernel
!    call lorentz_kernel(this%control%lld, kernel, 4.0d0)

      do n = 1, this%lattice%nrec ! Loop on self-consistent atoms
         ! Multiply the moments with the kernel
         do i = 1, 18
            this%recursion%mu_ng(i, i, :, n) = this%recursion%mu_n(i, i, :, n) * kernel(:)
         end do

         this%recursion%mu_ng(:, :, 2:size(kernel), n) = this%recursion%mu_ng(:, :, 2:size(kernel), n) * 2.0_RP

         do i = 1, size(kernel)
            write (400 + n, *) i, this%recursion%mu_n(5, 6, i, n)
         end do

         ! Calculate the Chebyshev polynomials
         call t_polynomial(size(w), size(kernel), wscale(:), polycheb)

         ! Calculate the density of states
         do i = 1, size(kernel)
            do l = 1, 18
               this%doscheb(l, :, n) = this%doscheb(l, :, n) + this%recursion%mu_ng(l, l, i, n) * polycheb(:, i - 1)
               do m = 1, 18
                  green(n, l, m, :) = green(n, l, m, :) + (-i_unit * (this%recursion%mu_ng(l, m, i, n) * exp(-i_unit * (i - 1) * acos(wscale(:)))))
               end do
            end do
         end do
         do l = 1, 18
            this%doscheb(l, :, n) = this%doscheb(l, :, n) / ((sqrt((a**2) - ((this%en%ene(:) - b)**2))) * pi)
            do m = 1, 18
               green(n, l, m, :) = green(n, l, m, :) / ((sqrt((a**2) - ((this%en%ene(:) - b)**2))))
            end do
         end do
         do l = 1, 18
            do i = 1, this%en%channels_ldos + 10
               if (isnan(this%doscheb(l, i, n))) this%doscheb(l, i, n) = 0.0D0
            end do
         end do

         do i = 1, this%en%channels_ldos + 10
            write (200 + n, '(8f16.4)') this%en%ene(i), (this%doscheb(1, i, n)), sum(this%doscheb(2:4, i, n)), sum(this%doscheb(5:9, i, n)), &
               (this%doscheb(10, i, n)), sum(this%doscheb(11:13, i, n)), sum(this%doscheb(14:18, i, n)), &
               sum(this%doscheb(1:18, i, n))
            write (125 + n, '(10f10.6)') this%en%ene(i), real(green(n, 1, 1, i)), real(green(n, 2, 2, i)), real(green(n, 3, 3, i)), real(green(n, 4, 4, i)), real(green(n, 5, 5, i)), &
                                                       &real(green(n, 5, 5, i)), real(green(n, 6, 6, i)), real(green(n, 8, 8, i)), aimag(green(n, 9, 9, i))
         end do
      end do  ! End loop on self-consistent atoms
   end subroutine chebyshev_dos_full

   ! Member functions
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the density of states
   !---------------------------------------------------------------------------
   subroutine density(this, tdens, ia, mdir)
      use mpi_mod
      class(dos) :: this
      ! Input
      integer, intent(in) :: mdir ! Direction index
      integer, intent(in) :: ia ! Atom type
      ! Output
      real(rp), dimension(18, this%en%channels_ldos + 10), intent(out) :: tdens
      ! Local variables
      integer :: k, l, ll1, nb, nbp1, nl, eidx, ll_in, ifail, npts, nw
      real(rp) :: a1, a2, dens, emax, emin, eps, err, e_shift, prefac
      

      integer, dimension(18) :: nb2
      
      
      integer, dimension(18) :: ll_fail
      real(rp), dimension(this%control%lld) :: aa, am, bb, bm, sqbb
      real(rp), dimension(10) :: edge, weight, width
      
      real(rp), dimension(18, this%control%lld) :: am2, bm2
      real(rp), dimension(18, 10) :: edge2, width2, weight2
      
      real(rp), dimension(2) :: bandedges
      real(rp), dimension(18) :: alpha_inf, beta_inf

      integer :: ia_glob

      eps = 1.0D-14
      err = 0.00001D0
      nbp1 = 2
      ll_in = this%control%lld
      ll_fail = ll_in
      nw = 10 * this%lattice%ntype * this%control%lld
      npts = this%en%channels_ldos + 10
      tdens = 0.0D0
      ia_glob = l2g_map(ia)

      do nl = 1, 18
         do l = 1, this%control%lld
            aa(l) = this%recursion%a(l, nl, ia, mdir)
            bb(l) = this%recursion%b2(l, nl, ia, mdir)!**2
            sqbb(l) = sqrt(bb(l))
         end do
         if (nl == 1 .or. nl == 10) bb = 1.025D0 * bb
         ll1 = this%control%lld
         am = 0.0D0; bm = 0.0D0; edge = 0.0D0; width = 0.0D0; weight = 0.0D0
         err = 0.0000001D0
         call this%recursion%bpOPT(this%control%lld, AA, sqBB, this%control%lld - 1, am(1), bm(1), ifail)
         if (nl == 1 .or. nl == 10) bm = 1.01D0 * bm
         alpha_inf(nl) = am(1)
         beta_inf(nl) = bm(1)
         nb = 1
         edge(1) = am(1) - 2.0D0 * bm(1)
         width(1) = 4.0D0 * bm(1)
         weight(1) = 1.0D0
         nb2(nl) = nb
         am = AA
         bm = BB
         do k = 1, nb
            a1 = edge(k)
            a2 = edge(k) + width(k)
            edge2(nl, k) = edge(k)
            width2(nl, k) = width(k)
            weight2(nl, k) = weight(k)
         end do
         if (nl == 1) then
            emin = a1
            emax = a2
         else
            emin = min(emin, a1)
            emax = max(emax, a2)
         end if
         if (nb > 0) then
            do l = 1, this%control%lld
               am2(nl, l) = am(l)
               bm2(nl, l) = bm(l)
            end do
         else
            do l = 1, this%control%lld
               am2(nl, l) = 0.0D0
               bm2(nl, l) = 0.0D0
            end do
         end if
      end do

      eps = 1.0D-10
      tdens = 0.0D0

      do eidx = 1, npts
         do nl = 1, 18
            nb = nb2(nl)
            if (nb > 0) then
               do l = 1, this%control%lld
                  aa(l) = this%recursion%a(l, nl, ia, mdir)
                  bb(l) = this%recursion%b2(l, nl, ia, mdir) !*0.5d0
                  am(l) = am2(nl, l)
                  bm(l) = bm2(nl, l) !*0.5d0
               end do
               edge = 0.0D0; width = 0.0D0; weight = 0.0D0
               do k = 1, nb
                  edge(k) = edge2(nl, k)
                  width(k) = width2(nl, k)
                  weight(k) = weight2(nl, k)
               end do
               e_shift = this%en%ene(eidx) / this%symbolic_atom(ia_glob)%potential%dw_l(nl) &
                         - 1.00D0 * this%symbolic_atom(ia_glob)%potential%cshi(nl)
               prefac = 1.0D0
               bandedges(1) = edge(1)!-0.1
               bandedges(2) = edge(1) + width(1)!+0.1
               ! only terminator 5 implemented
               dens = bprldos(e_shift, AA, BB, this%control%lld, bandedges)
            else
               dens = 0.0D0
            end if
            tdens(nl, eidx) = tdens(nl, eidx) + prefac * dens / this%symbolic_atom(ia_glob)%potential%dw_l(nl)
         end do
      end do

      do eidx = 1, npts
         write (300, *) this%en%ene(eidx), sum(tdens(1:9, eidx)), sum(tdens(10:18, eidx))
      end do
   end subroutine density

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the density of states
   !---------------------------------------------------------------------------
   real(rp) function bprldos(e, a, b2, ll, ei)
      ! Input
      integer, intent(in) :: LL
      real(rp), intent(in) :: E
      real(rp), intent(in) :: A(LL), B2(LL)
      real(rp), intent(in) :: ei(2)
      real(rp) :: DISC, BI2, AIT, RT
      complex(rp) :: etop, ebot, ea, eb, emid
      complex(rp) :: det, zoff, Qt
      integer :: l

      ebot = ei(1) !a(ll)-2*b2(ll)
      etop = ei(2) !a(ll)+2*b2(ll)
      !     ebot=a(ll)-2*b2(ll)
      !     etop=a(ll)+2*b2(ll)
      emid = 0.5D0 * (etop + ebot)
      ea = e - etop
      eb = e - ebot
      det = ea * eb
!    zoff = (0.0d0, 0.01d0)
      zoff = sqrt(det)
      Qt = (e - emid - zoff) * 0.5D0
      if (aimag(Qt) > 0.D0) then
         Qt = (e - emid + zoff) * 0.5D0
      end if
      do l = ll - 1, 1, -1
         Qt = B2(l) / (e - A(l) - Qt)
      end do
      bpRLDOS = -aimag(Qt) / PI
      return
   end function bprldos

   subroutine restore_to_default(this)
      class(dos), intent(inout) :: this
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('density_of_states.doscheb', this%doscheb, (/18, this%en%channels_ldos + 10, this%lattice%nrec/))
#else
      allocate (this%doscheb(18, this%en%channels_ldos + 10, this%lattice%nrec))
#endif

      this%doscheb(:, :, :) = 0.0D0
   end subroutine

end module density_of_states_mod
