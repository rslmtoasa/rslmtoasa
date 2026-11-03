! Matrix-free transport actions. No recursion support mask is used: every
! intermediate includes all sites, including two-hop paths and local hall rows.
module transport_operators_mod
   use precision_mod, only: rp
   implicit none
   private
   public :: transport_operators
   type :: transport_operators
      integer, pointer :: nn(:, :) => null(), iz(:) => null()
      real(rp), pointer :: cr(:, :) => null(), scale(:) => null()
      complex(rp), pointer :: h(:, :, :, :) => null(), h_local(:, :, :, :) => null()
      complex(rp), pointer :: overlap(:, :, :) => null(), enu(:, :, :) => null(), soc(:, :, :) => null()
      integer :: nmax = 0
      real(rp) :: alat = 1.0_rp
      logical :: hoh = .false.
   contains
      procedure :: hop
      procedure :: local_action
      procedure :: hamiltonian
      procedure :: velocity
      procedure :: response
   end type
contains
   ! Optional direction turns h into its bare velocity. Scale weights only
   ! the input slot, preserving the existing symmetric bond max(s_i,s_j).
   subroutine hop(this, psi, out, direction, weighted)
      class(transport_operators), intent(in) :: this
      complex(rp), intent(in) :: psi(:, :, :)
      complex(rp), intent(out) :: out(:, :, :)
      real(rp), intent(in), optional :: direction(3)
      logical, intent(in), optional :: weighted
      integer :: i, j, m, it
      complex(rp) :: block(size(psi, 1), size(psi, 1)), factor
      real(rp) :: distance
      out = (0.0_rp, 0.0_rp)
      !$omp parallel do default(shared) private(i,j,m,it,block,factor,distance)
      do i = 1, size(psi, 3)
         it = this%iz(i)
         do m = 1, this%nn(i, 1)
            j = i
            if (m > 1) j = this%nn(i, m)
            if (j == 0) cycle
            if (i <= this%nmax) then
               block = this%h_local(:, :, m, i)
            else
               block = this%h(:, :, m, it)
            end if
            if (present(direction)) then
               distance = dot_product(direction, this%cr(:, i)-this%cr(:, j))*this%alat
               factor = cmplx(0.0_rp, -distance, rp)
               if (present(weighted)) then
                  if (weighted) factor = factor*max(this%scale(it), this%scale(this%iz(j)))
               end if
               block = factor*block
            end if
            out(:, :, i) = out(:, :, i) + matmul(block, psi(:, :, j))
         end do
      end do
      !$omp end parallel do
   end subroutine

   subroutine local_action(this, matrix, psi, out)
      class(transport_operators), intent(in) :: this
      complex(rp), intent(in) :: matrix(:, :, :), psi(:, :, :)
      complex(rp), intent(out) :: out(:, :, :)
      integer :: i
      !$omp parallel do default(shared) private(i)
      do i = 1, size(psi, 3)
         out(:, :, i) = matmul(matrix(:, :, this%iz(i)), psi(:, :, i))
      end do
      !$omp end parallel do
   end subroutine

   ! H = h - h o h + E_nu + H_SOC. E_nu belongs only to the HOH representation.
   subroutine hamiltonian(this, psi, out, include_soc)
      class(transport_operators), intent(in) :: this
      complex(rp), intent(in) :: psi(:, :, :)
      complex(rp), intent(out) :: out(:, :, :)
      logical, intent(in), optional :: include_soc
      complex(rp), allocatable :: tmp(:, :, :), tmp2(:, :, :)
      logical :: with_soc
      with_soc = .true.
      if (present(include_soc)) with_soc = include_soc
      allocate(tmp, mold=psi)
      call this%hop(psi, out)
      if (this%hoh) then
         allocate(tmp2, mold=psi)
         call this%local_action(this%overlap, out, tmp)
         call this%hop(tmp, tmp2)
         out = out - tmp2
         call this%local_action(this%enu, psi, tmp)
         out = out + tmp
      end if
      if (with_soc) then
         call this%local_action(this%soc, psi, tmp)
         out = out + tmp
      end if
   end subroutine

   ! Full product rule V = v - v o h - h o v, without assuming [o,h]=0.
   subroutine velocity(this, psi, out, direction, weighted)
      class(transport_operators), intent(in) :: this
      complex(rp), intent(in) :: psi(:, :, :)
      complex(rp), intent(out) :: out(:, :, :)
      real(rp), intent(in) :: direction(3)
      logical, intent(in) :: weighted
      complex(rp), allocatable :: tmp(:, :, :), tmp2(:, :, :), correction(:, :, :)
      call this%hop(psi, out, direction, weighted)
      if (.not. this%hoh) return
      allocate(tmp, mold=psi)
      allocate(tmp2, mold=psi)
      allocate(correction, mold=psi)
      call this%local_action(this%overlap, out, tmp)
      call this%hop(tmp, correction)
      out = out - correction
      call this%hop(psi, tmp)
      call this%local_action(this%overlap, tmp, tmp2)
      call this%hop(tmp2, correction, direction, weighted)
      out = out - correction
   end subroutine

   ! q is the local S/L/Q/O matrix, in the same basis as h.
   ! spin_torque = (1/i)[S,H_without_SOC]; spin-independent terms commute
   ! with S, so this includes all spin-dependent HOH and E_nu contributions.
   subroutine response(this, kind, q, direction, weighted, psi, out)
      class(transport_operators), intent(in) :: this
      character(len=*), intent(in) :: kind
      complex(rp), intent(in) :: q(:, :), psi(:, :, :)
      real(rp), intent(in) :: direction(3)
      logical, intent(in) :: weighted
      complex(rp), intent(out) :: out(:, :, :)
      complex(rp), allocatable :: qpsi(:, :, :), tmp(:, :, :), other(:, :, :)
      integer :: i
      if (kind == 'charge') then
         call this%velocity(psi, out, direction, weighted)
         return
      end if
      allocate(qpsi, mold=psi)
      do i = 1, size(psi, 3)
         qpsi(:, :, i) = matmul(q, psi(:, :, i))
      end do
      if (kind == 'local') then
         out = qpsi
         return
      end if
      allocate(tmp, mold=psi)
      allocate(other, mold=psi)
      select case(kind)
      case('current')
         call this%velocity(psi, tmp, direction, weighted)
         call this%velocity(qpsi, other, direction, weighted)
      case('spin_torque')
         call this%hamiltonian(psi, tmp, .false.)
         call this%hamiltonian(qpsi, other, .false.)
      case('orbital_torque')
         call this%hamiltonian(psi, tmp)
         call this%hamiltonian(qpsi, other)
      case('soc_torque')
         call this%local_action(this%soc, psi, tmp)
         call this%local_action(this%soc, qpsi, other)
      case default
         error stop 'Invalid transport response kind'
      end select
      do i = 1, size(psi, 3)
         if (kind == 'current') then
            out(:, :, i) = 0.5_rp*(matmul(q, tmp(:, :, i)) + other(:, :, i))
         else
            out(:, :, i) = cmplx(0.0_rp, -1.0_rp, rp)*(matmul(q, tmp(:, :, i)) - other(:, :, i))
         end if
      end do
   end subroutine
end module
