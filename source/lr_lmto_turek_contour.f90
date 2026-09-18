!------------------------------------------------------------------------------
! DRESP-03TG-C -- native screened LMTO/Turek contour exchange.
!
! This module deliberately works on the path operator itself.  It does not
! obtain g from a finite-H resolvent and it does not use an eigenbasis.  The
! only representation-specific input is the live lattice, whose P and S
! objects are assembled by lr_lmto_turek_gf_mod.
!------------------------------------------------------------------------------
module lr_lmto_turek_contour_mod
   use lattice_mod, only: lattice
   use lr_lmto_turek_gf_mod, only: native_structure_constants, native_path_operator_from_structure, native_exchange_trace
   use math_mod, only: i_unit, pi
   use precision_mod, only: rp
   implicit none
   private

   type, public :: native_turek_contour_options
      integer :: contour_points = 64
      character(len=16) :: contour_shape = 'ellipse'
      real(rp) :: contour_margin = 0.25_rp
      real(rp) :: contour_height_fraction = 0.35_rp
      logical :: account_fermi_poles = .true.
   end type native_turek_contour_options

   type, public :: native_turek_contour_report
      integer :: contour_points = 0
      integer :: fermi_poles = 0
      integer :: k_points = 0
      integer :: q_points = 0
      real(rp) :: solve_seconds = 0.0_rp
      real(rp) :: contour_seconds = 0.0_rp
      real(rp) :: pole_seconds = 0.0_rp
   end type native_turek_contour_report

   public :: native_build_contour
   public :: native_exchange_q_contour
   public :: native_exchange_jij_contour
   public :: native_fourier_jq_to_jij
   public :: native_fourier_jij_to_jq
   public :: native_spin_site_block

contains

   !> Build the counter-clockwise ellipse used by the native contour route.
   !>
   !> The bounds are explicit because the native path operator has no
   !> eigenvalue API.  In production they are taken from the accepted
   !> reciprocal-H state; the native evaluator itself still only solves P-S.
   subroutine native_build_contour(energy_bounds, fermi, kT, options, nodes, weights, fermi_poles)
      real(rp), intent(in) :: energy_bounds(2), fermi, kT
      type(native_turek_contour_options), intent(in) :: options
      complex(rp), allocatable, intent(out) :: nodes(:), weights(:), fermi_poles(:)
      real(rp) :: lower, upper, center, semimajor, semiminor, height, theta, dtheta
      real(rp) :: pole_y, ellipse_value
      integer :: i, l, max_l, npoint, npole

      if (energy_bounds(2) <= energy_bounds(1)) error stop 'native contour: invalid energy bounds'
      if (options%contour_points < 8) error stop 'native contour: contour_points must be at least 8'
      if (trim(options%contour_shape) /= 'ellipse') error stop 'native contour: only ellipse is supported'
      if (options%contour_margin <= 0.0_rp .or. options%contour_height_fraction <= 0.0_rp) then
         error stop 'native contour: contour margin and height fraction must be positive'
      end if

      lower = energy_bounds(1); upper = energy_bounds(2)
      center = 0.5_rp*(lower+upper)
      semimajor = max(0.5_rp*(upper-lower)+options%contour_margin, options%contour_margin)
      height = max(options%contour_height_fraction*semimajor, 1.0e-8_rp)
      if (kT > 0.0_rp .and. .not. options%account_fermi_poles) height = min(height, 0.45_rp*pi*kT)
      semiminor = height

      npoint = options%contour_points
      allocate(nodes(npoint), weights(npoint))
      dtheta = 2.0_rp*pi/real(npoint,rp)
      do i = 1, npoint
         theta = dtheta*real(i-1,rp)
         nodes(i) = cmplx(center+semimajor*cos(theta), semiminor*sin(theta), rp)
         weights(i) = ((-semimajor*sin(theta)+i_unit*semiminor*cos(theta))*dtheta)/(2.0_rp*pi*i_unit)
      end do

      if (kT > 0.0_rp .and. options%account_fermi_poles) then
         max_l = max(2, int(ceiling(semiminor/(pi*kT)))+2)
         npole = 0
         do l = -max_l, max_l
            pole_y = pi*kT*real(2*l+1,rp)
            ellipse_value = ((fermi-center)/semimajor)**2 + (pole_y/semiminor)**2
            if (ellipse_value < 1.0_rp-1.0e-12_rp) npole = npole+1
         end do
         allocate(fermi_poles(npole))
         npole = 0
         do l = -max_l, max_l
            pole_y = pi*kT*real(2*l+1,rp)
            ellipse_value = ((fermi-center)/semimajor)**2 + (pole_y/semiminor)**2
            if (ellipse_value < 1.0_rp-1.0e-12_rp) then
               npole = npole+1
               fermi_poles(npole) = cmplx(fermi,pole_y,rp)
            end if
         end do
      else
         allocate(fermi_poles(0))
      end if
   end subroutine native_build_contour

   !> Evaluate the native Turek exchange matrix at every supplied q point.
   !>
   !> q is in the same direct reciprocal coordinates consumed by the live
   !> structure-factor backend.  The returned quantity is the absolute
   !> J_ij(q), with the convention
   !>
   !>   J_ij(q) = -1/4 Re int_C dz/(2*pi*i) sum_k w_k
   !>              f(z) Tr[DeltaP_i g^up_ij(k,z) DeltaP_j
   !>                 g^down_ji(k+q,z)].
   !>
   !> Fermi poles are removed from the contour integrand and added back as
   !> explicit residues.  This is the same finite-temperature identity used
   !> by lr_kl_contour, but every resolvent here is a native P-S solve.
   subroutine native_exchange_q_contour(lat, k_points, k_weights, q_points, fermi, kT, energy_bounds, options, jq, report)
      type(lattice), intent(inout) :: lat
      real(rp), intent(in) :: k_points(:, :), k_weights(:), q_points(:, :), fermi, kT, energy_bounds(2)
      type(native_turek_contour_options), intent(in) :: options
      real(rp), intent(out) :: jq(:, :, :)
      type(native_turek_contour_report), intent(out), optional :: report
      complex(rp), allocatable :: nodes(:), weights(:), poles(:)
      complex(rp), allocatable :: integral(:, :, :), smat_src(:, :, :), smat_end(:, :, :), sorb(:, :)
      complex(rp), allocatable :: pmat(:, :), gsrc(:, :), gend(:, :)
      complex(rp), allocatable :: delta(:, :, :), trace_value(:, :)
      complex(rp) :: coefficient
      real(rp) :: weight_sum
      integer :: nsite, norb, nlocal, nk, nq, nnode, npole, ik, iq
      integer :: nspin
      integer :: clock_start, clock_stop, clock_rate, contour_start, contour_stop, pole_start, pole_stop
      real(rp) :: solve_seconds, contour_seconds, pole_seconds

      nsite = lat%nrec
      if (nsite < 1) error stop 'native_exchange_q_contour: lattice has no reciprocal sites'
      if (size(k_points,1) /= 3 .or. size(q_points,1) /= 3) error stop 'native_exchange_q_contour: point shape mismatch'
      nk = size(k_points,2); nq = size(q_points,2)
      if (size(k_weights) /= nk .or. nk < 1 .or. nq < 1) error stop 'native_exchange_q_contour: mesh shape mismatch'
      weight_sum = sum(k_weights)
      if (weight_sum <= tiny(1.0_rp)) error stop 'native_exchange_q_contour: zero k-weight sum'
      norb = (lat%symbolic_atoms(lat%iz(lat%atlist(1)))%potential%lmax+1)**2
      nlocal = 2*norb; nspin = nlocal*nsite
      if (any(shape(jq) /= [nsite,nsite,nq])) error stop 'native_exchange_q_contour: output shape mismatch'

      call native_build_contour(energy_bounds, fermi, kT, options, nodes, weights, poles)
      nnode = size(nodes); npole = size(poles)
      allocate(integral(nsite,nsite,nq), smat_src(nspin,nspin,nk), smat_end(nspin,nspin,nk), &
               sorb(norb*nsite,norb*nsite), pmat(nspin,nspin), gsrc(nspin,nspin), gend(nspin,nspin), &
               delta(norb,norb,nsite), trace_value(nsite,nsite))
      integral = cmplx(0.0_rp,0.0_rp,rp)
      solve_seconds = 0.0_rp; contour_seconds = 0.0_rp; pole_seconds = 0.0_rp
      call system_clock(count_rate=clock_rate)

      ! S is energy independent. Cache the source mesh once, and rebuild only
      ! the endpoint mesh when q changes; the inner contour then consists only
      ! of native P(z)-S solves.
      do ik = 1, nk
         call native_structure_constants(lat,k_points(:,ik),sorb)
         call native_expand_spin_structure(sorb,nsite,norb,smat_src(:,:,ik))
      end do

      call system_clock(contour_start)
      do iq = 1, nq
         do ik = 1, nk
            call native_structure_constants(lat,k_points(:,ik)+q_points(:,iq),sorb)
            call native_expand_spin_structure(sorb,nsite,norb,smat_end(:,:,ik))
         end do
         do nnode = 1, size(nodes)
            do ik = 1, nk
               call system_clock(clock_start)
               call native_path_operator_from_structure(lat,nodes(nnode),smat_src(:,:,ik),pmat,gsrc)
               call native_path_operator_from_structure(lat,nodes(nnode),smat_end(:,:,ik),pmat,gend)
               call system_clock(clock_stop)
               solve_seconds = solve_seconds + elapsed_seconds(clock_start,clock_stop,clock_rate)
               call native_deltas(pmat,nsite,norb,delta)
               call native_trace_matrix(gsrc,gend,delta,norb,nsite,trace_value)
               coefficient = weights(nnode)
               if (kT > 0.0_rp) coefficient = coefficient*native_regularized_fermi(nodes(nnode),fermi,kT,poles)
               integral(:,:,iq) = integral(:,:,iq) + (k_weights(ik)/weight_sum)*coefficient*trace_value
            end do
         end do
      end do
      call system_clock(contour_stop)
      contour_seconds = elapsed_seconds(contour_start,contour_stop,clock_rate)

      ! The contour integrand has all enclosed Fermi poles removed.  The
      ! residue of f is -kT, so the physical contour adds +kT*R(z_pole).
      call system_clock(pole_start)
      do iq = 1, nq
         do ik = 1, nk
            call native_structure_constants(lat,k_points(:,ik)+q_points(:,iq),sorb)
            call native_expand_spin_structure(sorb,nsite,norb,smat_end(:,:,ik))
         end do
         do npole = 1, size(poles)
            do ik = 1, nk
               call system_clock(clock_start)
               call native_path_operator_from_structure(lat,poles(npole),smat_src(:,:,ik),pmat,gsrc)
               call native_path_operator_from_structure(lat,poles(npole),smat_end(:,:,ik),pmat,gend)
               call system_clock(clock_stop)
               solve_seconds = solve_seconds + elapsed_seconds(clock_start,clock_stop,clock_rate)
               call native_deltas(pmat,nsite,norb,delta)
               call native_trace_matrix(gsrc,gend,delta,norb,nsite,trace_value)
               integral(:,:,iq) = integral(:,:,iq) + (k_weights(ik)/weight_sum)*cmplx(kT,0.0_rp,rp)*trace_value
            end do
         end do
      end do
      call system_clock(pole_stop)
      pole_seconds = elapsed_seconds(pole_start,pole_stop,clock_rate)

      ! The contour weights already contain 1/(2*pi*i).  For the retarded
      ! real-axis convention Im[1/(E-e+i0)] = -pi*delta(E-e), the LKAG
      ! prefactor 1/(4*pi) therefore becomes -Re[I]/4 after the occupied
      ! contour is closed.  Taking Im(I) here would incorrectly return zero:
      ! the enclosed-pole residues are real.
      jq = -real(integral,rp)/4.0_rp
      if (present(report)) then
         report%contour_points = size(nodes)*nk*nq
         report%fermi_poles = size(poles)*nk*nq
         report%k_points = nk; report%q_points = nq
         report%solve_seconds = solve_seconds
         report%contour_seconds = contour_seconds
         report%pole_seconds = pole_seconds
      end if
      deallocate(nodes,weights,poles,integral,smat_src,smat_end,sorb,pmat,gsrc,gend,delta,trace_value)
   end subroutine native_exchange_q_contour

   !> Evaluate a common native q mesh and close it to real-space Jij.
   !>
   !> This wrapper intentionally performs no interpolation or fitted scale:
   !> `jq` is the native contour result and `jij` is its explicitly normalized
   !> discrete inverse Fourier transform.  For a complete commensurate mesh,
   !> the pair is an exact finite-mesh Fourier closure.
   subroutine native_exchange_jij_contour(lat, k_points, k_weights, q_points, real_space_vectors, fermi, kT, energy_bounds, &
      options, jq, jij, report)
      type(lattice), intent(inout) :: lat
      real(rp), intent(in) :: k_points(:, :), k_weights(:), q_points(:, :), real_space_vectors(:, :), fermi, kT, energy_bounds(2)
      type(native_turek_contour_options), intent(in) :: options
      real(rp), intent(out) :: jq(:, :, :), jij(:, :, :)
      type(native_turek_contour_report), intent(out), optional :: report
      type(native_turek_contour_report) :: local_report

      call native_exchange_q_contour(lat, k_points, k_weights, q_points, fermi, kT, energy_bounds, options, jq, local_report)
      call native_fourier_jq_to_jij(q_points, jq, real_space_vectors, jij)
      if (present(report)) report = local_report
   end subroutine native_exchange_jij_contour

   !> Convert q-space values to real-space coefficients.
   !>
   !> q_points and real_space_vectors are direct reciprocal/lattice
   !> coordinates.  For a complete uniform mesh this is an exact DFT; for a
   !> selected path the routine remains an explicit normalized quadrature and
   !> does not claim a shell-complete interaction.
   subroutine native_fourier_jq_to_jij(q_points, jq, real_space_vectors, jij)
      real(rp), intent(in) :: q_points(:, :), jq(:, :, :), real_space_vectors(:, :)
      real(rp), intent(out) :: jij(:, :, :)
      integer :: nq, nsite, nvec, iq, ir, ia, ja
      complex(rp) :: sum_value

      nq = size(q_points,2); nsite = size(jq,1); nvec = size(real_space_vectors,2)
      if (size(q_points,1) /= 3 .or. size(jq,2) /= nsite .or. size(jq,3) /= nq .or. &
          size(real_space_vectors,1) /= 3 .or. any(shape(jij) /= [nsite,nsite,nvec])) then
         error stop 'native_fourier_jq_to_jij: shape mismatch'
      end if
      do ir = 1, nvec
         do ia = 1, nsite
            do ja = 1, nsite
               sum_value = cmplx(0.0_rp,0.0_rp,rp)
               do iq = 1, nq
                  sum_value = sum_value + cmplx(jq(ia,ja,iq),0.0_rp,rp)* &
                     exp(i_unit*2.0_rp*pi*dot_product(q_points(:,iq),real_space_vectors(:,ir)))
               end do
               jij(ia,ja,ir) = real(sum_value,rp)/real(nq,rp)
            end do
         end do
      end do
   end subroutine native_fourier_jq_to_jij

   !> Reconstruct J_ij(q) from real-space coefficients using the inverse
   !> Fourier convention paired with native_fourier_jq_to_jij.
   subroutine native_fourier_jij_to_jq(q_points, real_space_vectors, jij, jq)
      real(rp), intent(in) :: q_points(:, :), real_space_vectors(:, :), jij(:, :, :)
      real(rp), intent(out) :: jq(:, :, :)
      integer :: nq, nsite, nvec, iq, ir, ia, ja
      complex(rp) :: sum_value

      nq = size(q_points,2); nsite = size(jij,1); nvec = size(real_space_vectors,2)
      if (size(q_points,1) /= 3 .or. size(jij,2) /= nsite .or. size(jij,3) /= nvec .or. &
          size(real_space_vectors,1) /= 3 .or. any(shape(jq) /= [nsite,nsite,nq])) then
         error stop 'native_fourier_jij_to_jq: shape mismatch'
      end if
      do iq = 1, nq
         do ia = 1, nsite
            do ja = 1, nsite
               sum_value = cmplx(0.0_rp,0.0_rp,rp)
               do ir = 1, nvec
                  sum_value = sum_value + cmplx(jij(ia,ja,ir),0.0_rp,rp)* &
                     exp(-i_unit*2.0_rp*pi*dot_product(q_points(:,iq),real_space_vectors(:,ir)))
               end do
               jq(ia,ja,iq) = real(sum_value,rp)
            end do
         end do
      end do
   end subroutine native_fourier_jij_to_jq

   !> Extract one site/spin block from the native site-major spin layout.
   subroutine native_spin_site_block(gmat, norb, site_i, spin_i, site_j, spin_j, block)
      complex(rp), intent(in) :: gmat(:, :)
      integer, intent(in) :: norb, site_i, spin_i, site_j, spin_j
      complex(rp), intent(out) :: block(:, :)
      integer :: nsite, i0, j0

      nsite = size(gmat,1)/(2*norb)
      if (norb < 1 .or. size(gmat,2) /= size(gmat,1) .or. size(block,1) /= norb .or. &
          size(block,2) /= norb .or. site_i < 1 .or. site_i > nsite .or. site_j < 1 .or. site_j > nsite .or. &
          (spin_i /= 1 .and. spin_i /= 2) .or. (spin_j /= 1 .and. spin_j /= 2)) then
         error stop 'native_spin_site_block: invalid shape or index'
      end if
      i0 = (site_i-1)*2*norb+(spin_i-1)*norb
      j0 = (site_j-1)*2*norb+(spin_j-1)*norb
      block = gmat(i0+1:i0+norb,j0+1:j0+norb)
   end subroutine native_spin_site_block

   subroutine native_expand_spin_structure(sorb, nsite, norb, smat)
      complex(rp), intent(in) :: sorb(:, :)
      integer, intent(in) :: nsite, norb
      complex(rp), intent(out) :: smat(:, :)
      integer :: isite, jsite, i0, j0

      if (any(shape(sorb) /= [norb*nsite,norb*nsite]) .or. any(shape(smat) /= [2*norb*nsite,2*norb*nsite])) then
         error stop 'native contour: structure expansion shape mismatch'
      end if
      smat = cmplx(0.0_rp,0.0_rp,rp)
      do isite = 1, nsite
         do jsite = 1, nsite
            i0 = (isite-1)*norb; j0 = (jsite-1)*norb
            smat((isite-1)*2*norb+1:(isite-1)*2*norb+norb, (jsite-1)*2*norb+1:(jsite-1)*2*norb+norb) = &
               sorb(i0+1:i0+norb,j0+1:j0+norb)
            smat((isite-1)*2*norb+norb+1:(isite-1)*2*norb+2*norb, (jsite-1)*2*norb+norb+1:(jsite-1)*2*norb+2*norb) = &
               sorb(i0+1:i0+norb,j0+1:j0+norb)
         end do
      end do
   end subroutine native_expand_spin_structure

   subroutine native_deltas(pmat, nsite, norb, delta)
      complex(rp), intent(in) :: pmat(:, :)
      integer, intent(in) :: nsite, norb
      complex(rp), intent(out) :: delta(:, :, :)
      integer :: site, p0
      if (any(shape(delta) /= [norb,norb,nsite]) .or. any(shape(pmat) /= [2*norb*nsite,2*norb*nsite])) then
         error stop 'native contour: delta shape mismatch'
      end if
      do site = 1, nsite
         p0 = (site-1)*2*norb
         delta(:,:,site) = pmat(p0+1:p0+norb,p0+1:p0+norb)- &
            pmat(p0+norb+1:p0+2*norb,p0+norb+1:p0+2*norb)
      end do
   end subroutine native_deltas

   subroutine native_trace_matrix(gsrc, gend, delta, norb, nsite, trace_value)
      complex(rp), intent(in) :: gsrc(:, :), gend(:, :), delta(:, :, :)
      integer, intent(in) :: norb, nsite
      complex(rp), intent(out) :: trace_value(:, :)
      complex(rp) :: gup(norb,norb), gdown(norb,norb)
      integer :: ia, ja

      if (any(shape(trace_value) /= [nsite,nsite])) error stop 'native contour: trace shape mismatch'
      do ia = 1, nsite
         do ja = 1, nsite
            call native_spin_site_block(gsrc,norb,ia,1,ja,1,gup)
            call native_spin_site_block(gend,norb,ja,2,ia,2,gdown)
            trace_value(ia,ja) = native_exchange_trace(delta(:,:,ia),delta(:,:,ja),gup,gdown)
         end do
      end do
   end subroutine native_trace_matrix

   pure complex(rp) function native_regularized_fermi(z, fermi, kT, poles) result(value)
      complex(rp), intent(in) :: z, poles(:)
      real(rp), intent(in) :: fermi, kT
      complex(rp) :: argument, reduced
      integer :: i

      argument = (z-fermi)/kT
      if (real(argument,rp) > 40.0_rp) then
         reduced = exp(-argument); value = reduced/(1.0_rp+reduced)
      else if (real(argument,rp) < -40.0_rp) then
         reduced = exp(argument); value = 1.0_rp/(1.0_rp+reduced)
      else
         value = 1.0_rp/(exp(argument)+1.0_rp)
      end if
      do i = 1, size(poles)
         value = value+kT/(z-poles(i))
      end do
   end function native_regularized_fermi

   pure function elapsed_seconds(start_count, end_count, rate) result(seconds)
      integer, intent(in) :: start_count, end_count, rate
      real(rp) :: seconds
      if (rate > 0) then
         seconds = real(end_count-start_count,rp)/real(rate,rp)
      else
         seconds = 0.0_rp
      end if
   end function elapsed_seconds

end module lr_lmto_turek_contour_mod
