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
   use lr_lmto_turek_gf_mod, only: native_structure_constants, native_path_operator_from_structure, native_exchange_trace, &
      native_screening_alpha
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
      ! If positive, choose the ellipse height to enclose this even number
      ! of Matsubara poles.  Zero retains the historical height-fraction
      ! prescription.
      integer :: target_fermi_poles = 0
   end type native_turek_contour_options

   type, public :: native_turek_contour_report
      integer :: contour_points = 0
      integer :: fermi_poles = 0
      integer :: k_points = 0
      integer :: q_points = 0
      real(rp) :: solve_seconds = 0.0_rp
      real(rp) :: contour_seconds = 0.0_rp
      real(rp) :: pole_seconds = 0.0_rp
      real(rp) :: native_max_ellipse_value = 0.0_rp
      logical :: native_bounds_verified = .false.
      integer :: native_spectral_poles = 0
   end type native_turek_contour_report

   public :: native_build_contour
   public :: native_exchange_q_contour
   public :: native_exchange_q_ordered_contour
   public :: native_exchange_jij_contour
   public :: native_exchange_pairs_contour
   public :: native_fourier_jq_to_jij
   public :: native_fourier_jij_to_jq
   public :: native_fourier_complex_jq_to_jij
   public :: native_fourier_complex_jij_to_jq
   public :: native_spin_site_block
   public :: native_complex_fermi
   public :: native_regularized_fermi
   public :: native_spectral_bounds
   public :: native_turek_static_reference
   public :: native_native_poles_from_structure

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

      call native_contour_axes(energy_bounds, fermi, kT, options, center, semimajor, semiminor)

      npoint = options%contour_points
      allocate(nodes(npoint), weights(npoint))
      dtheta = 2.0_rp*pi/real(npoint,rp)
      do i = 1, npoint
         theta = dtheta*real(i-1,rp)
         nodes(i) = cmplx(center+semimajor*cos(theta), semiminor*sin(theta), rp)
         weights(i) = ((-semimajor*sin(theta)+i_unit*semiminor*cos(theta))*dtheta)/(2.0_rp*pi*i_unit)
      end do

      if (kT > 0.0_rp .and. options%account_fermi_poles) then
         max_l = max(2, int(ceiling(semiminor/(pi*kT)))+4)
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

   subroutine native_contour_axes(energy_bounds, fermi, kT, options, center, semimajor, semiminor)
      real(rp), intent(in) :: energy_bounds(2), fermi, kT
      type(native_turek_contour_options), intent(in) :: options
      real(rp), intent(out) :: center, semimajor, semiminor
      real(rp) :: horizontal, scale, y_inside, y_outside
      integer :: npair

      center = 0.5_rp*sum(energy_bounds)
      semimajor = max(0.5_rp*(energy_bounds(2)-energy_bounds(1))+options%contour_margin, options%contour_margin)
      if (options%target_fermi_poles > 0) then
         if (.not. options%account_fermi_poles) error stop 'native contour: target poles requires pole accounting'
         if (kT <= 0.0_rp .or. mod(options%target_fermi_poles,2) /= 0) then
            error stop 'native contour: target_fermi_poles must be a positive even number at finite temperature'
         end if
         npair = options%target_fermi_poles/2
         horizontal = (fermi-center)/semimajor
         scale = sqrt(max(1.0_rp-horizontal*horizontal, 1.0e-14_rp))
         y_inside = pi*kT*real(2*npair-1,rp)
         y_outside = pi*kT*real(2*npair+1,rp)
         ! Place the ellipse midpoint between adjacent Matsubara pairs.
         semiminor = 0.5_rp*(y_inside+y_outside)/scale
      else
         semiminor = max(options%contour_height_fraction*semimajor, 1.0e-8_rp)
         if (kT > 0.0_rp .and. .not. options%account_fermi_poles) semiminor = min(semiminor, 0.45_rp*pi*kT)
      end if
   end subroutine native_contour_axes

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
   !> The pole-subtracted Fermi weight is evaluated on the contour and the
   !> explicit +kT residues are added back.  This is the finite-temperature
   !> identity certified by the scalar oracle; every resolvent here is still a
   !> native P-S solve.
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
         report = native_turek_contour_report()
         report%contour_points = size(nodes)*nk*nq
         report%fermi_poles = size(poles)*nk*nq
         report%k_points = nk; report%q_points = nq
         report%solve_seconds = solve_seconds
         report%contour_seconds = contour_seconds
         report%pole_seconds = pole_seconds
      end if
      deallocate(nodes,weights,poles,integral,smat_src,smat_end,sorb,pmat,gsrc,gend,delta,trace_value)
   end subroutine native_exchange_q_contour

   !> Evaluate both ordered reciprocal channels without collapsing their
   !> complex sublattice phases.  The scalar production wrapper above remains
   !> for compatibility; this is the independent ud/du route used by R1.
   subroutine native_exchange_q_ordered_contour(lat, k_points, k_weights, q_points, fermi, kT, energy_bounds, options, &
                                                 jq_ud, jq_du, report)
      type(lattice), intent(inout) :: lat
      real(rp), intent(in) :: k_points(:, :), k_weights(:), q_points(:, :), fermi, kT, energy_bounds(2)
      type(native_turek_contour_options), intent(in) :: options
      complex(rp), intent(out) :: jq_ud(:, :, :), jq_du(:, :, :)
      type(native_turek_contour_report), intent(out), optional :: report
      complex(rp), allocatable :: nodes(:), weights(:), poles(:)
      complex(rp), allocatable :: integral_ud(:, :, :), integral_du(:, :, :)
      complex(rp), allocatable :: smat_src(:, :, :), smat_end(:, :, :), sorb(:, :)
      complex(rp), allocatable :: pmat(:, :), gsrc(:, :), gend(:, :), delta(:, :, :)
      complex(rp), allocatable :: trace_ud(:, :), trace_du(:, :)
      complex(rp) :: coefficient
      real(rp) :: weight_sum
      integer :: nsite, norb, nk, nq, inode, ipole, ik, iq, nspin
      integer :: clock_start, clock_stop, clock_rate, contour_start, contour_stop, pole_start, pole_stop
      real(rp) :: solve_seconds, contour_seconds, pole_seconds

      nsite = lat%nrec; nk = size(k_points,2); nq = size(q_points,2)
      if (nsite < 1 .or. size(k_points,1) /= 3 .or. size(q_points,1) /= 3) then
         error stop 'native_exchange_q_ordered_contour: point shape mismatch'
      end if
      if (size(k_weights) /= nk .or. nk < 1 .or. nq < 1) error stop 'native_exchange_q_ordered_contour: mesh shape mismatch'
      weight_sum = sum(k_weights)
      if (weight_sum <= tiny(1.0_rp)) error stop 'native_exchange_q_ordered_contour: zero k-weight sum'
      norb = (lat%symbolic_atoms(lat%iz(lat%atlist(1)))%potential%lmax+1)**2
      nspin = 2*norb*nsite
      if (any(shape(jq_ud) /= [nsite,nsite,nq]) .or. any(shape(jq_du) /= [nsite,nsite,nq])) then
         error stop 'native_exchange_q_ordered_contour: output shape mismatch'
      end if

      call native_build_contour(energy_bounds, fermi, kT, options, nodes, weights, poles)
      allocate(integral_ud(nsite,nsite,nq), integral_du(nsite,nsite,nq), smat_src(nspin,nspin,nk), &
         smat_end(nspin,nspin,nk), sorb(norb*nsite,norb*nsite), pmat(nspin,nspin), gsrc(nspin,nspin), &
         gend(nspin,nspin), delta(norb,norb,nsite), trace_ud(nsite,nsite), trace_du(nsite,nsite))
      integral_ud = cmplx(0.0_rp,0.0_rp,rp); integral_du = integral_ud
      solve_seconds = 0.0_rp; contour_seconds = 0.0_rp; pole_seconds = 0.0_rp
      call system_clock(count_rate=clock_rate)

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
         do inode = 1, size(nodes)
            do ik = 1, nk
               call system_clock(clock_start)
               call native_path_operator_from_structure(lat,nodes(inode),smat_src(:,:,ik),pmat,gsrc)
               call native_path_operator_from_structure(lat,nodes(inode),smat_end(:,:,ik),pmat,gend)
               call system_clock(clock_stop)
               solve_seconds = solve_seconds + elapsed_seconds(clock_start,clock_stop,clock_rate)
               call native_deltas(pmat,nsite,norb,delta)
               call native_trace_matrix_ordered(gsrc,gend,delta,norb,nsite,trace_ud,trace_du)
               coefficient = weights(inode)
               if (kT > 0.0_rp) coefficient = coefficient*native_regularized_fermi(nodes(inode),fermi,kT,poles)
               integral_ud(:,:,iq) = integral_ud(:,:,iq) + (k_weights(ik)/weight_sum)*coefficient*trace_ud
               integral_du(:,:,iq) = integral_du(:,:,iq) + (k_weights(ik)/weight_sum)*coefficient*trace_du
            end do
         end do
      end do
      call system_clock(contour_stop)
      contour_seconds = elapsed_seconds(contour_start,contour_stop,clock_rate)

      call system_clock(pole_start)
      do iq = 1, nq
         do ik = 1, nk
            call native_structure_constants(lat,k_points(:,ik)+q_points(:,iq),sorb)
            call native_expand_spin_structure(sorb,nsite,norb,smat_end(:,:,ik))
         end do
         do ipole = 1, size(poles)
            do ik = 1, nk
               call system_clock(clock_start)
               call native_path_operator_from_structure(lat,poles(ipole),smat_src(:,:,ik),pmat,gsrc)
               call native_path_operator_from_structure(lat,poles(ipole),smat_end(:,:,ik),pmat,gend)
               call system_clock(clock_stop)
               solve_seconds = solve_seconds + elapsed_seconds(clock_start,clock_stop,clock_rate)
               call native_deltas(pmat,nsite,norb,delta)
               call native_trace_matrix_ordered(gsrc,gend,delta,norb,nsite,trace_ud,trace_du)
               integral_ud(:,:,iq) = integral_ud(:,:,iq) + (k_weights(ik)/weight_sum)*cmplx(kT,0.0_rp,rp)*trace_ud
               integral_du(:,:,iq) = integral_du(:,:,iq) + (k_weights(ik)/weight_sum)*cmplx(kT,0.0_rp,rp)*trace_du
            end do
         end do
      end do
      call system_clock(pole_stop)
      pole_seconds = elapsed_seconds(pole_start,pole_stop,clock_rate)

      jq_ud = -integral_ud/4.0_rp
      jq_du = -integral_du/4.0_rp
      if (present(report)) then
         report = native_turek_contour_report()
         report%contour_points = size(nodes)*nk*nq
         report%fermi_poles = size(poles)*nk*nq
         report%k_points = nk; report%q_points = nq
         report%solve_seconds = solve_seconds
         report%contour_seconds = contour_seconds
         report%pole_seconds = pole_seconds
      end if
      deallocate(nodes,weights,poles,integral_ud,integral_du,smat_src,smat_end,sorb,pmat,gsrc,gend,delta,trace_ud,trace_du)
   end subroutine native_exchange_q_ordered_contour

   !> Evaluate real-space ordered pairs directly from Fourier transformed
   !> native Green-function blocks.  This routine never forms J(q), so its
   !> result is independent of the reciprocal convolution route.
   subroutine native_exchange_pairs_contour(lat, k_points, k_weights, real_space_vectors, fermi, kT, energy_bounds, options, &
                                             jij_ud, jij_du, report)
      type(lattice), intent(inout) :: lat
      real(rp), intent(in) :: k_points(:, :), k_weights(:), real_space_vectors(:, :), fermi, kT, energy_bounds(2)
      type(native_turek_contour_options), intent(in) :: options
      complex(rp), intent(out) :: jij_ud(:, :, :), jij_du(:, :, :)
      type(native_turek_contour_report), intent(out), optional :: report
      complex(rp), allocatable :: nodes(:), weights(:), poles(:)
      complex(rp), allocatable :: integral_ud(:, :, :), integral_du(:, :, :)
      complex(rp), allocatable :: smat_src(:, :, :), sorb(:, :), pmat(:, :), gmat(:, :), delta(:, :, :)
      complex(rp), allocatable :: gup_k(:, :, :, :, :), gdown_k(:, :, :, :, :)
      complex(rp), allocatable :: gup_r(:, :, :, :, :), gup_minus_r(:, :, :, :, :)
      complex(rp), allocatable :: gdown_r(:, :, :, :, :), gdown_minus_r(:, :, :, :, :)
      complex(rp) :: coefficient, phase_minus, phase_plus
      complex(rp), allocatable :: block(:, :)
      real(rp) :: weight_sum
      integer :: nsite, norb, nk, nvec, nspin, inode, ipole, ik, ir, ia, ja
      integer :: clock_start, clock_stop, clock_rate, contour_start, contour_stop, pole_start, pole_stop
      real(rp) :: solve_seconds, contour_seconds, pole_seconds

      nsite = lat%nrec; nk = size(k_points,2); nvec = size(real_space_vectors,2)
      if (size(k_points,1) /= 3 .or. size(real_space_vectors,1) /= 3 .or. nk < 1 .or. nvec < 1) then
         error stop 'native_exchange_pairs_contour: point shape mismatch'
      end if
      if (size(k_weights) /= nk) error stop 'native_exchange_pairs_contour: weight shape mismatch'
      weight_sum = sum(k_weights)
      if (weight_sum <= tiny(1.0_rp)) error stop 'native_exchange_pairs_contour: zero k-weight sum'
      norb = (lat%symbolic_atoms(lat%iz(lat%atlist(1)))%potential%lmax+1)**2
      nspin = 2*norb*nsite
      if (any(shape(jij_ud) /= [nsite,nsite,nvec]) .or. any(shape(jij_du) /= [nsite,nsite,nvec])) then
         error stop 'native_exchange_pairs_contour: output shape mismatch'
      end if
      call native_build_contour(energy_bounds, fermi, kT, options, nodes, weights, poles)
      allocate(integral_ud(nsite,nsite,nvec), integral_du(nsite,nsite,nvec), smat_src(nspin,nspin,nk), &
         sorb(norb*nsite,norb*nsite), pmat(nspin,nspin), gmat(nspin,nspin), delta(norb,norb,nsite), &
         gup_k(norb,norb,nsite,nsite,nk), gdown_k(norb,norb,nsite,nsite,nk), &
         gup_r(norb,norb,nsite,nsite,nvec), gup_minus_r(norb,norb,nsite,nsite,nvec), &
         gdown_r(norb,norb,nsite,nsite,nvec), gdown_minus_r(norb,norb,nsite,nsite,nvec))
      allocate(block(norb,norb))
      integral_ud = cmplx(0.0_rp,0.0_rp,rp); integral_du = integral_ud
      solve_seconds = 0.0_rp; contour_seconds = 0.0_rp; pole_seconds = 0.0_rp
      call system_clock(count_rate=clock_rate)
      do ik = 1, nk
         call native_structure_constants(lat,k_points(:,ik),sorb)
         call native_expand_spin_structure(sorb,nsite,norb,smat_src(:,:,ik))
      end do

      call system_clock(contour_start)
      do inode = 1, size(nodes)
         do ik = 1, nk
            call system_clock(clock_start)
            call native_path_operator_from_structure(lat,nodes(inode),smat_src(:,:,ik),pmat,gmat)
            call system_clock(clock_stop)
            solve_seconds = solve_seconds + elapsed_seconds(clock_start,clock_stop,clock_rate)
            call native_deltas(pmat,nsite,norb,delta)
            do ia = 1, nsite
               do ja = 1, nsite
                  call native_spin_site_block(gmat,norb,ia,1,ja,1,block)
                  gup_k(:,:,ia,ja,ik) = block
                  call native_spin_site_block(gmat,norb,ia,2,ja,2,block)
                  gdown_k(:,:,ia,ja,ik) = block
               end do
            end do
         end do
         do ir = 1, nvec
            gup_r(:,:,:,:,ir) = cmplx(0.0_rp,0.0_rp,rp)
            gup_minus_r(:,:,:,:,ir) = cmplx(0.0_rp,0.0_rp,rp)
            gdown_r(:,:,:,:,ir) = cmplx(0.0_rp,0.0_rp,rp)
            gdown_minus_r(:,:,:,:,ir) = cmplx(0.0_rp,0.0_rp,rp)
            do ik = 1, nk
               phase_minus = exp(-i_unit*2.0_rp*pi*dot_product(k_points(:,ik),real_space_vectors(:,ir)))
               phase_plus = conjg(phase_minus)
               gup_r(:,:,:,:,ir) = gup_r(:,:,:,:,ir) + (k_weights(ik)/weight_sum)*phase_minus*gup_k(:,:,:,:,ik)
               gup_minus_r(:,:,:,:,ir) = gup_minus_r(:,:,:,:,ir) + (k_weights(ik)/weight_sum)*phase_plus*gup_k(:,:,:,:,ik)
               gdown_r(:,:,:,:,ir) = gdown_r(:,:,:,:,ir) + (k_weights(ik)/weight_sum)*phase_minus*gdown_k(:,:,:,:,ik)
               gdown_minus_r(:,:,:,:,ir) = gdown_minus_r(:,:,:,:,ir) + (k_weights(ik)/weight_sum)*phase_plus*gdown_k(:,:,:,:,ik)
            end do
            coefficient = weights(inode)
            if (kT > 0.0_rp) coefficient = coefficient*native_regularized_fermi(nodes(inode),fermi,kT,poles)
            do ia = 1, nsite
               do ja = 1, nsite
                  integral_ud(ia,ja,ir) = integral_ud(ia,ja,ir) + coefficient*native_exchange_trace(delta(:,:,ia),delta(:,:,ja), &
                     gup_r(:,:,ia,ja,ir),gdown_minus_r(:,:,ja,ia,ir))
                  integral_du(ia,ja,ir) = integral_du(ia,ja,ir) + coefficient*native_exchange_trace(delta(:,:,ia),delta(:,:,ja), &
                     gdown_r(:,:,ia,ja,ir),gup_minus_r(:,:,ja,ia,ir))
               end do
            end do
         end do
      end do
      call system_clock(contour_stop)
      contour_seconds = elapsed_seconds(contour_start,contour_stop,clock_rate)

      call system_clock(pole_start)
      do ipole = 1, size(poles)
         do ik = 1, nk
            call system_clock(clock_start)
            call native_path_operator_from_structure(lat,poles(ipole),smat_src(:,:,ik),pmat,gmat)
            call system_clock(clock_stop)
            solve_seconds = solve_seconds + elapsed_seconds(clock_start,clock_stop,clock_rate)
            call native_deltas(pmat,nsite,norb,delta)
            do ia = 1, nsite
               do ja = 1, nsite
                  call native_spin_site_block(gmat,norb,ia,1,ja,1,block)
                  gup_k(:,:,ia,ja,ik) = block
                  call native_spin_site_block(gmat,norb,ia,2,ja,2,block)
                  gdown_k(:,:,ia,ja,ik) = block
               end do
            end do
         end do
         do ir = 1, nvec
            gup_r(:,:,:,:,ir) = cmplx(0.0_rp,0.0_rp,rp)
            gup_minus_r(:,:,:,:,ir) = cmplx(0.0_rp,0.0_rp,rp)
            gdown_r(:,:,:,:,ir) = cmplx(0.0_rp,0.0_rp,rp)
            gdown_minus_r(:,:,:,:,ir) = cmplx(0.0_rp,0.0_rp,rp)
            do ik = 1, nk
               phase_minus = exp(-i_unit*2.0_rp*pi*dot_product(k_points(:,ik),real_space_vectors(:,ir)))
               phase_plus = conjg(phase_minus)
               gup_r(:,:,:,:,ir) = gup_r(:,:,:,:,ir) + (k_weights(ik)/weight_sum)*phase_minus*gup_k(:,:,:,:,ik)
               gup_minus_r(:,:,:,:,ir) = gup_minus_r(:,:,:,:,ir) + (k_weights(ik)/weight_sum)*phase_plus*gup_k(:,:,:,:,ik)
               gdown_r(:,:,:,:,ir) = gdown_r(:,:,:,:,ir) + (k_weights(ik)/weight_sum)*phase_minus*gdown_k(:,:,:,:,ik)
               gdown_minus_r(:,:,:,:,ir) = gdown_minus_r(:,:,:,:,ir) + (k_weights(ik)/weight_sum)*phase_plus*gdown_k(:,:,:,:,ik)
            end do
            do ia = 1, nsite
               do ja = 1, nsite
                  integral_ud(ia,ja,ir) = integral_ud(ia,ja,ir) + kT*native_exchange_trace(delta(:,:,ia),delta(:,:,ja), &
                     gup_r(:,:,ia,ja,ir),gdown_minus_r(:,:,ja,ia,ir))
                  integral_du(ia,ja,ir) = integral_du(ia,ja,ir) + kT*native_exchange_trace(delta(:,:,ia),delta(:,:,ja), &
                     gdown_r(:,:,ia,ja,ir),gup_minus_r(:,:,ja,ia,ir))
               end do
            end do
         end do
      end do
      call system_clock(pole_stop)
      pole_seconds = elapsed_seconds(pole_start,pole_stop,clock_rate)
      jij_ud = -integral_ud/4.0_rp; jij_du = -integral_du/4.0_rp
      if (present(report)) then
         report = native_turek_contour_report()
         report%contour_points = size(nodes)*nk
         report%fermi_poles = size(poles)*nk
         report%k_points = nk; report%q_points = nvec
         report%solve_seconds = solve_seconds
         report%contour_seconds = contour_seconds
         report%pole_seconds = pole_seconds
      end if
      deallocate(nodes,weights,poles,integral_ud,integral_du,smat_src,sorb,pmat,gmat,delta,gup_k,gdown_k,gup_r, &
         gup_minus_r,gdown_r,gdown_minus_r,block)
   end subroutine native_exchange_pairs_contour

   !> Evaluate a common native q mesh and close it to real-space Jij.
   !>
   !> This wrapper intentionally performs no interpolation or fitted scale:
   !> `jq` is the native contour result and `jij` is its explicitly normalized
   !> discrete inverse Fourier transform.  This is retained as an algebraic
   !> DFT-helper regression; it is not the independent native pair route.
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

   subroutine native_fourier_complex_jq_to_jij(q_points, jq, real_space_vectors, jij)
      real(rp), intent(in) :: q_points(:, :), real_space_vectors(:, :)
      complex(rp), intent(in) :: jq(:, :, :)
      complex(rp), intent(out) :: jij(:, :, :)
      integer :: nq, nsite, nvec, iq, ir, ia, ja

      nq = size(q_points,2); nsite = size(jq,1); nvec = size(real_space_vectors,2)
      if (size(q_points,1) /= 3 .or. size(jq,2) /= nsite .or. size(jq,3) /= nq .or. &
          size(real_space_vectors,1) /= 3 .or. any(shape(jij) /= [nsite,nsite,nvec])) then
         error stop 'native_fourier_complex_jq_to_jij: shape mismatch'
      end if
      do ir = 1, nvec
         do ia = 1, nsite
            do ja = 1, nsite
               jij(ia,ja,ir) = cmplx(0.0_rp,0.0_rp,rp)
               do iq = 1, nq
                  jij(ia,ja,ir) = jij(ia,ja,ir) + jq(ia,ja,iq)* &
                     exp(i_unit*2.0_rp*pi*dot_product(q_points(:,iq),real_space_vectors(:,ir)))/real(nq,rp)
               end do
            end do
         end do
      end do
   end subroutine native_fourier_complex_jq_to_jij

   subroutine native_fourier_complex_jij_to_jq(q_points, real_space_vectors, jij, jq)
      real(rp), intent(in) :: q_points(:, :), real_space_vectors(:, :)
      complex(rp), intent(in) :: jij(:, :, :)
      complex(rp), intent(out) :: jq(:, :, :)
      integer :: nq, nsite, nvec, iq, ir, ia, ja

      nq = size(q_points,2); nsite = size(jij,1); nvec = size(real_space_vectors,2)
      if (size(q_points,1) /= 3 .or. size(jij,2) /= nsite .or. size(jij,3) /= nvec .or. &
          size(real_space_vectors,1) /= 3 .or. any(shape(jq) /= [nsite,nsite,nq])) then
         error stop 'native_fourier_complex_jij_to_jq: shape mismatch'
      end if
      do iq = 1, nq
         do ia = 1, nsite
            do ja = 1, nsite
               jq(ia,ja,iq) = cmplx(0.0_rp,0.0_rp,rp)
               do ir = 1, nvec
                  jq(ia,ja,iq) = jq(ia,ja,iq) + jij(ia,ja,ir)* &
                     exp(-i_unit*2.0_rp*pi*dot_product(q_points(:,iq),real_space_vectors(:,ir)))
               end do
            end do
         end do
      end do
   end subroutine native_fourier_complex_jij_to_jq

   !> Determine bounds from the exact native P-S coefficient problem.
   !>
   !> With P0(z)=D(z-C), Q=qi-alpha, and Palpha=P0(I+Q P0)^-1,
   !> multiplying (Palpha-S) by (I+Q P0) gives the linear coefficient
   !> problem [(I-SQ)D] z = (I-SQ)DC+S.  Its eigenvalues are used only for
   !> contour selection/validation; all exchange resolvents remain direct
   !> P-S solves.
   subroutine native_spectral_bounds(lat, k_points, fermi, kT, options, energy_bounds, max_ellipse_value, all_inside, pole_count)
      type(lattice), intent(inout) :: lat
      real(rp), intent(in) :: k_points(:, :), fermi, kT
      type(native_turek_contour_options), intent(in) :: options
      real(rp), intent(out) :: energy_bounds(2), max_ellipse_value
      logical, intent(out) :: all_inside
      integer, intent(out) :: pole_count
      complex(rp), allocatable :: smat_orb(:, :), smat_spin(:, :), roots(:, :)
      real(rp) :: center, semimajor, semiminor, ellipse_value
      integer :: nsite, norb, nspin, nk, ik, ip

      if (size(k_points,1) /= 3 .or. size(k_points,2) < 1) error stop 'native_spectral_bounds: point shape mismatch'
      nsite = lat%nrec; nk = size(k_points,2)
      norb = (lat%symbolic_atoms(lat%iz(lat%atlist(1)))%potential%lmax+1)**2
      nspin = 2*norb*nsite
      allocate(smat_orb(norb*nsite,norb*nsite), smat_spin(nspin,nspin), roots(nspin,nk))
      do ik = 1, nk
         call native_structure_constants(lat,k_points(:,ik),smat_orb)
         call native_expand_spin_structure(smat_orb,nsite,norb,smat_spin)
         call native_native_poles_from_structure(lat,smat_spin,roots(:,ik))
      end do
      energy_bounds = [minval(real(roots,rp)),maxval(real(roots,rp))]
      if (energy_bounds(2)-energy_bounds(1) < 1.0e-10_rp) then
         energy_bounds = energy_bounds + [-0.5_rp,0.5_rp]
      end if
      call native_contour_axes(energy_bounds,fermi,kT,options,center,semimajor,semiminor)
      max_ellipse_value = -huge(1.0_rp)
      do ik = 1, nk
         do ip = 1, nspin
            ellipse_value = ((real(roots(ip,ik),rp)-center)/semimajor)**2 + &
               (aimag(roots(ip,ik))/semiminor)**2
            max_ellipse_value = max(max_ellipse_value,ellipse_value)
         end do
      end do
      pole_count = nspin*nk
      all_inside = max_ellipse_value < 1.0_rp-1.0e-10_rp
      deallocate(smat_orb,smat_spin,roots)
   end subroutine native_spectral_bounds

   !> Shared native static reference from accepted reciprocal-state metadata.
   !> Bounds/poles select and validate the contour; all exchange resolvents are
   !> still native P-S path-operator solves. Output arrays retain sublattice
   !> phases so finite-q consumers can compare representations without a scale.
   subroutine native_turek_static_reference(lat, k_points, k_weights, q_points, fermi, kT, options, &
      jq_ud, jq_du, jq_sym, delta_j, curvature, report)
      type(lattice), intent(inout) :: lat
      real(rp), intent(in) :: k_points(:, :), k_weights(:), q_points(:, :), fermi, kT
      type(native_turek_contour_options), intent(in) :: options
      complex(rp), intent(out) :: jq_ud(:, :, :), jq_du(:, :, :), jq_sym(:, :, :), delta_j(:, :, :), curvature(:, :, :)
      type(native_turek_contour_report), intent(out) :: report
      real(rp), allocatable :: native_points(:, :)
      real(rp) :: energy_bounds(2), max_ellipse_value
      integer :: nsite, nk, nq, iq, ik, gamma_index, pole_count
      logical :: bounds_verified

      nsite=lat%nrec; nk=size(k_points,2); nq=size(q_points,2)
      if (size(k_points,1)/=3 .or. size(q_points,1)/=3 .or. size(k_weights)/=nk .or. nk<1 .or. nq<1) then
         error stop 'native_turek_static_reference: input shape mismatch'
      end if
      if (any(shape(jq_ud)/=[nsite,nsite,nq]) .or. any(shape(jq_du)/=[nsite,nsite,nq]) .or. &
          any(shape(jq_sym)/=[nsite,nsite,nq]) .or. any(shape(delta_j)/=[nsite,nsite,nq]) .or. &
          any(shape(curvature)/=[nsite,nsite,nq])) error stop 'native_turek_static_reference: output shape mismatch'

      allocate(native_points(3,nk*nq))
      do iq=1,nq
         do ik=1,nk
            native_points(:,(iq-1)*nk+ik)=k_points(:,ik)+q_points(:,iq)
         end do
      end do
      call native_spectral_bounds(lat,native_points,fermi,kT,options,energy_bounds,max_ellipse_value,bounds_verified,pole_count)
      if (.not. bounds_verified) error stop 'native_turek_static_reference: native spectral poles are not inside contour'
      call native_exchange_q_ordered_contour(lat,k_points,k_weights,q_points,fermi,kT,energy_bounds,options,jq_ud,jq_du,report)
      report%native_max_ellipse_value=max_ellipse_value
      report%native_bounds_verified=bounds_verified
      report%native_spectral_poles=pole_count
      jq_sym=0.5_rp*(jq_ud+jq_du)
      gamma_index=0
      do iq=1,nq
         if (sqrt(sum(q_points(:,iq)**2))<=1.0e-12_rp) then
            gamma_index=iq
            exit
         end if
      end do
      if (gamma_index==0) error stop 'native_turek_static_reference: q path must contain Gamma for DeltaJ'
      delta_j=spread(jq_sym(:,:,gamma_index),3,nq)-jq_sym
      curvature=2.0_rp*delta_j
      deallocate(native_points)
   end subroutine native_turek_static_reference

   subroutine native_native_poles_from_structure(lat, smat, roots)
      type(lattice), intent(inout) :: lat
      complex(rp), intent(in) :: smat(:, :)
      complex(rp), intent(out) :: roots(:)
      complex(rp), allocatable :: mcoef(:, :), a(:, :), b(:, :), work(:), vl(:, :), vr(:, :), query(:)
      real(rp), allocatable :: rwork(:), diag_q(:), diag_d(:), diag_c(:), alpha(:)
      integer, allocatable :: ipiv(:)
      integer :: nsite, norb, n, lmax, site, spin, l, m, mls, idx, i, j, ntype, ia, it, info, lwork
      external :: zgesv, zgeev

      n = size(roots); nsite = lat%nrec
      norb = (lat%symbolic_atoms(lat%iz(lat%atlist(1)))%potential%lmax+1)**2
      lmax = lat%symbolic_atoms(lat%iz(lat%atlist(1)))%potential%lmax
      if (size(smat,1) /= n .or. size(smat,2) /= n .or. n /= 2*norb*nsite) then
         error stop 'native_native_poles_from_structure: shape mismatch'
      end if
      allocate(mcoef(n,n),a(n,n),b(n,n),work(n),vl(1,1),vr(1,1),query(1),rwork(2*n),ipiv(n), &
         diag_q(n),diag_d(n),diag_c(n),alpha(0:lmax))
      diag_q = 0.0_rp; diag_d = 0.0_rp; diag_c = 0.0_rp
      do site = 1, nsite
         ntype = lat%ib(site); ia = lat%atlist(ntype); it = lat%iz(ia)
         call native_screening_alpha(lat%symbolic_atoms(it),alpha)
         do spin = 1, 2
            do l = 0, lmax
               if (abs(lat%symbolic_atoms(it)%potential%dele(l,spin)) <= tiny(1.0_rp)) then
                  error stop 'native_native_poles_from_structure: zero potential width'
               end if
               do m = 1, 2*l+1
                  mls = l*l+m; idx = (site-1)*2*norb+(spin-1)*norb+mls
                  diag_q(idx) = lat%symbolic_atoms(it)%potential%qi(l,spin)-alpha(l)
                  diag_d(idx) = 1.0_rp/(lat%symbolic_atoms(it)%potential%dele(l,spin)**2)
                  diag_c(idx) = lat%symbolic_atoms(it)%potential%c(l,spin)+lat%symbolic_atoms(it)%potential%vmad
               end do
            end do
         end do
      end do
      ! Keep the generalized coefficient problem explicit.  Q and D are
      ! diagonal, but their side of the dense structure matrix matters:
      ! M = I - S Q, A = M D, B = M D C + S.
      mcoef = cmplx(0.0_rp,0.0_rp,rp)
      do i = 1, n
         mcoef(i,i) = cmplx(1.0_rp,0.0_rp,rp)
      end do
      do j = 1, n
         mcoef(:,j) = mcoef(:,j)-smat(:,j)*diag_q(j)
      end do
      a = mcoef
      do j = 1, n
         a(:,j) = a(:,j)*diag_d(j)
         b(:,j) = a(:,j)*diag_c(j)+smat(:,j)
      end do
      call zgesv(n,n,a,n,ipiv,b,n,info)
      if (info /= 0) error stop 'native_native_poles_from_structure: coefficient solve failed'
      call zgeev('N','N',n,b,n,roots,vl,1,vr,1,query,-1,rwork,info)
      if (info /= 0) error stop 'native_native_poles_from_structure: eigenvalue workspace query failed'
      lwork = max(1,int(real(query(1),rp)))
      deallocate(work)
      allocate(work(lwork))
      call zgeev('N','N',n,b,n,roots,vl,1,vr,1,work,lwork,rwork,info)
      if (info /= 0) error stop 'native_native_poles_from_structure: eigenvalue solve failed'
      deallocate(mcoef,a,b,work,vl,vr,query,rwork,ipiv,diag_q,diag_d,diag_c,alpha)
   end subroutine native_native_poles_from_structure

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

   subroutine native_trace_matrix_ordered(gsrc, gend, delta, norb, nsite, trace_ud, trace_du)
      complex(rp), intent(in) :: gsrc(:, :), gend(:, :), delta(:, :, :)
      integer, intent(in) :: norb, nsite
      complex(rp), intent(out) :: trace_ud(:, :), trace_du(:, :)
      complex(rp) :: gup_ij(norb,norb), gdown_ji(norb,norb), gdown_ij(norb,norb), gup_ji(norb,norb)
      integer :: ia, ja

      if (any(shape(trace_ud) /= [nsite,nsite]) .or. any(shape(trace_du) /= [nsite,nsite])) then
         error stop 'native contour: ordered trace shape mismatch'
      end if
      do ia = 1, nsite
         do ja = 1, nsite
            call native_spin_site_block(gsrc,norb,ia,1,ja,1,gup_ij)
            call native_spin_site_block(gend,norb,ja,2,ia,2,gdown_ji)
            call native_spin_site_block(gsrc,norb,ia,2,ja,2,gdown_ij)
            call native_spin_site_block(gend,norb,ja,1,ia,1,gup_ji)
            trace_ud(ia,ja) = native_exchange_trace(delta(:,:,ia),delta(:,:,ja),gup_ij,gdown_ji)
            trace_du(ia,ja) = native_exchange_trace(delta(:,:,ia),delta(:,:,ja),gdown_ij,gup_ji)
         end do
      end do
   end subroutine native_trace_matrix_ordered

   pure complex(rp) function native_complex_fermi(z, fermi, kT) result(value)
      complex(rp), intent(in) :: z
      real(rp), intent(in) :: fermi, kT
      complex(rp) :: argument, reduced

      if (kT <= 0.0_rp) then
         value = cmplx(merge(1.0_rp,0.0_rp,real(z,rp) < fermi),0.0_rp,rp)
         return
      end if
      argument = (z-fermi)/kT
      if (real(argument,rp) > 40.0_rp) then
         reduced = exp(-argument); value = reduced/(1.0_rp+reduced)
      else if (real(argument,rp) < -40.0_rp) then
         reduced = exp(argument); value = 1.0_rp/(1.0_rp+reduced)
      else
         value = 1.0_rp/(exp(argument)+1.0_rp)
      end if
   end function native_complex_fermi

   !> Pole-subtracted weight used by the alternative closed-contour identity.
   !> Its contour integral already contains the explicit Fermi-pole residues;
   !> callers must not add another pole sum to this weight.
   pure complex(rp) function native_regularized_fermi(z, fermi, kT, poles) result(value)
      complex(rp), intent(in) :: z, poles(:)
      real(rp), intent(in) :: fermi, kT
      integer :: i

      value = native_complex_fermi(z,fermi,kT)
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
