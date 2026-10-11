!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!> @brief spin_response radial_kernel_matrix on the developer oracle R1 (analytic kernel integrals).
!> @details Usage: test_spin_response_r1 <oracle directory>. Mesh, radial functions, m(r), dV(r) and the expected
!>          K_ll come from <dir>/r1_radial_kernel.nml, the tolerance r1_rel from <dir>/tolerances.nml.
!>          The Simpson weights are formed here from (a, b): dr/di = a (r + b), weights 1/3, 4/3, 2/3, ..., 1/3.
!>          Case 2 has m = 0 on a mesh point.
!------------------------------------------------------------------------------
program test_spin_response_r1
   use precision_mod, only: rp
   use spin_response_mod, only: radial_kernel_matrix
   use spin_response_tolerances_mod
   implicit none

   integer, parameter :: nmax = 4096
   character(len=512) :: dir, array_dims, definition
   integer :: nr, lmax, ncase, icase, iu, nch, i
   real(rp) :: a, b, r0, k_at_r0_limit, simpson_max_rel_error, err
   real(rp) :: r(nmax), gfac_up(8), gfac_dn(8), g_up(nmax), g_dn(nmax), s_up(nmax), s_dn(nmax), m_of_r(nmax), dv_of_r(nmax)
   real(rp) :: kernel_exact(64), kernel_simpson_reference(64)
   real(rp), allocatable :: wq(:), kmat(:, :), kref(:, :), m_case(:, :), dv_case(:, :)
   logical :: r0_on_mesh_point, failed
   namelist /r1_mesh/ nr, a, b, lmax, ncase, array_dims, definition, r, gfac_up, gfac_dn, g_up, g_dn, s_up, s_dn
   namelist /r1_case/ icase, r0, r0_on_mesh_point, k_at_r0_limit, m_of_r, dv_of_r
   namelist /r1_expected/ kernel_exact, kernel_simpson_reference, simpson_max_rel_error

   call get_command_argument(1, dir)
   if (len_trim(dir) == 0) error stop 'usage: test_spin_response_r1 <oracle directory>'
   call read_tolerance_file(trim(dir))
   call require('r1_rel', r1_rel)
   failed = .false.
   write (*, '(a)') 'test | measured | tolerance'

   open (newunit=iu, file=trim(dir)//'/r1_radial_kernel.nml', status='old', action='read')
   read (iu, nml=r1_mesh)
   nch = lmax + 1
   allocate (m_case(nr, ncase), dv_case(nr, ncase), kmat(nch, nch))
   do i = 1, ncase
      read (iu, nml=r1_case)
      m_case(:, i) = m_of_r(1:nr)
      dv_case(:, i) = dv_of_r(1:nr)
   end do
   read (iu, nml=r1_expected)
   close (iu)
   kref = reshape(kernel_exact(1:nch*nch), [nch, nch])
   wq = [(merge(1.0_rp, merge(4.0_rp, 2.0_rp, mod(i, 2) == 0), i == 1 .or. i == nr)/3.0_rp*a*(r(i) + b), i=1, nr)]

   do i = 1, ncase
      call radial_kernel_matrix(r(1:nr), wq, m_case(:, i), dv_case(:, i), reshape(g_up(1:nr*nch), [nr, nch]), &
                                reshape(g_dn(1:nr*nch), [nr, nch]), reshape(s_up(1:nr*nch), [nr, nch]), &
                                reshape(s_dn(1:nr*nch), [nr, nch]), spread(gfac_up(1:nch), 1, nr), &
                                spread(gfac_dn(1:nch), 1, nr), kmat)
      err = maxval(abs(kmat - kref)/abs(kref))
      ! catches: missing 1/(4 pi), large-component product instead of the SR metric, shifted channel index, no m = 0 guard
      write (*, '(a,i1,t58,es10.3,2x,es10.3,2x,a)') 'R1 K_ll vs analytic integral, case ', i, err, r1_rel, merge('ok  ', 'FAIL', err <= r1_rel)
      if (.not. err <= r1_rel) failed = .true.
   end do
   if (failed) error stop 1
   write (*, '(a)') 'RESULT: PASS'
end program test_spin_response_r1
