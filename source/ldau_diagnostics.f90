! Fixed-U trace setup and CSV serialization procedures.
submodule(self_mod) ldau_diagnostics
   implicit none
contains

   module procedure initialize_fixed_u_trace

      character(len=1024) :: trace_path
      character(len=128) :: trace_value
      character(len=128) :: root_run_id
      integer :: ios, unit, root_trace_enabled, root_selected_iteration, root_capture_only
      integer :: root_moments_only, root_ql_only
      real(rp) :: root_trace_lambda

      this%fixed_u_trace_enabled = .false.
      this%fixed_u_trace_filename = ''
      this%fixed_u_trace_run_id = 'run'
      this%fixed_u_trace_lambda = 1.0_rp
      this%fixed_u_trace_selected_iteration = 0
      this%fixed_u_trace_moments_only = .false.
      this%fixed_u_trace_ql_only = .false.

      ! Rank zero owns environment parsing and file creation.  Every rank
      ! reaches the broadcasts below, including the no-trace and file-open
      ! failure cases, so optional M1 paths cannot disagree across ranks.
      root_trace_enabled = 0
      root_selected_iteration = 0
      root_capture_only = 0
      root_moments_only = 0
      root_ql_only = 0
      root_trace_lambda = 1.0_rp
      root_run_id = 'run'
      if (rank == 0) then
         trace_path = ''
         trace_value = ''
         call get_environment_variable('RSLMTO_FIXED_U_TRACE', trace_path, status=ios)
         if (ios == 0 .and. len_trim(trace_path) > 0 .and. &
             this%control%nsp == 1 .and. this%control%lmax == 2 .and. &
             (.not. this%hamiltonian%hoh .or. this%hamiltonian%fixed_u_test_allow_hoh_w1_u0)) then
            ! Explicit trace opt-in also enables the M0 no-U reference.
            ! Ordinary no-U runs still take the untouched path when the
            ! environment variable is absent. Other Hubbard modes remain out.
            if (this%hamiltonian%hubbard_u_general_check) then
               if (trim(this%hamiltonian%hubbard_u_potential_form) == 'liechtenstein') root_trace_enabled = 1
            else if (.not. (this%hamiltonian%hubbard_u_sc_check .or. this%hamiltonian%hubbard_v_check) .and. &
                     (.not. this%hamiltonian%hubbard_u_impurity_check .or. &
                      (trim(this%control%calctype) == 'I' .and. this%control%nsp == 1 .and. &
                       this%control%lmax == 2 .and. &
                       trim(this%hamiltonian%hubbard_radial_moments) == 'spectral'))) then
               root_trace_enabled = 1
            end if
         end if

         if (root_trace_enabled == 1) then
            call get_environment_variable('RSLMTO_FIXED_U_RUN_ID', trace_value, status=ios)
            if (ios == 0 .and. len_trim(trace_value) > 0) root_run_id = trim(trace_value)
            call get_environment_variable('RSLMTO_FIXED_U_LAMBDA', trace_value, status=ios)
            if (ios == 0 .and. len_trim(trace_value) > 0) then
               read(trace_value, *, iostat=ios) root_trace_lambda
               if (ios /= 0) root_trace_enabled = 0
            end if
            ! Zero selects rolling-last mode: replace the trace file at each
            ! iteration start, leaving the last completed iteration on disk.
            call get_environment_variable('RSLMTO_FIXED_U_TRACE_ITERATION', trace_value, status=ios)
            if (ios == 0 .and. len_trim(trace_value) > 0) then
               read(trace_value, *, iostat=ios) root_selected_iteration
               if (ios /= 0 .or. root_selected_iteration < 0) root_trace_enabled = 0
            end if
            trace_value = ''
            call get_environment_variable('RSLMTO_FIXED_U_TRACE_CAPTURE_ONLY', trace_value, status=ios)
            if (ios == 0 .and. (trim(trace_value) == '1' .or. trim(trace_value) == 'true')) root_capture_only = 1
            trace_value = ''
            call get_environment_variable('RSLMTO_FIXED_U_TRACE_MOMENTS_ONLY', trace_value, status=ios)
            if (ios == 0 .and. (trim(trace_value) == '1' .or. trim(trace_value) == 'true')) then
               root_moments_only = 1
            end if
            trace_value = ''
            call get_environment_variable('RSLMTO_FIXED_U_TRACE_QL_ONLY', trace_value, status=ios)
            if (ios == 0 .and. (trim(trace_value) == '1' .or. trim(trace_value) == 'true')) root_ql_only = 1
         end if

         if (root_trace_enabled == 1) then
            open(newunit=unit, file=trim(trace_path), status='replace', action='write', iostat=ios)
            if (ios == 0) then
               write(unit, '(a)', iostat=ios) &
                  'kind,run_id,lambda,iteration,stage,site,type,l,spin,k,band,row,col,order,weight,energy,scalar,re,im'
               close(unit, iostat=ios)
            end if
            if (ios == 0) then
               this%fixed_u_trace_filename = trim(trace_path)
            else
               root_trace_enabled = 0
               call g_logger%warning('Could not create RSLMTO_FIXED_U_TRACE output; trace disabled.', __FILE__, __LINE__)
            end if
         end if
      end if

#ifdef USE_MPI
      call MPI_BCAST(root_trace_enabled, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(root_selected_iteration, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(root_capture_only, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(root_moments_only, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(root_ql_only, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(root_trace_lambda, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(root_run_id, len(root_run_id), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
#endif
      this%fixed_u_trace_enabled = root_trace_enabled == 1
      this%fixed_u_trace_capture_only = root_capture_only == 1
      this%fixed_u_trace_moments_only = root_moments_only == 1
      this%fixed_u_trace_ql_only = root_ql_only == 1
      this%fixed_u_trace_selected_iteration = root_selected_iteration
      this%fixed_u_trace_lambda = root_trace_lambda
      this%fixed_u_trace_run_id = root_run_id
   end procedure initialize_fixed_u_trace


   module procedure reset_fixed_u_trace

      integer :: unit, ios
      if (.not. this%fixed_u_trace_enabled .or. rank /= 0) return
      if (this%fixed_u_trace_ql_only) return
      open(newunit=unit, file=trim(this%fixed_u_trace_filename), status='replace', action='write', iostat=ios)
      if (ios /= 0) return
      write(unit, '(a)') 'kind,run_id,lambda,iteration,stage,site,type,l,spin,k,band,row,col,order,weight,energy,scalar,re,im'
      close(unit)
   end procedure reset_fixed_u_trace


   module procedure write_fixed_u_csv_row
      ! QL-only is the safe history mode used when diagnostic replay helpers
      ! are excluded.  Keep its artifact strictly limited to the accepted
      ! P_OUT occupations, including writes made by independent trace hooks.
      if (this%fixed_u_trace_ql_only .and. (trim(kind) /= 'ql' .or. trim(stage) /= 'P_OUT')) return
      write(unit, '(a,",",a,",",es24.16,",",i0,",",a,",",9(i0,","),es24.16,",",es24.16,",",es24.16,",",es24.16,",",es24.16)') &
         trim(kind), trim(this%fixed_u_trace_run_id), this%fixed_u_trace_lambda, iteration, trim(stage), &
         site, type_id, ell, spin, k, band, row, col, order, weight, energy_value, scalar, re, im
   end procedure write_fixed_u_csv_row



   module procedure write_fixed_u_trace
      integer :: unit, ios, ia, l, spin, i, j, site_type, plusbulk
      complex(rp) :: z, m1_trace, test_operator(5, 5)
      real(rp) :: center, q0, q2, test_q0_trace, test_trace_wp

      if (.not. this%fixed_u_trace_enabled .or. rank /= 0) return
      if (this%fixed_u_trace_selected_iteration > 0 .and. &
          this%fixed_u_trace_iteration /= this%fixed_u_trace_selected_iteration) return
      open(newunit=unit, file=trim(this%fixed_u_trace_filename), status='old', position='append', &
           action='write', iostat=ios)
      if (ios /= 0) return

      if (this%fixed_u_trace_ql_only) then
         if (trim(stage) == 'P_OUT') then
            do ia = 1, this%lattice%nrec
               plusbulk = this%lattice%nbulk + ia
               do l = 0, this%symbolic_atom(plusbulk)%potential%lmax
                  do spin = 1, 2
                     call this%write_fixed_u_csv_row(unit, 'ql', this%fixed_u_trace_iteration, 'P_OUT', ia - 1, &
                        this%lattice%iz(ia), l, spin, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, &
                        this%symbolic_atom(plusbulk)%potential%ql(1, l, spin), 0.0_rp, 0.0_rp)
                  end do
               end do
            end do
         end if
         close(unit)
         return
      end if

      test_operator = cmplx(0.0_rp, 0.0_rp, rp)
      if (this%fixed_u_trace_moments_only .and. this%hamiltonian%hubbard_test_operator_check) &
         call this%hamiltonian%get_hubbard_test_operator(test_operator)

      call this%write_fixed_u_csv_row(unit, 'ldm_availability', this%fixed_u_trace_iteration, trim(stage), &
         -1, -1, -1, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, &
         merge(1.0_rp, 0.0_rp, this%hamiltonian%hubbard_u_general_check), 0.0_rp, 0.0_rp)
      call this%write_fixed_u_csv_row(unit, 'vin_availability', this%fixed_u_trace_iteration, trim(stage), &
         -1, -1, -1, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, &
         merge(1.0_rp, 0.0_rp, allocated(this%hamiltonian%hubbard_u_pot)), 0.0_rp, 0.0_rp)

      do ia = 1, this%lattice%nrec
         plusbulk = this%lattice%nbulk + ia
         do l = 0, this%symbolic_atom(plusbulk)%potential%lmax
            do spin = 1, 2
               center = this%symbolic_atom(plusbulk)%potential%gravity_center(l + 1, spin) + &
                        this%symbolic_atom(plusbulk)%potential%vmad
               if (present(pl_diag_out)) then
                  if (l + 1 <= size(pl_diag_out, 2) .and. spin <= size(pl_diag_out, 3)) then
                     do i = 1, min(6, size(pl_diag_out, 4))
                        call this%write_fixed_u_csv_row(unit, 'pl_map_input', this%fixed_u_trace_iteration, &
                           'CALCULATE_PL', ia - 1, this%lattice%iz(ia), l, spin, -1, -1, -1, -1, i - 1, &
                           0.0_rp, 0.0_rp, pl_diag_out(ia, l + 1, spin, i), 0.0_rp, 0.0_rp)
                     end do
                  end if
               end if
               if (this%hamiltonian%hubbard_u_general_check) then
                  do i = 1, 2*l + 1
                     do j = 1, 2*l + 1
                        z = this%symbolic_atom(plusbulk)%potential%ldm(l + 1, spin, i, j)
                        call this%write_fixed_u_csv_row(unit, 'ldm', this%fixed_u_trace_iteration, trim(stage), ia - 1, &
                           this%lattice%iz(ia), l, spin, -1, -1, i - 1, j - 1, 0, 0.0_rp, 0.0_rp, 0.0_rp, &
                           real(z, rp), aimag(z))
                     end do
                  end do
               end if
               if (this%fixed_u_trace_moments_only .and. this%hamiltonian%hubbard_test_operator_check .and. &
                   trim(stage) == 'P_OUT' .and. l == 2) then
                  test_q0_trace = 0.0_rp
                  test_trace_wp = 0.0_rp
                  do i = 1, 5
                     do j = 1, 5
                        z = this%symbolic_atom(plusbulk)%potential%ldm(l + 1, spin, i, j)
                        call this%write_fixed_u_csv_row(unit, 'test_density_matrix', &
                           this%fixed_u_trace_iteration, 'P_OUT', ia - 1, this%lattice%iz(ia), l, spin - 1, &
                           -1, -1, i - 1, j - 1, 0, 0.0_rp, 0.0_rp, 0.0_rp, real(z, rp), aimag(z))
                        if (i == j) test_q0_trace = test_q0_trace + real(z, rp)
                        if (this%hamiltonian%hubbard_test_spin == 0 .or. &
                            this%hamiltonian%hubbard_test_spin == spin) then
                           test_trace_wp = test_trace_wp + real(test_operator(i, j) * &
                              this%symbolic_atom(plusbulk)%potential%ldm(l + 1, spin, j, i), rp)
                        end if
                     end do
                  end do
                  call this%write_fixed_u_csv_row(unit, 'test_q0_trace', this%fixed_u_trace_iteration, &
                     'P_OUT', ia - 1, this%lattice%iz(ia), l, spin - 1, -1, -1, -1, -1, 0, &
                     0.0_rp, 0.0_rp, test_q0_trace, 0.0_rp, 0.0_rp)
                  call this%write_fixed_u_csv_row(unit, 'test_trace_wp', this%fixed_u_trace_iteration, &
                     'P_OUT', ia - 1, this%lattice%iz(ia), l, spin - 1, -1, -1, -1, -1, 0, &
                     0.0_rp, 0.0_rp, test_trace_wp, 0.0_rp, 0.0_rp)
               end if
               if (present(m1_out)) then
                  if (l < size(m1_out, 2)) then
                     m1_trace = (0.0_rp, 0.0_rp)
                     do i = 1, min(2*l + 1, size(m1_out, 4))
                        do j = 1, min(2*l + 1, size(m1_out, 5))
                           z = m1_out(ia, l + 1, spin, i, j)
                           if (i == j) m1_trace = m1_trace + z
                           call this%write_fixed_u_csv_row(unit, 'ldm_m1', this%fixed_u_trace_iteration, 'M1_OUT', ia - 1, &
                              this%lattice%iz(ia), l, spin, -1, -1, i - 1, j - 1, 0, 0.0_rp, 0.0_rp, 0.0_rp, &
                              real(z, rp), aimag(z))
                        end do
                     end do
                     call this%write_fixed_u_csv_row(unit, 'm1_trace_consistency', &
                        this%fixed_u_trace_iteration, 'M1_OUT', ia - 1, this%lattice%iz(ia), l, spin, &
                        -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, &
                        real(m1_trace, rp) - this%symbolic_atom(plusbulk)%potential%ql(1, l, spin)*center, &
                        aimag(m1_trace), 0.0_rp)
                  end if
               end if
            end do
            do spin = 1, 2
               center = this%symbolic_atom(plusbulk)%potential%gravity_center(l + 1, spin) + &
                        this%symbolic_atom(plusbulk)%potential%vmad
               q0 = this%symbolic_atom(plusbulk)%potential%ql(1, l, spin)
               q2 = this%symbolic_atom(plusbulk)%potential%ql(3, l, spin)
               call this%write_fixed_u_csv_row(unit, 'ql', this%fixed_u_trace_iteration, trim(stage), ia - 1, &
                  this%lattice%iz(ia), l, spin, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, q0, 0.0_rp, 0.0_rp)
               call this%write_fixed_u_csv_row(unit, 'radial_moment', this%fixed_u_trace_iteration, trim(stage), ia - 1, &
                  this%lattice%iz(ia), l, spin, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, q0, 0.0_rp, 0.0_rp)
               call this%write_fixed_u_csv_row(unit, 'radial_moment', this%fixed_u_trace_iteration, trim(stage), ia - 1, &
                  this%lattice%iz(ia), l, spin, -1, -1, -1, -1, 1, 0.0_rp, center, q0*center, 0.0_rp, 0.0_rp)
               call this%write_fixed_u_csv_row(unit, 'radial_moment', this%fixed_u_trace_iteration, trim(stage), ia - 1, &
                  this%lattice%iz(ia), l, spin, -1, -1, -1, -1, 2, 0.0_rp, center**2 + q2, &
                  q0*(center**2 + q2), 0.0_rp, 0.0_rp)
               call this%write_fixed_u_csv_row(unit, 'radial_parameter', this%fixed_u_trace_iteration, trim(stage), ia - 1, &
                  this%lattice%iz(ia), l, spin, -1, -1, -1, -1, -1, 0.0_rp, &
                  this%symbolic_atom(plusbulk)%potential%gravity_center(l + 1, spin), 0.0_rp, 0.0_rp, 0.0_rp)
               call this%write_fixed_u_csv_row(unit, 'radial_parameter', this%fixed_u_trace_iteration, trim(stage), ia - 1, &
                  this%lattice%iz(ia), l, spin, -1, -1, -1, -1, 3, 0.0_rp, 0.0_rp, q2, 0.0_rp, 0.0_rp)
            end do
         end do
      end do

      if (allocated(this%hamiltonian%hubbard_u_pot)) then
         do site_type = 1, size(this%hamiltonian%hubbard_u_pot, 3)
            do i = 1, size(this%hamiltonian%hubbard_u_pot, 1)
               do j = 1, size(this%hamiltonian%hubbard_u_pot, 2)
                  z = this%hamiltonian%hubbard_u_pot(i, j, site_type)
                  call this%write_fixed_u_csv_row(unit, 'v_hub', this%fixed_u_trace_iteration, trim(stage), &
                     -1, site_type, -1, 0, -1, -1, i - 1, j - 1, 0, 0.0_rp, 0.0_rp, 0.0_rp, &
                     real(z, rp), aimag(z))
               end do
            end do
         end do
      end if
      if (allocated(this%hamiltonian%ee) .and. .not. this%fixed_u_trace_moments_only) then
         do site_type = 1, size(this%hamiltonian%ee, 4)
            do i = 1, size(this%hamiltonian%ee, 1)
               do j = 1, size(this%hamiltonian%ee, 2)
                  z = this%hamiltonian%ee(i, j, 1, site_type)
                  call this%write_fixed_u_csv_row(unit, 'h_total', this%fixed_u_trace_iteration, trim(stage), &
                     -1, site_type, -1, 0, -1, -1, i - 1, j - 1, 0, 0.0_rp, 0.0_rp, 0.0_rp, &
                     real(z, rp), aimag(z))
               end do
            end do
         end do
      end if
      if (stage == 'Pin_H' .and. allocated(this%hamiltonian%ee) .and. .not. this%fixed_u_trace_moments_only) &
         call this%write_fixed_u_rs_bonds(unit)
      if (stage == 'P_OUT' .and. .not. this%use_kspace) call this%write_fixed_u_rs_green(unit)
      close(unit)
   end procedure write_fixed_u_trace


   module procedure write_fixed_u_rs_bonds
      integer :: ntype, isite, ia, ja, jsite, iclu, ineigh, nr, i, j, coord
      complex(rp) :: z
      real(rp) :: r_direct(3)

      ! Export fractional direct vectors without invoking clusba: clusba
      ! updates lattice%sbarvec, so even a diagnostic trace must not call it.

      if (.not. allocated(this%lattice%irec) .or. .not. allocated(this%lattice%iz)) return
      if (size(this%lattice%irec) < this%lattice%nrec) return
      if (.not. allocated(this%lattice%atlist) .or. .not. allocated(this%lattice%nn)) return

      do ntype = 1, this%lattice%ntype
         if (ntype > size(this%lattice%atlist)) cycle
         ia = this%lattice%atlist(ntype)
         if (ia < 1 .or. ia > size(this%lattice%nn, 1)) cycle
         if (size(this%lattice%nn, 2) < 1) cycle
         nr = min(this%lattice%nn(ia, 1), size(this%lattice%nn, 2), size(this%hamiltonian%ee, 3))
         do ineigh = 1, nr
            do i = 1, size(this%hamiltonian%ee, 1)
               do j = 1, size(this%hamiltonian%ee, 2)
                  z = this%hamiltonian%ee(i, j, ineigh, ntype)
                  call this%write_fixed_u_csv_row(unit, 'rs_hopping', this%fixed_u_trace_iteration, 'RS_OPERATOR', &
                     -1, ntype - 1, -1, -1, -1, ineigh - 1, i - 1, j - 1, 0, 0.0_rp, 0.0_rp, 0.0_rp, &
                     real(z, rp), aimag(z))
               end do
            end do
         end do
      end do

      do isite = 1, this%lattice%nrec
         iclu = this%lattice%irec(isite)
         if (iclu < 1 .or. iclu > size(this%lattice%iz)) cycle
         ntype = this%lattice%iz(iclu)
         if (ntype < 1 .or. ntype > this%lattice%ntype) cycle
         if (ntype > size(this%lattice%atlist)) cycle
         ia = this%lattice%atlist(ntype)
         if (ia < 1 .or. ia > size(this%lattice%nn, 1)) cycle
         if (size(this%lattice%nn, 2) < 1) cycle
         nr = min(this%lattice%nn(ia, 1), size(this%lattice%nn, 2))
         do ineigh = 1, nr
            if (ineigh == 1) then
               jsite = isite
               r_direct = 0.0_rp
            else
               ja = this%lattice%nn(ia, ineigh)
               if (ja < 1 .or. ja > this%lattice%kk) cycle
               jsite = this%lattice%iz(ja)
               if (jsite < 1 .or. jsite > this%lattice%nrec) cycle
               if (this%lattice%a_cart_inv_ready) then
                  r_direct = matmul(this%lattice%a_cart_inv, &
                                    this%lattice%cr(:, ja) - this%lattice%cr(:, ia))
               else
                  r_direct = matmul(inverse_3x3(this%lattice%a), &
                                    this%lattice%cr(:, ja) - this%lattice%cr(:, ia))
               end if
            end if
            do coord = 1, 3
               call this%write_fixed_u_csv_row(unit, 'rs_bond_map', this%fixed_u_trace_iteration, 'RS_OPERATOR', &
                  isite - 1, ntype - 1, -1, -1, -1, ineigh - 1, -1, jsite - 1, coord - 1, 0.0_rp, &
                  r_direct(coord), 0.0_rp, 0.0_rp, 0.0_rp)
            end do
         end do
      end do
   end procedure write_fixed_u_rs_bonds


   module procedure write_fixed_u_radial_trace
      integer :: unit, ios, ir, spin

      if (.not. this%fixed_u_trace_enabled .or. rank /= 0) return
      if (this%fixed_u_trace_selected_iteration > 0 .and. &
          this%fixed_u_trace_iteration /= this%fixed_u_trace_selected_iteration) return
      open(newunit=unit, file=trim(this%fixed_u_trace_filename), status='old', position='append', &
           action='write', iostat=ios)
      if (ios /= 0) return
      do ir = 1, size(rofi)
         do spin = 1, size(v, 2)
            call this%write_fixed_u_csv_row(unit, 'radial_v', this%fixed_u_trace_iteration, 'POST_ATOMSC', &
               isite, this%lattice%iz(isite), -1, spin - 1, -1, -1, ir - 1, spin - 1, 0, 0.0_rp, &
               rofi(ir), v(ir, spin), 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'radial_density_in', this%fixed_u_trace_iteration, 'POST_ATOMSC', &
               isite, this%lattice%iz(isite), -1, spin - 1, -1, -1, ir - 1, spin - 1, 0, 0.0_rp, &
               rofi(ir), rho_in(ir, spin), 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'radial_density_out', this%fixed_u_trace_iteration, 'POST_ATOMSC', &
               isite, this%lattice%iz(isite), -1, spin - 1, -1, -1, ir - 1, spin - 1, 0, 0.0_rp, &
               rofi(ir), rho_out(ir, spin), 0.0_rp, 0.0_rp)
         end do
      end do
      call energy_row('sum_ev', atom%potential%SUMEV)
      call energy_row('sum_ec', atom%potential%SUMEC)
      call energy_row('ekin', atom%potential%EKIN)
      call energy_row('utot', atom%potential%UTOT)
      call energy_row('rhoeps', atom%potential%RHOEPS)
      call energy_row('etot_site', atom%potential%ETOT)
      close(unit)
   contains
      subroutine energy_row(label, energy_value)
         character(len=*), intent(in) :: label
         real(rp), intent(in) :: energy_value
         call this%write_fixed_u_csv_row(unit, label, this%fixed_u_trace_iteration, 'POST_ATOMSC', &
            isite, this%lattice%iz(isite), -1, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, energy_value, 0.0_rp, 0.0_rp)
      end subroutine energy_row
   end procedure write_fixed_u_radial_trace


   module procedure write_fixed_u_radial_replay_input
      integer :: unit, ios, ir, spin, ell, core_index

      if (.not. this%fixed_u_trace_enabled .or. rank /= 0) return
      if (this%fixed_u_trace_selected_iteration > 0 .and. &
          this%fixed_u_trace_iteration /= this%fixed_u_trace_selected_iteration) return
      open(newunit=unit, file=trim(this%fixed_u_trace_filename), status='old', position='append', &
           action='write', iostat=ios)
      if (ios /= 0) return
      call meta('radial_A', atom%A)
      call meta('radial_B', atom%B())
      call meta('radial_Z', atom%element%atomic_number)
      call meta('radial_tolerance', tol)
      call meta('radial_VHRMAX', vhrmax)
      call meta('radial_B_fsm', b_fsm)
      call meta('radial_txc', real(txc, rp))
      call meta('radial_hyperfine', merge(1.0_rp, 0.0_rp, hyperfine))
      call meta('radial_NR', real(size(rofi), rp))
      call meta('radial_lmax', real(atom%potential%lmax, rp))
      call meta('radial_nsp', real(size(v, 2), rp))
      call meta('radial_f_core', real(atom%element%f_core, rp))
      call meta('radial_replay_vmad', atom%potential%vmad)
      do ir = 1, size(rofi)
         do spin = 1, size(v, 2)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_mesh', this%fixed_u_trace_iteration, 'RADIAL_REPLAY', &
               this%fixed_u_trace_current_site - 1, -1, -1, spin - 1, -1, -1, ir - 1, spin - 1, 0, &
               0.0_rp, rofi(ir), 0.0_rp, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_v', this%fixed_u_trace_iteration, 'RADIAL_REPLAY', &
               this%fixed_u_trace_current_site - 1, -1, -1, spin - 1, -1, -1, ir - 1, spin - 1, 0, &
               0.0_rp, rofi(ir), v(ir, spin), 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_rho_in', this%fixed_u_trace_iteration, 'RADIAL_REPLAY', &
               this%fixed_u_trace_current_site - 1, -1, -1, spin - 1, -1, -1, ir - 1, spin - 1, 0, &
               0.0_rp, rofi(ir), rho_in(ir, spin), 0.0_rp, 0.0_rp)
         end do
      end do
      do ell = 0, atom%potential%lmax
         do spin = 1, size(pl, 2)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_pl', this%fixed_u_trace_iteration, 'RADIAL_REPLAY', &
               this%fixed_u_trace_current_site - 1, -1, ell, spin - 1, -1, -1, ell, spin - 1, 0, &
               0.0_rp, 0.0_rp, pl(ell + 1, spin), 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_pnu', this%fixed_u_trace_iteration, 'RADIAL_REPLAY', &
               this%fixed_u_trace_current_site - 1, -1, ell, spin - 1, -1, -1, ell, spin - 1, 0, &
               0.0_rp, 0.0_rp, atom%potential%pnu(ell, spin), 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_enu', this%fixed_u_trace_iteration, 'RADIAL_REPLAY', &
               this%fixed_u_trace_current_site - 1, -1, ell, spin - 1, -1, -1, ell, spin - 1, 0, &
               0.0_rp, 0.0_rp, atom%potential%enu(ell, spin), 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_gravity_center', this%fixed_u_trace_iteration, &
               'RADIAL_REPLAY', this%fixed_u_trace_current_site - 1, -1, ell, spin - 1, -1, -1, ell, spin - 1, 0, &
               0.0_rp, 0.0_rp, atom%potential%gravity_center(ell + 1, spin), 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_band_center', this%fixed_u_trace_iteration, &
               'RADIAL_REPLAY', this%fixed_u_trace_current_site - 1, -1, ell, spin - 1, -1, -1, ell, spin - 1, 0, &
               0.0_rp, 0.0_rp, atom%potential%gravity_center(ell + 1, spin) + atom%potential%vmad, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_c', this%fixed_u_trace_iteration, 'RADIAL_REPLAY', &
               this%fixed_u_trace_current_site - 1, -1, ell, spin - 1, -1, -1, ell, spin - 1, 0, &
               0.0_rp, 0.0_rp, atom%potential%c(ell, spin), 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_qpar', this%fixed_u_trace_iteration, 'RADIAL_REPLAY', &
               this%fixed_u_trace_current_site - 1, -1, ell, spin - 1, -1, -1, ell, spin - 1, 0, &
               0.0_rp, 0.0_rp, atom%potential%qpar(ell, spin), 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_srdel', this%fixed_u_trace_iteration, 'RADIAL_REPLAY', &
               this%fixed_u_trace_current_site - 1, -1, ell, spin - 1, -1, -1, ell, spin - 1, 0, &
               0.0_rp, 0.0_rp, atom%potential%srdel(ell, spin), 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_ql', this%fixed_u_trace_iteration, 'RADIAL_REPLAY', &
               this%fixed_u_trace_current_site - 1, -1, ell, spin - 1, -1, -1, ell, spin - 1, 0, &
               0.0_rp, 0.0_rp, ql(1, ell + 1, spin), 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_ql', this%fixed_u_trace_iteration, 'RADIAL_REPLAY', &
               this%fixed_u_trace_current_site - 1, -1, ell, spin - 1, -1, -1, ell, spin - 1, 1, &
               0.0_rp, 0.0_rp, ql(2, ell + 1, spin), 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_ql', this%fixed_u_trace_iteration, 'RADIAL_REPLAY', &
               this%fixed_u_trace_current_site - 1, -1, ell, spin - 1, -1, -1, ell, spin - 1, 2, &
               0.0_rp, 0.0_rp, ql(3, ell + 1, spin), 0.0_rp, 0.0_rp)
         end do
      end do
      if (size(ec_seed) < ncore) error stop 'fixed-U replay capture has incomplete core-energy seeds'
      do core_index = 1, ncore
         call this%write_fixed_u_csv_row(unit, 'radial_replay_ec_seed', this%fixed_u_trace_iteration, &
            'RADIAL_REPLAY', this%fixed_u_trace_current_site - 1, -1, -1, -1, -1, core_index - 1, &
            core_index - 1, -1, 0, 0.0_rp, 0.0_rp, ec_seed(core_index), 0.0_rp, 0.0_rp)
      end do
      close(unit)
   contains
      subroutine meta(kind, value)
         character(len=*), intent(in) :: kind
         real(rp), intent(in) :: value
         call this%write_fixed_u_csv_row(unit, kind, this%fixed_u_trace_iteration, 'RADIAL_REPLAY', &
            this%fixed_u_trace_current_site - 1, -1, -1, -1, -1, -1, -1, -1, 0, &
            0.0_rp, 0.0_rp, value, 0.0_rp, 0.0_rp)
      end subroutine meta
   end procedure write_fixed_u_radial_replay_input


   module procedure write_fixed_u_radial_replay_output
      integer :: unit, ios, ir, spin

      if (.not. this%fixed_u_trace_enabled .or. rank /= 0) return
      if (this%fixed_u_trace_selected_iteration > 0 .and. &
          this%fixed_u_trace_iteration /= this%fixed_u_trace_selected_iteration) return
      open(newunit=unit, file=trim(this%fixed_u_trace_filename), status='old', position='append', &
           action='write', iostat=ios)
      if (ios /= 0) return
      do ir = 1, size(rho, 1)
         do spin = 1, size(rho, 2)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_rho_out', this%fixed_u_trace_iteration, &
               'RADIAL_REPLAY', this%fixed_u_trace_current_site - 1, -1, -1, spin - 1, -1, -1, &
               ir - 1, spin - 1, 0, 0.0_rp, 0.0_rp, rho(ir, spin), 0.0_rp, 0.0_rp)
         end do
      end do
      do spin = 1, size(rho, 2)
         call this%write_fixed_u_csv_row(unit, 'radial_replay_sumec', this%fixed_u_trace_iteration, &
            'RADIAL_REPLAY', this%fixed_u_trace_current_site - 1, -1, -1, spin - 1, -1, -1, &
            -1, -1, 0, 0.0_rp, 0.0_rp, sumec(spin), 0.0_rp, 0.0_rp)
         call this%write_fixed_u_csv_row(unit, 'radial_replay_sumev', this%fixed_u_trace_iteration, &
            'RADIAL_REPLAY', this%fixed_u_trace_current_site - 1, -1, -1, spin - 1, -1, -1, &
            -1, -1, 0, 0.0_rp, 0.0_rp, sumev(spin), 0.0_rp, 0.0_rp)
      end do
      close(unit)
   end procedure write_fixed_u_radial_replay_output


   module procedure verify_fixed_u_radial_replay
      type(self) :: replay_state
      type(symbolic_atom) :: replay_atom
      real(rp), allocatable :: rho_replay(:, :), sumec_replay(:), sumev_replay(:), ec_replay(:), ev_replay(:)
      integer :: unit, ios
      real(rp) :: rho_error, sumec_error, sumev_error, ev_error, ec_error

      if (.not. this%fixed_u_trace_enabled .or. rank /= 0) return
      if (this%fixed_u_trace_selected_iteration > 0 .and. &
          this%fixed_u_trace_iteration /= this%fixed_u_trace_selected_iteration) return
      if (size(ec_seed) < ncore .or. size(ev_seed) < nval .or. size(ec_actual) < ncore .or. size(ev_actual) < nval) &
         error stop 'verify_fixed_u_radial_replay: inconsistent seed/output sizes'

      ! Repeat the unchanged production NEWRHO call against the exact accepted
      ! atom and input arrays.  The independent scratch object and atom copy
      ! keep this diagnostic replay from changing the live SCF state.
      replay_atom = atom
      replay_state%control => this%control
      allocate(replay_state%vtn(nr, 2), replay_state%vzt(nr, 2), &
               replay_state%fun2(nr, lmax + 1, 2), replay_state%phi_amp(nr, lmax + 1, 2))
      replay_state%vtn = 0.0_rp
      replay_state%vzt = 0.0_rp
      replay_state%fun2 = 0.0_rp
      replay_state%phi_amp = 0.0_rp
      allocate(rho_replay(nr, 2), sumec_replay(2), sumev_replay(2), &
               ec_replay(max(ncore, 1)), ev_replay(max(nval, 1)))
      rho_replay = 0.0_rp
      sumec_replay = 0.0_rp
      sumev_replay = 0.0_rp
      if (ncore > 0) ec_replay(1:ncore) = ec_seed(1:ncore)
      if (nval > 0) ev_replay(1:nval) = ev_seed(1:nval)
      call replay_state%newrho(replay_atom, z, lmax, a, b, nr, rofi, v, rho_replay, pl, ql, &
                               sumec_replay, sumev_replay, ec_replay, ev_replay, tol, 2, 0)
      rho_error = maxval(abs(rho_replay - rho_actual))
      sumec_error = maxval(abs(sumec_replay - sumec_actual))
      sumev_error = maxval(abs(sumev_replay - sumev_actual))
      ev_error = 0.0_rp
      ec_error = 0.0_rp
      if (nval > 0) ev_error = maxval(abs(ev_replay(1:nval) - ev_actual(1:nval)))
      if (ncore > 0) ec_error = maxval(abs(ec_replay(1:ncore) - ec_actual(1:ncore)))

      deallocate(replay_state%vtn, replay_state%vzt, replay_state%fun2, replay_state%phi_amp)

      open(newunit=unit, file=trim(this%fixed_u_trace_filename), status='old', position='append', &
           action='write', iostat=ios)
      if (ios == 0) then
         call result('rho_max_abs', rho_error, 0)
         call result('sumec_max_abs_ry', sumec_error, 1)
         call result('sumev_max_abs_ry', sumev_error, 2)
         call result('ec_max_abs_ry', ec_error, 3)
         call result('ev_max_abs_ry', ev_error, 4)
         close(unit)
      end if
      call this%replay_fixed_u_m0_moments(atom, z, lmax, a, b, nr, rofi, v, rho_in, rho_actual, &
         pl, ql, ec_seed, ev_seed, ncore, nval, tol, rofi_base, v_base, vhrmax, b_fsm, &
         rvh_actual, reps_actual, rmu_actual, vnucl_actual, sumec_actual, sumev_actual)
   contains
      subroutine result(label, value, result_index)
         character(len=*), intent(in) :: label
         real(rp), intent(in) :: value
         integer, intent(in) :: result_index
         call this%write_fixed_u_csv_row(unit, 'radial_replay_'//trim(label), this%fixed_u_trace_iteration, &
            'RADIAL_REPLAY', this%fixed_u_trace_current_site - 1, -1, -1, -1, -1, -1, -1, -1, &
            result_index, 0.0_rp, 0.0_rp, value, 0.0_rp, 0.0_rp)
      end subroutine result
   end procedure verify_fixed_u_radial_replay


   module procedure snapshot_fixed_u_m0_moments
      complex(rp), allocatable :: density_override(:, :, :, :, :), vout(:, :, :)
      integer :: ia, plusbulk, itype, iclu, spin, i, j, k, rowoff, norb_l, unit, ios
      complex(rp) :: p_ij, t_complex, vm1_complex, v2p_complex, v2_ij, p_trace, m1_trace
      real(rp), allocatable :: im_resid(:, :), p_trace_error(:, :), m1_trace_error(:, :)
      real(rp) :: q0, center, q2, t, local_q2, t_vout, vm1_vout, v2_vout
      real(rp) :: t_pt, vm1_pt, v2_pt
      logical :: p_available

      this%fixed_u_m0_snapshot_valid = .false.
      if (.not. this%fixed_u_trace_enabled .or. rank /= 0) return
      if (this%fixed_u_trace_selected_iteration > 0 .and. &
          this%fixed_u_trace_iteration /= this%fixed_u_trace_selected_iteration) return
      if (size(m1, 1) /= this%lattice%nrec .or. size(pl_diag, 1) /= this%lattice%nrec .or. &
          size(vin, 1) < nb .or. size(vin, 2) < nb) return
      if (lmax_basis < 2) return
      if (.not. allocated(this%lattice%irec) .or. .not. allocated(this%lattice%iz)) return
      if (size(this%lattice%irec) < this%lattice%nrec) return

      if (allocated(this%fixed_u_m0_q0)) deallocate(this%fixed_u_m0_q0)
      if (allocated(this%fixed_u_m0_center)) deallocate(this%fixed_u_m0_center)
      if (allocated(this%fixed_u_m0_q2)) deallocate(this%fixed_u_m0_q2)
      if (allocated(this%fixed_u_m0_local_center)) deallocate(this%fixed_u_m0_local_center)
      if (allocated(this%fixed_u_m0_local_q2)) deallocate(this%fixed_u_m0_local_q2)
      if (allocated(this%fixed_u_m0_trvp)) deallocate(this%fixed_u_m0_trvp)
      if (allocated(this%fixed_u_m0_vout_center)) deallocate(this%fixed_u_m0_vout_center)
      if (allocated(this%fixed_u_m0_vout_q2)) deallocate(this%fixed_u_m0_vout_q2)
      if (allocated(this%fixed_u_m0_pt_center)) deallocate(this%fixed_u_m0_pt_center)
      if (allocated(this%fixed_u_m0_pt_q2)) deallocate(this%fixed_u_m0_pt_q2)
      if (allocated(this%fixed_u_m0_t_vout)) deallocate(this%fixed_u_m0_t_vout)
      if (allocated(this%fixed_u_m0_t_pt)) deallocate(this%fixed_u_m0_t_pt)
      if (allocated(this%fixed_u_m0_ql)) deallocate(this%fixed_u_m0_ql)
      if (allocated(this%fixed_u_m0_pl_diag)) deallocate(this%fixed_u_m0_pl_diag)
      if (allocated(this%fixed_u_m0_vin)) deallocate(this%fixed_u_m0_vin)
      allocate(this%fixed_u_m0_q0(this%lattice%nrec, 2), this%fixed_u_m0_center(this%lattice%nrec, 2), &
               this%fixed_u_m0_q2(this%lattice%nrec, 2), this%fixed_u_m0_local_center(this%lattice%nrec, 2), &
               this%fixed_u_m0_local_q2(this%lattice%nrec, 2), this%fixed_u_m0_trvp(this%lattice%nrec, 2), &
               this%fixed_u_m0_vout_center(this%lattice%nrec, 2), this%fixed_u_m0_vout_q2(this%lattice%nrec, 2), &
               this%fixed_u_m0_pt_center(this%lattice%nrec, 2), this%fixed_u_m0_pt_q2(this%lattice%nrec, 2), &
               this%fixed_u_m0_t_vout(this%lattice%nrec, 2), this%fixed_u_m0_t_pt(this%lattice%nrec, 2), &
               this%fixed_u_m0_ql(this%lattice%nrec, 3, lmax_basis + 1, 2), &
               this%fixed_u_m0_pl_diag(this%lattice%nrec, lmax_basis + 1, 2, 6), &
               this%fixed_u_m0_vin(size(vin, 1), size(vin, 2), size(vin, 3)))
      this%fixed_u_m0_ql = 0.0_rp
      this%fixed_u_m0_pl_diag = pl_diag
      this%fixed_u_m0_vin = vin
      this%fixed_u_m0_vout_center = 0.0_rp
      this%fixed_u_m0_vout_q2 = 0.0_rp
      this%fixed_u_m0_pt_center = 0.0_rp
      this%fixed_u_m0_pt_q2 = 0.0_rp
      this%fixed_u_m0_t_vout = 0.0_rp
      this%fixed_u_m0_t_pt = 0.0_rp
      allocate(im_resid(this%lattice%nrec, 2))
      allocate(p_trace_error(this%lattice%nrec, 2), m1_trace_error(this%lattice%nrec, 2))
      im_resid = 0.0_rp
      p_trace_error = 0.0_rp
      m1_trace_error = 0.0_rp
      p_available = this%hamiltonian%hubbard_u_general_check
      if (p_available) then
         do ia = 1, this%lattice%nrec
            plusbulk = this%lattice%nbulk + ia
            if (plusbulk < 1 .or. plusbulk > size(this%symbolic_atom)) then
               p_available = .false.
               exit
            end if
            if (.not. allocated(this%symbolic_atom(plusbulk)%potential%ldm)) then
               p_available = .false.
               exit
            end if
            if (size(this%symbolic_atom(plusbulk)%potential%ldm, 1) < 3 .or. &
                size(this%symbolic_atom(plusbulk)%potential%ldm, 2) < 2 .or. &
                size(this%symbolic_atom(plusbulk)%potential%ldm, 3) < 5 .or. &
                size(this%symbolic_atom(plusbulk)%potential%ldm, 4) < 5) then
               p_available = .false.
               exit
            end if
         end do
      end if
      this%fixed_u_m0_p_available = p_available

      ! Evaluate V(Pout) without changing the live Hamiltonian.  This is a
      ! deliberately wrong-pairing control for the frozen radial replay.
      allocate(density_override(this%lattice%ntype, 4, 2, 7, 7), vout(nb, nb, this%lattice%ntype))
      density_override = (0.0_rp, 0.0_rp)
      vout = (0.0_rp, 0.0_rp)
      if (p_available) then
         do ia = 1, this%lattice%nrec
            plusbulk = this%lattice%nbulk + ia
            iclu = this%lattice%irec(ia)
            if (iclu < 1 .or. iclu > size(this%lattice%iz)) then
               deallocate(im_resid, p_trace_error, m1_trace_error)
               deallocate(density_override, vout)
               return
            end if
            itype = this%lattice%iz(iclu)
            if (itype < 1 .or. itype > this%lattice%ntype) then
               deallocate(im_resid, p_trace_error, m1_trace_error)
               deallocate(density_override, vout)
               return
            end if
            do spin = 1, 2
               do i = 1, min(5, size(this%symbolic_atom(plusbulk)%potential%ldm, 3))
                  do j = 1, min(5, size(this%symbolic_atom(plusbulk)%potential%ldm, 4))
                     density_override(itype, 3, spin, i, j) = &
                        this%symbolic_atom(plusbulk)%potential%ldm(3, spin, i, j)
                  end do
               end do
            end do
         end do
         call this%hamiltonian%calculate_hubbard_u_potential_general( &
            density_override=density_override, potential_out=vout)
      end if

      do ia = 1, this%lattice%nrec
         plusbulk = this%lattice%nbulk + ia
         iclu = this%lattice%irec(ia)
         if (iclu < 1 .or. iclu > size(this%lattice%iz)) then
            deallocate(im_resid, p_trace_error, m1_trace_error)
            deallocate(density_override, vout)
            return
         end if
         itype = this%lattice%iz(iclu)
         if (itype < 1 .or. itype > size(vin, 3)) return
         do spin = 1, 2
            do i = 1, 3
               this%fixed_u_m0_ql(ia, i, 1:lmax_basis + 1, spin) = &
                  this%symbolic_atom(plusbulk)%potential%ql(i, 0:lmax_basis, spin)
            end do
            q0 = this%symbolic_atom(plusbulk)%potential%ql(1, 2, spin)
            center = this%symbolic_atom(plusbulk)%potential%gravity_center(3, spin) + &
                     this%symbolic_atom(plusbulk)%potential%vmad
            q2 = this%symbolic_atom(plusbulk)%potential%ql(3, 2, spin)
            if (q0 <= tiny(1.0_rp)) then
               this%fixed_u_m0_snapshot_valid = .false.
               return
            end if
            this%fixed_u_m0_q0(ia, spin) = q0
            this%fixed_u_m0_center(ia, spin) = center
            this%fixed_u_m0_q2(ia, spin) = q2
            t_complex = (0.0_rp, 0.0_rp)
            vm1_complex = (0.0_rp, 0.0_rp)
            v2p_complex = (0.0_rp, 0.0_rp)
            p_trace = (0.0_rp, 0.0_rp)
            m1_trace = (0.0_rp, 0.0_rp)
            rowoff = 4 + spin_off*(spin - 1)
            norb_l = 5
            do i = 1, norb_l
               m1_trace = m1_trace + m1(ia, 3, spin, i, i)
            end do
            if (p_available) then
            do i = 1, norb_l
               do j = 1, norb_l
                  p_ij = this%symbolic_atom(plusbulk)%potential%ldm(3, spin, i, j)
                  if (i == j) p_trace = p_trace + p_ij
                  t_complex = t_complex + vin(rowoff + i, rowoff + j, itype)* &
                     this%symbolic_atom(plusbulk)%potential%ldm(3, spin, j, i)
                  vm1_complex = vm1_complex + vin(rowoff + i, rowoff + j, itype)* &
                     (m1(ia, 3, spin, j, i) - center*this%symbolic_atom(plusbulk)%potential%ldm(3, spin, j, i))
                  v2_ij = (0.0_rp, 0.0_rp)
                  do k = 1, norb_l
                     v2_ij = v2_ij + vin(rowoff + i, rowoff + k, itype)*vin(rowoff + k, rowoff + j, itype)
                  end do
                  v2p_complex = v2p_complex + v2_ij*this%symbolic_atom(plusbulk)%potential%ldm(3, spin, j, i)
               end do
            end do
            end if
            t = real(t_complex, rp)
            local_q2 = q2 - 2.0_rp*real(vm1_complex, rp) + real(v2p_complex, rp) - t*t/q0
            im_resid(ia, spin) = max(abs(aimag(t_complex)), abs(aimag(vm1_complex)), abs(aimag(v2p_complex)))
            p_trace_error(ia, spin) = real(p_trace, rp) - q0
            m1_trace_error(ia, spin) = real(m1_trace, rp) - q0*center
            this%fixed_u_m0_trvp(ia, spin) = t
            this%fixed_u_m0_local_center(ia, spin) = center - t/q0
            this%fixed_u_m0_local_q2(ia, spin) = local_q2

            t_vout = 0.0_rp
            vm1_vout = 0.0_rp
            v2_vout = 0.0_rp
            t_pt = 0.0_rp
            vm1_pt = 0.0_rp
            v2_pt = 0.0_rp
            if (p_available) then
            do i = 1, norb_l
               do j = 1, norb_l
                  t_vout = t_vout + real(vout(rowoff+i, rowoff+j, itype)* &
                     this%symbolic_atom(plusbulk)%potential%ldm(3, spin, j, i), rp)
                  vm1_vout = vm1_vout + real(vout(rowoff+i, rowoff+j, itype)* &
                     (m1(ia, 3, spin, j, i) - center*this%symbolic_atom(plusbulk)%potential%ldm(3, spin, j, i)), rp)
                  t_pt = t_pt + real(vin(rowoff+i, rowoff+j, itype)* &
                     this%symbolic_atom(plusbulk)%potential%ldm(3, spin, i, j), rp)
                  vm1_pt = vm1_pt + real(vin(rowoff+i, rowoff+j, itype)* &
                     (m1(ia, 3, spin, j, i) - center*this%symbolic_atom(plusbulk)%potential%ldm(3, spin, i, j)), rp)
                  v2_ij = (0.0_rp, 0.0_rp)
                  do k = 1, norb_l
                     v2_ij = v2_ij + vout(rowoff+i, rowoff+k, itype)*vout(rowoff+k, rowoff+j, itype)
                  end do
                  v2_vout = v2_vout + real(v2_ij*this%symbolic_atom(plusbulk)%potential%ldm(3, spin, j, i), rp)
                  v2_ij = (0.0_rp, 0.0_rp)
                  do k = 1, norb_l
                     v2_ij = v2_ij + vin(rowoff+i, rowoff+k, itype)*vin(rowoff+k, rowoff+j, itype)
                  end do
                  v2_pt = v2_pt + real(v2_ij*this%symbolic_atom(plusbulk)%potential%ldm(3, spin, i, j), rp)
               end do
            end do
            end if
            this%fixed_u_m0_vout_center(ia, spin) = center - t_vout/q0
            this%fixed_u_m0_vout_q2(ia, spin) = q2 - 2.0_rp*vm1_vout + v2_vout - t_vout*t_vout/q0
            this%fixed_u_m0_pt_center(ia, spin) = center - t_pt/q0
            this%fixed_u_m0_pt_q2(ia, spin) = q2 - 2.0_rp*vm1_pt + v2_pt - t_pt*t_pt/q0
            this%fixed_u_m0_t_vout(ia, spin) = t_vout
            this%fixed_u_m0_t_pt(ia, spin) = t_pt
         end do
      end do
      this%fixed_u_m0_snapshot_iteration = this%fixed_u_trace_iteration
      this%fixed_u_m0_snapshot_valid = .true.

      open(newunit=unit, file=trim(this%fixed_u_trace_filename), status='old', position='append', &
           action='write', iostat=ios)
      if (ios /= 0) return
      do ia = 1, this%lattice%nrec
         do spin = 1, 2
            itype = this%lattice%ib(ia)
            call this%write_fixed_u_csv_row(unit, 'm0_carrier_input', this%fixed_u_m0_snapshot_iteration, &
               'FROZEN_REPLAY', ia - 1, itype - 1, 2, spin, -1, -1, -1, -1, 0, 0.0_rp, &
               this%fixed_u_m0_center(ia, spin), this%fixed_u_m0_q0(ia, spin), &
               this%fixed_u_m0_q2(ia, spin), 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'm0_p_availability', this%fixed_u_m0_snapshot_iteration, &
               'FROZEN_REPLAY', ia - 1, itype - 1, 2, spin, -1, -1, -1, -1, 0, 0.0_rp, &
               merge(1.0_rp, 0.0_rp, p_available), 0.0_rp, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'm0_carrier_local', this%fixed_u_m0_snapshot_iteration, &
               'FROZEN_REPLAY', ia - 1, itype - 1, 2, spin, -1, -1, -1, -1, 1, &
               this%fixed_u_m0_trvp(ia, spin), this%fixed_u_m0_local_center(ia, spin), &
               this%fixed_u_m0_q0(ia, spin), this%fixed_u_m0_local_q2(ia, spin), 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'm0_carrier_control', this%fixed_u_m0_snapshot_iteration, &
               'VOUT_CONTROL', ia - 1, itype - 1, 2, spin, -1, -1, -1, -1, 0, 0.0_rp, &
               this%fixed_u_m0_t_vout(ia, spin), this%fixed_u_m0_vout_center(ia, spin), &
               this%fixed_u_m0_q0(ia, spin), this%fixed_u_m0_vout_q2(ia, spin))
            call this%write_fixed_u_csv_row(unit, 'm0_carrier_control', this%fixed_u_m0_snapshot_iteration, &
               'P_TRANSPOSE_CONTROL', ia - 1, itype - 1, 2, spin, -1, -1, -1, -1, 0, 0.0_rp, &
               this%fixed_u_m0_t_pt(ia, spin), this%fixed_u_m0_pt_center(ia, spin), &
               this%fixed_u_m0_q0(ia, spin), this%fixed_u_m0_pt_q2(ia, spin))
            call this%write_fixed_u_csv_row(unit, 'm0_carrier_identity', this%fixed_u_m0_snapshot_iteration, &
               'FROZEN_REPLAY', ia - 1, itype - 1, 2, spin, -1, -1, -1, -1, 2, &
               0.0_rp, 0.0_rp, im_resid(ia, spin), 0.0_rp, 0.0_rp)
            if (p_available) then
               call this%write_fixed_u_csv_row(unit, 'm0_trace_consistency', this%fixed_u_m0_snapshot_iteration, &
                  'FROZEN_REPLAY', ia - 1, itype - 1, 2, spin, -1, -1, -1, -1, 0, &
                  0.0_rp, 0.0_rp, p_trace_error(ia, spin), 0.0_rp, 0.0_rp)
            else
               call this%write_fixed_u_csv_row(unit, 'm0_p_trace_unavailable', &
                  this%fixed_u_m0_snapshot_iteration, 'FROZEN_REPLAY', ia - 1, itype - 1, 2, spin, &
                  -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp)
            end if
            call this%write_fixed_u_csv_row(unit, 'm0_trace_consistency', this%fixed_u_m0_snapshot_iteration, &
               'FROZEN_REPLAY', ia - 1, itype - 1, 2, spin, -1, -1, -1, -1, 1, &
               0.0_rp, 0.0_rp, m1_trace_error(ia, spin), 0.0_rp, 0.0_rp)
         end do
      end do
      close(unit)
      deallocate(im_resid, p_trace_error, m1_trace_error)
      deallocate(density_override, vout)
   end procedure snapshot_fixed_u_m0_moments


   module procedure replay_fixed_u_m0_moments
      type(bands) :: pl_bands
      type(symbolic_atom), allocatable, target :: pl_atoms(:)
      real(rp), allocatable :: ql_raw(:, :, :), ql_local(:, :, :), ql_vout(:, :, :), ql_pt(:, :, :)
      real(rp), allocatable :: rho_case(:, :), sumec_case(:), sumev_case(:)
      real(rp), allocatable :: ec_case(:), ev_case(:)
      real(rp), allocatable :: pl_raw(:, :), pl_local(:, :), pl_a11(:, :, :)
      real(rp), allocatable :: baseline_ev(:)
      integer :: ia, plusbulk, l, spin, unit, ios, row, iv, nshift
      real(rp) :: pl_reproduction_error, ql_vs_live_error, q2_min
      real(rp) :: baseline_sumev(2), baseline_sumec(2)
      real(rp) :: component_snapshot(10), actual_component(10), output_component(10)
      real(rp) :: actual_rhovh, actual_rhomu
      real(rp) :: actual_rhoeps, actual_zvnucl, actual_ekin, actual_utot, actual_etot
      logical :: have_raw_energy

      if (.not. this%fixed_u_m0_snapshot_valid) return
      if (this%fixed_u_m0_snapshot_iteration /= this%fixed_u_trace_iteration) return
      if (this%fixed_u_trace_current_site < 1 .or. this%fixed_u_trace_current_site > this%lattice%nrec) return
      if (.not. allocated(this%fixed_u_m0_pl_diag) .or. .not. allocated(this%fixed_u_m0_ql)) return
      if (size(ec_seed) < ncore .or. size(ev_seed) < nval) return

      allocate(pl_atoms(size(this%symbolic_atom)))
      pl_atoms = this%symbolic_atom
      pl_bands%symbolic_atom => pl_atoms
      pl_bands%lattice => this%lattice
      pl_bands%control => this%control
      call load_captured_pl_inputs(pl_atoms)
      call pl_bands%calculate_pl()
      allocate(pl_raw(lmax + 1, 2), pl_local(lmax + 1, 2))
      ia = this%fixed_u_trace_current_site
      plusbulk = this%lattice%nbulk + ia
      do spin = 1, 2
         do l = 0, lmax
            pl_raw(l + 1, spin) = pl_atoms(plusbulk)%potential%pl(l, spin)
         end do
      end do

      nshift = fixed_u_m0_a11_nshift
      allocate(pl_a11(lmax + 1, 2, nshift))
      do iv = 1, nshift
         pl_atoms = this%symbolic_atom
         call load_captured_pl_inputs(pl_atoms)
         do spin = 1, 2
            pl_atoms(plusbulk)%potential%gravity_center(3, spin) = &
               this%fixed_u_m0_pl_diag(ia, 3, spin, 1) - fixed_u_m0_a11_v(iv)
         end do
         call pl_bands%calculate_pl()
         do spin = 1, 2
            do l = 0, lmax
               pl_a11(l + 1, spin, iv) = pl_atoms(plusbulk)%potential%pl(l, spin)
            end do
         end do
      end do
      pl_reproduction_error = maxval(abs(pl_raw - this%fixed_u_m0_pl_diag(ia, 1:lmax + 1, 1:2, 6)))

      ! Reset the scratch atom to its independent copy of the accepted live
      ! state, then restore the original calculate_pl inputs.  Only d-shell
      ! gravity_center differs in this second call.
      pl_atoms = this%symbolic_atom
      call load_captured_pl_inputs(pl_atoms)
      do spin = 1, 2
         pl_atoms(plusbulk)%potential%gravity_center(3, spin) = &
            this%fixed_u_m0_local_center(ia, spin) - pl_atoms(plusbulk)%potential%vmad
      end do
      call pl_bands%calculate_pl()
      do spin = 1, 2
         do l = 0, lmax
            pl_local(l + 1, spin) = pl_atoms(plusbulk)%potential%pl(l, spin)
         end do
      end do

      allocate(ql_raw(3, lmax + 1, 2), ql_local(3, lmax + 1, 2), &
               ql_vout(3, lmax + 1, 2), ql_pt(3, lmax + 1, 2))
      allocate(baseline_ev(max(nval, 1)))
      baseline_ev = 0.0_rp
      have_raw_energy = .false.
      baseline_sumev = 0.0_rp
      baseline_sumec = 0.0_rp
      ql_raw = this%fixed_u_m0_ql(ia, :, 1:lmax + 1, :)
      ql_local = ql_raw
      ql_vout = ql_raw
      ql_pt = ql_raw
      do spin = 1, 2
         ql_local(3, 3, spin) = this%fixed_u_m0_local_q2(ia, spin)
         ql_vout(3, 3, spin) = this%fixed_u_m0_vout_q2(ia, spin)
         ql_pt(3, 3, spin) = this%fixed_u_m0_pt_q2(ia, spin)
      end do
      ql_vs_live_error = maxval(abs(ql_raw - ql_actual))

      open(newunit=unit, file=trim(this%fixed_u_trace_filename), status='old', position='append', &
           action='write', iostat=ios)
      if (ios == 0) then
         call this%write_fixed_u_csv_row(unit, 'm0_replay_match', this%fixed_u_trace_iteration, 'FROZEN_INPUT', &
            ia - 1, this%lattice%ib(ia) - 1, 2, -1, -1, -1, -1, -1, 0, 0.0_rp, &
            pl_reproduction_error, ql_vs_live_error, 0.0_rp, 0.0_rp)
         q2_min = minval(this%fixed_u_m0_local_q2(ia, :))
         call this%write_fixed_u_csv_row(unit, 'm0_replay_q2_min', this%fixed_u_trace_iteration, 'LOCAL_PL', &
            ia - 1, this%lattice%ib(ia) - 1, 2, -1, -1, -1, -1, -1, 0, 0.0_rp, q2_min, 0.0_rp, 0.0_rp, 0.0_rp)
         close(unit)
      end if

      ! Replay the exact PL/QL supplied to this live NEWRHO call before using
      ! the unmixed spectral q0/q2 values below.  This is the same-state
      ! control for the copied radial solver; RAW_PL remains a separate
      ! diagnostic because its moments are from the current spectral output.
      call run_replay_case('LIVE_STATE', pl_actual, ql_actual)
      call run_replay_case('RAW_PL', pl_raw, ql_raw)
      do iv = 1, nshift
         call run_replay_case(a11_case_name(iv), pl_a11(:, :, iv), ql_raw, fixed_u_m0_a11_v(iv), iv)
      end do
      if (minval(this%fixed_u_m0_local_q2(ia, :)) >= 0.0_rp) then
         call run_replay_case('LOCAL_PL', pl_local, ql_local)
      else
         open(newunit=unit, file=trim(this%fixed_u_trace_filename), status='old', position='append', &
              action='write', iostat=ios)
         if (ios == 0) then
            call this%write_fixed_u_csv_row(unit, 'm0_replay_skipped', this%fixed_u_trace_iteration, 'LOCAL_PL', &
               ia - 1, this%lattice%ib(ia) - 1, 2, -1, -1, -1, -1, -1, 0, 0.0_rp, &
               minval(this%fixed_u_m0_local_q2(ia, :)), 0.0_rp, 0.0_rp, 0.0_rp)
            close(unit)
         end if
      end if
      if (this%fixed_u_m0_p_available .and. minval(this%fixed_u_m0_vout_q2(ia, :)) >= 0.0_rp) &
         call run_replay_case('VOUT_CONTROL', pl_raw, ql_vout)
      if (this%fixed_u_m0_p_available .and. minval(this%fixed_u_m0_pt_q2(ia, :)) >= 0.0_rp) &
         call run_replay_case('P_TRANSPOSE_CONTROL', pl_raw, ql_pt)
      nullify(pl_bands%symbolic_atom, pl_bands%lattice)
      deallocate(pl_atoms, pl_raw, pl_local, pl_a11, ql_raw, ql_local, ql_vout, ql_pt, baseline_ev)
   contains
      function a11_case_name(shift_index) result(name)
         integer, intent(in) :: shift_index
         character(len=24) :: name
         write(name, '("A1_SHIFT_",i2.2)') shift_index
      end function a11_case_name

      function a11_stage_name(case_name) result(name)
         character(len=*), intent(in) :: case_name
         character(len=24) :: name
         name = case_name
      end function a11_stage_name

      logical function is_a11_shift_case(case_name)
         character(len=*), intent(in) :: case_name
         is_a11_shift_case = index(trim(case_name), 'A1_SHIFT_') == 1
      end function is_a11_shift_case

      subroutine load_captured_pl_inputs(atoms)
         type(symbolic_atom), intent(inout) :: atoms(:)
         integer :: site, pb, shell, channel
         do site = 1, this%lattice%nrec
            pb = this%lattice%nbulk + site
            do channel = 1, 2
               do shell = 0, lmax
                  atoms(pb)%potential%gravity_center(shell + 1, channel) = &
                     this%fixed_u_m0_pl_diag(site, shell + 1, channel, 1)
                  atoms(pb)%potential%c(shell, channel) = this%fixed_u_m0_pl_diag(site, shell + 1, channel, 2)
                  atoms(pb)%potential%qpar(shell, channel) = this%fixed_u_m0_pl_diag(site, shell + 1, channel, 3)
                  atoms(pb)%potential%srdel(shell, channel) = this%fixed_u_m0_pl_diag(site, shell + 1, channel, 4)
                  atoms(pb)%potential%pl(shell, channel) = this%fixed_u_m0_pl_diag(site, shell + 1, channel, 5)
               end do
            end do
         end do
      end subroutine load_captured_pl_inputs

      subroutine run_replay_case(case_name, pl_case, ql_case, a11_v, a11_index)
         character(len=*), intent(in) :: case_name
         real(rp), intent(in) :: pl_case(lmax + 1, 2), ql_case(3, lmax + 1, 2)
         real(rp), intent(in), optional :: a11_v
         integer, intent(in), optional :: a11_index
         type(self) :: replay_state
         type(symbolic_atom) :: replay_atom
      type(xc) :: replay_xc
      integer :: out_unit, out_ios, ir, channel, shell, moment, state_index, shift_index
      real(rp) :: center_input, shift_value, ev_raw_spin(2), ev_case_spin(2)
      real(rp) :: ev_raw_total, ev_case_total
      real(rp) :: charge(2), drdi, wgt
      logical :: a11_case
         shift_value = 0.0_rp
         shift_index = 0
         if (present(a11_v)) shift_value = a11_v
         if (present(a11_index)) shift_index = a11_index
         a11_case = trim(case_name) == 'RAW_PL' .or. is_a11_shift_case(case_name)
         allocate(rho_case(nr, 2), sumec_case(2), sumev_case(2), &
                  ec_case(max(ncore, 1)), ev_case(max(nval, 1)))
         allocate(replay_state%vtn(nr, 2), replay_state%vzt(nr, 2), &
                  replay_state%fun2(nr, lmax + 1, 2), replay_state%phi_amp(nr, lmax + 1, 2))
         replay_state%control => this%control
         replay_state%magnetic_seed_active = this%magnetic_seed_active
         replay_state%magnetic_seed_field = this%magnetic_seed_field
         replay_state%vtn = 0.0_rp
         replay_state%vzt = 0.0_rp
         replay_state%fun2 = 0.0_rp
         replay_state%phi_amp = 0.0_rp
         rho_case = 0.0_rp
         sumec_case = 0.0_rp
         sumev_case = 0.0_rp
         if (ncore > 0) ec_case(1:ncore) = ec_seed(1:ncore)
         if (nval > 0) ev_case(1:nval) = ev_seed(1:nval)
         replay_atom = atom
         call replay_state%newrho(replay_atom, z, lmax, a, b, nr, rofi, v, rho_case, pl_case, ql_case, &
                                  sumec_case, sumev_case, ec_case, ev_case, tol, 2, 0)
         if (trim(case_name) == 'RAW_PL') then
            baseline_sumev = sumev_case
            baseline_sumec = sumec_case
            if (nval > 0) baseline_ev(1:nval) = ev_case(1:nval)
            have_raw_energy = .true.
         end if
         replay_xc = xc(this%lattice%control)
         charge = 0.0_rp
         do channel = 1, 2
            do ir = 1, nr
               wgt = 2.0_rp*(mod(ir + 1, 2) + 1)/3.0_rp
               if (ir == 1 .or. ir == nr) wgt = 1.0_rp/3.0_rp
               drdi = a*(rofi_base(ir) + b)
               charge(channel) = charge(channel) + wgt*drdi*rho_case(ir, channel)
            end do
         end do
         open(newunit=out_unit, file=trim(this%fixed_u_trace_filename), status='old', position='append', &
              action='write', iostat=out_ios)
         if (out_ios == 0) then
            ! Recompute both component sets on copied arrays and write them
            ! only after the trace stream is open. ATOMSC evaluates the
            ! actual electrostatic/XC terms on rho_in before NEWRHO; the
            ! output-rho set is a separately labelled diagnostic functional.
            call evaluate_frozen_components(this, replay_state, replay_xc, atom, z, a, b, nr, rofi_base, &
               v_base, vhrmax, b_fsm, rho_in, sumec_case, sumev_case, out_unit, ia, &
               'ATOMSC_INPUT', component_snapshot)
            ! Use the VNUCL saved by the actual ATOMSC call-site, not an
            ! inferred boundary value, for the direct production comparison.
            actual_rhovh = sum(rvh_actual)
            actual_rhomu = sum(rmu_actual)
            actual_rhoeps = sum(reps_actual)
            actual_zvnucl = -z*vnucl_actual
            actual_ekin = sum(sumev_actual)+sum(sumec_actual)-actual_rhovh-actual_rhomu
            actual_utot = 0.5_rp*(actual_rhovh+actual_zvnucl)
            actual_etot = actual_ekin+actual_utot+actual_rhoeps
            actual_component = [sum(sumev_actual), sum(sumec_actual), actual_rhovh, actual_rhomu, &
                                vnucl_actual, actual_zvnucl, actual_ekin, actual_utot, actual_rhoeps, actual_etot]
            call evaluate_frozen_components(this, replay_state, replay_xc, atom, z, a, b, nr, rofi_base, &
               v_base, vhrmax, b_fsm, rho_case, sumec_case, sumev_case, out_unit, ia, &
               'OUTPUT_RHO', output_component)
            ! Emit the exact ATOMSC site-energy reconstruction residual. This
            ! is a same-state replay check, not an energy correction.
            do shell = 0, 9
               call this%write_fixed_u_csv_row(out_unit, 'm0_replay_component_delta', &
                  this%fixed_u_trace_iteration, trim(case_name), ia-1, this%lattice%ib(ia)-1, 2, -1, &
                  -1, -1, -1, -1, shell, 0.0_rp, 0.0_rp, component_snapshot(shell+1)-actual_component(shell+1), &
                  component_snapshot(shell+1), actual_component(shell+1))
            end do
            call this%write_fixed_u_csv_row(out_unit, 'm0_replay_control', this%fixed_u_trace_iteration, &
               trim(case_name), ia - 1, this%lattice%ib(ia) - 1, 2, -1, -1, -1, -1, -1, 0, &
               0.0_rp, a, z, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(out_unit, 'm0_replay_control', this%fixed_u_trace_iteration, &
               trim(case_name), ia - 1, this%lattice%ib(ia) - 1, 2, -1, -1, -1, -1, -1, 1, &
               0.0_rp, b, tol, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(out_unit, 'm0_replay_control', this%fixed_u_trace_iteration, &
               trim(case_name), ia - 1, this%lattice%ib(ia) - 1, 2, -1, -1, -1, -1, -1, 2, &
               0.0_rp, real(nr, rp), real(ncore, rp), real(nval, rp), 0.0_rp)
            call this%write_fixed_u_csv_row(out_unit, 'm0_replay_control', this%fixed_u_trace_iteration, &
               trim(case_name), ia - 1, this%lattice%ib(ia) - 1, 2, -1, -1, -1, -1, -1, 3, &
               0.0_rp, shift_value, this%fixed_u_m0_trvp(ia, 1), &
               this%fixed_u_m0_trvp(ia, 2), 0.0_rp)
            call this%write_fixed_u_csv_row(out_unit, 'm0_replay_control', this%fixed_u_trace_iteration, &
               trim(case_name), ia - 1, this%lattice%ib(ia) - 1, 2, -1, -1, -1, -1, -1, 4, &
               0.0_rp, vhrmax, 0.0_rp, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(out_unit, 'm0_replay_control', this%fixed_u_trace_iteration, &
               trim(case_name), ia - 1, this%lattice%ib(ia) - 1, 2, -1, -1, -1, -1, -1, 5, &
               0.0_rp, b_fsm, 0.0_rp, 0.0_rp, 0.0_rp)
            ! Serialize the exact PL tuple passed to the production scratch
            ! calculate_pl call and its returned PL for every shell/spin.
            do channel = 1, 2
               do shell = 0, lmax
                  call this%write_fixed_u_csv_row(out_unit, 'm0_replay_radial_parameter', &
                     this%fixed_u_trace_iteration, trim(case_name), ia - 1, this%lattice%ib(ia) - 1, &
                     shell, channel, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, &
                     atom%potential%c(shell, channel), 0.0_rp, 0.0_rp)
                  call this%write_fixed_u_csv_row(out_unit, 'm0_replay_radial_parameter', &
                     this%fixed_u_trace_iteration, trim(case_name), ia - 1, this%lattice%ib(ia) - 1, &
                     shell, channel, -1, -1, -1, -1, 1, 0.0_rp, 0.0_rp, &
                     atom%potential%qpar(shell, channel), 0.0_rp, 0.0_rp)
                  call this%write_fixed_u_csv_row(out_unit, 'm0_replay_radial_parameter', &
                     this%fixed_u_trace_iteration, trim(case_name), ia - 1, this%lattice%ib(ia) - 1, &
                     shell, channel, -1, -1, -1, -1, 2, 0.0_rp, 0.0_rp, &
                     atom%potential%srdel(shell, channel), 0.0_rp, 0.0_rp)
                  center_input = this%fixed_u_m0_pl_diag(ia, shell + 1, channel, 1)
                  if (shell == 2) then
                     select case (trim(case_name))
                     case ('LOCAL_PL')
                        center_input = this%fixed_u_m0_local_center(ia, channel) - atom%potential%vmad
                     case ('VOUT_CONTROL')
                        center_input = this%fixed_u_m0_vout_center(ia, channel) - atom%potential%vmad
                     case ('P_TRANSPOSE_CONTROL')
                        center_input = this%fixed_u_m0_pt_center(ia, channel) - atom%potential%vmad
                     end select
                     if (is_a11_shift_case(case_name)) &
                        center_input = this%fixed_u_m0_pl_diag(ia, shell + 1, channel, 1) - shift_value
                  end if
                  do moment = 0, 4
                     select case (moment)
                     case (0)
                        center_input = this%fixed_u_m0_pl_diag(ia, shell + 1, channel, 1)
                        if (shell == 2) then
                           select case (trim(case_name))
                           case ('LOCAL_PL')
                              center_input = this%fixed_u_m0_local_center(ia, channel) - atom%potential%vmad
                           case ('VOUT_CONTROL')
                              center_input = this%fixed_u_m0_vout_center(ia, channel) - atom%potential%vmad
                           case ('P_TRANSPOSE_CONTROL')
                              center_input = this%fixed_u_m0_pt_center(ia, channel) - atom%potential%vmad
                           end select
                           if (is_a11_shift_case(case_name)) &
                              center_input = this%fixed_u_m0_pl_diag(ia, shell + 1, channel, 1) - shift_value
                        end if
                     case (1)
                        center_input = this%fixed_u_m0_pl_diag(ia, shell + 1, channel, 2)
                     case (2)
                        center_input = this%fixed_u_m0_pl_diag(ia, shell + 1, channel, 3)
                     case (3)
                        center_input = this%fixed_u_m0_pl_diag(ia, shell + 1, channel, 4)
                     case (4)
                        center_input = this%fixed_u_m0_pl_diag(ia, shell + 1, channel, 5)
                     end select
                     call this%write_fixed_u_csv_row(out_unit, 'm0_replay_pl_input', this%fixed_u_trace_iteration, &
                        trim(case_name), ia - 1, this%lattice%ib(ia) - 1, shell, channel, -1, -1, -1, -1, moment, &
                        0.0_rp, 0.0_rp, center_input, 0.0_rp, 0.0_rp)
                  end do
                  call this%write_fixed_u_csv_row(out_unit, 'm0_replay_pl_return', this%fixed_u_trace_iteration, &
                     trim(case_name), ia - 1, this%lattice%ib(ia) - 1, shell, channel, -1, -1, -1, -1, 0, &
                     0.0_rp, 0.0_rp, pl_case(shell + 1, channel), 0.0_rp, 0.0_rp)
                  do moment = 1, 3
                     call this%write_fixed_u_csv_row(out_unit, 'm0_replay_ql_input', this%fixed_u_trace_iteration, &
                        trim(case_name), ia - 1, this%lattice%ib(ia) - 1, shell, channel, -1, -1, -1, -1, moment - 1, &
                        0.0_rp, 0.0_rp, ql_case(moment, shell + 1, channel), 0.0_rp, 0.0_rp)
                  end do
               end do
            end do
            ! Each case records the same immutable radial solver inputs. The
            ! independent checker compares these rows exactly across cases.
            do ir = 1, nr
               call this%write_fixed_u_csv_row(out_unit, 'm0_replay_rofi', this%fixed_u_trace_iteration, &
                  trim(case_name), ia - 1, this%lattice%ib(ia) - 1, -1, -1, -1, -1, ir - 1, -1, 0, &
                  0.0_rp, 0.0_rp, rofi(ir), 0.0_rp, 0.0_rp)
               do channel = 1, 2
                  call this%write_fixed_u_csv_row(out_unit, 'm0_replay_v_input', this%fixed_u_trace_iteration, &
                     trim(case_name), ia - 1, this%lattice%ib(ia) - 1, -1, channel, -1, -1, ir - 1, -1, 0, &
                     0.0_rp, 0.0_rp, v(ir, channel), 0.0_rp, 0.0_rp)
                  call this%write_fixed_u_csv_row(out_unit, 'm0_replay_rho_input', this%fixed_u_trace_iteration, &
                     trim(case_name), ia - 1, this%lattice%ib(ia) - 1, -1, channel, -1, -1, ir - 1, -1, 0, &
                     0.0_rp, 0.0_rp, rho_in(ir, channel), 0.0_rp, 0.0_rp)
               end do
            end do
            do row = 1, ncore
               call this%write_fixed_u_csv_row(out_unit, 'm0_replay_ec_seed', this%fixed_u_trace_iteration, &
                  trim(case_name), ia - 1, this%lattice%ib(ia) - 1, -1, -1, -1, -1, row - 1, -1, 0, &
                  0.0_rp, 0.0_rp, ec_seed(row), 0.0_rp, 0.0_rp)
            end do
            do row = 1, nval
               call this%write_fixed_u_csv_row(out_unit, 'm0_replay_ev_seed', this%fixed_u_trace_iteration, &
                  trim(case_name), ia - 1, this%lattice%ib(ia) - 1, -1, -1, -1, -1, row - 1, -1, 0, &
                  0.0_rp, 0.0_rp, ev_seed(row), 0.0_rp, 0.0_rp)
            end do
            do channel = 1, 2
               call this%write_fixed_u_csv_row(out_unit, 'm0_replay_energy', this%fixed_u_trace_iteration, &
                  trim(case_name), ia - 1, this%lattice%ib(ia) - 1, 2, channel, -1, -1, -1, -1, 0, &
                  0.0_rp, sumec_case(channel), sumev_case(channel), 0.0_rp, 0.0_rp)
               if (have_raw_energy) then
                  call this%write_fixed_u_csv_row(out_unit, 'm0_replay_delta', this%fixed_u_trace_iteration, &
                     trim(case_name), ia - 1, this%lattice%ib(ia) - 1, 2, channel, -1, -1, -1, -1, 0, &
                     0.0_rp, sumec_case(channel)-baseline_sumec(channel), &
                     sumev_case(channel)-baseline_sumev(channel), 0.0_rp, 0.0_rp)
               end if
            end do
            if (a11_case) then
               ev_raw_spin = 0.0_rp
               ev_case_spin = 0.0_rp
               ev_raw_total = 0.0_rp
               ev_case_total = 0.0_rp
               state_index = 0
               do channel = 1, 2
                  do shell = 0, lmax
                     if (ql_case(1, shell + 1, channel) < 1.0e-5_rp) cycle
                     state_index = state_index + 1
                     if (state_index > nval) error stop 'A1 trace: NEWRHO EV ordering exceeds nval'
                     ev_raw_spin(channel) = ev_raw_spin(channel) + baseline_ev(state_index)
                     ev_case_spin(channel) = ev_case_spin(channel) + ev_case(state_index)
                     call this%write_fixed_u_csv_row(out_unit, 'm0_a11_ev', this%fixed_u_trace_iteration, &
                        a11_stage_name(case_name), ia-1, this%lattice%ib(ia)-1, shell, channel, -1, &
                        state_index-1, -1, -1, shift_index, shift_value, baseline_ev(state_index), &
                        ev_case(state_index), ev_case(state_index)-baseline_ev(state_index), 0.0_rp)
                  end do
                  call this%write_fixed_u_csv_row(out_unit, 'm0_a11_ev_sum', this%fixed_u_trace_iteration, &
                     a11_stage_name(case_name), ia-1, this%lattice%ib(ia)-1, 2, channel, -1, -1, -1, -1, &
                     shift_index, shift_value, ev_raw_spin(channel), ev_case_spin(channel), &
                     ev_case_spin(channel)-ev_raw_spin(channel), 0.0_rp)
                  ev_raw_total = ev_raw_total + ev_raw_spin(channel)
                  ev_case_total = ev_case_total + ev_case_spin(channel)
               end do
               call this%write_fixed_u_csv_row(out_unit, 'm0_a11_ev_sum', this%fixed_u_trace_iteration, &
                  a11_stage_name(case_name), ia-1, this%lattice%ib(ia)-1, 2, -1, -1, -1, -1, -1, &
                  shift_index, shift_value, ev_raw_total, ev_case_total, ev_case_total-ev_raw_total, 0.0_rp)
               do channel = 1, 2
                  call this%write_fixed_u_csv_row(out_unit, 'm0_a11_case', this%fixed_u_trace_iteration, &
                     a11_stage_name(case_name), ia-1, this%lattice%ib(ia)-1, 2, channel, -1, -1, -1, -1, &
                     shift_index, shift_value, baseline_sumev(channel), sumev_case(channel), &
                     this%fixed_u_m0_q0(ia, channel), this%fixed_u_m0_ql(ia, 3, 3, channel))
                  call this%write_fixed_u_csv_row(out_unit, 'm0_a11_center', this%fixed_u_trace_iteration, &
                     a11_stage_name(case_name), ia-1, this%lattice%ib(ia)-1, 2, channel, -1, -1, -1, -1, &
                     shift_index, shift_value, this%fixed_u_m0_pl_diag(ia, 3, channel, 1), &
                     this%fixed_u_m0_pl_diag(ia, 3, channel, 1)-shift_value, &
                     pl_raw(3, channel), pl_case(3, channel))
               end do
               call this%write_fixed_u_csv_row(out_unit, 'm0_a11_case', this%fixed_u_trace_iteration, &
                  a11_stage_name(case_name), ia-1, this%lattice%ib(ia)-1, 2, -1, -1, -1, -1, -1, &
                  shift_index, shift_value, sum(baseline_sumev), sum(sumev_case), &
                  sum(this%fixed_u_m0_q0(ia, :)), sum(this%fixed_u_m0_ql(ia, 3, 3, :)))
               call this%write_fixed_u_csv_row(out_unit, 'm0_a11_complete', this%fixed_u_trace_iteration, &
                  a11_stage_name(case_name), ia-1, this%lattice%ib(ia)-1, 2, -1, -1, -1, -1, -1, &
                  shift_index, shift_value, real(nval, rp), 1.0_rp, 0.0_rp, 0.0_rp)
            end if
            ! NEWRHO returns rho/SUMEC/SUMEV, but the full ATOMSC ETOT
            ! components are recomputed independently below from copied arrays.
            call this%write_fixed_u_csv_row(out_unit, 'm0_replay_radial_charge', &
               this%fixed_u_trace_iteration, trim(case_name), ia - 1, this%lattice%ib(ia)-1, &
               -1, 1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, charge(1), 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(out_unit, 'm0_replay_radial_charge', &
               this%fixed_u_trace_iteration, trim(case_name), ia - 1, this%lattice%ib(ia)-1, &
               -1, 2, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, charge(2), 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(out_unit, 'm0_replay_radial_charge', &
               this%fixed_u_trace_iteration, trim(case_name), ia - 1, this%lattice%ib(ia)-1, &
               -1, 0, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, charge(1)+charge(2), &
               charge(1)-charge(2), 0.0_rp)
            do ir = 1, nr
               do channel = 1, 2
                  call this%write_fixed_u_csv_row(out_unit, 'm0_replay_rho', this%fixed_u_trace_iteration, &
                     trim(case_name), ia - 1, this%lattice%ib(ia) - 1, -1, channel, -1, -1, ir - 1, -1, 0, &
                     0.0_rp, 0.0_rp, rho_case(ir, channel), 0.0_rp, 0.0_rp)
               end do
            end do
            close(out_unit)
         end if
         deallocate(replay_state%vtn, replay_state%vzt, replay_state%fun2, replay_state%phi_amp)
         deallocate(rho_case, sumec_case, sumev_case, ec_case, ev_case)
      end subroutine run_replay_case

      subroutine evaluate_frozen_components(this, replay_state, replay_xc, atom, z, a, b, nr, rofi_base, &
                                            v_base, vhrmax, b_fsm, rho_sample, sumec_case, sumev_case, &
                                            out_unit, site, component_stage, components)
         class(self), intent(inout) :: this, replay_state
         type(xc), intent(in) :: replay_xc
         class(symbolic_atom), intent(in) :: atom
         real(rp), intent(in) :: z, a, b, rofi_base(nr), v_base(nr,2), vhrmax, b_fsm
         integer, intent(in) :: nr, out_unit, site
         character(len=*), intent(in) :: component_stage
         real(rp), intent(in) :: rho_sample(nr,2), sumec_case(2), sumev_case(2)
         real(rp), intent(out) :: components(10)
         integer :: isp, comp
         real(rp) :: rhoeps_local, rhomu_local, rhovh_local, utot_local, ekin_local, etot_local
         real(rp) :: vnucl_local, zvnucl_local, vsum_local, value_local
         real(rp), allocatable :: rho_work(:,:), v_work(:,:), rofi_work(:)
         real(rp) :: rvh_local(2), rho0_local(2), reps_local(2), rmu_local(2)

         allocate(rho_work(nr,2), v_work(nr,2), rofi_work(nr))
         rho_work = rho_sample
         v_work = v_base
         rofi_work = rofi_base
         rvh_local = 0.0_rp
         rho0_local = 0.0_rp
         reps_local = 0.0_rp
         rmu_local = 0.0_rp
         vsum_local = 0.0_rp
         call POISS0(z, a, b, rofi_work, rho_work, nr, vhrmax, v_work, rvh_local, vsum_local, 2)
         vnucl_local = v_work(1,1)
         zvnucl_local = -z*vnucl_local
         call replay_state%VXC0SP(replay_xc, z, a, b, rofi_work, rho_sample, nr, v_work, &
                                  rho0_local, reps_local, rmu_local, 2, b_fsm)
         rhovh_local = sum(rvh_local)
         rhoeps_local = sum(reps_local)
         rhomu_local = sum(rmu_local)
         utot_local = 0.5_rp*(rhovh_local + zvnucl_local)
         ekin_local = sum(sumev_case) + sum(sumec_case) - rhovh_local - rhomu_local
         etot_local = ekin_local + utot_local + rhoeps_local
         components = [sum(sumev_case), sum(sumec_case), rhovh_local, rhomu_local, vnucl_local, &
                       zvnucl_local, ekin_local, utot_local, rhoeps_local, etot_local]
         do comp = 0, 9
            value_local = components(comp+1)
            call this%write_fixed_u_csv_row(out_unit, 'm0_replay_energy_component', &
               this%fixed_u_trace_iteration, trim(component_stage), site-1, this%lattice%ib(site)-1, &
               2, -1, -1, -1, -1, -1, comp, 0.0_rp, 0.0_rp, value_local, 0.0_rp, 0.0_rp)
         end do
         do isp = 1, 2
            call this%write_fixed_u_csv_row(out_unit, 'm0_replay_energy_component', &
               this%fixed_u_trace_iteration, trim(component_stage), site-1, this%lattice%ib(site)-1, &
               2, isp, -1, -1, -1, -1, 10, 0.0_rp, 0.0_rp, rvh_local(isp), 0.0_rp, 0.0_rp)
         end do
         deallocate(rho_work, v_work, rofi_work)
      end subroutine evaluate_frozen_components
   end procedure replay_fixed_u_m0_moments


   module procedure update_fixed_u_residuals
      complex(rp), allocatable :: density_override(:, :, :, :, :), vout(:, :, :)
      complex(rp), allocatable :: p_in(:, :, :, :), p_out(:, :, :, :), v_in(:, :, :, :), v_out(:, :, :, :)
      integer, allocatable :: site_for_type(:), rec_type(:)
      integer :: nsite, ntype, block, isite, iclu, itype, il, ispin, mmax, i, j, flat, orbital_i, orbital_j, worst_block
      integer :: worst_l, worst_site, unit, ios
      real(rp) :: rp_abs, rp_rel, rv_abs, rv_rel, antiherm_p, antiherm_v, u_value, j_value
      real(rp) :: rp_target, rv_target
      logical :: active, physical

      this%fixed_u_residual_active = .false.
      this%fixed_u_residual_converged = .false.
      this%fixed_u_residual_physical = .false.
      if (.not. this%hamiltonian%hubbard_u_general_check .or. &
          trim(this%hamiltonian%hubbard_u_potential_form) /= 'liechtenstein' .or. &
          this%hamiltonian%hubbard_u_sc_check .or. this%hamiltonian%hubbard_v_check) then
         this%mix%fixed_u_stable_updates = 0
         return
      end if

      nsite = this%lattice%nrec
      ntype = this%lattice%ntype
      if (nsite < 1 .or. ntype < 1) return
      if (.not. allocated(this%lattice%irec)) error stop 'active fixed-U residual: recursion-site map is not allocated'
      if (.not. allocated(this%lattice%iz)) error stop 'active fixed-U residual: cluster type map is not allocated'
      if (.not. allocated(this%lattice%symbolic_atoms)) &
         error stop 'active fixed-U residual: symbolic atom table is not allocated'
      if (size(this%lattice%irec) < nsite) error stop 'active fixed-U residual: recursion-site map is too short'
      if (size(this%lattice%symbolic_atoms) < ntype) &
         error stop 'active fixed-U residual: symbolic atom table is too short'
      allocate(site_for_type(ntype), rec_type(nsite))
      site_for_type = 0
      do isite = 1, nsite
         iclu = this%lattice%irec(isite)
         if (iclu < 1) error stop 'active fixed-U residual: recursion site has no cluster mapping'
         if (iclu > size(this%lattice%iz)) error stop 'active fixed-U residual: cluster mapping is out of range'
         itype = this%lattice%iz(iclu)
         if (itype < 1 .or. itype > ntype) error stop 'active fixed-U residual: mapped atom type is out of range'
         rec_type(isite) = itype
         if (site_for_type(itype) == 0) site_for_type(itype) = isite
      end do

      active = this%hamiltonian%hubbard_v_check .or. this%hamiltonian%hubbard_u_sc_check .or. &
               trim(this%hamiltonian%hubbard_u_potential_form) == 'acbn0'
      do itype = 1, ntype
         do il = 1, min(4, size(this%mix%ldm_new, 2))
            u_value = 0.0_rp
            j_value = 0.0_rp
            if (allocated(this%lattice%symbolic_atoms(itype)%potential%hubbard_u)) then
               if (il <= size(this%lattice%symbolic_atoms(itype)%potential%hubbard_u)) &
                  u_value = this%lattice%symbolic_atoms(itype)%potential%hubbard_u(il)
            end if
            if (allocated(this%lattice%symbolic_atoms(itype)%potential%hubbard_j)) then
               if (il <= size(this%lattice%symbolic_atoms(itype)%potential%hubbard_j)) &
                  j_value = this%lattice%symbolic_atoms(itype)%potential%hubbard_j(il)
            end if
            if (abs(u_value) > 1.0e-10_rp .or. abs(j_value) > 1.0e-10_rp) active = .true.
         end do
      end do
      if (.not. active) then
         this%mix%fixed_u_stable_updates = 0
         return
      end if
      this%fixed_u_residual_active = .true.
      if (.not. allocated(this%hamiltonian%hubbard_u_pot)) &
         error stop 'active fixed-U residual: Hubbard potential is not allocated'

      allocate(density_override(ntype, 4, 2, 7, 7), vout(nb, nb, ntype))
      allocate(p_in(nsite*4, 2, 7, 7), p_out(nsite*4, 2, 7, 7), &
               v_in(nsite*4, 2, 7, 7), v_out(nsite*4, 2, 7, 7))
      density_override = 0.0_rp
      vout = 0.0_rp
      p_in = 0.0_rp
      p_out = 0.0_rp
      v_in = 0.0_rp
      v_out = 0.0_rp

      ! The production Hubbard operator is type-local. Use the same first
      ! site-of-type density convention for its read-only Pout evaluation;
      ! residuals below remain site-resolved using each site's own P matrices.
      do itype = 1, ntype
         isite = site_for_type(itype)
         if (isite == 0) cycle
         ! Only traverse angular-momentum channels represented by this build's
         ! allocated LDM.  The fixed four-shell potential storage can exceed
         ! the active basis (e.g. a d-only build has ldm_new(:,1:3,...)).
         do il = 1, min(4, size(this%mix%ldm_new, 2))
            mmax = 2*il - 1
            do ispin = 1, 2
               do i = 1, mmax
                  do j = 1, mmax
                     flat = (i - 1)*mmax + j
                     density_override(itype, il, ispin, i, j) = this%mix%ldm_new(isite, il, ispin, flat)
                  end do
               end do
            end do
         end do
      end do
      call this%hamiltonian%calculate_hubbard_u_potential_general( &
         density_override=density_override, potential_out=vout)

      do isite = 1, nsite
         itype = rec_type(isite)
         do il = 1, min(4, size(this%mix%ldm_new, 2))
            u_value = 0.0_rp
            j_value = 0.0_rp
            if (allocated(this%lattice%symbolic_atoms(itype)%potential%hubbard_u)) then
               if (il <= size(this%lattice%symbolic_atoms(itype)%potential%hubbard_u)) &
                  u_value = this%lattice%symbolic_atoms(itype)%potential%hubbard_u(il)
            end if
            if (allocated(this%lattice%symbolic_atoms(itype)%potential%hubbard_j)) then
               if (il <= size(this%lattice%symbolic_atoms(itype)%potential%hubbard_j)) &
                  j_value = this%lattice%symbolic_atoms(itype)%potential%hubbard_j(il)
            end if
            if (abs(u_value) <= 1.0e-10_rp .and. abs(j_value) <= 1.0e-10_rp) cycle
            block = (il - 1)*nsite + isite
            mmax = 2*il - 1
            do ispin = 1, 2
               do i = 1, mmax
                  do j = 1, mmax
                     flat = (i - 1)*mmax + j
                     p_in(block, ispin, i, j) = this%mix%ldm_old(isite, il, ispin, flat)
                     p_out(block, ispin, i, j) = this%mix%ldm_new(isite, il, ispin, flat)
                     orbital_i = (il - 1)**2 + i + (ispin - 1)*spin_off
                     orbital_j = (il - 1)**2 + j + (ispin - 1)*spin_off
                     v_in(block, ispin, i, j) = this%hamiltonian%hubbard_u_pot(orbital_i, orbital_j, itype)
                     v_out(block, ispin, i, j) = vout(orbital_i, orbital_j, itype)
                  end do
               end do
            end do
         end do
      end do

      call fixed_u_residual_metrics(p_in, p_out, v_in, v_out, rp_abs, rp_rel, rv_abs, rv_rel, &
                                    physical, antiherm_p, antiherm_v, worst_block)
      if (worst_block > 0) then
         worst_site = mod(worst_block - 1, nsite) + 1
         worst_l = (worst_block - 1)/nsite
      else
         worst_site = 0
         worst_l = -1
      end if
      this%fixed_u_rp_abs = rp_abs
      this%fixed_u_rp_rel = rp_rel
      this%fixed_u_rv_abs = rv_abs
      this%fixed_u_rv_rel = rv_rel
      this%fixed_u_antiherm_p = antiherm_p
      this%fixed_u_antiherm_v = antiherm_v
      this%fixed_u_worst_site = worst_site
      this%fixed_u_worst_l = worst_l
      this%fixed_u_residual_physical = physical
      call fixed_u_r05_targets(this%use_kspace, this%r05_p_tolerance, this%r05_v_tolerance, &
         this%hamiltonian%fixed_u_test_tighten_r05, rp_target, rv_target)
      call this%mix%advance_fixed_u_convergence(rp_abs, rv_abs, physical, this%fixed_u_residual_converged, &
                                                rp_target, rv_target)
      if (rank == 0) then
         call g_logger%info('FIXED_U_RESIDUAL rp_abs='//fmt('es12.4',rp_abs)//' rp_rel='//fmt('es12.4',rp_rel)// &
            ' rv_abs_Ry='//fmt('es12.4',rv_abs)//' rv_rel='//fmt('es12.4',rv_rel)// &
            ' target_rp='//fmt('es12.4',rp_target)//' target_rv_Ry='//fmt('es12.4',rv_target)// &
            ' required_streak='//int2str(fixed_u_required_stable_updates)// &
            ' antiherm_P='//fmt('es12.4',antiherm_p)//' antiherm_V='//fmt('es12.4',antiherm_v)// &
            ' physical='//trim(log2str(physical))//' streak='//int2str(this%mix%fixed_u_stable_updates)// &
            ' worst_site='//int2str(worst_site)//' l='//int2str(worst_l), __FILE__, __LINE__)
      end if
      if (this%fixed_u_trace_enabled .and. rank == 0) then
         open(newunit=unit, file=trim(this%fixed_u_trace_filename), status='old', position='append', &
              action='write', iostat=ios)
         if (ios == 0) then
            call this%write_fixed_u_csv_row(unit, 'fixed_u_rp_abs', this%fixed_u_trace_iteration, 'SCF_CONVERGENCE', &
               worst_site - 1, -1, worst_l, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, rp_abs, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'fixed_u_rp_rel', this%fixed_u_trace_iteration, 'SCF_CONVERGENCE', &
               worst_site - 1, -1, worst_l, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, rp_rel, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'fixed_u_rv_abs_ry', this%fixed_u_trace_iteration, 'SCF_CONVERGENCE', &
               worst_site - 1, -1, worst_l, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, rv_abs, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'fixed_u_rv_rel', this%fixed_u_trace_iteration, 'SCF_CONVERGENCE', &
               worst_site - 1, -1, worst_l, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, rv_rel, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'fixed_u_rp_target', this%fixed_u_trace_iteration, 'SCF_CONVERGENCE', &
               worst_site - 1, -1, worst_l, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, rp_target, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'fixed_u_rv_target_ry', this%fixed_u_trace_iteration, 'SCF_CONVERGENCE', &
               worst_site - 1, -1, worst_l, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, rv_target, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'fixed_u_antiherm_p', this%fixed_u_trace_iteration, 'SCF_CONVERGENCE', &
               worst_site - 1, -1, worst_l, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, antiherm_p, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'fixed_u_antiherm_v', this%fixed_u_trace_iteration, 'SCF_CONVERGENCE', &
               worst_site - 1, -1, worst_l, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, antiherm_v, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'fixed_u_physical', this%fixed_u_trace_iteration, 'SCF_CONVERGENCE', &
               worst_site - 1, -1, worst_l, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, merge(1.0_rp,0.0_rp,physical), &
               0.0_rp, 0.0_rp)
            close(unit)
         end if
      end if
      deallocate(site_for_type, rec_type, density_override, vout, p_in, p_out, v_in, v_out)
   end procedure update_fixed_u_residuals


   module procedure write_fixed_u_total_energy
      integer :: unit, ios
      real(rp) :: value

      if (.not. this%fixed_u_trace_enabled .or. rank /= 0) return
      if (this%fixed_u_trace_selected_iteration > 0 .and. &
          this%fixed_u_trace_iteration /= this%fixed_u_trace_selected_iteration) return
      open(newunit=unit, file=trim(this%fixed_u_trace_filename), status='old', position='append', &
           action='write', iostat=ios)
      if (ios /= 0) return
      value = sum(this%symbolic_atom(:)%potential%sumev)
      call total_row('sum_ev', value)
      value = sum(this%symbolic_atom(:)%potential%sumec)
      call total_row('sum_ec', value)
      value = sum(this%symbolic_atom(:)%potential%ekin)
      call total_row('ekin', value)
      value = sum(this%symbolic_atom(:)%potential%utot)
      call total_row('utot', value)
      value = sum(this%symbolic_atom(:)%potential%rhoeps)
      call total_row('rhoeps', value)
      value = sum(this%symbolic_atom(:)%potential%etot)
      call total_row('etot_total', value)
      close(unit)
   contains
      subroutine total_row(label, energy_value)
         character(len=*), intent(in) :: label
         real(rp), intent(in) :: energy_value
         call this%write_fixed_u_csv_row(unit, label, this%fixed_u_trace_iteration, 'POST_SCF_RADIAL', &
            -1, -1, -1, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, energy_value, 0.0_rp, 0.0_rp)
      end subroutine total_row
   end procedure write_fixed_u_total_energy



   module procedure write_fixed_u_kspace_trace
      integer :: unit, ios, ik, ig, ib, i, j, n, nsite, nper_spin, base, m, spin, ell
      real(rp) :: kt, occ
      complex(rp) :: z

      if (.not. this%fixed_u_trace_enabled) return
      if (rank /= 0) return
      if (this%fixed_u_trace_ql_only) return
      if (this%fixed_u_trace_moments_only) return
      if (this%fixed_u_trace_selected_iteration > 0 .and. &
          this%fixed_u_trace_iteration /= this%fixed_u_trace_selected_iteration) then
         call g_logger%warning('Fixed-U k-space trace skipped: selected iteration mismatch.', __FILE__, __LINE__)
         return
      end if
      if (.not. allocated(this%reciprocal_scf_cache)) then
         call g_logger%warning('Fixed-U k-space trace skipped: reciprocal cache is not allocated.', __FILE__, __LINE__)
         return
      end if
      if (capture_spectrum .and. (.not. allocated(this%reciprocal_scf_cache%eigenvalues) .or. &
          .not. allocated(this%reciprocal_scf_cache%eigenvectors))) then
         call g_logger%warning('Fixed-U k-space spectrum trace skipped: eigenpairs are not allocated.', __FILE__, __LINE__)
         return
      end if
      open(newunit=unit, file=trim(this%fixed_u_trace_filename), status='old', position='append', &
           action='write', iostat=ios)
      if (ios /= 0) return
      kt = max(this%reciprocal_scf_cache%temperature*kB_Ry_per_K, 1.0e-10_rp)
      if (allocated(this%reciprocal_scf_cache%hk_bulk)) then
         n = size(this%reciprocal_scf_cache%hk_bulk, 1)
      else if (allocated(this%reciprocal_scf_cache%eigenvectors)) then
         n = size(this%reciprocal_scf_cache%eigenvectors, 1)
      else
         call g_logger%warning('Fixed-U k-space trace skipped: neither Hamiltonian nor eigenvectors are allocated.', &
                               __FILE__, __LINE__)
         close(unit)
         return
      end if
      nsite = this%lattice%nrec
      nper_spin = n/(2*nsite)
      if (capture_spectrum) then
         call this%write_fixed_u_csv_row(unit, 'fermi_level', this%fixed_u_trace_iteration, 'EIGENPAIRS', &
            -1, -1, -1, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, &
            this%reciprocal_scf_cache%fermi_level, 0.0_rp, 0.0_rp)
         call this%write_fixed_u_csv_row(unit, 'temperature_k', this%fixed_u_trace_iteration, 'EIGENPAIRS', &
            -1, -1, -1, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, &
            this%reciprocal_scf_cache%temperature, 0.0_rp, 0.0_rp)
         call this%write_fixed_u_csv_row(unit, 'temperature_ry', this%fixed_u_trace_iteration, 'EIGENPAIRS', &
            -1, -1, -1, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, kt, 0.0_rp, 0.0_rp)
         do ik = 1, size(this%reciprocal_scf_cache%eigenvalues, 2)
            ig = k_global_index(ik)
            do ib = 1, size(this%reciprocal_scf_cache%eigenvalues, 1)
               occ = fermi_dirac_occupation(this%reciprocal_scf_cache%eigenvalues(ib, ik), &
                                            this%reciprocal_scf_cache%fermi_level, kt)
               call this%write_fixed_u_csv_row(unit, 'eigenvalue', this%fixed_u_trace_iteration, 'EIGENPAIRS', &
                  -1, -1, -1, -1, ig - 1, ib - 1, -1, -1, 0, 0.0_rp, &
                  this%reciprocal_scf_cache%eigenvalues(ib, ik), 0.0_rp, 0.0_rp, 0.0_rp)
               call this%write_fixed_u_csv_row(unit, 'occupation', this%fixed_u_trace_iteration, 'EIGENPAIRS', &
                  -1, -1, -1, -1, ig - 1, ib - 1, -1, -1, 0, 0.0_rp, 0.0_rp, occ, 0.0_rp, 0.0_rp)
               do i = 1, size(this%reciprocal_scf_cache%eigenvectors, 1)
                  z = this%reciprocal_scf_cache%eigenvectors(i, ib, ik)
                  call this%write_fixed_u_csv_row(unit, 'eigenvector', this%fixed_u_trace_iteration, 'EIGENPAIRS', &
                     -1, -1, -1, -1, ig - 1, ib - 1, i - 1, -1, 0, 0.0_rp, 0.0_rp, 0.0_rp, &
                     real(z, rp), aimag(z))
               end do
            end do
         end do
      else
      do ik = 1, size(this%reciprocal_scf_cache%hk_bulk, 3)
         ig = ik
         if (allocated(this%reciprocal_scf_cache%k_l2g_map)) then
            if (ik <= size(this%reciprocal_scf_cache%k_l2g_map)) ig = this%reciprocal_scf_cache%k_l2g_map(ik)
         end if
         call write_k_header(ik)
         do i = 1, n
            do j = 1, n
               z = this%reciprocal_scf_cache%hk_bulk(i, j, ik)
               call this%write_fixed_u_csv_row(unit, 'h_total', this%fixed_u_trace_iteration, 'EIGENPAIRS', &
                  -1, -1, -1, -1, ig - 1, -1, i - 1, j - 1, 0, 0.0_rp, 0.0_rp, 0.0_rp, real(z, rp), aimag(z))
               if (allocated(this%reciprocal_scf_cache%sk_overlap)) then
                  z = this%reciprocal_scf_cache%sk_overlap(i, j, ik)
                  call this%write_fixed_u_csv_row(unit, 's_overlap', this%fixed_u_trace_iteration, 'EIGENPAIRS', &
                     -1, -1, -1, -1, ig - 1, -1, i - 1, j - 1, 0, 0.0_rp, 0.0_rp, 0.0_rp, &
                     real(z, rp), aimag(z))
               end if
            end do
         end do
      end do
      end if
      ! Basis ordering is site-major, then spin-major, then l,m; for the
      ! supported scalar-collinear spd layout, d orbitals are local indices 5:9.
      if (capture_spectrum) then
         do i = 1, nsite
            do spin = 1, 2
               base = (i - 1)*(2*nper_spin) + (spin - 1)*nper_spin
               do ell = 0, int(sqrt(real(nper_spin, rp))) - 1
                  do m = 0, 2*ell
                     call this%write_fixed_u_csv_row(unit, 'basis_map', this%fixed_u_trace_iteration, 'EIGENPAIRS', &
                        i - 1, this%lattice%iz(i), ell, spin - 1, -1, -1, base + ell*ell + m, m - ell, 0, &
                        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp)
                  end do
               end do
            end do
         end do
      end if
      close(unit)
   contains
      subroutine write_k_header(k_local)
         integer, intent(in) :: k_local
         integer :: kg
         kg = k_global_index(k_local)
         if (allocated(this%reciprocal_scf_cache%k_points) .and. allocated(this%reciprocal_scf_cache%k_weights)) then
            do i = 1, 3
               call krow('k_coordinate', kg, i - 1, -1, 0.0_rp, 0.0_rp, &
                         this%reciprocal_scf_cache%k_points(i, kg), 'fractional')
            end do
            call this%write_fixed_u_csv_row(unit, 'k_weight', this%fixed_u_trace_iteration, 'EIGENPAIRS', &
               -1, -1, -1, -1, kg - 1, -1, -1, -1, 0, &
               this%reciprocal_scf_cache%k_weights(kg), 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp)
         end if
      end subroutine write_k_header
      integer function k_global_index(k_local)
         integer, intent(in) :: k_local
         k_global_index = k_local
         if (allocated(this%reciprocal_scf_cache%k_l2g_map)) then
            if (k_local <= size(this%reciprocal_scf_cache%k_l2g_map)) &
               k_global_index = this%reciprocal_scf_cache%k_l2g_map(k_local)
         end if
      end function k_global_index
      subroutine krow(kind, k_index, row, col, re, im, value, name)
         integer, intent(in) :: k_index, row, col
         character(len=*), intent(in) :: kind, name
         real(rp), intent(in) :: re, im, value
         call this%write_fixed_u_csv_row(unit, trim(kind), this%fixed_u_trace_iteration, 'EIGENPAIRS', &
            -1, -1, -1, -1, k_index - 1, -1, row, col, 0, 0.0_rp, 0.0_rp, value, re, im)
      end subroutine krow
   end procedure write_fixed_u_kspace_trace


   module procedure write_fixed_u_rs_green
      integer :: ia, iloc, ie, i, j, npts, nmat, orbital, ell, m, spin
      real(rp), allocatable :: quadrature_weight(:)
      real(rp) :: h, tail
      complex(rp) :: z

      if (this%fixed_u_trace_moments_only) return
      if (.not. allocated(this%green%g0) .or. .not. allocated(this%en%ene)) return
      npts = this%bands%nv1
      if (npts < 3 .or. size(this%en%ene) < npts + 2) return
      nmat = size(this%green%g0, 1)
      allocate(quadrature_weight(size(this%en%ene)))
      quadrature_weight = 0.0_rp
      h = this%en%edel
      do i = 2, npts - 1, 2
         quadrature_weight(i - 1) = quadrature_weight(i - 1) + h/3.0_rp
         quadrature_weight(i) = quadrature_weight(i) + 4.0_rp*h/3.0_rp
         quadrature_weight(i + 1) = quadrature_weight(i + 1) + h/3.0_rp
      end do
      if (this%bands%e1 /= this%en%fermi) then
         tail = (this%en%fermi - this%bands%e1)/6.0_rp
         quadrature_weight(npts) = quadrature_weight(npts) + tail
         quadrature_weight(npts + 1) = quadrature_weight(npts + 1) + 4.0_rp*tail
         quadrature_weight(npts + 2) = quadrature_weight(npts + 2) + tail
      end if

      do ia = start_atom, end_atom
         iloc = g2l_map(ia)
         if (iloc < 1 .or. iloc > size(this%green%g0, 4)) cycle
         do i = 0, nmat - 1
            orbital = mod(i, spin_off)
            spin = i/spin_off
            ell = int(sqrt(real(orbital, rp)))
            m = orbital - ell*ell - ell
            call this%write_fixed_u_csv_row(unit, 'basis_map', this%fixed_u_trace_iteration, 'RS_GREEN', ia - 1, &
               this%lattice%iz(ia), ell, spin, -1, -1, i, m, 0, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp)
         end do
         do ie = 1, npts + 2
            do i = 1, nmat
               do j = 1, nmat
                  z = this%green%g0(i, j, ie, iloc)
                  call this%write_fixed_u_csv_row(unit, 'green_function', this%fixed_u_trace_iteration, 'RS_GREEN', &
                     ia - 1, this%lattice%iz(ia), -1, -1, ie - 1, -1, i - 1, j - 1, 0, &
                     quadrature_weight(ie), this%en%ene(ie), 0.0_rp, real(z, rp), aimag(z))
               end do
            end do
            call this%write_fixed_u_csv_row(unit, 'fermi_step', this%fixed_u_trace_iteration, 'RS_GREEN', &
               ia - 1, this%lattice%iz(ia), -1, -1, ie - 1, -1, -1, -1, 0, quadrature_weight(ie), &
               this%en%ene(ie), merge(1.0_rp, 0.0_rp, this%en%ene(ie) <= this%en%fermi), 0.0_rp, 0.0_rp)
         end do
      end do
      deallocate(quadrature_weight)
   end procedure write_fixed_u_rs_green

   module procedure write_hubbard_run_summary
      integer :: itype, l
      logical :: active, adopted
      character(len=4096) :: line

      active = .false.
      do itype = 1, this%lattice%ntype
         if (allocated(this%lattice%symbolic_atoms(itype)%potential%hubbard_u)) &
            active = active .or. any(abs(this%lattice%symbolic_atoms(itype)%potential%hubbard_u) > 1.0e-10_rp)
         if (allocated(this%lattice%symbolic_atoms(itype)%potential%hubbard_j)) &
            active = active .or. any(abs(this%lattice%symbolic_atoms(itype)%potential%hubbard_j) > 1.0e-10_rp)
      end do
      if (.not. active .or. rank /= 0) return

      adopted = this%is_eligible_local_hubbard() .and. this%hamiltonian%hubbard_radial_moments == 'local'
      line = 'HUBBARD_RUN form='//trim(this%hamiltonian%hubbard_u_potential_form)// &
         ' radial_mode='//trim(this%hamiltonian%hubbard_radial_moments)// &
         ' radial_source='//trim(radial_mode_source)// &
         ' radial_reason='//trim(radial_mode_reason)// &
         ' energy_routing='//merge('adopted', 'off    ', adopted)//' provenance=PASS'
      do itype = 1, this%lattice%ntype
         if (.not. allocated(this%lattice%symbolic_atoms(itype)%potential%hubbard_u) .or. &
             .not. allocated(this%lattice%symbolic_atoms(itype)%potential%hubbard_j)) cycle
         line = trim(line)//' type='//fmt('i0', itype)//' label='// &
            trim(this%lattice%symbolic_atoms(itype)%element%symbol)//' U_eV='
         do l = 1, min(size(this%lattice%symbolic_atoms(itype)%potential%hubbard_u), &
                       size(this%lattice%symbolic_atoms(itype)%potential%hubbard_j))
            if (l > 1) line = trim(line)//','
            line = trim(line)//trim(adjustl(fmt('f10.6', this%lattice%symbolic_atoms(itype)%potential%hubbard_u(l)*ry2ev)))// &
               '/'//trim(adjustl(fmt('f10.6', this%lattice%symbolic_atoms(itype)%potential%hubbard_j(l)*ry2ev)))
         end do
      end do
      call g_logger%info(trim(line), __FILE__, __LINE__)
   end procedure write_hubbard_run_summary


end submodule ldau_diagnostics
