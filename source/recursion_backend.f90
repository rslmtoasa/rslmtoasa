module recursion_backend_mod

   use iso_c_binding, only: c_ptr, c_null_ptr, c_char, c_int, c_double, c_double_complex, c_associated, c_null_char
   use precision_mod, only: rp
   use hamiltonian_mod
   use control_mod
   use lattice_mod
   use block_sparse_operator_mod
   use logger_mod, only: g_logger
   use string_mod, only: lower, log2str, int2str, real2str
   use basis_mod, only: nb
   implicit none

   private

   type, bind(C), public :: recursion_backend_c_handle
      type(c_ptr) :: value
   end type recursion_backend_c_handle

   interface
      function rslmto_backend_plugin_open(path, errbuf, errbuf_len) bind(C)
         import :: c_ptr, c_char, c_int
         type(c_ptr) :: rslmto_backend_plugin_open
         character(kind=c_char), intent(in) :: path(*)
         character(kind=c_char), intent(out) :: errbuf(*)
         integer(c_int), value :: errbuf_len
      end function rslmto_backend_plugin_open

      subroutine rslmto_backend_plugin_close(runtime_handle) bind(C)
         import :: c_ptr
         type(c_ptr), value :: runtime_handle
      end subroutine rslmto_backend_plugin_close

      function rslmto_backend_plugin_create_backend(runtime_handle, precision_mode, block_dim, n_sites, library_name, errbuf, errbuf_len) bind(C)
         import :: c_ptr, c_char, c_int
         type(c_ptr) :: rslmto_backend_plugin_create_backend
         type(c_ptr), value :: runtime_handle
         character(kind=c_char), intent(in) :: precision_mode(*)
         integer(c_int), value :: block_dim, n_sites
         character(kind=c_char), intent(in) :: library_name(*)
         character(kind=c_char), intent(out) :: errbuf(*)
         integer(c_int), value :: errbuf_len
      end function rslmto_backend_plugin_create_backend

      function rslmto_backend_plugin_destroy_backend(runtime_handle, backend_handle, errbuf, errbuf_len) bind(C)
         import :: c_ptr, c_char, c_int
         integer(c_int) :: rslmto_backend_plugin_destroy_backend
         type(c_ptr), value :: runtime_handle, backend_handle
         character(kind=c_char), intent(out) :: errbuf(*)
         integer(c_int), value :: errbuf_len
      end function rslmto_backend_plugin_destroy_backend

      function rslmto_backend_plugin_upload_operator(runtime_handle, backend_handle, op_name, block_dim, n_sites, nnzb, row_ptr, col_ind, site_types, blocks, errbuf, errbuf_len) bind(C)
         import :: c_ptr, c_char, c_int, c_double_complex
         integer(c_int) :: rslmto_backend_plugin_upload_operator
         type(c_ptr), value :: runtime_handle, backend_handle
         character(kind=c_char), intent(in) :: op_name(*)
         integer(c_int), value :: block_dim, n_sites, nnzb, errbuf_len
         integer(c_int), intent(in) :: row_ptr(*), col_ind(*), site_types(*)
         complex(c_double_complex), intent(in) :: blocks(block_dim, block_dim, nnzb)
         character(kind=c_char), intent(out) :: errbuf(*)
      end function rslmto_backend_plugin_upload_operator

      function rslmto_backend_plugin_apply_operator(runtime_handle, backend_handle, op_name, trans_mode, block_dim, n_sites, psi_in, active_in, psi_out, active_out, errbuf, errbuf_len) bind(C)
         import :: c_ptr, c_char, c_int, c_double_complex
         integer(c_int) :: rslmto_backend_plugin_apply_operator
         type(c_ptr), value :: runtime_handle, backend_handle
         character(kind=c_char), intent(in) :: op_name(*)
         character(kind=c_char), value :: trans_mode
         integer(c_int), value :: block_dim, n_sites, errbuf_len
         complex(c_double_complex), intent(in) :: psi_in(block_dim, block_dim, n_sites)
         integer(c_int), intent(in) :: active_in(*)
         complex(c_double_complex), intent(out) :: psi_out(block_dim, block_dim, n_sites)
         integer(c_int), intent(out) :: active_out(*)
         character(kind=c_char), intent(out) :: errbuf(*)
      end function rslmto_backend_plugin_apply_operator

      function rslmto_backend_plugin_run_scalar_lanczos(runtime_handle, backend_handle, llmax, n_targets, target_sites, hoh_enabled, a_out, b2_out, errbuf, errbuf_len) bind(C)
         import :: c_ptr, c_char, c_int, c_double
         integer(c_int) :: rslmto_backend_plugin_run_scalar_lanczos
         type(c_ptr), value :: runtime_handle, backend_handle
         integer(c_int), value :: llmax, n_targets, hoh_enabled, errbuf_len
         integer(c_int), intent(in) :: target_sites(*)
         real(c_double), intent(out) :: a_out(*)
         real(c_double), intent(out) :: b2_out(*)
         character(kind=c_char), intent(out) :: errbuf(*)
      end function rslmto_backend_plugin_run_scalar_lanczos

      function rslmto_backend_plugin_run_block_lanczos(runtime_handle, backend_handle, llmax, n_workflows, site_i, site_j, variant_id, hoh_enabled, a_b_out, b2_b_out, a_diag_out, b2_diag_out, errbuf, errbuf_len) bind(C)
         import :: c_ptr, c_char, c_int, c_double, c_double_complex
         integer(c_int) :: rslmto_backend_plugin_run_block_lanczos
         type(c_ptr), value :: runtime_handle, backend_handle
         integer(c_int), value :: llmax, n_workflows, hoh_enabled, errbuf_len
         integer(c_int), intent(in) :: site_i(*), site_j(*), variant_id(*)
         complex(c_double_complex), intent(out) :: a_b_out(*)
         complex(c_double_complex), intent(out) :: b2_b_out(*)
         real(c_double), intent(out) :: a_diag_out(*)
         real(c_double), intent(out) :: b2_diag_out(*)
         character(kind=c_char), intent(out) :: errbuf(*)
      end function rslmto_backend_plugin_run_block_lanczos

      function rslmto_backend_plugin_run_block_chebyshev(runtime_handle, backend_handle, llmax, n_workflows, site_i, site_j, variant_id, hoh_enabled, a_scale, b_shift, mu_out, errbuf, errbuf_len) bind(C)
         import :: c_ptr, c_char, c_int, c_double, c_double_complex
         integer(c_int) :: rslmto_backend_plugin_run_block_chebyshev
         type(c_ptr), value :: runtime_handle, backend_handle
         integer(c_int), value :: llmax, n_workflows, hoh_enabled, errbuf_len
         integer(c_int), intent(in) :: site_i(*), site_j(*), variant_id(*)
         real(c_double), value :: a_scale, b_shift
         complex(c_double_complex), intent(out) :: mu_out(*)
         character(kind=c_char), intent(out) :: errbuf(*)
      end function rslmto_backend_plugin_run_block_chebyshev

      function rslmto_backend_plugin_run_transport_stochastic(runtime_handle, backend_handle, llmax, loop_over, calc_mode, hoh_enabled, type_sites, a_scale, b_shift, mu_nm_out, errbuf, errbuf_len) bind(C)
         import :: c_ptr, c_char, c_int, c_double, c_double_complex
         integer(c_int) :: rslmto_backend_plugin_run_transport_stochastic
         type(c_ptr), value :: runtime_handle, backend_handle
         integer(c_int), value :: llmax, loop_over, calc_mode, hoh_enabled, errbuf_len
         integer(c_int), intent(in) :: type_sites(*)
         real(c_double), value :: a_scale, b_shift
         complex(c_double_complex), intent(out) :: mu_nm_out(*)
         character(kind=c_char), intent(out) :: errbuf(*)
      end function rslmto_backend_plugin_run_transport_stochastic

      subroutine rslmto_backend_plugin_name(runtime_handle, namebuf, namebuf_len) bind(C)
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: runtime_handle
         character(kind=c_char), intent(out) :: namebuf(*)
         integer(c_int), value :: namebuf_len
      end subroutine rslmto_backend_plugin_name

      subroutine rslmto_backend_plugin_capabilities(runtime_handle, capbuf, capbuf_len) bind(C)
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: runtime_handle
         character(kind=c_char), intent(out) :: capbuf(*)
         integer(c_int), value :: capbuf_len
      end subroutine rslmto_backend_plugin_capabilities
   end interface

   type, public :: recursion_backend_context
      class(hamiltonian), pointer :: hamiltonian => null()
      class(control), pointer :: control => null()
      class(lattice), pointer :: lattice => null()
      type(block_sparse_operator) :: h_plain
      type(block_sparse_operator) :: h_full
      type(block_sparse_operator) :: h_core
      type(block_sparse_operator) :: h_overlap
      type(block_sparse_operator) :: h_shift
      type(block_sparse_operator) :: v_a
      type(block_sparse_operator) :: v_b
      type(block_sparse_operator) :: vo_a
      type(block_sparse_operator) :: vo_b
      type(c_ptr) :: plugin_runtime = c_null_ptr
      type(recursion_backend_c_handle) :: plugin_handle = recursion_backend_c_handle(c_null_ptr)
      character(len=32) :: backend_kind = 'cpu_reference'
      character(len=256) :: backend_plugin = ''
      character(len=128) :: backend_library = ''
      character(len=128) :: plugin_name = ''
      character(len=256) :: plugin_capabilities = ''
      character(len=32) :: precision_mode = 'complex_fp64'
      logical :: export_enabled = .false.
      logical :: validate_roundtrip = .false.
      logical :: initialized = .false.
      logical :: exports_written = .false.
      logical :: validation_done = .false.
      logical :: plugin_operators_uploaded = .false.
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: sync_from_hamiltonian
      procedure :: apply_hamiltonian
      procedure :: apply_hamiltonian_hoh
      procedure :: apply_velocity
      procedure :: apply_velocity_hoh
      procedure :: validate_against_legacy
      procedure :: ensure_plugin_ready
      procedure :: upload_plugin_operator
      procedure :: apply_plugin_operator
      procedure :: plugin_supports
      procedure :: log_backend_dispatch
      procedure :: run_scalar_lanczos
      procedure :: run_block_lanczos
      procedure :: run_block_chebyshev
      procedure :: run_transport_stochastic
   end type recursion_backend_context

contains

   subroutine initialize(this, hamiltonian_obj, control_obj, lattice_obj)
      class(recursion_backend_context), intent(inout) :: this
      class(hamiltonian), target, intent(inout) :: hamiltonian_obj
      class(control), target, intent(in) :: control_obj
      class(lattice), target, intent(in) :: lattice_obj

      this%hamiltonian => hamiltonian_obj
      this%control => control_obj
      this%lattice => lattice_obj
      this%backend_kind = lower(trim(control_obj%recur_backend))
      this%backend_plugin = trim(control_obj%recur_backend_plugin)
      this%backend_library = trim(control_obj%recur_backend_library)
      this%precision_mode = lower(trim(control_obj%recur_precision))
      this%export_enabled = control_obj%export_hamiltonian
      this%validate_roundtrip = control_obj%validate_backend_roundtrip
      this%plugin_operators_uploaded = .false.
      this%initialized = .true.

      if (this%backend_kind == 'external_plugin') call this%ensure_plugin_ready()
   end subroutine initialize

   subroutine finalize(this)
      class(recursion_backend_context), intent(inout) :: this

      ! Treat plugin state as process-lifetime for now. Explicit teardown of
      ! dynamically loaded backends has been unstable at shutdown, especially
      ! once threaded plugin code has been exercised.
      this%plugin_handle%value = c_null_ptr
      if (c_associated(this%plugin_runtime)) then
         call rslmto_backend_plugin_close(this%plugin_runtime)
         this%plugin_runtime = c_null_ptr
      end if
      this%plugin_name = ''
      this%plugin_capabilities = ''
      this%plugin_operators_uploaded = .false.
   end subroutine finalize

   subroutine sync_from_hamiltonian(this)
      class(recursion_backend_context), intent(inout) :: this
      character(len=:), allocatable :: base

      if (.not. associated(this%hamiltonian)) return

      if (.not. this%h_plain%is_ready()) call this%hamiltonian%export_block_operator('hamiltonian', this%h_plain)
      if (.not. this%h_full%is_ready()) call this%hamiltonian%export_block_operator('hamiltonian_hoh', this%h_full)
      if (.not. this%h_core%is_ready()) call this%hamiltonian%export_block_operator('hamiltonian_core', this%h_core)
      if (.not. this%h_overlap%is_ready()) call this%hamiltonian%export_block_operator('hamiltonian_overlap', this%h_overlap)
      if (.not. this%h_shift%is_ready()) call this%hamiltonian%export_block_operator('hamiltonian_shift', this%h_shift)

      if (this%backend_kind == 'external_plugin' .and. .not. this%plugin_operators_uploaded) then
         call this%ensure_plugin_ready()
         call this%upload_plugin_operator('hamiltonian', this%h_plain)
         call this%upload_plugin_operator('hamiltonian_hoh', this%h_full)
         call this%upload_plugin_operator('hamiltonian_core', this%h_core)
         call this%upload_plugin_operator('hamiltonian_overlap', this%h_overlap)
         call this%upload_plugin_operator('hamiltonian_shift', this%h_shift)
         this%plugin_operators_uploaded = .true.
      end if

      if (this%export_enabled .and. .not. this%exports_written) then
         call this%hamiltonian%block_to_sparse()
         base = trim(this%control%export_hamiltonian_path)
         if (len_trim(base) > 0) then
            call this%h_plain%write_text(trim(base)//'_ham.dat')
            call this%h_full%write_text(trim(base)//'_ham_hoh.dat')
         end if
         this%exports_written = .true.
      end if

      if (this%validate_roundtrip .and. .not. this%validation_done) then
         this%validation_done = .true.
         call this%validate_against_legacy()
      end if
   end subroutine sync_from_hamiltonian

   subroutine apply_hamiltonian(this, psi_in, psi_out, active_in, active_out, a, b)
      class(recursion_backend_context), intent(inout) :: this
      complex(rp), intent(in) :: psi_in(:, :, :)
      complex(rp), intent(out) :: psi_out(:, :, :)
      integer, intent(in) :: active_in(:)
      integer, intent(out) :: active_out(:)
      real(rp), intent(in) :: a, b

      if (this%backend_kind == 'cpu_reference') then
         call legacy_apply_hamiltonian(this%hamiltonian, this%lattice, psi_in, psi_out, active_in, active_out)
         psi_out(:, :, :) = (psi_out(:, :, :) - b*psi_in(:, :, :))/a
         return
      end if

      if (this%backend_kind == 'external_plugin') then
         call this%sync_from_hamiltonian()
         call this%apply_plugin_operator('hamiltonian', 'n', psi_in, psi_out, active_in, active_out)
         psi_out(:, :, :) = (psi_out(:, :, :) - b*psi_in(:, :, :))/a
         return
      end if

      call this%sync_from_hamiltonian()
      call apply_block_operator_cpu(this%h_plain, psi_in, psi_out, active_in, active_out, 'n')
      psi_out(:, :, :) = (psi_out(:, :, :) - b*psi_in(:, :, :))/a
   end subroutine apply_hamiltonian

   subroutine apply_hamiltonian_hoh(this, psi_in, psi_out, active_in, active_out, a, b)
      class(recursion_backend_context), intent(inout) :: this
      complex(rp), intent(in) :: psi_in(:, :, :)
      complex(rp), intent(out) :: psi_out(:, :, :)
      integer, intent(in) :: active_in(:)
      integer, intent(out) :: active_out(:)
      real(rp), intent(in) :: a, b
      complex(rp), allocatable :: tmp_core(:, :, :), tmp_overlap(:, :, :), tmp_shift(:, :, :)
      integer, allocatable :: active_core(:), active_overlap(:), active_shift(:)
      integer :: i

      if (this%backend_kind == 'cpu_reference') then
         call legacy_apply_hamiltonian_hoh(this%hamiltonian, this%lattice, psi_in, psi_out, active_in, active_out)
         psi_out(:, :, :) = (psi_out(:, :, :) - b*psi_in(:, :, :))/a
         return
      end if

      if (this%backend_kind == 'external_plugin') then
         call this%sync_from_hamiltonian()
         allocate(tmp_core(nb, nb, size(psi_in, 3)), tmp_overlap(nb, nb, size(psi_in, 3)), tmp_shift(nb, nb, size(psi_in, 3)))
         allocate(active_core(size(active_in)), active_overlap(size(active_in)), active_shift(size(active_in)))
         call this%apply_plugin_operator('hamiltonian_core', 'n', psi_in, tmp_core, active_in, active_core)
         call this%apply_plugin_operator('hamiltonian_overlap', 'n', tmp_core, tmp_overlap, active_core, active_overlap)
         call this%apply_plugin_operator('hamiltonian_shift', 'n', psi_in, tmp_shift, active_in, active_shift)
         psi_out(:, :, :) = tmp_core(:, :, :) - tmp_overlap(:, :, :) + tmp_shift(:, :, :)
         active_out(:) = 0
         do i = 1, size(active_out)
            if (active_core(i) /= 0 .or. active_overlap(i) /= 0 .or. active_shift(i) /= 0) active_out(i) = 1
         end do
         psi_out(:, :, :) = (psi_out(:, :, :) - b*psi_in(:, :, :))/a
         deallocate(tmp_core, tmp_overlap, tmp_shift, active_core, active_overlap, active_shift)
         return
      end if

      call this%sync_from_hamiltonian()

      allocate(tmp_core(nb, nb, size(psi_in, 3)), tmp_overlap(nb, nb, size(psi_in, 3)), tmp_shift(nb, nb, size(psi_in, 3)))
      allocate(active_core(size(active_in)), active_overlap(size(active_in)), active_shift(size(active_in)))

      call apply_block_operator_cpu(this%h_core, psi_in, tmp_core, active_in, active_core, 'n')
      call apply_block_operator_cpu(this%h_overlap, tmp_core, tmp_overlap, active_core, active_overlap, 'n')
      call apply_block_operator_cpu(this%h_shift, psi_in, tmp_shift, active_in, active_shift, 'n')

      psi_out(:, :, :) = tmp_core(:, :, :) - tmp_overlap(:, :, :) + tmp_shift(:, :, :)
      active_out(:) = 0
      do i = 1, size(active_out)
         if (active_core(i) /= 0 .or. active_overlap(i) /= 0 .or. active_shift(i) /= 0) active_out(i) = 1
      end do
      psi_out(:, :, :) = (psi_out(:, :, :) - b*psi_in(:, :, :))/a

      deallocate(tmp_core, tmp_overlap, tmp_shift, active_core, active_overlap, active_shift)
   end subroutine apply_hamiltonian_hoh

   subroutine apply_velocity(this, operator_name, v_op, psi_in, psi_out, active_in, active_out, trans_mode)
      class(recursion_backend_context), intent(inout) :: this
      character(len=*), intent(in) :: operator_name
      complex(rp), intent(in) :: v_op(:, :, :, :)
      complex(rp), intent(in) :: psi_in(:, :, :)
      complex(rp), intent(out) :: psi_out(:, :, :)
      integer, intent(in) :: active_in(:)
      integer, intent(out) :: active_out(:)
      character(len=1), intent(in) :: trans_mode
      type(block_sparse_operator) :: op

      if (this%backend_kind == 'external_plugin') then
         call this%hamiltonian%export_velocity_operator(operator_name, v_op, op)
         call this%upload_plugin_operator(operator_name, op)
         call this%apply_plugin_operator(operator_name, trans_mode, psi_in, psi_out, active_in, active_out)
         call op%clear()
         return
      end if

      call this%hamiltonian%export_velocity_operator(operator_name, v_op, op)
      call apply_block_operator_cpu(op, psi_in, psi_out, active_in, active_out, trans_mode)
      call op%clear()
   end subroutine apply_velocity

   subroutine apply_velocity_hoh(this, velocity_name, overlap_name, v_op, vo_op, psi_in, psi_out, active_in, active_out)
      class(recursion_backend_context), intent(inout) :: this
      character(len=*), intent(in) :: velocity_name, overlap_name
      complex(rp), intent(in) :: v_op(:, :, :, :)
      complex(rp), intent(in) :: vo_op(:, :, :, :)
      complex(rp), intent(in) :: psi_in(:, :, :)
      complex(rp), intent(out) :: psi_out(:, :, :)
      integer, intent(in) :: active_in(:)
      integer, intent(out) :: active_out(:)
      type(block_sparse_operator) :: op_v, op_vo
      complex(rp), allocatable :: tmp_h(:,:,:), tmp_v(:,:,:), tmp_vo(:,:,:)
      integer, allocatable :: tmp_active(:), tmp_active2(:)

      if (this%backend_kind == 'external_plugin') then
         call this%sync_from_hamiltonian()
         call this%hamiltonian%export_velocity_operator(velocity_name, v_op, op_v)
         call this%hamiltonian%export_velocity_operator(overlap_name, vo_op, op_vo)
         call this%upload_plugin_operator(velocity_name, op_v)
         call this%upload_plugin_operator(overlap_name, op_vo)
         allocate(tmp_h(nb, nb, size(psi_in, 3)), tmp_v(nb, nb, size(psi_in, 3)), tmp_vo(nb, nb, size(psi_in, 3)))
         allocate(tmp_active(size(active_in)), tmp_active2(size(active_in)))
         call this%apply_plugin_operator(velocity_name, 'n', psi_in, tmp_v, active_in, tmp_active)
         call this%apply_plugin_operator('hamiltonian_core', 'n', psi_in, tmp_h, active_in, tmp_active2)
         call this%apply_plugin_operator(overlap_name, 'n', tmp_h, tmp_vo, tmp_active2, active_out)
         psi_out(:, :, :) = tmp_v(:, :, :) - tmp_vo(:, :, :)
         deallocate(tmp_h, tmp_v, tmp_vo, tmp_active, tmp_active2)
         call op_v%clear()
         call op_vo%clear()
         return
      end if

      call this%sync_from_hamiltonian()
      call this%hamiltonian%export_velocity_operator(velocity_name, v_op, op_v)
      call this%hamiltonian%export_velocity_operator(overlap_name, vo_op, op_vo)

      allocate(tmp_h(nb, nb, size(psi_in, 3)), tmp_v(nb, nb, size(psi_in, 3)), tmp_vo(nb, nb, size(psi_in, 3)))
      allocate(tmp_active(size(active_in)), tmp_active2(size(active_in)))

      call apply_block_operator_cpu(op_v, psi_in, tmp_v, active_in, tmp_active, 'n')
      call apply_block_operator_cpu(this%h_core, psi_in, tmp_h, active_in, tmp_active2, 'n')
      call apply_block_operator_cpu(op_vo, tmp_h, tmp_vo, tmp_active2, active_out, 'n')

      psi_out(:, :, :) = tmp_v(:, :, :) - tmp_vo(:, :, :)

      deallocate(tmp_h, tmp_v, tmp_vo, tmp_active, tmp_active2)
      call op_v%clear()
      call op_vo%clear()
   end subroutine apply_velocity_hoh

   subroutine validate_against_legacy(this)
      class(recursion_backend_context), intent(inout) :: this
      complex(rp), allocatable :: psi_in(:, :, :), backend_out(:, :, :), legacy_out(:, :, :)
      integer, allocatable :: active_in(:), active_out(:), legacy_active(:)
      real(rp) :: diff_plain, diff_full

      if (.not. associated(this%lattice)) return

      allocate(psi_in(nb, nb, this%lattice%kk), backend_out(nb, nb, this%lattice%kk), legacy_out(nb, nb, this%lattice%kk))
      allocate(active_in(this%lattice%kk), active_out(this%lattice%kk), legacy_active(this%lattice%kk))

      psi_in(:, :, :) = (0.0_rp, 0.0_rp)
      active_in(:) = 0
      active_in(1:min(2, this%lattice%kk)) = 1
      psi_in(:, :, 1) = 0.0_rp
      psi_in(1:min(nb, size(psi_in, 1)), 1:min(nb, size(psi_in, 2)), 1) = eye_block(min(nb, size(psi_in, 1)))
      if (this%lattice%kk > 1) psi_in(:, :, 2) = cmplx(0.125_rp, -0.0625_rp, kind=rp)

      call this%apply_hamiltonian(psi_in, backend_out, active_in, active_out, 1.0_rp, 0.0_rp)
      call legacy_apply_hamiltonian(this%hamiltonian, this%lattice, psi_in, legacy_out, active_in, legacy_active)
      diff_plain = maxval(abs(backend_out - legacy_out))

      call this%apply_hamiltonian_hoh(psi_in, backend_out, active_in, active_out, 1.0_rp, 0.0_rp)
      call legacy_apply_hamiltonian_hoh(this%hamiltonian, this%lattice, psi_in, legacy_out, active_in, legacy_active)
      diff_full = maxval(abs(backend_out - legacy_out))

      if (diff_plain > 1.0e-10_rp .or. diff_full > 1.0e-10_rp) then
         call g_logger%fatal('Backend roundtrip validation failed: plain='//trim(real2str(diff_plain))// &
                             ' full='//trim(real2str(diff_full)), __FILE__, __LINE__)
      else
         call g_logger%info('Backend roundtrip validation passed: plain='//trim(real2str(diff_plain))// &
                            ' full='//trim(real2str(diff_full)), __FILE__, __LINE__)
      end if

      deallocate(psi_in, backend_out, legacy_out, active_in, active_out, legacy_active)
   end subroutine validate_against_legacy

   subroutine ensure_plugin_ready(this)
      class(recursion_backend_context), intent(inout) :: this
      character(kind=c_char) :: errbuf(512), namebuf(256), capbuf(512)
      character(kind=c_char), allocatable :: c_path(:), c_precision(:), c_library(:)
      character(len=:), allocatable :: plugin_path, library_name

      if (this%backend_kind /= 'external_plugin') return
      if (c_associated(this%plugin_runtime) .and. c_associated(this%plugin_handle%value)) return
      plugin_path = trim(this%backend_plugin)
      library_name = trim(this%backend_library)

      if (len_trim(plugin_path) == 0 .and. looks_like_shared_library_path(trim(this%backend_library))) then
         plugin_path = trim(this%backend_library)
         library_name = 'builtin'
      end if
      if (len_trim(plugin_path) == 0) then
         call g_logger%fatal('external_plugin backend requires control%recur_backend_plugin to point to a shared library', __FILE__, __LINE__)
      end if
      if (len_trim(library_name) == 0) then
         library_name = 'builtin'
      end if

      call to_c_string(trim(plugin_path), c_path)
      call to_c_string(trim(this%precision_mode), c_precision)
      call to_c_string(trim(library_name), c_library)

      this%plugin_runtime = rslmto_backend_plugin_open(c_path, errbuf, int(size(errbuf), c_int))
      if (.not. c_associated(this%plugin_runtime)) then
         call g_logger%fatal('Failed to open backend plugin: '//trim(from_c_buffer(errbuf)), __FILE__, __LINE__)
      end if

      this%plugin_handle%value = rslmto_backend_plugin_create_backend(this%plugin_runtime, c_precision, int(nb, c_int), int(this%lattice%kk, c_int), c_library, errbuf, int(size(errbuf), c_int))
      if (.not. c_associated(this%plugin_handle%value)) then
         call g_logger%fatal('Failed to create backend plugin instance: '//trim(from_c_buffer(errbuf)), __FILE__, __LINE__)
      end if

      call rslmto_backend_plugin_name(this%plugin_runtime, namebuf, int(size(namebuf), c_int))
      call rslmto_backend_plugin_capabilities(this%plugin_runtime, capbuf, int(size(capbuf), c_int))
      this%plugin_name = trim(from_c_buffer(namebuf))
      this%plugin_capabilities = trim(from_c_buffer(capbuf))
      this%backend_plugin = trim(plugin_path)
      this%backend_library = trim(library_name)
      call g_logger%info('Loaded backend plugin: '//trim(this%plugin_name)//' from '//trim(this%backend_plugin), __FILE__, __LINE__)
      call g_logger%info('Backend plugin library selector: '//trim(this%backend_library), __FILE__, __LINE__)
      call g_logger%info('Backend plugin capabilities: '//trim(this%plugin_capabilities), __FILE__, __LINE__)
   end subroutine ensure_plugin_ready

   subroutine upload_plugin_operator(this, operator_name, op)
      class(recursion_backend_context), intent(inout) :: this
      character(len=*), intent(in) :: operator_name
      type(block_sparse_operator), intent(in) :: op
      character(kind=c_char) :: errbuf(512)
      character(kind=c_char), allocatable :: c_name(:)
      integer(c_int), allocatable :: row_ptr_c(:), col_ind_c(:), site_types_c(:)
      integer(c_int) :: ierr

      if (this%backend_kind /= 'external_plugin') return
      call this%ensure_plugin_ready()
      call g_logger%info('Uploading backend operator '//trim(operator_name)//' with nnzb='//trim(int2str(op%nnzb))// &
                         ' block_dim='//trim(int2str(op%block_dim))//' n_sites='//trim(int2str(op%n_sites)), __FILE__, __LINE__)
      call to_c_string(trim(operator_name), c_name)

      allocate(row_ptr_c(size(op%row_ptr)), col_ind_c(size(op%col_ind)), site_types_c(size(op%site_types)))
      row_ptr_c(:) = int(op%row_ptr(:), c_int)
      col_ind_c(:) = int(op%col_ind(:), c_int)
      site_types_c(:) = int(op%site_types(:), c_int)

      ierr = rslmto_backend_plugin_upload_operator(this%plugin_runtime, this%plugin_handle%value, c_name, int(op%block_dim, c_int), int(op%n_sites, c_int), int(op%nnzb, c_int), row_ptr_c, col_ind_c, site_types_c, op%blocks, errbuf, int(size(errbuf), c_int))
      if (ierr /= 0) then
         call g_logger%fatal('Backend plugin upload failed for operator '//trim(operator_name)//': '//trim(from_c_buffer(errbuf)), __FILE__, __LINE__)
      end if

      deallocate(row_ptr_c, col_ind_c, site_types_c)
   end subroutine upload_plugin_operator

   subroutine apply_plugin_operator(this, operator_name, trans_mode, psi_in, psi_out, active_in, active_out)
      class(recursion_backend_context), intent(inout) :: this
      character(len=*), intent(in) :: operator_name
      character(len=1), intent(in) :: trans_mode
      complex(rp), intent(in) :: psi_in(:, :, :)
      complex(rp), intent(out) :: psi_out(:, :, :)
      integer, intent(in) :: active_in(:)
      integer, intent(out) :: active_out(:)
      character(kind=c_char) :: errbuf(512)
      character(kind=c_char), allocatable :: c_name(:)
      integer(c_int), allocatable :: active_in_c(:), active_out_c(:)
      integer(c_int) :: ierr

      call this%ensure_plugin_ready()
      call to_c_string(trim(operator_name), c_name)
      allocate(active_in_c(size(active_in)), active_out_c(size(active_out)))
      active_in_c(:) = int(active_in(:), c_int)
      active_out_c(:) = 0_c_int

      ierr = rslmto_backend_plugin_apply_operator(this%plugin_runtime, this%plugin_handle%value, c_name, char(iachar(trans_mode), kind=c_char), &
                                                 int(nb, c_int), int(size(psi_in, 3), c_int), psi_in, active_in_c, psi_out, active_out_c, &
                                                 errbuf, int(size(errbuf), c_int))
      if (ierr /= 0) then
         call g_logger%fatal('Backend plugin apply failed for operator '//trim(operator_name)//': '//trim(from_c_buffer(errbuf)), __FILE__, __LINE__)
      end if

      active_out(:) = int(active_out_c(:))
      deallocate(active_in_c, active_out_c)
   end subroutine apply_plugin_operator

   logical function plugin_supports(this, capability_name)
      class(recursion_backend_context), intent(in) :: this
      character(len=*), intent(in) :: capability_name
      character(len=:), allocatable :: caps, needle

      caps = lower(trim(this%plugin_capabilities))
      needle = lower(trim(capability_name))
      plugin_supports = len_trim(needle) > 0 .and. index(caps, needle) > 0
   end function plugin_supports

   subroutine log_backend_dispatch(this, kernel_name, status, detail)
      class(recursion_backend_context), intent(in) :: this
      character(len=*), intent(in) :: kernel_name, status, detail
      character(len=:), allocatable :: message

      message = 'Backend dispatch ['//trim(kernel_name)//'] '//trim(status)
      if (len_trim(detail) > 0) message = trim(message)//': '//trim(detail)

      select case (lower(trim(status)))
      case ('fallback')
         call g_logger%warning(trim(message), __FILE__, __LINE__)
      case default
         call g_logger%info(trim(message), __FILE__, __LINE__)
      end select
   end subroutine log_backend_dispatch

   subroutine run_scalar_lanczos(this, llmax, target_sites, hoh_enabled, a_out, b2_out, handled)
      class(recursion_backend_context), intent(inout) :: this
      integer, intent(in) :: llmax
      integer, intent(in) :: target_sites(:)
      logical, intent(in) :: hoh_enabled
      real(rp), intent(out) :: a_out(:, :, :)
      real(rp), intent(out) :: b2_out(:, :, :)
      logical, intent(out) :: handled
      character(kind=c_char) :: errbuf(512)
      integer(c_int), allocatable :: target_sites_c(:)
      integer(c_int) :: ierr

      handled = .false.
      if (this%backend_kind /= 'external_plugin') return

      call this%sync_from_hamiltonian()
      if (.not. this%plugin_supports('scalar_lanczos')) then
         call this%log_backend_dispatch('scalar_lanczos', 'fallback', 'plugin lacks scalar_lanczos capability')
         return
      end if

      allocate(target_sites_c(size(target_sites)))
      target_sites_c(:) = int(target_sites(:), c_int)

      ierr = rslmto_backend_plugin_run_scalar_lanczos(this%plugin_runtime, this%plugin_handle%value, int(llmax, c_int), int(size(target_sites), c_int), &
                                                      target_sites_c, merge(1_c_int, 0_c_int, hoh_enabled), a_out, b2_out, errbuf, int(size(errbuf), c_int))
      deallocate(target_sites_c)
      if (ierr /= 0) then
         call this%log_backend_dispatch('scalar_lanczos', 'fallback', 'plugin call failed: '//trim(from_c_buffer(errbuf)))
         return
      end if

      handled = .true.
      call this%log_backend_dispatch('scalar_lanczos', 'plugin', trim(this%plugin_name))
   end subroutine run_scalar_lanczos

   subroutine run_block_lanczos(this, llmax, site_i, site_j, variant_id, hoh_enabled, a_b_out, b2_b_out, a_diag_out, b2_diag_out, handled)
      class(recursion_backend_context), intent(inout) :: this
      integer, intent(in) :: llmax
      integer, intent(in) :: site_i(:), site_j(:), variant_id(:)
      logical, intent(in) :: hoh_enabled
      complex(rp), intent(out) :: a_b_out(:, :, :, :)
      complex(rp), intent(out) :: b2_b_out(:, :, :, :)
      real(rp), intent(out) :: a_diag_out(:, :, :)
      real(rp), intent(out) :: b2_diag_out(:, :, :)
      logical, intent(out) :: handled
      character(kind=c_char) :: errbuf(512)
      integer(c_int), allocatable :: site_i_c(:), site_j_c(:), variant_id_c(:)
      integer(c_int) :: ierr

      handled = .false.
      if (this%backend_kind /= 'external_plugin') return

      call this%sync_from_hamiltonian()
      if (.not. this%plugin_supports('block_lanczos')) then
         call this%log_backend_dispatch('block_lanczos', 'fallback', 'plugin lacks block_lanczos capability')
         return
      end if

      allocate(site_i_c(size(site_i)), site_j_c(size(site_j)), variant_id_c(size(variant_id)))
      site_i_c(:) = int(site_i(:), c_int)
      site_j_c(:) = int(site_j(:), c_int)
      variant_id_c(:) = int(variant_id(:), c_int)

      ierr = rslmto_backend_plugin_run_block_lanczos(this%plugin_runtime, this%plugin_handle%value, int(llmax, c_int), int(size(site_i), c_int), &
                                                     site_i_c, site_j_c, variant_id_c, merge(1_c_int, 0_c_int, hoh_enabled), &
                                                     a_b_out, b2_b_out, a_diag_out, b2_diag_out, errbuf, int(size(errbuf), c_int))
      deallocate(site_i_c, site_j_c, variant_id_c)
      if (ierr /= 0) then
         call this%log_backend_dispatch('block_lanczos', 'fallback', 'plugin call failed: '//trim(from_c_buffer(errbuf)))
         return
      end if

      handled = .true.
      call this%log_backend_dispatch('block_lanczos', 'plugin', trim(this%plugin_name))
   end subroutine run_block_lanczos

   subroutine run_block_chebyshev(this, llmax, site_i, site_j, variant_id, hoh_enabled, a_scale, b_shift, mu_out, handled)
      class(recursion_backend_context), intent(inout) :: this
      integer, intent(in) :: llmax
      integer, intent(in) :: site_i(:), site_j(:), variant_id(:)
      logical, intent(in) :: hoh_enabled
      real(rp), intent(in) :: a_scale, b_shift
      complex(rp), intent(out) :: mu_out(:, :, :, :)
      logical, intent(out) :: handled
      character(kind=c_char) :: errbuf(512)
      integer(c_int), allocatable :: site_i_c(:), site_j_c(:), variant_id_c(:)
      integer(c_int) :: ierr

      handled = .false.
      if (this%backend_kind /= 'external_plugin') return

      call this%sync_from_hamiltonian()
      if (.not. this%plugin_supports('block_chebyshev')) then
         call this%log_backend_dispatch('block_chebyshev', 'fallback', 'plugin lacks block_chebyshev capability')
         return
      end if

      allocate(site_i_c(size(site_i)), site_j_c(size(site_j)), variant_id_c(size(variant_id)))
      site_i_c(:) = int(site_i(:), c_int)
      site_j_c(:) = int(site_j(:), c_int)
      variant_id_c(:) = int(variant_id(:), c_int)

      ierr = rslmto_backend_plugin_run_block_chebyshev(this%plugin_runtime, this%plugin_handle%value, int(llmax, c_int), int(size(site_i), c_int), &
                                                       site_i_c, site_j_c, variant_id_c, merge(1_c_int, 0_c_int, hoh_enabled), &
                                                       real(a_scale, c_double), real(b_shift, c_double), mu_out, errbuf, int(size(errbuf), c_int))
      deallocate(site_i_c, site_j_c, variant_id_c)
      if (ierr /= 0) then
         call this%log_backend_dispatch('block_chebyshev', 'fallback', 'plugin call failed: '//trim(from_c_buffer(errbuf)))
         return
      end if

      handled = .true.
      call this%log_backend_dispatch('block_chebyshev', 'plugin', trim(this%plugin_name))
   end subroutine run_block_chebyshev

   subroutine run_transport_stochastic(this, llmax, calc_mode, type_sites, hoh_enabled, a_scale, b_shift, v_a, vo_a, v_b, vo_b, mu_nm_out, handled)
      class(recursion_backend_context), intent(inout) :: this
      integer, intent(in) :: llmax, calc_mode
      integer, intent(in) :: type_sites(:)
      logical, intent(in) :: hoh_enabled
      real(rp), intent(in) :: a_scale, b_shift
      complex(rp), intent(in) :: v_a(:, :, :, :), vo_a(:, :, :, :), v_b(:, :, :, :), vo_b(:, :, :, :)
      complex(rp), intent(out) :: mu_nm_out(:, :, :, :, :)
      logical, intent(out) :: handled
      character(kind=c_char) :: errbuf(512)
      integer(c_int), allocatable :: type_sites_c(:)
      integer(c_int) :: ierr
      type(block_sparse_operator) :: op_va, op_voa, op_vb, op_vob

      handled = .false.
      if (this%backend_kind /= 'external_plugin') return

      call this%sync_from_hamiltonian()
      if (.not. this%plugin_supports('transport')) then
         call this%log_backend_dispatch('transport_stochastic', 'fallback', 'plugin lacks transport capability')
         return
      end if

      call this%hamiltonian%export_velocity_operator('transport_velocity_a', v_a, op_va)
      call this%hamiltonian%export_velocity_operator('transport_velocity_overlap_a', vo_a, op_voa)
      call this%hamiltonian%export_velocity_operator('transport_velocity_b', v_b, op_vb)
      call this%hamiltonian%export_velocity_operator('transport_velocity_overlap_b', vo_b, op_vob)

      call this%upload_plugin_operator('transport_velocity_a', op_va)
      call this%upload_plugin_operator('transport_velocity_overlap_a', op_voa)
      call this%upload_plugin_operator('transport_velocity_b', op_vb)
      call this%upload_plugin_operator('transport_velocity_overlap_b', op_vob)

      allocate(type_sites_c(size(type_sites)))
      type_sites_c(:) = int(type_sites(:), c_int)
      ierr = rslmto_backend_plugin_run_transport_stochastic(this%plugin_runtime, this%plugin_handle%value, int(llmax, c_int), int(size(type_sites), c_int), &
                                                            int(calc_mode, c_int), merge(1_c_int, 0_c_int, hoh_enabled), type_sites_c, &
                                                            real(a_scale, c_double), real(b_shift, c_double), mu_nm_out, errbuf, int(size(errbuf), c_int))
      deallocate(type_sites_c)
      call op_va%clear()
      call op_voa%clear()
      call op_vb%clear()
      call op_vob%clear()
      if (ierr /= 0) then
         call this%log_backend_dispatch('transport_stochastic', 'fallback', 'plugin call failed: '//trim(from_c_buffer(errbuf)))
         return
      end if

      handled = .true.
      call this%log_backend_dispatch('transport_stochastic', 'plugin', trim(this%plugin_name))
   end subroutine run_transport_stochastic

   subroutine apply_block_operator_cpu(op, psi_in, psi_out, active_in, active_out, trans_mode)
      type(block_sparse_operator), intent(in) :: op
      complex(rp), intent(in) :: psi_in(:, :, :)
      complex(rp), intent(out) :: psi_out(:, :, :)
      integer, intent(in) :: active_in(:)
      integer, intent(out) :: active_out(:)
      character(len=1), intent(in) :: trans_mode
      integer :: row, entry, col
      logical :: touched
      complex(rp), parameter :: one = (1.0_rp, 0.0_rp)

      psi_out(:, :, :) = (0.0_rp, 0.0_rp)
      active_out(:) = 0

      do row = 1, op%n_sites
         touched = .false.
         do entry = op%row_ptr(row), op%row_ptr(row + 1) - 1
            col = op%col_ind(entry)
            if (active_in(col) == 0) cycle
            touched = .true.
            select case (lower(trans_mode))
            case ('c')
               call zgemm('c', 'n', nb, nb, nb, one, op%blocks(:, :, entry), nb, psi_in(:, :, col), nb, one, psi_out(:, :, row), nb)
            case default
               call zgemm('n', 'n', nb, nb, nb, one, op%blocks(:, :, entry), nb, psi_in(:, :, col), nb, one, psi_out(:, :, row), nb)
            end select
         end do
         if (touched) active_out(row) = 1
      end do
   end subroutine apply_block_operator_cpu

   subroutine legacy_apply_hamiltonian(hamiltonian_obj, lattice_obj, psi_in, psi_out, active_in, active_out)
      class(hamiltonian), intent(in) :: hamiltonian_obj
      class(lattice), intent(in) :: lattice_obj
      complex(rp), intent(in) :: psi_in(:, :, :)
      complex(rp), intent(out) :: psi_out(:, :, :)
      integer, intent(in) :: active_in(:)
      integer, intent(out) :: active_out(:)
      integer :: k, ih, nr, ineigh, nnmap, nlimplus1
      complex(rp) :: locham(nb, nb)

      psi_out(:, :, :) = (0.0_rp, 0.0_rp)
      active_out(:) = 0
      nlimplus1 = lattice_obj%nmax + 1

      if (lattice_obj%nmax /= 0) then
         do k = 1, lattice_obj%nmax
            ih = lattice_obj%iz(k)
            nr = lattice_obj%nn(k, 1)
            if (active_in(k) /= 0) then
               locham = hamiltonian_obj%hall(:, :, 1, k) + hamiltonian_obj%lsham(:, :, ih)
               psi_out(:, :, k) = psi_out(:, :, k) + matmul(locham, psi_in(:, :, k))
               active_out(k) = 1
            end if
            if (nr >= 2) then
               do ineigh = 2, nr
                  nnmap = lattice_obj%nn(k, ineigh)
                  if (nnmap /= 0) then
                     if (active_in(nnmap) /= 0) then
                        psi_out(:, :, k) = psi_out(:, :, k) + matmul(hamiltonian_obj%hall(:, :, ineigh, k), psi_in(:, :, nnmap))
                        active_out(k) = 1
                     end if
                  end if
               end do
            end if
         end do
      end if

      do k = nlimplus1, lattice_obj%kk
         ih = lattice_obj%iz(k)
         nr = lattice_obj%nn(k, 1)
         if (active_in(k) /= 0) then
            locham = hamiltonian_obj%ee(:, :, 1, ih) + hamiltonian_obj%lsham(:, :, ih)
            psi_out(:, :, k) = psi_out(:, :, k) + matmul(locham, psi_in(:, :, k))
            active_out(k) = 1
         end if
         if (nr >= 2) then
            do ineigh = 2, nr
               nnmap = lattice_obj%nn(k, ineigh)
               if (nnmap /= 0) then
                  if (active_in(nnmap) /= 0) then
                     psi_out(:, :, k) = psi_out(:, :, k) + matmul(hamiltonian_obj%ee(:, :, ineigh, ih), psi_in(:, :, nnmap))
                     active_out(k) = 1
                  end if
               end if
            end do
         end if
      end do
   end subroutine legacy_apply_hamiltonian

   subroutine legacy_apply_hamiltonian_hoh(hamiltonian_obj, lattice_obj, psi_in, psi_out, active_in, active_out)
      class(hamiltonian), intent(in) :: hamiltonian_obj
      class(lattice), intent(in) :: lattice_obj
      complex(rp), intent(in) :: psi_in(:, :, :)
      complex(rp), intent(out) :: psi_out(:, :, :)
      integer, intent(in) :: active_in(:)
      integer, intent(out) :: active_out(:)
      integer :: k, ih, nr, ineigh, nnmap, nlimplus1
      complex(rp), allocatable :: psi_core(:, :, :), hoh(:, :, :), soc(:, :, :), enu(:, :, :)
      integer, allocatable :: active_core(:)

      allocate(psi_core(nb, nb, lattice_obj%kk), hoh(nb, nb, lattice_obj%kk), soc(nb, nb, lattice_obj%kk), enu(nb, nb, lattice_obj%kk))
      allocate(active_core(size(active_in)))

      psi_core(:, :, :) = (0.0_rp, 0.0_rp)
      hoh(:, :, :) = (0.0_rp, 0.0_rp)
      soc(:, :, :) = (0.0_rp, 0.0_rp)
      enu(:, :, :) = (0.0_rp, 0.0_rp)
      psi_out(:, :, :) = (0.0_rp, 0.0_rp)
      active_out(:) = 0
      active_core(:) = 0
      nlimplus1 = lattice_obj%nmax + 1

      if (lattice_obj%nmax /= 0) then
         do k = 1, lattice_obj%nmax
            ih = lattice_obj%iz(k)
            nr = lattice_obj%nn(k, 1)
            if (active_in(k) /= 0) then
               psi_core(:, :, k) = psi_core(:, :, k) + matmul(hamiltonian_obj%hall(:, :, 1, k), psi_in(:, :, k))
               soc(:, :, k) = soc(:, :, k) + matmul(hamiltonian_obj%lsham(:, :, ih), psi_in(:, :, k))
               enu(:, :, k) = enu(:, :, k) + matmul(hamiltonian_obj%enim(:, :, ih), psi_in(:, :, k))
               active_core(k) = 1
            end if
            if (nr >= 2) then
               do ineigh = 2, nr
                  nnmap = lattice_obj%nn(k, ineigh)
                  if (nnmap /= 0) then
                     if (active_in(nnmap) /= 0) then
                        psi_core(:, :, k) = psi_core(:, :, k) + matmul(hamiltonian_obj%hall(:, :, ineigh, k), psi_in(:, :, nnmap))
                        active_core(k) = 1
                     end if
                  end if
               end do
            end if
         end do
      end if

      do k = nlimplus1, lattice_obj%kk
         ih = lattice_obj%iz(k)
         nr = lattice_obj%nn(k, 1)
         if (active_in(k) /= 0) then
            psi_core(:, :, k) = psi_core(:, :, k) + matmul(hamiltonian_obj%ee(:, :, 1, ih), psi_in(:, :, k))
            soc(:, :, k) = soc(:, :, k) + matmul(hamiltonian_obj%lsham(:, :, ih), psi_in(:, :, k))
            enu(:, :, k) = enu(:, :, k) + matmul(hamiltonian_obj%enim(:, :, ih), psi_in(:, :, k))
            active_core(k) = 1
         end if
         if (nr >= 2) then
            do ineigh = 2, nr
               nnmap = lattice_obj%nn(k, ineigh)
               if (nnmap /= 0) then
                  if (active_in(nnmap) /= 0) then
                     psi_core(:, :, k) = psi_core(:, :, k) + matmul(hamiltonian_obj%ee(:, :, ineigh, ih), psi_in(:, :, nnmap))
                     active_core(k) = 1
                  end if
               end if
            end do
         end if
      end do

      if (lattice_obj%nmax /= 0) then
         do k = 1, lattice_obj%nmax
            nr = lattice_obj%nn(k, 1)
            if (active_core(k) /= 0) then
               hoh(:, :, k) = hoh(:, :, k) + matmul(hamiltonian_obj%hallo(:, :, 1, k), psi_core(:, :, k))
               active_out(k) = 1
            end if
            if (nr >= 2) then
               do ineigh = 2, nr
                  nnmap = lattice_obj%nn(k, ineigh)
                  if (nnmap /= 0) then
                     if (active_core(nnmap) /= 0) then
                        hoh(:, :, k) = hoh(:, :, k) + matmul(hamiltonian_obj%hallo(:, :, ineigh, k), psi_core(:, :, nnmap))
                        active_out(k) = 1
                     end if
                  end if
               end do
            end if
         end do
      end if

      do k = nlimplus1, lattice_obj%kk
         ih = lattice_obj%iz(k)
         nr = lattice_obj%nn(k, 1)
         if (active_core(k) /= 0) then
            hoh(:, :, k) = hoh(:, :, k) + matmul(hamiltonian_obj%eeo(:, :, 1, ih), psi_core(:, :, k))
            active_out(k) = 1
         end if
         if (nr >= 2) then
            do ineigh = 2, nr
               nnmap = lattice_obj%nn(k, ineigh)
               if (nnmap /= 0) then
                  if (active_core(nnmap) /= 0) then
                     hoh(:, :, k) = hoh(:, :, k) + matmul(hamiltonian_obj%eeo(:, :, ineigh, ih), psi_core(:, :, nnmap))
                     active_out(k) = 1
                  end if
               end if
            end do
         end if
      end do

      psi_out(:, :, :) = psi_core(:, :, :) - hoh(:, :, :) + enu(:, :, :) + soc(:, :, :)

      deallocate(psi_core, hoh, soc, enu, active_core)
   end subroutine legacy_apply_hamiltonian_hoh

   pure function eye_block(n) result(block)
      integer, intent(in) :: n
      complex(rp) :: block(n, n)
      integer :: i

      block(:, :) = (0.0_rp, 0.0_rp)
      do i = 1, n
         block(i, i) = (1.0_rp, 0.0_rp)
      end do
   end function eye_block

   subroutine to_c_string(input, output)
      character(len=*), intent(in) :: input
      character(kind=c_char), allocatable, intent(out) :: output(:)
      integer :: i, n

      n = len_trim(input)
      allocate(output(n + 1))
      do i = 1, n
         output(i) = input(i:i)
      end do
      output(n + 1) = c_null_char
   end subroutine to_c_string

   function from_c_buffer(buffer) result(text)
      character(kind=c_char), intent(in) :: buffer(:)
      character(len=:), allocatable :: text
      integer :: i, n

      n = 0
      do i = 1, size(buffer)
         if (buffer(i) == c_null_char) exit
         n = n + 1
      end do
      allocate(character(len=n) :: text)
      do i = 1, n
         text(i:i) = buffer(i)
      end do
   end function from_c_buffer

   logical function looks_like_shared_library_path(value)
      character(len=*), intent(in) :: value
      character(len=:), allocatable :: lowered

      lowered = lower(trim(value))
      looks_like_shared_library_path = index(lowered, '.so') > 0 .or. index(lowered, '.dylib') > 0 .or. index(lowered, '.dll') > 0
   end function looks_like_shared_library_path

end module recursion_backend_mod
