submodule(recursion_mod) recursion_core

contains

   !> @brief Construct a recursion object from Hamiltonian, energy, and sparse helpers.
   !> @details The recursion object is the real-space spectral kernel used by SCF,
   !>          DOS, exchange, conductivity, and orbital-moment workflows. It stores
   !>          pointers to the upstream objects and allocates the coefficient and
   !>          moment arrays sized by the active lattice/control state.
   !> @param[in] hamiltonian_obj Hamiltonian object providing real-space blocks.
   !> @param[in] energy_obj Energy mesh used by DOS and Green reconstruction.
   !> @param[in] sparse_obj Sparse algebra helper for legacy kernels.
   !> @return Initialized recursion object.
   module function constructor(hamiltonian_obj, energy_obj, sparse_obj) result(obj)
      type(recursion) :: obj
      type(hamiltonian), target, intent(in) :: hamiltonian_obj
      type(energy), target, intent(in) :: energy_obj
      type(sparse), intent(in) :: sparse_obj

      obj%hamiltonian => hamiltonian_obj
      allocate(obj%sparse_obj)
      obj%sparse_obj = sparse_obj
      obj%lattice => hamiltonian_obj%charge%lattice
      obj%en => energy_obj
      obj%control => hamiltonian_obj%charge%lattice%control

      call obj%restore_to_default()
   end function constructor

   !> @brief Test whether the CUDA recursion plugin was requested and compiled in.
   !> @details Centralizes the build-time and namelist checks before GPU-specific
   !>          recursion paths attempt to acquire or upload a backend context.
   !> @param[in] this Recursion object whose control flags select plugin usage.
   !> @return True when GPU dispatch is allowed by both build and input state.
   module logical function gpu_plugin_enabled(this)
      class(recursion), intent(in) :: this

      gpu_plugin_enabled = this%control%gpu_plugin
   end function gpu_plugin_enabled

   !> @brief Test whether the selected GPU backend expects periodic/k-space data.
   !> @details Separates real-space BSR/convolution recursion backends from GPU
   !>          backends intended for periodic representations so callers can avoid
   !>          uploading incompatible Hamiltonian layouts.
   !> @param[in] this Recursion object whose control flags encode the GPU backend.
   !> @return True for periodic GPU backend selectors.
   module logical function gpu_periodic_backend(this)
      class(recursion), intent(in) :: this
      integer :: backend

      backend = decode_gpu_backend(this%control%gpu_backend)
      gpu_periodic_backend = (backend == gpu_backend_fft .or. backend == gpu_backend_conv)
   end function gpu_periodic_backend

   !> @brief Decide whether a named recursion feature may run on the GPU plugin.
   !> @details Applies feature-level guards such as scalar-only spin support,
   !>          optional HOH support, backend compatibility, and plugin availability.
   !>          Callers use this before choosing CUDA paths for Lanczos, block
   !>          Lanczos, Chebyshev, stochastic, and intersite workflows.
   !> @param[in] this Recursion object with control/Hamiltonian state.
   !> @param[in] feature Human-readable feature name for logs and diagnostics.
   !> @param[in] require_nsp1 Optional guard for scalar-relativistic-only kernels.
   !> @param[in] allow_hoh Optional guard declaring whether HOH is supported.
   !> @return True when the GPU plugin can be used for the requested feature.
   module logical function gpu_plugin_ready(this, feature, require_nsp1, allow_hoh)
      class(recursion), intent(in) :: this
      character(len=*), intent(in) :: feature
      logical, intent(in), optional :: require_nsp1
      logical, intent(in), optional :: allow_hoh
      integer :: backend
      logical :: nsp1_only, hoh_ok

      gpu_plugin_ready = .false.
      if (.not. gpu_plugin_enabled(this)) return
      nsp1_only = .false.
      if (present(require_nsp1)) nsp1_only = require_nsp1
      hoh_ok = .false.
      if (present(allow_hoh)) hoh_ok = allow_hoh

      if (.not. rsrec_cuda_plugin_compiled()) then
         call g_logger%warning('gpu_plugin requested for '//trim(feature)// &
            ', but the executable was built without ENABLE_CUDA_PLUGIN. '// &
            'Falling back to current recursion path.', __FILE__, __LINE__)
         return
      end if

      ! hoh is supported on the GPU only via the two-sweep block-ELL path:
      ! the Chebyshev moment kernels and the block-Lanczos (Haydock) kernel
      ! (allow_hoh from those call sites). The structured / FFT periodic backend
      ! assumes an ee-only Hamiltonian, so hoh is rejected there regardless.
      if (this%hamiltonian%hoh .and. (.not. hoh_ok .or. gpu_periodic_backend(this))) then
         call g_logger%warning('gpu_plugin does not support hoh in '// &
            trim(feature)//'. Falling back to current recursion path.', __FILE__, __LINE__)
         return
      end if

      if (this%hamiltonian%local_axis) then
         call g_logger%warning('gpu_plugin does not support local-axis '// &
            'rotation in '//trim(feature)//'. Falling back to current recursion path.', &
            __FILE__, __LINE__)
         return
      end if

      if (this%hamiltonian%ccor_2c .and. this%hamiltonian%hoh) then
         call g_logger%warning('gpu_plugin does not support additive ccor_2c with hoh in '// &
            trim(feature)//'. Falling back to current recursion path.', __FILE__, __LINE__)
         return
      end if

      if (nsp1_only .and. this%control%nsp /= 1) then
         call g_logger%warning('gpu_plugin scalar Lanczos path is restricted '// &
            'to nsp=1. Falling back to current recursion path.', __FILE__, __LINE__)
         return
      end if

      backend = decode_gpu_backend(this%control%gpu_backend)
      if (backend < 0) then
         call g_logger%warning('Unknown gpu_backend='//trim(this%control%gpu_backend)// &
            '. Falling back to current recursion path.', __FILE__, __LINE__)
         return
      end if
      if (gpu_periodic_backend(this)) then
         if (.not. this%lattice%pbc) then
            call g_logger%warning('gpu_backend='//trim(this%control%gpu_backend)// &
               ' requires lattice%pbc=.true.. Falling back to current recursion path.', &
               __FILE__, __LINE__)
            return
         end if
         if (.not. this%lattice%b1 .or. .not. this%lattice%b2 .or. .not. this%lattice%b3) then
            call g_logger%warning('gpu_backend='//trim(this%control%gpu_backend)// &
               ' requires full periodic directions b1/b2/b3. Falling back to current recursion path.', &
               __FILE__, __LINE__)
            return
         end if
         if (this%lattice%nmax /= 0) then
            call g_logger%warning('gpu_backend='//trim(this%control%gpu_backend)// &
               ' is limited to ee-only periodic Hamiltonians in v1. Falling back to current recursion path.', &
               __FILE__, __LINE__)
            return
         end if
      end if

      gpu_plugin_ready = .true.
   end function gpu_plugin_ready

   !> @brief Ensure the CUDA backend has the current real-space Hamiltonian uploaded.
   !> @details Acquires the shared plugin context and transfers the Hamiltonian in
   !>          the layout required by the selected real-space GPU backend. The GPU
   !>          context owns its own fingerprinting, so repeated calls with unchanged
   !>          inputs are cheap.
   !> @param[inout] this Recursion object; sets this%gpu_backend when available.
   !> @note Fatal diagnostics are raised by the backend when a requested upload is unsupported.
   module subroutine gpu_plugin_upload_hamiltonian(this)
      class(recursion), intent(inout) :: this
      integer :: backend

      if (this%hamiltonian%ccor_2c .and. .not. this%hamiltonian%hoh) call ensure_ccor_operator_blocks(this)
      backend = decode_gpu_backend(this%control%gpu_backend)
      this%gpu_backend => get_gpu_context(this%lattice%kk, nb, size(this%lattice%nn, 2), &
         this%lattice%ntype, this%lattice%nmax)
      if (this%lattice%nmax > 0) then
         if (this%hamiltonian%hoh) then
            call this%gpu_backend%set_hamiltonian(this%hamiltonian%ee, &
               this%hamiltonian%hall, this%hamiltonian%lsham, this%lattice%nn, &
               this%lattice%iz, this%lattice%nmax, eeo=this%hamiltonian%eeo, &
               hallo=this%hamiltonian%hallo, enim=this%hamiltonian%enim)
         else if (this%hamiltonian%ccor_2c) then
            call this%gpu_backend%set_hamiltonian(this%ee_ccor_work, &
               this%hall_ccor_work, this%hamiltonian%lsham, this%lattice%nn, &
               this%lattice%iz, this%lattice%nmax)
         else
            call this%gpu_backend%set_hamiltonian(this%hamiltonian%ee, &
               this%hamiltonian%hall, this%hamiltonian%lsham, this%lattice%nn, &
               this%lattice%iz, this%lattice%nmax)
         end if
      else
         if (this%hamiltonian%hoh) then
            call this%gpu_backend%set_hamiltonian(this%hamiltonian%ee, &
               lsham=this%hamiltonian%lsham, nn=this%lattice%nn, &
               iz=this%lattice%iz, nmax=this%lattice%nmax, &
               eeo=this%hamiltonian%eeo, enim=this%hamiltonian%enim)
         else if (this%hamiltonian%ccor_2c) then
            call this%gpu_backend%set_hamiltonian(this%ee_ccor_work, &
               lsham=this%hamiltonian%lsham, nn=this%lattice%nn, &
               iz=this%lattice%iz, nmax=this%lattice%nmax)
         else
            call this%gpu_backend%set_hamiltonian(this%hamiltonian%ee, &
               lsham=this%hamiltonian%lsham, nn=this%lattice%nn, &
               iz=this%lattice%iz, nmax=this%lattice%nmax)
         end if
      end if
      call this%gpu_backend%set_periodic_lattice(this%lattice%pbc, this%lattice%n1, &
         this%lattice%n2, this%lattice%n3, this%lattice%a, this%lattice%crd)
      call this%gpu_backend%set_backend(backend)
      if (backend == gpu_backend_bsr) then
         if (this%hamiltonian%ccor_2c .and. .not. this%hamiltonian%hoh) then
            call this%sparse_obj%build_bsr_from_operator(this%hall_ccor_work, this%ee_ccor_work)
         end if
         call this%gpu_backend%upload_bsr(this%sparse_obj)
      end if
   end subroutine gpu_plugin_upload_hamiltonian

   !> @brief Build cached first-order ccor_2c Hamiltonian blocks for fast kernels.
   !> @details Combines the ordinary and ccor_2c real-space operator blocks into
   !>          work arrays used by no-HOH Chebyshev, block-Lanczos, and transport
   !>          kernels that expect one effective hopping operator.
   !> @param[inout] this Recursion object; populates ee_ccor_work and hall_ccor_work.
   !> @note The cached arrays are reused across moment calls until the recursion object is reset.
   module subroutine ensure_ccor_operator_blocks(this)
      class(recursion), intent(inout) :: this

      if (.not. this%hamiltonian%ccor_2c) return
      if (this%hamiltonian%hoh) then
         call g_logger%fatal('ensure_ccor_operator_blocks called for hoh+ccor_2c; use the exact separate-additive path instead.', __FILE__, __LINE__)
      end if

      if (.not. allocated(this%ee_ccor_work)) then
         allocate(this%ee_ccor_work(size(this%hamiltonian%ee, 1), size(this%hamiltonian%ee, 2), &
                                   size(this%hamiltonian%ee, 3), size(this%hamiltonian%ee, 4)))
      else if (any(shape(this%ee_ccor_work) /= shape(this%hamiltonian%ee))) then
         deallocate(this%ee_ccor_work)
         allocate(this%ee_ccor_work(size(this%hamiltonian%ee, 1), size(this%hamiltonian%ee, 2), &
                                   size(this%hamiltonian%ee, 3), size(this%hamiltonian%ee, 4)))
      end if
      this%ee_ccor_work = this%hamiltonian%ee + this%hamiltonian%eecc

      if (.not. allocated(this%hall_ccor_work)) then
         allocate(this%hall_ccor_work(size(this%hamiltonian%hall, 1), size(this%hamiltonian%hall, 2), &
                                     size(this%hamiltonian%hall, 3), size(this%hamiltonian%hall, 4)))
      else if (any(shape(this%hall_ccor_work) /= shape(this%hamiltonian%hall))) then
         deallocate(this%hall_ccor_work)
         allocate(this%hall_ccor_work(size(this%hamiltonian%hall, 1), size(this%hamiltonian%hall, 2), &
                                     size(this%hamiltonian%hall, 3), size(this%hamiltonian%hall, 4)))
      end if
      this%hall_ccor_work = this%hamiltonian%hall + this%hamiltonian%hallcc
   end subroutine ensure_ccor_operator_blocks

   !> @brief Run the optimized no-HOH block-Lanczos kernel for one starting block.
   !> @details Wraps haydock_fast for the common block-recursion path where the
   !>          Hamiltonian is represented by a single hopping sweep. The caller has
   !>          initialized psi_b and copies atemp_b/b2temp_b into public coefficient
   !>          arrays after this routine returns.
   !> @param[inout] this Recursion object with current starting block in psi_b.
   !> @param[in] llmax Number of block-Lanczos steps.
   !> @param[in] use_ccor Use cached ccor_2c operator blocks instead of raw blocks.
   module subroutine block_lanczos_fast_nohoh(this, llmax, use_ccor)
      class(recursion), intent(inout) :: this
      integer, intent(in) :: llmax
      logical, intent(in) :: use_ccor

      if (use_ccor) then
         associate(ee_op => this%ee_ccor_work, hall_op => this%hall_ccor_work)
            call block_lanczos_fast(this%psi_b, llmax, this%atemp_b, this%b2temp_b, &
               ee_op, hall_op, this%hamiltonian%lsham, &
               this%lattice%nn, this%lattice%iz, this%lattice%kk, nb, size(this%lattice%nn, 2), &
               this%lattice%ntype, this%lattice%nmax, &
               trim(this%control%cheb_backend) == 'fast', .false., &
               ee_op, hall_op, this%hamiltonian%lsham)
         end associate
      else
         associate(ee_op => this%hamiltonian%ee, hall_op => this%hamiltonian%hall)
            call block_lanczos_fast(this%psi_b, llmax, this%atemp_b, this%b2temp_b, &
               ee_op, hall_op, this%hamiltonian%lsham, &
               this%lattice%nn, this%lattice%iz, this%lattice%kk, nb, size(this%lattice%nn, 2), &
               this%lattice%ntype, this%lattice%nmax, &
               trim(this%control%cheb_backend) == 'fast', .false., &
               ee_op, hall_op, this%hamiltonian%lsham)
         end associate
      end if
   end subroutine block_lanczos_fast_nohoh

   !> @brief Finalize a recursion object.
   !> @details Releases allocatable coefficient, wavefunction, moment, work, and
   !>          helper arrays owned by the object when it leaves scope.
   !> @param[inout] this Recursion object being finalized.
   module subroutine destructor(this)
      type(recursion) :: this
      nullify (this%gpu_backend)
#ifdef USE_SAFE_ALLOC
      if (allocated(this%a)) call g_safe_alloc%deallocate('recursion.a', this%a)
      if (allocated(this%b2)) call g_safe_alloc%deallocate('recursion.b2', this%b2)
      if (allocated(this%izero)) call g_safe_alloc%deallocate('recursion.izero', this%izero)
      if (allocated(this%izeroll)) call g_safe_alloc%deallocate('recursion.izeroll', this%izeroll)
      if (allocated(this%irlist)) call g_safe_alloc%deallocate('recursion.irlist', this%irlist)
      if (allocated(this%psi)) call g_safe_alloc%deallocate('recursion.psi', this%psi)
      if (allocated(this%psi1)) call g_safe_alloc%deallocate('recursion.psi1', this%psi1)
      if (allocated(this%psi2)) call g_safe_alloc%deallocate('recursion.psi2', this%psi2)
      if (allocated(this%ee_ccor_work)) deallocate (this%ee_ccor_work)
      if (allocated(this%hall_ccor_work)) deallocate (this%hall_ccor_work)
      if (allocated(this%pmn)) call g_safe_alloc%deallocate('recursion.pmn', this%pmn)
      if (allocated(this%v)) call g_safe_alloc%deallocate('recursion.v', this%v)
      if (allocated(this%mu_n)) call g_safe_alloc%deallocate('recursion.mu_n', this%mu_n)
      if (allocated(this%mu_ng)) call g_safe_alloc%deallocate('recursion.mu_ng', this%mu_ng)
      if (allocated(this%atemp)) call g_safe_alloc%deallocate('recursion.atemp', this%atemp)
      if (allocated(this%b2temp)) call g_safe_alloc%deallocate('recursion.b2temp', this%b2temp)
      if (allocated(this%a_b)) call g_safe_alloc%deallocate('recursion.a_b', this%a_b)
      if (allocated(this%b2_b)) call g_safe_alloc%deallocate('recursion.b2_b', this%b2_b)
      if (allocated(this%psi_b)) call g_safe_alloc%deallocate('recursion.psi_b', this%psi_b)
      if (allocated(this%hpsi)) call g_safe_alloc%deallocate('recursion.hpsi', this%hpsi)
      if (allocated(this%hohpsi)) call g_safe_alloc%deallocate('recursion.hohpsi', this%hohpsi)
      if (allocated(this%enupsi)) call g_safe_alloc%deallocate('recursion.enupsi', this%enupsi)
      if (allocated(this%socpsi)) call g_safe_alloc%deallocate('recursion.socpsi', this%socpsi)
      if (allocated(this%atemp_b)) call g_safe_alloc%deallocate('recursion.atemp_b', this%atemp_b)
      if (allocated(this%b2temp_b)) call g_safe_alloc%deallocate('recursion.b2temp_b', this%b2temp_b)
      if (allocated(this%pmn_b)) call g_safe_alloc%deallocate('recursion.pmn_b', this%pmn_b)
      if (allocated(this%t_h)) call g_safe_alloc%deallocate('recursion.t_h', this%t_h)
      if (allocated(this%mu_nm_stochastic)) call g_safe_alloc%deallocate('recursion.mu_nm_stochastic', this%mu_nm_stochastic)
      if (allocated(this%nzero)) call g_safe_alloc%deallocate('recursion.nzero', this%nzero)
      if (allocated(this%mzero)) call g_safe_alloc%deallocate('recursion.mzero', this%mzero)
      if (allocated(this%ndum)) call g_safe_alloc%deallocate('recursion.ndum', this%ndum)
      if (allocated(this%mdum)) call g_safe_alloc%deallocate('recursion.mdum', this%mdum)
#else
      if (allocated(this%a)) deallocate (this%a)
      if (allocated(this%b2)) deallocate (this%b2)
      if (allocated(this%izero)) deallocate (this%izero)
      if (allocated(this%izeroll)) deallocate (this%izeroll)
      if (allocated(this%irlist)) deallocate (this%irlist)
      if (allocated(this%psi)) deallocate (this%psi)
      if (allocated(this%psi1)) deallocate (this%psi1)
      if (allocated(this%psi2)) deallocate (this%psi2)
      if (allocated(this%ee_ccor_work)) deallocate (this%ee_ccor_work)
      if (allocated(this%hall_ccor_work)) deallocate (this%hall_ccor_work)
      if (allocated(this%pmn)) deallocate (this%pmn)
      if (allocated(this%v)) deallocate (this%v)
      if (allocated(this%mu_n)) deallocate (this%mu_n)
      if (allocated(this%mu_ng)) deallocate (this%mu_ng)
      if (allocated(this%atemp)) deallocate (this%atemp)
      if (allocated(this%b2temp)) deallocate (this%b2temp)
      if (allocated(this%a_b)) deallocate (this%a_b)
      if (allocated(this%b2_b)) deallocate (this%b2_b)
      if (allocated(this%psi_b)) deallocate (this%psi_b)
      if (allocated(this%hpsi)) deallocate (this%hpsi)
      if (allocated(this%hohpsi)) deallocate (this%hohpsi)
      if (allocated(this%enupsi)) deallocate (this%enupsi)
      if (allocated(this%socpsi)) deallocate (this%socpsi)
      if (allocated(this%atemp_b)) deallocate (this%atemp_b)
      if (allocated(this%b2temp_b)) deallocate (this%b2temp_b)
      if (allocated(this%pmn_b)) deallocate (this%pmn_b)
      if (allocated(this%t_h)) deallocate(this%t_h)
      if (allocated(this%mu_nm_stochastic)) deallocate(this%mu_nm_stochastic)
      if (allocated(this%nzero)) deallocate (this%nzero)
      if (allocated(this%mzero)) deallocate (this%mzero)
      if (allocated(this%ndum)) deallocate (this%ndum)
      if (allocated(this%mdum)) deallocate (this%mdum)
#endif
   end subroutine destructor

   !> @brief Reset recursion-owned pointers, buffers, and backend state.
   !> @details Restores the object to a reusable baseline before construction or
   !>          teardown. A full reset also releases large coefficient, moment, and
   !>          wavefunction arrays so a new lattice/control shape can be installed.
   !> @param[inout] this Recursion object to reset.
   !> @param[in] full Optional flag requesting deallocation of all owned arrays.
   !> @note Clears the CUDA backend pointer association but does not own the plugin context.
   module subroutine restore_to_default(this, full)
      class(recursion) :: this
      logical, intent(in), optional :: full
      integer :: lmax

      lmax = lmax_basis

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('recursion.a', this%a, (/max(this%lattice%control%llsp, this%lattice%control%lld), nb, this%lattice%nrec, 3/))
      call g_safe_alloc%allocate('recursion.atemp', this%atemp, max(this%lattice%control%llsp, this%lattice%control%lld))
      call g_safe_alloc%allocate('recursion.b2temp', this%b2temp, max(this%lattice%control%llsp, this%lattice%control%lld))
      call g_safe_alloc%allocate('recursion.b2', this%b2, (/max(this%lattice%control%llsp, this%lattice%control%lld), nb, this%lattice%nrec, 3/))
      allocate (this%izero(0:this%lattice%kk), this%izeroll(0:this%lattice%kk, this%lattice%control%lld + 1), this%idum(0:this%lattice%kk))
      call g_safe_alloc%report_allocate('recursion.izero', this%izero)
      call g_safe_alloc%report_allocate('recursion.izeroll', this%izeroll)
      call g_safe_alloc%report_allocate('recursion.idum', this%idum)
      allocate (this%irlist(0:this%lattice%kk))
      call g_safe_alloc%report_allocate('recursion.irlist', this%irlist)

      call g_safe_alloc%allocate('recursion.psi', this%psi, (/nb, this%lattice%kk/))
      call g_safe_alloc%allocate('recursion.pmn', this%pmn, (/nb, this%lattice%kk/))
      call g_safe_alloc%allocate('recursion.psi1', this%psi1, (/nb, nb, this%lattice%kk/))
      call g_safe_alloc%allocate('recursion.psi2', this%psi2, (/nb, nb, this%lattice%kk/))
      call g_safe_alloc%allocate('recursion.psi0', this%psi0, (/nb, nb, this%lattice%kk/))
      call g_safe_alloc%allocate('recursion.v', this%v, (/nb, this%lattice%kk/))
      if ((this%lattice%njij == 0) .and. (this%lattice%njijk == 0)) then
         call g_safe_alloc%allocate('recursion.a_b', this%a_b, (/nb, nb, this%control%lld, this%lattice%nrec/))
         call g_safe_alloc%allocate('recursion.b2_b', this%b2_b, (/nb, nb, this%control%lld,, this%lattice%nrec/))
         call g_safe_alloc%allocate('recursion.mu_n', this%mu_n, (/nb, nb, (2*this%lattice%control%lld) + 2, this%lattice%nrec/))
         call g_safe_alloc%allocate('recursion.mu_ng', this%mu_ng, (/nb, nb, (2*this%lattice%control%lld) + 2, this%lattice%nrec/))
      else
         call g_safe_alloc%allocate('recursion.a_b', this%a_b, (/nb, nb, this%control%lld, this%lattice%njij*4/))
         call g_safe_alloc%allocate('recursion.b2_b', this%b2_b, (/nb, nb, this%control%lld, this%lattice%njij*4/))
         call g_safe_alloc%allocate('recursion.mu_n', this%mu_n, (/nb, nb, (2*this%lattice%control%lld) + 2, this%lattice%njij*4/))
         call g_safe_alloc%allocate('recursion.mu_ng', this%mu_ng, (/nb, nb, (2*this%lattice%control%lld) + 2, this%lattice%njij*4/))
      end if
      call g_safe_alloc%allocate('recursion.psi_b', this%psi_b, (/nb, nb, this%lattice%kk/))
      call g_safe_alloc%allocate('recursion.hpsi', this%hpsi, (/nb, nb, this%lattice%kk/))
      call g_safe_alloc%allocate('recursion.hohpsi', this%hohpsi, (/nb, nb, this%lattice%kk/))
      call g_safe_alloc%allocate('recursion.enupsi', this%enupsi, (/nb, nb, this%lattice%kk/))
      call g_safe_alloc%allocate('recursion.socpsi', this%socpsi, (/nb, nb, this%lattice%kk/))
      call g_safe_alloc%allocate('recursion.atemp_b', this%atemp_b, (/nb, nb, this%control%lld/))
      call g_safe_alloc%allocate('recursion.b2temp_b', this%b2temp_b, (/nb, nb, this%control%lld/))
      call g_safe_alloc%allocate('recursion.pmn_b', this%pmn_b, (/nb, nb, this%lattice%kk/))
      call g_safe_alloc%allocate('recursion.cheb_mom_temp', this%cheb_mom_temp, (/nb, nb/))
#else
      allocate (this%a(max(this%lattice%control%llsp, this%lattice%control%lld),&
                          &nb, this%lattice%nrec, 3))
      allocate (this%atemp(max(this%lattice%control%llsp, this%lattice%control%lld)))
      allocate (this%b2temp(max(this%lattice%control%llsp, this%lattice%control%lld)))
      allocate (this%b2(max(this%lattice%control%llsp, this%lattice%control%lld),&
                          &nb, this%lattice%nrec, 3))
      allocate (this%izero(0:this%lattice%kk), this%idum(0:this%lattice%kk), this%izeroll(0:this%lattice%kk, this%lattice%control%lld + 1))
      allocate (this%irlist(0:this%lattice%kk))

      allocate (this%psi(nb, this%lattice%kk), this%pmn(nb, this%lattice%kk))
      allocate (this%psi1(nb, nb, this%lattice%kk), this%psi2(nb, nb, this%lattice%kk), this%psi0(nb, nb, this%lattice%kk))
      allocate (this%v(nb, this%lattice%kk))
      if (this%lattice%njij == 0) then
         allocate (this%a_b(nb, nb, this%control%lld, this%lattice%nrec))
         allocate (this%b2_b(nb, nb, this%control%lld, this%lattice%nrec))
         allocate (this%mu_n(nb, nb, (2*this%lattice%control%lld) + 2, this%lattice%nrec), &
                   this%mu_ng(nb, nb, (2*this%lattice%control%lld) + 2, this%lattice%nrec))
      else
         allocate (this%a_b(nb, nb, this%control%lld, this%lattice%njij*4))
         allocate (this%b2_b(nb, nb, this%control%lld, this%lattice%njij*4))
         allocate (this%mu_n(nb, nb, (2*this%lattice%control%lld) + 2, this%lattice%njij*4), &
                   this%mu_ng(nb, nb, (2*this%lattice%control%lld) + 2, this%lattice%njij*4))
      end if
      allocate (this%psi_b(nb, nb, this%lattice%kk))
      allocate (this%hpsi(nb, nb, this%lattice%kk))
      allocate (this%hohpsi(nb, nb, this%lattice%kk))
      allocate (this%enupsi(nb, nb, this%lattice%kk))
      allocate (this%socpsi(nb, nb, this%lattice%kk))
      allocate (this%atemp_b(nb, nb, this%control%lld))
      allocate (this%b2temp_b(nb, nb, this%control%lld))
      allocate (this%pmn_b(nb, nb, this%lattice%kk))
      allocate (this%cheb_mom_temp(nb, nb))
#endif
      this%v(:, :) = 0.0d0
      this%psi(:, :) = 0.0d0
      this%psi1(:, :, :) = 0.0d0
      this%psi2(:, :, :) = 0.0d0
      this%pmn(:, :) = 0.0d0
      this%mu_n(:, :, :, :) = 0.0d0
      this%mu_ng(:, :, :, :) = 0.0d0
      this%izero(:) = 0
      this%idum(:) = 0
      this%izeroll(:, :) = 0
      this%irlist(:) = 0
      this%a(:, :, :, :) = 0.0d0
      this%b2(:, :, :, :) = 0.0d0
      this%atemp(:) = 0.0d0
      this%b2temp(:) = 0.0d0
      this%a_b(:, :, :, :) = 0.0d0
      this%b2_b(:, :, :, :) = 0.0d0
      this%psi_b(:, :, :) = 0.0d0
      this%hpsi(:, :, :) = 0.0d0
      this%hohpsi(:, :, :) = 0.0d0
      this%enupsi(:, :, :) = 0.0d0
      this%socpsi(:, :, :) = 0.0d0
      this%atemp_b(:, :, :) = 0.0d0
      this%b2temp_b(:, :, :) = 0.0d0
      this%pmn_b(:, :, :) = 0.0d0
      this%cheb_mom_temp(:, :) = 0.0d0
      if (present(full)) then
         if (full) then
            if (associated(this%hamiltonian)) call this%hamiltonian%restore_to_default()
            if (associated(this%lattice)) call this%lattice%restore_to_default()
            if (associated(this%en)) call this%en%restore_to_default()
         end if
      end if

   end subroutine restore_to_default

end submodule recursion_core
