# Add library
set(cmake_archive_output_directory lib)
message(STATUS "The project name is: ${PROJECT_NAME}")
add_library(${PROJECT_NAME} STATIC
    array.f90
    basis.f90
    bands.f90
    calculation.f90
    calculation_preprocessing.f90
    calculation_reciprocal.f90
    charge.f90
    charge_interface.f90
    charge_madelung_2d.f90
    chebyshev_fast.f90
    haydock_fast.f90
    bsf_kernel.f90
    conductivity.f90
    control.f90
    density_of_states.f90
    dyson_kernel.f90
    electrostatics_multipole.f90
    element.f90
    energy.f90
    exchange.f90
    exchange_q.f90
    exchange_dynamics.f90
    face.F90
    frozen_magnon.f90
    globals.f90
    green.f90
    green_lifecycle.f90
    green_chebyshev.f90
    green_block.f90
    green_lanczos.f90
    magnetic_representation.f90
    gbt_structure.f90
    lmto_magnetic_tangent.f90
    lr_kl_hessian.f90
    lr_kl_contour.f90
    lmto_radial_augmentation.f90
    pauli_ground_state_projection.f90
    spin_density.f90
    hamiltonian.f90
    hamiltonian_build.f90
    hamiltonian_ccor.f90
    hamiltonian_paoflow_io.f90
    hamiltonian_hubbard.f90
    lattice.f90
    lehmann_kernel.f90
    moment_kernel.f90
    lattice_lifecycle.f90
    lattice_cluster.f90
    lattice_strux.f90
    lattice_print.f90
    # lists.f90
    logger.f90
    math.f90
    mix.f90
    mpi.f90
    namelist_generator.f90
    os.f90
    potential.f90
    radial_ground_state.f90
    precision.f90
    kpm_profile.f90
    scf_benchmark_profile.f90
    kpoint_workset.f90
    reciprocal.f90
    reciprocal_backend.f90
    reciprocal_lifecycle.f90
    reciprocal_fourier.f90
    reciprocal_bands.f90
    reciprocal_export.f90
    reciprocal_occupations.f90
    reciprocal_dos.f90
    reciprocal_projection.f90
    reciprocal_green.f90
    reciprocal_bsf.f90
    reciprocal_moments.f90
    reciprocal_spin_density.f90
    region_registry.f90
    recursion.f90
    recursion_core.f90
    recursion_haydock.f90
    recursion_chebyshev.f90
    recursion_transport.f90
    response_angular_basis.f90
    response_basis_mapping.f90
    lr_response_space.f90
    lr_pauli_transition_vertex.f90
    lr_lmto_product_response_basis.f90
    lr_projected_site_spin.f90
    lr_compact_static_interaction.f90
    lr_gf_endpoint_augmentation.f90
    lr_ks_susceptibility.f90
    lr_gf_susceptibility.f90
    lr_product_gf_susceptibility.f90
    lr_projected_reciprocal_chi0.f90
    lr_projected_interacting_response.f90
    lr_projected_juelich_interaction.f90
    lr_dresp06a_bridge.f90
    lr_ward_mode_analysis.f90
    lr_kl_static_bridge.f90
    lr_lmto_turek_gf.f90
    lr_lmto_turek_contour.f90
    lr_rs_gf_susceptibility.f90
    lr_dresp07_radial_oracle.f90
    lr_exact_ks_ward.f90
    lr_dresp08_native_mapping.f90
    lr_dresp08_bridge.f90
    lr_dresp09_field_insertion.f90
    lr_dresp07_bridge.f90
    tddft_native_rsgf_provider.f90
    lr_alsda_kernel.f90
    lr_goldstone_sumrule.f90
    lr_goldstone_correction.f90
    tddft_production_driver.f90
    tddft_dyson.f90
    sparse.f90
    report.f90
    rsrec_cuda_plugin.f90
    safe_alloc.f90
    sigma_provider.f90
    self.f90
    self_reciprocal.f90
    bounds.f90
    spglib_interface.f90
    spin_dynamics.f90
    string.f90
    symbolic_atom.f90
    symmetry.f90
    timer.f90
    vacuum_lead.f90
    xc_radial.f90
    xc.f90
)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    # ifx 202x can ICE in the legacy FORMAT handling of the 2D Madelung
    # implementation. Keep -O0 local to that isolated translation unit.
    set_source_files_properties(charge_madelung_2d.f90 PROPERTIES COMPILE_OPTIONS "-O0")
endif()

if(ENABLE_MKL_KERNELS)
    add_library(rslmto_mkl_kernels INTERFACE)
    target_compile_definitions(rslmto_mkl_kernels INTERFACE USE_MKL_BATCH USE_MKL_SPARSE)
    if(DEFINED ENV{MKLROOT})
        target_include_directories(rslmto_mkl_kernels INTERFACE "$ENV{MKLROOT}/include")
        if(EXISTS "$ENV{MKLROOT}/include/intel64/lp64")
            target_include_directories(rslmto_mkl_kernels INTERFACE "$ENV{MKLROOT}/include/intel64/lp64")
        endif()
        if(EXISTS "$ENV{MKLROOT}/include/intel64/ilp64")
            target_include_directories(rslmto_mkl_kernels INTERFACE "$ENV{MKLROOT}/include/intel64/ilp64")
        endif()
    endif()
    target_link_libraries(${PROJECT_NAME} PRIVATE rslmto_mkl_kernels)
    if(DEFINED ENV{MKLROOT} AND EXISTS "$ENV{MKLROOT}/include/mkl_spblas.f90")
        target_sources(${PROJECT_NAME} PRIVATE "$ENV{MKLROOT}/include/mkl_spblas.f90")
    endif()
endif()

if(ENABLE_CUDA_PLUGIN)
    add_library(rslmto_cuda_plugin_options INTERFACE)
    target_compile_features(rslmto_cuda_plugin_options INTERFACE cxx_std_11)
    target_link_libraries(rslmto_cuda_plugin_options INTERFACE
        CUDA::cudart CUDA::cusparse CUDA::cufft CUDA::cublas)
    if(ENABLE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(rslmto_cuda_plugin_options INTERFACE OpenMP::OpenMP_CXX)
    endif()

    # rsrec_cuda is a backend dispatcher that provides rsrec_cuda_* interface
    # for backward compatibility with rsrec_cuda_plugin.f90; internally dispatches to
    # rsrec_cpu.cpp (CPU) or rsrec_gpu.cu (GPU) based on backend selection.
    # Both rsrec_gpu.cu (CUDA) and rsrec.cpp (wrapper) are compiled into this target.
    add_library(rsrec_cuda STATIC
        cuda/rsrec_cuda.cpp
        cuda/rsrec_gpu.cu
    )
    target_link_libraries(rsrec_cuda PUBLIC rslmto_cuda_plugin_options)
    target_link_libraries(${PROJECT_NAME} PRIVATE rsrec_cuda)
    target_compile_definitions(${PROJECT_NAME} PRIVATE USE_CUDA_PLUGIN)
endif()

if(ENABLE_CUDA_RECIPROCAL)
    add_library(rslmto_reciprocal_cuda STATIC
        cuda/reciprocal_cuda.cpp
        cuda/reciprocal_lehmann.cu
    )
    target_compile_features(rslmto_reciprocal_cuda PUBLIC cxx_std_11)
    target_link_libraries(rslmto_reciprocal_cuda PUBLIC CUDA::cudart CUDA::cusolver)
    target_compile_definitions(${PROJECT_NAME} PRIVATE USE_CUDA_RECIPROCAL)
    target_link_libraries(${PROJECT_NAME} PRIVATE rslmto_reciprocal_cuda)
endif()

add_subdirectory(./include_codes/abspinlib)
target_link_libraries(${PROJECT_NAME} PRIVATE abspinlib)
add_subdirectory(./include_codes/strux_lib)
target_link_libraries(${PROJECT_NAME} PRIVATE strux_lib)
