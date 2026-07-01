# Add library
set(cmake_archive_output_directory lib)
message(STATUS "The project name is: ${PROJECT_NAME}")
add_library(${PROJECT_NAME} STATIC
    array.f90
    basis.f90
    bands.f90
    calculation.f90
    charge.f90
    chebyshev_fast.f90
    haydock_fast.f90
    conductivity.f90
    control.f90
    density_of_states.f90
    element.f90
    energy.f90
    exchange.f90
    face.F90
    globals.f90
    green.f90
    #hamiltonian_base_mod.f90
    #hamiltonian_cpu_mod.f90
    hamiltonian.f90
    lattice.f90
    # lists.f90
    logger.f90
    math.f90
    mix.f90
    mpi.f90
    namelist_generator.f90
    os.f90
    potential.f90
    precision.f90
    reciprocal.f90
    recursion.f90
    cuda/recursion_gpu_mod.f90
    sparse.f90
    report.f90
    rsrec_cuda_plugin.f90
    rsk_cuda_plugin.f90
    safe_alloc.f90
    self.f90
    bounds.f90
    spglib_interface.f90
    spin_dynamics.f90
    string.f90
    symbolic_atom.f90
    symmetry.f90
    timer.f90
    xc.f90
)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    # ifx 202x can ICE in charge.f90 around legacy FORMAT handling in release
    # builds. Keep the workaround local to this translation unit.
    set_source_files_properties(charge.f90 PROPERTIES COMPILE_OPTIONS "-O0")
endif()

if(ENABLE_MKL_BATCH OR ENABLE_MKL_SPARSE)
    if(DEFINED ENV{MKLROOT})
        target_include_directories(${PROJECT_NAME} PRIVATE "$ENV{MKLROOT}/include")
        if(EXISTS "$ENV{MKLROOT}/include/intel64/lp64")
            target_include_directories(${PROJECT_NAME} PRIVATE "$ENV{MKLROOT}/include/intel64/lp64")
        endif()
        if(EXISTS "$ENV{MKLROOT}/include/intel64/ilp64")
            target_include_directories(${PROJECT_NAME} PRIVATE "$ENV{MKLROOT}/include/intel64/ilp64")
        endif()
    endif()
endif()

if(ENABLE_MKL_BATCH)
    target_sources(${PROJECT_NAME} PRIVATE chebyshev_mkl_batch.cpp)
    target_compile_definitions(${PROJECT_NAME} PRIVATE USE_MKL_BATCH)
endif()

if(ENABLE_MKL_SPARSE)
    target_compile_definitions(${PROJECT_NAME} PRIVATE USE_MKL_SPARSE)
    if(DEFINED ENV{MKLROOT} AND EXISTS "$ENV{MKLROOT}/include/mkl_spblas.f90")
        target_sources(${PROJECT_NAME} PRIVATE "$ENV{MKLROOT}/include/mkl_spblas.f90")
    endif()
endif()

if(NOT ENABLE_CUDA_PLUGIN)
    # cuda/recursion_gpu_mod.f90 is compiled unconditionally and binds the
    # rsrec_* C interface. The real kernels (cuda/rsrec_gpu.cu) are only built
    # with the CUDA plugin, so on plain CPU builds we link a no-op stub that
    # makes the symbols resolve and fails cleanly if the gpu_plugin path is
    # ever requested at runtime. No compute kernels are pulled in.
    target_sources(${PROJECT_NAME} PRIVATE cuda/rsrec_stub.c cuda/rsk_stub.c)
endif()

if(ENABLE_CUDA_PLUGIN)
    # rsrec_cuda is a backend dispatcher that provides rsrec_cuda_* interface
    # for backward compatibility with rsrec_cuda_plugin.f90; internally dispatches to
    # rsrec_cpu.cpp (CPU) or rsrec_gpu.cu (GPU) based on backend selection.
    # Both rsrec_gpu.cu (CUDA) and rsrec.cpp (wrapper) are compiled into this target.
    add_library(rsrec_cuda STATIC
        cuda/rsrec_cuda.cpp
        cuda/rsrec_gpu.cu
        cuda/rsk_cuda.cpp
        cuda/rsk_gpu.cu
    )
    target_compile_features(rsrec_cuda PRIVATE cxx_std_11)
    # Link against rsrec_cpp and CUDA runtime libraries
    target_link_libraries(rsrec_cuda CUDA::cudart CUDA::cusparse CUDA::cufft CUDA::cublas)
    if(ENABLE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(rsrec_cuda OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(${PROJECT_NAME} rsrec_cuda)
    target_compile_definitions(${PROJECT_NAME} PRIVATE USE_CUDA_PLUGIN)
endif()

add_subdirectory(./include_codes/abspinlib)
target_link_libraries(${PROJECT_NAME} abspinlib)
add_subdirectory(./include_codes/strux_lib)
target_link_libraries(${PROJECT_NAME} strux_lib)
