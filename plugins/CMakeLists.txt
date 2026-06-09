set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/plugins)
file(MAKE_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

include(CheckLanguage)

set(ONEMKL_ROOT "" CACHE PATH "Root directory of an Intel oneMKL installation")

add_library(rslmto_cpu_sparse_plugin SHARED
    cpu_sparse_plugin.c
)

target_include_directories(rslmto_cpu_sparse_plugin PRIVATE
    ${CMAKE_SOURCE_DIR}/source
)

set_target_properties(rslmto_cpu_sparse_plugin PROPERTIES
    OUTPUT_NAME rslmto_cpu_sparse_plugin
    PREFIX "lib"
)

add_library(rslmto_openmp_sparse_plugin SHARED
    openmp_sparse_plugin.c
)

target_include_directories(rslmto_openmp_sparse_plugin PRIVATE
    ${CMAKE_SOURCE_DIR}/source
)

if(ENABLE_OPENMP)
  find_package(OpenMP QUIET COMPONENTS C)
  if(OpenMP_C_FOUND)
    target_link_libraries(rslmto_openmp_sparse_plugin PRIVATE OpenMP::OpenMP_C)
  endif()
endif()

if(BUILD_CUDA_PLUGIN)
  check_language(CUDA)
  if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    find_package(CUDAToolkit QUIET)

    if(CUDAToolkit_FOUND)
      add_library(rslmto_cuda_sparse_plugin SHARED
        cuda_sparse_plugin.cu
      )

      target_link_libraries(rslmto_cuda_sparse_plugin PRIVATE
        CUDA::cudart
        CUDA::cusparse
        CUDA::cublas
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
      )

      target_include_directories(rslmto_cuda_sparse_plugin PRIVATE
        ${CMAKE_SOURCE_DIR}/source
      )

      set_target_properties(rslmto_cuda_sparse_plugin PROPERTIES
        OUTPUT_NAME rslmto_cuda_sparse_plugin
        PREFIX "lib"
        CUDA_STANDARD 14
        CUDA_STANDARD_REQUIRED ON
      )

      target_compile_definitions(rslmto_cuda_sparse_plugin PRIVATE
        RSLMTO_CUDA_PLUGIN=1
      )

      message(STATUS "Building CUDA sparse plugin with compiler: ${CMAKE_CUDA_COMPILER}")
    else()
      message(STATUS "CUDA sparse plugin disabled: CUDAToolkit not found.")
    endif()
  else()
    message(STATUS "CUDA sparse plugin disabled: no CUDA compiler detected.")
  endif()
endif()

set_target_properties(rslmto_openmp_sparse_plugin PROPERTIES
    OUTPUT_NAME rslmto_openmp_sparse_plugin
    PREFIX "lib"
)

if(BUILD_ONEMKL_PLUGIN)
  set(_onemkl_hints)
  if(ONEMKL_ROOT)
    list(APPEND _onemkl_hints "${ONEMKL_ROOT}")
  endif()
  if(DEFINED ENV{MKLROOT})
    list(APPEND _onemkl_hints "$ENV{MKLROOT}")
  endif()
  list(APPEND _onemkl_hints
    /opt/intel/oneapi/mkl/latest
    /opt/intel/oneapi
    /opt/intel/oneapi/compiler/latest
    /opt/intel/oneapi/mkl
    /usr/local/oneapi/mkl/latest
  )

  find_path(ONEMKL_INCLUDE_DIR
    NAMES mkl_spblas.h
    HINTS ${_onemkl_hints}
    PATH_SUFFIXES include
  )

  find_library(ONEMKL_RT_LIBRARY
    NAMES mkl_rt
    HINTS ${_onemkl_hints}
    PATH_SUFFIXES lib lib64 lib/intel64 lib/intel64_lin lib/intel64_win
  )

  if(ONEMKL_INCLUDE_DIR AND ONEMKL_RT_LIBRARY)
    add_library(rslmto_onemkl_sparse_plugin SHARED
      onemkl_sparse_plugin.c
    )

    target_include_directories(rslmto_onemkl_sparse_plugin PRIVATE
      ${CMAKE_SOURCE_DIR}/source
      ${ONEMKL_INCLUDE_DIR}
    )

    target_link_libraries(rslmto_onemkl_sparse_plugin PRIVATE
      ${ONEMKL_RT_LIBRARY}
      Threads::Threads
      ${CMAKE_DL_LIBS}
    )

    set_target_properties(rslmto_onemkl_sparse_plugin PROPERTIES
      OUTPUT_NAME rslmto_onemkl_sparse_plugin
      PREFIX "lib"
    )

    target_compile_definitions(rslmto_onemkl_sparse_plugin PRIVATE
      RSLMTO_ONEMKL_PLUGIN=1
    )

    message(STATUS "Building oneMKL sparse plugin with include dir: ${ONEMKL_INCLUDE_DIR}")
    message(STATUS "Building oneMKL sparse plugin with runtime: ${ONEMKL_RT_LIBRARY}")
  else()
    message(STATUS "oneMKL sparse plugin disabled: could not find mkl_spblas.h and mkl_rt. Set ONEMKL_ROOT or MKLROOT to enable it.")
  endif()
endif()
