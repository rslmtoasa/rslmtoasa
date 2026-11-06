# if(COMMAND cmake_policy)
#   cmake_policy(SET CMP0003 NEW)
#   cmake_policy(SET CMP0076 NEW)
# endif()
# 
# cmake_policy(SET CMP0011 NEW)
# cmake_policy(SET CMP0012 NEW)
# cmake_policy(SET CMP0048 NEW)
# cmake_policy(SET CMP0054 NEW)


cmake_minimum_required(VERSION 3.19)

set(PROJECT_NAME rslmto)
project(${PROJECT_NAME})

# Find all the .in files
# file(GLOB inFiles RELATIVE "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/src/*.in" "${PROJECT_SOURCE_DIR}/test/*.in")
# foreach(infileName ${inFiles})
#   string(REGEX REPLACE ".in\$" "" outfileName "${infileName}")
#   configure_file(${PROJECT_SOURCE_DIR}/${infileName} ${PROJECT_SOURCE_DIR}/${outfileName})
# endforeach(infileName)

set(SRC ${CMAKE_SOURCE_DIR}/source)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules/")

# Enable Fortran language
enable_language(Fortran)
message(STATUS "Using Fortran compiler ${CMAKE_Fortran_COMPILER_ID} Fortran ${CMAKE_Fortran_COMPILER_VERSION}")

# Fortran module and output directories
set(fortran_module_directory "${CMAKE_BINARY_DIR}/modules" CACHE PATH "Directory for fortran module files")
file(MAKE_DIRECTORY ${fortran_module_directory})
set(CMAKE_Fortran_MODULE_DIRECTORY ${fortran_module_directory})
set(fortran_output_directory "${CMAKE_BINARY_DIR}/bin" CACHE PATH "Directory for executable file")
file(MAKE_DIRECTORY ${fortran_output_directory})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${fortran_output_directory})

##############################################################################################
# Set default CMake options
##############################################################################################
option(ENABLE_OPENMP "Use OpenMP for parallelization" ON)
option(ENABLE_MPI "Use MPI for parallelization" OFF)
option(ENABLE_LIBXC "Use libXC for exchange-correlation functionals" OFF)
option(ENABLE_CUDA_PLUGIN "Build the optional CUDA recursion plugin backend" OFF)
option(ENABLE_CUDA_RECIPROCAL "Build the optional CUDA reciprocal backend skeleton" OFF)
option(ENABLE_MKL_KERNELS "Build optional MKL Chebyshev kernels (cgemm_batch and inspector-executor sparse)" OFF)
option(ENABLE_FLUSH "Flush printouts" OFF)
option(COLOR "Color printouts" ON)
option(RUN_REG_TESTS "Run regression tests" OFF)
option(ENABLE_SPGLIB "Enable spglib for symmetry operations and k-path generation" ON)
if (APPLE)
  option(ENABLE_FUSED_RECIPROCAL "Fuse reciprocal H(k) assembly and eigensolution" OFF)
else()
  option(ENABLE_FUSED_RECIPROCAL "Fuse reciprocal H(k) assembly and eigensolution" ON)
endif()
option(RUN_EXAMPLE_TESTS "Run fast example-based SCF tests" OFF)
option(RUN_UNIT_TESTS "Build and run standalone Fortran unit tests" OFF)
set(EXAMPLE_REF_ABS_TOL "1e-4" CACHE STRING "Absolute tolerance for example reference comparisons")
set(EXAMPLE_REF_REL_TOL "1e-5" CACHE STRING "Relative tolerance for example reference comparisons")
mark_as_advanced(ENABLE_FLUSH fortran_module_directory fortran_output_directory CMAKE_OSX_*)
if (APPLE)
  mark_as_advanced(CMAKE_OSX_ARCHITECTURES CMAKE_OSX_DEPLOYMENT_TARGET CMAKE_OSX_SYSROOT)
endif(APPLE)

if (ENABLE_MKL_KERNELS OR ENABLE_CUDA_PLUGIN OR ENABLE_CUDA_RECIPROCAL)
  enable_language(CXX)
  message(STATUS "Using C++ compiler ${CMAKE_CXX_COMPILER_ID} C++ ${CMAKE_CXX_COMPILER_VERSION}")
endif()

if (ENABLE_CUDA_PLUGIN OR ENABLE_CUDA_RECIPROCAL)
  # Let a CUDA build target the GPU present while it is configured.  This
  # avoids carrying a stale architecture from a different machine (for
  # example, sm_75 code falling back to PTX JIT on an Ada GPU).  Do not
  # replace an explicit -DCMAKE_CUDA_ARCHITECTURES=... or CUDAARCHS choice:
  # those are required for cross-compiles, headless build hosts, and
  # distributable multi-architecture binaries.
  #
  # CMake added the `native` architecture selector in 3.24.  Retain the
  # project's 3.19 minimum for CPU and older-CMake CUDA builds; on those
  # versions the CUDA compiler's default remains in effect and we explain how
  # to select an architecture explicitly.
  if (NOT DEFINED CMAKE_CUDA_ARCHITECTURES AND "$ENV{CUDAARCHS}" STREQUAL "")
    if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.24)
      set(CMAKE_CUDA_ARCHITECTURES native CACHE STRING
          "CUDA architectures to build for; use native to detect the configuration host GPU")
      message(STATUS "CUDA architectures: native (detected from the configuration host GPU)")
    else()
      message(STATUS
        "CMake ${CMAKE_VERSION} cannot auto-detect CUDA architectures; "
        "set -DCMAKE_CUDA_ARCHITECTURES=<NN> when configuring CUDA builds")
    endif()
  endif()

  enable_language(CUDA)
  message(STATUS "Using CUDA compiler ${CMAKE_CUDA_COMPILER_ID} CUDA ${CMAKE_CUDA_COMPILER_VERSION}")
  message(STATUS "CUDA architectures: ${CMAKE_CUDA_ARCHITECTURES}")
  set(_cuda_toolkit_components cudart)
  if (ENABLE_CUDA_PLUGIN)
    list(APPEND _cuda_toolkit_components cusparse cufft cublas)
  endif()
  if (ENABLE_CUDA_RECIPROCAL)
    list(APPEND _cuda_toolkit_components cusolver)
  endif()
  find_package(CUDAToolkit REQUIRED COMPONENTS ${_cuda_toolkit_components})
  message(STATUS "CUDA toolkit root: ${CUDAToolkit_LIBRARY_ROOT}")
endif()

##############################################################################################
# Determine and set compiler-dependent flags
##############################################################################################
set(LOCAL_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")
include(${LOCAL_MODULE_PATH}/SetFortranFlags.cmake)
set(CMAKE_VERBOSE_MAKEFILE off)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

if (COLOR)
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DCOLOR")
endif(COLOR)

if (ENABLE_FLUSH)
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DEnable_flush")
endif(ENABLE_FLUSH)

if (ENABLE_OPENMP)
  find_package(OpenMP)
  message(STATUS "OpenMP libraries: ${OpenMP_Fortran_LIBRARIES}")
  if (OpenMP_Fortran_FOUND)
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DOpenMP_Fortran_FOUND")
  endif()
endif(ENABLE_OPENMP)

if(ENABLE_MPI)
    # The user may select an MPI wrapper before configuration, but never mutate
    # CMAKE_Fortran_COMPILER after Fortran has been enabled. MPI::MPI_Fortran
    # carries the correct compiler options, module paths, and link libraries.
    find_package(MPI REQUIRED COMPONENTS Fortran)
endif(ENABLE_MPI)

# The source uses default-width BLAS/LAPACK integer arguments (and the MKL kernels use c_int).
# Avoid silently changing all Fortran integers (including MPI arguments) to eight bytes.
if(BLA_VENDOR MATCHES "ilp" OR BLA_SIZEOF_INTEGER STREQUAL "8")
    message(FATAL_ERROR "Use LP64 BLAS/LAPACK, e.g. -DBLA_VENDOR=Intel10_64lp")
endif()
if(BLA_VENDOR MATCHES "^Intel" AND DEFINED ENV{MKL_IL})
    if(NOT "$ENV{MKL_IL}" STREQUAL "lp64")
        message(FATAL_ERROR "This build requires MKL_IL=lp64 (or unset MKL_IL)")
    endif()
endif()

if(ENABLE_SPGLIB)
    # Try to find spglib using pkg-config first
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
        pkg_check_modules(SPGLIB QUIET spglib)
    endif()
    
    # If pkg-config didn't find it, try manual search
    if(NOT SPGLIB_FOUND)
        find_path(SPGLIB_INCLUDE_DIR
            NAMES spglib.h
            HINTS /usr/local/include /usr/include /opt/local/include
            PATH_SUFFIXES spglib
        )
        
        find_library(SPGLIB_LIBRARY
            NAMES spglib symspg
            HINTS /usr/local/lib /usr/lib /opt/local/lib
        )
        
        if(SPGLIB_INCLUDE_DIR AND SPGLIB_LIBRARY)
            set(SPGLIB_FOUND TRUE)
            set(SPGLIB_LIBRARIES ${SPGLIB_LIBRARY})
            set(SPGLIB_INCLUDE_DIRS ${SPGLIB_INCLUDE_DIR})
        endif()
    endif()
    
    if(SPGLIB_FOUND)
    # If pkg-config reported specific library file variables (pkgcfg_lib_*), prefer those
    # They typically contain full paths like /usr/local/lib/libsymspg.dylib
    if(DEFINED pkgcfg_lib_SPGLIB_symspg)
      set(SPGLIB_LIBRARIES ${pkgcfg_lib_SPGLIB_symspg})
    endif()

    # If SPGLIB_LDFLAGS contains -L flags, add them to the link directories so the linker finds the library
    if(DEFINED SPGLIB_LDFLAGS)
      string(REGEX MATCHALL "-L[^; ]+" _spglib_ldflags_list "${SPGLIB_LDFLAGS}")
      foreach(_lf ${_spglib_ldflags_list})
        string(SUBSTRING ${_lf} 2 -1 _libdir)
        if(EXISTS ${_libdir})
          list(APPEND SPGLIB_LIBRARY_DIRS ${_libdir})
        endif()
      endforeach()
    endif()

    message(STATUS "Found spglib: ${SPGLIB_LIBRARIES}")
    add_definitions(-DUSE_SPGLIB)
    include_directories(${SPGLIB_INCLUDE_DIRS})
    if(SPGLIB_LIBRARY_DIRS)
      link_directories(${SPGLIB_LIBRARY_DIRS})
    endif()
    else()
        message(WARNING "spglib requested but not found. Install spglib or set ENABLE_SPGLIB=OFF")
        set(ENABLE_SPGLIB OFF)
    endif()
endif(ENABLE_SPGLIB)

if(ENABLE_LIBXC)
    # Find libXC using pkg-config (preferred method for portability)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(LIBXC REQUIRED IMPORTED_TARGET libxc>=5.0)
    pkg_check_modules(LIBXCF03 REQUIRED IMPORTED_TARGET libxcf03>=5.0)
    
    if(LIBXC_FOUND AND LIBXCF03_FOUND)
        message(STATUS "libXC found via pkg-config")
        message(STATUS "libXC version: ${LIBXC_VERSION}")
        message(STATUS "libXC includes: ${LIBXC_INCLUDE_DIRS}")
        message(STATUS "libXC libraries: ${LIBXC_LIBRARIES}")
        message(STATUS "libXC library paths: ${LIBXC_LIBRARY_DIRS}")
        message(STATUS "libXCF03 libraries: ${LIBXCF03_LIBRARIES}")
        
        # Set up compiler flags
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DHAVE_LIBXC")
        
        # Create aliases for cleaner linking
        add_library(LibXC::LibXC ALIAS PkgConfig::LIBXC)
        add_library(LibXC::LibXCF03 ALIAS PkgConfig::LIBXCF03)
    else()
        message(FATAL_ERROR "libXC or libXCF03 not found via pkg-config. Please install libXC or set PKG_CONFIG_PATH.")
    endif()
endif(ENABLE_LIBXC)

find_package(BLAS REQUIRED)
if(BLAS_FOUND)
  message(STATUS "BLAS libraries: ${BLAS_LIBRARIES}")
endif()

find_package(LAPACK REQUIRED)
if(LAPACK_FOUND)
  message(STATUS "LAPACK libraries: ${LAPACK_LIBRARIES}")
endif()

if(BLAS95_FOUND)
  message(STATUS "BLAS Fortran 95 interface Found")
  message(STATUS "BLAS95 libraries: ${BLAS95_LIBRARIES}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DBLAS95_FOUND")
endif()
if(LAPACK95_FOUND)
  message(STATUS "LAPACK Fortran 95 interface Found")
  message(STATUS "LAPACK95 libraries: ${LAPACK95_LIBRARIES}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DLAPACK95_FOUND")
endif()

# FindBLAS/FindLAPACK supply the selected MKL link dependencies for both Intel (ifort) and
# IntelLLVM (ifx). Do not override them with hard-coded archive paths, global -i8 flags, or
# CMAKE_EXE_LINKER_FLAGS. The MKL include paths of the optional kernels are set in source/CMakeLists.txt.
if (BLAS_FOUND AND BLAS95_FOUND AND BLA_VENDOR MATCHES "^Intel")
  set(HAVE_INTEL_MKL TRUE)
endif()

if (ENABLE_MKL_KERNELS)
  if (NOT HAVE_INTEL_MKL AND NOT DEFINED ENV{MKLROOT})
    message(FATAL_ERROR "ENABLE_MKL_KERNELS requires Intel oneMKL. Set BLA_VENDOR to an Intel MKL vendor or define MKLROOT.")
  endif()
  message(STATUS "MKL Sparse Inspector-Executor Chebyshev kernels: ON")
  message(STATUS "MKL cgemm_batch Chebyshev kernels: ON")
endif()

##############################################################################################
# Determine and set compiler-dependent flags
##############################################################################################
# Add main executable
set(PROJECT_BINARY ${PROJECT_NAME}.x)
add_executable(${PROJECT_BINARY} ${SRC}/main.f90)

# Add the directory of the sources
add_subdirectory(${SRC})

if(ENABLE_MPI AND MPI_Fortran_FOUND)
    target_compile_definitions(${PROJECT_NAME} PUBLIC USE_MPI)
    target_link_libraries(${PROJECT_NAME} PUBLIC MPI::MPI_Fortran)
endif()

if(NOT ENABLE_FUSED_RECIPROCAL)
    target_compile_definitions(${PROJECT_NAME} PUBLIC RSLMTO_DISABLE_FUSED_RECIPROCAL)
endif()

# Linking with Libraries
set(COMMON_LIBRARIES ${PROJECT_NAME} ${BLAS95_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK95_LIBRARIES} Threads::Threads ${CMAKE_DL_LIBS})

# Add spglib if enabled
if(ENABLE_SPGLIB AND SPGLIB_FOUND)
  list(APPEND COMMON_LIBRARIES ${SPGLIB_LIBRARIES})
endif()

# Add libXC libraries if enabled
if(ENABLE_LIBXC AND LIBXC_FOUND)
    list(APPEND COMMON_LIBRARIES PkgConfig::LIBXC PkgConfig::LIBXCF03)
    # Add include directories to the library target for Fortran modules
    target_include_directories(${PROJECT_NAME} PRIVATE ${LIBXC_INCLUDE_DIRS})
endif()

if(ENABLE_OPENMP AND ENABLE_MPI AND OpenMP_Fortran_FOUND AND MPI_Fortran_FOUND)
    target_link_libraries(${PROJECT_BINARY} ${COMMON_LIBRARIES} MPI::MPI_Fortran OpenMP::OpenMP_Fortran)
elseif(ENABLE_OPENMP AND OpenMP_Fortran_FOUND)
    target_link_libraries(${PROJECT_BINARY} ${COMMON_LIBRARIES} OpenMP::OpenMP_Fortran)
elseif(ENABLE_MPI AND MPI_Fortran_FOUND)
    target_link_libraries(${PROJECT_BINARY} ${COMMON_LIBRARIES} MPI::MPI_Fortran)
else()
    target_link_libraries(${PROJECT_BINARY} ${COMMON_LIBRARIES})
endif()

##############################################################################################
# Print current flags
##############################################################################################
message(STATUS "--------Summary of compiler options---------")
message(STATUS "| Compiler ID : ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "| Build type  : ${CMAKE_BUILD_TYPE}")
message(STATUS "| OpenMP      : ${ENABLE_OPENMP}")
message(STATUS "| MPI         : ${ENABLE_MPI}")
message(STATUS "| libXC       : ${ENABLE_LIBXC}")
message(STATUS "| CUDA plugin : ${ENABLE_CUDA_PLUGIN}")
message(STATUS "| CUDA reciprocal : ${ENABLE_CUDA_RECIPROCAL}")
message(STATUS "--------------------------------------------")

##############################################################################################
# Tell where to install the executable
##############################################################################################
IF(WIN32)
   SET(CMAKE_INSTALL_PREFIX "C:/Program Files")
ELSE()
   SET(CMAKE_INSTALL_PREFIX /usr/local)
ENDIF(WIN32)
# Convert to native path (optional, but useful if you do further path manipulations)
file(TO_NATIVE_PATH "${CMAKE_INSTALL_PREFIX}" CMAKE_INSTALL_PREFIX_NATIVE)
message(STATUS "Install prefix (native): ${CMAKE_INSTALL_PREFIX_NATIVE}")

SET_TARGET_PROPERTIES(${PROJECT_BINARY} PROPERTIES LINKER_LANGUAGE "Fortran")
INSTALL(TARGETS ${PROJECT_BINARY} RUNTIME DESTINATION bin)

##############################################################################################
# Tests
##############################################################################################
include(${CMAKE_SOURCE_DIR}/tests/register_tests.cmake)

## ADD Cpack
##############################################################################################
# Include CPack
INCLUDE(${LOCAL_MODULE_PATH}/git-watcher.cmake)
# The response driver writes the same build identifier into production output
# provenance.  The global definition above is evaluated after the library
# target is created, so apply it explicitly to the library as well.
target_compile_definitions(${PROJECT_NAME} PRIVATE ${gitlabel})
message(STATUS "Git version info : ${gitversion}.")
# Set package properties
set(CPACK_PACKAGE_NAME "rslmto")
set(CPACK_PACKAGE_VERSION "0.9.0")
#set(CPACK_PACKAGE_VERSION "${gitversion}")
message(STATUS "CPACK_PACKAGE_VERSION: ${CPACK_PACKAGE_VERSION}")
set(CPACK_PACKAGE_CONTACT "ramon@rslmto.org")
set(CPACK_PACKAGE_VENDOR "The RS-LMTO-ASA collective")
set(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}/cpack/icon.jpg")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/README.md")
set(CPACK_RESOURCE_FILE_WELCOME "${CMAKE_SOURCE_DIR}/cpack/WELCOME.txt")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/cpack/README.txt")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "RS-LMTO-ASA is a real-space linearized muffin-tin orbital method for electronic structure calculations.")
set(CPACK_PACKAGE_HOMEPAGE_URL "https://github.com/rslmto/rslmto")
#set(CPACK_GENERATOR "STGZ;DEB;RPM")
# DEB specifics
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "${CPACK_PACKAGE_CONTACT}")
set (CPACK_DEBIAN_PACKAGE_RELEASE ${GIT_HASH})
set (CPACK_DEBIAN_PACKAGE_GENERATE_SHLIBS ON)
set (CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
set (CPACK_PACKAGE_NAME "${CMAKE_PROJECT_NAME}")
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "${CMAKE_PROJECT_NAME} - ${CMAKE_PROJECT_DESCRIPTION}")
set (CPACK_DEBIAN_PACKAGE_DEBUG ON)

if(WIN32)
   include(InstallRequiredSystemLibraries)

    # Include NSIS for creating Windows installers
    set(CPACK_GENERATOR "NSIS64")

    # Include the custom NSIS script
    # set(CPACK_NSIS_INCLUDE "${CMAKE_SOURCE_DIR}/cpack/CustomInstallNSIS.nsh")

    # Additional NSIS settings if needed
    set(CPACK_NSIS_MODIFY_PATH ON)
    set(CPACK_NSIS_ENABLE_UNINSTALL_BEFORE_INSTALL ON)
    set(CPACK_NSIS_EXTRA_INSTALL_COMMANDS "Call AddToPath")

    # WIX installer
    # set(CPACK_WIX_UPGRADE_GUID "E25F58C1-EC45-422A-9700-7254DB19547B")
    # WIX installer
    set(CPACK_WIX_LICENSE_RTF "${CMAKE_SOURCE_DIR}/cpack/LICENSE.rtf")
    set(CPACK_WIX_UI_BANNER "${CMAKE_SOURCE_DIR}/cpack/banner.jpg")
    set(CPACK_WIX_UI_DIALOG "${CMAKE_SOURCE_DIR}/cpack/dialog.jpg")
    set(CPACK_WIX_UPGRADE_GUID "02cf4cf7-fa5b-43de-a8ff-b39098d35db3")
    set(CPACK_WIX_PATCH_FILE "${CMAKE_SOURCE_DIR}/cpack/wixpath.xml")
    # set(CPACK_WIX_EXTRA_SOURCES "${CMAKE_SOURCE_DIR}/PathEnvironment.wxs")
   set(CPACK_WIX_VERSION "4")
   endif()

include(CPack)
