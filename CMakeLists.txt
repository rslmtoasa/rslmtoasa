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
  enable_language(CUDA)
  message(STATUS "Using CUDA compiler ${CMAKE_CUDA_COMPILER_ID} CUDA ${CMAKE_CUDA_COMPILER_VERSION}")
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

# MKL machine dependent settings
if (BLAS_FOUND AND BLAS95_FOUND)
  if (BLA_VENDOR MATCHES "^Intel")
    set(HAVE_INTEL_MKL TRUE)
    if(NOT DEFINED ENV{MKL_IL})
      set(ENV{MKL_IL} ilp64)
    elseif(NOT $ENV{MKL_IL} STREQUAL lp64 OR NOT $ENV{MKL_IL} STREQUAL ilp64)
      message(FATAL_ERROR "Invalid option for MKL_IL")
    endif()

    if(NOT EXISTS "$ENV{MKLROOT}/include")
      set(ENV{MKLROOT} "$ENV{MKLROOT}/mkl")
    endif()

    if ($ENV{MKL_IL} STREQUAL lp64)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64 -I$ENV{MKLROOT}/include -I$ENV{MKLROOT}/include/intel64/lp64")
    endif()

    if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      if ($ENV{MKL_IL} STREQUAL ilp64)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64 -fdefault-integer-8 -I$ENV{MKLROOT}/include -I$ENV{MKLROOT}/include/intel64/ilp64")
      endif()
    endif()

    if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
      if ($ENV{MKL_IL} STREQUAL ilp64)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64 -i8 -I$ENV{MKLROOT}/include -I$ENV{MKLROOT}/include/intel64/ilp64")
      endif()
      if (APPLE)
        set(MKL_LIBRARY_DIRS "$ENV{MKLROOT}/lib")
        find_library(MKL_CORE mkl_core.a HINTS ${MKL_LIBRARY_DIRS})
        find_library(MKL_COMPILER mkl_intel_$ENV{MKL_IL}.a HINTS ${MKL_LIBRARY_DIRS})
        find_library(MKL_THREADING mkl_intel_thread.a HINTS ${MKL_LIBRARY_DIRS})
      endif()
      if (UNIX AND NOT APPLE)
        set(MKL_LIBRARY_DIRS "$ENV{MKLROOT}/lib/intel64")
        find_library(MKL_CORE libmkl_core.a HINTS ${MKL_LIBRARY_DIRS})
        find_library(MKL_COMPILER libmkl_intel_$ENV{MKL_IL}.a HINTS ${MKL_LIBRARY_DIRS})
        find_library(MKL_THREADING libmkl_intel_thread.a HINTS ${MKL_LIBRARY_DIRS})
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-intel -Wl,--start-group ${MKL_CORE} ${MKL_COMPILER} ${MKL_THREADING} -Wl,--end-group")
      endif()
    endif()
  endif()
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
set(COMMON_LIBRARIES ${PROJECT_NAME} ${BLAS95_LIBRARIES} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES} Threads::Threads ${CMAKE_DL_LIBS})

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
include(CTest)
#CMAKE_DEPENDENT_OPTION(RUN_REGRESSION_TESTS "Run regression tests"  OFF "BUILD_TESTING" OFF)
###############################################################################################
## Add tests directory
###############################################################################################
set(TestBinary ${CMAKE_BINARY_DIR}/bin/${PROJECT_BINARY})
set(TestDIR ${CMAKE_SOURCE_DIR}/tests)

# Read optional semantic labels from one manifest case.  Missing or empty
# labels are valid so older cases do not need to be changed.  Keep duplicate
# removal at the registration site so explicit labels cannot duplicate a
# derived label such as "scf" or "backend".
function(read_manifest_case_labels MANIFEST_JSON CASE_INDEX OUT_VAR)
  set(_labels "")
  string(JSON _label_count ERROR_VARIABLE _label_error LENGTH
    "${MANIFEST_JSON}" "cases" ${CASE_INDEX} "labels")
  if (NOT _label_error)
    if (_label_count GREATER 0)
      math(EXPR _label_last "${_label_count} - 1")
      foreach(_label_index RANGE ${_label_last})
        string(JSON _label GET
          "${MANIFEST_JSON}" "cases" ${CASE_INDEX} "labels" ${_label_index})
        list(APPEND _labels "${_label}")
      endforeach()
    endif()
  endif()
  set(${OUT_VAR} "${_labels}" PARENT_SCOPE)
endfunction()

# Resolve the effective structure-constant backend for one manifest case.  The
# case-local value wins over the optional manifest default; absent values retain
# the production default (legacy) and are labelled accordingly.
function(read_manifest_case_strux_backend MANIFEST_JSON CASE_INDEX OUT_VAR)
  set(_backend "legacy")
  string(JSON _default_backend ERROR_VARIABLE _default_error GET
    "${MANIFEST_JSON}" "defaults" "namelists" "lattice" "strux_backend")
  if (NOT _default_error)
    set(_backend "${_default_backend}")
  endif()
  string(JSON _case_backend ERROR_VARIABLE _case_error GET
    "${MANIFEST_JSON}" "cases" ${CASE_INDEX} "namelists" "lattice" "strux_backend")
  if (NOT _case_error)
    set(_backend "${_case_backend}")
  endif()
  set(${OUT_VAR} "${_backend}" PARENT_SCOPE)
endfunction()

if (RUN_REG_TESTS)
set(Python3_FIND_VIRTUALENV FIRST)
find_package(Python3 REQUIRED COMPONENTS Interpreter)

set(ThisTestDIR ${TestDIR}/regression/bccFe_lanczos)
add_test(NAME Lanczos   COMMAND  ${ThisTestDIR}/oneliner.sh ${TestBinary}  WORKING_DIRECTORY ${ThisTestDIR})
set_tests_properties(Lanczos PROPERTIES LABELS "functional;regression;baseline;metallic;rs;lanczos")
set(ThisTestDIR ${TestDIR}/regression/bccFe_block)
add_test(NAME Block     COMMAND  ${ThisTestDIR}/oneliner.sh ${TestBinary}  WORKING_DIRECTORY ${ThisTestDIR})
set_tests_properties(Block PROPERTIES LABELS "functional;regression;baseline;metallic;rs;block")
set(ThisTestDIR ${TestDIR}/regression/bccFe_chebyshev)
add_test(NAME Chebyshev COMMAND  ${ThisTestDIR}/oneliner.sh ${TestBinary}  WORKING_DIRECTORY ${ThisTestDIR})
set_tests_properties(Chebyshev PROPERTIES LABELS "functional;regression;baseline;metallic;rs;chebyshev")

# The WP8/WP9 decks and runners remain under tests/regression as historical
# scientific evidence, but their active CTest registrations were deliberately
# removed. Re-enable them only through a separate validation decision.

set(REGRESSION_CASES_JSON "${TestDIR}/regression/cases.json")
file(READ "${REGRESSION_CASES_JSON}" _regression_json)
string(JSON _regression_n LENGTH "${_regression_json}" "cases")
math(EXPR _regression_last "${_regression_n} - 1")

foreach(_i RANGE ${_regression_last})
  string(JSON _regression_name GET "${_regression_json}" "cases" ${_i} "name")
  string(JSON _regression_timeout GET "${_regression_json}" "cases" ${_i} "timeout")
  read_manifest_case_labels("${_regression_json}" ${_i} _regression_case_labels)
  read_manifest_case_strux_backend("${_regression_json}" ${_i} _regression_strux_backend)
  set(_regression_labels "functional;regression;backend")
  list(APPEND _regression_labels ${_regression_case_labels})
  if (_regression_strux_backend STREQUAL "strux_lib")
    list(APPEND _regression_labels "strux")
  elseif (_regression_strux_backend STREQUAL "legacy")
    list(APPEND _regression_labels "legacy_strux")
  endif()
  list(REMOVE_DUPLICATES _regression_labels)
  add_test(NAME "Regression_${_regression_name}"
    COMMAND "${Python3_EXECUTABLE}" "${TestDIR}/regression/run_matrix.py"
      "--binary" "${TestBinary}"
      "--cases-json" "${REGRESSION_CASES_JSON}"
      "--case-name" "${_regression_name}"
      "--scratch-root" "${CMAKE_BINARY_DIR}/Testing/regression"
      "--references" "${TestDIR}/regression/references"
  )
  set_tests_properties("Regression_${_regression_name}" PROPERTIES
    LABELS "${_regression_labels}" TIMEOUT "${_regression_timeout}")
endforeach()

# B5.2 route-agnostic estimator triads (recursion vs lehmann vs dyson). Each case
# runs one observable three ways and checks the documented agreement envelopes
# (docs/dev/route_agnostic_estimators.md).
set(TRIAD_CASES_JSON "${TestDIR}/regression/triad_cases.json")
file(READ "${TRIAD_CASES_JSON}" _triad_json)
string(JSON _triad_n LENGTH "${_triad_json}" "cases")
math(EXPR _triad_last "${_triad_n} - 1")
foreach(_i RANGE ${_triad_last})
  string(JSON _triad_name GET "${_triad_json}" "cases" ${_i} "name")
  string(JSON _triad_timeout GET "${_triad_json}" "cases" ${_i} "timeout")
  add_test(NAME "Triad_${_triad_name}"
    COMMAND "${Python3_EXECUTABLE}" "${TestDIR}/regression/run_triad.py"
      "--binary" "${TestBinary}"
      "--cases-json" "${TRIAD_CASES_JSON}"
      "--case-name" "${_triad_name}"
      "--scratch-root" "${CMAKE_BINARY_DIR}/Testing/triad"
      "--references" "${TestDIR}/regression/references_triad"
  )
  set_tests_properties("Triad_${_triad_name}" PROPERTIES
    LABELS "functional;regression;triad" TIMEOUT "${_triad_timeout}")
endforeach()

# Small no-binary fixture: this checks report parsing and campaign policy only;
# it is tooling, not scientific regression evidence.
add_test(NAME TddftValidationCampaign
  COMMAND "${Python3_EXECUTABLE}" "${TestDIR}/regression/tddft_validation/test_validation.py")
set_tests_properties(TddftValidationCampaign PROPERTIES LABELS "tooling;tddft" TIMEOUT 30)

  add_test(NAME Val04LdaUPhysics
    COMMAND "${Python3_EXECUTABLE}" "${TestDIR}/validation/val04_lda_u.py"
      "--binary" "${TestBinary}"
      "--scratch-root" "${CMAKE_BINARY_DIR}/Testing/validation/val04_lda_u")
  set_tests_properties(Val04LdaUPhysics PROPERTIES
    LABELS "validation;lda_u;magnetic;kspace" TIMEOUT 300)

  add_test(NAME Val05GreenConvergence
    COMMAND "${Python3_EXECUTABLE}" "${TestDIR}/validation/val05_green_convergence.py"
      "--binary" "${TestBinary}"
      "--scratch-root" "${CMAKE_BINARY_DIR}/Testing/validation/val05_green_convergence")
  set_tests_properties(Val05GreenConvergence PROPERTIES
    LABELS "validation;functional;green;lehmann;dyson;kspace;rs;convergence" TIMEOUT 1800)

  add_test(NAME Val07ExchangeFormulations
    COMMAND "${Python3_EXECUTABLE}" "${TestDIR}/validation/val07_exchange_formulations.py"
      "--binary" "${TestBinary}"
      "--scratch-root" "${CMAKE_BINARY_DIR}/Testing/validation/val07_exchange_formulations")
  set_tests_properties(Val07ExchangeFormulations PROPERTIES
    LABELS "validation;functional;exchange;magnetic;tensor;multisite" TIMEOUT 600)

  add_test(NAME Val08DampingInertia
    COMMAND "${Python3_EXECUTABLE}" "${TestDIR}/validation/val08_damping_inertia.py"
      "--binary" "${TestBinary}"
      "--scratch-root" "${CMAKE_BINARY_DIR}/Testing/validation/val08_damping_inertia")
  set_tests_properties(Val08DampingInertia PROPERTIES
    LABELS "validation;functional;exchange;magnetic;damping;inertia;convergence" TIMEOUT 1800)

  add_test(NAME Val09KuboBastinTransport
    COMMAND "${Python3_EXECUTABLE}" "${TestDIR}/validation/val09_kubo_bastin_transport.py"
      "--binary" "${TestBinary}"
      "--scratch-root" "${CMAKE_BINARY_DIR}/Testing/validation/val09_kubo_bastin_transport")
  set_tests_properties(Val09KuboBastinTransport PROPERTIES
    LABELS "validation;functional;conductivity;transport;kpm;spin_transport;orbital_transport;convergence;tensor" TIMEOUT 1800)

  add_test(NAME Val12LmtoFieldsTorques
    COMMAND "${Python3_EXECUTABLE}" "${TestDIR}/validation/val12_lmto_fields_torques.py"
      "--binary" "${TestBinary}"
      "--scratch-root" "${CMAKE_BINARY_DIR}/Testing/validation/val12_lmto_fields_torques")
  set_tests_properties(Val12LmtoFieldsTorques PROPERTIES
    LABELS "validation;functional;magnetic;spin_dynamics;lmto;torque;constraints;soc" TIMEOUT 900)

  add_test(NAME Val13AbInitioSpinDynamics
    COMMAND "${Python3_EXECUTABLE}" "${TestDIR}/validation/val13_abinitio_spin_dynamics.py"
      "--binary" "${TestBinary}"
      "--scratch-root" "${CMAKE_BINARY_DIR}/Testing/validation/val13_abinitio_spin_dynamics")
  set_tests_properties(Val13AbInitioSpinDynamics PROPERTIES
    LABELS "validation;functional;magnetic;spin_dynamics;abinitio;trajectory;convergence" TIMEOUT 1800)

  add_test(NAME Val15MultilayerVacuumElectrostatics
    COMMAND "${Python3_EXECUTABLE}" "${TestDIR}/validation/val15_multilayer_vacuum.py"
      "--binary" "${TestBinary}"
      "--scratch-root" "${CMAKE_BINARY_DIR}/Testing/validation/val15_multilayer_vacuum")
  set_tests_properties(Val15MultilayerVacuumElectrostatics PROPERTIES
    LABELS "validation;functional;electrostatics;vacuum;green;convergence" TIMEOUT 1800)

  add_test(NAME Val16GbtCommensurateSupercells
    COMMAND "${Python3_EXECUTABLE}" "${TestDIR}/validation/val16_gbt_supercell.py"
      "--binary" "${TestBinary}"
      "--scratch-root" "${CMAKE_BINARY_DIR}/Testing/validation/val16_gbt_supercell")
  set_tests_properties(Val16GbtCommensurateSupercells PROPERTIES
    LABELS "validation;functional;magnetic;supercell;commensurate;convergence"
    TIMEOUT 14400)

  add_test(NAME Val17GbtHarmonicGoldstone
    COMMAND "${Python3_EXECUTABLE}" "${TestDIR}/validation/val17_gbt_harmonic_goldstone.py"
      "--binary" "${TestBinary}"
      "--scratch-root" "${CMAKE_BINARY_DIR}/Testing/validation/val17_gbt_harmonic_goldstone")
  set_tests_properties(Val17GbtHarmonicGoldstone PROPERTIES
    LABELS "validation;functional;magnetic;gbt;kspace;goldstone;convergence"
    TIMEOUT 14400)
endif (RUN_REG_TESTS)

##############################################################################################
# Standalone Fortran unit tests (known-answer and focused production-path tests)
##############################################################################################
if (RUN_UNIT_TESTS)
  find_package(Python3 COMPONENTS Interpreter REQUIRED)
  set(_standalone_unit_targets "")

  # Each test is a small program linked against the project library. Link the
  # same conditional MPI/OpenMP libraries as the main binary.
  function(add_fortran_unit_test _name _source)
    add_executable(${_name} ${_source})
    if(ENABLE_OPENMP AND ENABLE_MPI AND OpenMP_Fortran_FOUND AND MPI_Fortran_FOUND)
      target_link_libraries(${_name} ${COMMON_LIBRARIES} MPI::MPI_Fortran OpenMP::OpenMP_Fortran)
    elseif(ENABLE_OPENMP AND OpenMP_Fortran_FOUND)
      target_link_libraries(${_name} ${COMMON_LIBRARIES} OpenMP::OpenMP_Fortran)
    elseif(ENABLE_MPI AND MPI_Fortran_FOUND)
      target_link_libraries(${_name} ${COMMON_LIBRARIES} MPI::MPI_Fortran)
    else()
      target_link_libraries(${_name} ${COMMON_LIBRARIES})
    endif()
    set_target_properties(${_name} PROPERTIES LINKER_LANGUAGE "Fortran")
    add_test(NAME ${_name} COMMAND ${_name})
    set_tests_properties(${_name} PROPERTIES LABELS "unit")
    list(APPEND _standalone_unit_targets ${_name})
    set(_standalone_unit_targets "${_standalone_unit_targets}" PARENT_SCOPE)
  endfunction()

  add_fortran_unit_test(UnitBasisDimensions ${TestDIR}/unit/test_basis_dimensions.f90)
  set_tests_properties(UnitBasisDimensions PROPERTIES LABELS "unit;basis")
  add_fortran_unit_test(UnitGreenLifecycle ${TestDIR}/unit/test_green_lifecycle.f90)
  set_tests_properties(UnitGreenLifecycle PROPERTIES LABELS "unit;green;lifecycle")
  if(ENABLE_LIBXC)
    add_fortran_unit_test(UnitLibxcXcBaseline ${TestDIR}/unit/test_libxc_xc_baseline.f90)
    set_tests_properties(UnitLibxcXcBaseline PROPERTIES LABELS "unit;libxc;xc")
  endif()
  add_fortran_unit_test(UnitLdaUHamiltonian ${TestDIR}/unit/test_lda_u_hamiltonian.f90)
  set_tests_properties(UnitLdaUHamiltonian PROPERTIES LABELS "unit;lda_u;hamiltonian")
  add_fortran_unit_test(UnitLehmannChain ${TestDIR}/unit/test_lehmann_chain.f90)
  add_fortran_unit_test(UnitDysonEquivalence ${TestDIR}/unit/test_dyson_equivalence.f90)
  add_fortran_unit_test(UnitGammaSupercell ${TestDIR}/unit/test_gamma_supercell.f90)
  add_fortran_unit_test(UnitBsfSumRule ${TestDIR}/unit/test_bsf_sumrule.f90)
  add_fortran_unit_test(UnitMomentKernel ${TestDIR}/unit/test_moment_kernel.f90)
  add_fortran_unit_test(UnitKpmTransport ${TestDIR}/unit/test_kpm_transport.f90)
  set_tests_properties(UnitKpmTransport PROPERTIES
    LABELS "unit;conductivity;transport;kpm;spin_transport;orbital_transport")
  add_fortran_unit_test(UnitMadl2dExhGuard ${TestDIR}/unit/test_madl2d_exh_guard.f90)
  add_fortran_unit_test(UnitMadelungConventions ${TestDIR}/unit/test_madelung_conventions.f90)
  add_fortran_unit_test(UnitRegionRegistry ${TestDIR}/unit/test_region_registry.f90)
  add_fortran_unit_test(UnitInterfaceAlignmentOracle ${TestDIR}/unit/test_interface_alignment_oracle.f90)
  add_fortran_unit_test(UnitVacuumLead ${TestDIR}/unit/test_vacuum_lead.f90)
  add_fortran_unit_test(UnitVacuumRegionWiring ${TestDIR}/unit/test_vacuum_region_wiring.f90)
  add_fortran_unit_test(UnitAlignmentSolver ${TestDIR}/unit/test_alignment_solver.f90)
  add_fortran_unit_test(UnitKspaceOccupations ${TestDIR}/unit/test_kspace_occupations.f90)
  add_fortran_unit_test(UnitParallelContext ${TestDIR}/unit/test_parallel_context.f90)
  add_fortran_unit_test(UnitKpointWorkset ${TestDIR}/unit/test_kpoint_workset.f90)
  set_tests_properties(UnitParallelContext PROPERTIES LABELS "unit;mpi;parallel_context")
  add_fortran_unit_test(UnitArbitraryKEigenpairs ${TestDIR}/unit/test_arbitrary_k_eigenpairs.f90)
  add_test(NAME Acc04ArbitraryKSource
    COMMAND ${Python3_EXECUTABLE} ${TestDIR}/unit/test_acc04_arbitrary_k_source.py)
  set_tests_properties(Acc04ArbitraryKSource PROPERTIES LABELS "unit;cuda;reciprocal;arbitrary-k;source-contract")
  add_fortran_unit_test(UnitMagneticRepresentation ${TestDIR}/unit/test_magnetic_representation.f90)
  set_tests_properties(UnitMagneticRepresentation PROPERTIES LABELS "unit;gbt;wp3")
  add_fortran_unit_test(UnitGbtStructure ${TestDIR}/unit/test_gbt_structure.f90)
  set_tests_properties(UnitGbtStructure PROPERTIES LABELS "unit;gbt;wp4")
  add_fortran_unit_test(UnitStructureConstantsBackends ${TestDIR}/unit/test_structure_constants_backends.f90)
  target_compile_definitions(UnitStructureConstantsBackends PRIVATE
    TEST_STRUCTURE_CONSTANT_DATABASE="${CMAKE_SOURCE_DIR}/tests/scf/cases/bulk/bccFe")
  file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/Testing/structure_constants")
  set_tests_properties(UnitStructureConstantsBackends PROPERTIES
    LABELS "unit;structure_constants;strux;legacy_strux"
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/Testing/structure_constants")
  add_fortran_unit_test(UnitLmtoMagneticTangents ${TestDIR}/unit/test_lmto_magnetic_tangents.f90)
  set_tests_properties(UnitLmtoMagneticTangents PROPERTIES LABELS "unit;tddft;ward;lmto;tangent")
  add_fortran_unit_test(UnitLmtoPairPotential ${TestDIR}/unit/test_lmto_pair_potential.f90)
  set_tests_properties(UnitLmtoPairPotential PROPERTIES LABELS "unit;tddft;ward;lmto;pair-potential")
  add_fortran_unit_test(UnitSpinDensity ${TestDIR}/unit/test_spin_density.f90)
  set_tests_properties(UnitSpinDensity PROPERTIES LABELS "unit;gbt;wp7;density")
  add_fortran_unit_test(UnitAbspinlibIntegrators ${TestDIR}/validation/test_abspinlib_integrators.f90)
  set_tests_properties(UnitAbspinlibIntegrators PROPERTIES
    LABELS "validation;sd;spin_integrators;physics;convergence")
  add_fortran_unit_test(UnitConstrainingField ${TestDIR}/unit/test_constraining_field.f90)
  set_tests_properties(UnitConstrainingField PROPERTIES LABELS "unit;constraints;physics")
  add_test(NAME UnitConstrainingFieldSource
    COMMAND ${Python3_EXECUTABLE} ${TestDIR}/unit/test_constraining_field_source.py)
  set_tests_properties(UnitConstrainingFieldSource PROPERTIES LABELS "unit;constraints;source-contract")
  add_test(NAME CudaPluginSurface
    COMMAND ${Python3_EXECUTABLE} ${TestDIR}/cuda/test_plugin_surface.py)
  set_tests_properties(CudaPluginSurface PROPERTIES LABELS "unit;cuda;gpu;coverage;source-contract")
  if(ENABLE_CUDA_RECIPROCAL)
    add_executable(ReciprocalCudaArbitraryK ${TestDIR}/unit/test_reciprocal_cuda_arbitrary_k.f90)
    target_link_libraries(ReciprocalCudaArbitraryK PRIVATE ${COMMON_LIBRARIES} rslmto_reciprocal_cuda)
    set_target_properties(ReciprocalCudaArbitraryK PROPERTIES LINKER_LANGUAGE "Fortran")
    add_test(NAME ReciprocalCudaArbitraryK COMMAND ReciprocalCudaArbitraryK)
    set_tests_properties(ReciprocalCudaArbitraryK PROPERTIES
      LABELS "unit;cuda;reciprocal;arbitrary-k;integration" SKIP_RETURN_CODE 77)
    add_executable(ReciprocalCudaLifecycle ${TestDIR}/unit/test_reciprocal_cuda_lifecycle.cpp)
    target_include_directories(ReciprocalCudaLifecycle PRIVATE ${SRC}/cuda)
    target_link_libraries(ReciprocalCudaLifecycle PRIVATE rslmto_reciprocal_cuda)
    add_test(NAME ReciprocalCudaLifecycle COMMAND ReciprocalCudaLifecycle)
    set_tests_properties(ReciprocalCudaLifecycle PROPERTIES
      LABELS "unit;cuda;reciprocal;lifecycle" SKIP_RETURN_CODE 77)
    add_executable(ReciprocalCudaEigensolver ${TestDIR}/unit/test_reciprocal_cuda_eigensolver.cpp)
    target_include_directories(ReciprocalCudaEigensolver PRIVATE ${SRC}/cuda)
    target_link_libraries(ReciprocalCudaEigensolver PRIVATE rslmto_reciprocal_cuda)
    add_test(NAME ReciprocalCudaEigensolver COMMAND ReciprocalCudaEigensolver)
    set_tests_properties(ReciprocalCudaEigensolver PROPERTIES
      LABELS "unit;cuda;reciprocal;eigensolver" SKIP_RETURN_CODE 77)
  endif()
  add_fortran_unit_test(UnitResponseConventions ${TestDIR}/unit/test_response_conventions.f90)
  set_tests_properties(UnitResponseConventions PROPERTIES LABELS "unit;tddft;conventions")
  add_fortran_unit_test(UnitTddftWardConventions ${TestDIR}/unit/test_tddft_ward_conventions.f90)
  set_tests_properties(UnitTddftWardConventions PROPERTIES LABELS "unit;tddft;ward;conventions")
  add_fortran_unit_test(UnitResponseVertices ${TestDIR}/unit/test_response_vertices.f90)
  set_tests_properties(UnitResponseVertices PROPERTIES LABELS "unit;tddft;vertices")
  add_fortran_unit_test(UnitTddftChiKS ${TestDIR}/unit/test_tddft_chi_ks.f90)
  set_tests_properties(UnitTddftChiKS PROPERTIES LABELS "unit;tddft;chi_ks;stoner")
  add_fortran_unit_test(UnitTddftTransitionWorkspace ${TestDIR}/unit/test_tddft_transition_workspace.f90)
  set_tests_properties(UnitTddftTransitionWorkspace PROPERTIES LABELS "unit;tddft;workspace")
  add_fortran_unit_test(UnitTddftDirectXi ${TestDIR}/unit/test_tddft_direct_xi.f90)
  set_tests_properties(UnitTddftDirectXi PROPERTIES LABELS "unit;tddft;xi;ward")
  add_fortran_unit_test(UnitTddftCpuProfile ${TestDIR}/unit/test_tddft_cpu_profile.f90)
  set_tests_properties(UnitTddftCpuProfile PROPERTIES LABELS "performance;tddft;profile")
  add_fortran_unit_test(UnitTddftGreenChiKS ${TestDIR}/unit/test_tddft_green_chi0.f90)
  set_tests_properties(UnitTddftGreenChiKS PROPERTIES LABELS "unit;tddft;chi_ks;green")
  add_fortran_unit_test(UnitTddftGoldstone ${TestDIR}/unit/test_tddft_goldstone.f90)
  set_tests_properties(UnitTddftGoldstone PROPERTIES LABELS "unit;tddft;goldstone;alsda")
  add_fortran_unit_test(UnitTddftDysonModes ${TestDIR}/unit/test_tddft_dyson_modes.f90)
  set_tests_properties(UnitTddftDysonModes PROPERTIES LABELS "unit;tddft;dyson;modes;damping")
  add_fortran_unit_test(UnitTddftLongitudinal ${TestDIR}/unit/test_tddft_longitudinal.f90)
  set_tests_properties(UnitTddftLongitudinal PROPERTIES LABELS "unit;tddft;longitudinal;relaxation")
  add_fortran_unit_test(UnitTddftFourComponent ${TestDIR}/unit/test_tddft_four_component.f90)
  set_tests_properties(UnitTddftFourComponent PROPERTIES LABELS "unit;tddft;four-component;soc;noncollinear")
  add_fortran_unit_test(UnitTddftConfig ${TestDIR}/unit/test_tddft_config.f90)
  set_tests_properties(UnitTddftConfig PROPERTIES LABELS "unit;tddft;input;dispatch")
  add_test(NAME TddftCrossMilestoneEquivalence
    COMMAND ${Python3_EXECUTABLE} ${TestDIR}/unit/test_tddft_cross_milestone_equivalence.py
      --binary-dir ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  set_tests_properties(TddftCrossMilestoneEquivalence PROPERTIES
    LABELS "unit;tddft;equivalence;TDDFT-03-10" TIMEOUT 120)
  add_test(NAME UnitTddftDispatch
    COMMAND ${Python3_EXECUTABLE} ${TestDIR}/unit/test_tddft_dispatch.py)
  set_tests_properties(UnitTddftDispatch PROPERTIES LABELS "tooling;tddft;dispatch")
  set_tests_properties(UnitKspaceOccupations PROPERTIES LABELS "unit;gbt;wp2;wp5")
  if(ENABLE_MPI AND MPI_Fortran_FOUND)
    # Homebrew Open MPI on macOS runners may report only one available slot,
    # while the parallel-context unit intentionally exercises four ranks.
    # Pass the launcher override in front of the test executable so it is
    # consumed by mpiexec rather than forwarded as an application argument.
    set(_rslmto_mpi_unit_preflags ${MPIEXEC_PREFLAGS})
    if(APPLE)
      list(APPEND _rslmto_mpi_unit_preflags --map-by :OVERSUBSCRIBE)
    endif()
    add_test(NAME UnitKspaceOccupations_mpi
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${_rslmto_mpi_unit_preflags}
              $<TARGET_FILE:UnitKspaceOccupations> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(UnitKspaceOccupations_mpi PROPERTIES
      LABELS "unit;gbt;wp2;mpi" PROCESSORS 2)
    foreach(_parallel_context_ranks IN ITEMS 1 2 4)
      add_test(NAME "UnitParallelContext_mpi_${_parallel_context_ranks}"
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${_parallel_context_ranks} ${_rslmto_mpi_unit_preflags}
                $<TARGET_FILE:UnitParallelContext> ${MPIEXEC_POSTFLAGS})
      set_tests_properties("UnitParallelContext_mpi_${_parallel_context_ranks}" PROPERTIES
        LABELS "unit;mpi;parallel_context" PROCESSORS ${_parallel_context_ranks})
    endforeach()
  endif()
  add_test(NAME UnitGbtOracles
    COMMAND ${Python3_EXECUTABLE} ${TestDIR}/unit/test_gbt_oracles.py)
  set_tests_properties(UnitGbtOracles PROPERTIES LABELS "unit;gbt")
  add_test(NAME UnitGbtWp5SourceContract
    COMMAND ${Python3_EXECUTABLE} ${TestDIR}/unit/test_gbt_wp5_source_contract.py)
  set_tests_properties(UnitGbtWp5SourceContract PROPERTIES LABELS "tooling;gbt;wp5")
  add_test(NAME UnitGbtWp10SourceContract
    COMMAND ${Python3_EXECUTABLE} ${TestDIR}/unit/test_gbt_wp10_source_contract.py)
  set_tests_properties(UnitGbtWp10SourceContract PROPERTIES LABELS "tooling;gbt;wp10")
  add_test(NAME UnitGbtWp6Hoh
    COMMAND ${Python3_EXECUTABLE} ${TestDIR}/unit/test_gbt_wp6_hoh.py)
  set_tests_properties(UnitGbtWp6Hoh PROPERTIES LABELS "unit;gbt;wp6;hoh;overlap")
  add_test(NAME UnitGbtWp6Ccor
    COMMAND ${Python3_EXECUTABLE} ${TestDIR}/unit/test_gbt_wp6_ccor.py)
  set_tests_properties(UnitGbtWp6Ccor PROPERTIES LABELS "unit;gbt;wp6;ccor")
  add_test(NAME UnitGbtWp6cVelocity
    COMMAND ${Python3_EXECUTABLE} ${TestDIR}/unit/test_gbt_wp6c_velocity.py)
  set_tests_properties(UnitGbtWp6cVelocity PROPERTIES LABELS "unit;gbt;wp6;velocity;torque")
  # WP6 integration: keeps tests/gbt_wp6_matrix/feature_matrix.json in step with
  # the guards actually present in source/, and owns the oracles that span the
  # WP6a/WP6b/WP6c slice boundaries (onsite fold, CCOR+HOH composition).
  add_test(NAME UnitGbtWp6Matrix
    COMMAND ${Python3_EXECUTABLE} ${TestDIR}/unit/test_gbt_wp6_matrix.py)
  set_tests_properties(UnitGbtWp6Matrix PROPERTIES LABELS "unit;gbt;wp6;matrix")
  add_test(NAME UnitManifestDefaults
    COMMAND ${Python3_EXECUTABLE} ${TestDIR}/unit/test_manifest_defaults.py)
  set_tests_properties(UnitManifestDefaults PROPERTIES LABELS "tooling;manifest")
  add_test(NAME BenchmarkHarnessSchema
    COMMAND ${Python3_EXECUTABLE} ${TestDIR}/benchmarks/test_benchmark_harness.py)
  set_tests_properties(BenchmarkHarnessSchema PROPERTIES LABELS "tooling;performance;benchmark")
  # This test also inspects the madl2d source directly, so that the corrected
  # overflow guard cannot silently regress (see the test's header).
  target_compile_definitions(UnitMadl2dExhGuard PRIVATE
    MADL2D_SOURCE="${CMAKE_CURRENT_SOURCE_DIR}/source/charge_madelung_2d.f90")
  add_custom_target(unit_tests DEPENDS ${_standalone_unit_targets})
endif (RUN_UNIT_TESTS)

if (RUN_EXAMPLE_TESTS)

# Prefer an active venv, then the system Python. Users can override via
# -DEXAMPLE_PYTHON_EXECUTABLE=/path/to/venv/python3 without reconfiguring.
set(Python3_FIND_VIRTUALENV FIRST)
find_package(Python3 REQUIRED COMPONENTS Interpreter)
set(EXAMPLE_PYTHON_EXECUTABLE "${Python3_EXECUTABLE}" CACHE FILEPATH
    "Python interpreter for example tests (must have f90nml installed)")

# Read cases.json and register CTest entries for one example test suite.
# Args: CASES_JSON  SUITE_LABEL  SCRATCH_ROOT  REFS_DIR
function(register_example_tests CASES_JSON SUITE SCRATCH_ROOT REFS_DIR)
  file(READ "${CASES_JSON}" _json)
  string(JSON _n LENGTH "${_json}" "cases")
  math(EXPR _last "${_n} - 1")

  foreach(_i RANGE ${_last})
    string(JSON _name    GET "${_json}" "cases" ${_i} "name")
    string(JSON _recur   GET "${_json}" "cases" ${_i} "namelists" "control" "recur")
    string(JSON _timeout GET "${_json}" "cases" ${_i} "timeout")

    # Cases requiring an optional CMake feature (e.g. ENABLE_MKL_KERNELS) are
    # skipped entirely when that feature is off, mirroring the regression
    # matrix's requires_cmake_option gate in run_matrix.py.
    string(JSON _req_opt ERROR_VARIABLE _req_opt_err GET "${_json}" "cases" ${_i} "requires_cmake_option")
    if (NOT _req_opt_err AND NOT ${_req_opt})
      continue()
    endif()

    # Cases marked "quick": true get an extra label so CI can select a small
    # representative subset for pull_request runs (see tests.yml) without
    # needing a second cases.json.
    string(JSON _quick ERROR_VARIABLE _quick_err GET "${_json}" "cases" ${_i} "quick")
    set(_quick_label "")
    if (NOT _quick_err AND ${_quick})
      set(_quick_label ";quick")
    endif()
    read_manifest_case_labels("${_json}" ${_i} _case_labels)
    read_manifest_case_strux_backend("${_json}" ${_i} _case_strux_backend)
    if (_case_strux_backend STREQUAL "strux_lib")
      list(APPEND _case_labels "strux")
    elseif (_case_strux_backend STREQUAL "legacy")
      list(APPEND _case_labels "legacy_strux")
    endif()

    set(_cmd
      "${EXAMPLE_PYTHON_EXECUTABLE}" "${TestDIR}/run_test.py"
      "--binary"       "${TestBinary}"
      "--cases-json"   "${CASES_JSON}"
      "--case-name"    "${_name}"
      "--scratch-root" "${SCRATCH_ROOT}"
      "--compare-ref"  "${REFS_DIR}"
      "--abs-tol"      "${EXAMPLE_REF_ABS_TOL}"
      "--rel-tol"      "${EXAMPLE_REF_REL_TOL}"
    )

    # Cases listing "launch_modes": ["serial", "mpi"] register a serial/MPI
    # pair that shares one reference, so a rank-decomposition bug shows up as
    # a mismatch against the same values the serial run is held to.
    string(JSON _n_modes ERROR_VARIABLE _modes_err LENGTH "${_json}" "cases" ${_i} "launch_modes")

    if (_modes_err)
      set(_default_cmd ${_cmd})
      if (ENABLE_MPI)
        list(APPEND _default_cmd "--mpi-enabled")
      endif()
      add_test(NAME "${_name}" COMMAND ${_default_cmd})
      set(_labels "functional;example;${SUITE};${_recur}${_quick_label}")
      list(APPEND _labels ${_case_labels})
      list(REMOVE_DUPLICATES _labels)
      set_tests_properties("${_name}" PROPERTIES
        LABELS "${_labels}" TIMEOUT "${_timeout}")
    else()
      string(JSON _mpi_procs GET "${_json}" "cases" ${_i} "mpi_procs")

      set(_serial_cmd ${_cmd} "--force-serial")
      add_test(NAME "${_name}" COMMAND ${_serial_cmd})
      set(_labels "functional;example;${SUITE};${_recur};launch_serial${_quick_label}")
      list(APPEND _labels ${_case_labels})
      list(REMOVE_DUPLICATES _labels)
      set_tests_properties("${_name}" PROPERTIES
        LABELS "${_labels}" TIMEOUT "${_timeout}")

      if (ENABLE_MPI)
        string(JSON _n_mpi_ranks ERROR_VARIABLE _mpi_ranks_err LENGTH "${_json}" "cases" ${_i} "mpi_rank_counts")
        if (_mpi_ranks_err)
          set(_n_mpi_ranks 1)
          set(_mpi_rank_0 "${_mpi_procs}")
        endif()
        math(EXPR _last_mpi_rank "${_n_mpi_ranks} - 1")
        foreach(_mpi_rank_index RANGE ${_last_mpi_rank})
          if (NOT _mpi_ranks_err)
            string(JSON _mpi_rank GET "${_json}" "cases" ${_i} "mpi_rank_counts" ${_mpi_rank_index})
          else()
            set(_mpi_rank "${_mpi_rank_0}")
          endif()
          if (_mpi_rank EQUAL _mpi_procs)
            set(_mpi_suffix "_mpi")
          else()
            set(_mpi_suffix "_mpi_${_mpi_rank}")
          endif()
          set(_mpi_test_name "${_name}${_mpi_suffix}")
          set(_mpi_cmd ${_cmd} "--mpi-enabled" "--mpi-procs" "${_mpi_rank}" "--scratch-name" "${_mpi_test_name}")
          add_test(NAME "${_mpi_test_name}" COMMAND ${_mpi_cmd})
          set(_labels "functional;example;${SUITE};${_recur};launch_mpi")
          list(APPEND _labels ${_case_labels})
          list(REMOVE_DUPLICATES _labels)
          set_tests_properties("${_mpi_test_name}" PROPERTIES
            LABELS "${_labels}" TIMEOUT "${_timeout}"
            PROCESSORS "${_mpi_rank}")
        endforeach()
      endif()
    endif()
  endforeach()
endfunction()

register_example_tests(
  "${TestDIR}/scf/cases.json"     scf
  "${CMAKE_BINARY_DIR}/Testing/scf"
  "${TestDIR}/scf/references"
)

register_example_tests(
  "${TestDIR}/postproc/cases.json"     postproc
  "${CMAKE_BINARY_DIR}/Testing/postproc"
  "${TestDIR}/postproc/references"
)

add_test(
  NAME Example_si_chebyshev_kspace_dos_equivalence
  COMMAND
    "${EXAMPLE_PYTHON_EXECUTABLE}"
    "${TestDIR}/scf/test_si_dos_equivalence.py"
    --binary "${TestBinary}"
    --case-root "${TestDIR}/scf/cases/bulk/diamondSi"
    --scratch-root "${CMAKE_BINARY_DIR}/Testing/scf/si_dos_equivalence"
)
set_tests_properties(Example_si_chebyshev_kspace_dos_equivalence PROPERTIES
  LABELS "functional;example;scf;equivalence;rs;kspace;chebyshev;basis_sp;nonmagnetic;strux;quick"
  TIMEOUT 300
)

if (ENABLE_LIBXC)
  add_test(
    NAME Val03LibxcScfEquivalence
    COMMAND
      "${EXAMPLE_PYTHON_EXECUTABLE}"
      "${TestDIR}/validation/val03_libxc_scf.py"
      --binary "${TestBinary}"
      --scratch-root "${CMAKE_BINARY_DIR}/Testing/validation/val03_libxc_scf"
      --kspace-only
  )
  set_tests_properties(Val03LibxcScfEquivalence PROPERTIES
    LABELS "validation;functional;scf;libxc;equivalence;nonmagnetic;magnetic;strux"
    TIMEOUT 300
  )
endif()

endif (RUN_EXAMPLE_TESTS)
 
##############################################################################################
## ADD Cpack
##############################################################################################
# Include CPack
INCLUDE(${LOCAL_MODULE_PATH}/git-watcher.cmake)
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
