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
option(ENABLE_FLUSH "Flush printouts" OFF)
option(COLOR "Color printouts" ON)
option(RUN_REG_TESTS "Run regression tests" OFF)
option(RUN_EXAMPLE_TESTS "Run fast example-based SCF tests" OFF)
set(EXAMPLE_REF_ABS_TOL "1e-6" CACHE STRING "Absolute tolerance for example reference comparisons")
set(EXAMPLE_REF_REL_TOL "1e-6" CACHE STRING "Relative tolerance for example reference comparisons")
mark_as_advanced(ENABLE_FLUSH fortran_module_directory fortran_output_directory CMAKE_OSX_*)
if (APPLE)
  mark_as_advanced(CMAKE_OSX_ARCHITECTURES CMAKE_OSX_DEPLOYMENT_TARGET CMAKE_OSX_SYSROOT)
endif(APPLE)

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

# Keep the compiler chosen at configure time. For Intel oneAPI, use ifx
# with -DMPI_Fortran_COMPILER=mpiifx, or select mpiifx as the compiler.
# Do not set I_MPI_F90 to a wrapper: that can make it invoke itself.
if(ENABLE_MPI)
    find_package(MPI REQUIRED COMPONENTS Fortran)
    add_compile_definitions(USE_MPI)
endif()

# The source uses default-width BLAS/LAPACK integer arguments. Avoid silently
# changing all Fortran integers (including MPI arguments) to eight bytes.
if(BLA_VENDOR MATCHES "ilp" OR BLA_SIZEOF_INTEGER STREQUAL "8")
    message(FATAL_ERROR "Use LP64 BLAS/LAPACK, e.g. -DBLA_VENDOR=Intel10_64lp")
endif()
if(BLA_VENDOR MATCHES "^Intel" AND DEFINED ENV{MKL_IL})
    if(NOT "$ENV{MKL_IL}" STREQUAL "lp64")
        message(FATAL_ERROR "This build requires MKL_IL=lp64 (or unset MKL_IL)")
    endif()
endif()

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

# FindBLAS/FindLAPACK supply the selected MKL link dependencies for both
# Intel (ifort) and IntelLLVM (ifx). Do not override them with hard-coded
# archive paths, global -i8 flags, or CMAKE_EXE_LINKER_FLAGS.

##############################################################################################
# Determine and set compiler-dependent flags
##############################################################################################
# Add main executable
set(PROJECT_BINARY ${PROJECT_NAME}.x)
add_executable(${PROJECT_BINARY} ${SRC}/main.f90)

# Add the directory of the sources
add_subdirectory(${SRC})

# Apply MPI usage requirements to the library that compiles mpi.f90.
# The plain signature matches source/CMakeLists.txt and propagates these
# requirements to the main executable and the existing test executables.
if(ENABLE_MPI)
    target_link_libraries(${PROJECT_NAME} MPI::MPI_Fortran)
endif()

# Linking with Libraries
if(ENABLE_OPENMP AND OpenMP_Fortran_FOUND)
    target_link_libraries(${PROJECT_BINARY} ${PROJECT_NAME} ${BLAS95_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK95_LIBRARIES} OpenMP::OpenMP_Fortran Threads::Threads ${CMAKE_DL_LIBS})
else()
    target_link_libraries(${PROJECT_BINARY} ${PROJECT_NAME} ${BLAS95_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK95_LIBRARIES} Threads::Threads ${CMAKE_DL_LIBS})
endif()

##############################################################################################
# Print current flags
##############################################################################################
message(STATUS "--------Summary of compiler options---------")
message(STATUS "| Compiler ID : ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "| OpenMP      : ${ENABLE_OPENMP}")
message(STATUS "| MPI         : ${ENABLE_MPI}")
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
# Dense-reference fixture runner: link the same production modules and libraries.
# Built explicitly by tests/regression/test_hoh_transport.py, not by default.
if(BUILD_TESTING)
  add_executable(transport_operator_test EXCLUDE_FROM_ALL tests/regression/transport_driver.f90)
  target_include_directories(transport_operator_test PRIVATE ${fortran_module_directory})
  get_target_property(transport_test_libraries ${PROJECT_BINARY} LINK_LIBRARIES)
  target_link_libraries(transport_operator_test PRIVATE ${transport_test_libraries})
  add_executable(response_units_test EXCLUDE_FROM_ALL tests/regression/response_units_driver.f90)
  target_include_directories(response_units_test PRIVATE ${fortran_module_directory})
  target_link_libraries(response_units_test PRIVATE ${transport_test_libraries})
endif()

#CMAKE_DEPENDENT_OPTION(RUN_REGRESSION_TESTS "Run regression tests"  OFF "BUILD_TESTING" OFF)
###############################################################################################
## Add tests directory
###############################################################################################
set(TestBinary ${CMAKE_BINARY_DIR}/bin/${PROJECT_BINARY})
set(TestDIR ${CMAKE_SOURCE_DIR}/tests)

if (RUN_REG_TESTS)
set(ThisTestDIR ${TestDIR}/regression/bccFe_lanczos)
add_test(NAME Lanczos   COMMAND  ${ThisTestDIR}/oneliner.sh ${TestBinary}  WORKING_DIRECTORY ${ThisTestDIR})
set(ThisTestDIR ${TestDIR}/regression/bccFe_block)
add_test(NAME Block     COMMAND  ${ThisTestDIR}/oneliner.sh ${TestBinary}  WORKING_DIRECTORY ${ThisTestDIR})
set(ThisTestDIR ${TestDIR}/regression/bccFe_chebyshev)
add_test(NAME Chebyshev COMMAND  ${ThisTestDIR}/oneliner.sh ${TestBinary}  WORKING_DIRECTORY ${ThisTestDIR})
endif (RUN_REG_TESTS)

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

    set(_cmd
      "${EXAMPLE_PYTHON_EXECUTABLE}" "${TestDIR}/run_test.py"
      "--binary"       "${TestBinary}"
      "--cases-json"   "${CASES_JSON}"
      "--case-name"    "${_name}"
      "--scratch-root" "${SCRATCH_ROOT}"
    )

    list(APPEND _cmd
      "--compare-ref" "${REFS_DIR}"
      "--abs-tol"     "${EXAMPLE_REF_ABS_TOL}"
      "--rel-tol"     "${EXAMPLE_REF_REL_TOL}"
    )

    add_test(NAME "${_name}" COMMAND ${_cmd})
    set_tests_properties("${_name}" PROPERTIES
      LABELS "example;${SUITE};${_recur}" TIMEOUT "${_timeout}")
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