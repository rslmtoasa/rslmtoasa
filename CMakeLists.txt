# if(COMMAND cmake_policy)
#   cmake_policy(SET CMP0003 NEW)
#   cmake_policy(SET CMP0076 NEW)
# endif()
# 
# cmake_policy(SET CMP0011 NEW)
# cmake_policy(SET CMP0012 NEW)
# cmake_policy(SET CMP0048 NEW)
# cmake_policy(SET CMP0054 NEW)


cmake_minimum_required(VERSION 3.5)

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
set(fortran_module_directory "${CMAKE_BINARY_DIR}/modules" CACHE PATH "directory for fortran module files")
file(MAKE_DIRECTORY ${fortran_module_directory})
set(CMAKE_Fortran_MODULE_DIRECTORY ${fortran_module_directory})
set(fortran_output_directory "${CMAKE_BINARY_DIR}/bin" CACHE PATH "directory for executable file")
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

if(ENABLE_MPI)
    find_package(MPI REQUIRED)
    if(MPI_Fortran_FOUND)
        set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
        include_directories(${MPI_Fortran_INCLUDE_PATH})
        add_definitions(${MPI_Fortran_COMPILE_FLAGS})
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DUSE_MPI")
    endif()
endif(ENABLE_MPI)

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

# Add main executable
set(PROJECT_BINARY ${PROJECT_NAME}.x)
add_executable(${PROJECT_BINARY} ${SRC}/main.f90)

# Add the directory of the sources
add_subdirectory(${SRC})

# Linking with Libraries
if(ENABLE_OPENMP AND ENABLE_MPI AND OpenMP_Fortran_FOUND AND MPI_Fortran_FOUND)
    target_link_libraries(${PROJECT_BINARY} ${PROJECT_NAME} ${BLAS95_LIBRARIES} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES} ${MPI_Fortran_LIBRARIES} OpenMP::OpenMP_Fortran Threads::Threads ${CMAKE_DL_LIBS})
elseif(ENABLE_OPENMP AND OpenMP_Fortran_FOUND)
    target_link_libraries(${PROJECT_BINARY} ${PROJECT_NAME} ${BLAS95_LIBRARIES} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES} OpenMP::OpenMP_Fortran Threads::Threads ${CMAKE_DL_LIBS})
elseif(ENABLE_MPI AND MPI_Fortran_FOUND)
    target_link_libraries(${PROJECT_BINARY} ${PROJECT_NAME} ${BLAS95_LIBRARIES} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES} ${MPI_Fortran_LIBRARIES} Threads::Threads ${CMAKE_DL_LIBS})
else()
    target_link_libraries(${PROJECT_BINARY} ${PROJECT_NAME} ${BLAS95_LIBRARIES} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES} Threads::Threads ${CMAKE_DL_LIBS})
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
# Tests
##############################################################################################
include(CTest)
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
 
##############################################################################################
## ADD Cpack
##############################################################################################
# Include CPack
include(CPack)
INCLUDE(${LOCAL_MODULE_PATH}/git-watcher.cmake)
message(STATUS "Git version info : ${gitversion}.")
# Set package properties
set(CPACK_PACKAGE_NAME "rslmto")
set(CPACK_PACKAGE_VERSION "0.9.0")
#set(CPACK_PACKAGE_VERSION "${gitversion}")
message(STATUS "CPACK_PACKAGE_VERSION: ${CPACK_PACKAGE_VERSION}")
set(CPACK_PACKAGE_CONTACT "ramond@rslmto.org")
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