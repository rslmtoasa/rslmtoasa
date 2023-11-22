set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# Add library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY lib)
#set(CMAKE_POSITION_INDEPENDENT_CODE True)
add_library(${PROJECT_NAME} STATIC
    array.f90
    bands.f90
    calculation.f90
    charge.f90
    control.f90
    density_of_states.f90
    element.f90
    energy.f90
    exchange.f90
    face.F90
    globals.f90
    green.f90
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
    recursion.f90
    report.f90
    safe_alloc.f90
    self.f90
    spin_dynamics.f90
    string.f90
    symbolic_atom.f90
    timer.f90
    xc.f90
    include_codes/abspinlib/depondt.f90
    include_codes/abspinlib/constants.f90
    include_codes/abspinlib/mtprng.f90
    include_codes/abspinlib/stdtypes.f90
    include_codes/abspinlib/randomnumbers.f90
    include_codes/abspinlib/constants.f90
    include_codes/abspinlib/parameters.f90
)

# Add main executable
add_executable(
  ${PROJECT_NAME}.x
  main.f90
)

if(USE_SAFE_ALLOC)
  target_compile_definitions(${PROJECT_NAME}.x PUBLIC USE_SAFE_ALLOC)
endif()

if (COLOR)
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DCOLOR")
endif(COLOR)

## Enable the check bounds
#if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
#  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -ffree-line-length-0 -g -fbacktrace -cpp")
#endif()

if (ENABLE_FLUSH)
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DEnable_flush")
endif(ENABLE_FLUSH)

if (ENABLE_OPENMP)
  find_package(OpenMP)
  message(STATUS "OpenMP libraries: ${OpenMP_Fortran_LIBRARIES}")
  if (OpenMP_Fortran_FOUND)
    # set_config_option(HAVE_OPENMP "OpenMP")
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
endif()

set(BLA_F95 OFF) #put on if Intel LAPACK 95 interface is needed

find_package(BLAS REQUIRED)
  if(BLAS_FOUND)
  message(STATUS "BLAS libraries: ${BLAS_LIBRARIES}")
  # set_config_option(HAVE_BLAS "BLAS")
endif()

find_package(LAPACK REQUIRED)
  if(LAPACK_FOUND)
  message(STATUS "LAPACK libraries: ${LAPACK_LIBRARIES}")
  # set_config_option(HAVE_LAPACK "LAPACK")
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

# MKL machine dependant stuffs
# check https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor to comply with
if (BLAS_FOUND AND BLAS95_FOUND)
  if (BLA_VENDOR MATCHES "^Intel")
  if(NOT DEFINED ENV{MKL_IL})
    set(ENV{MKL_IL} ilp64)
  elseif(NOT $ENV{MKL_IL} STREQUAL lp64 OR NOT $ENV{MKL_IL} STREQUAL ilp64)
  message(FATAL_ERROR "Unvalid option for MKL_IL")
endif()

# test existance of directories
if(NOT EXISTS "$ENV{MKLROOT}/include")
  set(ENV{MKLROOT} "$ENV{MKLROOT}/mkl")
endif()

# 32-bits integer
if ($ENV{MKL_IL} STREQUAL lp64)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64 -I$ENV{MKLROOT}/include -I$ENV{MKLROOT}/include/intel64/lp64")
endif()

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
# 64-bits integer
  if ($ENV{MKL_IL} STREQUAL ilp64)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64 -fdefault-integer-8 -I$ENV{MKLROOT}/include -I$ENV{MKLROOT}/include/intel64/ilp64")
  endif()
endif()

if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
# 64-bits integer
if ($ENV{MKL_IL} STREQUAL ilp64)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64 -i8 -I$ENV{MKLROOT}/include -I$ENV{MKLROOT}/include/intel64/ilp64")
endif()
if (APPLE)
set(MKL_LIBRARY_DIRS "$ENV{MKLROOT}/lib")
find_library(MKL_CORE      mkl_core.a               HINTS ${MKL_LIBRARY_DIRS})
find_library(MKL_COMPILER  mkl_intel_$ENV{MKL_IL}.a HINTS ${MKL_LIBRARY_DIRS})
find_library(MKL_THREADING mkl_intel_thread.a       HINTS ${MKL_LIBRARY_DIRS})
endif()
if (UNIX AND NOT APPLE)
set(MKL_LIBRARY_DIRS "$ENV{MKLROOT}/lib/intel64")
find_library(MKL_CORE      libmkl_core.a               HINTS ${MKL_LIBRARY_DIRS})
find_library(MKL_COMPILER  libmkl_intel_$ENV{MKL_IL}.a HINTS ${MKL_LIBRARY_DIRS})
find_library(MKL_THREADING libmkl_intel_thread.a       HINTS ${MKL_LIBRARY_DIRS})
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-intel -Wl, --start-group ${MKL_CORE} ${MKL_COMPILER} ${MKL_THREADING} -Wl, --end-group")
endif()
endif()
endif()
endif()
#if(ENABLE_OPENMP)
#  target_link_libraries(${PROJECT_NAME}.x ${PROJECT_NAME} ${BLAS95_LIBRARIES} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES} OpenMP::OpenMP_Fortran Threads::Threads ${CMAKE_DL_LIBS})
  #target_link_libraries(bands.x ${PROJECT_NAME} ${BLAS95_LIBRARIES} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES} OpenMP::OpenMP_Fortran Threads::Threads ${CMAKE_DL_LIBS})
  #target_link_libraries(pdos.x ${PROJECT_NAME} ${BLAS95_LIBRARIES} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES} OpenMP::OpenMP_Fortran Threads::Threads ${CMAKE_DL_LIBS})
#else(ENABLE_OPENMP)
#  target_link_libraries(${PROJECT_NAME}.x ${PROJECT_NAME} ${BLAS95_LIBRARIES} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES} Threads::Threads ${CMAKE_DL_LIBS})
  #target_link_libraries(bands.x ${PROJECT_NAME} ${BLAS95_LIBRARIES} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES} Threads::Threads ${CMAKE_DL_LIBS})
  #target_link_libraries(pdoss.x ${PROJECT_NAME} ${BLAS95_LIBRARIES} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES} Threads::Threads ${CMAKE_DL_LIBS})
#endif(ENABLE_OPENMP)

# Linking with Libraries

# Both OpenMP and MPI
if(ENABLE_OPENMP AND ENABLE_MPI AND OpenMP_Fortran_FOUND AND MPI_Fortran_FOUND)
    target_link_libraries(${PROJECT_NAME}.x ${PROJECT_NAME} ${BLAS95_LIBRARIES} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES} ${MPI_Fortran_LIBRARIES} OpenMP::OpenMP_Fortran Threads::Threads ${CMAKE_DL_LIBS})

# Only OpenMP
elseif(ENABLE_OPENMP AND OpenMP_Fortran_FOUND)
    target_link_libraries(${PROJECT_NAME}.x ${PROJECT_NAME} ${BLAS95_LIBRARIES} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES} OpenMP::OpenMP_Fortran Threads::Threads ${CMAKE_DL_LIBS})

# Only MPI
elseif(ENABLE_MPI AND MPI_Fortran_FOUND)
    target_link_libraries(${PROJECT_NAME}.x ${PROJECT_NAME} ${BLAS95_LIBRARIES} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES} ${MPI_Fortran_LIBRARIES} Threads::Threads ${CMAKE_DL_LIBS})

# Neither OpenMP nor MPI
else()
    target_link_libraries(${PROJECT_NAME}.x ${PROJECT_NAME} ${BLAS95_LIBRARIES} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES} Threads::Threads ${CMAKE_DL_LIBS})
endif()

