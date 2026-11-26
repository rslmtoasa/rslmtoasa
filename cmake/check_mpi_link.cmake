# Usage: cmake -DBINARY=<executable> -P check_mpi_link.cmake
# Fails if the executable loads libmpi.so with more than one soname (two MPI installations mixed).
execute_process(COMMAND ldd ${BINARY} OUTPUT_VARIABLE _ldd RESULT_VARIABLE _rc ERROR_QUIET)
if(NOT _rc EQUAL 0)
    return()
endif()
string(REGEX MATCHALL "libmpi\\.so\\.[0-9]+ => [^ ]+" _libs "${_ldd}")
list(REMOVE_DUPLICATES _libs)
set(_sonames "")
foreach(_l IN LISTS _libs)
    string(REGEX REPLACE " =>.*" "" _s "${_l}")
    list(APPEND _sonames ${_s})
endforeach()
list(REMOVE_DUPLICATES _sonames)
list(LENGTH _sonames _n)
if(_n GREATER 1)
    list(JOIN _libs "\n  " _found)
    message(FATAL_ERROR "${BINARY} loads libmpi from more than one MPI installation:\n  ${_found}\n"
        "Configure and build in an environment with a single MPI (unset LIBRARY_PATH and LD_LIBRARY_PATH entries of "
        "other MPIs, e.g. oneAPI setvars), with -DMPI_Fortran_COMPILER and -DMPIEXEC_EXECUTABLE from the same installation.")
endif()
