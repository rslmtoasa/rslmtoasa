######################################################
# Determine and set the Fortran compiler flags we want 
######################################################

message(STATUS "Setting compiler flags.")
####################################################################
# Make sure that the default build type is Release if not specified.
####################################################################
INCLUDE(${LOCAL_MODULE_PATH}/SetCompileFlag.cmake)

set(RSLMTO_BUILD_TYPES Debug Release RelWithDebInfo Test)
set(RSLMTO_BUILD_TYPE_HELP
    "Choose the type of build, options are: Debug, Release, RelWithDebInfo, Test.")

if(CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_CONFIGURATION_TYPES "${RSLMTO_BUILD_TYPES}" CACHE STRING
        "${RSLMTO_BUILD_TYPE_HELP}" FORCE)
else()
    string(TOUPPER "${CMAKE_BUILD_TYPE}" BT)

    if(NOT BT)
        set(CMAKE_BUILD_TYPE Release CACHE STRING "${RSLMTO_BUILD_TYPE_HELP}" FORCE)
        message(STATUS "CMAKE_BUILD_TYPE not given, defaulting to Release")
    elseif(BT STREQUAL "RELEASE")
        set(CMAKE_BUILD_TYPE Release CACHE STRING "${RSLMTO_BUILD_TYPE_HELP}" FORCE)
    elseif(BT STREQUAL "DEBUG")
        set(CMAKE_BUILD_TYPE Debug CACHE STRING "${RSLMTO_BUILD_TYPE_HELP}" FORCE)
    elseif(BT STREQUAL "RELWITHDEBINFO")
        set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "${RSLMTO_BUILD_TYPE_HELP}" FORCE)
    elseif(BT STREQUAL "RELWITHDEB" OR
           BT STREQUAL "RELEASEWITHDEB" OR
           BT STREQUAL "RELEASWITHDEB" OR
           BT STREQUAL "RELEASEWITHDEBINFO" OR
           BT STREQUAL "RELEASEWITHDEBUGINFO")
        set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "${RSLMTO_BUILD_TYPE_HELP}" FORCE)
        message(STATUS "Normalized CMAKE_BUILD_TYPE=${BT} to RelWithDebInfo")
    elseif(BT STREQUAL "TEST" OR BT STREQUAL "TESTING")
        set(CMAKE_BUILD_TYPE Test CACHE STRING "${RSLMTO_BUILD_TYPE_HELP}" FORCE)
        if(BT STREQUAL "TESTING")
            message(STATUS "Normalized legacy CMAKE_BUILD_TYPE=TESTING to Test")
        endif()
    else()
        message(FATAL_ERROR
            "CMAKE_BUILD_TYPE='${CMAKE_BUILD_TYPE}' is not valid. "
            "Choices are: Debug, Release, RelWithDebInfo, Test.")
    endif()

    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${RSLMTO_BUILD_TYPES})
endif()

#########################################################
# If the compiler flags have already been set, return now
#########################################################
#
IF(CMAKE_Fortran_FLAGS_RELEASE AND CMAKE_Fortran_FLAGS_RELWITHDEBINFO AND CMAKE_Fortran_FLAGS_TEST AND CMAKE_Fortran_FLAGS_DEBUG)
   #unset(CMAKE_Fortran_FLAGS CACHE)
    RETURN ()
ENDIF(CMAKE_Fortran_FLAGS_RELEASE AND CMAKE_Fortran_FLAGS_RELWITHDEBINFO AND CMAKE_Fortran_FLAGS_TEST AND CMAKE_Fortran_FLAGS_DEBUG)

########################################################################
# Determine the appropriate flags for this compiler for each build type.
# For each option type, a list of possible flags is given that work
# for various compilers.  The first flag that works is chosen.
# If none of the flags work, nothing is added (unless the REQUIRED 
# flag is given in the call).  This way unknown compiles are supported.
#######################################################################

#####################
### GENERAL FLAGS ###
#####################


## Don't add underscores in symbols for C-compatability
#SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
#                 Fortran "-fno-underscoring")

# No limits on line-lengths with GNU
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                 Fortran "-ffree-line-length-0")
#SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
#                 Fortran "-std=legacy")

message(STATUS "CMAKE_Fortran_FLAGS_RELEASE        : ${CMAKE_Fortran_FLAGS_RELEASE}")
message(STATUS "CMAKE_Fortran_FLAGS_RELWITHDEBINFO : ${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}")
message(STATUS "CMAKE_Fortran_FLAGS_TEST           : ${CMAKE_Fortran_FLAGS_TEST}")
message(STATUS "CMAKE_Fortran_FLAGS_DEBUG          : ${CMAKE_Fortran_FLAGS_DEBUG}")
message(STATUS "CMAKE_Fortran_FLAGS                : ${CMAKE_Fortran_FLAGS}")

if(USE_VSL)
   SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                 Fortran "-fno-range-check")
endif()

# Ensure that preprocessor flags are invoked
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                 Fortran "-cpp"
                         "-fpp")

# Optimize for the host's architecture (default: native -march=native).
# For portable/CI builds where the binary may run on a different machine,
# pass -DENABLE_MARCH_NATIVE=OFF to use generic baseline (-mtune=generic).
option(ENABLE_MARCH_NATIVE "Optimize for native CPU architecture (-march=native)" ON)
IF(APPLE)
    SET(GNUNATIVE "-mtune=native")
ELSEIF(NOT ENABLE_MARCH_NATIVE)
    # Portable/CI: generic baseline
    SET(GNUNATIVE "-mtune=generic")
ELSE()
    # Default: native optimization
    SET(GNUNATIVE "-march=native")
ENDIF()
# Optimize for the host's architecture
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                 Fortran "-xHost"        # Intel
                         ${GNUNATIVE}    # GNU
                         "-ta=host"      # Portland Group
                         "/QxHost"       # Intel Windows
                )

###################
### DEBUG FLAGS ###
###################

# NOTE: debugging symbols (-g or /debug:full) are already on by default

# Disable optimizations
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                           Fortran REQUIRED "-O0" # All compilers not on Windows
                                            "/Od" # Intel Windows
                )

# Turn on all warnings 
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-warn all" # Intel
                         "-Minform=warn" #Portland
                         "-Wall"     # GNU
                          "/warn:all" # Intel Windows
                )

# Traceback
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-traceback"   # Intel/Portland Group
                         "-fbacktrace"  # GNU (gfortran)
                         "-ftrace=full" # GNU (g95)
                          "/traceback"   # Intel Windows
                )

SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran "-traceback"   # Intel/Portland Group
                         "-fbacktrace -g"  # GNU (gfortran)
                         "-ftrace=full" # GNU (g95)
                          "/traceback"   # Intel Windows
                )

SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran "-g"
                )

# Check array bounds
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-check all"   # Intel
                         "-fcheck=all"  # GNU (bounds + do + mem + pointer + recursion)
                         "-Mbounds"     # Portland Group
                          "/check:all"   # Intel Windows
                )

#####################
### TEST FLAGS ###
#####################

# Optimizations
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_TEST "${CMAKE_Fortran_FLAGS_TEST}"
                 Fortran REQUIRED "-O2" # All compilers not on Windows
                                  "/O2" # Intel Windows
                )
# Profiling
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_TEST "${CMAKE_Fortran_FLAGS_TEST}"
                 Fortran "-fprofile-arcs" # Profiling
                )

# Coverage
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_TEST "${CMAKE_Fortran_FLAGS_TEST}"
                 Fortran "-ftest-coverage" # Profiling
                )

# Preserve compatibility with existing build trees or scripts that still use
# -DCMAKE_BUILD_TYPE=TESTING.
set(CMAKE_Fortran_FLAGS_TESTING "${CMAKE_Fortran_FLAGS_TEST}" CACHE STRING
    "Legacy alias for CMAKE_Fortran_FLAGS_TEST" FORCE)

#############################
### RELWITHDEBINFO FLAGS  ###
#############################

# Optimized build with symbols, intended for CLI profilers and debuggers.
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}"
                 Fortran REQUIRED "-O2" # All compilers not on Windows
                                  "/O2" # Intel Windows
                )

SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}"
                 Fortran "-g"
                         "/debug:full"
                )

SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}"
                 Fortran "-traceback"   # Intel/Portland Group
                         "-fbacktrace"  # GNU (gfortran)
                         "-ftrace=full" # GNU (g95)
                          "/traceback"  # Intel Windows
                )

#####################
### RELEASE FLAGS ###
#####################

# NOTE: agressive optimizations (-O3) are already turned on by default
# For GNU, force conservative optimization in RELEASE to match stable runtime
# behavior observed in DEBUG for strux/SPDF workflows.
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                   Fortran REQUIRED "-O0")
endif()

# Optimizations
#SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
#                 Fortran REQUIRED "-fastsse" # All compilers not on Windows
#                                  "-Ofast" # All compilers not on Windows
#                                  "-O3" # All compilers not on Windows
#                                  "/O3" # Intel Windows
#                )


# Unroll/inline are disabled for GNU release builds due numerical-instability
# regressions in LMTO47 screening. Keep them for non-GNU compilers.
if (NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                   Fortran "-unroll"        # Intel
                           "-Munroll"       # Portland Group
                           "/unroll"        # Intel Windows
                  )

  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                   Fortran "-inline"            # Intel
                           "-Minline"           # Portland Group
                           "/Qinline"           # Intel Windows
                  )
endif()

             ## Interprocedural (link-time) optimizations
             #SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
             #                 Fortran "-ipo"     # Intel
             #                         "/Qipo"    # Intel Windows
             #                         "-Mipa=fast"    # Portland Group
             #)

# Single-file optimizations. ifx/IntelLLVM no longer supports -ip.
if (NOT CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                   Fortran "-ip"  # classic Intel
                   "-Mnoipa"    # Portland
                   "/Qip" # Intel Windows
                  )
endif()

# Vectorize code
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran "-qopt-report0"  # Intel
                         "-Mvect"        # Portland Group
                         "/Qvec-report0" # Intel Windows
                )

             
mark_as_advanced(CMAKE_Fortran_FLAGS_TEST CMAKE_Fortran_FLAGS_TESTING)
list(REMOVE_DUPLICATES CMAKE_Fortran_FLAGS_RELEASE)
list(REMOVE_DUPLICATES CMAKE_Fortran_FLAGS_RELWITHDEBINFO)
list(REMOVE_DUPLICATES CMAKE_Fortran_FLAGS_TEST)
list(REMOVE_DUPLICATES CMAKE_Fortran_FLAGS_DEBUG)
list(REMOVE_DUPLICATES CMAKE_Fortran_FLAGS)
message(STATUS "CMAKE_Fortran_FLAGS_RELEASE        : ${CMAKE_Fortran_FLAGS_RELEASE}")
message(STATUS "CMAKE_Fortran_FLAGS_RELWITHDEBINFO : ${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}")
message(STATUS "CMAKE_Fortran_FLAGS_TEST           : ${CMAKE_Fortran_FLAGS_TEST}")
message(STATUS "CMAKE_Fortran_FLAGS_DEBUG          : ${CMAKE_Fortran_FLAGS_DEBUG}")
message(STATUS "CMAKE_Fortran_FLAGS                : ${CMAKE_Fortran_FLAGS}")
