if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif()

cmake_policy(SET CMP0011 NEW)
cmake_policy(SET CMP0012 NEW)
cmake_policy(SET CMP0048 NEW)
cmake_policy(SET CMP0054 NEW)

cmake_minimum_required(VERSION 3.5)

set(PROJECT_NAME rslmto)
project(${PROJECT_NAME})

# Find all the .in files
file(GLOB inFiles RELATIVE "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/src/*.in" "${PROJECT_SOURCE_DIR}/test/*.in")
foreach(infileName ${inFiles})
  string(REGEX REPLACE ".in\$" "" outfileName "${infileName}")
  #message(STATUS "Stamp ${infileName} => ${outfileName}")
  configure_file(${PROJECT_SOURCE_DIR}/${infileName} ${PROJECT_SOURCE_DIR}/${outfileName})
endforeach(infileName)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules/")

# Enable Fortran language
enable_language(Fortran)
# Force the use of the local version of the compilers
#find_program(CMAKE_Fortran_COMPILER NAMES $ENV{FC} gfortran PATHS ENV PATH NO_DEFAULT_PATH)
message(STATUS "Using Fortran compiler ${CMAKE_Fortran_COMPILER_ID} Fortran ${CMAKE_Fortran_COMPILER_VERSION}")

# Fortran module and output directories
set(fortran_module_directory "${CMAKE_BINARY_DIR}/modules" CACHE PATH "directory for fortran module files")
file(MAKE_DIRECTORY ${fortran_module_directory})
set(CMAKE_Fortran_MODULE_DIRECTORY ${fortran_module_directory})
set(fortran_output_directory "${CMAKE_BINARY_DIR}/bin" CACHE PATH "directory for executable file")
file(MAKE_DIRECTORY ${fortran_output_directory})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${fortran_output_directory})

##############################################################################################
# Determine and set compiler-dependent flags
##############################################################################################
set(LOCAL_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")
include(${LOCAL_MODULE_PATH}/SetFortranFlags.cmake)
##############################################################################################
# # Enable the preprocessor globally
# set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")

# Generate position independent code (PIC) suitable for use in shared libraries
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")
set(CMAKE_VERBOSE_MAKEFILE off)

# add the directory of the sources
add_subdirectory(source)
