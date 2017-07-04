cmake_minimum_required(VERSION 2.6)
project(celltrack)

enable_language (Fortran)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-loops -fforce-addr -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fbounds-check -fbacktrace -Wall -Wextra -ggdb -fcheck=all")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
else ()
  message (FATAL_ERROR "Compiler gfortran not found! Stop!")
endif ()

FIND_PATH(CDI_INCLUDE
  NAMES
  cdi.inc
  PATHS
  $ENV{CDI_INCLUDE_PATH}
  /usr/include
  $ENV{HOME}/include
  DOC "cdi include directory (must contain cdi.inc)"
)

FIND_LIBRARY(CDI_LIB
  NAMES
  libcdi.so
  PATHS
  $ENV{CDI_LIB_PATH}
  /usr/lib
  /usr/lib64
  $ENV{HOME}/lib
  $ENV{HOME}/lib64
  DOC "cdi library"
)

add_subdirectory(src)
