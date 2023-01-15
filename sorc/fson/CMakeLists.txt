#=======================================================================

# CMAKE: sorc/fson/CMakeLists.txt

# Author: Henry R. Winterbottom

# Email: henry.winterbottom@noaa.gov

# This program is free software: you can redistribute it and/or modify
# it under the terms of the respective public license published by the
# Free Software Foundation and included with the repository within
# which this application is contained.

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#=======================================================================

message (STATUS "Building with Fortran compiler ${CMAKE_Fortran_COMPILER}.")

if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
   set (CMAKE_Fortran_FLAGS_RELEASE -g -fcheck=all -fbacktrace -Wall)
   set (CMAKE_Fortran_FLAGS_DEBUG "-O0")
endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")

file (GLOB f90_objs *.f90)
set (F90_OBJS ${f90_objs})
add_library (fson SHARED ${f90_objs})
