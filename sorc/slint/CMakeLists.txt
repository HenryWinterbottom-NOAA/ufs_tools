#=======================================================================

#$$$ CMAKE DOCUMENTATION BLOCK

# Copyright (C) 2021 Henry R. Winterbottom

# Email: henry.winterbottom@noaa.gov

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see
# <http://www.gnu.org/licenses/>.

# Review the README, within the top-level directory, which provides
# relevant instructions and (any) references cited by algorithms
# within this software suite.

#=======================================================================

cmake_minimum_required (VERSION 3.18)
message (STATUS "Building with Fortran compiler ${CMAKE_Fortran_COMPILER}.")

if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
   set (CMAKE_Fortran_FLAGS_RELEASE -g -fcheck=all -fbacktrace -Wall)
   set (CMAKE_Fortran_FLAGS_DEBUG "-O0")
endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")

file (GLOB f90_objs *.F90)
set (F90_OBJS ${f90_objs})
add_library (slint SHARED ${f90_objs})
