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
project (CONTAINERS_STACK)

enable_language (Fortran)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CONTAINERS_STACK_SOURCE_DIR}/build/lib)
set (CMAKE_Fortran_MODULE_DIRECTORY ${CONTAINERS_STACK_SOURCE_DIR}/build/include)

add_subdirectory (fson)
add_subdirectory (gridprojs)
add_subdirectory (slint)
