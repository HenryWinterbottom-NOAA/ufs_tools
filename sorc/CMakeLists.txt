#=======================================================================

# CMAKE: sorc/CMakeLists.txt

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

cmake_minimum_required (VERSION 2.0)
project (UFS_TOOLS_STACK)

enable_language (Fortran)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${UFS_TOOLS_STACK_SOURCE_DIR}/build/lib)
set (CMAKE_Fortran_MODULE_DIRECTORY ${UFS_TOOLS_STACK_SOURCE_DIR}/build/include)

add_subdirectory (fson)
add_subdirectory (gridprojs)
add_subdirectory (slint)
