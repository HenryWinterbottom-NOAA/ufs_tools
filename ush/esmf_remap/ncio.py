# =========================================================================

# Module: ush/esmf_remap/ncio.py

# Author: Henry R. Winterbottom

# Email: henry.winterbottom@noaa.gov

# This program is free software: you can redistribute it and/or modify
# it under the terms of the respective public license published by the
# Free Software Foundation and included with the repository within
# which this application is contained.

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

# =========================================================================

"""
Module
------

    ncio.py

Description
-----------

    This module contains functions to compute/define remapping
    attributes for either source or destination grids collected from
    Network Common Data Format (netCDF) formatted files.

Functions
---------

    define_vararray(srcgrid_info_obj)

        This function defines and returns an array of zeros of the
        same dimension as the input file variable grids; this array
        can then be used to compute the remapped variable and
        subsequently produce the ESMF formatted file containing the
        remapping attributes.

    input_grids(grid_obj)

        This function defines a grid object to be used by xESMF to
        compute the respective ESMF remapping attributes; this method
        is identical for both the source and destination grid
        configurations.

Requirements
------------

- ufs_pytils; https://github.com/HenryWinterbottom-NOAA/ufs_pyutils

Author(s)
---------

    Henry R. Winterbottom; 09 February 2023

History
-------

    2023-02-09: Henry Winterbottom -- Initial implementation.

"""

# ----

import numpy
import xarray
from exceptions import ESMFRemapError
from ioapps import netcdf4_interface
from tools import parser_interface
from utils.logger_interface import Logger

# ----

logger = Logger()

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


def define_vararray(srcgrid_info_obj: object) -> numpy.array:
    """
    Description
    -----------

    This function defines and returns an array of zeros of the same
    dimension as the input file variable grids; this array can then be
    used to compute the remapped variable and subsequently produce the
    ESMF formatted file containing the remapping attributes.

    Parameters
    ----------

    srcgrid_info_obj: object

        A Python object containing the respective grid (e.g., source
        or destination) attributes.

    Returns
    -------

    vararray: array-type

        A Python array of zeros of the same dimension as the input
        file variable grids.

    Raises
    ------

    ESMFRemapError:

        * raised if the source grid attributes, collected from the
          experiment configuration file, does not contain a required
          dimension attribute.

    """

    # Initialize and build the variable array; proceed accordingly.
    vararray_obj = parser_interface.object_define()
    array_attrs_list = ["ncxdim_name", "ncydim_name"]

    # Collect the array attributes from the experiment configuration
    # file.
    for array_attr in array_attrs_list:
        value = parser_interface.object_getattr(
            object_in=srcgrid_info_obj, key=array_attr, force=True
        )
        if value is None:
            msg = (
                "The source grid attributes does not contain the "
                f"attribute {array_attr}. Aborting!!!"
            )
            raise ESMFRemapError(msg=msg)

        # Define the array attributes.
        ncdim = netcdf4_interface.ncreaddim(
            ncfile=srcgrid_info_obj.ncfile, ncdimname=value
        )
        vararray_obj = parser_interface.object_setattr(
            object_in=vararray_obj, key=array_attr, value=ncdim
        )

    msg = (
        f"Defining an empty array of dimension ({vararray_obj.ncydim_name}, "
        f"{vararray_obj.ncxdim_name})."
    )
    logger.info(msg=msg)

    # Define an array of zeros.
    vararray = numpy.zeros(
        [vararray_obj.ncydim_name, vararray_obj.ncxdim_name])

    return vararray


# ----


def input_grids(grid_obj: object) -> object:
    """
    Description
    -----------

    This function defines a grid object to be used by xESMF to compute
    the respective ESMF remapping attributes; this method is identical
    for both the source and destination grid configurations.

    Parameters
    ----------

    grid_obj: object

        A Python object containing the respective grid (e.g., source
        or destination) attributes.

    Returns
    -------

    grid: object

        A Python xarray object containing the respective grid (e.g.,
        source or destination) netCDF file attributes and the
        respective netCDF geographical coordinate variable names in
        accordance with the xESMF application expections.

    Raises
    ------

    ESMFRemapError:

        * raised if the netCDF-formatted file path, containing the
          respective grid configuration and attributes, is not defined
          within the experiment configuration.

        * raised if a required netCDF grid coordinate cannot be
          determined from the experiment configuration.

    """

    # Define the netCDF-formatted file path.
    ncfile = parser_interface.object_getattr(
        object_in=grid_obj, key="ncfile", force=True
    )
    if ncfile is None:
        msg = (
            "The grid input netCDF file is either NoneType or undefined. Aborting!!!"
        )
        raise ESMFRemapError(msg=msg)

    # Collect the attributes from the netCDF-formatted file path and
    # define the xarray object; proceed accordingly.
    msg = f"Reading netCDF formatted file {ncfile}."
    logger.info(msg=msg)

    # Define the xarray grid attributes; proceed accordingly.
    grid = xarray.open_dataset(ncfile)
    grid_ncattrs_dict = {}

    # Build the coordinate attributes for the xarray object.
    coords_dict = {"nclat": "lat", "nclon": "lon"}

    for coord in coords_dict:
        coord_name = parser_interface.object_getattr(
            object_in=grid_obj, key=coord, force=True
        )
        if coord_name is None:
            msg = (
                f"The netCDF grid coordinate {coord} could not be "
                "determined from the experiment configuration. Aborting!!!"
            )
            raise ESMFRemapError(msg=msg)

        # Update the xarray grid attribute.
        grid_ncattrs_dict[coord_name] = parser_interface.dict_key_value(
            dict_in=coords_dict, key=coord, no_split=True
        )

    # Update the xarray object and close the open netCDF-formatted
    # file path.
    grid = grid.rename(**grid_ncattrs_dict)
    grid.close()

    msg = (
        f"The netCDF longitudinal coordinate dimension is {grid.lon.shape}."
    )
    logger.info(msg=msg)

    msg = f"The netCDF latitudinal coordinate dimension is {grid.lat.shape}."
    logger.info(msg=msg)

    return grid
