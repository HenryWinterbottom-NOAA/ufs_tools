# =========================================================================

# Module: ush/remapper/remapio/xarray_interface.py

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

    xarray_interface.py

Description
-----------

    This module contains interfaces to the Python xarray library; this
    interface requires xarray version 0.16.2 or earlier.

Functions
---------

    dataset(ncfile, varobj_list, unlimitdim = None)

        This function defines a xarray dataset object and writes the
        respective variable objects (within the input varobj_list) to
        the specified netCDF formatted file.

    open(ncfile)

        This function opens a netCDF file and returns a Python object
        containing the contents of the respective netCDF file.

    read(ncfile, ncvarname)

        This function parses a netCDF file and returns the attributes
        for the specified netCDF variable.

    varobj(varval, coords, dims, ncvarname)

        This function defines an xarray DataArray object in accordance
        with the specified arguments.

    write(ncfile, var_obj, var_arr)

        This function writes an array of values to an existing variable
        within the specified netCDF file.

Requirements
------------

- ufs_pytils; https://github.com/HenryWinterbottom-NOAA/ufs_pyutils

- xarray; https://github.com/pydata/xarray

Author(s)
---------

    Henry R. Winterbottom; 12 February 2023

History
-------

    2023-02-12: Henry Winterbottom -- Initial implementation.

"""

# ----

# pylint: disable=broad-except
# pylint: disable=redefined-builtin
# pylint: disable=redefined-outer-name

# ----

from typing import Dict, List

import numpy
import xarray
from exceptions import RemapperError
from tools import parser_interface
from utils.logger_interface import Logger

# ----

logger = Logger()

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


def dataset(ncfile: str, varobj_list: List, unlimitdim: str = None) -> None:
    """
    Description
    -----------

    This function defines a xarray dataset object and writes the
    respective variable objects (within the input varobj_list) to the
    specified netCDF formatted file.

    Parameters
    ----------

    ncfile: str

         A Python string specifying the path to the netCDF formatted
         file to be written.

    varobj_list: list

         A Python list of variable objects to be written to the netCDF
         formatted file.

    Keywords
    --------

    unlimitdim: str, optional

         A Python string specifying the coordinate dimension which is
         to be 'unlimited'.

    """
    dataset = xarray.merge(varobj_list)

    kwargs = {}
    if unlimitdim is not None:
        kwargs["unlimited_dims"] = unlimitdim

    try:
        dataset.to_netcdf(ncfile, engine="scipy", **kwargs)

    except Exception:
        msg = "Unable to use xarray engine scipy; trying netCDF4."
        logger.warn(msg=msg)
        dataset.to_netcdf(ncfile, engine="netcdf4", **kwargs)

    dataset.close()


# ----


def open(ncfile: str) -> object:
    """
    Description
    -----------

    This function opens a netCDF file and returns a Python object
    containing the contents of the respective netCDF file.

    Parameters
    ----------

    ncfile: str

        A Python string specifying the path to the netCDF formatted
        file to be opened.

    Returns
    -------

    ncfile_obj: object

        A Python object containing the contents of the input netCDF
        file.

    """

    # Open the Python object containing the netCDF-formatted file
    # contents; proceed accordingly.
    try:
        ncfile_obj = xarray.open_dataset(ncfile, engine="scipy")

    except Exception:
        msg = "Unable to use xarray engine scipy; trying netCDF4."
        logger.warn(msg=msg)
        ncfile_obj = xarray.open_dataset(ncfile, engine="netcdf4")

    return ncfile_obj


# ----


def read(ncfile: str, ncvarname: str) -> object:
    """
    Description
    -----------

    This function parses a netCDF file and returns the attributes for
    the specified netCDF variable.

    Parameters
    ----------

    ncfile: str

        A Python string specifying the path to the netCDF formatted
        file to be opened.

    ncvarname: str

        A Python string specifying the netCDF variable name.

    Returns
    -------

    ncvar_obj: object

        A Python object containing the specified netCDF variable
        attributes.

    """

    # Open the Python object containing the netCDF-formatted file
    # contents; proceed accordingly.
    try:
        ncfile_obj = xarray.open_dataset(ncfile, engine="scipy")

    except Exception:
        msg = "Unable to use xarray engine scipy; trying netCDF4."
        logger.warn(msg=msg)
        ncfile_obj = xarray.open_dataset(ncfile, engine="netcdf4")

    # Read the attributes for the netCDF variable specified upon
    # entry; proceed accordingly.
    ncvar_obj = parser_interface.object_getattr(
        object_in=ncfile_obj, key=ncvarname, force=True
    )
    if ncvar_obj is None:
        msg = (
            f"The netCDF variable {ncvarname} could not be found in netCDF "
            f"file {ncfile}. Aborting!!!"
        )
        raise RemapperError(msg=msg)

    ncfile_obj.close()

    return ncvar_obj


# ----


def varobj(varval: numpy.array, coords: Dict, dims: List, ncvarname: str) -> object:
    """
    Description
    -----------

    This function defines an xarray DataArray object in accordance
    with the specified arguments.

    Parameters
    ----------

    varval: array-type

        A Python array-type variable containing the value array for
        the xarray DataArray object.

    coords: dict

        A Python dictionary containing the coordinate dimension key
        and value pairs for the respective variable.

    dims: list

        A Python list containing the coordinate dimension names for
        the respective variable.

    ncvarname: str

        A Python string specifying the netCDF variable name.

    Returns
    -------

    var_obj: object

        A Python xarray DataArray object.

    """

    # Define the xarray object.
    xarray_obj = xarray.DataArray(varval, coords=coords, dims=dims)
    var_obj = xarray_obj.to_dataset(name=ncvarname)
    xarray_obj.close()

    return var_obj


# ----


def write(ncfile: str, var_obj: object, var_arr: numpy.array) -> None:
    """
    Description
    -----------

    This function writes an array of values to an existing variable
    within the specified netCDF file.

    Parameters
    ----------

    ncfile: str

        A Python string specifying the path to the netCDF formatted
        file to be opened.

    var_obj: object

        A Python xarray DataArray object.

    var_arr: array-type

        A Python array-type variable containing the value array to be
        written to the existing variable within the netCDF file.

    """

    # Open the Python object containing the netCDF-formatted file
    # contents; proceed accordingly.
    try:
        dataset = xarray.open_dataset(ncfile, engine="scipy")

    except Exception:
        msg = "Unable to use xarray engine scipy; trying netCDF4."
        logger.warn(msg=msg)
        dataset = xarray.open_dataset(ncfile, engine="netcdf4")

    # Collect the attributes for the respective netCDF variable.
    ncvarname = parser_interface.object_getattr(
        object_in=var_obj, key="ncvarname", force=True
    )

    if ncvarname is None:
        msg = (
            "The netCDF variable name could not be determined "
            "from the input variable attributes. Aborting!!!"
        )
        raise RemapperError(msg=msg)

    msg = f"Writing variable {ncvarname} to {ncfile}."
    logger.info(msg=msg)

    #    ncdims = dataset[ncvarname].dims
    #    nccords = dataset[ncvarname].coords

    # Write the variable values to the netCDF-formatted file path.
    dataset[ncvarname][:] = var_arr[:]
    dataset.to_netcdf(ncfile, engine="scipy")
    dataset.close()
