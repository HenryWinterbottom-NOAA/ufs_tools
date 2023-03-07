# =========================================================================

# Module: ush/tcdiags/interp/__init__.py

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

    __init__.py

Description
-----------

    This module functions to handle various interpolation types.

Functions
---------

    interp_ll2ra(varin, lats, lons, lat_0, lon_0, max_radius, drho,
                 dphi)

        This function interpolates a 2-dimensional variable, defined
        on a Cartesian type grid, to a polar projection grid defined
        by the attributes for the polar projection specified upon
        entry.

    interp_vertical(varin, zarr, levs)

        This method interpolates a 3-dimensional variable to specified
        vertical levels.

Requirements
------------

- ufs_pytils; https://github.com/HenryWinterbottom-NOAA/ufs_pyutils

Author(s)
---------

    Henry R. Winterbottom; 07 March 2023

History
-------

    2023-03-07: Henry Winterbottom -- Initial implementation.

"""

# ----

# pylint: disable=invalid-name
# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals

# ----

from typing import List

import numpy
from exceptions import TCDiagsError
from tcdiags.geomets import haversine
from tools import parser_interface
from utils.logger_interface import Logger
from wrf import interplevel

# ----

# Define all available functions.
__all__ = ["interp_ll2ra", "interp_vertical"]

# ----

logger = Logger()

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


def interp_ll2ra(
    varin: numpy.array,
    lats: numpy.array,
    lons: numpy.array,
    lat_0: float,
    lon_0: float,
    max_radius: float,
    drho: float,
    dphi: float,
) -> object:
    """
    Description
    -----------

    This function interpolates a 2-dimensional variable, defined on a
    Cartesian type grid, to a polar projection grid defined by the
    attributes for the polar projection specified upon entry.

    Parameters
    ----------

    varin: array-type

        A Python array-type variable containing the 2-dimensional
        variable defined along the respective Cartesian grid.

    lats: array-type

        A Python array-type variable containing the 2-dimensional grid
        of latitude coordinate values; the coordinate values are
        assumed order south to north; units are degrees.

    lons: array-type

        A Python array-type variable containing the 2-dimensional grid
        of longitude coordinate values; the coordinate values are
        assumed to be within in the range[-180.0 to 180.0]; units are
        degrees.

    lat_0: float

        A Python float value defining the reference latitude
        coordinate value from which the polar grid projection will be
        defined; units are degrees.

    lon_0: float

        A Python float value defining the reference longitude
        coordinate value from which the polar grid projection will be
        defined; the coordinate values are assumed to be within in the
        range[-180.0 to 180.0]; units are degrees;

    max_radius: float

        A Python float value defining the maximum radial distance for
        which to define the polar grid projection; units are meters.

    drho: float

        A Python float variable defining the radial distance interval
        for the polar projection; units are meters.

    dphi: float

        A Python float value defining the aximuthal interval for the
        polar projection; units are degrees.

    Returns
    -------

    varout_obj: object

        A Python object containing the interpolated variable as well
        as the attributes of the polar projection.

    """

    # Initialize the coordinate arrays.
    varin = numpy.ravel(varin)
    lats = numpy.ravel(lats)
    lons = numpy.ravel(lons)
    dphi = numpy.radians(dphi)

    # Compute the radial distance relative to the specified
    # geographical coordinate location.
    fix = (lat_0, lon_0)
    rho = numpy.array(
        [haversine(fix, (lats[idx], lons[idx])) for idx in range(len(lats))]
    )

    xx = numpy.array([haversine(fix, (lat_0, lons[idx])) for idx in range(len(lats))])
    xx = numpy.where(lons < lon_0, -1.0 * xx, xx)

    yy = numpy.array([haversine(fix, (lats[idx], lon_0)) for idx in range(len(lats))])
    yy = numpy.where(lats < lat_0, -1.0 * yy, yy)

    phi = numpy.arctan2(yy, xx)

    # Interpolate the Cartesian grid to the defined polar coordinate
    # grid.
    radial = numpy.arange(0.0, (max_radius + drho), drho)
    azimuth = numpy.arange(-1.0 * numpy.pi, (numpy.pi + dphi), dphi)

    msg = (
        f"Defining polar projection grid of resolution {drho} meters "
        f"and {dphi} degrees centered at longitude coordinate {lon_0} "
        f"and latitude coordinate location {lat_0}."
    )
    logger.info(msg=msg)

    # Interpolate the variable defined on the Cartesian grid to the
    # established the radial and azimuthal angle coordinates; proceed
    # accordingly.
    var = []

    # Define the radial coordinate.
    for ridx in enumerate(radial):
        radii = radial[ridx[0]]

        # Define the azimuthal coordinate.
        for aidx in enumerate(azimuth):
            theta = azimuth[aidx[0]]
            array = []

            # Determine all (if any) input variable values within the
            # radii and azimuth (theta) interval; proceed accordingly.
            idxs = numpy.ndarray.tolist(numpy.where(rho >= radii)[0])
            for idx in idxs:
                if rho[idx] < (radii + drho):
                    if (phi[idx] >= theta) and (phi[idx] < (theta + dphi)):
                        array.append(varin[idx])

            if len(array) == 0:
                var.append(numpy.nan)

            else:
                var.append(numpy.nanmean(array))

    # Interpolate in order to fill and missing (i.e., NaN) values.
    interp_var = numpy.array(var)

    check = numpy.logical_not(numpy.isnan(interp_var))
    xp = check.ravel().nonzero()[0]
    fp = interp_var[numpy.logical_not(numpy.isnan(interp_var))]
    x = numpy.isnan(interp_var).ravel().nonzero()[0]

    interp_var[numpy.isnan(var)] = numpy.interp(x, xp, fp)
    (nrho, nphi) = [len(radial), len(azimuth)]

    # Define and build the output variable object attributes.
    varout_dict = {
        "azimuth": azimuth,
        "dphi": dphi,
        "drho": drho,
        "lat_0": lat_0,
        "lon_0": lon_0,
        "max_radius": max_radius,
        "nphi": nphi,
        "nrho": nrho,
        "radial": radial,
        "varout": numpy.array(interp_var).reshape((nrho, nphi)),
    }

    msg = f"Output variable has radial dimension {nrho} and azimuthal dimension {nphi}."
    logger.info(msg=msg)

    varout_obj = parser_interface.object_define()

    for varout_attr in varout_dict:
        value = parser_interface.dict_key_value(
            dict_in=varout_dict, key=varout_attr, no_split=True
        )
        varout_obj = parser_interface.object_setattr(
            object_in=varout_obj, key=varout_attr, value=value
        )

    return varout_obj


# ----


def interp_vertical(varin: numpy.array, zarr: numpy.array, levs: List) -> numpy.array:
    """
    Description
    -----------

    This method interpolates a 3-dimensional variable to specified
    vertical levels.

    Parameters
    ----------

    varin: array-type

        A Python array for the 3-dimensional variable to be
        interpolated.

    zarr: array-type

        A Python array for the vertical level type; this array must be
        of the same dimension as varin.

    levs: list

        A Python list of levels to which to interpolate; the units of
        this list must be identical to the units of the zarr array.

    Returns
    -------

    varout: array-type

        A Python array containing the 3-dimensional variable
        interpolated to the specified vertical levels.

    Raises
    ------

    TCDiagsError:

        * raised if an exception is encountered during the vertical
          interpolation.

    """

    # Interpolate the 3-dimensional variable specified upon input to
    # the specified vertical-type levels.
    try:
        varout = interplevel(varin, zarr, levs)

    except Exception as errmsg:
        msg = f"The vertical interpolation failed with error {errmsg}. Aborting!!!"
        raise TCDiagsError(msg=msg) from errmsg

    return varout
