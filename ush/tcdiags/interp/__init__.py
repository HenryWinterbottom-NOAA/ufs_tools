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



"""

# ----

from math import asin, cos, radians, sin
from typing import List

# from scipy.spatial import KDTree

import numpy
from tcdiags.geomets import haversine
from exceptions import TCDiagsError
from utils.logger_interface import Logger
# from wrf import interplevel

# ----

logger = Logger()

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


def interp_ll2ra(varin: numpy.array, lats: numpy.array,
                 lons: numpy.array, lon_0: float, lat_0: float,
                 max_radius: float, drho: float, dphi: float) -> numpy.array:
    """ """

    # Initialize the coordinate arrays.
    varin = numpy.ravel(varin)
    lats = numpy.ravel(lats)
    lons = numpy.ravel(lons)

    # Compute the radial distance relative to the specified
    # geographical coordinate location.
    fix = (lat_0, lon_0)
    rho = numpy.array([haversine(fix, (lats[idx], lons[idx]))
                       for idx in range(len(lats))])

    xx = numpy.array([haversine(fix, (lat_0, lons[idx]))
                      for idx in range(len(lats))])
    xx = numpy.where(lons < lon_0, -1.0*xx, xx)

    yy = numpy.array([haversine(fix, (lats[idx], lon_0))
                      for idx in range(len(lats))])
    yy = numpy.where(lats < lat_0, -1.0*yy, yy)

    phi = numpy.arctan2(yy, xx)

    # Define the grid attributes relative to the specified
    # geographical coordinate location.
    xlocs = numpy.where(numpy.logical_and(rho >= 0.0, rho <= max_radius))[0]
    xlats = numpy.array([lats[idx] for idx in xlocs])
    xlons = numpy.array([lons[idx] for idx in xlocs])
    xrho = numpy.array([rho[idx] for idx in xlocs])
    xphi = numpy.array([phi[idx] for idx in xlocs])
    xvar = numpy.array([varin[idx] for idx in xlocs])


#    xx = numpy.array([haversine(fix, (lat_0, lons[idx]))
#                      for idx in xlocs])
#    xx = numpy.where(xlons < lon_0, -1.0*xx, xx)

#    yy = numpy.array([haversine(fix, (lats[idx], lon_0))
#                      for idx in xlocs])
#    yy = numpy.where(xlats < lat_0, -1.0*yy, yy)

    # Compute the azimuthal angle relative to the specified
    # geographical coordinate location.
#    xphi = numpy.arctan2(yy, xx)
    dphi = numpy.radians(dphi)

    # Interpolate the Cartesian grid to the defined polar coordinate
    # grid.
    radius = numpy.arange(0.0, (max_radius + drho), drho)
    azimuth = numpy.arange(-1.0*numpy.pi, (numpy.pi+dphi), dphi)

    var = []

    for ridx in range(len(radius)):
        radii = radius[ridx]

        for aidx in range(len(azimuth)):
            theta = azimuth[aidx]
            array = []

            idxs = numpy.ndarray.tolist(numpy.where(rho >= radii)[0])
            for idx in idxs:

                if rho[idx] < (radii + drho):
                    if (phi[idx] >= theta) and (phi[idx] < (theta+dphi)):
                        array.append(varin[idx])

            if len(array) == 0:
                var.append(numpy.nan)

            else:
                var.append(numpy.nanmean(array))

    var = numpy.array(var)

    check = numpy.logical_not(numpy.isnan(var))

    xp = check.ravel().nonzero()[0]
    fp = var[numpy.logical_not(numpy.isnan(var))]

    x = numpy.isnan(var).ravel().nonzero()[0]

    var[numpy.isnan(var)] = numpy.interp(x, xp, fp)

    (nrho, nphi) = [len(radius), len(azimuth)]

    varout = numpy.array(var).reshape((nrho, nphi))

    print(varout.min(), varout.max())

    return (radius, azimuth, varout)


# ----


def interp_vertical(varin: numpy.array, pres: numpy.array,
                    levs: List) -> numpy.array:
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
        raise TCDiagsError(msg=msg)

    return varout
