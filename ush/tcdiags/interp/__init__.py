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

import numpy
from wrf import interplevel

from typing import List

from exceptions import TCDiagsError
from utils.logger_interface import Logger

# ----

logger = Logger()

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

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
