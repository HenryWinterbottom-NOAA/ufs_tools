# =========================================================================

# Module: ush/tcdiags/atmos/heights.py

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

    heights.py

Description
-----------

    This module contains various height profile computational methods.

Functions
---------

    height_from_pressure(inputs_obj)

        This function computes the geometric height profile from the
        pressure profile array.

Requirements
------------

- metpy; https://unidata.github.io/MetPy/latest/index.html

- ufs_pytils; https://github.com/HenryWinterbottom-NOAA/ufs_pyutils

Author(s)
---------

    Henry R. Winterbottom; 09 March 2023

History
-------

    2023-03-09: Henry Winterbottom -- Initial implementation.

"""

# ----

from metpy.calc import pressure_to_height_std
from metpy.units import units
from tools import parser_interface
from utils.logger_interface import Logger

# ----

# Define all available functions.
__all__ = ["height_from_pressure"]

# ----

logger = Logger()

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


def height_from_pressure(inputs_obj: object) -> object:
    """
    Description
    -----------

    This function computes the geometric height profile from the
    pressure profile array.

    Parameters
    ----------

    inputs_obj: object

        A Python object containing, at minimum, the pressure profile
        from which the geometric heights will be computed.

    Returns
    -------

    inputs_obj: object

        A Python object updated to contain the geometric heights
        profile.

    """

    # Compute the geometric height profile using the pressure profile.
    msg = (
        "Computing the geometric height array of dimension " f"{inputs_obj.pres.shape}."
    )
    logger.info(msg=msg)

    hght = units.Quantity(pressure_to_height_std(pressure=inputs_obj.pres), "meter")

    # Update the input variable object.
    inputs_obj = parser_interface.object_setattr(
        object_in=inputs_obj, key="hght", value=hght
    )

    return inputs_obj
