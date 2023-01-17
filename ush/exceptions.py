# =========================================================================

# Module: ush/exceptions.py

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

    exceptions.py

Description
-----------

    This module loads the exceptions package.

Classes
-------

    InnovStatsError()

        This is the base-class for exceptions encountered within the
        ush/innov_stats module; it is a sub-class of Error.


Author(s)
---------

    Henry R. Winterbottom; 15 January 2023

History
-------

    2023-01-15: Henry Winterbottom -- Initial implementation.

"""

# ----

from utils.error_interface import Error

# ----

# Define all available attributes.
__all__ = [
    "InnovStatsError"
]

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


class InnovStatsError(Error):
    """
    Description
    -----------

    This is the base-class for exceptions encountered within the
    ush/innov_stats module; it is a sub-class of Error.

    """
