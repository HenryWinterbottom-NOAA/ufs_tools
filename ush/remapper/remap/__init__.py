# =========================================================================

# Module: ush/remapper/remap/__init__.py

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

    This module loads the remap package.

Classes
-------

    Remap()

        This is the base-class object for all remapping tasks and
        applications.

Requirements
------------

- ufs_pytils; https://github.com/HenryWinterbottom-NOAA/ufs_pyutils

Author(s)
---------

    Henry R. Winterbottom; 12 February 2023

History
-------

    2023-02-12: Henry Winterbottom -- Initial implementation.

"""

# ----

from dataclasses import dataclass

from exceptions import RemapperError
from utils.logger_interface import Logger

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


@dataclass
class Remap:
    """
    Description
    -----------

    This is the base-class object for all remapping tasks and
    applications.

    """

    def __init__(self):
        """
        Description
        -----------

        Creates a new Remap object.

        """

        # Define the base-class attributes.
        self.logger = Logger()
