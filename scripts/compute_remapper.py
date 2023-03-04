# =========================================================================

# Script: scripts/compute_remapper.py

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
Script
------

    compute_remapper.py

Description
-----------

    This script is the driver script for remapper application.

Classes
-------

    ComputeRemapper(options_obj)

        This is the base-class object for all remapper applications.

Functions
---------

    main()

        This is the driver-level function to invoke the tasks within
        this script.

Usage
-----

    user@host:$ python compute_remapper.py --yaml_file /path/to/yaml_file

Parameters
----------

    yaml_file: str

        A Python string specifying the path to the YAML-formatted
        configuration file for the respective remapper application.

        --yaml_file /path/to/yaml/file or -yaml_file /path/to/yaml/file

Requirements
------------

- ufs_pytils; https://github.com/HenryWinterbottom-NOAA/ufs_pyutils

Author(s)
---------

    Henry R. Winterbottom; 03 March 2023

History
-------

    2023-03-03: Henry Winterbottom -- Initial implementation.

"""

# ----

import os
import time
from dataclasses import dataclass

from remapper import Remapper
from utils.arguments_interface import Arguments
from utils.logger_interface import Logger

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----

# Specify whether to evaluate the format for the respective parameter
# values.
EVAL_SCHEMA = True

# Define the schema attributes.
CLS_SCHEMA = {"yaml_file": str}

# ----


@dataclass
class ComputeRemapper:
    """
    Description
    -----------

    This is the base-class object for all remapper applications.

    Parameters
    ----------

    options_obj: object

        A Python object containing the command line argument
        attributes.

    """

    def __init__(self, options_obj: object):
        """
        Description
        -----------

        Creates a new ComputeRemapper object.

        """

        # Define the base-class attributes.
        self.options_obj = options_obj
        self.remapper = Remapper(options_obj=self.options_obj)

    def run(self) -> None:
        """
        Description
        -----------

        This method performs the following tasks:

        (1) Executes the remapper application to interpolate (i.e.,
            remap) variable fields between grid projections and/or
            grid-spacing resolutions.

        """

        # Remap the variable fields accordingly.
        self.remapper.run()


# ----


def main():
    """
    Description
    -----------

    This is the driver-level function to invoke the tasks within this
    script.

    """

    # Collect the command line arguments.
    script_name = os.path.basename(__file__)
    start_time = time.time()
    msg = f"Beginning application {script_name}."
    Logger().info(msg=msg)
    options_obj = Arguments().run(eval_schema=EVAL_SCHEMA, cls_schema=CLS_SCHEMA)

    # Launch the task.
    task = ComputeRemapper(options_obj=options_obj)
    task.run()

    stop_time = time.time()
    msg = f"Completed application {script_name}."
    Logger().info(msg=msg)
    total_time = stop_time - start_time
    msg = f"Total Elapsed Time: {total_time} seconds."
    Logger().info(msg=msg)


# ----


if __name__ == "__main__":
    main()
