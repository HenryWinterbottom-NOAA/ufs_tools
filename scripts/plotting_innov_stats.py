# =========================================================================

# Script: scripts/compute_innov_stats.py

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

    innov_stats.py

Description
-----------

    This script contains classes and functions to plot innovation
    statistic diagnostics as specified within the user experiment
    configuration.

Classes
-------

    InnovStats()

        This is the base-class object for all supported innovation
        statistic diagnostics.

Functions
---------

    error(msg)

        This function is the exception handler for the respective
        script.

    main()

        This is the driver-level function to invoke the tasks within
        this script.

Author(s)
---------

    Henry R. Winterbottom; 29 January 2023

History
-------

    2023-01-29: Henry Winterbottom -- Initial implementation.

"""

# ----

import os
import time
from dataclasses import dataclass

from tools import fileio_interface
from tools import parser_interface

from exceptions import InnovStatsPlotError

from innov_stats.plot.apps import gsi_atmos
#from innov_stats.plot.apps import soca_ice
#from innov_stats.plot.apps import soca_ocean

from confs.yaml_interface import YAML

from utils.arguments_interface import Arguments
from utils.error_interface import msg_except_handle
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
CLS_SCHEMA = {
    "yaml_file": str,
    "app": str,
    "plot": str
}

# Define the allowable applications.
INNOV_STATS_DICT = {'gsi_atmos_omfprofile': gsi_atmos.GSIAtmosOMFProfile,
                    'gsi_atmos_omftimeseries': gsi_atmos.GSIAtmosOMFTimeseries}
#                    'soca_ice_omftimeseries': soca_ice.SOCAIceOMFTimeseries,
#                    'soca_ocean_omfprofile': soca_ocean.SOCAOceanOMFProfile,
#                    'soca_ocean_omftimeseries': soca_ocean.SOCAOceanOMFTimeseries
#                    }

# ----


class InnovStats:
    """
    Description
    -----------

    This is the base-class object for all supported innovation
    statistic diagnostics.


    """

    def __init__(self, options_obj: object):
        """
        Description
        -----------

        Creates a new InnovStats object.

        """

        # Define the base-class attributes.
        self.options_obj = options_obj
        self.basedir = os.getcwd()

        # Parse the user experiment configuration and proceed
        # accordingly.
        self.yaml_file = self.options_obj.yaml_file
        self.yaml_dict = YAML().read_yaml(yaml_file=self.yaml_file)

        # Define the plotting application.
        plot_app = f"{self.options_obj.app}_{self.options_obj.plot}"
        self.innov_stats_plot = parser_interface.dict_key_value(
            dict_in=INNOV_STATS_DICT, key=plot_app, no_split=True)

    def run(self):
        """
        Description
        -----------

        This method performs the following tasks:

        (1) Determines and defines the innovation statistics
            application based on the user experiment configuration
            attributes.

        (2) Launches the respective innovation statistics application.

        """

        # Launch the respective innovation statistics plotting
        # application.
        app = self.innov_stats_plot(
            options_obj=self.options_obj, basedir=self.basedir)
        app.run()

# ----


@msg_except_handle(InnovStatsPlotError)
def error(msg: str) -> None:
    """
    Description
    -----------

    This function is the exception handler for the respective module.

    Parameters
    ----------

    msg: str

        A Python string containing a message to accompany the
        exception.

    """

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
    task = InnovStats(options_obj=options_obj)
    task.run()

    stop_time = time.time()
    msg = f"Completed application {script_name}."
    Logger().info(msg=msg)
    total_time = stop_time - start_time
    msg = f"Total Elapsed Time: {total_time} seconds."
    Logger().info(msg=msg)

# ----


if __name__ == '__main__':
    main()
