# =========================================================================

# Script: scripts/run_innov_stats.py

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

    run_innov_stats.py

Description
-----------

    This script is the driver script for the innovation statistics
    diagnostics.

Classes
-------

    InnovStats()

        This is the base-class for all supported innovation statistic
        diagnostics.

Functions
---------

    error(msg)

        This function is the exception handler for the respective
        script.

    main()

        This is the driver-level function to invoke the tasks within
        this script.

Usage
-----

    user@host:$ python run_innov_stats --<yaml_file> --<cycle> --<app> \
                    [--expt_name]

Parameters
----------

    yaml_file: str

        A Python string specifying the path to the YAML-formatted
        configuration file for the respective application type (see
        app)

        Enter the parameter value as:

        --yaml_file /path/to/yaml/file or -yaml_file /path/to/yaml/file

    cycle: str

        A Python string specifying the respective forecast cycle; this
        string must be formatted as %Y%m%d%H%M%S assuming the POSIX
        convention; enter the parameter value as follows for a
        forecast cycle beginning 0000 UTC 01 January 2000:

        --cycle 20000101000000 or -cycle 20000101000000

    app: str

        A Python string specifying the application from which to
        compute the innovation statistics; the following are the
        allowable application arguments:

        gsi_atmos: Gridpoint Statistical Interpolation (GSI)
                   atmosphere observation innovation statistics and
                   timeseries.

        soca_ice: Sea-ice and Ocean Coupled Analysis (SOCA) ice
                  observation innovation timeseries.

        soca_ocean: Sea-ice and Ocean Coupled Analysis (SOCA) ocean
                    observation innovation statistics and timeseries.

        Enter the parameter value as:

        --app gsi_atmos|soca_ice|soca_ocean or

        -app gsi_atmos|soca_ice|soca_ocean

    expt_name: str, optional

        A Python string specifying an (unique) name for the respective
        experiment; for an experiment named SPAM, enter the parameter
        value as follows.

        --expt_name SPAM or -expt_name SPAM


Author(s)
---------

    Henry R. Winterbottom; 16 January 2023

History
-------

    2023-01-16: Henry Winterbottom -- Initial implementation.

"""

# ----

# pylint: disable=broad-except
# pylint: disable=unused-argument

# ----

import os
import time
from dataclasses import dataclass

from exceptions import InnovStatsError
from innov_stats import gsi_atmos, soca_ice, soca_ocean
from schema import Optional, Or
from tools import datetime_interface, parser_interface
from utils import timestamp_interface
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
    "cycle": Or(str, int),
    "app": str,
    Optional("expt_name"): Or(str, None),
}

# Define the allowable applications.
INNOV_STATS_DICT = {
    "gsi_atmos": gsi_atmos.GSIAtmos,
    "soca_ice": soca_ice.SOCAIce,
    "soca_ocean": soca_ocean.SOCAOcean,
}

# ----


@dataclass
class InnovStats:
    """
    Description
    -----------

    This is the base-class for all supported innovation statistic
    diagnostics.

    Parameters
    ----------

    options_obj: object

        A Python object containing the command line argument
        attributes.

    Raises
    ------

    InnovStatsError:

        * raised if the user experiment configuration does not specify
          the supported innovation statistics application.

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
        self.cycle = self.options_obj.cycle
        self.app = self.options_obj.app

        # Check that the analysis cycle timestamp format is valid;
        # proceed accordingly.
        try:
            datetime_interface.datestrupdate(
                datestr=self.cycle,
                in_frmttyp=timestamp_interface.GLOBAL,
                out_frmttyp=timestamp_interface.GLOBAL,
            )

        except Exception:
            msg = (
                "The formatting of the command line variable cycle is not of "
                f"format {timestamp_interface.GLOBAL}; received {self.cycle} "
                "upon entry. Aborting!!!"
            )
            error(msg=msg)

        # Define the base-class object for the respective appliction
        # for which to compute the innovation statistics; proceed
        # accordingly.
        innov_stats = parser_interface.dict_key_value(
            dict_in=INNOV_STATS_DICT, key=self.app.lower(), force=True, no_split=True
        )
        if innov_stats is None:
            msg = (
                f"The innovation statistics application {self.app} is not "
                "supported. Aborting!!!"
            )
            error(msg=msg)

        self.innov_stats = innov_stats(
            options_obj=self.options_obj, basedir=self.basedir
        )

    def run(self):
        """
        Description
        -----------

        This method performs the following tasks.

        (1) Defines the base-class object corresponding to the
            respective innovation statistics application.

        (2) Executes the respective innovation statistics application.

        Raises
        ------

        InnovStatsError:

            * raised if the innovation statistic application specified
              within the user experiment configuration is not
              supported.

        """

        # Compute the respective innovation statistics.
        self.innov_stats.run()


# ----


@msg_except_handle(InnovStatsError)
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


if __name__ == "__main__":
    main()
