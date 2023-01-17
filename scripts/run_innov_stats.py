# =========================================================================

# $$$ SCRIPT DOCUMENTATION BLOCK

# UFS-RNR-containers :: innov_stats/py/innov_stats/innov_stats.py

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

    InnovStatsError(msg)

        This is the base-class for all exceptions; it is a sub-class
        of Error.

Functions
---------

    main()

        This is the driver-level function to invoke the tasks within
        this script.

Author(s)
---------

    Henry R. Winterbottom; 16 January 2023

History
-------

    2023-01-16: Henry Winterbottom -- Initial implementation.

"""

# ----

import inspect
import os
import time

from confs.yaml_interface import YAML
from innov_stats import gsi_atmos
from innov_stats import soca_ice
from innov_stats import soca_ocean
from tools import fileio_interface
from tools import parser_interface
from utils.arguments_interface import Arguments
from utils.error_interface import Error, msg_except_handle
from utils.logger_interface import Logger

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


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
        self.basedir = os.getcwd()
        self.yaml_file = options_obj.yaml_file

        # Parse the experiment configuration file.
        yaml_dict = YAML().read_yaml(yaml_file=self.yaml_file)

        # Define the allowable applications; proceed accordingly.
        self.innov_stats_dict = {
            'gsi_atmos': gsi_atmos.GSIAtmos,
            'soca_ice': soca_ice.SOCAIce,
            'soca_ocean': soca_ocean.SOCAOcean}

        self.innov_stats_app = parser_interface.dict_key_value(
            dict_in=yaml_dict, key='innov_stats_app',
            force=True, no_split=True)
        if self.innov_stats_app is None:
            msg = ('The innovation statistics application could '
                   'not be determined from the user experiment '
                   'configuration. Aborting!!!')
            error(msg=msg)

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
        innov_stats = parser_interface.dict_key_value(
            dict_in=self.innov_stats_dict, key=self.innov_stats_app,
            force=True, no_split=True)
        if innov_stats is None:
            msg = ('The innovation statistics application {0} is not '
                   'supported. Aborting!!!'.format(self.innov_stats_app))
            error(msg=msg)

        app = innov_stats(yaml_file=self.yaml_file, basedir=self.basedir)
        app.run()

# ----


class InnovStatsError(Error):
    """
    Description
    -----------

    This is the base-class for all exceptions; it is a sub-class of
    Error.

    """

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

    # Define the schema attributes.
    cls_schema = {"yaml_file": str}

    # Collect the command line arguments.
    script_name = os.path.basename(__file__)
    start_time = time.time()
    msg = f'Beginning application {script_name}.'
    Logger().info(msg=msg)
    options_obj = Arguments().run(eval_schema=True, cls_schema=cls_schema)

    # Launch the task.
    task = InnovStats(options_obj=options_obj)
    task.run()

    stop_time = time.time()
    msg = f'Completed application {script_name}.'
    Logger().info(msg=msg)
    total_time = stop_time - start_time
    msg = f'Total Elapsed Time: {total_time} seconds.'
    Logger().info(msg=msg)

# ----


if __name__ == '__main__':
    main()
