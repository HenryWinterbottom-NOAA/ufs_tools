# =========================================================================

# Script: wrappers/wrapper_esmf_remap.py

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

    wrapper_esmf_remap.py

Description
-----------

    This script is the wrapper script for Earth System Modeling
    Framework (ESMF) remapping attribute computations and definitions
    provided by scripts/compute_esmf_remap.py

Classes
-------

    WrapperESMFRemap(options_obj)

        This is the base-class object for the Earth System Modeling
        Framework (ESMF) remapping application; it acts as wrapper and
        is governed by the contents of the command line attribute
        yaml_file.

Functions
---------

    main()

        This is the driver-level function to invoke the tasks within
        this script.

Usage
-----

    user@host: $ python wrapper_esmf_remap.py --yaml_file --yaml_template --script_path

Parameters
----------

    yaml_file: str

        A Python string specifying the path to the YAML-formatted
        configuration file for the ESMF remap application.

        --yaml_file / path/to/yaml/file or -yaml_file / path/to/yaml/file

    yaml_file: str

        A Python string specifying the path to the YAML-formatted
        template file; this is typically found beneath
        parm/compute_remap.

        --yaml_template /path/to/yaml_template  or -yaml_template /path/to/yaml_template

    script_path: str

        A Python string specifying the path to the
        compute_esmf_remap.py script; this is typically found beneath
        scripts/.

        --script_path /path/to/compute_esmf_remap.py or -script_path /path/to/compute_esmf_remap.py

Requirements
------------

- ufs_pytils; https: // github.com/HenryWinterbottom-NOAA/ufs_pyutils

Author(s)
---------

    Henry R. Winterbottom; 10 February 2023

History
-------

    2023-02-10: Henry Winterbottom - - Initial implementation.

"""

# ----

# pylint: disable=too-many-instance-attributes

# ----

import os
import time
from dataclasses import dataclass
from typing import Dict

from confs.yaml_interface import YAML
from exceptions import ESMFRemapError
from execute import subprocess_interface
from tools import parser_interface, system_interface
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
CLS_SCHEMA = {"script_path": str, "yaml_file": str, "yaml_template": str}

# ----


@dataclass
class WrapperESMFRemap:
    """
    Description
    -----------

    This is the base-class object for the Earth System Modeling
    Framework (ESMF) remapping application; it acts as wrapper and is
    governed by the contents of the command line attribute yaml_file.

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

        Creates a new WrapperESMFRemap object.

        """

        # Define the base-class attributes.
        self.logger = Logger()
        self.options_obj = options_obj
        self.yaml_template = self.options_obj.yaml_template
        self.yaml_file = self.options_obj.yaml_file
        self.yaml_dict = YAML().read_yaml(yaml_file=self.yaml_file)

        # Define the application configuration attributes.
        self.grid_type_dict = {"destination": "DST_", "source": "SRC_"}

        # Define the attributes for the ESMF remapping application
        # script.
        [self.exe, self.job_type] = [
            system_interface.get_app_path(app="python"), "app"]

    def __defenv__(self, yaml_key: str, yaml_dict: Dict) -> None:
        """
        Description
        -----------

        This method prepares the run-time environment for the ESMF
        remapping application.

        Parameters
        ----------

        yaml_key: str

            A Python string specifying the YAML-formatted file key
            containing the ESMF remap application attributes.

        yaml_dict: dict

            A Python dictionary containing the YAML-formatted file
            attributes for the respective ESMF remapping application.

        Raises
        ------

        ESMFRemapError:

            * raised if the a mandatory attribute cannot be determined
              from the ESMF remapping application attributes.

            * raised if either the destination or source grid
              attributes can be determined from the experiment
              configuration for the respective ESMF remapping
              application.

        """

        # Define the run-time environment interpolation type and
        # output netCDF-formatted file.
        for item in ["interp_type", "output_netCDF"]:

            value = parser_interface.dict_key_value(
                dict_in=yaml_dict, key=item, force=True, no_split=True
            )
            if value is None:
                msg = (
                    f"The mandatory remapping attribute {item} could not "
                    f"be determined for application {yaml_key} in "
                    f"YAML-formatted configuration file {self.yaml_file}. "
                    "Aborting!!!"
                )
                raise ESMFRemapError(msg=msg)

            parser_interface.enviro_set(
                envvar=f"{item.upper()}", value=str(value))

        # Update the run-time environment for the respective remapping
        # application.
        for grid_type in self.grid_type_dict:

            # Collect the attributes for the respective grid type;
            # proceed accordingly.
            grid_dict = parser_interface.dict_key_value(
                dict_in=yaml_dict, key=grid_type, force=True, no_split=True
            )
            if grid_dict is None:
                msg = (
                    f"The grid attributes for grid type {grid_type} for "
                    f"application {yaml_key} could not be determined "
                    f"from YAML-formatted configuration file {self.yaml_file}. "
                    "Aborting!!!"
                )
                raise ESMFRemapError(msg=msg)

            # Define the prefix string for the respective grid type
            # run-time environment variables.
            prefix = parser_interface.dict_key_value(
                dict_in=self.grid_type_dict, key=grid_type, no_split=True
            )

            # Define the run-time environment for the respective grid
            # type.
            for item in grid_dict:

                # Define the run-time environment.
                envvar = f"{prefix}{item.upper()}"
                value = parser_interface.dict_key_value(
                    dict_in=grid_dict, key=item, no_split=True
                )

                parser_interface.enviro_set(envvar=envvar, value=str(value))

    def build_inputs(self, yaml_key: str) -> None:
        """
        Description
        -----------

        This method builds the inputs for the ESMF remap application;
        this includes defining the run-time environment and providing
        the respective application configuration.

        Parameters
        ----------

        yaml_key: str

            A Python string specifying the YAML-formatted file key
            containing the ESMF remap application attributes.

        Raises
        ------

        ESMFRemapError:

            * raised if specified ESMF remapping application
              attributes cannot be determined from the experiment
              configuration file.

        """

        # Collect the attributes for the respective remapping
        # application from the configuration file.
        msg = f"Configuration remapping application {yaml_key}."
        self.logger.info(msg=msg)

        yaml_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key=yaml_key, force=True, no_split=True
        )
        if yaml_dict is None:
            msg = (
                f"The attributes for ESMF remap application {yaml_key} "
                "could not be determined from YAML-formatted configuration "
                f"file {self.yaml_file}. Aborting!!!"
            )
            raise ESMFRemapError(msg=msg)

        # Define the run-time environment and write the respective
        # configuration to standard out.
        self.__defenv__(yaml_key=yaml_key, yaml_dict=yaml_dict)
        YAML().dict_to_yaml(yaml_dict=yaml_dict, level="warn", nspace=2)

    def esmf_remap(self, yaml_key: str) -> None:
        """
        Description
        -----------

        This method executes the ESMF remapping application for the
        configuration defined by the yaml_key parameter value upon
        entry.

        Parameters
        ----------

        yaml_key: str

            A Python string specifying the YAML-formatted file key
            containing the ESMF remap application attributes.

        """

        # Run the remapping application.
        msg = f"Launching remapping application {yaml_key}."
        self.logger.info(msg=msg)

        args = [f"{self.options_obj.script_path}", "--yaml_file", self.yaml_template]
        errlog = f"{yaml_key}.err"
        outlog = f"{yaml_key}.out"

        subprocess_interface.run(
            exe=self.exe,
            job_type=self.job_type,
            args=args,
            errlog=errlog,
            outlog=outlog,
        )

        msg = f"Completed remapping application {yaml_key}."
        self.logger.info(msg=msg)

    def run(self) -> None:
        """
        Description
        -----------

        This method performs the following tasks:

        (1) Looping through each ESMF remapping application, defines
            the run-time environment.

        (2) Looping through each ESMF remapping application, executes
            the ESMF remap application for the respective remapping
            configuration.

        """

        # Loop through the YAML-keys; the YAML-keys define the
        # different ESMF remappings to generate.
        for yaml_key in self.yaml_dict:

            # Define the run-time environment.
            self.build_inputs(yaml_key=yaml_key)

            # Launch the respective remapping configuration.
            self.esmf_remap(yaml_key=yaml_key)


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
    task = WrapperESMFRemap(options_obj=options_obj)
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
