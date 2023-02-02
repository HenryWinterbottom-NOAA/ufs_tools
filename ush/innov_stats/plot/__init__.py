# =========================================================================

# Module: ush/innov_stats/plot/__init__.py

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

    This module contains the base-class and methods specific to
    innovation statistics applications.

Classes
-------

    InnovStats(options_obj, basedir)

        This is the base-class object for all innovation statistics
        applications.

Functions
---------

    error(msg)

        This function is the exception handler for the respective
        module.

Author(s)
---------

    Henry R. Winterbottom; 29 January 2023

History
-------

    2023-01-30: Henry Winterbottom -- Initial implementation.

"""

# ----

# pylint: disable=broad-except
# pylint: disable=too-many-instance-attributes
# pylint: disable=unused-argument

# ----

from dataclasses import dataclass

from confs.yaml_interface import YAML
from exceptions import InnovStatsPlotError
from tools import datetime_interface, fileio_interface, parser_interface
from utils import timestamp_interface
from utils.error_interface import msg_except_handle
from utils.logger_interface import Logger

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


@dataclass
class InnovStats:
    """
    Description
    -----------

    This is the base-class object for all innovation statistics
    applications.

    Parameters
    ----------

    options_obj: object

        A Python object containing the command line argument
        attributes.

    basedir: str

        A Python string specifying the path to the directory path in
        which the application is being executed.

    """

    def __init__(self, options_obj: object, basedir: str) -> None:
        """
        Description
        -----------

        Creates a new InnovStats object.

        """

        # Define the base-class attributes.
        self.logger = Logger()
        self.options_obj = options_obj
        self.basedir = basedir
        self.yaml_dict = YAML().read_yaml(yaml_file=self.options_obj.yaml_file)

        # Define the respective (i.e., supported) application
        # variables.
        self.gsi_atmos_variables = ["spechumid", "temperature", "uvwind"]
        self.soca_ice_variables = ["thickness"]
        self.soca_ocean_variables = ["salinity", "temperature"]

        # Define the diagnostics to be computed.
        self.stats_vars_list = ["bias", "count", "rmsd"]

        # Collect the respective innovations statistics application
        # attributes.
        self.regions_obj = parser_interface.object_define()
        self.get_regions()

        self.cycles_obj = parser_interface.object_define()
        self.get_cycles()

        self.expts_obj = parser_interface.object_define()
        self.get_experiments()

        self.vars_obj = parser_interface.object_define()
        self.get_variables()

    def get_cycles(self) -> None:
        """
        Description
        -----------

        This method collects the forecast cycle information from the
        user experiment configuration and defines the base-class
        attribute cycle_obj.

        Raises
        ------

        InnovStatsPlotError:

            * raised if the cycles attribute cannot be determined from
              the user experiment configuration.

            * raised if an attribute specifying the forecast cycle
              configuration cannot be determined from the user
              experiment configuration.

        """

        # Initialize the base-class cycles_obj attribute and collect
        # the analysis cycles information from the user experiment
        # configuration.
        cycles_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key="cycles", force=True, no_split=True
        )
        if cycles_dict is None:
            msg = (
                "The attribute cycles could not be determined from "
                "the user experiment configuration. Aborting!!!"
            )
            error(msg=msg)

        # Construct the base-class attribute cycles_obj using the
        # information specified in the user experiment configuration.
        cycle_attrs_list = ["cycle_interval", "cycle_start", "cycle_stop"]
        for cycle_attr in cycle_attrs_list:
            value = parser_interface.dict_key_value(
                dict_in=cycles_dict, key=cycle_attr, force=True, no_split=True
            )

            if value is None:
                msg = (
                    f"The cycles attribute {cycle_attr} could not be determined "
                    "from the user experiment configuration. "
                    "Aborting!!!"
                )
                error(msg=msg)

            self.cycles_obj = parser_interface.object_setattr(
                object_in=self.cycles_obj, key=cycle_attr, value=value
            )

    def get_experiments(self) -> None:
        """
        Description
        -----------

        This method collects the information for the experiments for
        which to computed the innovation statistics from the user
        experiment configuration and defines the base-class attribute
        expts_obj.

        Raises
        ------

        InnovStatsPlotError:

            * raised if the user experiment configuration does not
              contain the attribute experiments.

        """

        # Initialize the base-class expts_obj attribute and collect
        # the experiment information from the user experiment
        # configuration.

        expts_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key="experiments", force=True, no_split=True
        )
        if expts_dict is None:
            msg = (
                "The attribute experiments could not be determined "
                "from the user experiment configuration. Aborting!!!"
            )
            error(msg=msg)

        # Construct the base-class attribute expts_obj using the
        # information specified in the user experiment configuration.
        for expt in expts_dict.keys():
            value = parser_interface.dict_key_value(
                dict_in=expts_dict, key=expt, no_split=True
            )
            self.expts_obj = parser_interface.object_setattr(
                object_in=self.expts_obj, key=expt, value=value
            )

    def get_regions(self) -> None:
        """
        Description
        -----------

        This method collects the information for the analysis regions
        to be plotted from the user experiment configuration and
        defines the base-class attribute regions_obj.

        Raises
        ------

        InnovStatsPlotError:

            * raised if the user experiment configuration does not
              contain the attribute regions.

        """

        # Initialize the base-class regions_obj attribute and collect
        # the regions information from the user experiment
        # configuration.
        regions_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key="regions", force=True, no_split=True
        )

        # Construct the base-class attribute regions_obj accordingly.
        if regions_dict is None:
            self.regions_obj = None
            msg = (
                "The attribute regions is not specified or could not "
                "be determined from the user experiment configuration."
            )
            self.logger.warn(msg=msg)

        if regions_dict is not None:
            for region in regions_dict.keys():
                value = parser_interface.dict_key_value(
                    dict_in=regions_dict, key=region, no_split=True
                )

                self.regions_obj = parser_interface.object_setattr(
                    object_in=self.regions_obj, key=region, value=value
                )

    def get_times(self, expt_name: str) -> list:
        """
        Description
        -----------

        This method collects and returns the valid analysis time
        timestamp values valid for the respective experiment name
        provided upon entry.

        Parameters
        ----------

        expt_name: str

            A Python string specifying an experiment name defined
            within the user experiment configuration.

        Returns
        -------

        timestamp_list: list

            A Python list containing the valid timestamp values for
            the respective experiment.

        Raises
        ------

        InnovStatsPlotError:

            * raised if the attributes for the specified experiment
              name cannot be determined from the user experiment
              configuration.

            * raised if an exception is encountered while defining the
              analysis starting and stopping timestamp values.

        """

        # Collect the experiment attributes from the user experiment
        # configuration using the experiment name (expt_name) provided
        # upon entry.
        expt_dict = parser_interface.object_getattr(
            object_in=self.expts_obj, key=expt_name, force=True
        )
        if expt_dict is None:
            msg = (
                f"The attributes for experiment {expt_name} could not "
                "be determined from the user experiment "
                "configuration. Aborting!!!"
            )
            error(msg=msg)

        # Define the start and stop analysis times using the
        # attributes within the base-class attribute cycles_obj.
        try:
            timestamp = datetime_interface.datestrupdate(
                datestr=self.cycles_obj.cycle_start,
                in_frmttyp=timestamp_interface.GENERAL,
                out_frmttyp=timestamp_interface.GLOBAL,
                offset_seconds=self.cycles_obj.cycle_interval,
            )

            last_timestamp = datetime_interface.datestrupdate(
                datestr=self.cycles_obj.cycle_stop,
                in_frmttyp=timestamp_interface.GENERAL,
                out_frmttyp=timestamp_interface.GLOBAL,
            )

        except Exception as errmsg:
            msg = (
                "Defining the initial and final timestamps failed with error "
                f"{errmsg}. Aborting!!!"
            )
            error(msg=msg)

        # Create a list of valid timestamps for the respective
        # experiment name provided upon entry.
        timestamp_list = []

        msg = f"The valid timestamps for experiment {expt_name} are the following:\n\n"
        while timestamp <= last_timestamp:
            timestamp_list.append(timestamp)
            timestamp = datetime_interface.datestrupdate(
                datestr=timestamp,
                in_frmttyp=timestamp_interface.GLOBAL,
                out_frmttyp=timestamp_interface.GLOBAL,
                offset_seconds=self.cycles_obj.cycle_interval,
            )

            msg = msg + f"{timestamp}\n"
        self.logger.warn(msg=msg)

        return timestamp_list

    def get_variables(self) -> None:
        """
        Description
        -----------

        This method collects the information for the analysis
        variables to be plotted from the user experiment configuration
        and defines the base-class attribute vars_obj.

        Raises
        ------

        InnovStatsPlotError:

            * raised if the user experiment configuration does not
              contain the attribute variables.

        """

        # Initialize the base-class vars_obj attribute and collect the
        # variables information from the user experiment
        # configuration.

        vars_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key="variables", force=True, no_split=True
        )
        if vars_dict is None:
            msg = (
                "The attribute variables could not be determined "
                "from the user experiment configuration. Aborting!!!"
            )
            error(msg=msg)

        # Construct the base-class attribute expts_obj using the
        # information specified in the user experiment configuration.
        for var in vars_dict.keys():
            value = parser_interface.dict_key_value(
                dict_in=vars_dict, key=var, no_split=True
            )

            self.vars_obj = parser_interface.object_setattr(
                object_in=self.vars_obj, key=var, value=value
            )


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
