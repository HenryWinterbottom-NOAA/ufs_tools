# =========================================================================

# Module: ush/innov_stats/plot/omf_profile.py

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

    omf_profile.py

Description
-----------


    This module contains the base-class object and methods specific to
    temporal mean innovation statistic profile applications.

Classes
-------

    OMFProfile(options_obj, basedir)

         This is the base-class object for profile-type innovation
         statistics calculations; it is a sub-class of InnovStats.

Author(s)
---------

    Henry R. Winterbottom; 30 January 2023

History
-------

    2023-01-30: Henry Winterbottom -- Initial implementation.

"""

# ----

# pylint: disable=too-many-arguments
# pylint: disable=too-many-branches
# pylint: disable=too-many-locals

# ----

import os

import numpy
from ioapps import netcdf4_interface
from tools import datetime_interface, fileio_interface, parser_interface
from utils import timestamp_interface

from innov_stats.plot import InnovStats, error
from innov_stats.plot.plots import Plots

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


class OMFProfile(InnovStats):
    """
    Description
    -----------

    This is the base-class object for profile-type innovation
    statistics calculations; it is a sub-class of InnovStats.

    Parameters
    ----------

    options_obj: object

        A Python object containing the command line argument
        attributes.

    basedir: str

        A Python string specifying the path to the directory path in
        which the application is being executed.

    Keywords
    --------

    is_gsi_atmos: bool, optional

        A Python boolean valued variables specifying whether the
        methods relative to the GSI atmosphere application are to be
        implemented.

    is_soca_ice: bool, optional

        A Python boolean valued variables specifying whether the
        methods relative to the SOCA ice application are to be
        implemented.

    is_soca_ocean: bool, optional

        A Python boolean valued variables specifying whether the
        methods relative to the SOCA ocean application are to be
        implemented.

    Raises
    ------

    InnovStatsPlotError:

        * raised if the argument(s) pertaining the specified
          innovation statistics method are not valid.

    """

    def __init__(
        self,
        options_obj: object,
        basedir: str,
        is_gsi_atmos: bool = False,
        is_soca_ice: bool = False,
        is_soca_ocean: bool = False,
    ):
        """
        Description
        -----------

        Creates a new OMFProfile object.

        """

        # Define the base-class attributes.
        super().__init__(options_obj=options_obj, basedir=basedir)
        apps_dict = {
            is_gsi_atmos: self.gsi_atmos_variables,
            is_soca_ice: self.soca_ice_variables,
            is_soca_ocean: self.soca_ocean_variables,
        }

        # Check that the constructor arguments are valid; proceed
        # accordingly.
        if not any(apps_dict.keys()):
            msg = (
                "The supported innovation statistics application "
                "has not been specified. Aborting!!!"
            )
            error(msg=msg)

        for (app, _) in apps_dict.items():
            if app:
                self.variables = parser_interface.dict_key_value(
                    dict_in=apps_dict, key=app, no_split=True
                )

        msg = (
            "The analysis variable(s) for which to collect the "
            f"innovation statistics is (are) {self.variables}."
        )
        self.logger.info(msg=msg)

        if is_gsi_atmos:
            self.vars_list = self.gsi_atmos_variables
            self.vcoord_name = "plevs"

        if is_soca_ice:
            self.vars_list = self.soca_ice_variable

        if is_soca_ocean:
            self.vars_list = self.soca_ocean_variables
            self.vcoord_name = "depth"

        # Define the plotter object.
        self.plotter = Plots(options_obj=options_obj,
                             stats_list=self.stats_vars_list,
                             vars_list=self.vars_list)

    def compute(self, filelist_obj_dict: dict) -> object:
        """
        Description
        -----------

        This method computes the temporal mean analysis variable
        innovation statistics as a function of user specified
        experiment(s) and region(s) of interest.

        Parameters
        ----------

        filelist_obj_dict: dict

            A Python dictionary containing the file attributes for the
            respective experiments and analysis variables.

        Returns
        -------

        innov_stats_obj: object

            A Python object containing the temporal mean analysis
            variable innovation statistics as a function of user
            specified experiment(s) and region(s) of interest.

        Raises
        ------

        InnovStatsPlotError:

            * raised if the dictionary containing file attributes for
              the respective experiments and analysis variables upon
              entry does not contain the attributes for a specified
              experiment.

            * raised if the datapath attribute cannot be determined
              for a specified experiment.

            * raised if the innvation statistic file names for a
              specified experiment and analysis variable canno be
              determined from the dictionary containing file
              attributes for the respective experiments and analysis
              variables upon entry.

        """

        # Define the respective experiments and regions.
        expts_list = list(filelist_obj_dict.keys())

        # Loop through the respective experiments and proceed with the
        # temporal mean innovation statistic computation as a function
        # region and variable.
        innov_stats_obj = parser_interface.object_define()
        for expt in expts_list:
            innov_stats_obj = parser_interface.object_setattr(
                object_in=innov_stats_obj, key=expt, value=None
            )

            expt_dict = parser_interface.object_getattr(
                object_in=self.expts_obj, key=expt, force=True
            )

            if expt_dict is None:
                msg = (
                    "The input filelist dictionary (filelist_obj_dict) "
                    f"does not contain the attributes for experiment {expt} "
                    "upon entry. Aborting!!!"
                )
                error(msg=msg)

            datapath=parser_interface.dict_key_value(
                dict_in=expt_dict, key="datapath", force=True, no_split=True
            )

            if datapath is None:
                msg=(
                    f"The datapath attribute for experiment {expt} could "
                    "not be determined from the user experiment "
                    "configuration. Aborting!!!."
                )
                error(msg=msg)

            # Loop through the respective analysis variables and
            # compute the respective innovation statistics attributes.
            region_dict={}
            for variable in self.variables:

                # Collect the base filenames for the respective
                # experiment and analysis variable.
                filelist = parser_interface.object_getattr(
                    object_in=filelist_obj_dict[expt], key=variable, force=True
                )

                if filelist is None:
                    msg = (
                        "The innovation statistics base filenames could "
                        f"not be determined for experiment {expt} and analysis "
                        f"variable {variable}. Aborting!!!"
                    )
                    error(msg=msg)

                # Define and collect the netCDF attributes for the
                # respective experiment and analysis variables.
                ncfile = os.path.join(datapath, filelist[0])
                exist = fileio_interface.fileexist(path=ncfile)
                if not exist:
                    ncfile = filelist[0]

                # Define the levels attributes.
                nlevs = netcdf4_interface.ncreaddim(
                    ncfile=ncfile, ncdimname=self.vcoord_name
                )
                levs = netcdf4_interface.ncreadvar(
                    ncfile=ncfile, ncvarname=self.vcoord_name
                )

                # Loop through the specified regions and
                # compute/define the respective innovation statistic
                # attributes.
                region_dict[variable] = {}
                for region in vars(self.regions_obj):
                    (stats_dict, stats_mean_dict) = [{} for i in range(2)]
                    for stat_var in self.stats_vars_list:
                        stats_dict[stat_var] = numpy.zeros(
                            (len(filelist), nlevs), dtype=numpy.float64
                        )

                    for (i, filename) in enumerate(filelist):
                        ncfile = os.path.join(datapath, filename)
                        exist = fileio_interface.fileexist(path=ncfile)
                        if not exist:
                            ncfile = filename

                        for stat_var in self.stats_vars_list:
                            ncvarname = f"{stat_var}_{region}"
                            stats_dict[stat_var][i, :] = netcdf4_interface.ncreadvar(
                                ncfile=ncfile, ncvarname=ncvarname
                            )

                    # Compute the respective temporal mean innovation
                    # statistics.
                    for stat_var in self.stats_vars_list:
                        msg = (
                            f"Computing innovation statistic {stat_var} for variable "
                            f"{variable} of region {region} for experiment {expt}."
                        )
                        self.logger.info(msg=msg)

                        stats_mean_dict[stat_var] = numpy.nanmean(
                            stats_dict[stat_var], axis=0
                        )
                    region_dict[variable][region] = stats_mean_dict

                innov_stats_obj = parser_interface.object_setattr(
                    object_in=innov_stats_obj, key=expt, value=region_dict
                )

            # Append level attribute to innovation statistics object.
            innov_stats_obj = parser_interface.object_setattr(
                object_in=innov_stats_obj, key="levels", value=levs
            )

        return innov_stats_obj

    def get_common_filelist(self, filelist_obj_dict_in: dict) -> dict:
        """
        Description
        -----------

        This method determines the common file path basenames for the
        respective experiments and analysis variables.

        Parameters
        ----------

        filelist_obj_dict_in: dict

            A Python dictionary containing the file attributes for the
            respective experiments and analysis variables.

        Returns
        -------

        filelist_obj_dict_out: dict

            A Python dictionary containing the common file attributes
            for the respective experiments and analysis variables.

        Raises
        ------

        InnovStatsPlotError:

            * raised if the filelist for an experiment cannot be
              determines from the Python dictionary passed upon entry.

        """

        # Create a deep copy the input filelist dictionary.
        filelist_obj_dict_out = parser_interface.object_deepcopy(
            object_in=filelist_obj_dict_in
        )

        # Loop through each variable and experiment, determine the
        # common file basenames, and update the output filelist
        # dictionary.
        expt_names_list = filelist_obj_dict_in.keys()
        for varname in self.variables:
            filelist_all = []

            for expt_name in expt_names_list:
                expt_obj = parser_interface.dict_key_value(
                    dict_in=filelist_obj_dict_in,
                    key=expt_name,
                    force=True,
                    no_split=True,
                )

                if expt_obj is None:
                    msg = (
                        f"The filelist for experiment {expt_name} could not "
                        "be determined from the passed dictionary "
                        "(filelist_obj_dict_in) upon entry. Aborting!!!"
                    )
                    error(msg=msg)

                filelist = parser_interface.object_getattr(
                    object_in=expt_obj, key=varname
                )

                filenames = [os.path.basename(filename)
                             for filename in filelist]
                filelist_all.append(filenames)

            # Determine the common filename list.
            common_filelist = list(set.intersection(*map(set, filelist_all)))
            msg = (
                f"The common files for variable {varname} collected from "
                f"{len(expt_names_list)} experiments are the following.\n"
            )

            for common_file in sorted(common_filelist):
                msg = msg + (f"{common_file}\n")
            self.logger.info(msg=msg)

            # Update the output filelist dictionary.
            for expt_name in expt_names_list:
                expt_obj = parser_interface.object_setattr(
                    object_in=expt_obj, key=varname, value=common_filelist
                )

                filelist_obj_dict_out[expt_name] = expt_obj

        return filelist_obj_dict_out

    def get_omf_filelist(self, expt_name: str) -> object:
        """
        Description
        -----------

        This method compiles a list of innovation statistics files
        available to the respective experiment for each specified
        (i.e., allowable) analysis variable.

        Parameters
        ----------

        expt_name: str

            A Python string specifying an experiment name defined
            within the user experiment configuration.

        Returns
        -------

        filelist_obj: object

            A Python object containing the list of existing innovation
            statistics files available to the respective experiment
            for each specified (i.e., allocate) analysis variable.

        Raises
        ------

        InnovStatsPlotError:

            * raised if the attributes for the respective experiment,
              denoted by expt_name upon entry, cannot be determined
              from the user experiment configuration.

            * raised if the datapath cannot be determined for the
              respective experiment, denoted by expt_name upon entry,
              from the user experiment configuration.

            * raised if the variables attribute cannot be determined
              for the respective experiment, denote by expt_name upon
              entry, from the user experiment configuration.

            * raised if the attributes for an expected analysis
              variable for the respective experiment, denoted by
              expt_name upon entry, cannot be determined from the user
              experiment configuration.

            * raised if a required attribute for the file path
              constructor cannot be determined for the expected
              analysis variable for the respective experiment, denoted
              by expt_name upon entry, cannot be determined from the
              user experiment configuration.

        """

        # Collect a list of valid analysis times.
        timestamp_list = self.get_times(expt_name=expt_name)

        # Collect the data path attribute for the respective
        # experiment.
        expt_dict = parser_interface.object_getattr(
            object_in=self.expts_obj, key=expt_name, force=True
        )
        if expt_dict is None:
            msg = (
                f"The attributes for experiment {expt_name} could not be "
                "determined from the user experiment configuration. "
                "Aborting!!!."
            )
            error(msg=msg)

        datapath = parser_interface.dict_key_value(
            dict_in=expt_dict, key="datapath", force=True, no_split=True
        )

        if datapath is None:
            msg = (
                "The attribute datapath could not be determined for "
                f"experiment {expt_name}. Aborting!!!"
            )
            error(msg=msg)

        # Collect the variables attribute for the respective
        # experiment.
        vars_dict = parser_interface.dict_key_value(
            dict_in=expt_dict, key="variables", force=True
        )

        if vars_dict is None:
            msg = (
                "The variables attribute could not be determined "
                f"for experiment {expt_name} within the user experiment "
                "configuration. Aborting!!!"
            )
            error(msg=msg)

        # Collect the file naming attributes from the user experiment
        # configuration for the respective experiment(s) and
        # variables.
        filelist_obj = parser_interface.object_define()
        filelist_dict = {"filename": None, "filename_offset_seconds": 0}

        for var in vars_dict.keys():
            var_dict = parser_interface.dict_key_value(
                dict_in=vars_dict, key=var)

            if var_dict is None:
                msg = (
                    f"The variable attributes for variable {var} of "
                    f"experiment {expt_name} could not be determined from "
                    "the user experiment configuration. Aborting!!!"
                )
                error(msg=msg)

            for item in filelist_dict:
                value = parser_interface.dict_key_value(
                    dict_in=var_dict, key=item, force=True, no_split=True
                )
                if value is None:
                    value = parser_interface.dict_key_value(
                        dict_in=filelist_dict, key=item, no_split=True
                    )

                    if value is None:
                        msg = (
                            f"The attribute {item} for experiment {expt_name} "
                            f"and variable {var} cannot be NoneType. "
                            "Aborting!!!"
                        )
                        error(msg=msg)

                    if value is not None:
                        filelist_dict[item] = value

                if value is not None:
                    filelist_dict[item] = value

            # Build the file list for the respective experiment and
            # variables; append only existing file paths.
            filelist = []
            (filename, offset_seconds) = (
                filelist_dict["filename"],
                filelist_dict["filename_offset_seconds"],
            )

            filepath = os.path.join(datapath, filename)
            for timestamp in timestamp_list:
                filename = datetime_interface.datestrupdate(
                    datestr=timestamp,
                    in_frmttyp=timestamp_interface.GLOBAL,
                    out_frmttyp=filepath,
                    offset_seconds=offset_seconds,
                )

                exist = fileio_interface.fileexist(path=filename)
                if exist:
                    msg = f"File {filename} exists for experiment {expt_name} variable {var}."
                    self.logger.info(msg=msg)

                    filelist.append(filename)

                if not exist:
                    msg = (
                        f"File {filename} does not exist for experiment {expt_name} "
                        f"variable {var} and will not be used for creation of "
                        "the diagnostics plots."
                    )
                    self.logger.warn(msg=msg)

            filelist_obj = parser_interface.object_setattr(
                object_in=filelist_obj, key=var, value=filelist
            )

        return filelist_obj
