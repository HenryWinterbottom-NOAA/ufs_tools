# =========================================================================

# $$$ MODULE DOCUMENTATION BLOCK

# UFS-RNR-analysis :: innov_stats/py/apps/soca_ocean.py

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

    soca_ocean.py

Description
-----------

    This module contains the base-class object and methods specific to
    SOCA ocean innovation statistics applications.

Classes
-------

    SOCAOceanOMFProfile(yaml_file, basedir)

        This is the base-class object for all SOCA ocean OMF profile
        plots; it is a sub-class of OMFProfile.

    SOCAOceanOMFTimeseries(yaml_file, basedir)

        This is the base-class object for all SOCA ocean OMF
        time-series plots; it is a sub-class of OMFTimeseries.

    SOCAOceanOMFError(msg)

         This is the base-class for all exceptions; it is a sub-class
         of OMFProfileError.

Author(s)
---------

    Henry R. Winterbottom; 19 August 2022

History
-------

    2022-08-19: Henry Winterbottom -- Initial implementation.

"""

# ----

from innov_stats.omf_profile import OMFProfile
from innov_stats.omf_profile import OMFProfileError
from innov_stats.omf_timeseries import OMFTimeseries
from innov_stats.omf_timeseries import OMFTimeseriesError
from plots.omf_profile import OMFProfilePlots
from plots.omf_timeseries import OMFTimeseriesPlots
from tools import parser_interface

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


class SOCAOceanOMFProfile(OMFProfile):
    """
    Description
    -----------

    This is the base-class object for all SOCA ocean OMF plots; it is
    a sub-class of OMFProfile.

    Parameters
    ----------

    yaml_file: str

        A Python string specifying the path to the user experiment
        configuration.

    basedir: str

        A Python string specifying the path to the directory path in
        which the application is being executed.

    """

    def __init__(self, yaml_file, basedir):
        """
        Description
        -----------

        Creates a new SOCAOceanOMFProfile object.

        """

        # Define the base-class attibutes.
        super(SOCAOceanOMFProfile, self).__init__(yaml_file=yaml_file,
                                                  basedir=basedir,
                                                  is_soca_ocean=True)
        self.omf_plots = OMFProfilePlots(yaml_file=yaml_file,
                                         basedir=basedir,
                                         stats_list=self.stats_vars_list,
                                         vars_list=self.soca_ocean_variables,
                                         is_soca_ocean=True)

    def run(self):
        """
        Description
        -----------

        This method performs the following tasks:

        (1) Defines the regions of  interest specified within the user
            experiment configuration.

        (2) Collects the innovation statistics files for the
            respective experiments and computes the temporal mean.

        (3) Plots the temporal mean innovation statistics.

        """

        # Collect the regions of interest specified within the user
        # experiment configuration.
        self.get_regions()

        # Collect the innovation statistics files for the respective
        # experiments specified within the user experiment
        # configuration.
        filelist_obj_dict = dict()
        for expt_name in vars(self.expts_obj):
            filelist_obj = self.get_omf_filelist(expt_name=expt_name)
            filelist_obj_dict[expt_name] = filelist_obj
        if len(filelist_obj_dict.keys()) > 1:

            # If more the innovation statistics are being compared
            # amongst different experiments, define a list of common
            # innovation statistics files (i.e., those valid at the
            # same analysis times).
            filelist_obj_dict = self.get_common_filelist(
                filelist_obj_dict_in=filelist_obj_dict)

        # Compute the temporal mean innovation statistics for all user
        # specified experiments.
        innov_stats_obj = self.compute(filelist_obj_dict=filelist_obj_dict)

        # Plot the temporal mean innovation statistics for all user
        # specified experiments.
        self.omf_plots.plot(regions_obj=self.regions_obj,
                            innov_stats_obj=innov_stats_obj)

# ----


class SOCAOceanOMFTimeseries(OMFTimeseries):
    """
    Description
    -----------

    This is the base-class object for all SOCA ocean OMF time-series
    plots; it is a sub-class of OMFTimeseries.

    Parameters
    ----------

    yaml_file: str

        A Python string specifying the path to the user experiment
        configuration.

    basedir: str

        A Python string specifying the path to the directory path in
        which the application is being executed.

    """

    def __init__(self, yaml_file, basedir):
        """
        Description
        -----------

        Creates a new SOCAOceanOMFTimeseries object.

        """

        # Define the base-class attibutes.
        super(SOCAOceanOMFTimeseries, self).__init__(yaml_file=yaml_file,
                                                     basedir=basedir,
                                                     is_soca_ocean=True)

        self.omf_plots = OMFTimeseriesPlots(yaml_file=yaml_file,
                                            basedir=basedir,
                                            stats_list=self.stats_vars_list,
                                            vars_list=self.soca_ocean_variables,
                                            is_soca_ocean=True)

    def run(self):
        """
        Description
        -----------

        This method performs the following tasks:

        (1) Defines the regions of interest specified within the user
            experiment configuration.

        (2) Collects the innovation statistics files for the
            respective experiment.

        (3) Plots the innovation statistics time-series.

        """

        # Loop through each experiment and build the time-series for
        # each respective region and variable.
        innov_stats_obj = parser_interface.object_define()
        for expt_name in vars(self.expts_obj):
            expt_attrs = parser_interface.object_getattr(
                object_in=self.expts_obj, key=expt_name, force=True)
            if expt_attrs is None:
                msg = ('The attributes for experiment {0} could not be '
                       'determined from the user experiment configuration. '
                       'Aborting!!!'.format(expt_name))
                raise OMFTimeseriesError(msg=msg)

            # Collect the time-series information from the SQLite3
            # database file for the respective experiment.
            innov_stats_obj = parser_interface.object_setattr(
                object_in=innov_stats_obj, key=expt_name, value=self.get_timeseries(
                    expt_name=expt_name))

        # Plot the innovation statistic timeseries for all user
        # specified experiments.
        self.omf_plots.plot(innov_stats_obj=innov_stats_obj,
                            times_list=self.times_list)
