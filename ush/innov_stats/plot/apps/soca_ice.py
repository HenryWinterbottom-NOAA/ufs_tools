# =========================================================================

# $$$ MODULE DOCUMENTATION BLOCK

# UFS-RNR-analysis :: innov_stats/py/apps/soca_ice.py

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

    soca_ice.py

Description
-----------

    This module contains the base-class object and methods specific to
    SOCA ice innovation statistics applications.

Classes
-------

    SOCAOIceOMFTimeseries(yaml_file, basedir)

        This is the base-class object for all SOCA ice OMF time-series
        plots; it is a sub-class of OMFTimeseries.

    SOCAOceanOMFError(msg)

         This is the base-class for all exceptions; it is a sub-class
         of OMFProfileError.

Author(s)
---------

    Henry R. Winterbottom; 30 August 2022

History
-------

    2022-08-30: Henry Winterbottom -- Initial implementation.

"""

# ----

from innov_stats.omf_timeseries import OMFTimeseries
from innov_stats.omf_timeseries import OMFTimeseriesError
from plots.omf_timeseries import OMFTimeseriesPlots
from tools import parser_interface

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


class SOCAIceOMFTimeseries(OMFTimeseries):
    """
    Description
    -----------

    This is the base-class object for all SOCA ice OMF plots; it is a
    sub-class of OMFProfile.

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

        Creates a new SOCAIceOMFTimeseries object.

        """

        # Define the base-class attributes.
        super(SOCAIceOMFTimeseries, self).__init__(yaml_file=yaml_file,
                                                   basedir=basedir,
                                                   is_soca_ice=True)

        self.omf_plots = OMFTimeseriesPlots(yaml_file=yaml_file,
                                            basedir=basedir,
                                            stats_list=self.stats_vars_list,
                                            vars_list=self.soca_ice_variables,
                                            is_soca_ice=True)

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
