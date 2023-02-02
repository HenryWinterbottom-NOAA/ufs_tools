# =========================================================================

# $$$ MODULE DOCUMENTATION BLOCK

# UFS-RNR-analysis :: innov_stats/py/innov_stats/omf_timeseries.py

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

    omf_timeseries.py

Description
-----------

    This module contains the base-class object and methods specific to
    innovation statistic time-series applications.

Classes
-------

    OMFTimeseries(yaml_file, basedir)

         This is the base-class object for innovation statistics
         time-series calculations; it is a sub-class of InnovStats.

    InnovStatsPlotError(msg)

         This is the base-class for all exceptions; it is a sub-class
         of InnovStatsError.

Author(s)
---------

    Henry R. Winterbottom; 27 August 2022

History
-------

    2022-08-27: Henry Winterbottom -- Initial implementation.

"""

# ----

import os
import numpy

from innov_stats.plot import InnovStats, error
from ioapps import sqlite3_interface
from mycolorpy import colorlist
from tools import datetime_interface
from tools import fileio_interface
from tools import parser_interface

from innov_stats.plot.plots import Plots

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


class OMFTimeseries(InnovStats):
    """
    Description
    -----------

    This is the base-class object for innovation statistics
    time-series calculations; it is a sub-class of InnovStats.

    Parameters
    ----------

    yaml_file: str

        A Python string specifying the path to the user experiment
        configuration.

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

    def __init__(self, options_obj, basedir, is_gsi_atmos=False,
                 is_soca_ice=False, is_soca_ocean=False):
        """
        Description
        -----------

        Creates a new OMFTimeseries object.

        """

        # Define the base-class attributes.
        super().__init__(options_obj=options_obj, basedir=basedir)
        apps_dict = {is_gsi_atmos: self.gsi_atmos_variables,
                     is_soca_ice: self.soca_ice_variables,
                     is_soca_ocean: self.soca_ocean_variables
                     }
        self.get_regions()

        # Define the database table column name attributes.
        if is_gsi_atmos:
            self.column_str = 'hPa_'
            self.column_scale = 1.0
            self.vars_list = self.gsi_atmos_variables
        if is_soca_ice:
            self.column_str = 'surface_'
            self.column_scale = 1.0
            self.vars_list = self.soca_ice_variables
        if is_soca_ocean:
            self.column_str = 'meters_'
            self.column_scale = 1000.0
            self.vars_list = self.soca_ocean_variables

        # Check that the constructor arguments are valid; proceed
        # accordingly.
        if not any(apps_dict.keys()):
            msg = ('The supported innovation statistics application '
                   'has not been specified. Aborting!!!')
            error(msg=msg)
        for app in apps_dict.keys():
            if app:
                self.variables = apps_dict[app]

        # Define the plotter object.
        self.plotter = Plots(options_obj=options_obj,
                             stats_list=self.stats_vars_list,
                             vars_list=self.vars_list)

    def build_colorlevs(self, cmap, nlevs):
        """
        Description
        -----------

        This method builds a list of colors normalized with respect to
        the specified valid matplotlib colormap and specified the
        number of levels upon entry.

        Parameters
        ----------

        cmap: str

            A Python string specifying a valid matplotlib color map.

        nlevs: int

            A Python integer identifying the total number of color
            levels.

        Returns
        -------

        colorlevs: list

            A Python list of normalized color values.

        Raises
        ------

        PlotsError:

            * raised if an exception is encountered while attempting
              to define the respective color levels.

        """

        # Define the color levels relative to the specified
        # attributes.
        try:
            data_arr = numpy.linspace(0, 1.0, nlevs)
            colorlevs = colorlist.gen_color_normalized(
                cmap=cmap, data_arr=data_arr)
        except Exception as errmsg:
            msg = ('The normalized color levels generation failed '
                   f'with error {errmsg}. Aborting!!!')
            error(msg=msg)

        return colorlevs

    def get_levels(self, sql_path, table_name):
        """
        Description
        -----------

        This method returns the vertical level values for the
        respective application; the levels are determined from the
        respective database table name within the SQLite3 database
        path.

        Parameters
        ----------

        sql_path: str

            A Python string specifying the path to the SQLite3
            database file path.

        table_name: str

            A Python string specifying SQLite3 database table name.

        Returns
        -------

        levels: list

            A Python list of vertical level values; these values are
            scaled by the attributes corresponding to the respective
            application.

        """

        # Collect the relevant column names from the SQLite3 database
        # table.
        columns = sqlite3_interface.read_columns(
            path=sql_path, table_name=table_name)
        values = [column for column in columns if column != 'CYCLE']

        # Scale the vertical level values accordingly.
        levels = list()
        for value in values:
            level = float(value.replace(self.column_str, ''))/self.column_scale
            levels.append(level)
        return levels

    def get_stats_regions(self, expt_dict, variable):
        """
        Description
        -----------

        This method defines a list of regions for which to plot the
        innovation statistics profiles.

        Parameters
        ----------

        expt_dict: dict

            A Python dictionary containing the experiment innovation
            statistics attributes.

        variable: str

            A Python string specifying the innovation statistics
            variable.

        Returns
        -------

        regions_list: list

            A Python list speciing the respective regions for the
            innovation statistics variable specified upon entry.

        Raises
        ------

        PlotsError:

            * raised if the regions for which to plot the respective
              variable innnovation statistics cannot be determined
              from the base-class attribute expt_obj upon entry.

        """

        # Define the list of regions for the respective innovation
        # statistics variable.
        regions_list = list(parser_interface.dict_key_value(
            dict_in=expt_dict, key=variable, force=True,
            no_split=True).keys())
        if regions_list is None:
            msg = ('The regions for which to plot the respective '
                   'innovation statistics could not be determined from '
                   'the base-class attribute expt_obj upon entry. '
                   'Aborting!!!')
            raise PlotsError(msg=msg)

        return regions_list

    def get_timeseries(self, expt_name):
        """
        Description
        -----------

        This method parses the SQLite3 database file for the
        respective experiment and compiles a Python dictionary
        containing the time-series attributes as a function of region,
        statistic, and level for all times within the temporal window
        specified in the user experiment configuration; any requested
        time-stamps that do not exist within the specified temporal
        window are set to NoneType.

        Parameters
        ----------

        expt_name: str

            A Python string specifying the experiment name; this
            string is used to identify the respective SQLite3 database
            file(s) and build the Python dictionary accordingly.

        Returns
        -------

        timeseries_dict: dict

            A Python dictionary containing the respective experiment
            time-series attributes.

        Raises
        ------

        InnovStatsPlotError:

            * raised if the variable attributes could not be
              determined from the user experiment configuration for
              the respective experiment.

            * raised if the SQLite3 database directory file path
              cannot be determined from the user experiment
              configuration for the respective experiment.

            * raised if the SQLite3 database basename file path cannot
              be from the user experiment configuration for the
              respective experiment.

            * raised if the SQLite3 database path for the respective
              experiment does not exist.

        """

        # Collect the experiment attributes from the user experiment
        # configuration file and proceed accordingly.
        expt_dict = parser_interface.object_getattr(
            object_in=self.expts_obj, key=expt_name, force=True)
        vars_dict = parser_interface.dict_key_value(
            dict_in=expt_dict, key='variables', force=True, no_split=True)
        if vars_dict is None:
            msg = ('The variable attributes could not be determined for '
                   'experiment {0}. Aborting!!!'.format(expt_name))
            error(msg=msg)

        # Define the basename path for the variable SQLite3 database
        # files for the respective experiment.
        datapath = parser_interface.dict_key_value(
            dict_in=expt_dict, key='datapath', force=True, no_split=True)
        if datapath is None:
            msg = ('The datapath attribute for experiment {0} could not '
                   'be determined from the user experiment configuration. '
                   'Aborting!!!'.format(expt_name))
            error(msg=msg)

        # Define the list of timestamps for the respective experiment.
        self.times_list = self.get_times(expt_name=expt_name)

        # Loop through each variable and build the time-series object
        # for the respective experiment; the time-series object
        # contains the time-series statistics as a function of the
        # regions specified within the user experiment configuration.
        timeseries_dict = dict()
        for variable in vars_dict.keys():

            # Define the SQLite3 filename path for the respective
            # variable and confirm that the file path exists.
            timeseries_dict[variable] = dict()
            filename = parser_interface.dict_key_value(
                dict_in=vars_dict[variable], key='filename', force=True,
                no_split=True)
            if filename is None:
                msg = ('The basename for the SQLite3 database file for '
                       'variable {0} could not be determined from the user '
                       'experiment configuration. Aborting!!!'.format(variable))
                error(msg=msg)
            sql_path = os.path.join(datapath, filename)
            exist = fileio_interface.fileexist(path=sql_path)
            if not exist:
                msg = ('The SQLite3 database file path {0} does not exist. '
                       'Aborting!!!'.format(sql_path))
                error(msg=msg)

            # Parse the SQLite3 database file path and collect the
            # timeseries attributes for the respective experiment and
            # respective regions.
            for region in vars(self.regions_obj):
                for stats_var in self.stats_vars_list:
                    table_name = '{0}_{1}'.format(stats_var, region)
                    timeseries_dict[variable][table_name] = dict()
                    table_dict = sqlite3_interface.read_table(
                        path=sql_path, table_name=table_name)

                    levels = self.get_levels(
                        sql_path=sql_path, table_name=table_name)

                    # Loop through all levels and proceed accordingly.
                    for (lidx, level) in enumerate(levels):

                        timeseries_dict[variable][table_name][level] = list()
                        # Loop through the list of requested times and
                        # proceed accordingly.
                        for (tidx, time) in enumerate(self.times_list):

                            # Find all timestamps corresponding to the
                            # requested time; build the generator
                            # function and proceed accordingly.
                            try:

                                value = (row for row in
                                         list(table_dict.values()) if int(time) ==
                                         int(row[0]))

                                timeseries_dict[variable][table_name][level].append(
                                    next(value)[lidx+1])
                            except StopIteration:
                                timeseries_dict[variable][table_name][level].append(
                                    numpy.nan)
        return timeseries_dict
