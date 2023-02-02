# =========================================================================

# $$$ MODULE DOCUMENTATION BLOCK

# UFS-RNR-analysis :: apps/innov_stats/py/plots/omf_timeseries.py

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

    This module contains the base-class for all OMF timeseries
    plotting.

Classes
-------

    OMFTimeseriesPlots(yaml_file, basedir, vars_list=None, 
                       stats_list=None, is_gsi_atmos=False, 
                       is_soca_ice=False, is_soca_ocean=False)

        This is the base-class object for all OMF timeseries plots; it
        is a sub-class of Plots.

Author(s)
---------

    Henry R. Winterbottom; 27 August 2022

History
-------

    2022-08-27: Henry Winterbottom -- Initial implementation.

"""

# ----

import numpy
import re

from itertools import product
from innov_stats.plotting import Plots
#from plots import PlotsError
from tools import parser_interface

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


class OMFTimeseriesPlots(Plots):
    """
    Description
    -----------

    This is the base-class object for all OMF timeseries plots; it is
    a sub-class of Plots.

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

    vars_list: list, optional

        A Python list of variables for which to collect the user
        specified attributes.

    stats_list: list, optional

        A Python list of statistics specified within the user
        experiment configuration for each respective variable (see
        vars_list).

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

    PlotsError:

        * raised if the argument(s) pertaining the specified
          innovation statistics plots method are not valid.

    """

    def __init__(self, yaml_file, basedir, vars_list=None,
                 stats_list=None, is_gsi_atmos=False, is_soca_ice=False,
                 is_soca_ocean=False):
        """
        Description
        -----------

        Creates a new OMFTimeseriesPlots object.

        """

        # Define the base-class attributes.
        super(OMFTimeseriesPlots, self).__init__(yaml_file=yaml_file,
                                                 basedir=basedir,
                                                 vars_list=vars_list,
                                                 stats_list=stats_list)

        self.expt_plot_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key='experiments', force=True,
            no_split=True)
        if self.expt_plot_dict is None:
            msg = ('The respective experiments plotting attributes '
                   'could not be determined from the user experiment '
                   'configuration. Aborting!!!')
            raise PlotsError(msg=msg)

        # Check that the variable attributes within the user
        # experiment configuration are valid.
        self.var_plot_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key='variables', force=True,
            no_split=True)
        if self.var_plot_dict is None:
            msg = ('The respective variables plotting attributes '
                   'could not be determined from the user experiment '
                   'configuration. Aborting!!!')
            raise PlotsError(msg=msg)

        # Check that the constructor arguments are valid; proceed
        # accordingly.
        if not any([is_gsi_atmos, is_soca_ice, is_soca_ocean]):
            msg = ('The supported innovation statistics application '
                   'has not been specified. Aborting!!!')
            raise PlotsError(msg=msg)
        (self.is_gsi_atmos, self.is_soca_ice,
         self.is_soca_ocean) = (is_gsi_atmos, is_soca_ice, is_soca_ocean)

    def plot(self, innov_stats_obj, times_list):
        """
        Description
        -----------

        This method plots the OMF timeseries as a function of region
        and a respective analysis variable's innovation statistic(s).

        Parameters
        ----------

        innov_stats_obj: obj

            A Python object containing the respective computed
            temporal mean innovation statistics for each experiment
            and respective region(s).

        times_list: list

            A Python list containing the times for which to construct
            the time-series plot(s).

        Raises
        ------

        PlotsError:

            * raised if the levels for which to plot the respective
              innovation statistic(s) has not been specified in the
              user experiment configuration.

            * raised if a specified innovation statistic has not be
              defined for a respective analysis variable and
              experiment.

            * raised if an innovation statistic timeseries cannot be
              determined for a specified experiment.

            * raised if an exception is encountered while plotting the
              respective innovation statistic.

        """

        # Collect the levels attributes from the user experiment
        # configuration.
        levels = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key='levels', force=True,
            no_split=True)
        if levels is None:
            msg = ('The user experiment configuration has not specified '
                   'levels for which to plot the respective innovation '
                   'statistic timeseries. Aborting!!!')
            raise PlotsError(msg=msg)

        if levels is not None:
            try:
                levels_list = list(numpy.float_(levels.split(',')))
            except AttributeError:
                levels_list = [levels
                               ]
            except ValueError:
                pass
            msg = ('Timeseries will be plotted for the following levels.\n')
            for level in levels_list:
                msg = msg + '{0}\n'.format(level)
            self.logger.warn(msg=msg)

        # Define the color attributes for the respective levels of the
        # time-series.
        cmap = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key='colormap', force=True, no_split=True)
        if cmap is None:
            cmap = 'viridis'
            msg = ('The user experiment does not specify the color map '
                   'attribute; assigning the colormap as {0}.'.
                   format(cmap))
            self.logger.warn(msg=msg)
        colorlevs = self.build_colorlevs(cmap=cmap, nlevs=len(levels_list))

        # Collect the innovation statistics attributes.
        experiments_list = self._get_experiments(
            innov_stats_obj=innov_stats_obj)
        expt_dict = parser_interface.object_getattr(object_in=innov_stats_obj,
                                                    key=experiments_list[0])
        stats_regions_list = self._get_regions(expt_dict=expt_dict,
                                               variable=self.vars_list[0])

        # Loop through each innovation statistic and region and
        # proceed accordingly.
        for stat_region in stats_regions_list:

            # Define the respective statistics variable and region
            # values.
            stat = stat_region.split('_')[0]
            region = re.sub('{0}_'.format(stat), '', stat_region)

            # Loop through each analysis variable and proceed
            # accordingly.
            for variable in self.vars_list:

                # Initialilze the figure attributes.
                (plot_obj, fig_obj, ax_obj) = self.build_figure()
                ax_obj.spines['top'].set_visible(False)
                ax_obj.spines['right'].set_visible(False)

                # Define the plot attributes for the respective
                # variable from the user experiment configuration and
                # proceed accordingly.
                var_plot_dict = parser_interface.dict_key_value(
                    dict_in=self.var_plot_dict[variable], key=stat, force=True,
                    no_split=True)
                if var_plot_dict is None:
                    msg = ('The plot attributes for variable {0} and statistics '
                           'region variable {1} could not be determined from the '
                           'user experiment configuration; the variable will not '
                           'be plotted.'.format(variable, stat_region))
                    self.logger.warn(msg=msg)

                if var_plot_dict is not None:

                    # Define the plotting object attributes and
                    # proceed accordingly.
                    keys_dict = {'xmin': 0, 'xmax': len(times_list)}
                    for key in keys_dict.keys():
                        value = parser_interface.dict_key_value(
                            dict_in=keys_dict, key=key)
                        if 'axes' in var_plot_dict.keys():
                            var_plot_dict['axes'][key] = value
                        if 'hlines' in var_plot_dict.keys():
                            var_plot_dict['hlines'][key] = value

                    # Loop through each experiment and plot the
                    # computed innovation statistics accordingly.
                    for experiment in experiments_list:

                        # Collect the vertical levels available for
                        # the respective experiment innovation
                        # statistic.
                        expt_dict = parser_interface.object_getattr(
                            object_in=innov_stats_obj, key=experiment, force=True)
                        if expt_dict is None:
                            msg = ('The innovation statistics for experiment {0} could '
                                   'not be determined from the user experiment '
                                   'configuration. Aborting!!!'.format(experiment))
                            raise PlotsError(msg=msg)

                        var_dict = parser_interface.dict_key_value(
                            dict_in=expt_dict, key=variable, force=True, no_split=True)
                        if var_dict is None:
                            msg = ('The {0} variable attributes could not be determined '
                                   'for experiment {1}. Aborting!!!'.format(variable, experiment))
                            raise PlotsError(msg=msg)

                        expt_stat_region = parser_interface.dict_key_value(
                            dict_in=var_dict, key=stat_region, force=True, no_split=True)
                        if expt_stat_region is None:
                            msg = ('The attributes could not be determined for innovation '
                                   'statistic region {0} variable {1} for experiment {2}. '
                                   'Aborting!!!'.format(stat_region, variable, experiment))
                            raise PlotsError(msg=msg)

                        expt_levels = list(expt_stat_region.keys())

                        # Loop through each specified vertical level;
                        # check whether the level already exists for
                        # the respective application; proceed
                        # accordingly.
                        for (cidx, level) in enumerate(levels_list):

                            # If level is define, collect the index
                            # and plot the values.
                            if level in expt_levels:

                                # Define the respective analysis
                                # variable time-series.
                                var_ts = parser_interface.dict_key_value(
                                    dict_in=expt_stat_region, key=level, force=True)
                                if var_ts is None:
                                    msg = ('The innovation statistic region {0} time-series '
                                           'for variable {1} could not be determined for '
                                           'experiment {2}. Aborting!!!'.format(stat_region,
                                                                                variable,
                                                                                experiment))
                                    raise PlotsError(msg=msg)

                            # If the respective level is not defined
                            # for the respective experiment, attempt
                            # to interpolate to the respective level;
                            # if unable to interpolate, assign as
                            # numpy.nan.
                            if level not in expt_levels:

                                # Find the bounding levels relative to
                                # the respective level.
                                level_list = [level]
                                order = sorted(product(expt_levels, level_list), key=lambda t:
                                               abs(t[0] - t[1]))
                                (minval, maxval) = (order[0][0], order[1][0])

                                # Compute the weights relative to the
                                # respective level.
                                (w1, w2) = ((level - minval)/(maxval - minval),
                                            (maxval - level)/(maxval - minval))

                                # Define the interpolated time-series.
                                var_ts1 = w1*numpy.asarray(parser_interface.dict_key_value(
                                    dict_in=expt_stat_region, key=minval))
                                var_ts2 = w2*numpy.asarray(parser_interface.dict_key_value(
                                    dict_in=expt_stat_region, key=maxval))
                                var_ts = var_ts1 + var_ts2

                                # Check that the level is valid;
                                # proceed accordingly.
                                if (level > minval) or (level < maxval):
                                    msg = ('The level {0} is in valid since it is not '
                                           'within the interval ({1}, {2}); setting to '
                                           '{3}.'.format(level, min(expt_levels),
                                                         max(expt_levels), numpy.nan))
                                    self.logger.warn(msg=msg)
                                    var_ts[:] = numpy.nan

                            # Plot the time-series for the respective
                            # level.
                            color = colorlevs[cidx]
                            kwargs = parser_interface.dict_key_value(
                                dict_in=var_plot_dict, key='lines', force=True)
                            if kwargs is None:
                                kwargs = dict()
                            try:
                                plot_obj.plot(var_ts, color=color,
                                              label=level, **kwargs)
                            except TypeError:
                                plot_obj.plot(var_ts, label=level, **kwargs)

                    # Customize the figure aspects based on the user
                    # experiment configuration and proceed
                    # accordingly.

                    plot_regions_title_dict = parser_interface.dict_key_value(
                        dict_in=self.yaml_dict['regions'][region], key='title',
                        force=True)

                    # Build the figure attributes for the respective
                    # timeseries variable and statistic.
                    self.build_ylabel(ax_obj=ax_obj, **var_plot_dict)
                    self.build_legend(plot_obj=plot_obj, **var_plot_dict)
                    var_plot_dict['axes']['xmin'] = 0
                    var_plot_dict['axes']['xmax'] = len(times_list)
                    xint = parser_interface.dict_key_value(
                        dict_in=var_plot_dict['axes'], key='xint', force=True)
                    if xint is None:
                        xint = 1
                    var_plot_dict['axes']['xlabels'] = times_list[::xint]
                    self.build_axes(plot_obj=plot_obj, ax_obj=ax_obj,
                                    **var_plot_dict)
                    self.build_hlines(plot_obj=plot_obj, **var_plot_dict)

                    # Build the figure title.
                    if plot_regions_title_dict is not None:
                        self.build_title(plot_obj, **plot_regions_title_dict)

                    # Build the output figure name and create the
                    # figure.
                    savename = '{0}.{1}.{2}.png'.format(region, stat, variable)
                    self.create_image(plot_obj=plot_obj, savename=savename,
                                      **var_plot_dict)
