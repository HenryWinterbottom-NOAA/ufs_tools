# =========================================================================

# Module: ush/innov_stats/plot/plots/omf_timeseries_plot.py

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

    omf_timeseries_plot.py

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

from tools import parser_interface
from typing import Tuple

from innov_stats.plot import error
from innov_stats.plot.omf_timeseries import OMFTimeseries

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


class OMFTimeseriesPlots(OMFTimeseries):
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

    def __init__(self, options_obj, basedir, is_gsi_atmos=False, is_soca_ice=False,
                 is_soca_ocean=False):
        """
        Description
        -----------

        Creates a new OMFTimeseriesPlots object.

        """

        # Define the base-class attributes.
        super().__init__(options_obj=options_obj,
                         basedir=basedir,
                         is_gsi_atmos=is_gsi_atmos,
                         is_soca_ice=is_soca_ice,
                         is_soca_ocean=is_soca_ocean)

        # Collect the experiment attributes within the user experiment
        # configuration.
        self.expt_plot_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key='experiments', force=True,
            no_split=True)
        if self.expt_plot_dict is None:
            msg = ('The respective experiments plotting attributes '
                   'could not be determined from the user experiment '
                   'configuration. Aborting!!!')
            error(msg=msg)

        # Check that the variable attributes within the user
        # experiment configuration are valid.
        self.var_plot_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key='variables', force=True,
            no_split=True)
        if self.var_plot_dict is None:
            msg = ('The respective variables plotting attributes '
                   'could not be determined from the user experiment '
                   'configuration. Aborting!!!')
            error(msg=msg)

        # Check that the constructor arguments are valid; proceed
        # accordingly.
        if not any([is_gsi_atmos, is_soca_ice, is_soca_ocean]):
            msg = ('The supported innovation statistics application '
                   'has not been specified. Aborting!!!')
            error(msg=msg)

        (self.is_gsi_atmos, self.is_soca_ice,
         self.is_soca_ocean) = (is_gsi_atmos, is_soca_ice, is_soca_ocean)

    def __plotstat__(
        self, innov_stats_obj: object, region: str, stat: str, variable: str
    ) -> Tuple[object, object]:
        """
        Description
        -----------

        This method plots the computed regional variable innovation
        statistic for all defined experiments.

        Parameters
        ----------

        innov_stats_obj: object

            A Python object containing the respective computed
            temporal mean innovation statistics for each experiment
            and respective region(s).

        region: str

            A Python string specifiying the region for which to
            plot the respective variable innovations statistics.

        stat: str

            A Python string specifying the innovation statistic to be
            plotted.

        variable: str

            A Python string specifying the variable for which to plot
            the respective innovation statistic.

        Returns
        -------

        plot_obj: object

            A Python matplotlib pyplot object.

        ax_obj: object

            A Python matplotlib pyplot axis-level subplot object

        Raises
        ------

        InnovStatsPlotError:

            * raised if an analysis variable specified in the user
              experiment configuration has not been computed for the
              respective experiment.

            * raised if a region specified in the user experiment
              configuration has not been defined for a specified
              analysis variable and region.

            * raised if a specified innovation statistic has not be
              defined for a respective analysis variable and
              experiment.

            * raised if an exception is encountered while plotting the
              respective innovation statistic.

        """

        # Collect the list of experiments from the experiment
        # configuration.
        experiments_list = list(vars(self.expts_obj))

        # Define the plotting object attributes and proceed
        # accordingly.
        (plot_obj, _, ax_obj) = self.plotter.build_figure()

        # Loop through the respective experiments.
        for experiment in experiments_list:

            # Collect the respective experiment plotting keyword
            # arguments.
            plot_kwargs = parser_interface.dict_key_value(
                dict_in=self.expt_plot_dict[experiment],
                key="plotting",
                force=True,
            )

            # Collect all attributes for the respective region,
            # variable, experiment, and innovation statistic.
            expt_dict = parser_interface.object_getattr(
                object_in=innov_stats_obj, key=experiment
            )

            var_dict = parser_interface.dict_key_value(
                dict_in=expt_dict, key=variable, force=True, no_split=True
            )
            if var_dict is None:
                msg = (
                    f"The variable {variable} has not been defined for "
                    f"experiment {experiment}. Aborting!!!"
                )
                error(msg=msg)

            region_dict = parser_interface.dict_key_value(
                dict_in=var_dict, key=region, force=True, no_split=True
            )
            if region_dict is None:
                msg = (
                    f"The region {region} has not be defined for variable "
                    f"{variable} for experiment {experiment}. Aborting!!!"
                )
                error(msg=msg)

            statvals = parser_interface.dict_key_value(
                dict_in=region_dict, key=stat, force=True, no_split=True
            )
            if statvals is None:
                msg = (
                    f"The statistic {stat} values could not be determined "
                    f"for variable {variable} experiment {experiment}. "
                    "Aborting!!!"
                )
                error(msg=msg)

            # Plot the respective innovation statistic for the
            # respective region, variable, and experiment.
            if plot_kwargs is None:
                msg = (
                    f"Plotting arguments for innovation statistic {stat} "
                    f"for variable {variable} of experiment {experiment} "
                    "could not be determined from the user experiment "
                    "configuration; plotting will not be performed."
                )
                self.logger.warn(msg=msg)
                break

            # Plot the respective innovation statistics
            # for the respective experiment and region.
            try:
                plot_obj.plot(statvals, innov_stats_obj.levels, **plot_kwargs)

            except Exception as errmsg:
                msg = (
                    f"The plotting of innovation statistic {stat} for experiment "
                    f"{experiment} variable {variable} in region {region} "
                    f"failed with error {errmsg}. Aborting!!!"
                )

                error(msg=msg)

        return (plot_obj, ax_obj)

    def plot(self, innov_stats_obj: object) -> None:
        """
        Description
        -----------

        This method plots the OMF profiles as a function of region and
        a respective analysis variable's innovation statistic(s).

        Parameters
        ----------

        innov_stats_obj: object

            A Python object containing the respective computed
            temporal mean innovation statistics for each experiment
            and respective region(s).

        """

        # Collect the innovations statistics attributes; proceed
        # accordingly.
        experiments_list = list(vars(self.expts_obj))
        regions_list = list(vars(self.regions_obj))

        if "levels" in experiments_list:
            experiments_list.remove("levels")

        # Loop through the respective region(s), variable(s), and
        # experiment(s) for which the innovation statistics have been
        # computed and proceed accordingly.
        for region in regions_list:
            for variable in self.vars_list:

                # Define attributes for the respective innovation
                # statistics; if the plotting attributes for the
                # respective experiment cannot be determined, the
                # keyword arguments for the matplotlib.pyplot
                # object will be set to NoneType.
                for stat in self.stats_vars_list:

                    (plot_obj, ax_obj) = self.__plotstat__(
                        innov_stats_obj=innov_stats_obj,
                        region=region,
                        stat=stat,
                        variable=variable,
                    )
                    plot_obj = self.__plotattrs__(
                        ax_obj=ax_obj,
                        plot_obj=plot_obj,
                        region=region,
                        stat=stat,
                        variable=variable,
                    )
                    self.__plotsave__(
                        plot_obj=plot_obj, region=region, stat=stat, variable=variable
                    )
        

#    def plot(self, innov_stats_obj, times_list):
#        """
#        Description
#        -----------###

#        This method plots the OMF timeseries as a function of region
#        and a respective analysis variable's innovation statistic(s).#

#        Parameters
#        ----------

#        innov_stats_obj: obj

#            A Python object containing the respective computed
#            temporal mean innovation statistics for each experiment
#            and respective region(s).

#        times_list: list

#            A Python list containing the times for which to construct
#            the time-series plot(s).

#        Raises
#        ------

#        PlotsError:

#            * raised if the levels for which to plot the respective
#              innovation statistic(s) has not been specified in the
#              user experiment configuration.

#            * raised if a specified innovation statistic has not be
#              defined for a respective analysis variable and
#              experiment.

#            * raised if an innovation statistic timeseries cannot be
#              determined for a specified experiment.

#            * raised if an exception is encountered while plotting the
#              respective innovation statistic.

#        """

        # Collect the levels attributes from the user experiment
        # configuration.
#        levels = parser_interface.dict_key_value(
#            dict_in=self.yaml_dict, key='levels', force=True,
#            no_split=True)
#        if levels is None:
#            msg = ('The user experiment configuration has not specified '
#                   'levels for which to plot the respective innovation '
#                   'statistic timeseries. Aborting!!!')
#            error(msg=msg)###

#        if levels is not None:
#            try:
#                levels_list = list(numpy.float_(levels.split(',')))

#            except AttributeError:
#                levels_list = [levels
                               ]

#            except ValueError:
#                pass

#            msg = ('Timeseries will be plotted for the following levels.\n')
#            for level in levels_list:
#                msg = msg + '{0}\n'.format(level)
#            self.logger.warn(msg=msg)

        # Define the color attributes for the respective levels of the
        # time-series.
#        cmap = parser_interface.dict_key_value(
#            dict_in=self.yaml_dict, key='colormap', force=True, no_split=True)
#        if cmap is None:
#            cmap = 'viridis'
#            msg = ('The user experiment does not specify the color map '
#                   'attribute; assigning the colormap as {0}.'.
#                   format(cmap))
#            self.logger.warn(msg=msg)
#        colorlevs = self.build_colorlevs(cmap=cmap, nlevs=len(levels_list))

        # Collect the innovation statistics attributes.
#        experiments_list = list(self.expt_plot_dict)
#        expt_dict = parser_interface.object_getattr(object_in=innov_stats_obj,
#                                                    key=experiments_list[0])

#        print(self.regions_obj)
#        quit()

#        stats_regions_list = self.get_regions()  # (expt_dict=expt_dict,
        # variable=list(self.var_plot_dict))  # self.vars_list[0])

        # Loop through each innovation statistic and region and
        # proceed accordingly.
#        for stat_region in vars(self.regions_obj):  # stats_regions_list:

            # Define the respective statistics variable and region
            # values.
#            stat = stat_region.split('_')[0]
#            region = re.sub('{0}_'.format(stat), '', stat_region)

            # Loop through each analysis variable and proceed
            # accordingly.
#            for variable in vars(self.vars_obj):  # self.vars_list:

                # Initialilze the figure attributes.
#                (plot_obj, _, ax_obj) = self.plotter.build_figure()
#                ax_obj.spines['top'].set_visible(False)
#                ax_obj.spines['right'].set_visible(False)

                # Define the plot attributes for the respective
                # variable from the user experiment configuration and
                # proceed accordingly.

#                print(variable)
#                quit()

#                var_plot_dict = parser_interface.dict_key_value(
#                    dict_in=self.var_plot_dict[variable], key=stat, force=True,
#                    no_split=True)

#                if var_plot_dict is None:
#                    msg = ('The plot attributes for variable {0} and statistics '
#                           'region variable {1} could not be determined from the '
#                           'user experiment configuration; the variable will not '
#                           'be plotted.'.format(variable, stat_region))
#                    self.logger.warn(msg=msg)

#                if var_plot_dict is not None:

                    # Define the plotting object attributes and
                    # proceed accordingly.
#                    keys_dict = {'xmin': 0, 'xmax': len(times_list)}
#                    for key in keys_dict.keys():
#                        value = parser_interface.dict_key_value(
#                            dict_in=keys_dict, key=key)
#                        if 'axes' in var_plot_dict.keys():
#                            var_plot_dict['axes'][key] = value
#                        if 'hlines' in var_plot_dict.keys():
#                            var_plot_dict['hlines'][key] = value

                    # Loop through each experiment and plot the
                    # computed innovation statistics accordingly.
#                    for experiment in experiments_list:

                        # Collect the vertical levels available for
                        # the respective experiment innovation
                        # statistic.
#                        expt_dict = parser_interface.object_getattr(
#                            object_in=innov_stats_obj, key=experiment, force=True)
#                        if expt_dict is None:
#                            msg = ('The innovation statistics for experiment {0} could '
#                                   'not be determined from the user experiment '
#                                   'configuration. Aborting!!!'.format(experiment))
#                            error(msg=msg)

#                        var_dict = parser_interface.dict_key_value(
#                            dict_in=expt_dict, key=variable, force=True, no_split=True)
#                        if var_dict is None:
#                            msg = ('The {0} variable attributes could not be determined '
#                                   'for experiment {1}. Aborting!!!'.format(variable, experiment))
#                            error(msg=msg)

#                        expt_stat_region = parser_interface.dict_key_value(
#                            dict_in=var_dict, key=stat_region, force=True, no_split=True)
#                        if expt_stat_region is None:
#                            msg = ('The attributes could not be determined for innovation '
#                                   'statistic region {0} variable {1} for experiment {2}. '
#                                   'Aborting!!!'.format(stat_region, variable, experiment))
#                            error(msg=msg)

#                        expt_levels = list(expt_stat_region.keys())

                        # Loop through each specified vertical level;
                        # check whether the level already exists for
                        # the respective application; proceed
                        # accordingly.
#                        for (cidx, level) in enumerate(levels_list):

                            # If level is define, collect the index
                            # and plot the values.
#                            if level in expt_levels:

                                # Define the respective analysis
                                # variable time-series.
#                                var_ts = parser_interface.dict_key_value(
#                                    dict_in=expt_stat_region, key=level, force=True)
#                                if var_ts is None:
#                                    msg = ('The innovation statistic region {0} time-series '
#                                           'for variable {1} could not be determined for '
#                                           'experiment {2}. Aborting!!!'.format(stat_region,
#                                                                                variable,
#                                                                                experiment))
#                                    error(msg=msg)

                            # If the respective level is not defined
                            # for the respective experiment, attempt
                            # to interpolate to the respective level;
                            # if unable to interpolate, assign as
                            # numpy.nan.
#                            if level not in expt_levels:

                                # Find the bounding levels relative to
                                # the respective level.
#                                level_list = [level]
#                                order = sorted(product(expt_levels, level_list), key=lambda t:
#                                               abs(t[0] - t[1]))
#                                (minval, maxval) = (order[0][0], order[1][0])

                                # Compute the weights relative to the
                                # respective level.
#                                (w1, w2) = ((level - minval)/(maxval - minval),
#                                            (maxval - level)/(maxval - minval))

                                # Define the interpolated time-series.
#                                var_ts1 = w1*numpy.asarray(parser_interface.dict_key_value(
#                                    dict_in=expt_stat_region, key=minval))
#                                var_ts2 = w2*numpy.asarray(parser_interface.dict_key_value(
#                                    dict_in=expt_stat_region, key=maxval))
#                                var_ts = var_ts1 + var_ts2

                                # Check that the level is valid;
                                # proceed accordingly.
#                                if (level > minval) or (level < maxval):
#                                    msg = ('The level {0} is in valid since it is not '
#                                           'within the interval ({1}, {2}); setting to '
#                                           '{3}.'.format(level, min(expt_levels),
#                                                         max(expt_levels), numpy.nan))
#                                    self.logger.warn(msg=msg)
#                                    var_ts[:] = numpy.nan

                            # Plot the time-series for the respective
                            # level.#
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
