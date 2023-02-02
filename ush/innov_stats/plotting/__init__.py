# =========================================================================

# Module: ush/innov_stats/plotting/__init__.py

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

    This module contains the base-class for all plotting applications.

Classes
-------

    Plots(yaml_file, basedir, vars_list=None, stats_list=None)

        This is the base-class object for all plotting applications.

Functions
---------

    error(msg)

        This function is the exception handler for the respective
        module.

Author(s)
---------

    Henry R. Winterbottom; 18 August 2022

History
-------

    2022-08-18: Henry Winterbottom -- Initial implementation.

"""

# ----

from dataclasses import dataclass

import matplotlib.pyplot
import mycolorpy
import numpy

from mycolorpy import colorlist


from tools import fileio_interface
from tools import parser_interface
from utils.error_interface import msg_except_handle
from utils.logger_interface import Logger

from confs.yaml_interface import YAML
from exceptions import InnovStatsPlottingError

from tools import datetime_interface

from utils import timestamp_interface

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


class Plots:
    """
    Description
    -----------

    This is the base-class object for all plotting applications.

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

    Raises
    ------

    InnovStatsPlottingError:

        * raised if either vars_list or stats_list are NoneType upon
          entry.

    """

    def __init__(self, options_obj: object, basedir: str,
                 vars_list: list = None, stats_list: list = None):
        """
        Description
        -----------

        Creates a new Plots object.

        """

        # Define the base-class attributes.
        self.options_obj = options_obj
        self.basedir = basedir
        self.yaml_file = self.options_obj.yaml_file
        self.logger = Logger()
        self.yaml_dict = YAML().read_yaml(yaml_file=self.yaml_file)
        self.pyplot = matplotlib.pyplot

        # Check the attributes specified upon entry; proceed
        # accordingly.
        (self.stats_list, self.vars_list) = (
            stats_list, vars_list)

        if any([self.stats_list, self.vars_list]) is None:
            msg = ('The attributes stats_list and vars_list must either '
                   'both be NoneType or both be not NoneType; received '
                   f'the following upon entry:\n'
                   f'\t stats_list = {self.stats_list}\n'
                   f'\t vars_list = {self.vars_list}\n'
                   '\Aborting!!!')
            error(msg=msg)

        self.cycles_obj = parser_interface.object_define()
        self.plots_obj = parser_interface.object_define()
        self.build_plotsinfo()

        self.get_cycles()

    def _get_experiments(self, innov_stats_obj):
        """
        Description
        -----------

        This method parses the innov_stats_obj upon entry and defines
        and returns a list of the experiments contained within.

        Parameters
        ----------

        innov_stats_obj: object

            A Python object containing the innovation statistics
            attributes.

        Returns
        -------

        experiments_list: list

            A Python list containing a list of the experiments
            contained within the innov_stats_obj attribute upon entry.

        """

        # Define and return the list of experiments.
        experiments_list = list(vars(innov_stats_obj))

        return experiments_list

    def _get_regions(self, expt_dict, variable):
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

        InnovStatsPlottingError:

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

    def build_axes(self, plot_obj, ax_obj, **kwargs):
        """
        Description
        -----------

        This method parses the keyword arguments dictionary to
        determine the configuration for of the axes for the respective
        figure.

        Parameters
        ----------

        plot_obj: obj

            A Python matplotlib pyplot object.

        ax_obj:

            A Python matplotlib pyplot subplot object

        Keywords
        --------

        kwargs: dict, optional

            All keyword arguments to modify the default figure axes
            should be based via the dictionary upon entry.

        Raises
        ------

        InnovStatsPlottingError:

            * raised if the axes attribute is NoneType.

        Notes
        -----

        A list of compliant keyword arguments may be found here:

        https://matplotlib.org/stable/api/axes_api.html#matplotlib-axes

        Additional compliant arguments are as follows.

        xint: float/int

            A Python float or integer value specifying the interval at
            which to compute the figure x-axis tick marks; if present
            the x-axis tick marks will be computed within this method.

        yint: float/int

            A Python float or integer value specifying the interval at
            which to compute the figure y-axis tick marks; if present
            the y-axis tick marks will be computed within this method.

        """

        # Parse the keyword argument dictionary; proceed accordingly.
        axes_obj = parser_interface.object_define()
        if kwargs:
            if 'axes' in kwargs.keys():

                # Define the axes attributes contained within the key
                # word argument dictionary.
                axes_dict = parser_interface.dict_key_value(
                    dict_in=kwargs, key='axes', force=True)
                if axes_dict is None:
                    msg = ('The axes attributes cannot be NoneType; check '
                           'the user experiment configuration is correctly '
                           'formatted. Aborting!!!')
                    raise PlotsError(msg=msg)

                # Build the axes according to the available
                # attributes.
                for key in axes_dict.keys():

                    # Define the available axes attributes.
                    value = parser_interface.dict_key_value(
                        dict_in=axes_dict, key=key, no_split=True)
                    axes_obj = parser_interface.object_setattr(
                        object_in=axes_obj, key=key, value=value)

                # Build the x-axis attributes accordingly.
                if any(parser_interface.object_hasattr(
                        object_in=axes_obj, key=key) for key in
                       ['xmin', 'xmax']):
                    if 'xbuffer' in vars(axes_obj) and \
                       axes_obj.xbuffer:
                        plot_obj.gca().set_xlim([axes_obj.xmin - 1.0,
                                                 axes_obj.xmax + 1.0])
                    else:
                        plot_obj.gca().set_xlim([axes_obj.xmin,
                                                 axes_obj.xmax])
                    if 'xint' in axes_dict.keys():
                        xticks = numpy.arange(axes_obj.xmin,
                                              axes_obj.xmax + 1.e-6,
                                              axes_obj.xint)
                        plot_obj.xticks(xticks)
                        if 'xlabels' in axes_dict.keys():
                            xlabels_kwargs = dict()
                            for attr in ['xfontsize', 'xrotation']:
                                if attr in vars(axes_obj):
                                    xlabels_kwargs[attr[1:]] = \
                                        parser_interface.object_getattr(
                                            object_in=axes_obj, key=attr)
                            plot_obj.xlim([0, len(axes_obj.xlabels) + 1])
                            xlabels = [i for i in axes_obj.xlabels]
                            ax_obj.set_xticks(numpy.arange(0, axes_obj.xmax,
                                                           axes_obj.xint))
                            ax_obj.set_xticklabels(xlabels, **xlabels_kwargs)

                # Build the y-axis attributes accordingly.
                if all(parser_interface.object_hasattr(
                        object_in=axes_obj, key=key) for key in
                       ['ymin', 'ymax']):
                    if 'ybuffer' in vars(axes_obj) and \
                       axes_obj.ybuffer:
                        plot_obj.gca().set_ylim([axes_obj.ymin - 1.0,
                                                 axes_obj.ymax + 1.0])
                    else:
                        plot_obj.gca().set_ylim([axes_obj.ymin,
                                                 axes_obj.ymax])
                    if 'yint' in axes_dict.keys():
                        yticks = numpy.arange(axes_obj.ymin,
                                              axes_obj.ymax + 1.e-6,
                                              axes_obj.yint)
                        plot_obj.yticks(yticks)
                        if 'ylabels' in axes_dict.keys():
                            ylabels_kwargs = dict()
                            for attr in ['yfontsize', 'yrotation']:
                                if attr in vars(axes_obj):
                                    ylabels_kwargs[attr[1:]] = \
                                        parser_interface.object_getattr(
                                            object_in=axes_obj, key=attr)
                            plot_obj.ylim([0, len(axes_obj.ylabels) + 1])
                            ylabels = [i for i in axes_obj.ylabels]
                            ax_obj.set_yticks(numpy.arange(0, axes_obj.ymax,
                                                           axes_obj.yint))
                            ax_obj.set_yticklabels(ylabels, **ylabels_kwargs)

            else:

                # Proceed without build the axis object.
                msg = ('The axis attributes have not been specified; the '
                       'figure axes will not be modified.')

        if not kwargs:

            # Proceed without build the axis object.
            msg = ('The axis attributes have not been specified; the '
                   'figure axes will not be modified.')
            self.logger.warn(msg=msg)

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

        InnovStatsPlottingError:

            * raised if an exception is encountered while attempting
              to define the respective color levels.

        """

        # Define the color levels relative to the specified
        # attributes.
        try:
            data_arr = numpy.linspace(0, 1.0, nlevs)
            colorlevs = colorlist.gen_color_normalized(
                cmap=cmap, data_arr=data_arr)
        except Exception as error:
            msg = ('The normalized color levels generation failed '
                   'with error {0}. Aborting!!!'.format(error))
            raise PlotsError(msg=msg)

        return colorlevs

    def build_figure(self, suppress_clf=True):
        """
        Description
        -----------

        This method initializes the attributes for the respective
        matplotlib pyplot object and defines the base-class attributes
        plot_obj, fig_obj, and ax_obj.

        Keywords
        --------

        suppress_clf: bool, optional

            A Python boolean variable specifying whether to suppress
            the closing and clearing of a matplotlib.pyplot figure
            object prior to defining a new one.

        Returns
        -------

        plot_obj: obj

            A Python matplotlib pyplot object.

        fig_obj:

            A Python matplotlib pyplot figure object.

        ax_obj:

            A Python matplotlib pyplot subplot object.

        """

        # Initialize the matplotlib pyplot object and define the
        # respective attributes and proceed accordingly.
        msg = ('Initializing new pyplot objects.')
        self.logger.warn(msg=msg)
        plot_obj = matplotlib.pyplot
        if suppress_clf:
            try:
                plot_obj.close('all')
                plot_obj.clf()
            except Exception:
                pass
        fig_obj = plot_obj.figure()
        ax_obj = plot_obj.subplot()

        return (plot_obj, fig_obj, ax_obj)

    def build_hlines(self, plot_obj, **kwargs):
        """
        Description
        -----------

        This method parses the keyword arguments dictionary to
        determine the configuration for the horizontal lines (hlines)
        for the respective figure.

        Parameters
        ----------

        plot_obj: obj

            A Python matplotlib pyplot object.

        Keywords
        --------

        kwargs: dict, optional

            All keyword arguments to modify the default hlines should
            be based via the dictionary upon entry.

        Notes
        -----

        A list of compliant keyword arguments may be found here:

        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.hlines.html#matplotlib-pyplot-hlines

        """

        # Parse the keyword argument dictionary; proceed accordingly.
        if kwargs is not None:
            hlines_kwargs = parser_interface.dict_key_value(
                dict_in=kwargs, key='hlines', force=True, no_split=True)
            if hlines_kwargs is None:
                msg = ('No horizontal line attributes were found in the user '
                       'experiment configuration; horizontal lines will not be '
                       'plotted.')
                self.logger.warn(msg=msg)

            if hlines_kwargs is not None:
                plot_obj.hlines(**hlines_kwargs)

    def build_legend(self, plot_obj, **kwargs):
        """
        Description
        -----------

        This method parses the keyword arguments dictionary to
        determine the configuration for the legend for the respective
        figure.

        Parameters
        ----------

        plot_obj: obj

            A Python matplotlib pyplot object.

        Keywords
        --------

        kwargs: dict, optional

            All keyword arguments to modify the default figure legend
            should be based via the dictionary upon entry.

        Notes
        -----

        A list of compliant keyword arguments may be found here:

        https://matplotlib.org/stable/api/legend_api.html#module-matplotlib.legend

        """

        # Parse the keyword argument dictionary; proceed accordingly.
        legend_kwargs = parser_interface.dict_key_value(
            dict_in=kwargs, key='legend', force=True, no_split=True)
        if legend_kwargs is None:
            msg = ('No legend attributes were found in the user experiment '
                   'configuration; the legend will not be created.')
            self.logger.warn(msg=msg)

        if legend_kwargs is not None:
            plot_obj.legend(**legend_kwargs)

    def build_plotsinfo(self):
        """
        Description
        -----------

        This method builds the base-class attribute by collecting all
        plotting information for the respective variables and
        respective statistics specified upon entry by vars_list and
        stats_list, respectively.

        Raises
        ------

        InnovStatsPlottingError:

            * raised if the variables attribute cannot be determined
              from the user experiment configuration.

        """

        # Collect the variable attributes from the user experiment
        # configuration.
        variables_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key='variables', force=True,
            no_split=True)
        if variables_dict is None:
            msg = ('The variables attribute could not be determined '
                   'from the user experiment configuration. '
                   'Aborting!!!')
            raise PlotsError(msg=msg)

        # Define the attributes for all variables to be plotted; if
        # not specified reset to the defined default value(s).

        plots_attrs_dict = {'create_images': False}
        for plots_attr in plots_attrs_dict.keys():
            value = parser_interface.dict_key_value(
                dict_in=plots_attrs_dict, key=plots_attr, force=True,
                no_split=True)
            if value is None:
                value = parser_interface.dict_key_value(
                    dict_in=plots_attrs_dict, key=plots_attr, no_split=True)
                msg = ('The plotting attribute {0} has not been '
                       'specified and will be set to {1}.'.format(
                           plots_attr, value))
                self.logger.warn(msg=msg)

            # Loop through the specified variables and statistics
            # values and build the base-class attribute plots_obj.
            for var in self.vars_list:
                stat_obj = parser_interface.object_define()
                for stat in self.stats_list:
                    value = parser_interface.dict_key_value(
                        dict_in=variables_dict, key=stat, force=True,
                        no_split=True)
                    if value is None:
                        msg = ('The statistics attribute {0} for variable '
                               '{1} has not been specified within the user '
                               'experiment configuration; it will be reset '
                               'to NoneType.'.format(stat, var))
                        self.logger.warn(msg=msg)
                        stat_obj = None
                    if value is not None:
                        stat_obj = parser_interface.object_setattr(
                            object_in=stat_obj, key=stat, value=value)
                self.plots_obj = parser_interface.object_setattr(
                    object_in=self.plots_obj, key=var, value=stat_obj)

    def build_title(self, plot_obj, **kwargs):
        """
        Description
        -----------

        This method parses the keyword arguments dictionary to
        determine the configuration for the title for the respective
        figure.

        Parameters
        ----------

        plot_obj: obj

            A Python matplotlib pyplot object.

        Keywords
        --------

        kwargs: dict, optional

            All keyword arguments to modify the default figure title
            should be based via the dictionary upon entry.

        Notes
        -----

        A list of compliant keyword arguments may be found here:

        https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.title.html#matplotlib-pyplot-title

        """

        # Parse the keyword argument dictionary; proceed accordingly.
        if kwargs is not None:
            try:
                plot_obj.title(**kwargs)
            except TypeError:
                try:
                    kwargs['s'] = kwargs.pop('label')

                except KeyError:
                    kwargs['label'] = kwargs.pop('s')

                plot_obj.title(**kwargs)

    def build_varinfo(self, variable):
        """
        Description
        -----------

        This method returns the attributes for the respective analysis
        variable contained within the user experiment configuration.

        Parameters
        ----------

        variable: str

            A Python string specifying the respective analysis
            variable for which to collect the attributes within the
            user experiment configuration.

        Returns
        -------

        varinfo_dict: dict

            A Python dictionary containing the respective analysis
            variable attributes collected from the user experiment
            configuration.

        Raises
        ------

        InnovStatsPlottingError:

            * raised if the user experiment configuration does not
              contain the 'variables' attribute.

            * raised if the attributes for the respective analysis
              variable cannot be determined from the user experiment
              configuration.

        """

        # Parse the user experiment configuration for the respective
        # analysis variable and proceed accordingly.
        variable_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key='variables', force=True,
            no_split=True)
        if variable_dict is None:
            msg = ('The variables attribute could not be determined from '
                   'the user experiment configuration. Aborting!!!')
            raise PlotsError(msg=msg)

        varinfo_dict = parser_interface.dict_key_value(
            dict_in=variable_dict, key=variable, force=True, no_split=True)
        if varinfo_dict is None:
            msg = ('The attributes for variable {0} could not be determined '
                   'from the user experiment configuration. Aborting!!!'.
                   format(variable))
            raise PlotsError(msg=msg)

        return varinfo_dict

    def build_vlines(self, plot_obj, **kwargs):
        """
        Description
        -----------

        This method parses the keyword arguments dictionary to
        determine the configuration for the vertical lines (vlines)
        for the respective figure.

        Parameters
        ----------

        plot_obj: obj

            A Python matplotlib pyplot object.

        Keywords
        --------

        kwargs: dict, optional

            All keyword arguments to modify the default vlines should
            be based via the dictionary upon entry.

        Notes
        -----

        A list of compliant keyword arguments may be found here:

        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.vlines.html#matplotlib-pyplot-vlines

        """

        # Parse the keyword argument dictionary; proceed accordingly.
        if kwargs is not None:
            vlines_kwargs = parser_interface.dict_key_value(
                dict_in=kwargs, key='vlines', force=True, no_split=True)
            if vlines_kwargs is None:
                msg = ('No vertical line attributes were found in the user '
                       'experiment configuration; vertical lines will not be '
                       'plotted.')
                self.logger.warn(msg=msg)

            if vlines_kwargs is not None:
                plot_obj.vlines(**vlines_kwargs)

    def build_xlabel(self, ax_obj, **kwargs):
        """
        Description
        -----------

        This method parses the keyword arguments dictionary to
        determine the configuration for the x-label of the respective
        figure.

        Parameters
        ----------

        ax_obj:

            A Python matplotlib pyplot subplot object.

        Keywords
        --------

        kwargs: dict, optional

            All keyword arguments to modify the default figure x-label
            should be based via the dictionary upon entry.

        Notes
        -----

        A list of compliant keyword arguments may be found here:

        https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.set_xlabel.html#matplotlib-axes-axes-set-xlabel

        """

        # Parse the keyword argument dictionary; proceed accordingly.
        labels_kwargs = parser_interface.dict_key_value(
            dict_in=kwargs, key='xlabel', force=True, no_split=True)
        if labels_kwargs is None:
            msg = ('No x-label attributes were found in the user experiment '
                   'experiment configuration; x-labels will not be created.')
            self.logger.warn(msg=msg)

        if labels_kwargs is not None:
            ax_obj.set_xlabel(**labels_kwargs)

    def build_ylabel(self, ax_obj, **kwargs):
        """
        Description
        -----------

        This method parses the keyword arguments dictionary to
        determine the configuration for the y-label of the respective
        figure.

        Parameters
        ----------

        ax_obj:

            A Python matplotlib pyplot subplot object.

        Keywords
        --------

        kwargs: dict, optional

            All keyword arguments to modify the default figure y-label
            should be based via the dictionary upon entry.

        Notes
        -----

        A list of compliant keyword arguments may be found here:

        https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.set_xlabel.html#matplotlib-axes-axes-set-ylabel

        """

        # Parse the keyword argument dictionary; proceed accordingly.
        labels_kwargs = parser_interface.dict_key_value(
            dict_in=kwargs, key='ylabel', force=True, no_split=True)
        if labels_kwargs is None:
            msg = ('No y-label attributes were found in the user experiment '
                   'experiment configuration; y-labels will not be created.')
            self.logger.warn(msg=msg)

        if labels_kwargs is not None:
            ax_obj.set_ylabel(**labels_kwargs)

    def compute(self, filelist_obj_dict, regions_obj):
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

        innov_stats_obj: obj

            A Python object containing the temporal mean analysis
            variable innovation statistics as a function of user
            specified experiment(s) and region(s) of interest.

        Raises
        ------

        OMFProfileError:

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
        (expts_list, regions_list) = (list(filelist_obj_dict.keys()),
                                      list(vars(regions_obj).keys()))

        # Loop through the respective experiments and proceed with the
        # temporal mean innovation statistic computation as a function
        # region and variable.
        innov_stats_obj = parser_interface.object_define()
        for expt in expts_list:
            innov_stats_obj = parser_interface.object_setattr(
                object_in=innov_stats_obj, key=expt, value=None)
            expt_dict = parser_interface.object_getattr(
                object_in=self.expts_obj, key=expt, force=True)
            if expt_dict is None:
                msg = ('The input filelist dictionary (filelist_obj_dict) '
                       'does not contain the attributes for experiment {0} '
                       'upon entry. Aborting!!!'.format(expt))
                raise OMFProfileError(msg=msg)
            datapath = parser_interface.dict_key_value(
                dict_in=expt_dict, key='datapath', force=True, no_split=True)
            if datapath is None:
                msg = ('The datapath attribute for experiment {0} could '
                       'not be determined from the user experiment '
                       'configuration. Aborting!!!.'.format(expt))
                raise OMFProfileError(msg=msg)

            # Loop through the respective analysis variables and
            # compute the respective innovation statistics attributes.
            region_dict = dict()
            for variable in self.variables:

                # Collect the base filenames for the respective
                # experiment and analysis variable.
                filelist = parser_interface.object_getattr(
                    object_in=filelist_obj_dict[expt], key=variable, force=True)
                if filelist is None:
                    msg = ('The innovation statistics base filenames could '
                           'not be determined for experiment {0} and analysis '
                           'variable {1}. Aborting!!!'.format(expt, variable))
                    raise OMFProfileError(msg=msg)

                # Define and collect the netCDF attributes for the
                # respective experiment and analysis variables.
                ncfile = os.path.join(datapath, filelist[0])
                exist = fileio_interface.fileexist(path=ncfile)
                if not exist:
                    ncfile = filelist[0]

                # Define the levels attributes.
                nlevs = netcdf4_interface.ncreaddim(
                    ncfile=ncfile, ncdimname=self.vcoord_name)
                levs = netcdf4_interface.ncreadvar(ncfile=ncfile,
                                                   ncvarname=self.vcoord_name)

                # Loop through the specified regions and
                # compute/define the respective innovation statistic
                # attributes.
                region_dict[variable] = dict()
                for region in vars(self.regions_obj):
                    (stats_dict, stats_mean_dict) = [dict() for i in range(2)]
                    for stat_var in self.stats_vars_list:
                        stats_dict[stat_var] = numpy.zeros((len(filelist), nlevs),
                                                           dtype=numpy.float64)
                    for (i, filename) in enumerate(filelist):
                        ncfile = os.path.join(datapath, filename)
                        exist = fileio_interface.fileexist(path=ncfile)
                        if not exist:
                            ncfile = filename
                        for stat_var in self.stats_vars_list:
                            ncvarname = '{0}_{1}'.format(stat_var, region)
                            stats_dict[stat_var][i, :] = netcdf4_interface.ncreadvar(
                                ncfile=ncfile, ncvarname=ncvarname)

                    # Compute the respective temporal mean innovation
                    # statistics.
                    for stat_var in self.stats_vars_list:
                        msg = ('Computing innovation statistic {0} for variable '
                               '{1} of region {2} for experiment {3}.'.
                               format(stat_var, variable, region, expt))
                        self.logger.info(msg=msg)
                        stats_mean_dict[stat_var] = numpy.nanmean(
                            stats_dict[stat_var], axis=0)
                    region_dict[variable][region] = stats_mean_dict
                innov_stats_obj = parser_interface.object_setattr(
                    object_in=innov_stats_obj, key=expt, value=region_dict)

            # Append level attribute to innovation statistics object.
            innov_stats_obj = parser_interface.object_setattr(
                object_in=innov_stats_obj, key='levels', value=levs)

        return innov_stats_obj

    def create_image(self, plot_obj, savename, **kwargs):
        """
        Description
        -----------

        This method creates an external file containing the respective
        figure.

        Parameters
        ----------

        plot_obj: obj

            A Python matplotlib pyplot object.

        savename: str

            A Python string specifying the path to the figure image to
            be created.

        Keywords
        --------

        kwargs: dict, optional

            All keyword arguments attributes for the figure creation
            should be based via the dictionary upon entry.

        Notes
        -----

        A list of compliant keyword arguments may be found here:

        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.savefig.html#matplotlib-pyplot-savefig

        """

        # Parse the keyword argument dictionary; proceed accordingly.
        msg = ('Creating figure {0}.'.format(savename))
        self.logger.info(msg=msg)
        savefig_kwargs = parser_interface.dict_key_value(
            dict_in=kwargs, key='savefig', force=True, no_split=True)

        if savefig_kwargs is None:
            savefig_kwargs = dict()
        plot_obj.savefig(fname=savename, **savefig_kwargs)

    def get_cycles(self):
        """
        Description
        -----------

        This method collects the forecast cycle information from the
        user experiment configuration and defines the base-class
        attribute cycle_obj.

        Raises
        ------

        InnovStatsError:

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
            dict_in=self.yaml_dict, key='cycles', force=True,
            no_split=True)

        if cycles_dict is None:
            msg = ('The attribute cycles could not be determined from '
                   'the user experiment configuration. Aborting!!!')
            raise InnovStatsError(msg=msg)

        # Construct the base-class attribute cycles_obj using the
        # information specified in the user experiment configuration.
        cycle_attrs_list = ['cycle_interval',
                            'cycle_start',
                            'cycle_stop'
                            ]
        for cycle_attr in cycle_attrs_list:
            value = parser_interface.dict_key_value(
                dict_in=cycles_dict, key=cycle_attr, force=True,
                no_split=True)
            if value is None:
                msg = ('The cycles attribute {0} could not be determined '
                       'from the user experiment configuration. '
                       'Aborting!!!'.format(cycle_attr))
                raise InnovStatsError(msg=msg)

            self.cycles_obj = parser_interface.object_setattr(
                object_in=self.cycles_obj, key=cycle_attr, value=value)

    def get_experiments(self):
        """
        Description
        -----------

        This method collects the information for the experiments for
        which to computed the innovation statistics from the user
        experiment configuration and defines the base-class attribute
        expts_obj.

        Raises
        ------

        InnovStatsError:

            * raised if the user experiment configuration does not
              contain the attribute experiments.

        """

        # Initialize the base-class expts_obj attribute and collect
        # the experiment information from the user experiment
        # configuration.
        self.expts_obj = parser_interface.object_define()
        expts_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key='experiments', force=True,
            no_split=True)
        if expts_dict is None:
            msg = ('The attribute experiments could not be determined '
                   'from the user experiment configuration. Aborting!!!')
            raise InnovStatsError(msg=msg)

        # Construct the base-class attribute expts_obj using the
        # information specified in the user experiment configuration.
        for expt in expts_dict.keys():
            value = parser_interface.dict_key_value(
                dict_in=expts_dict, key=expt, no_split=True)
            self.expts_obj = parser_interface.object_setattr(
                object_in=self.expts_obj, key=expt, value=value)

    def get_regions(self):
        """
        Description
        -----------

        This method collects the information for the analysis regions
        to be plotted from the user experiment configuration and
        defines the base-class attribute regions_obj.

        Raises
        ------

        InnovStatsError:

            * raised if the user experiment configuration does not
              contain the attribute regions.

        """

        # Initialize the base-class regions_obj attribute and collect
        # the regions information from the user experiment
        # configuration.
        regions_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key='regions', force=True,
            no_split=True)

        # Construct the base-class attribute regions_obj accordingly.
        if regions_dict is None:
            self.regions_obj = None
            msg = ('The attribute regions is not specified or could not '
                   'be determined from the user experiment configuration.')
            self.logger.warn(msg=msg)
        if regions_dict is not None:
            self.regions_obj = parser_interface.object_define()
            for region in regions_dict.keys():
                value = parser_interface.dict_key_value(
                    dict_in=regions_dict, key=region, no_split=True)
                self.regions_obj = parser_interface.object_setattr(
                    object_in=self.regions_obj, key=region, value=value)

    def get_times(self, expt_name: str) -> None:
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

        InnovStatsError:

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
            object_in=self.expts_obj, key=expt_name, force=True)

        if expt_dict is None:
            msg = ('The attributes for experiment {0} could not '
                   'be determined from the user experiment '
                   'configuration. Aborting!!!'.format(expt_name))
            error(msg=msg)

        # Define the start and stop analysis times using the
        # attributes within the base-class attribute cycles_obj.
        try:
            timestamp = datetime_interface.datestrupdate(
                datestr=self.cycles_obj.cycle_start,
                in_frmttyp=self.cyclestr_infrmttyp,
                out_frmttyp=self.cyclestrfrmt,
                offset_seconds=self.cycles_obj.cycle_interval)
            last_timestamp = datetime_interface.datestrupdate(
                datestr=self.cycles_obj.cycle_stop,
                in_frmttyp=self.cyclestr_infrmttyp,
                out_frmttyp=self.cyclestrfrmt)

        except Exception as errmsg:
            msg = ('Defining the initial and final timestamps failed with error '
                   f'{errmsg}. Aborting!!!')
            error(msg=msg)

        # Create a list of valid timestamps for the respective
        # experiment name provided upon entry.
        timestamp_list = list()
        msg = 'The valid timestamps for experiment {0} are the following:\n\n'.\
            format(expt_name)
        while timestamp <= last_timestamp:
            timestamp_list.append(timestamp)
            timestamp = datetime_interface.datestrupdate(
                datestr=timestamp, in_frmttyp=self.cyclestrfrmt,
                out_frmttyp=self.cyclestrfrmt, offset_seconds=self.cycles_obj.cycle_interval)
            msg = msg + '{0}\n'.format(timestamp)
        self.logger.warn(msg=msg)

        return timestamp_list

# ----


@msg_except_handle(InnovStatsPlottingError)
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
