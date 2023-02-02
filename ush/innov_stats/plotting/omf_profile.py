# =========================================================================

# $$$ MODULE DOCUMENTATION BLOCK

# UFS-RNR-analysis :: apps/innov_stats/py/plots/omf_profile.py

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

    This module contains the base-class for all OMF profile plotting.

Classes
-------

    OMFProfilePlots(options_obj, basedir, vars_list=None, stats_list=None)

        This is the base-class object for all OMF profile plots; it is
        a sub-class of Plots.

Author(s)
---------

    Henry R. Winterbottom; 29 January 2023

History
-------

    2023-01-29: Henry Winterbottom -- Initial implementation.

"""

# ----

from innov_stats.plotting import Plots, error
from tools import parser_interface

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


class OMFProfilePlots(Plots):
    """
    Description
    -----------

    This is the base-class object for all OMF profile plots; it is a
    sub-class of Plots.

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

    is_gsi: bool, optional

        A Python boolean valued variables specifying whether the
        methods relative to the GSI application are to be implemented.

    is_soca: bool, optional

        A Python boolean valued variables specifying whether the
        methods relative to the SOCA application are to be
        implemented.

    Raises
    ------

    InnovStatsPlottingError:

        * raised if the argument(s) pertaining the specified
          innovation statistics plots method are not valid.

    """

    def __init__(self, options_obj: object, basedir: str,
                 vars_list: list = None,
                 stats_list: list = None, is_gsi_atmos: bool = False,
                 is_soca_ice: bool = False, is_soca_ocean: bool = False):
        """
        Description
        -----------

        Creates a new OMFProfilePlots object.

        """

        # Define the base-class attributes.
        super().__init__(options_obj=options_obj, basedir=basedir,
                         vars_list=vars_list, stats_list=stats_list)

        if is_soca_ice:
            msg = ("OMF profiles are currently not supported for SOCA "
                   "ice applications. Aborting!!!")
            error(msg=msg)

        # Define the experiment attributes.
        self.expt_plot_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key='experiments', force=True,
            no_split=True)
        if self.expt_plot_dict is None:
            msg = ('The respective experiment plotting attributes '
                   'could not be determined from the user experiment '
                   'configuration. Aborting!!!')
            error(msg=msg)

    def __plot_attrs(self, plot_obj: object, fig_obj: object, ax_obj: object):
        """ """

        # Customize the figure aspects based on the user experiment
        # configuration and proceed accordingly.
        plot_attrs_dict = parser_interface.dict_key_value(
            dict_in=varinfo_dict, key=stat, force=True)

        # Build the figure axes and legend.
        if plot_attrs_dict is not None:
            self.build_axes(plot_obj=plot_obj,
                            ax_obj=ax_obj, **plot_attrs_dict)
        ax_obj.spines['top'].set_visible(False)
        ax_obj.spines['right'].set_visible(False)
        if any([self.is_gsi_atmos, self.is_soca_ocean]):
            plot_obj.gca().invert_yaxis()

        self.build_xlabel(ax_obj=ax_obj, **plot_attrs_dict)
        self.build_ylabel(ax_obj=ax_obj, **plot_attrs_dict)
        self.build_vlines(plot_obj=plot_obj, **plot_attrs_dict)

        self.build_legend(plot_obj=plot_obj, **plot_attrs_dict)

        # Build the figure title.
        plot_regions_title_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict['regions'][region], key='title', force=True)
        if plot_regions_title_dict is not None:
            self.build_title(plot_obj, **plot_regions_title_dict)

    def __plot_expt__(self, plot_obj: object, innov_stats_obj: object, region: str,
                      variable: str, experiment: str) -> object:
        """
        """

        # Collect the respective experiment plotting keyword
        # arguments.
        plot_kwargs = parser_interface.dict_key_value(
            dict_in=self.expt_plot_dict[experiment],
            key='plotting', force=True)

        # Collect all attributes for the respective region, variable,
        # experiment, and innovation statistic; proceed accordingly.
        expt_dict = parser_interface.object_getattr(
            object_in=innov_stats_obj, key=experiment)
        var_dict = parser_interface.dict_key_value(
            dict_in=expt_dict, key=variable, force=True,
            no_split=True)

        if var_dict is None:
            msg = (f'The variable {variable} has not been define for '
                   f'experiment {experiment}. Aborting!!!')
            error(msg=msg)

        # Collect the region attributes.
        region_dict = parser_interface.dict_key_value(
            dict_in=var_dict, key=region, force=True, no_split=True)
        if region_dict is None:
            msg = ('The region {region} has not be defined for variable '
                   '{variable} for experiment {experiment}. Aborting!!!')
            error(msg=msg)

        statvals = parser_interface.dict_key_value(
            dict_in=region_dict, key=stat, force=True, no_split=True)
        if statvals is None:
            msg = (f'The statistic {stat} values could not be determined '
                   f'for variable {variable} experiment {experiment}. Aborting!!!')
            error(msg=msg)

        # Plot the respective innovation statistic for the respective
        # region, variable, and experiment.
        if plot_kwargs is None:
            msg = (f'Plotting arguments for innovation statistic {stat} '
                   f'for variable {variable} of experiment {experiment} could not be '
                   'determined from the user experiment configuration; '
                   'plotting will not be performed.')
            self.logger.warn(msg=msg)

        # Plot the respective innovation statistics for the respective
        # experiment and region.
        try:
            plot_obj.plot(
                statvals, innov_stats_obj.levels, **plot_kwargs)
            stat_varinfo_dict = parser_interface.dict_key_value(
                dict_in=varinfo_dict, key=stat, force=True)

        except Exception as errmsg:
            msg = (f'The plotting of innovation statistic {stat} for experiment '
                   f'{experiment} variable {variable} in region {region} failed '
                   'with error {errmsg}. Aborting!!!')
            error(msg=msg)

        return plot_obj

    def __plot_stats__(self, innov_stats_obj: object, region: str, variable: str) -> object:
        """ """

        # Define the plotting object attributes and proceed
        # accordingly.
        (plot_obj, fig_obj, ax_obj) = self.build_figure()
        varinfo_dict = self.build_varinfo(variable=variable)

        # Loop through the respective experiments.
        for experiment in experiments_list:

            # Plot the innovations statistics profile for the
            # respective experiment.
            plot_obj = self.__plot_expt__(plot_obj=plot_obj,
                                          innov_stats_obj=innov_stat_obj, region=region,
                                          variable=variable, experiment=experiment)

        return (plot_obj, fig_obj, ax_obj)

    def get_omf_filelist(self, expt_name):
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

        filelist_obj: obj

            A Python object containing the list of existing innovation
            statistics files available to the respective experiment
            for each specified (i.e., allocate) analysis variable.

        Raises
        ------

        OMFProfileError:

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
            object_in=self.expts_obj, key=expt_name, force=True)
        if expt_dict is None:
            msg = ('The attributes for experiment {0} could not be '
                   'determined from the user experiment configuration. '
                   'Aborting!!!.'.format(expt_name))
            raise OMFProfileError(msg=msg)
        datapath = parser_interface.dict_key_value(
            dict_in=expt_dict, key='datapath', force=True, no_split=True)

        if datapath is None:
            msg = ('The attribute datapath could not be determined for '
                   'experiment {0}. Aborting!!!'.format(expt_name))
            raise OMFProfileError(msg=msg)

        # Collect the variables attribute for the respective
        # experiment.
        vars_dict = parser_interface.dict_key_value(
            dict_in=expt_dict, key='variables', force=True)
        if vars_dict is None:
            msg = ('The variables attribute could not be determined '
                   'for experiment {0} within the user experiment '
                   'configuration. Aborting!!!'.format(expt_name))
            raise OMFProfileError(msg=msg)

        # Collect the file naming attributes from the user experiment
        # configuration for the respective experiment(s) and
        # variables.
        filelist_obj = parser_interface.object_define()
        filelist_dict = {'filename': None, 'filename_offset_seconds':
                         0}
        for var in vars_dict.keys():
            var_dict = parser_interface.dict_key_value(
                dict_in=vars_dict, key=var)
            if var_dict is None:
                msg = ('The variable attributes for variable {0} of '
                       'experiment {1} could not be determined from '
                       'the user experiment configuration. Aborting!!!'.
                       format(var, expt_name))
                raise OMFProfileError(msg=msg)
            for item in filelist_dict.keys():
                value = parser_interface.dict_key_value(
                    dict_in=var_dict, key=item, force=True, no_split=True)
                if value is None:
                    value = parser_interface.dict_key_value(
                        dict_in=filelist_dict, key=item, no_split=True)

                    if value is None:
                        msg = ('The attribute {0} for experiment {1} '
                               'and variable {2} cannot be NoneType. '
                               'Aborting!!!')
                        raise OMFProfileError(msg=msg)
                    if value is not None:
                        filelist_dict[item] = value
                if value is not None:
                    filelist_dict[item] = value

            # Build the file list for the respective experiment and
            # variables; append only existing file paths.
            filelist = list()
            (filename, offset_seconds) = (filelist_dict['filename'],
                                          filelist_dict['filename_offset_seconds'])
            filepath = os.path.join(datapath, filename)
            for timestamp in timestamp_list:
                filename = datetime_interface.datestrupdate(
                    datestr=timestamp, in_frmttyp=self.cyclestrfrmt,
                    out_frmttyp=filepath, offset_seconds=offset_seconds)
                exist = fileio_interface.fileexist(path=filename)
                if exist:
                    msg = ('File {0} exists for experiment {1} variable {2}.'.
                           format(filename, expt_name, var))
                    self.logger.info(msg=msg)
                    filelist.append(filename)
                if not exist:
                    msg = ('File {0} does not exist for experiment {1} '
                           'variable {2} and will not be used for creation of '
                           'the diagnostics plots.'.format(filename, expt_name,
                                                           var))
                    self.logger.warn(msg=msg)
            filelist_obj = parser_interface.object_setattr(object_in=filelist_obj,
                                                           key=var, value=filelist)
        return filelist_obj

    def plot(self, regions_obj, innov_stats_obj):
        """
        Description
        -----------

        This method plots the OMF profiles as a function of region and
        a respective analysis variable's innovation statistic(s).

        Parameters
        ----------

        regions_obj: obj

            A Python object containing the analysis region(s) to be
            plotted from the user experiment configuration.

        innov_stats_obj: obj

            A Python object containing the respective computed
            temporal mean innovation statistics for each experiment
            and respective region(s).

        Raises
        ------

        InnovStatsPlottingError:

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

        # Collect the innovations statistics attributes.
        experiments_list = self._get_experiments(
            innov_stats_obj=innov_stats_obj)
        if 'levels' in experiments_list:
            experiments_list.remove('levels')

        expt_dict = parser_interface.object_getattr(
            object_in=innov_stats_obj, key=experiments_list[0])
        regions_list = self.get_regions(expt_dict=expt_dict,
                                        variable=self.vars_list[0])

        print(expt_dict)
        quit()

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
                for stat in self.stats_list:

                    # Define the plotting object attributes and
                    # proceed accordingly.
                    (plot_obj, fig_obj, ax_obj) = \
                        self.__plot_stats__(
                            innov_stats_obj=innov_stats_obj, variable=variable)
                    varinfo_dict = self.build_varinfo(variable=variable)

                    # Build the output figure name and create the
                    # figure.
                    savename = f'{region}.{stat}.{variable}.png'
                    try:
                        self.create_image(plot_obj=plot_obj, savename=savename,
                                          dpi=500)

                    except TypeError:
                        self.create_image(plot_obj=plot_obj, savename=savename)
