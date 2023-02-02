# =========================================================================

# Module: ush/innov_stats/plot/plots/omf_profile_plot.py

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

    omf_profile_plot.py

Description
-----------

    This module contains the base-class for OMF profile plotting for
    the respective supported applications.

Classes
-------

    OMFProfilePlots(options_obj, basedir, vars_list=None, stats_list=None,
                    is_gsi_atmos=False, is_soca_ice=False,
                    is_soca_ocean=False)

        This is the base-class object for all OMF profile plots; it is
        a sub-class of OMFProfile.

Author(s)
---------

    Henry R. Winterbottom; 31 January 2023

History
-------

    2023-01-31: Henry Winterbottom -- Initial implementation.

"""

# ----

# pylint: disable=broad-except
# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals

# ----

from typing import Tuple

from innov_stats.plot import error
from innov_stats.plot.omf_profile import OMFProfile
from tools import parser_interface

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


class OMFProfilePlots(OMFProfile):
    """
    Description
    -----------

    This is the base-class object for all OMF profile plots; it is a
    sub-class of OMFProfile.

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
          innovation statistics plots method are not valid.

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

        Creates a new OMFProfilePlots object.

        """

        # Define the base-class attributes.
        super().__init__(
            options_obj=options_obj,
            basedir=basedir,
            is_gsi_atmos=is_gsi_atmos,
            is_soca_ice=is_soca_ice,
            is_soca_ocean=is_soca_ocean,
        )

        # Check that the constructor arguments are valid; proceed
        # accordingly.
        if not any([is_gsi_atmos, is_soca_ice, is_soca_ocean]):
            msg = (
                "The supported innovation statistics application "
                "has not been specified. Aborting!!!"
            )
            error(msg=msg)

        (self.is_gsi_atmos, self.is_soca_ice, self.is_soca_ocean) = (
            is_gsi_atmos,
            is_soca_ice,
            is_soca_ocean,
        )

        # Define the respective experiment plotting attributes.
        self.expt_plot_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key="experiments", force=True, no_split=True
        )
        if self.expt_plot_dict is None:
            msg = (
                "The respective experiment plotting attributes "
                "could not be determined from the user experiment "
                "configuration. Aborting!!!"
            )
            error(msg=msg)

    def __plotattrs__(
        self, plot_obj: object, ax_obj: object, region: str, stat: str, variable: str
    ) -> object:
        """
        Description
        -----------

        This method plots the figure attributes; this method should be
        called following the plotting of the respective innovation
        statistic (see __plotstat__).

        Parameters
        ----------

        plot_obj: object

            A Python matplotlib pyplot object containing the
            respective figure attributes.

        ax_obj: object

            A Python matplotlib pyplot axis-level subplot object

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

            A Python matplotlib pyplot object containing the updated
            figure attributes.

        """

        # Collect the plotting attributes for the respective variable.
        varinfo_dict = self.plotter.build_varinfo(variable=variable)

        # Customize the figure aspects based on the user
        # experiment configuration and proceed
        # accordingly.
        plot_attrs_dict = parser_interface.dict_key_value(
            dict_in=varinfo_dict, key=stat, force=True
        )
        plot_regions_title_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict["regions"][region],
            key="title",
            force=True,
        )

        # Build the figure axes and legend.
        if plot_attrs_dict is not None:
            self.plotter.build_axes(
                plot_obj=plot_obj, ax_obj=ax_obj, **plot_attrs_dict)
            ax_obj.spines["top"].set_visible(False)
            ax_obj.spines["right"].set_visible(False)

            if any([self.is_gsi_atmos, self.is_soca_ocean]):
                plot_obj.gca().invert_yaxis()
                self.plotter.build_legend(plot_obj=plot_obj, **plot_attrs_dict)
                self.plotter.build_xlabel(ax_obj=ax_obj, **plot_attrs_dict)
                self.plotter.build_ylabel(ax_obj=ax_obj, **plot_attrs_dict)
                self.plotter.build_vlines(plot_obj=plot_obj, **plot_attrs_dict)

        # Build the figure title.
        if plot_regions_title_dict is not None:
            self.plotter.build_title(plot_obj, **plot_regions_title_dict)

        return plot_obj

    def __plotsave__(
        self, plot_obj: object, region: str, stat: str, variable: str
    ) -> None:
        """
        Description
        -----------

        This method saves the attributes of the Python pyplot object
        to an external Portable Network Graphics (PNG) formatted
        image.

        Parameters
        ----------

        plot_obj: object

            A Python matplotlib pyplot object containing the
            respective figure attributes.

        region: str

            A Python string specifiying the region for which to
            plot the respective variable innovations statistics.

        stat: str

            A Python string specifying the innovation statistic to be
            plotted.

        variable: str

            A Python string specifying the variable for which to plot
            the respective innovation statistic.

        """

        # Define the PNG formatted image file name.
        savename = f"{region}.{stat}.{variable}.png"

        # Create the PNG formatted image file accordingly.
        try:
            self.plotter.create_image(
                plot_obj=plot_obj, savename=savename, dpi=500)

        except TypeError:
            self.plotter.create_image(plot_obj=plot_obj, savename=savename)

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
