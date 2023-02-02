# =========================================================================

# Module: ush/innov_stats/plot/plots/__init__.py

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

Author(s)
---------

    Henry R. Winterbottom; 31 January 2023

History
-------

    2023-01-31: Henry Winterbottom -- Initial implementation.

"""

# ----

# pylint: disable=broad-except
# pylint: disable=too-many-lines
# pylint: disable=too-many-nested-blocks

# ----

from typing import Tuple

import matplotlib.pyplot
import mycolorpy
import numpy
from confs.yaml_interface import YAML
from innov_stats.plot import error
from mycolorpy import colorlist
from tools import fileio_interface, parser_interface
from utils.error_interface import Error
from utils.logger_interface import Logger

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

    options_obj: object

        A Python object containing the command line argument
        attributes.

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

    InnovStatsPlotError:

        * raised if either vars_list or stats_list are NoneType upon
          entry.

    """

    def __init__(
        self, options_obj: object, vars_list: list = None, stats_list: list = None
    ):
        """
        Description
        -----------

        Creates a new Plots object.

        """

        # Define the base-class attributes.
        self.logger = Logger()
        self.options_obj = options_obj
        self.yaml_dict = YAML().read_yaml(yaml_file=self.options_obj.yaml_file)

        (self.stats_list, self.vars_list) = (stats_list, vars_list)
        if any([self.stats_list, self.vars_list]) is None:
            msg = (
                "The attributes stats_list and vars_list must "
                "not be NoneType on entry; received the following "
                f"upon entry:\n stats_list={self.stats_list}\n "
                f"vars_list={self.vars_list}\n Aborting!!!"
            )
            error(msg=msg)

        self.build_plotsinfo()
        self.pyplot = matplotlib.pyplot

    def __xaxes__(
        self, plot_obj: object, ax_obj: object, axes_obj: object, axes_dict: dict
    ) -> Tuple[object, object]:
        """
        Description
        -----------

        This method builds the x-axis attributes accordingly.

        Parameters
        ----------

        plot_obj: object

            A Python matplotlib pyplot object.

        ax_obj: object

            A Python matplotlib pyplot subplot object.

        axes_obj: object

            A Python object containing defined x-axis attributes; see
            method build_axes.

        axes_dict: dict

            A Python dictionary containing experiment configuration
            the x-axis attributes; see method build_axes.

        Returns
        -------

        plot_obj: object

            An updated Python matplotlib pyplot object.

        ax_obj: object

            An updated Python matplotlib pyplot axis-level subplot
            object.

        """

        # Build the x-axis attributes accordingly.
        if any(
            parser_interface.object_hasattr(object_in=axes_obj, key=key)
            for key in ["xmin", "xmax"]
        ):

            if "xbuffer" in vars(axes_obj) and axes_obj.xbuffer:
                plot_obj.gca().set_xlim(
                    [axes_obj.xmin - 1.0, axes_obj.xmax + 1.0])

            else:
                plot_obj.gca().set_xlim([axes_obj.xmin, axes_obj.xmax])

                if "xint" in axes_dict.keys():
                    xticks = numpy.arange(
                        axes_obj.xmin, axes_obj.xmax + 1.0e-6, axes_obj.xint
                    )
                    plot_obj.xticks(xticks)

                    if "xlabels" in axes_dict.keys():
                        xlabels_kwargs = {}
                        for attr in ["xfontsize", "xrotation"]:
                            if attr in vars(axes_obj):
                                xlabels_kwargs[
                                    attr[1:]
                                ] = parser_interface.object_getattr(
                                    object_in=axes_obj, key=attr
                                )
                        plot_obj.xlim([0, len(axes_obj.xlabels) + 1])
                        xlabels = list(axes_obj.xlabels)
                        ax_obj.set_xticks(numpy.arange(
                            0, axes_obj.xmax, axes_obj.xint))
                        ax_obj.set_xticklabels(xlabels, **xlabels_kwargs)

        return (plot_obj, ax_obj)

    def __yaxes__(
        self, plot_obj: object, ax_obj: object, axes_obj: object, axes_dict: dict
    ) -> Tuple[object, object]:
        """
        Description
        -----------

        This method builds the y-axis attributes accordingly.

        Parameters
        ----------

        plot_obj: object

            A Python matplotlib pyplot object.

        ax_obj: object

            A Python matplotlib pyplot subplot object.

        axes_obj: object

            A Python object containing defined y-axis attributes; see
            method build_axes.

        axes_dict: dict

            A Python dictionary containing experiment configuration
            the y-axis attributes; see method build_axes.

        Returns
        -------

        plot_obj: object

            An updated Python matplotlib pyplot object.

        ax_obj: object

            An updated Python matplotlib pyplot axis-level subplot
            object.

        """

        # Build the y-axis attributes accordingly.
        if all(
            parser_interface.object_hasattr(object_in=axes_obj, key=key)
            for key in ["ymin", "ymax"]
        ):

            if "ybuffer" in vars(axes_obj) and axes_obj.ybuffer:
                plot_obj.gca().set_ylim([axes_obj.ymin - 1.0, axes_obj.ymax + 1.0])

            else:
                plot_obj.gca().set_ylim([axes_obj.ymin, axes_obj.ymax])
            if "yint" in axes_dict.keys():
                yticks = numpy.arange(
                    axes_obj.ymin, axes_obj.ymax + 1.0e-6, axes_obj.yint
                )
                plot_obj.yticks(yticks)
                if "ylabels" in axes_dict.keys():
                    ylabels_kwargs = {}
                    for attr in ["yfontsize", "yrotation"]:
                        if attr in vars(axes_obj):
                            ylabels_kwargs[attr[1:]] = parser_interface.object_getattr(
                                object_in=axes_obj, key=attr
                            )
                    plot_obj.ylim([0, len(axes_obj.ylabels) + 1])
                    ylabels = list(axes_obj.ylabels)
                    ax_obj.set_yticks(numpy.arange(0, axes_obj.ymax, axes_obj.yint))
                    ax_obj.set_yticklabels(ylabels, **ylabels_kwargs)

        return (plot_obj, ax_obj)

    def build_axes(self, plot_obj: object, ax_obj: object, **kwargs: dict) -> None:
        """
        Description
        -----------

        This method parses the keyword arguments dictionary to
        determine the configuration for of the axes for the respective
        figure.

        Parameters
        ----------

        plot_obj: object

            A Python matplotlib pyplot object.

        ax_obj: object

            A Python matplotlib pyplot subplot object.

        Keywords
        --------

        kwargs: dict, optional

            All keyword arguments to modify the default figure axes
            should be based via the dictionary upon entry.

        Raises
        ------

        InnovStatsPlotError:

            * raised if the axes attribute is NoneType.

        Notes
        -----

        A list of compliant keyword arguments may be found here:

        https://tinyurl.com/pyplot-axes

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
            if "axes" in kwargs:

                # Define the axes attributes contained within the key
                # word argument dictionary.
                axes_dict = parser_interface.dict_key_value(
                    dict_in=kwargs, key="axes", force=True
                )

                if axes_dict is None:
                    msg = (
                        "The axes attributes cannot be NoneType; check "
                        "the user experiment configuration is correctly "
                        "formatted. Aborting!!!"
                    )
                    error(msg=msg)

                # Build the axes according to the available
                # attributes.
                for key in axes_dict.keys():

                    # Define the available axes attributes.
                    value = parser_interface.dict_key_value(
                        dict_in=axes_dict, key=key, no_split=True
                    )
                    axes_obj = parser_interface.object_setattr(
                        object_in=axes_obj, key=key, value=value
                    )

                # Build the axes attributes accordingly.
                (plot_obj, ax_obj) = self.__xaxes__(
                    plot_obj=plot_obj,
                    ax_obj=ax_obj,
                    axes_obj=axes_obj,
                    axes_dict=axes_dict,
                )
                (plot_obj, ax_obj) = self.__yaxes__(
                    plot_obj=plot_obj,
                    ax_obj=ax_obj,
                    axes_obj=axes_obj,
                    axes_dict=axes_dict,
                )

            else:

                # Proceed without build the axis object.
                msg = (
                    "The axis attributes have not been specified; the "
                    "figure axes will not be modified."
                )

        if not kwargs:

            # Proceed without build the axis object.
            msg = (
                "The axis attributes have not been specified; the "
                "figure axes will not be modified."
            )
            self.logger.warn(msg=msg)

    def build_colorlevs(self, cmap: str, nlevs: int) -> list:
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

        InnovStatsPlotError:

            * raised if an exception is encountered while attempting
              to define the respective color levels.

        """

        # Define the color levels relative to the specified
        # attributes.
        try:
            data_arr = numpy.linspace(0, 1.0, nlevs)
            colorlevs = colorlist.gen_color_normalized(cmap=cmap, data_arr=data_arr)

        except Exception as errmsg:
            msg = (
                "The normalized color levels generation failed "
                f"with error {errmsg}. Aborting!!!"
            )
            error(msg=msg)

        return colorlevs

    def build_figure(self, suppress_clf: bool = True) -> Tuple[object, object, object]:
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

        plot_obj: object

            A Python matplotlib pyplot object.

        fig_obj: object

            A Python matplotlib pyplot figure-level object.

        ax_obj: object

            A Python matplotlib pyplot axis-level subplot object.

        """

        # Initialize the matplotlib pyplot object and define the
        # respective attributes and proceed accordingly.
        msg = "Initializing new pyplot objects."
        self.logger.warn(msg=msg)

        plot_obj = matplotlib.pyplot
        if suppress_clf:
            try:
                plot_obj.close("all")
                plot_obj.clf()

            except Exception:
                pass

        fig_obj = plot_obj.figure()
        ax_obj = plot_obj.subplot()

        return (plot_obj, fig_obj, ax_obj)

    def build_hlines(self, plot_obj: object, **kwargs: dict) -> None:
        """
        Description
        -----------

        This method parses the keyword arguments dictionary to
        determine the configuration for the horizontal lines (hlines)
        for the respective figure.

        Parameters
        ----------

        plot_obj: object

            A Python matplotlib pyplot object.

        Keywords
        --------

        kwargs: dict, optional

            All keyword arguments to modify the default hlines should
            be based via the dictionary upon entry.

        Notes
        -----

        A list of compliant keyword arguments may be found here:

        https://tinyurl.com/pyplot-hlines

        """

        # Parse the keyword argument dictionary; proceed accordingly.
        if kwargs is not None:
            hlines_kwargs = parser_interface.dict_key_value(
                dict_in=kwargs, key="hlines", force=True, no_split=True
            )

            if hlines_kwargs is None:
                msg = (
                    "No horizontal line attributes were found in the user "
                    "experiment configuration; horizontal lines will not be "
                    "plotted."
                )
                self.logger.warn(msg=msg)

            if hlines_kwargs is not None:
                plot_obj.hlines(**hlines_kwargs)

    def build_legend(self, plot_obj: object, **kwargs: dict) -> None:
        """
        Description
        -----------

        This method parses the keyword arguments dictionary to
        determine the configuration for the legend for the respective
        figure.

        Parameters
        ----------

        plot_obj: object

            A Python matplotlib pyplot object.

        Keywords
        --------

        kwargs: dict, optional

            All keyword arguments to modify the default figure legend
            should be based via the dictionary upon entry.

        Notes
        -----

        A list of compliant keyword arguments may be found here:

        https://tinyurl.com/pyplot-legend

        """

        # Parse the keyword argument dictionary; proceed accordingly.
        legend_kwargs = parser_interface.dict_key_value(
            dict_in=kwargs, key="legend", force=True, no_split=True
        )

        if legend_kwargs is None:
            msg = (
                "No legend attributes were found in the user experiment "
                "configuration; the legend will not be created."
            )
            self.logger.warn(msg=msg)

        if legend_kwargs is not None:
            plot_obj.legend(**legend_kwargs)

    def build_plotsinfo(self) -> None:
        """
        Description
        -----------

        This method builds the base-class attribute by collecting all
        plotting information for the respective variables and
        respective statistics specified upon entry by vars_list and
        stats_list, respectively.

        Raises
        ------

        InnovStatsPlotError:

            * raised if the variables attribute cannot be determined
              from the user experiment configuration.

        """

        # Collect the variable attributes from the user experiment
        # configuration.
        variables_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key="variables", force=True, no_split=True
        )
        if variables_dict is None:
            msg = (
                "The variables attribute could not be determined "
                "from the user experiment configuration. "
                "Aborting!!!"
            )
            error(msg=msg)

        # Define the attributes for all variables to be plotted; if
        # not specified reset to the defined default value(s).
        self.plots_obj = parser_interface.object_define()
        plots_attrs_dict = {"create_images": False}

        for plots_attr in plots_attrs_dict:
            value = parser_interface.dict_key_value(
                dict_in=plots_attrs_dict, key=plots_attr, force=True, no_split=True
            )
            if value is None:
                value = parser_interface.dict_key_value(
                    dict_in=plots_attrs_dict, key=plots_attr, no_split=True
                )
                msg = (
                    f"The plotting attribute {plots_attr} has not been "
                    f"specified and will be set to {value}."
                )
                self.logger.warn(msg=msg)

            # Loop through the specified variables and statistics
            # values and build the base-class attribute plots_obj.
            for var in self.vars_list:
                stat_obj = parser_interface.object_define()

                for stat in self.stats_list:
                    value = parser_interface.dict_key_value(
                        dict_in=variables_dict, key=stat, force=True, no_split=True
                    )
                    if value is None:
                        msg = (
                            f"The statistics attribute {stat} for variable "
                            f"{var} has not been specified within the user "
                            "experiment configuration; it will be reset "
                            "to NoneType."
                        )
                        self.logger.warn(msg=msg)
                        stat_obj = None

                    if value is not None:
                        stat_obj = parser_interface.object_setattr(
                            object_in=stat_obj, key=stat, value=value
                        )

                self.plots_obj = parser_interface.object_setattr(
                    object_in=self.plots_obj, key=var, value=stat_obj
                )

    def build_title(self, plot_obj: object, **kwargs: dict) -> None:
        """
        Description
        -----------

        This method parses the keyword arguments dictionary to
        determine the configuration for the title for the respective
        figure.

        Parameters
        ----------

        plot_obj: object

            A Python matplotlib pyplot object.

        Keywords
        --------

        kwargs: dict, optional

            All keyword arguments to modify the default figure title
            should be based via the dictionary upon entry.

        Notes
        -----

        A list of compliant keyword arguments may be found here:

        https://tinyurl.com/pyplot-title

        """

        # Parse the keyword argument dictionary; proceed accordingly.
        if kwargs is not None:
            try:
                plot_obj.title(**kwargs)

            except TypeError:
                try:
                    kwargs["s"] = kwargs.pop("label")

                except KeyError:
                    kwargs["label"] = kwargs.pop("s")

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

        InnovStatsPlotError:

            * raised if the user experiment configuration does not
              contain the 'variables' attribute.

            * raised if the attributes for the respective analysis
              variable cannot be determined from the user experiment
              configuration.

        """

        # Parse the user experiment configuration for the respective
        # analysis variable and proceed accordingly.
        variable_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key="variables", force=True, no_split=True
        )
        if variable_dict is None:
            msg = (
                "The variables attribute could not be determined from "
                "the user experiment configuration. Aborting!!!"
            )
            error(msg=msg)

        varinfo_dict = parser_interface.dict_key_value(
            dict_in=variable_dict, key=variable, force=True, no_split=True
        )
        if varinfo_dict is None:
            msg = (
                f"The attributes for variable {variable} could not be determined "
                "from the user experiment configuration. Aborting!!!"
            )
            error(msg=msg)

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

        https://tinyurl.com/pyplot-vlines

        """

        # Parse the keyword argument dictionary; proceed accordingly.
        if kwargs is not None:
            vlines_kwargs = parser_interface.dict_key_value(
                dict_in=kwargs, key="vlines", force=True, no_split=True
            )
            if vlines_kwargs is None:
                msg = (
                    "No vertical line attributes were found in the user "
                    "experiment configuration; vertical lines will not be "
                    "plotted."
                )
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

        https://tinyurl.com/pyplot-xlabel

        """

        # Parse the keyword argument dictionary; proceed accordingly.
        labels_kwargs = parser_interface.dict_key_value(
            dict_in=kwargs, key="xlabel", force=True, no_split=True
        )
        if labels_kwargs is None:
            msg = (
                "No x-label attributes were found in the user experiment "
                "experiment configuration; x-labels will not be created."
            )
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

        https://tinyurl.com/pyplot-ylabel

        """

        # Parse the keyword argument dictionary; proceed accordingly.
        labels_kwargs = parser_interface.dict_key_value(
            dict_in=kwargs, key="ylabel", force=True, no_split=True
        )
        if labels_kwargs is None:
            msg = (
                "No y-label attributes were found in the user experiment "
                "experiment configuration; y-labels will not be created."
            )
            self.logger.warn(msg=msg)

        if labels_kwargs is not None:
            ax_obj.set_ylabel(**labels_kwargs)

    def create_image(self, plot_obj, savename, **kwargs):
        """
        Description
        -----------

        This method creates an external file containing the respective
        figure.

        Parameters
        ----------

        plot_obj: object

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

        https://tinyurl.com/pyplot-savefig

        """

        # Parse the keyword argument dictionary; proceed accordingly.
        msg = f"Creating figure {savename}."
        self.logger.info(msg=msg)

        savefig_kwargs = parser_interface.dict_key_value(
            dict_in=kwargs, key="savefig", force=True, no_split=True
        )
        if savefig_kwargs is None:
            savefig_kwargs = {}

        plot_obj.savefig(fname=savename, **savefig_kwargs)
