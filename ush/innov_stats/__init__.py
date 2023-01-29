# =========================================================================

# Module: ush/innov_stats/__init__.py

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

    This module contains the base-class module for all innovation
    statistics calculations and netCDF formatted file creations.

Classes
-------

    InnovStats(yaml_file, basedir, cyclestrfrmt=None,
               ncconcat_filename=None)

        This is the base-class object for all innovation statistics
        sub-classes.

Functions
---------

    error(msg)

        This function is the exception handler for the respective
        module.

Author(s)
---------

    Henry R. Winterbottom; 15 January 2023

History
-------

    2023-01-15: Henry Winterbottom -- Initial implementation.

"""

# ----

# pylint: disable=broad-except
# pylint: disable=too-many-instance-attributes
# pylint: disable=too-many-lines
# pylint: disable=too-many-locals
# pylint: disable=undefined-loop-variable
# pylint: disable=unnecessary-dict-index-lookup
# pylint: disable=unused-argument

# ----

import os
import statistics
import sys

import numpy
import tabulate
from confs.yaml_interface import YAML
from exceptions import InnovStatsError
from ioapps import netcdf4_interface, sqlite3_interface
from tools import datetime_interface, fileio_interface, parser_interface
from utils import timestamp_interface
from utils.error_interface import msg_except_handle
from utils.logger_interface import Logger

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


class InnovStats:
    """
    Description
    -----------

    This is the base-class object for all innovation statistics
    calculations and output file generation.

    Parameters
    ----------

    options_obj: object

        A Python object containing the command line argument
        attributes.

    basedir: str

        A Python string specifying the top-level of the working
        directory tree.

    """

    def __init__(self: object, options_obj: object, basedir=str):
        """
        Description
        -----------

        Creates a new InnovStats object.

        """

        # Define the base-class attributes.
        # self = cls
        self.logger = Logger()
        self.options_obj = options_obj
        self.basedir = basedir
        self.cycle = self.options_obj.cycle

        self.exptname = parser_interface.object_getattr(
            object_in=self.options_obj, key="expt_name", force=True
        )
        if self.exptname is None:
            self.exptname = "UNKNOWN"
        self.yaml_dict = YAML().read_yaml(yaml_file=self.options_obj.yaml_file)

        # Define the innovation statistic types.
        self.stats_type_list = ["bias", "count", "rmsd"]

        # Collect the relevant information from the experiment
        # configuration.
        self.diagsinfo_obj = parser_interface.object_define()
        self.regions_obj = self.get_regions()
        self.diagsvars = self.get_diagsinfo()
        self.levels_obj = parser_interface.object_define()

    def build_database(
        self, column_frmt: str, levels_list: list = None, column_scale: float = 1.0
    ) -> None:
        """
        Description
        -----------

        This method initializes the SQLite3 database for each
        diagnostics variable; if the database already exists and a
        corresponding database table exists it will not be
        overwritten.

        Parameters
        ----------

        column_frmt: str

            A Python string specifying the format for the table columm
            names; it is specific to the sub-class instance of the
            base-class.

        Keywords
        --------

        levels_list: list, optional

            A Python list of levels; if NoneType on entry a single
            level is assumed.

        column_scale: float, optional

            A Python float value specifying the scaling value for the
            table column name creation/construction.

        Raises
        ------

        InnovStatsError:

            * raised if experiment configuration does not specify the
              database_sql_path attribute for the respective
              diagnostics variable.

        """

        # Loop through each diagnostics variable and proceed
        # accordingly.
        for diagsinfo in vars(self.diagsinfo_obj):

            # Define the SQLite3 database attributes and proceed
            # accordingly.
            diagsinfo_dict = parser_interface.object_getattr(
                object_in=self.diagsinfo_obj, key=diagsinfo
            )

            sql_database = parser_interface.dict_key_value(
                dict_in=diagsinfo_dict,
                key="database_sql_path",
                force=True,
                no_split=True,
            )
            if sql_database is None:
                msg = (
                    "The SQLite3 database path, database_sql_path, "
                    f"for variable {diagsinfo} has not be specified in the "
                    "experiment configuration. Aborting!!!"
                )
                error(msg=msg)

            msg = (
                f"The SQLite3 database table path for variable {diagsinfo} is "
                f"{sql_database}."
            )
            self.logger.info(msg=msg)

            # Check that the directory tree for the SQLite3 database
            # file exists; proceed accordingly.
            path = os.path.dirname(sql_database)
            exist = fileio_interface.fileexist(path=path)
            if not exist:
                msg = (
                    f"The directory path {path} does not exist; an attempt "
                    "will be made to create it."
                )
                self.logger.warn(msg=msg)
                fileio_interface.makedirs(path=path)

            # Check the status of the SQLite3 database file and
            # proceed accordingly.
            path = sql_database
            exist = fileio_interface.fileexist(path=path)
            if exist:
                msg = (
                    f"The SQLite3 database path {path} exists and will be " "appended."
                )
                self.logger.info(msg=msg)

            if not exist:
                msg = (
                    f"The SQLite3 database path {path} does not exist and it "
                    "will be initialized."
                )
                self.logger.warn(msg=msg)

            # Build the SQLite3 database tables for the respective
            # forecast cycle, statistic, and region of interest;
            # proceed accordingly.
            for stats_type in self.stats_type_list:
                for region in vars(self.regions_obj).keys():
                    table_name = f"{stats_type}_{region}"
                    table_dict = {}
                    table_dict["CYCLE"] = "INT"
                    if levels_list is None:
                        table_dict[column_frmt] = "REAL"

                    if levels_list is not None:
                        for level in levels_list:
                            table_dict[column_frmt %
                                       int(level * column_scale)] = "REAL"

                    # Create the SQLite3 database table.
                    sqlite3_interface.create_table(
                        path=path, table_name=table_name, table_dict=table_dict
                    )

    def build_ncdims(self, dimsdict: dict, ncvar_obj: object) -> object:
        """
        Description
        -----------

        This method defines the dimension attributes for the netCDF
        formatted file to be created.

        Parameters
        ----------

        dimsdict: dict

            A Python dictionary containing the netCDF variable key and
            value pairs; the Python dictionary keys are the respective
            netCDF coordinate dimension names while the respective
            Python key values are the netCDF coordinate dimension
            values.

        ncvar_obj: object

            A Python object containing the variable attributes for the
            netCDF-formatted file.

        Returns
        -------

        ncvar_obj: object

            A Python object containing the updated variable attributes
            for the netCDF-formatted file.

        Raises
        ------

        InnovStatsError:

            * raised if an error is encountered while defining the
              respective netCDF formatted file coordinate dimension
              attribute.

        """

        # Define the netCDF dimension attributes for the output
        # netCDF-formatted file.
        try:
            for key in dimsdict.keys():
                ncvar_dict = {
                    key: {
                        "varname": key,
                        "dims": key,
                        "type": dimsdict[key]["type"],
                        "values": dimsdict[key]["values"],
                    }
                }
                ncvar_obj = self.update_ncvar(
                    ncvar_dict=ncvar_dict, ncvar_obj=ncvar_obj
                )

        except Exception as errmsg:
            msg = (
                "Building the netCDF variable dimensions failed with "
                f"error {errmsg}. Aborting!!!"
            )
            error(msg=msg)

        return ncvar_obj

    def get_depths(self):
        """
        Description
        -----------

        This method collects the layer information from the user
        experiment configuration; the base-class attributes nlevs and
        levels_obj is defined and contains the following:

        layer_bottom: the bottom interface of a given layer; this is
                      collected from the user experiment
                      configuration.

        layer_top: the top interface of a given layer; this is
                   collected from the user experiment configuration.

        layer_mean: this is the mean value of the layer_bottom and
                    layer_top value for a given layer; it is computed
                    within this method.

        Raises
        ------

        InnovStatsError:

            * raised if the depth levels information cannot be
              determined from the user experiment configuration.

        """

        # Define the depth profile information from the user
        # experiment configuration and proceed accordingly.
        depth_levels_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key='depth_levels',
            force=True)

        if depth_levels_dict is None:
            msg = ('The depth levels information could not be '
                   'determined from the user experiment configuration. '
                   'Aborting!!!')
            error(msg=msg)

        # Define the base-class attributes levels_obj; this namespace
        # contains both the layer bounding regions (i.e., top and
        # bottom layer interfaces for the respective level) as well as
        # the mean value for the respective interval (i.e., the layer
        # mean -- the middle of the layer).
        levels_attr_dict = {'depth_bottom': 'layer_bottom', 'depth_top':
                            'layer_top'}
        for (levels_attr, _) in levels_attr_dict.items():
            value = parser_interface.dict_key_value(
                dict_in=depth_levels_dict, key=levels_attr,
                force=True)

            if value is None:
                msg = ('The attribute {0} could not be determined from '
                       'the user experiment configuration. Aborting!!!'.
                       format(levels_attr))
                error(msg=msg)

            levels = [float(level) for level in value]

            self.levels_obj = parser_interface.object_setattr(
                object_in=self.levels_obj, key=levels_attr_dict[levels_attr],
                value=levels)

        value = [statistics.mean(k) for k in zip(self.levels_obj.layer_bottom,
                                                 self.levels_obj.layer_top)]
        self.levels_obj = parser_interface.object_setattr(
            object_in=self.levels_obj, key='layer_mean', value=value)

    def get_diagsinfo(self) -> None:
        """
        Description
        -----------

        This method defines the base-class attributes diagsinfo_obj
        and diagsvars which contain the respective innovation
        statistics variable diagnostics information and a list of
        diagnostics variable names respectively.

        Raises
        ------

        InnovStatsError:

            * raised if the attribute diagsinfo cannot be determined
              from the experiment configuration file.

        """

        # Collect and define the diagnostics attributes from the user
        # experiment configuration file.
        diagsinfo = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key="diagsinfo", force=True
        )
        if diagsinfo is None:
            msg = (
                "The observation diagnostics attribute diagsinfo "
                "could not be determined from the experiment "
                "configuration. Aborting!!!"
            )
            error(msg=msg)

        diagsvars = list(diagsinfo.keys())

        for diagsvar in diagsvars:
            value = parser_interface.dict_key_value(
                dict_in=diagsinfo, key=diagsvar, no_split=True
            )
            self.diagsinfo_obj = parser_interface.object_setattr(
                object_in=self.diagsinfo_obj, key=diagsvar, value=value
            )

        return diagsvars

    def get_innovinfo(self, vardict: dict, addinfo_list: list = None) -> object:
        """
        Description
        -----------

        This method returns the innovation statistics information
        using the respective variable information collected from the
        experiment configuration (vardict); this class, by default,
        defines the minimum innovation attribute (i.e., nclat and
        nclon); additional innovation statistics information to be
        collected from the experiment configuration may be passed via
        the keyword argument addinfo_list.

        Parameters
        ----------

        vardict: dict

            A Python dictionary containing the variable information
            collected from the experiment configuration file.

        Keywords
        --------

        addinfo_list: list, optional

            A Python list containing the additional innovatrion
            statistics information be collected from the experiment
            configuration file.

        Returns
        -------

        innovinfo_obj: object

            A Python object containing the innovation statistics
            information for the requested variable.

        Raises
        ------

        InnovStatsError:

            * raised if an error is encountered while appending the
              items in the keyword argument addinfo_list.

            * raised if a specified netCDF attribute could not be
              determined for the requested variable.

        """

        # Initialize the innovation statistics object.
        innovinfo_obj = parser_interface.object_define()
        innovinfo_attr_list = ["nclat", "nclon"]

        if addinfo_list is not None:
            try:
                for item in addinfo_list:
                    innovinfo_attr_list.append(item)

            except Exception as errmsg:
                msg = (
                    "The appending of the addinfo_list items failed with "
                    f"error {errmsg}. Aborting!!!"
                )
                error(msg=msg)

        # Collect the innovation statistics information from the user
        # experiment configuration.
        for innovinfo_attr in innovinfo_attr_list:
            ncinfo_dict = self.get_ncinfo(vardict=vardict)
            value = parser_interface.dict_key_value(
                dict_in=ncinfo_dict, key=innovinfo_attr, force=True, no_split=True
            )
            if value is None:
                msg = (
                    f"The netCDF variable attribute {innovinfo_attr} could not be "
                    "deterimined from the experiment configuration. "
                    "Aborting!!!"
                )
                error(msg=msg)

            innovinfo_obj = parser_interface.object_setattr(
                object_in=innovinfo_obj, key=innovinfo_attr, value=value
            )

        return innovinfo_obj

    def get_ncfilename(self, vardict: dict, ncdim: str = None) -> str:
        """
        Description
        -----------

        This method defines the netCDF filename path for the
        respective variable; if multiple netCDF file paths for a given
        variable are provided in the experiment configuration this
        method will concatenate the respective netCDF files into a
        single netCDF formatted file (specified by the base-class
        attribute ncconcat_filename) and along the netCDF coordinate
        dimension specified by ncdim on entry.

        Parameters
        ----------

        vardict: dict

            A Python dictionary containing the respective variable
            attributes collected from the experiment configuration.

        Keywords
        --------

        ncdim: str, optional

            A Python string specifying the netCDF coordinate dimension
            variable along which to concatenate multiple netCDF
            formatted files.

        Returns
        -------

        ncfilename: str

            A Python string specifying the netCDF formatted filename
            path for the respective variable.

        Raises
        ------

        InnovStatsError:

            * raised if the netCDF formatted file could not be
              determined from the experiment configuration.

            * raised if the netCDF coordinate dimension variable name
              along which to concatenate multople netCDF formatted
              file paths is not specified upon entry; this may only be
              raised if multiple input netCDF formatted file names are
              specified in the experiment configuration.

            * raised if base-class constructor variable
              ncconcat_filename has not be specified upon base-class
              object instantation; this may only be raised if multiple
              input netCDF formatted file names are specified in the
              experiment configuration.

        """

        # Collect the list of netCDF files for the respective
        # diagnostics variable as specified within the experiment
        # configuration.
        ncfilelist = parser_interface.dict_key_value(
            dict_in=vardict, key="ncinfile", force=True, no_split=True
        ).split(",")
        if ncfilelist is None:
            msg = (
                "The netCDF formatted file could not be determined "
                "from the experiment configuration. Aborting!!!"
            )
            error(msg=msg)

        # Build the netCDF file list for the respective diagnostics
        # variable.
        offset_seconds = parser_interface.dict_key_value(
            dict_in=vardict, key="offset_seconds", force=True
        )
        if offset_seconds is None:
            offset_seconds = 0

        ncfilenames = []
        for ncfilename in ncfilelist:
            ncfile = datetime_interface.datestrupdate(
                datestr=self.cycle,
                in_frmttyp=timestamp_interface.GLOBAL,
                out_frmttyp=ncfilename.strip(),
                offset_seconds=offset_seconds,
            )

            path = os.path.join(self.basedir, ncfile)
            exist = fileio_interface.fileexist(path=path)
            if exist:
                msg = f"File {path} exists and will be processed."
                self.logger.info(msg=msg)
                ncfilenames.append(path)

            if not exist:
                msg = f"File {path} does not exist and will not be processed."
                self.logger.warn(msg=msg)

        # Check the number of netCDF files available for the
        # respective diagnostics variable; proceed accordingly.
        if len(ncfilenames) > 1:
            if ncdim is None:
                msg = (
                    "In order to concatenate netCDF files the netCDF "
                    "dimension axis which to concatenate along must be "
                    "specified; got NoneType on entry. Aborting!!!"
                )
                error(msg=msg)

            ncfilename = os.path.join(self.basedir, "ncconcat.nc")
            netcdf4_interface.ncconcat(
                ncfilelist=ncfilenames, ncfile=ncfilename, ncdim=ncdim
            )

        else:
            ncfilename = datetime_interface.datestrupdate(
                datestr=self.cycle,
                in_frmttyp=timestamp_interface.GLOBAL,
                out_frmttyp=ncfilename,
                offset_seconds=offset_seconds,
            ).strip()
            path = os.path.join(self.basedir, ncfilename)
            exist = fileio_interface.fileexist(path=path)

            if not exist:
                msg = (
                    f"The netCDF file path {path} does not exist; the file "
                    "will not be processed."
                )
                self.logger.warn(msg=msg)
                ncfilename = None

        return ncfilename

    def get_ncinfo(self, vardict: dict) -> dict:
        """
        Description
        -----------

        This method collects the netCDF file variable attributes for
        the respective variable in the experiment configuration.

        Parameters
        ----------

        vardict: dict

            A Python dictionary containing the respective variable
            attributes collected from the experiment configuration.

        Returns
        -------

        ncinfo_dict: dict

            A Python dictionary containing the netCDF attributes for
            the respective variable.

        Raises
        ------

        InnovStatsError:

            * raised if the ncinfo attribute cannot be determined from
              the experiment configuration.

        """

        # Collect the netCDF attributes for the respective diagnostics
        # variable.
        ncinfo_dict = parser_interface.dict_key_value(
            dict_in=vardict, key="ncinfo", force=True
        )

        if ncinfo_dict is None:
            msg = (
                "The netCDF information could not be determined from "
                "the provided dictionary. Aborting!!!"
            )
            error(msg=msg)

        return ncinfo_dict

    def get_obslocs(
        self, ncfilename: str, innovinfo_obj: object, region: str
    ) -> numpy.array:
        """
        Description
        -----------

        This method defines a boolean array specifying the locations
        of observations valid for the respective (specified) region.

        Parameters
        ----------

        ncfilename: str

            A Python string specifying the netCDF filename path.

        innovinfo_obj: object

            A Python object containing the innovation statistics
            information for the requested variable.

        region: str

            A Python string specifying the region for which to
            determine the valid observation locations; this value must
            correspond to a value within the experiment configuration.

        Returns
        -------

        obslocs: array-type

            A Python array type variable containing the boolean values
            specifying the locations of observations valid for the
            respective (specified) region.

        Raises
        ------

        InnovStatsError:

            * raised if an error is encountered while determining the
              locations of valid observations.

        """

        # Read the netCDF-formatted latitude and longitude arrays from
        # the respective netCDF-formatted file and bin the latitude
        # and longitude coordinate array values relative to the region
        # of interest.
        try:

            # Define the regions of interest in accordance with the
            # experiment configuration file; proceed accordingly.
            region_attrs_dict = parser_interface.object_getattr(
                object_in=self.regions_obj, key=region
            )

            regions_obj = parser_interface.object_define()
            for region_attr in region_attrs_dict.keys():
                value = parser_interface.dict_key_value(
                    dict_in=region_attrs_dict, key=region_attr, no_split=True
                )
                regions_obj = parser_interface.object_setattr(
                    object_in=regions_obj, key=region_attr, value=value
                )

            # Collect the geographical coordinatte values and update
            # accordingly.
            nclats = netcdf4_interface.ncreadvar(
                ncfile=ncfilename, ncvarname=innovinfo_obj.nclat
            )

            nclons = netcdf4_interface.ncreadvar(
                ncfile=ncfilename, ncvarname=innovinfo_obj.nclon
            )
            if nclons.min() < 0.0:
                nclons = nclons + 180.0

            # Define the valid observation locations.
            obslocs = numpy.logical_and(
                numpy.logical_and(
                    nclats >= regions_obj.lat_min, nclats <= regions_obj.lat_max
                ),
                nclons >= regions_obj.lon_min,
                nclons <= regions_obj.lon_max,
            )

        except Exception as errmsg:
            msg = (
                "The determination of valid observation locations for "
                f"region {region} failed with error {errmsg}. Aborting!!!"
            )
            error(msg=msg)

        return obslocs

    def get_pressures(self):
        """
        Description
        -----------

        This method collects the layer information from the user
        experiment configuration; the base-class attributes nlevs and
        levels_obj is defined and contains the following:

        layer_bottom: the bottom interface of a given layer; this is
                      collected from the user experiment
                      configuration.

        layer_top: the top interface of a given layer; this is
                   collected from the user experiment configuration.

        layer_mean: this is the mean value of the layer_bottom and
                    layer_top value for a given layer; it is computed
                    within this method.

        Raises
        ------

        InnovStatsError:

            * raised if the pressure levels information cannot be
              determined from the user experiment configuration.

        Notes
        -----

        * The orientation of the layers for the atmosphere column is
          flipped (i.e., the layer_top value is the bottom of the
          layer and layer_bottom is the top of the layer) due to the
          orientation of the atmosphere pressure (isobaric) profile.

        """

        # Define the pressure-level information from the user
        # experiment configuration and proceed accordingly.
        pressure_levels_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key="pressure_levels", force=True
        )

        if pressure_levels_dict is None:
            msg = (
                "The pressure levels information could not be "
                "determined from the user experiment configuration. "
                "Aborting!!!"
            )
            error(msg=msg)

        # Define the base-class attributes levels_obj; this namespace
        # contains both the layer bounding regions (i.e., top and
        # bottom layer interfaces for the respective level) as well as
        # the mean value for the respective interval (i.e., the layer
        # mean -- the middle of the layer).
        levels_attr_dict = {"plevs_bottom": "layer_top",
                            "plevs_top": "layer_bottom"}

        for (levels_attr, _) in levels_attr_dict.items():

            value = parser_interface.dict_key_value(
                dict_in=pressure_levels_dict, key=levels_attr, force=True
            )
            if value is None:
                msg = (
                    f"The attribute {levels_attr} could not be determined from "
                    "the user experiment configuration. Aborting!!!"
                )
                error(msg=msg)

            levels = [float(level) for level in value]

            self.levels_obj = parser_interface.object_setattr(
                object_in=self.levels_obj,
                key=levels_attr_dict[levels_attr],
                value=levels,
            )

        value = [
            statistics.mean(k)
            for k in zip(self.levels_obj.layer_bottom, self.levels_obj.layer_top)
        ]

        self.levels_obj = parser_interface.object_setattr(
            object_in=self.levels_obj, key="layer_mean", value=value
        )

    def get_regions(self) -> object:
        """
        Description
        -----------

        This method collects the regions of interest from the user
        experiment configuration and defines the base-class attribute
        regions_obj and returns a table to the user stdout containing
        the attributes for the respective regions.

        Returns
        -------

        regions_obj: object

            A Python object containing the region attributes from the
            experiment configuration.

        Raises
        ------

        InnovStatsError:

            * raised if the regions attribute cannot be determined
              from the experiment configuration.

        """

        # Initialize the base-class object and define the default
        # values for the bounding geographical coordinate values for
        # the respective regions of interest specified within the user
        # experiment configuration.
        regions_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key="regions", force=True, no_split=True
        )
        if regions_dict is None:
            msg = (
                "The regions attribute could not be determined from "
                "the experiment configuration. Aborting!!!"
            )
            error(msg=msg)

        region_default_dict = {
            "lat_min": -90.0,
            "lat_max": 90.0,
            "lon_min": 0.0,
            "lon_max": 360.0,
        }
        regions_obj = parser_interface.object_define()

        # Loop through each region of interest and proceed
        # accordingly; compile the respective row of values for the
        # respective region.
        table_list = []

        for region in regions_dict.keys():
            region_obj = parser_interface.object_define()
            region_obj = parser_interface.object_setattr(
                object_in=region_obj, key="region", value=region
            )

            for (region_info, _) in region_default_dict.items():
                value = parser_interface.dict_key_value(
                    dict_in=regions_dict[region], key=region_info, force=True
                )
                if value is None:
                    value = parser_interface.dict_key_value(
                        dict_in=region_default_dict, key=region_info
                    )
                region_obj = parser_interface.object_setattr(
                    object_in=region_obj, key=region_info, value=value
                )

            row = [
                region_obj.region,
                region_obj.lat_min,
                region_obj.lat_max,
                region_obj.lon_min,
                region_obj.lon_max,
            ]
            table_list.append(row)

            region_dict = {
                "lat_min": region_obj.lat_min,
                "lat_max": region_obj.lat_max,
                "lon_min": region_obj.lon_min,
                "lon_max": region_obj.lon_max,
            }

            regions_obj = parser_interface.object_setattr(
                object_in=regions_obj, key=region, value=region_dict
            )

        # Print table for the experiment configuration regions of
        # interest and the respective region latitude and longitude
        # threshold values.
        headers = [
            "Region Name",
            "Minimum Latitude",
            "Maximum Latitude",
            "Minimum Longitude",
            "Maximum Longitude",
        ]

        kwargs = {"tablefmt": "fancy_grid", "numalign": "center", "stralign": "center"}

        table_obj = tabulate.tabulate(table_list, headers, **kwargs)
        sys.stdout.write(table_obj + "\n")

        return regions_obj

    def update_ncvar(self, ncvar_dict: dict, ncvar_obj: object) -> None:
        """
        Description
        -----------

        This method updates the base-class netCDF variable object
        (ncvar_obj).

        Parameters
        ----------

        ncvar_dict: dict

            A Python dictionary containing the respective netCDF
            variable attributes.

        """

        # Loop through all netCDF variable attributes and update the
        # base-class netCDF variable object.
        for key in ncvar_dict.keys():
            dict_in = {}
            for item in ncvar_dict[key].keys():
                value = parser_interface.dict_key_value(
                    dict_in=ncvar_dict[key], key=item, no_split=True
                )
                dict_in[item] = value

            # Update the netCDF variable object.
            ncvar_obj = parser_interface.object_setattr(
                object_in=ncvar_obj, key=key, value=dict_in
            )

        return ncvar_obj

    def write_database(
        self, vardict: dict, variable: str, column_frmt: str, column_scale: float
    ) -> None:
        """
        Description
        -----------

        This method updates a SQLite3 database file for the respective
        analysis date and variable as a function of region and
        statistic; if the analysis date already exists within the
        respective SQLite3 database file, it will be removed and
        updated with the most recently computed values.

        Parameters
        ----------

        vardict: dict

            A Python dictionary containing the respective variable
            attributes and diagnostics information.

        variable: str

            A Python string specifying the variable name.

        column_frmt: str

            A Python string specifying the column title string format.

        column_scale: float

            A Python string specifying the scaling value to be applied
            to the definition of the respective column titles.

        """

        # Initialize the SQLite3 database file attributes.
        offset_seconds = parser_interface.dict_key_value(
            dict_in=vardict, key="offset_seconds", force=True
        )
        if offset_seconds is None:
            offset_seconds = 0
        analdate = int(
            datetime_interface.datestrupdate(
                datestr=self.cycle,
                in_frmttyp=timestamp_interface.GLOBAL,
                out_frmttyp=timestamp_interface.GLOBAL,
                offset_seconds=offset_seconds,
            )
        )

        # Define the local attributes for the SQLite3 database file
        # creation.
        diagsinfo_dict = parser_interface.object_getattr(
            object_in=self.diagsinfo_obj, key=variable
        )
        database_sql_path = parser_interface.dict_key_value(
            dict_in=diagsinfo_dict, key="database_sql_path", no_split=True
        )
        ncoutfile = parser_interface.dict_key_value(
            dict_in=diagsinfo_dict, key="ncoutfile", no_split=True
        )
        ncoutfile = datetime_interface.datestrupdate(
            datestr=str(analdate),
            in_frmttyp=timestamp_interface.GLOBAL,
            out_frmttyp=ncoutfile,
        )

        # Loop through all regions specified within the user
        # experiment configuration and proceed accordingly.
        msg = f"Writing to SQLite3 database file {database_sql_path}."
        self.logger.info(msg=msg)
        for region in vars(self.regions_obj):

            for stats_type in self.stats_type_list:
                row_dict = {}
                row_dict["CYCLE"] = analdate
                ncvarname = f"{stats_type}_{region}"
                ncvalues = netcdf4_interface.ncreadvar(
                    ncfile=ncoutfile, ncvarname=ncvarname
                )
                table_name = f"{stats_type}_{region}"
                for (idx, lev) in enumerate(self.levels_obj.layer_mean):
                    if not numpy.isnan(ncvalues[idx]):
                        row_dict[column_frmt % int(lev * column_scale)] = numpy.float(
                            ncvalues[idx]
                        )
                contents = sqlite3_interface.read_table(
                    path=database_sql_path, table_name=table_name
                )

                # If the the respective analysis date previously
                # exists within the SQLite3 database, remove the row
                # (in order to prevent duplicates) and proceed
                # accordingly.
                for (row, _) in contents.items():
                    table_row = parser_interface.dict_key_value(
                        dict_in=contents, key=row, no_split=True
                    )
                    table_analdate = int(table_row[0])
                    if table_analdate == analdate:
                        msg = (
                            f"Found analysis date {table_name} in table {analdate}; "
                            "removing row and updating with new values."
                        )
                        self.logger.warn(msg=msg)

                        # Define the table removal criteria; proceed
                        # accordingly.
                        rmcond = f"CYCLE={analdate}"
                        sqlite3_interface.delete_row(
                            path=database_sql_path, table_name=table_name, rmcond=rmcond
                        )

                # Write/update the SQLite3 database table.
                sqlite3_interface.write_table(
                    path=database_sql_path, table_name=table_name, row_dict=row_dict
                )

    def write_ncout(
        self, vardict: dict, variable: str, ncdim_obj: object, ncvar_obj: object
    ) -> None:
        """
        Description
        -----------

        This method writes a netCDF formatted file containing the
        regions and variables of interest defined within the user
        experiment configuration.

        Parameters
        ----------

        vardict: dict

            A Python dictionary containing the respective variable
            attributes and diagnostics information.

        variable: str

            A Python string specifying the variable name.

        ncdim_obj: object

            A Python object containing the netCDF formatted file
            dimension attributes.

        ncvar_obj: object

            A Python object containing the netCDF formatted file
            variable attributes.

        Raises
        ------

        InnovStatsError:

            * raised if an error is encountered while creating the
              output netCDF formatted file for the respective
              variable.

        """

        # Define the netCDF file to contain the innovation statistics
        # for the respective specified statistic(s) (i.e., bias, rmse,
        # etc.,) and proceed accordingly.
        try:
            ncoutfile = parser_interface.dict_key_value(
                dict_in=vardict, key="ncoutfile", force=True, no_split=True
            )
            if ncoutfile is None:
                msg = (
                    "The output netCDF formatted file attribute (ncoutfile) "
                    f"for variable {variable} could not be determined from the user "
                    "experiment configuration. Aborting!!!"
                )
                error(msg=msg)

            offset_seconds = parser_interface.dict_key_value(
                dict_in=vardict, key="offset_seconds", force=True
            )
            if offset_seconds is None:
                offset_seconds = 0

            nc_cycle_attr_frmt = timestamp_interface.GENERAL
            expt_cycle = datetime_interface.datestrupdate(
                datestr=self.cycle,
                in_frmttyp=timestamp_interface.GLOBAL,
                out_frmttyp=nc_cycle_attr_frmt,
            )

            anly_cycle = datetime_interface.datestrupdate(
                datestr=self.cycle,
                in_frmttyp=timestamp_interface.GLOBAL,
                out_frmttyp=nc_cycle_attr_frmt,
                offset_seconds=offset_seconds,
            )

            glbattrs_dict = {
                "analysis_time": anly_cycle,
                "experiment_cycle": expt_cycle,
                "_FillValue": numpy.nan,
                "experiment_name": self.exptname,
                "variable": variable,
            }

            # Write the netCDF-formatted file containing the
            # innovation statistics.
            ncoutfile = datetime_interface.datestrupdate(
                datestr=self.cycle,
                in_frmttyp=timestamp_interface.GLOBAL,
                out_frmttyp=ncoutfile,
                offset_seconds=offset_seconds,
            )

            ncfile = os.path.join(self.basedir, ncoutfile)
            msg = (
                f"Creating innovation statistics file {ncfile} for variable {variable}."
            )
            self.logger.info(msg=msg)

            # Update the netCDF file with the respective variable.
            netcdf4_interface.ncwrite(
                ncfile=ncfile,
                ncdim_obj=ncdim_obj,
                ncvar_obj=ncvar_obj,
                ncfrmt="NETCDF3_64BIT_DATA",
                glbattrs_dict=glbattrs_dict,
            )

        except Exception as errmsg:
            msg = (
                f"The writing of the netCDF formatted file for variable {variable} "
                f"failed with error {errmsg}."
            )
            error(msg=msg)


# ----


@msg_except_handle(InnovStatsError)
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
