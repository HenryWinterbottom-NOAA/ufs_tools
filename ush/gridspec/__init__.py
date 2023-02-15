# =========================================================================

# Module: ush/gridspec/__init__.py

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

    This module contains the base-class for all reduced grid
    definitions determined by a defined supergrid structure specified
    upon entry.

Classes
-------

    GridSpec(options_obj)

        This is the base-class for all reduced-grid definitions using
        a supergrid structure specified upon input.

Requirements
------------

- ufs_pytils; https://github.com/HenryWinterbottom-NOAA/ufs_pyutils

Author(s)
---------

    Henry R. Winterbottom; 08 February 2023

History
-------

    2023-02-08: Henry Winterbottom -- Initial implementation.

"""

# ----

# pylint: disable=too-many-instance-attributes

# ----

import os
from dataclasses import dataclass
from typing import Tuple, Union

import numpy
from confs.yaml_interface import YAML
from exceptions import GridSpecError
from ioapps import netcdf4_interface
from tools import datetime_interface, fileio_interface, parser_interface
from utils import timestamp_interface
from utils.logger_interface import Logger

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


@dataclass
class GridSpec:
    """
    Description
    -----------

    This is the base-class for all reduced-grid definitions using a
    supergrid structure specified upon input.

    Parameters
    ----------

    options_obj: object

        A Python object containing the command line argument
        attributes.

    """

    def __init__(self, options_obj: object):
        """
        Description
        -----------

        Creates a new GridSpec object.

        """

        # Define base-class attributes.
        self.options_obj = options_obj
        self.logger = Logger()
        self.yaml_file = self.options_obj.yaml_file

        self.yaml_dict = YAML().read_yaml(yaml_file=self.yaml_file)

        self.grid_attr_list = ["angle", "latitude", "longitude"]

        self.reduce_grid_obj = parser_interface.object_define()
        self.reduce_grid_list = [
            "angle",
            "qlat",
            "qlon",
            "tlat",
            "tlon",
            "ulat",
            "ulon",
            "vlat",
            "vlon",
        ]

        (self.grids_obj, self.ncdim_obj, self.ncvar_obj) = [
            parser_interface.object_define() for idx in range(3)
        ]

        # Define the format/orientation of the resulting grid
        # projection; proceed accordingly.
        self.is_tripolar = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key="is_tripolar", force=True
        )
        self.is_tripolar = self.is_tripolar if self.is_tripolar else False

        self.is_wrap_lons = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key="is_wrap_lons", force=True
        )
        self.is_wrap_lons = self.is_wrap_lons if self.is_wrap_lons else False

    def read_ncfile(self) -> None:
        """
        Description
        -----------

        This method parses the specified netCDF-formatted file path
        containing the supergrid structure; the supergrid structure is
        encapsulated within the base-class object grids_obj.

        Raises
        ------

        GridSpecError:

            * raised if the specified netCDF attribute can not be
              determined for the grid attribute defined within the
              experiment configuration file.

            * raised if a specified mandatory attribute cannot be
              determined from the experiment configuration file; see
              base-class attribute grid_attr_list.

            * raised if the specified netCDF-formatted file path does
              not exist upon entry.

        """

        # Define the netCDF attributes for the grid variables.
        nc_attr_list = ["ncfile", "ncvarname", "ncxdim", "ncydim"]

        for grid_attr in self.grid_attr_list:

            # Collect the YAML attributes for the respective grid
            # attribute; proceed accordingly.
            grid_dict = parser_interface.dict_key_value(
                dict_in=self.yaml_dict, key=grid_attr, force=True
            )
            if grid_dict is None:
                msg = (
                    f"The attributes for grid attribute {grid_attr} could not "
                    f"be determined from configurtion file {self.yaml_file}. "
                    "Aborting!!!"
                )
                raise GridSpecError(msg=msg)

            # Collect the netCDF attributes; proceed accordingly.
            nc_obj = parser_interface.object_define()

            for nc_attr in nc_attr_list:
                value = parser_interface.dict_key_value(
                    dict_in=grid_dict, key=nc_attr, force=True, no_split=True
                )
                if value is None:
                    msg = (
                        f"netCDF attribute {nc_attr} could not be determined for "
                        f"grid {grid_attr} in configurtion file {self.yaml_file}. "
                        "Aborting!!!"
                    )
                    raise GridSpecError(msg=msg)

                nc_obj = parser_interface.object_setattr(
                    object_in=nc_obj, key=nc_attr, value=value
                )

            # Check that the netCDF-formatted file path exists;
            # proceed accordingly.
            fileexist = fileio_interface.fileexist(path=nc_obj.ncfile)
            if not fileexist:
                msg = (
                    f"The netCDF-formatted file path {nc_obj.ncfile} does not exist. "
                    "Aborting!!!"
                )
                raise GridSpecError(msg=msg)

            # Collect the netCDF variable values.
            msg = (
                f"Collecting grid variable {grid_attr} attributes from netCDF-formatted "
                f"file path {nc_obj.ncfile}."
            )
            self.logger.info(msg=msg)
            ncvalues = netcdf4_interface.ncreadvar(
                ncfile=nc_obj.ncfile, ncvarname=nc_obj.ncvarname
            )
            nc_obj = parser_interface.object_setattr(
                object_in=nc_obj, key="ncvalues", value=ncvalues
            )

            self.grids_obj = parser_interface.object_setattr(
                object_in=self.grids_obj, key=grid_attr, value=nc_obj
            )

    def update_ncvar(self, ncvar_dict: dict, ncvar_obj: object) -> object:
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

        ncvar_obj: object

            A Python object containing the netCDF variable object;
            this object is the accumulated netCDF variable attributes
            to be written to the specified netCDF-formatted output
            file path.

        Returns
        -------

        ncvar_obj: object

            A Python object containing the updated netCDF variable
            object; this object is the accumulated netCDF variable
            attributes to be written to the specified netCDF-formatted
            output file path.

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

    def wrap_lons(
        self, lons: numpy.array, angle: numpy.array = None
    ) -> Tuple[numpy.array, Union[numpy.array, None]]:
        """
        Description
        -----------

        This method will order the longitude coordinate array,
        provided upon entry, within the range [0,
        numpy.degrees(2*numpy.pi)].

        Parameters
        ----------

        lons: numpy.array

            A Python numpy.array containing longitude values.

        Keywords
        --------

        angle: numpy.array

            ----

        Returns
        -------

        lons: numpy.array

            A Python numpy.array containing the ordered longitude
            values.

        """

        # Define the longitude scaling value.
        scale = 2.0 * numpy.degrees(numpy.pi)

        # Update the longitude values accordingly.
        lons = numpy.where((lons > scale), (scale - lons), lons)
        lons = numpy.where((lons < 0.0), (lons + scale), lons)

        if angle is not None:
            angle = -1.0 * angle

        return (lons, angle)

    def write_ncfile(self) -> None:
        """
        Description
        -----------

        This method writes the netCDF-formatted output file path using
        the attributes contained within the base-class attributes
        ncdim_obj and ncvar_obj; if a netCDF-formatted output file
        path has not been specified the reduced grid definitions to
        the file gridspec.nc in the run-time directory tree path.

        """

        # Define the netCDF-formatted output path; proceed
        # accordingly.
        ncfile = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key="output_netcdf", force=True, no_split=True
        )
        if ncfile is None:

            # Define a generic netCDF-formatted file output path.
            ncfile = os.path.join(os.getcwd(), "gridspec.nc")

            msg = (
                f"The experiment configuration file {self.yaml_file} does not "
                f"specify an output file name; setting to {ncfile}."
            )
            self.logger.warn(msg=msg)

        # Define the global attributes for the netCDF-formatted file.
        glbattrs_dict = {
            "_FillValue": numpy.nan,
            "date": datetime_interface.current_date(timestamp_interface.INFO),
            "created": str(parser_interface.enviro_get(envvar="HOSTNAME")),
        }

        # Write the netCDF-formatted file.
        msg = f"Writing reduced to path {ncfile}."
        self.logger.info(msg=msg)
        netcdf4_interface.ncwrite(
            ncfile=ncfile,
            ncdim_obj=self.ncdim_obj,
            ncvar_obj=self.ncvar_obj,
            ncfrmt="NETCDF4",
            glbattrs_dict=glbattrs_dict,
        )

    def run(self) -> None:
        """
        Description
        -----------

        This method is generic and used to run the methods for the
        respective calling class (i.e., GridSpec sub-class)
        applications.

        """
