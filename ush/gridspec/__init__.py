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



"""

# ----

import os
import numpy

from confs.yaml_interface import YAML

from exceptions import GridSpecError

from tools import fileio_interface
from tools import parser_interface
from ioapps import netcdf4_interface

from utils.logger_interface import Logger

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


class GridSpec:
    """

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

        self.grid_attr_list = ["latitude", "longitude", "mask", "topography"]

        self.reduce_grid_obj = parser_interface.object_define()
        self.reduce_grid_list = ["qlat", "qlon", "tlat", "tlon", "ulat", "ulon",
                                 "vlat", "vlon"
                                 ]

        (self.grids_obj, self.ncdim_obj, self.ncvar_obj) = \
            [parser_interface.object_define() for idx in range(3)]

    def read_ncfile(self):
        """ """

        # Define the netCDF attributes for the grid variables.
        nc_attr_list = ["ncfile", "ncvarname", "ncxdim", "ncydim"]

        for grid_attr in self.grid_attr_list:

            # Collect the YAML attributes for the respective grid
            # attribute; proceed accordingly.
            grid_dict = parser_interface.dict_key_value(
                dict_in=self.yaml_dict, key=grid_attr, force=True)
            if grid_dict is None:
                msg = (f"The attributes for grid attribute {grid_attr} could not "
                       f"be determined from configurtion file {self.yaml_file}. "
                       "Aborting!!!")
                raise GridSpecError(msg=msg)

            # Collect the netCDF attributes; proceed accordingly.
            nc_obj = parser_interface.object_define()

            for nc_attr in nc_attr_list:
                value = parser_interface.dict_key_value(
                    dict_in=grid_dict, key=nc_attr, force=True, no_split=True)
                if value is None:
                    msg = (f"netCDF attribute {nc_attr} could not be determined for "
                           f"grid {grid_attr} in configurtion file {self.yaml_file}. "
                           "Aborting!!!"
                           )
                    raise GridSpecError(msg=msg)

                nc_obj = parser_interface.object_setattr(
                    object_in=nc_obj, key=nc_attr, value=value)

            # Check that the netCDF-formatted file path exists;
            # proceed accordingly.
            fileexist = fileio_interface.fileexist(path=nc_obj.ncfile)
            if not fileexist:
                msg = (f"The netCDF-formatted file path {nc_obj.ncfile} does not exist. "
                       "Aborting!!!")
                raise GridSpecError(msg=msg)

            # Collect the netCDF variable values.
            msg = (f"Collecting grid variable {grid_attr} attributes from netCDF-formatted "
                   f"file path {nc_obj.ncfile}."
                   )
            self.logger.info(msg=msg)
            ncvalues = netcdf4_interface.ncreadvar(ncfile=nc_obj.ncfile,
                                                   ncvarname=nc_obj.ncvarname)
            nc_obj = parser_interface.object_setattr(
                object_in=nc_obj, key="ncvalues", value=ncvalues)

            self.grids_obj = parser_interface.object_setattr(
                object_in=self.grids_obj, key=grid_attr, value=nc_obj)

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

    def write_ncfile(self) -> None:
        """ """

        # Define the netCDF-formatted output path; proceed
        # accordingly.
        ncfile = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key="output_netcdf", force=True,
            no_split=True)
        if ncfile is None:

            # Define a generic netCDF-formatted file output path.
            ncfile = os.path.join(os.getcwd(), "gridspec.nc")

            msg = (f"The experiment configuration file {self.yaml_file} does not "
                   f"specify an output file name; setting to {ncfile}."
                   )
            self.logger.warn()

        # Write the netCDF-formatted file.
        msg = (f"Writing reduced to path {ncfile}.")
        self.logger.info(msg=msg)
        netcdf4_interface.ncwrite(ncfile=ncfile, ncdim_obj=self.ncdim_obj,
                                  ncvar_obj=self.ncvar_obj, ncfrmt="NETCDF4")

    def run(self) -> None:
        """

        """
