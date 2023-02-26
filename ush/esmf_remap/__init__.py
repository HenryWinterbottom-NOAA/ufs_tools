# =========================================================================

# Module: ush/esmf_remap/__init__.py

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

    This module contains the base-class object for all Earth System
    Modeling Framework (ESMF) supported source to destination grid
    interpolation type attributes.

Classes
-------

    ESMFRemap(options_obj)

        This is the base-class object for all Earth System Modeling
        Framework (ESMF) remapping attributes and remapping
        coefficient file generation.

Requirements
------------

- ufs_pytils; https://github.com/HenryWinterbottom-NOAA/ufs_pyutils

- xesmf; https://github.com/JiaweiZhuang/xESMF

Author(s)
---------

    Henry R. Winterbottom; 09 February 2023

History
-------

    2023-02-09: Henry Winterbottom -- Initial implementation.

"""

# ----

# pylint: disable=too-many-instance-attributes

# ----

import os
from dataclasses import dataclass

import xesmf
from confs.yaml_interface import YAML
from exceptions import ESMFRemapError
from tools import fileio_interface, parser_interface
from utils.logger_interface import Logger

from esmf_remap import ncio

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


@dataclass
class ESMFRemap:
    """
    Descriptions
    ------------

    This is the base-class object for all Earth System Modeling
    Framework (ESMF) remapping attributes and remapping coefficient
    file generation.

    Parameters
    ----------

    options_obj: object

        A Python object containing the attributes collect via the
        command line from the application driver script.

    Raises
    ------

    ESMFRemapError:

        * raised if the YAML-formatted configuration file does not
          exist upon entry.

    """

    def __init__(self, options_obj: object):
        """
        Description
        -----------

        Creates a new ESMFRemap object.

        """

        # Define the base-class attributes.
        self.logger = Logger()
        self.options_obj = options_obj
        self.yaml_file = self.options_obj.yaml_file

        # Define the experiment configuration; proceed accordingly.
        exist = fileio_interface.fileexist(path=self.yaml_file)
        if not exist:
            msg = (
                f"The YAML-formatted file {self.yaml_file} does not exist. "
                "Aborting!!!"
            )
            raise ESMFRemapError(msg=msg)

        self.yaml_dict = YAML().read_yaml(yaml_file=self.yaml_file)

        # Define the configuration attributes for the respective grids
        # and remapping application.
        self.grid_info_attr_list = [
            "is_grib",
            "gribfile",
            "is_netcdf",
            "ncfile",
            "ncxdim_name",
            "ncydim_name",
            "nclat",
            "nclon",
        ]

        self.remap_attrs_list = ["interp_type", "output_netCDF"]

        # Define the respective source and destination grid
        # attributes.
        self.srcgrid_info_obj = self.grid_info(grid_name="source")
        self.dstgrid_info_obj = self.grid_info(grid_name="destination")

        # Define the destination grid.
        if self.dstgrid_info_obj.is_grib:
            pass  # for now

        if self.dstgrid_info_obj.is_netcdf:
            self.dstgrid_obj = ncio.input_grids(grid_obj=self.dstgrid_info_obj)

        # Define the source grid.
        if self.srcgrid_info_obj.is_grib:
            pass  # for now

        if self.srcgrid_info_obj.is_netcdf:
            self.srcgrid_obj = ncio.input_grids(grid_obj=self.srcgrid_info_obj)

    def compute_remap(self) -> None:
        """
        Description
        -----------

        This method computes the ESMF remapping attributes and writes
        and external netCDF formatted file containing the remapping
        attributes for the interpolation type specified in the
        experiment configuration.

        Raises
        ------

        ESMFRemapError:

            * raised if a mandatory remapping attribute cannot be
              determined from the YAML-formatted experiment
              configuration file upon entry.

        """

        # Collecting the remapping attributes from the experiment
        # configuration.
        remap_obj = parser_interface.object_define()

        for remap_attr in self.remap_attrs_list:
            value = parser_interface.dict_key_value(
                dict_in=self.yaml_dict, key=remap_attr, force=True, no_split=True
            )
            if value is None:
                msg = (
                    f"The remapping attribute {remap_attr} could not be "
                    f"determined from the configuration file {self.yaml_file}. "
                    "Aborting!!!"
                )
                raise ESMFRemapError(msg=msg)

            remap_obj = parser_interface.object_setattr(
                object_in=remap_obj, key=remap_attr, value=value
            )

        # Define the empty arrays in accordance with the grid filename
        # path format.
        if self.srcgrid_info_obj.is_grib:
            pass  # for now

        if self.srcgrid_info_obj.is_netcdf:
            vararray = ncio.define_vararray(
                srcgrid_info_obj=self.srcgrid_info_obj)

        # Compute the ESMF regridding attributes.
        remapper = xesmf.Regridder(
            ds_in=self.srcgrid_obj,
            ds_out=self.dstgrid_obj,
            method=remap_obj.interp_type,
            periodic=True,
            reuse_weights=False,
            filename=remap_obj.output_netCDF,
        )
        remapper(vararray)

    def grid_info(self, grid_name: str) -> object:
        """
        Description
        -----------

        This method parses the YAML-formatted experiment configuration
        file path and collects the respective grid-type (i.e., source
        or destination) grid configuration attributes and builds a
        corresponding Python object; the returned object attributes
        specify either WMO GRIB formatted file path or the netCDF
        variable names for the respective geographical coordinate
        variables (e.g., latitude and longitude).

        Parameters
        ----------

        grid_name: str

            A Python string specifying the grid (e.g., source or
            destination) to be collected from the user experiment
            configuration.

        Returns
        -------

        grid_obj: object

            A Python object containing the respective grid (e.g.,
            source or destination) attributes.

        Raises
        ------

        ESMFRemapError:

            * raised if a grid name within the experiment
              configuration file is not supported; should be either
              source or destination for the respective YAML keys.

            * raised if the experiment configuration files does not
              contain the YAML key grids.

            * raised if the attributes for a grid type cannot be
              determined from the experiment configuration upon entry.

        """

        # Define the grid attributes; proceed accordingly.
        matching = [grid for grid in [
            "destination", "source"] if grid_name in grid]
        if not matching:
            msg = (
                f"The specified grid type {grid_name} is not supported. " "Aborting!!!"
            )
            raise ESMFRemapError(msg=msg)

        msg = (
            f"Collecting {grid_name} grid information from experiment "
            f"configurtion file {self.yaml_file}."
        )
        self.logger.info(msg=msg)

        # Collect the respective grid type attributes.
        grids_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key="grids", force=True, no_split=True
        )
        if grids_dict is None:
            msg = (
                f"The configuration file {self.yaml_file} does not contain the "
                "attribute grids. Aborting!!!"
            )
            raise ESMFRemapError(msg=msg)

        grid_dict = parser_interface.dict_key_value(
            dict_in=grids_dict, key=grid_name, force=True, no_split=True
        )
        if grid_dict is None:
            msg = (
                f"The grid type {grid_name} could not be determined from the "
                f"grid attributes in configuration file {self.yaml_file}. "
                "Aborting!!!"
            )
            raise ESMFRemapError(msg=msg)

        # Build the respective grid-type object.
        grid_obj = parser_interface.object_define()

        for grid_info_attr in self.grid_info_attr_list:
            value = parser_interface.dict_key_value(
                dict_in=grid_dict, key=grid_info_attr, force=True, no_split=True
            )
            grid_obj = parser_interface.object_setattr(
                object_in=grid_obj, key=grid_info_attr, value=value
            )

        return grid_obj

    def run(self) -> None:
        """
        Description
        -----------

        This method performs the following tasks:

        (1) Computes and defines the remapping attributes for the
            source grid to destination grid interpolation type.

        """

        # Compute and define the ESMF remapping attribute netCDF file
        # path.
        self.compute_remap()
