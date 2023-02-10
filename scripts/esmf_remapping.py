# =========================================================================

# Script: scripts/esmf_remapping.py

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

"""Script
------

    esmf_remapping.py

Description
-----------

    This script is the driver script for Earth System Modeling
    Framework (ESMF) remapping coefficient value computations.

Classes
-------

   ESMFRemapping()

       This is the base-class object to compute the remapping
       attributes required to interpolate (i.e., remap) a variable
       defined on a source grid to a destination grid projection.

   ESMFRemappingError(msg)

       This is the base-class for all exceptions; it is a sub-class of
       Error.

Methods
-------

   main()

      This is the driver-level method to invoke the tasks within this
      script.

Author(s)
---------

   Henry R. Winterbottom; 01 December 2021

History
-------

   2021-12-01: Henry Winterbottom -- Initial implementation.

"""

# ----

from dataclasses import dataclass
import os
import time

from confs.yaml_interface import YAML

import gribio
import ncio

import warnings
import xesmf

from tools import parser_interface
from exceptions import ESMFRemappingError
from utils.arguments_interface import Arguments
from utils.logger_interface import Logger

warnings.filterwarnings("ignore")

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----

# Specify whether to evaluate the format for the respective parameter
# values.
EVAL_SCHEMA = True

# Define the schema attributes.
CLS_SCHEMA = {
    "yaml_file": str
}

# ----


@dataclass
class ESMFRemapping:
    """
    Description
    -----------

    This is the base-class object to compute the remapping attributes
    required to interpolate (i.e., remap) a variable defined on a
    source grid to a destination grid projection.

    """

    def __init__(self, options_obj: object):
        """
        Description
        -----------

        Creates a new ESMFRemapping object.

        """

        # Define the base-class attributes.
        self.logger = Logger()
        self.options_obj = options_obj
        self.yaml_file = self.options_obj.yaml_file
        self.yaml_dict = YAML().read_yaml(yaml_file=self.yaml_file)

        # Initialize the source and destination grid attributes from
        # the experiment configuration file.
        self.srcgrid_info_obj = self.getgrid_info(grid_name="source")
        self.dstgrid_info_obj = self.getgrid_info(grid_name="destination")
        self.input_grids()

    def compute_remap(self) -> None:
        """
        Description
        -----------

        This method computes the ESMF remapping attributes and writes
        and external netCDF formatted file containing the remapping
        attributes for the interpolation type specified in the user
        experiment configuration.

        """

        remap_obj = tools.parser_interface.object_define()
        remap_attrs_list = ['interp_type', 'output_netCDF']
        for remap_attr in remap_attrs_list:
            kwargs = {'dict_in': self.yaml_dict, 'key': remap_attr, 'force':
                      True, 'no_split': True}
            value = tools.parser_interface.dict_key_value(**kwargs)
            if value is None:
                msg = ('The attribute %s could not be determined from the '
                       'user experiment configuration. Aborting!!!' %
                       remap_attr)
                raise ESMFRemappingError(msg=msg)
            kwargs = {'object_in': remap_obj, 'key': remap_attr, 'value':
                      value}
            remap_obj = tools.parser_interface.object_setattr(**kwargs)

        if self.srcgrid_info_obj.is_grib:
            kwargs = {'srcgrid_obj': self.srcgrid_obj}
            vararray = gribio.define_vararray(**kwargs)
        if self.srcgrid_info_obj.is_netcdf:
            kwargs = {'srcgrid_info_obj': self.srcgrid_info_obj}
            vararray = ncio.define_vararray(**kwargs)
        msg = ('Computing the ESMF remapping attributes for interpolation type '
               '%s and writing to netCDF file path %s.' % (remap_obj.interp_type,
                                                           remap_obj.output_netCDF))
        self.logger.info(msg=msg)
        kwargs = {'ds_in': self.srcgrid_obj, 'ds_out': self.dstgrid_obj,
                  'method': remap_obj.interp_type, 'periodic': True,
                  'reuse_weights': False, 'filename': remap_obj.output_netCDF}
        remapper = xesmf.Regridder(**kwargs)
        remapper(vararray)

    def getgrid_info(self, grid_name: str) -> object:
        """
        Description
        -----------

        This method parses the user experiment configuration and
        collects the input and output grid configuration attributes
        and defines the Python objects for the respective source and
        destination grids; the returned object attributes specify
        either WMO GRIB formatted file path or the netCDF variable
        names for the respective geographical coordinate variables
        (e.g., latitude and longitude) in accordance with the
        respective Arakawa grid staggerings as well as the netCDF
        formatted file paths for the respective source and desination
        grids.

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

        """

        # Collect the attributes for the respective source and
        # destination grids from the experiment configuration.
        grid_info_attr_list = ["grib_file",
                               'is_grib',
                               'is_netcdf',
                               'ncfile',
                               'ncxdim_name',
                               'ncydim_name',
                               'nclat',
                               'nclon'
                               ]

        msg = (f'Collecting the {grid_name} grid information from the user '
               'experiment configuration.')
        self.logger.info(msg=msg)

        # Check that the grid type is supported; proceed accordingly.
        matching = [s for s in ['source', 'destination'] if grid_name in s]
        if not matching:
            msg = (f'The specified grid type {grid_name} is not supported or cannot '
                   'be determined from the user experiment configuration. '
                   'Aborting!!!' % grid_name)
            raise ESMFRemappingError(msg=msg)

        # Build the respective grid attributes; proceed accordingly.
        grid_obj = parser_interface.object_define()

        grids_dict = tools.parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key="grids", force=True)
        if grids_dict is None:
            msg = ('The grids attribute could not be determined from the '
                   f'user experiment configuration in YAML file {self.yaml_file}. Aborting!!!'
                   )
            raise ESMFRemappingError(msg=msg)

        grid_dict = parser_interface.dict_key_value(
            dict_in=grids_dict, key=grid_name, force=True)

        if grid_dict is None:
            msg = (f'The grid type {grid_name} could not be determined from the user '
                   'experiment grid attributes. Aborting!!!'
                   )
            raise ESMFRemappingError(msg=msg)

        for grid_info_attr in grid_info_attr_list:
            value = parser_interface.dict_key_value(
                dict_in=grid_dict, key=grid_info_attr, force=True, no_split=True)
            grid_obj = parser_interface.object_setattr(
                object_in=grid_obj, key=grid_info_attr, value=value)

        return grid_obj

    def input_grids(self):
        """
        Description
        -----------

        This method defines the attributes for the respective source
        and destination grids and updates the corresponding base-class
        objects.

        """

        # Collect and define the source grid attributes.
        if self.srcgrid_info_obj.is_grib:
            self.srcgrid_obj = gribio.input_grids(
                grid_obj=self.srcgrid_info_obj)

        if self.srcgrid_info_obj.is_netcdf:
            self.srcgrid_obj = ncio.input_grids(grid_obj=self.srcgrid_info_obj)

        # Collect and define the destination grid attributes.
        if self.dstgrid_info_obj.is_grib:
            self.dstgrid_obj = gribio.input_grids(
                grid_obj=self.dstgrid_info_obj)

        if self.dstgrid_info_obj.is_netcdf:
            self.dstgrid_obj = ncio.input_grids(
                grid_obj=self.dstgrid_info_obj)

    def run(self):
        """
        Description
        -----------

        This method performs the following tasks:

        (1) Collects the source and destination grid information from
            the user experiment configuration.

        (2) Defines the base-class objects for the respective source
            and destination grids.

        (3) Computes the remapping attributes for the respective
            interpolation type using the xESMF application.

        (4) Writes the computed remapping attributes to the user
            specified external netCDF formatted file.

        """
        msg = ('Container application beginning.')
        self.logger.info(msg=msg)

        self.compute_remap()
        msg = ('Container application completed.')
        self.logger.info(msg=msg)

# ----


def main():
    """
    Description
    -----------

    This is the driver-level function to invoke the tasks within this
    script.

    """

    # Collect the command line arguments.
    script_name = os.path.basename(__file__)
    start_time = time.time()
    msg = f"Beginning application {script_name}."
    Logger().info(msg=msg)
    options_obj = Arguments().run(eval_schema=EVAL_SCHEMA, cls_schema=CLS_SCHEMA)

    # Launch the task.
    task = ESMFRemapping(options_obj=options_obj)
    task.run()

    stop_time = time.time()
    msg = f"Completed application {script_name}."
    Logger().info(msg=msg)
    total_time = stop_time - start_time
    msg = f"Total Elapsed Time: {total_time} seconds."
    Logger().info(msg=msg)


# ----

if __name__ == '__main__':
    main()
