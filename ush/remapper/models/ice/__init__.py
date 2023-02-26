# =========================================================================

# Module: ush/remapper/models/ice/__init__.py

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

    This module contains the base-class object for all supported ice
    analysis CICE remappings; it is a sub-class of Models.

Classes
-------

    Ice():

        This is the base-class object for all supported ice analysis
        CICE remappings; it is a sub-class of Models.

Author(s)
---------

    Henry R. Winterbottom; 25 February 2023

History
-------

    2023-02-25: Henry Winterbottom -- Initial implementation.

"""

# ----

# pylint: disable=invalid-name
# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals

# ----

from typing import List

import numpy
from exceptions import RemapperError
from remapper.models import Models
from remapper.remapio import xarray_interface
from tools import parser_interface

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


class Ice(Models):
    """
    Description
    -----------

    This is the base-class object for all supported ice analysis CICE
    remappings; it is a sub-class of Models.

    """

#    def __init__(self: Models):
#        """
#        Description
#        -----------

#        Creates a new Ice object.

#        """

#        # Define the base-class attributes.
#        super().__init__()

    def build_cfmetadata(
        self: Models, grid_obj: object, dstgrid_obj: object, nlevs: int
    ) -> object:
        """
        Description
        -----------

        This method defines the dimension variable attributes for the
        output netCDF-formatted file assuming the CF metadata
        convention.

        Parameters
        ----------

        grid_obj: object

            A Python object containing the contents of a xarray
            formatted dataset.

        dstgrid_obj: object

            A Python object containing the destination grid attributes
            collected from the experiment configuration.

        nlevs: int

            A Python integer value specifying the total number of
            unstaggered vertical levels for the MOM6 destination grid.

        Returns
        -------

        dims_obj: object

            A Python object containing the dimension variable
            attributes for the output netCDF-formatted file (assuming
            the CF metadata convention).

        Raises
        ------

        RemapperError:

            * raised if the mass-grid is not defined within the
              experiment configuration.

        """

        # Define the attributes for the respective grid types.
        (ncat, time) = (numpy.zeros([nlevs]), numpy.zeros([1]))
        ncat = list(numpy.arange(0, nlevs, 1))

        mass_grid_dict = parser_interface.object_getattr(
            object_in=dstgrid_obj, key="mass", force=True
        )

        if mass_grid_dict is None:
            msg = (
                "The destination mass grid geographical location "
                "attributes could not be determined from the "
                "experiment configuration. Aborting!!!"
            )
            raise RemapperError(msg=msg)

        # Define the mass grid type attributes; proceed accordingly.
        nclat = parser_interface.dict_key_value(
            dict_in=mass_grid_dict, key="nclat", force=True, no_split=True
        )
        nclon = parser_interface.dict_key_value(
            dict_in=mass_grid_dict, key="nclon", force=True, no_split=True
        )

        if any(item is None for item in [nclat, nclon]):
            msg = (
                "A coordinate variable (e.g., nclat or nclon)  could not be "
                "determined for the mass variable grid. Aborting!!!"
            )
            raise RemapperError(msg=msg)

        latm = (
            parser_interface.object_getattr(object_in=grid_obj, key=nclat)
            .mean(axis=1)
            .values
        )
        lonm = (
            parser_interface.object_getattr(object_in=grid_obj, key=nclon)
            .mean(axis=0)
            .values
        )

        # Build the Python object containing the dimension attributes.
        dims_attrs_dict = {
            "latm": latm,
            "lonm": lonm,
            "ncat": ncat,
            "time": time,
        }

        dims_obj = parser_interface.object_define()

        for dims_attr in dims_attrs_dict:
            value = parser_interface.dict_key_value(
                dict_in=dims_attrs_dict, key=dims_attr, no_split=True
            )
            dims_obj = parser_interface.object_setattr(
                object_in=dims_obj, key=dims_attr, value=value
            )

        return dims_obj

    def build_ncoutput(
        self: Models,
        dstgrid_obj: object,
        varinfo_obj: object,
        nlevs: int,
        variable_list: List,
        output_netcdf: str,
    ) -> None:
        """
        Description
        -----------

        This method builds the external netCDF formatted file to
        contain the interpolated input variables.

        Parameters
        ----------

        dstgrid_obj: object

            A Python object containing the destination grid attributes
            collected from the experiment configuration.

        landmask_obj: object

            A Python object containing the interpolated source grid
            landmask, the destination grid landmask, and the
            destination grid dimensions.

        varinfo_obj: object

            A Python object containing the variable attributes for the
            variables to be interpolated (regridded) from a source
            grid projection to a destination grid projection.

        nlevs: int

            A Python integer value specifying the total number of
            unstaggered vertical levels for the CICE destination grid.

        variable_list: list

            A Python list containing the list of specified variables
            to be interpolated (from the experiment configuration).

        output_netcdf: str

            A Python string specifying the path to the external netCDF
            formatted file to contain the CICE initial conditions.

        """

        # Initialize the external netCDF-formatted file to contain the
        # remapped variables.
        msg = f"Preparing output file {output_netcdf} for interpolated CICE variables."
        self.logger.info(msg=msg)

        ncfile = parser_interface.object_getattr(
            object_in=dstgrid_obj, key="grid_ncfile", force=True
        )
        if ncfile is None:
            msg = (
                "The destination grid netCDF formatted file path "
                "could not be determined from the experiment "
                "configuration. Aborting!!!"
            )
            raise RemapperError(msg=msg)

        grid_obj = xarray_interface.open(ncfile=ncfile)
        dims_obj = self.build_cfmetadata(
            grid_obj=grid_obj, dstgrid_obj=dstgrid_obj, nlevs=nlevs
        )
        grid_obj.close()

        varobj_list = []
        for variable in variable_list:
            varinfo_dict = parser_interface.object_getattr(
                object_in=varinfo_obj, key=variable, force=True
            )
            if varinfo_dict is None:
                msg = (
                    "The experiment configuration does not specify "
                    f"the attributes for variable {variable} and/or could not be "
                    "determined from the experiment configuration. Aborting!!!"
                )
                raise RemapperError(msg=msg)

            # Define the mass-grid variable attributes; proceed
            # accordingly.
            if varinfo_dict["grid_stagger"].lower() == "mass":

                coords = {
                    "nj": (["nj"], dims_obj.latm),
                    "ni": (["ni"], dims_obj.lonm),
                    "Time": (["Time"], dims_obj.time),
                }

                (nx, ny) = (len(dims_obj.lonm), len(dims_obj.latm))

                # Build the respective variable array; proceed
                # accordingly.
                try:
                    coords.update(
                        {
                            varinfo_dict["zdim_name"]: (
                                varinfo_dict["zdim_name"],
                                dims_obj.ncat,
                            )
                        }
                    )
                    dims = ["Time", "ncat", "nj", "ni"]

                    msg = (
                        f"Build netCDF variable {variable} of (x,y,z) dimension "
                        f"({nx}, {ny}, {nlevs})."
                    )
                    varval = numpy.zeros([1, nlevs, ny, nx])

                except KeyError:
                    varinfo_dict["zdim_name"] = None
                    dims = ["Time", "nj", "ni"]

                    msg = (
                        f"Build netCDF variable {variable} of (x,y) dimension "
                        f"({nx}, {ny})."
                    )
                    varval = numpy.zeros([1, ny, nx])

                self.logger.info(msg=msg)

                ncvarname = parser_interface.dict_key_value(
                    dict_in=varinfo_dict, key="dst_ncvarname", no_split=True
                )
                varobj = xarray_interface.varobj(
                    varval=varval, coords=coords, dims=dims, ncvarname=ncvarname
                )
                varobj_list.append(varobj)

        # Build the netCDF-formatted output file.
        xarray_interface.dataset(
            ncfile=output_netcdf, varobj_list=varobj_list, unlimitdim="Time"
        )
