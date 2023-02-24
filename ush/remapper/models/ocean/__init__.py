# =========================================================================

# Module: ush/remapper/models/ocean/__init__.py

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

    This module contains the base-class object for all supported ocean
    analysis MOM6 remappings; it is a sub-class of Models.

Classes
-------

    Ocean():

        This is the base-class object for all supported ocean analysis
        MOM6 remappings; it is a sub-class of Models.

Author(s)
---------

    Henry R. Winterbottom; 23 February 2023

History
-------

    2023-02-23: Henry Winterbottom -- Initial implementation.

"""

# ----

# pylint: disable=too-many-arguments

# ----

from typing import Dict, List, Tuple

import numpy
from exceptions import RemapperError
from ioapps import netcdf4_interface
from remapper.models import Models
from remapper.remapio import xarray_interface
from tools import fileio_interface, parser_interface

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


class Ocean(Models):
    """
    Description
    -----------

    This is the base-class object for all supported ocean analysis
    MOM6 remappings; it is a sub-class of Models.

    """

    def __init__(self: Models):
        """
        Description
        -----------

        Creates a new Ocean object.

        """

        # Define the base-class attributes.
        super().__init__()

        self.bathy_adjust_default_dict = {"angstrom_z": 1.0e-10, "htolerance": 1.0e-10}

    def bathy_adjust(
        self: Models, bathy_obj: object, landmask_obj: object, output_netcdf: str
    ) -> None:
        """
        Description
        -----------

        This method adjusts the destination grid in accordance with
        the destination grid bathymetry and writes both the updated
        sea-surface height (sfc) and thickness profile (h) to the
        output netCDF file (base-class attribute output_netcdf); this
        method follows from the MOM6 MOM_state_initialization module
        adjust adjustetatofitbathymetry subroutine.

        Parameters
        ----------

        bathy_obj: object

            A Python object containing the MOM6 destination grid
            bathymetry attributes.

        landmask_obj: object

            A Python object containing the interpolated source grid
            landmask, the destination grid landmask, and the
            destination grid dimensions.

        output_netcdf: str

            A Python string specifying the path to the external netCDF
            formatted file to contain the MOM6 initial conditions.

        """

        # Update the sea-surface height variable using the land/sea
        # mask and update the netCDF-formatted output file.
        msg = "Applying the land/sea mask to the sea-surface height variable."
        self.logger.info(msg=msg)

        sfc = netcdf4_interface.ncreadvar(
            ncfile=output_netcdf, ncvarname="sfc", squeeze=True, axis=0
        )
        sfc = numpy.where(landmask_obj.dstgrid_mask[0, :, :] > 0, sfc, 0.0)
        var_obj = self.build_varobj(ncvarname="sfc")

        xarray_interface.write(ncfile=output_netcdf, var_obj=var_obj, var_arr=sfc)

        # Adjust the interpolated variable thickness profile with
        # respect to the destination grid bathymetry.
        msg = (
            "Adjusting the thickness profile based on the destination grid "
            "bathymetry."
        )
        self.logger.info(msg=msg)

        thck = netcdf4_interface.ncreadvar(
            ncfile=output_netcdf, ncvarname="h", squeeze=True, axis=0
        )
        eta = numpy.zeros([(landmask_obj.nlevs + 1), landmask_obj.nj, landmask_obj.ni])
        eta[0, :, :] = sfc[:, :]

        # The following code block is adopted from the MOM6
        # MOM_state_initialization module adjust
        # adjustetatofitbathymetry subroutine
        for level in range(0, landmask_obj.nlevs):
            eta[(level + 1), :, :] = eta[level, :, :] - thck[level, :, :]

        for level in range(0, (landmask_obj.nlevs + 1)):
            eta[level, :, :] = numpy.where(
                -1.0 * eta[level, :, :]
                > (bathy_obj.depth[:, :] + bathy_obj.htolerance),
                -1.0 * bathy_obj.depth[:, :],
                eta[level, :, :],
            )

        for level in range((landmask_obj.nlevs - 1), 0, -1):
            thck[level, :, :] = numpy.where(
                eta[level, :, :] < (eta[(level + 1), :, :] + bathy_obj.angstrom_z),
                bathy_obj.angstrom_z,
                (eta[level, :, :] - eta[(level + 1), :, :]),
            )

        for xidx in range(landmask_obj.ni):
            for yidx in range(landmask_obj.nj):
                if -eta[landmask_obj.nlevs, yidx, xidx] < (
                    bathy_obj.depth[yidx, xidx] - bathy_obj.htolerance
                ):
                    if eta[0, yidx, xidx] <= eta[landmask_obj.nlevs, yidx, xidx]:
                        thck[:, yidx, xidx] = (
                            eta[1, yidx, xidx] + bathy_obj.depth[yidx, xidx]
                        ) / (float(landmask_obj.nlevs))
                    else:
                        dilate = (eta[0, yidx, xidx] + bathy_obj.depth[yidx, xidx]) / (
                            eta[0, yidx, xidx] - eta[landmask_obj.nlevs, yidx, xidx]
                        )
                        thck[:, yidx, xidx] = thck[:, yidx, xidx] * dilate

        # Apply the land mask to the updated/adjusted thickness
        # profile and update the netCDF-formatted output file.
        msg = "Applying land mask to updated thickness profile."
        self.logger.info(msg=msg)
        thck = numpy.where(landmask_obj.dstgrid_mask > 0, thck, 0.0)

        msg = f"Writing updated thickness profile to netCDF file {output_netcdf}."
        self.logger.info(msg=msg)
        var_obj = self.build_varobj(ncvarname="h")

        xarray_interface.write(ncfile=output_netcdf, var_obj=var_obj, var_arr=thck)

    def bathy_edits(
        self: Models, varinfo_obj: object, depth: numpy.array
    ) -> numpy.array:
        """
        Description
        -----------

        This method edits (e.g., updates) if a bathymetry grid if the
        bathymetry edits file (bathy_edits_file) is defined within the
        experiment configuration.

        Parameters
        ----------

        varinfo_obj: object

            A Python object containing the variable attributes for the
            variables to be interpolated (regridded) from a source
            grid projection to a destination grid projection.

        depth: array-type

            A Python array-type variable containing the bathymetry
            values.

        Returns
        -------

        depth: array-type

            A Python array-type variable containing the updated
            bathymetry values in accordance with the experiment
            configuration; if the bathymetry edits file has not be
            specified in the experiment configuration, this variable
            is equal to the depth array-type values specified upon
            entry.

        Raises
        ------

        RemapperError:

            * raised if the bathymetry edits file path within the
              experiment configuration does not exist.

        """

        # Collect the bathymetry edits attributes; proceed
        # accordingly.
        bathy_edits_file = parser_interface.object_getattr(
            object_in=varinfo_obj, key="bathy_edits_file", force=True
        )
        if bathy_edits_file is None:
            msg = (
                "The user experiment configuration file does not specify a "
                "netCDF file path containing the bathymetry grid updates (i.e., "
                "edits; bathy_edits_file); no bathymetry edits will be made; this may"
                "cause the application to not perform as expected."
            )
            self.logger.warn(msg=msg)

        if bathy_edits_file is not None:
            exist = fileio_interface.fileexist(path=bathy_edits_file)
            if not exist:
                msg = (
                    "The netCDF-formatted file containing the bathymetry edits file "
                    f"{bathy_edits_file} does not exist. Aborting!!!"
                )
                raise RemapperError(msg=msg)

            msg = f"Reading bathymetry edits from file {bathy_edits_file}."
            self.logger.info(msg=msg)

            # Define a Python object containing the bathymetry edits
            # attributes.
            bathy_edits_obj = parser_interface.object_define()
            bathy_edits_attr_list = ["iEdit", "jEdit", "nEdits", "zEdit"]

            for bathy_edits_attr in bathy_edits_attr_list:
                ncvar_obj = xarray_interface.read(
                    ncfile=bathy_edits_file, ncvarname=bathy_edits_attr
                )
                bathy_edits_obj = parser_interface.object_setattr(
                    object_in=bathy_edits_obj, key=bathy_edits_attr, value=ncvar_obj
                )

            # Edit the bathymetry values accordingly.
            for idx in bathy_edits_obj.nEdits.values:
                xcoord = bathy_edits_obj.iEdit[idx].values
                ycoord = bathy_edits_obj.jEdit[idx].values

                input_bathy = depth[ycoord, xcoord]
                depth[ycoord, xcoord] = bathy_edits_obj.zEdit[idx].values

                msg = (
                    f"Resetting input bathymetry coordinate ({ycoord},{xcoord}) "
                    f"value of {input_bathy} to a value of {depth[ycoord, xcoord]}."
                )
                self.logger.warn(msg=msg)

        # Make global bathymetry changes (if applicable).
        depth[:, :] = numpy.where(depth < 9.5, 9.5, depth)
        depth[:, :] = numpy.where(depth > 6500.0, 6500.0, depth)

        return depth

    def build_bathy(self: Models, varinfo_obj: object) -> object:
        """
        Description
        -----------

        This method builds the base-class attribute bathymetry object
        (bathy_obj) containing the MOM6 desination grid
        depth/bathymetry attributes.

        Parameters
        ----------

        varinfo_obj: object

            A Python object containing the variable attributes for the
            variables to be interpolated (regridded) from a source
            grid projection to a destination grid projection.

        Returns
        -------

        bathy_obj: object

            A Python object containing the MOM6 destination grid
            bathymetry attributes.

        Raises
        ------

        RemapperError:

            * raised if the bathymetry file path has not been
              specified within the experiment configuration.

        """

        # Collect the bathymetry attributes.
        bathy_file = parser_interface.object_getattr(
            object_in=varinfo_obj, key="bathy_file", force=True
        )
        if bathy_file is None:
            msg = (
                "The netCDF formatted file containing the bathymetry grid "
                "(bathy_file) could not be determined from the experiment "
                "configuration. Aborting!!!"
            )
            raise RemapperError(msg=msg)

        msg = f"Reading bathymetry from file {bathy_file}."
        self.logger.info(msg=msg)
        depth_obj = xarray_interface.read(ncfile=bathy_file, ncvarname="depth")
        depth = depth_obj.values

        # Edit the bathymetry accordingly.
        depth = self.bathy_edits(varinfo_obj=varinfo_obj, depth=depth)

        # Define a Python object containing the bathymetry attributes.
        bathy_obj = parser_interface.object_define()
        bathy_obj = parser_interface.object_setattr(
            object_in=bathy_obj, key="depth", value=depth
        )

        for bathy_adjust_default in self.bathy_adjust_default_dict:
            value = parser_interface.dict_key_value(
                dict_in=self.bathy_adjust_default_dict,
                key=bathy_adjust_default,
                no_split=True,
            )
            bathy_obj = parser_interface.object_setattr(
                object_in=bathy_obj, key=bathy_adjust_default, value=value
            )

        return bathy_obj

    def rotate_currents(
        self,
        grid_obj: object,
        ucurr: numpy.array,
        vcurr: numpy.array,
        nlevs: int,
        grid_rel: bool = False,
        earth_rel: bool = False,
    ) -> Tuple[numpy.array, numpy.array]:
        """
        Description
        -----------
        This method rotates the ocean current velocity components to
        either an Earth relative projection (from a grid relative
        projection) or to a grid relative projection (from an Earth
        relative projection); the angular component required to rotate
        to the respective projection is collected from the netCDF file
        attribute (ncfile) within the grid_obj argument.

        Parameters
        ----------

        grid_obj; object

            A Python object containing the contents of a xarray
            formatted dataset.

        ucurr: array-type

            A Python 3-dimensional array (z,y,x) containing the
            zonal-current velocity interpolated to respective grid
            mass points.

        vcurr: array-type

            A Python 3-dimensional array (z,y,x) containing the
            meridional-current velocity interpolated to respective
            grid mass points.

        Keywords
        --------

        grid_rel: bool, optional

            A Python boolean valued variable specifying whether to
            rotate the input variables from an Earth relative
            projection to a grid relative projection.

        earth_rel: bool, optional

            A Python boolean valued variable specifying whether to
            rotate the input variables from a grid relative projection
            to an Earth relative projection.

        Returns
        -------

        urot: array-type

            A Python 3-dimensional array (z,y,x) containing the\
            rotated zonal-current velocity.

        vrot: array-type

            A Python 3-dimensional array (z,y,x) containing the
            rotated meridional-current velocity.

        Raises
        ------

        RemapperError:

            * raised if the netCDF-formatted file path containing the
              respective grid projection rotation angle has not been
              defined within the experiment configuration.

            * raised if the netCDF-formatted file path containing the
              respective grid projection rotation angle does not
              exist.

        """

        # Collect the rotation angle for the respective grid
        # projection.
        ncfile = parser_interface.object_getattr(
            object_in=grid_obj, key="grid_ncfile", force=True
        )
        if ncfile is None:
            msg = (
                "The netCDF file containing the rotation angle could not be "
                "be determined from the experiment configuration. Aborting!!!"
            )
            raise RemapperError(msg=msg)

        exist = fileio_interface.fileexist(path=ncfile)
        if not exist:
            msg = f"The netCDF-formatted file path {ncfile} does not exist. Aborting!!!"
            raise RemapperError(msg=msg)

        anglet = netcdf4_interface.ncreadvar(ncfile=ncfile, ncvarname="angle")

        # Rotate the current velocity vector components accordingly.
        if earth_rel:
            msg = (
                "Rotating currents from grid relative to Earth relative projection "
                "values."
            )
            self.logger.info(msg=msg)

            urot = numpy.zeros(numpy.shape(ucurr))
            vrot = numpy.zeros(numpy.shape(vcurr))
            for level in range(0, nlevs):
                urot[level, :, :] = ucurr[level, :, :] * numpy.cos(
                    anglet[:, :]
                ) + vcurr[level, :, :] * numpy.sin(anglet[:, :])
                vrot[level, :, :] = vcurr[level, :, :] * numpy.cos(
                    anglet[:, :]
                ) - ucurr[level, :, :] * numpy.sin(anglet[:, :])

        if grid_rel:
            msg = (
                "Rotating currents from Earth relative to grid relative "
                "projection values."
            )
            self.logger.info(msg=msg)

            urot = numpy.zeros(numpy.shape(ucurr))
            vrot = numpy.zeros(numpy.shape(vcurr))
            for level in range(0, nlevs):
                urot[level, :, :] = ucurr[level, :, :] * numpy.cos(
                    anglet[:, :]
                ) - vcurr[level, :, :] * numpy.sin(anglet)
                vrot[level, :, :] = vcurr[level, :, :] * numpy.cos(
                    anglet[:, :]
                ) + ucurr[level, :, :] * numpy.sin(anglet[:, :])

        return (urot, vrot)
