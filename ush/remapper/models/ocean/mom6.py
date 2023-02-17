# =========================================================================

# Module: ush/remapper/models/ocean/mom6.py

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

    mom6.py

Description
-----------

    This module contains classes and methods to interpolate/regrid an
    external analysis variable(s) to a MOM6 (destination) grid
    projection using the user experiment configuration specified
    interpolation scheme.

Classes
-------

    MOM6(srcgrid_obj, dstgrid_obj, remap_obj, varinfo_obj)

        This is the base-class object for interpolating (regridding)
        user specified variables from a source grid projection to the
        Modular Ocean Model Version 6 (MOM6) forecast model grid
        projection; it is a sub-class of Ocean.

Author(s)
---------

    Henry R. Winterbottom; 17 February 2023

History
-------

    2023-02-17: Henry Winterbottom -- Initial implementation.

"""

# ----

# pylint: disable=too-many-instance-attributes

# ----

from typing import Tuple

import numpy
from exceptions import RemapperError
from ioapps import netcdf4_interface
from remapper.models.ocean import Ocean
from remapper.remap import maskupdate
from remapper.remapio import xarray_interface
from tools import parser_interface

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


class MOM6(Ocean):
    """
    Description
    -----------

    This is the base-class object for interpolating (regridding) user
    specified variables from a source grid projection to the Modular
    Ocean Model Version 6 (MOM6) forecast model grid projection; it is
    a sub-class of Ocean

    Parameters
    ----------

    dstgrid_obj: object

        A Python object containing the destination grid attributes
        collected from the user experiment configuration.

    srcgrid_obj: object

        A Python object containing the source grid attributes
        collected from the user experiment configuration.

    remap_obj: object

        A Python object containing the interpolation (e.g.,
        regridding) attributes collected from the user experiment
        configuration.

    varinfo_obj: object

        A Python object containing the variable attributes for the
        variables to be interpolated (regridded) from a source grid
        projection to a destination grid projection.

    """

    def __init__(
        self: Ocean,
        srcgrid_obj: object,
        dstgrid_obj: object,
        remap_obj: object,
        varinfo_obj: object,
    ):
        """
        Description
        -----------

        Creates a new MOM6 object.

        """

        # Define the base-class attributes.
        super().__init__()

        [self.dstgrid_obj, self.remap_obj, self.srcgrid_obj, self.varinfo_obj] = [
            dstgrid_obj,
            remap_obj,
            srcgrid_obj,
            varinfo_obj,
        ]
        self.interp_scheme = self.remap_obj.interp_scheme.lower()
        self.nlevs = self.remap_obj.nlevs

        self.output_netcdf = parser_interface.object_getattr(
            object_in=self.dstgrid_obj, key="output_netcdf"
        )

        self.variable_list = parser_interface.object_getattr(
            object_in=self.varinfo_obj, key="variable_list", force=True
        )
        if self.variable_list is None:
            msg = (
                "The formatted list containing the list of variables "
                "within the user experiment configuration which to "
                "be interpolated (regridded) could not be determined "
                "from the user experiment configuration. Aborting!!!"
            )
            raise RemapperError(msg=msg)

        # Initialize all regridding/remapping objects.
        self.maskupdate = maskupdate.MaskUpdate()

        self.remap_app = self.get_interpscheme(interp_scheme=self.interp_scheme)

        self.bathy_obj = self.build_bathy(varinfo_obj=self.varinfo_obj)

        self.landmask_obj = self.build_mom6_landmask(
            srcgrid_obj=self.srcgrid_obj,
            dstgrid_obj=self.dstgrid_obj,
            remap_obj=self.remap_obj,
            nlevs=self.nlevs,
        )

        self.build_mom6_ncoutput(
            dstgrid_obj=self.dstgrid_obj,
            varinfo_obj=self.varinfo_obj,
            nlevs=self.nlevs,
            variable_list=self.variable_list,
            output_netcdf=self.output_netcdf,
        )

    def get_uvvars(self: Ocean) -> object:
        """
        Description
        -----------

        This method collects the velocity vector components from the
        respective netCDF-formatted file(s) containing the input
        fields.

        Returns
        -------

        dst_var_obj: object

            A Python object containing the attributes for the
            respective velocity vector components to be remapped.

        """

        # Define the velocity vector components.
        dst_var_obj = parser_interface.object_define()
        dst_ncvarname_list = ["u", "v"]

        # Loop through each variable and remap accordingly.
        for variable in self.variable_list:

            # Define the attributes for the respective variable;
            # proceed accordingly.
            varinfo_dict = parser_interface.object_getattr(
                object_in=self.varinfo_obj, key=variable, force=True
            )
            if varinfo_dict is None:
                msg = (
                    f"The attributes for variable {variable} could not be "
                    "determined from the user experiment configuration. "
                    "Aborting!!!"
                )
                raise RemapperError(msg=msg)

            dst_ncvarname = parser_interface.dict_key_value(
                dict_in=varinfo_dict, key="dst_ncvarname", force=True, no_split=True
            )

            # Build a Python dictionary containing the attributes for
            # the respective momentum (e.g., current velocity)
            # variable(s) to be remapped.
            if dst_ncvarname in dst_ncvarname_list:

                # Read the variable from the netCDF-formatted input
                # file.
                dst_vardict = {}
                msg = f"Reading MOM6 momentum variable {dst_ncvarname}."
                self.logger.info(msg=msg)

                (ncattr_obj, _) = self.build_mom6_varattrobjs(varinfo_dict=varinfo_dict)
                ncvar = netcdf4_interface.ncreadvar(
                    ncfile=ncattr_obj.ncfilename,
                    ncvarname=ncattr_obj.src_ncvarname,
                    squeeze=True,
                    axis=0,
                )
                dst_vardict["ncvar"] = ncvar

                # Define the grid-type staggering for the respective
                # variable.
                grid_stagger = parser_interface.dict_key_value(
                    dict_in=varinfo_dict, key="grid_stagger", force=True, no_split=True
                )
                if grid_stagger is None:
                    msg = (
                        "The grid staggering attributes for MOM6 momentum "
                        f"variable {dst_ncvarname} could not be determined from the user "
                        "experiment configuration. Aborting!!!"
                    )
                    raise RemapperError(msg=msg)

                dst_vardict["grid_stagger"] = grid_stagger
                dst_var_obj = parser_interface.object_setattr(
                    object_in=dst_var_obj, key=dst_ncvarname, value=dst_vardict
                )

        return dst_var_obj

    def mass2uv(
        self: Ocean, ucurr: numpy.array, vcurr: numpy.array, grid_obj: object
    ) -> Tuple[numpy.array, numpy.array]:
        """
        Description
        -----------

        This method remaps the velocity vector components from the
        mass grid locations to the respective staggered grid
        locations.

        Parameters
        ----------

        ucurr: array-type

            A Python array-type variable containing the zonal velocity
            vector component defined a the mass grid locations.

        vcurr: array-type

            A Python array-type variable containing the meridional
            velocity vector component defined a the mass grid
            locations.

        grid_obj: object

            A Python object containing the contents of a xarray
            formatted dataset.

        Returns
        -------

        ucurr: array-type

            A Python array type containing the zonal velocity vector
            component remapped to the respective staggered-grid
            locations.

        vcurr: array-type

            A Python array type containing the meridional velocity
            vector component remapped to the respective staggered-grid
            locations.

        """

        # Remap the velocity vector components from the mass grid to
        # the respective staggered grid locations.
        msg = (
            "Regridding MOM6 current velocity vector components from mass variable "
            "locations to the respective staggered grid points locations."
        )
        self.logger.info(msg=msg)

        ucurr = self.remap_app(
            invar=ucurr,
            interp_type="bilinear",
            srcgrid_stagger="mass",
            dstgrid_stagger="uvel",
            dstgrid_obj=grid_obj,
            srcgrid_obj=grid_obj,
            remap_obj=self.remap_obj,
            reuse_weights=True,
            dst2dst=True,
        )

        vcurr = self.remap_app(
            invar=vcurr,
            interp_type="bilinear",
            srcgrid_stagger="mass",
            dstgrid_stagger="vvel",
            dstgrid_obj=grid_obj,
            srcgrid_obj=grid_obj,
            remap_obj=self.remap_obj,
            reuse_weights=True,
            dst2dst=True,
        )

        return (ucurr, vcurr)

    def regrid_massvars(self: Ocean) -> None:
        """
        Description
        -----------

        This method interpolates (regrids) the user-specified
        variables defined on the respective Arakawa grid
        mass-coordinate from the source grid projection to the
        destination grid projection using the user-specified
        interpolation scheme.

        Raises
        ------

        RemapperError:

            * raised if the attributes for a variable to be remapped
              cannot be determined from the experiment configuration
              file.

        """

        # Regrid the mass-variables; proceed accordingly.
        msg = "Regridding variables defined on the Arakawa mass grid."
        self.logger.info(msg=msg)

        # Loop through each variable and remap accordingly.
        for variable in self.variable_list:

            # Define the attributes for the respective variable;
            # proceed accordingly.
            varinfo_dict = parser_interface.object_getattr(
                object_in=self.varinfo_obj, key=variable, force=True
            )
            if varinfo_dict is None:
                msg = (
                    f"The attributes for variable {variable} could not be "
                    "determined from the experiment configuration file. "
                    "Aborting!!!"
                )
                raise RemapperError(msg=msg)

            grid_stagger = parser_interface.dict_key_value(
                dict_in=varinfo_dict, key="grid_stagger", force=True, no_split=True
            )

            # If the grid-staggering is for a variable defined at the
            # respective grid-type mass locations, proceed
            # accordingly.
            if grid_stagger.lower() == "mass":
                msg = f"Preparing to regrid variable {variable}."
                self.logger.info(msg=msg)

                (ncattr_obj, varattr_obj) = self.build_mom6_varattrobjs(
                    varinfo_dict=varinfo_dict
                )
                msg = (
                    f"Regridding variable {variable} using interpolation "
                    f"type {varattr_obj.interp_type}."
                )
                self.logger.info(msg=msg)

                # Read the variable from the netCDF-formatted input
                # file.
                ncvar = netcdf4_interface.ncreadvar(
                    ncfile=ncattr_obj.ncfilename,
                    ncvarname=ncattr_obj.src_ncvarname,
                    squeeze=True,
                    axis=0,
                )

                invar = numpy.where(ncvar >= 1.0e10, 0.0, ncvar)

                # Remap the source variable to the destination grid
                # projection.
                outvar = self.remap_app(
                    invar=invar,
                    interp_type=varattr_obj.interp_type,
                    srcgrid_stagger=varattr_obj.grid_stagger,
                    dstgrid_stagger=varattr_obj.grid_stagger,
                    dstgrid_obj=self.dstgrid_obj,
                    srcgrid_obj=self.srcgrid_obj,
                    remap_obj=self.remap_obj,
                )

                # Update the remapped variable relative to the
                # destination grid landmask.
                msg = f"Preparing to update variable {variable} relative to the specified landmask."
                self.logger.info(msg=msg)

                outvar = self.maskupdate.run(
                    invar=outvar,
                    mask_obj=self.landmask_obj,
                    zdim_name=varattr_obj.zdim_name,
                    **self.gridfill_kwargs,
                )

                # Write the interpolated variable to the output
                # netCDF-formatted file path.
                msg = (
                    f"Writing interpolated variable {variable} to output netCDF file "
                    f"{self.dstgrid_obj.output_netcdf}."
                )
                self.logger.info(msg=msg)

                var_obj = self.build_varobj(ncvarname=ncattr_obj.dst_ncvarname)
                xarray_interface.write(
                    ncfile=self.output_netcdf, var_obj=var_obj, var_arr=outvar
                )

    def regrid_momentumvars(self: Ocean) -> None:
        """
        Description
        -----------

        This method regrids the MOM6 momentum variables as follows:

        (1) The zonal- and meridional-current velocites are retrieved
            from the respective netCDF formatted file(s) for the
            source grid projection.

        (2) The zonal- and meridional-current velocities are
            interpolated to the source grid projection mass coordinate
            locations.

        (3) The interpolated zonal- and meridional-current velocities
            are rotated from the grid relative projection to an Earth
            relative projection (if necessary).

        (4) The rotated (if necessary) zonal- and meridional-current
            velocities are interpolated to the destination grid
            projection mass coordinate locations.

        (5) The interpolated zonal- and meridional-current velocities,
            now defined on the destination grid projection, are
            rotated if necessary from an Earth relative projection to
            a grid relative projection.

        (6) The interpolated and rotated (if necessary) zonal- and
            meridional current velocities are then interpolated to
            their respective positions along the Arakawa-C grid
            staggering.

        Note that the respective current velocity components must only
        be rotated if they are defined at locations other than the
        mass coordinate variables. The resulting interpolated
        (regridded) momentum variables are updated within the output
        netCDF file.

        """

        # Collect the momentum (e.g., current velocity) variable
        # attributes.
        dst_var_obj = self.get_uvvars()

        # Remap the respective source grid momentum variables to the
        # destination grid mass point locations.
        (umass, vmass) = self.uv2mass(dst_var_obj=dst_var_obj)

        # Rotate the velocity vector components to Earth relative
        # components.
        (ucurr, vcurr) = self.rotate_currents(
            grid_obj=self.srcgrid_obj, ucurr=umass, vcurr=vmass, earth_rel=True
        )

        # Interpolate the velocity vector components, defined at the
        # source grid mass points, to the destination grid mass
        # points.
        umass = self.remap_app(
            invar=ucurr,
            interp_type="bilinear",
            srcgrid_stagger="mass",
            dstgrid_stagger="mass",
            dstgrid_obj=self.dstgrid_obj,
            srcgrid_obj=self.srcgrid_obj,
            remap_obj=self.remap_obj,
            reuse_weights=True,
        )
        vmass = self.remap_app(
            invar=vcurr,
            interp_type="bilinear",
            srcgrid_stagger="mass",
            dstgrid_stagger="mass",
            dstgrid_obj=self.dstgrid_obj,
            srcgrid_obj=self.srcgrid_obj,
            remap_obj=self.remap_obj,
            reuse_weights=True,
        )

        # Rotate the velocity vector components to grid relative
        # components.
        (ucurr, vcurr) = self.rotate_currents(
            grid_obj=self.dstgrid_obj, ucurr=umass, vcurr=vmass, grid_rel=True
        )

        # Interpolate the destination grid current velocity vector
        # components from the mass point locations to the respective
        # staggered locations.
        (ucurr, vcurr) = self.mass2uv(
            ucurr=ucurr, vcurr=vcurr, grid_obj=self.dstgrid_obj
        )

        # Update the netCDF-formatted file path.
        msg = (
            "Updating the current velocity vector components in "
            f"{self.dstgrid_obj.output_netcdf}."
        )
        self.logger.info(msg=msg)

        var_obj = self.build_varobj(ncvarname="u")
        xarray_interface.write(
            ncfile=self.dstgrid_obj.output_netcdf, var_obj=var_obj, var_arr=ucurr
        )

        var_obj = self.build_varobj(ncvarname="v")
        xarray_interface.write(
            ncfile=self.dstgrid_obj.output_netcdf, var_obj=var_obj, var_arr=vcurr
        )

    def uv2mass(self: Ocean, dst_var_obj: object) -> Tuple[numpy.array, numpy.array]:
        """
        Description
        -----------

        This method remaps the velocity vector components from the
        respective staggered-grid locations to the mass-grid
        locations.

        Parameters
        ----------

        dst_var_obj: object

            A Python object containing the destination variable
            attributes.

        Returns
        -------

        umass: array-type

            A Python array-type variable containing the zonal velocity
            vector component, defined at the respective staggered-grid
            locations, remapped to the mass-grid locations.

        vmass: array-type

            A Python array-type variable containing the zonal velocity
            vector component, defined at the respective staggered-grid
            locations, remapped to the mass-grid locations.

        """

        # For each velocity vector component, proceed accordingly to
        # remap the component from the respective staggered-grid
        # location to the mass-grid locations.
        for dst_ncvarname in vars(dst_var_obj):

            # Collect the grid attributes for the respective velocity
            # vector component.
            dst_vardict = parser_interface.object_getattr(
                object_in=dst_var_obj, key=dst_ncvarname
            )
            grid_stagger = parser_interface.dict_key_value(
                dict_in=dst_vardict, key="grid_stagger", no_split=True
            )

            # Define the respective velocity component values at the
            # respective staggered-grid locations.
            ncvar = parser_interface.dict_key_value(dict_in=dst_vardict, key="ncvar")
            invar = numpy.where(ncvar >= 1.0e10, 0.0, ncvar)

            # Remap the respective velocity component values from the
            # respective staggered-grid locations to the mass-grid
            # locations.
            msg = f"Regridding MOM6 momentum variable {dst_ncvarname.lower()}."
            self.logger.info(msg=msg)
            outvar = self.remap_app(
                invar=invar,
                interp_type="bilinear",
                srcgrid_stagger=grid_stagger,
                dstgrid_stagger="mass",
                dstgrid_obj=self.srcgrid_obj,
                srcgrid_obj=self.srcgrid_obj,
                remap_obj=self.remap_obj,
                reuse_weights=True,
                src2src=True,
            )

            # Define the velocity vector component arrays accordingly.
            if dst_ncvarname.lower() == "u":
                umass = outvar

            if dst_ncvarname.lower() == "v":
                vmass = outvar

        return (umass, vmass)

    def run(self: Ocean):
        """
        Description
        -----------

        This method performs the following tasks:

        (1) Establishes the interpolation scheme to be applied (from
            the user experiment configuration) and defines the
            base-class object to be used to interpolate the variables
            specified in the user experiment configuration.

        (2) Defines the base-class objects containing the user
            experiment configuration attributes.

        (3) Builds the output netCDF formatted file to contain the
            interpolated variables.

        (4) Interpolates (regrids) the user specified variables to the
            specified destination grid projection.

        """

        # Remap the mass variables accordingly.
        self.regrid_massvars()

        # Adjust the sea-surface heights and thickness profile
        # relative to the destination grid bathymetry.
        self.bathy_adjust(
            bathy_obj=self.bathy_obj,
            landmask_obj=self.landmask_obj,
            output_netcdf=self.output_netcdf,
        )

        # Remap the respective momentum variables accordingly.
        self.regrid_momentumvars()
