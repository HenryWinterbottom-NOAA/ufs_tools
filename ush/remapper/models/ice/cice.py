# =========================================================================

# Module: ush/remapper/models/ice/cice.py

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


# ----

from typing import Dict, Tuple

import numpy
from exceptions import RemapperError

from remapper.models.ice import Ice
from remapper.remap import maskupdate
from remapper.remapio import xarray_interface
from ioapps import netcdf4_interface
from tools import fileio_interface, parser_interface

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


class CICE(Ice):
    """

    """

    def __init__(
        self: Ice,
        srcgrid_obj: object,
        dstgrid_obj: object,
        remap_obj: object,
        varinfo_obj: object,
    ):
        """
        Description
        -----------

        Creates a new CICE object.

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

        self.remap_app = self.get_interpscheme(
            interp_scheme=self.interp_scheme)

        self.landmask_obj = self.build_landmask(
            srcgrid_obj=self.srcgrid_obj,
            dstgrid_obj=self.dstgrid_obj,
            remap_obj=self.remap_obj,
            nlevs=self.nlevs,
        )

        self.build_ncoutput(
            dstgrid_obj=self.dstgrid_obj,
            varinfo_obj=self.varinfo_obj,
            nlevs=self.nlevs,
            variable_list=self.variable_list,
            output_netcdf=self.output_netcdf,
        )

    def build_landmask(self: Ice, srcgrid_obj: object, dstgrid_obj: object,
                       remap_obj: object, nlevs: int) -> object:
        """
        Description
        -----------

        This method builds the base-class landmask object; the input
        variable landmask is interpolated to the output variable
        projection; both the interpolated input and the output grid
        landmasks are projected into 3-dimensions; the interpolated
        input and output variable landmask are defined as
        input_interp_mask and output_mask, respectively, within the
        base-class landmask object.

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

        nlevs: int

            A Python integer value specifying the total number of
            ice levels.

        Returns
        -------

        landmask_obj: object

            A Python object containing the interpolated source grid
            landmask, the destination grid landmask, and the
            destination grid dimensions.

        """

        # Define the source and destination grid masks.
        ncfile = parser_interface.object_getattr(
            object_in=srcgrid_obj, key="topo_ncfile", force=True)
        if ncfile is None:
            msg = (
                "The source grid topography file (topo_ncfile) could not be "
                "determined from the experiment configuration. Aborting!!!"
            )
            raise RemapperError(msg=msg)

        # Build the CICE land/sea mask.
        msg = ('Interpolating the CICE land/sea mask.')
        self.logger.info(msg=msg)

        msg = f"Collecting the source grid mask from netCDF-formatted file {ncfile}."
        self.logger.info(msg=msg)
        srcgrid_mask_obj = xarray_interface.read(
            ncfile=ncfile, ncvarname="wet")

        ncfile = parser_interface.object_getattr(
            object_in=dstgrid_obj, key="topo_ncfile", force=True
        )
        if ncfile is None:
            msg = (
                "The destination grid topography file (topo_ncfile) could not be "
                "determined from the experiment configuration. Aborting!!!"
            )
            raise RemapperError(msg=msg)

        msg = (
            f"Collecting the destination grid mask from netCDF-formatted file {ncfile}."
        )
        self.logger.info(msg=msg)
        dstgrid_mask_obj = xarray_interface.read(
            ncfile=ncfile, ncvarname="wet")

        # Update the respective source and destination grid masks
        # accordingly.
        msg = "Resetting/updating values within source grid land mask."
        self.logger.info(msg=msg)
        srcgrid_mask = srcgrid_mask_obj.values
        srcgrid_mask = numpy.where(srcgrid_mask > 0.0, 1.0, 0.0)

        msg = "Resetting/updating values within destination grid land mask."
        self.logger.info(msg=msg)
        dstgrid_mask = dstgrid_mask_obj.values
        dstgrid_mask = numpy.where(dstgrid_mask > 0.0, 1.0, 0.0)

        # Interpolate the source grid mask to the destination grid
        # projection.
        msg = "Interpolating the source grid mask to the destination grid."
        self.logger.info(msg=msg)
        interp_mask = self.remap_app(
            invar=srcgrid_mask,
            interp_type="bilinear",
            srcgrid_stagger="mass",
            dstgrid_stagger="mass",
            srcgrid_obj=srcgrid_obj,
            dstgrid_obj=dstgrid_obj,
            remap_obj=remap_obj,
        )

        interp_mask = numpy.where(interp_mask < 0.99, 0.0, 1.0)
        (ni, nj, _) = (len(interp_mask[0, :]), len(interp_mask[:, 0]), nlevs)

        # Build the landmask object.
        landmask_obj = parser_interface.object_define()
        landmask_list = ["dstgrid_mask", "interp_mask"]

        for item in landmask_list:
            landmask = numpy.zeros([nlevs, nj, ni])
            landmask_obj = parser_interface.object_setattr(
                object_in=landmask_obj, key=item, value=landmask
            )

        landmask_obj.dstgrid_mask[:, ...] = dstgrid_mask[:, :]
        landmask_obj.interp_mask[:, ...] = interp_mask[:, :]
        (ni, nj, _) = (len(interp_mask[0, :]), len(interp_mask[:, 0]), nlevs)

        # Define the destination grid landmask dimensions.
        grid_dim_list = ["ni", "nj", "nlevs"]

        for grid_dim in grid_dim_list:
            landmask_obj = parser_interface.object_setattr(
                object_in=landmask_obj, key=grid_dim, value=eval(grid_dim)
            )

        return landmask_obj

    def build_varattrobjs(self: Ice, varinfo_dict: Dict) -> Tuple[object, object]:
        """
        Description
        -----------

        This method collects the attributes for the CICE variables to
        be interpolated (regridded) from the experiment configuration
        file.

        Parameters
        ----------

        varinfo_dict: dict

            A Python dictionary containing the experiment
            configuration CICE variable attributes.

        Returns
        -------

        ncattr_obj: object

            A Python object containing the netCDF file attributes
            collected from the experiment configuration.

        varattr_obj: object

            A Python object containing the CICE variable attributes
            collected from the experiment configuration.

        Raises
        ------

        RemapperError:

            * raised if a netCDF attribute cannot be determined from
              the experiment configuration.

            * raised if a specified netCDF-formatted file path does
              not exist.

            * raised if a variable attribute cannot be determined from
              the experiment configuration.

        """

        # Define the Python objects containing the netCDF-formatted
        # file and CICE variable attributes; proceed accordingly.
        (ncattr_obj, varattr_obj) = [
            parser_interface.object_define() for i in range(2)]

        # Build the netCDF-formatted file attributes.
        ncattrs_list = ["ncfilename", "src_ncvarname", "dst_ncvarname"]

        for ncattr in ncattrs_list:
            value = parser_interface.dict_key_value(
                dict_in=varinfo_dict, key=ncattr, force=True, no_split=True
            )
            if value is None:
                msg = (
                    f"The netCDF attribute {ncattr} could not be determined from "
                    "the experiment configuration. Aborting!!!"
                )
                raise RemapperError(msg=msg)

            ncattr_obj = parser_interface.object_setattr(
                object_in=ncattr_obj, key=ncattr, value=value
            )
            exist = fileio_interface.fileexist(path=ncattr_obj.ncfilename)
            if not exist:
                msg = f"The netCDF file {ncattr_obj.ncfilename} does not exist. Aborting!!!"
                raise RemapperError(msg=msg)

        # Build the CICE variable attributes.
        varattr_list = ["grid_stagger",
                        "interp_type", "xdim_name", "ydim_name"]

        for varattr in varattr_list:
            value = parser_interface.dict_key_value(
                dict_in=varinfo_dict, key=varattr, force=True, no_split=True
            )
            if value is None:
                msg = (
                    f"The variable attribute {varattr} could not be determined for "
                    "variable from the experiment configuration. Aborting!!!"
                )
                raise RemapperError(msg=msg)

            varattr_obj = parser_interface.object_setattr(
                object_in=varattr_obj, key=varattr, value=value
            )

            value = parser_interface.dict_key_value(
                dict_in=varinfo_dict, key="zdim_name", force=True, no_split=True
            )
            varattr_obj = parser_interface.object_setattr(
                object_in=varattr_obj, key="zdim_name", value=value
            )

        return (ncattr_obj, varattr_obj)

    def regrid_massvars(self: Ice) -> None:
        """
        Description
        -----------

        This method interpolates (regrids) the specified variables
        defined on the respective Arakawa grid mass-coordinate from
        the source grid projection to the destination grid projection
        using the specified interpolation scheme.

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

                (ncattr_obj, varattr_obj) = self.build_varattrobjs(
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
                    ncvarname=ncattr_obj.src_ncvarname
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
                outvar = outvar[0, ...]
                outvar[numpy.isnan(outvar)] = 0.0

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

    def run(self: Ice):
        """

        """

        self.regrid_massvars()
