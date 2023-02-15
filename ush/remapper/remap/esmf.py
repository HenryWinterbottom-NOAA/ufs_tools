# =========================================================================

# Module: ush/remapper/remap/esmf.py

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

    esmf.py

Description
-----------

    This module contains classes and methods to interpolate source
    grid variables, defined on the respective source grid projection,
    to a destination grid projection via ESMF -- by means of xesmf.

Classes
-------

    ESMF(dstgrid_obj,srcgrid_obj,remap_obj)

        This is the base-class object for all ESMF remapping
        applications via xesmf; it is a sub-class of Remap.

Requirements
------------

- ufs_pytils; https://github.com/HenryWinterbottom-NOAA/ufs_pyutils

- xesmf; https://github.com/pangeo-data/xESMF

Author(s)
---------

    Henry R. Winterbottom; 12 February 2023

History
-------

    2023-02-12: Henry Winterbottom -- Initial implementation.

"""

# ----

# pylint: disable=too-many-arguments
# pylint: disable=too-many-function-args

# ----

import numpy
import xesmf
from exceptions import RemapperError
from remapper.remapio import xarray_interface
from tools import fileio_interface, parser_interface

from remapper.remap import Remap

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


class ESMF(Remap):
    """
    Description
    -----------

    This is the base-class object for all ESMF remapping applications
    via xesmf; it is a sub-class of Remap.

    Parameters
    ----------

    dstgrid_obj: object

        A Python object containing the destination grid attributes
        collected from the user experiment configuration.

    srcgrid_obj: object

        A Python object containing the source grid attributes
        collected from the user experiment configuration.

    remap_obj: object

        A Python object containing the remapping (e.g., regridding)
        attributes collected from the user experiment configuration.

    """

    def __init__(
        self: Remap, dstgrid_obj: object, srcgrid_obj: object, remap_obj: object
    ):
        """
        Description
        -----------

        Creates a new ESMF object.

        """

        # Define the base-class attributes.
        super().__init__(self)
        self.dstgrid_obj = dstgrid_obj
        self.remap_obj = remap_obj
        self.srcgrid_obj = srcgrid_obj

    def get_grid(self: Remap, grid_obj: object, grid_stagger: str):
        """
        Description
        -----------

        This method defines the input data objects (e.g., datasets)
        for the xesmf remapping application.

        Parameters
        ----------

        grid_obj: object

            A Python object containing the grid (either source or
            destination) attributes collected from the user experiment
            configuration.

        grid_stagger: str

            A Python string specifying the grid staggering for grid
            variable; may be either mass, uvel, or vvel.

        Returns
        -------

        grid: object

            A Pyton object specifying the respective source or
            destination grid input data object (e.g., dataset) for the
            xesmf remapping application.

        """

        # Collect the netCDF-formatted file and respective variable
        # array attributes; proceed accordingly.
        ncfile = parser_interface.object_getattr(
            object_in=grid_obj, key="ncfile", force=True
        )
        if ncfile is None:
            msg = (
                "The grid input netCDF file is either NoneType or "
                "undefined. Aborting!!!"
            )
            raise RemapperError(msg=msg)

        grid = xarray_interface.open(ncfile=ncfile)

        # Define the respective grid attributes.
        stagger_coords_dict = {"nclat": "lat", "nclon": "lon"}

        stagger_attrs_dict = parser_interface.object_getattr(
            object_in=grid_obj, key=grid_stagger, force=True
        )
        if stagger_attrs_dict is None:
            msg = (
                "The grid attributes for the mass point locations "
                "could not be determined from the user experiment "
                "configuration. Aborting!!!"
            )
            raise RemapperError(msg=msg)

        grid_ncattrs_dict = {}

        for stagger_coord in stagger_coords_dict:
            stagger_coord_name = parser_interface.dict_key_value(
                dict_in=stagger_attrs_dict, key=stagger_coord, force=True, no_split=True
            )
            if stagger_coord_name is None:
                msg = (
                    f"The grid attribute {stagger_coord} could not be determined "
                    f"for {grid_stagger} staggered variables. Aborting!!!"
                )
                raise RemapperError(msg=msg)

            value = parser_interface.dict_key_value(
                dict_in=stagger_coords_dict, key=stagger_coord, no_split=True
            )
            grid_ncattrs_dict[stagger_coord_name] = value

        # Update the respective grid accordingly.
        grid = grid.rename(**grid_ncattrs_dict)
        grid.close()

        return grid

    def get_remapwghtsfile(
        self: Remap,
        srcgrid_stagger: str,
        dstgrid_stagger: str,
        interp_type: str,
        src2src: bool,
        dst2dst: bool,
    ) -> str:
        """
        Description
        -----------

        This method defines the netCDF file path containing the
        remapping weights for the respective source and destination
        variable grid-staggerings as a function of the interpolation
        type (e.g., bilinear or nearest-neighbor).

        Parameters
        ----------

        srcgrid_stagger: str

            A Python string specifying the respective grid staggering
            for source grid variable; may be either mass, uvel, or
            vvel.

        dstgrid_stagger: str

            A Python string specifying the respective grid staggering
            for destination grid variable; may be either mass, uvel,
            or vvel.

        interp_type: str

            A Python string specifying the interpolation type;
            supported values are either bilinear or nrstnghbr
            (nearest-neighbor).

        src2src: bool

            A Python boolean valued variable specifying whether the
            remapping weights are used to interpolate source grid
            values from one grid-staggering to another.

        dst2dst: bool

            A Python boolean valued variable specifying whether the
            remapping weights are used to interpolate destination grid
            values from one grid-staggering to another.

        Returns
        -------

        remap_ncfile: str

            A Python string specifying the path to the netCDF file
            containing the remapping weights as a function of the
            interpolation type and source and destination grid
            staggerings.

        Raises
        ------

        RemapperError:

            * raised if the netCDF-formatted file path containing the
              remapping attributes has not been defined.

            * raised if the netCDF-formatted file path containing the
              remapping attributes does not exist.

        """

        # Define the remapping type in accordance with the parameter
        # values upon entry.
        if src2src:
            remap_type = f"src{srcgrid_stagger}2src{dstgrid_stagger}_{interp_type}"
        elif dst2dst:
            remap_type = f"dst{srcgrid_stagger}2dst{dstgrid_stagger}_{interp_type}"

        else:
            remap_type = f"src{srcgrid_stagger}2dst{dstgrid_stagger}_{interp_type}"

        # Define the netCDF-formatted file path containing the
        # remapping attributes; proceed accordingly.
        remap_ncfile = parser_interface.object_getattr(
            object_in=self.remap_obj, key=remap_type, force=True
        )

        if remap_ncfile is None:
            msg = (
                f"The remap file {remap_ncfile} has not been defined "
                "within the user experiment configuration. Aborting!!!"
            )
            raise RemapperError(msg=msg)

        exist = fileio_interface.fileexist(path=remap_ncfile)
        if not exist:
            msg = (
                f"The netCDF file {remap_ncfile} containing the ESMF remapping "
                "weights does not exist and/or could not be located. "
                "Aborting!!!"
            )
            raise RemapperError(msg=msg)

        return remap_ncfile

    def remap(
        self: Remap,
        invar: numpy.array,
        srcgrid: object,
        dstgrid: object,
        interp_type: str,
        remap_ncfile: str,
        reuse_weights: bool,
    ) -> numpy.array:
        """
        Description
        -----------

        This method performs the source grid variable remapping to the
        destination grid projection via ESMF.

        Parameters
        ----------

        invar: array-type

            A Python array-type variable containing the input (source
            grid) variable.

        srcgrid: object

            A Python object containing the source grid attributes.

        dstgrid: object

            A Python object containing the destination grid
            attributes.

        interp_type: str

            A Python string specifying the interpolation type;
            supported values are either bilinear or nrstnghbr
            (nearest-neighbor).

        remap_ncfile: str

            A Python string specifying the path to the netCDF file
            containing the remapping weights as a function of the
            interpolation type and source and destination grid
            staggerings.

        reuse_weights: bool

            A Python boolean valued variable specifying whether to use
            an existing netCDF file path containing the ESMF remapping
            attributes; if False, the ESMF remapping attributes will
            be written to a netCDF formatted file specified by the
            value of remap_ncfile (above).

        Returns
        -------

        outvar: array-type

            A Python array type variable containing the input variable
            that has been interpolated/regridded to the destination
            grid projection.

        """

        # Define the remapping method and coefficients attributes for
        # the respective variable; proceed accordingly.
        if reuse_weights:
            msg = (
                "Remapping variable using ESMF remapping attributes "
                f"in file {remap_ncfile}."
            )

        if not reuse_weights:
            msg = (
                "No remapping attributes provided; the regridder "
                "will attempt to compute them and write them to file "
                f"{remap_ncfile}."
            )

        self.logger.info(msg=msg)

        # Remap the variable array.
        remapper = xesmf.Regridder(
            ds_in=srcgrid,
            ds_out=dstgrid,
            method=interp_type,
            periodic=True,
            reuse_weights=reuse_weights,
            filename=remap_ncfile,
        )
        outvar = remapper(invar)

        return outvar

    def run(
        self: Remap,
        interp_type: str,
        srcgrid_stagger: str,
        dstgrid_stagger: str,
        invar: numpy.array,
        reuse_weights: bool = True,
        src2src: bool = False,
        dst2dst: bool = False,
    ) -> numpy.array:
        """
        Description
        -----------

        This method performs the following tasks:

        (1) Defines the source grid data object (e.g., dataset) for
            the xesmf remapping application.

        (2) Defines the destination grid data object (e.g., dataset)
            for the xesmf remapping application.

        (3) Defines the path for the netCDF file containing the ESMF
            remapping weights as function of interpolation type (e.g.,
            bilinear or nearest-neighbor).

        (4) Computes the remapping of the source grid variable from
            the source grid projection to the destination grid
            projection via xesmf.

        Parameters
        ----------

        interp_type: str

            A Python string specifying the interpolation type;
            supported values are either bilinear or nrstnghbr
            (nearest-neighbor).

        srcgrid_stagger: str

            A Python string specifying the respective grid staggering
            for source grid variable; may be either mass, uvel, or
            vvel.

        dstgrid_stagger: str

            A Python string specifying the respective grid staggering
            for destination grid variable; may be either mass, uvel,
            or vvel.

        invar: array-type

            A Python array type variable containing the source grid
            variable defined on the source grid projection.

        Keywords
        --------

        reuse_weights: bool, optional

            A Python boolean valued variable specifying whether to use
            an existing netCDF file path containing the ESMF remapping
            attributes; if False, the ESMF remapping attributes will
            be written to a netCDF formatted file specified by the
            user experiment configuration.

        src2src: bool, optional

            A Python boolean valued variable specifying whether the
            remapping weights are used to interpolate source grid
            values from one grid-staggering to another.

        dst2dst: bool, optional

            A Python boolean valued variable specifying whether the
            remapping weights are used to interpolate destination grid
            values from one grid-staggering to another.

        Returns
        -------

        outvar: array-type

            A Python array type variable containing the input variable
            that has been interpolated/regridded to the destination
            grid projection.

        """

        # Define the respective grid-type attributes.
        dstgrid = self.get_grid(grid_obj=self.dstgrid_obj,
                                grid_stagger=dstgrid_stagger)
        srcgrid = self.get_grid(grid_obj=self.srcgrid_obj,
                                grid_stagger=srcgrid_stagger)

        # Define the file path to contain (or containing) the
        # remapping attributes.
        if reuse_weights:
            remap_ncfile = self.get_remapwghtsfile(
                srcgrid_stagger=srcgrid_stagger,
                dstgrid_stagger=dstgrid_stagger,
                interp_type=interp_type,
                src2src=src2src,
                dst2dst=dst2dst,
            )

        if not reuse_weights:
            remap_ncfile = f"{interp_type}.remap.nc"

        # Remap the respective variable array accordingly.
        outvar = self.remap(
            invar=invar,
            srcgrid=srcgrid,
            dstgrid=dstgrid,
            interp_type=interp_type,
            remap_ncfile=remap_ncfile,
            reuse_weights=reuse_weights,
        )

        return outvar
