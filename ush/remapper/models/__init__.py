# =========================================================================

# Module: ush/remapper/models/__init__.py

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

    This module loads the models package.

Classes
-------

    Models()

        This is the base-class object for all supported forecast model
        initial condition remapping applications.

Requirements
------------

- ufs_pytils; https://github.com/HenryWinterbottom-NOAA/ufs_pyutils

Author(s)
---------

    Henry R. Winterbottom; 25 October 2021

History
-------

    2023-02-12: Henry Winterbottom -- Initial implementation.

"""

# ----

# pylint: disable=too-many-arguments

# ----

import warnings
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy
import xarray
from exceptions import RemapperError
from ioapps import netcdf4_interface
from remapper.remap import esmf, maskupdate
from remapper.remapio import xarray_interface
from tools import fileio_interface, parser_interface
from utils.logger_interface import Logger

warnings.filterwarnings("ignore")

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


@dataclass
class Models:
    """
    Description
    -----------

    This is the base-class object for all supported forecast model
    initial condition remapping applications.

    """

    def __init__(self):
        """
        Description
        -----------

        Creates a new Models object.

        """

        # Define the base-class attributes.
        self.logger = Logger()
        self.maskupdate = maskupdate.MaskUpdate()

        self.gridfill_kwargs = {
            "eps": 1.0e-2,
            "relax": 0.6,
            "itermax": 1500,
            "initzonal": True,
            "cyclic": True,
        }

    def build_varobj(self, ncvarname: str) -> object:
        """
        Description
        -----------

        This method defines a Python object to be used by the xarray
        library/package to define/update the respective netCDF output
        file/variables.

        Parameters
        ----------

        ncvarname: str

            A Python string specifying the netCDF variable name to be
            defined and/or written within the output netCDF file.

        Returns
        -------

        var_obj: object

            A Pyton object containing the netCDF variable attributes
            to be used by xarray library/package to define/update the
            respective netCDF output file/variables.

        """

        # Define the Python object containing the netCDF variable
        # attributes.
        var_obj = parser_interface.object_define()
        var_obj = parser_interface.object_setattr(
            object_in=var_obj, key="ncvarname", value=ncvarname
        )

        return var_obj

    def esmf_remap(
        self,
        invar: numpy.array,
        interp_type: str,
        srcgrid_stagger: str,
        dstgrid_stagger: str,
        srcgrid_obj: object,
        dstgrid_obj: object,
        remap_obj: object,
        reuse_weights: bool = True,
        src2src: bool = False,
        dst2dst: bool = False,
    ) -> numpy.array:
        """
        Description
        -----------

        This method remaps a source grid variable, defined on the
        source grid projection, to a destination grid projection using
        the ESMF algorithm (via xesmf).

        Parameters
        ----------

        invar: array-type

            A Python array-type variable containing the input (source
            grid) variable.

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

            A Python array-type variable containing the output
            (destination grid) variable.

        """

        # Define the remap attributes and the remapper object.
        esmf_remapper = esmf.ESMF(
            srcgrid_obj=srcgrid_obj, dstgrid_obj=dstgrid_obj, remap_obj=remap_obj
        )

        # Remap the respective input variable array.
        outvar = esmf_remapper.run(
            interp_type=interp_type,
            srcgrid_stagger=srcgrid_stagger,
            dstgrid_stagger=dstgrid_stagger,
            invar=invar,
            reuse_weights=reuse_weights,
            src2src=src2src,
            dst2dst=dst2dst,
        )

        return outvar

    def get_interpscheme(self, interp_scheme: str) -> object:
        """
        Description
        -----------

        This method defines the base-class attribute remap_app which
        defines the base-class module to be used to interpolate the
        respective variables; the method is determined by the user
        experiment configuration value for interp_scheme.

        Parameters
        ----------

        interp_scheme: str

            A Python string specifying the interpolation scheme to be
            applied; the supported option is (currently) esmf only.

        Returns
        -------

        remap_app: object

            A Python object specifying the interpolation
            application/methodlogy to be used for the
            remapping/regridding.

        Raises
        ------

        RemapperError:

            * raised when a non-supported scheme/application is
              requested.

        """

        # Define the remapping scheme/application; proceed
        # accordingly.
        if interp_scheme == "esmf":
            remap_app = self.esmf_remap

        if interp_scheme == "slint":
            msg = "The SLINT interpolation scheme is not yet supported. Aborting!!!"
            raise RemapperError(msg=msg)

        return remap_app
