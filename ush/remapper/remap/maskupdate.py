# =========================================================================

# Module: ush/remapper/remap/maskupdate.py

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

    maskupdate.py

Description
-----------

    This module contains classes and methods to update an interpolated
    variable grid relative to the destination grid projection using
    the associated land/sea masks.

Classes
-------

    MaskUpdate()

        This is the base-class object for updating all interpolated
        variables relative to the destination grid land/sea mask; it
        is a sub-class of Remap.

Requirements
------------

- gridfill; https://github.com/ajdawson/gridfill

- ufs_pytils; https://github.com/HenryWinterbottom-NOAA/ufs_pyutils

Author(s)
---------

   Henry R. Winterbottom; 12 February 2023

History
-------

   2023-02-12: Henry Winterbottom -- Initial implementation.

"""

# ----

# pylint: disable=too-many-function-args

# ----

import gridfill
import numpy

from remap import Remap

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


class MaskUpdate(Remap):
    """
    Description
    -----------

    This is the base-class object for updating all interpolated
    variables relative to the destination grid land/sea mask; it is a
    sub-class of Remap.

    """

    def __init__(self):
        """
        Description
        -----------

        Creates a new MaskUpdate object.

        """

        # Define the base-class attributes.
        super().__init__(self)

    def grid_fill(
        self: Remap,
        invar: numpy.array,
        mask_obj: object,
        zdim_name: str,
        **gridfill_kwargs
    ) -> numpy.array:
        """
        Description
        -----------

        This method uses the gridfill library/application to update
        the interpolated variable grid relative to the destination
        grid land/sea masks.

        Parameters
        ----------

        invar: array-type

            A Python array-type variable containing the interpolated
            variable grid.

        mask_obj: object

            A Python object containing the destination grid land/sea
            mask attributes.

        zdim_name: str

            A Python string specifying the z-coordinate netCDF
            dimension name; if NoneType this method assumes that the
            interpolated variable grid is of 2-dimensions.

        Keywords
        --------

        gridfill_kwargs: dict

            A Python dictionary of keyword arguments for the gridfill
            library/application.

        Returns
        -------

        outvar: array-type

            A Python array-type variable containing the interpolated
            variable grid that has been updated relative to the
            destination grid land/sea masks.

        """

        # Apply the gridfill application to the 2-dimensional
        # variable.
        if zdim_name is None:
            msg = "Updating 2-dimensional variable using 2-dimensional mask."
            self.logger.info(msg)

            try:
                mask = numpy.ma.masked_where(
                    mask_obj.interp_mask[:, :] == 0, invar)

            except IndexError:
                mask = numpy.ma.masked_where(
                    mask_obj.interp_mask[0, :, :] == 0, invar)

            (filled, _) = gridfill.fill(mask, 1, 0, **gridfill_kwargs)
            outvar = filled.reshape([1, mask_obj.nj, mask_obj.ni])

        # Apply the gridfill application to the 3-dimensional
        # variable.
        if zdim_name is not None:
            msg = "Updating 3-dimensional variable using 3-dimensional mask."
            self.logger.info(msg)

            mask = numpy.ma.masked_where(mask_obj.interp_mask == 0, invar)
            (filled, _) = gridfill.fill(mask, 2, 1, **gridfill_kwargs)
            filled = numpy.where(mask_obj.dstgrid_mask == 0.0, 0.0, filled)

            outvar = filled.reshape(
                [1, mask_obj.nlevs, mask_obj.nj, mask_obj.ni])

        return outvar

    def run(
        self: Remap,
        invar: numpy.array,
        mask_obj: object,
        zdim_name: str,
        **gridfill_kwargs
    ) -> numpy.array:
        """
        Description
        -----------

        This method performs the following tasks:

        (1) Updates the interpolated variable relative to the
            destination grid land/sea masks.

        Parameters
        ----------

        invar: array-type

            A Python array-type variable containing the interpolated
            variable grid.

        mask_obj: object

            A Python object containing the destination grid land/sea
            mask attributes.

        zdim_name: str

            A Python string specifying the z-coordinate netCDF
            dimension name; if NoneType this method assumes that the
            interpolated variable grid is of 2-dimensions.

        gridfill_kwargs: dict

            A Python dictionary of keyword arguments for the gridfill
            library/application.

        Returns
        -------

        outvar: array-type

            A Python array-type variable containing the interpolated
            variable grid that has been updated relative to the
            destination grid land/sea masks.

        """

        # Update the interpolated variable relative to the destination
        # grid land/sea masks.
        outvar = self.grid_fill(
            invar=invar, mask_obj=mask_obj, zdim_name=zdim_name**gridfill_kwargs
        )

        return outvar
