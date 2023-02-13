# =========================================================================

# Module: ush/esmf_remap/cyclic.py

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

    cyclic.py

Description
-----------

    This module contains classes and methods to update input variable
    arrays with various cyclic coordinate values.

Classes
-------

    Cyclic()

        This is the base-class object for updating all input arrays
        with cyclic variable columns and subsequently update all
        missing data values; it is a sub-class of Remap.

Requirements
------------

- cartopy; https://github.com/SciTools/cartopy

- ufs_pytils; https://github.com/HenryWinterbottom-NOAA/ufs_pyutils

Author(s)
---------

   Henry R. Winterbottom; 23 March 2022

History
-------

   2022-03-23: Henry Winterbottom -- Initial implementation.

"""

# ----

# pylint: disable=too-many-function-args

# ----

import numpy
import scipy.ndimage
from cartopy.util import add_cyclic_point
from exceptions import RemapperError

from remap import Remap

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


class Cyclic(Remap):
    """
    Description
    -----------

    This is the base-class object for updating all input arrays with
    cyclic variable columns and subsequently update all missing data
    values; it is a sub-class of Remap.

    """

    def __init__(self: Remap):
        """
        Description
        -----------

        Creates a new Cyclic object.

        """

        # Define the base-class attributes.
        super().__init__(self)

    def add_cyclic_loncoord(
        self: Remap, invar: numpy.array, nlevs: int = None
    ) -> numpy.array:
        """
        Description
        -----------

        This method adds a right-most longitude variable column to the
        input variable on entry.

        Parameters
        ----------

        invar: array-type

            A Python array-type variable containing the input variable
            for which to append a right-most longitude column.

        Keywords
        --------

        nlevs: int, optional

            A Python integer specifying the total number of vertical
            levels; if NoneType on entry, a 2-dimensional variable is
            assumed.

        Returns
        -------

        outvar: array-type

            A Python array-type variable containing the input variable
            with the appended right-most longitude column.

        Raises
        ------

        RemapperError:

            * raised if an exception is encountered while
              computing/defining the cyclic longitude.

        """

        # Initialize the output variable array.
        outvar = invar

        # Define the right-most longitude coordinate dimension column
        # for a 2-dimensional grid provided upon entry.
        if nlevs is None:
            msg = (
                "Adding a cyclic longitude coordinate value for "
                "2-dimensional input variable."
            )
            self.logger.info(msg=msg)

            try:
                tmparr = self.fill(data=add_cyclic_point(invar[:, :], axis=1))
                outvar[:, :] = tmparr[:, 0:-1]

            except Exception as errmsg:
                msg = (
                    "Cyclic longitude coordinate value computation "
                    f"failed with error {errmsg}. Aborting!!!"
                )
                raise RemapperError(msg=msg) from errmsg

        # Define the right-most longitude coordinate dimension column
        # for a 3-dimensional grid provided upon entry.
        if nlevs is not None:
            try:
                for lev in range(nlevs):
                    msg = (
                        "Adding a cyclic longitude coordinate value for "
                        f"3-dimensional input variable level {lev+1}."
                    )
                    self.logger.info(msg=msg)

                    tmparr = self.fill(data=add_cyclic_point(invar[lev, :, :], axis=1))
                    outvar[lev, :, :] = tmparr[:, 0:-1]

            except Exception as errmsg:
                msg = (
                    "Cyclic longitude coordinate value computation "
                    f"failed with error {errmsg}. Aborting!!!"
                )
                raise RemapperError(msg=msg) from errmsg

        return outvar

    def fill(
        self: Remap, data: numpy.array, invalid: numpy.array = None
    ) -> numpy.array:
        """
        Description
        -----------

        This method replaces the value of invalid 'data' cells
        (indicated by 'invalid') by the value of the nearest valid
        data cell

        Parameters
        ----------

        data: array-type

            An array-type variable of any dimension.

        Keywords
        --------

        invalid: array-type, optional

            A binary array of same shape as 'data'; true cells set
            where data value should be replaced; if None (default),
            use: invalid = np.isnan(data)

        Returns
        -------

        data: array-type

            An array-type variable containing the a filled array
            corresponding to the data variable on input.

        """

        # Define the grid of invalid (i.e., NaN) variable values;
        # proceed accordingly.
        if invalid is None:
            invalid = numpy.isnan(data)

        # Update the variable array values accordingly.
        ind = scipy.ndimage.distance_transform_edt(
            invalid, return_distances=False, return_indices=True
        )

        return data[tuple(ind)]
