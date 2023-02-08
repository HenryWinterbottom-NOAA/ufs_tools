# =========================================================================

# Module: ush/gridspec/arakawa_c.py

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

    arakawa_c.py

Description
-----------

    This module contains the base-class for all Arakawa-C type reduced
    grid projections from a specified supergrid structure.

Classes
-------

    ArakawaC(options_obj)

        This is the base-class object for all Arakawa-C type reduced
        grid projections generated from a specified supergrid
        structure; it is a sub-class of GridSpec.

Requirements
------------

- ufs_pytils; https://github.com/HenryWinterbottom-NOAA/ufs_pyutils

Author(s)
---------

    Henry R. Winterbottom; 08 February 2023

History
-------

    2023-02-08: Henry Winterbottom -- Initial implementation.

"""

# ----

# pylint: disable=eval-used
# pylint: disable=invalid-name
# pylint: disable=unused-variable

# ----

from exceptions import GridSpecError
from tools import parser_interface

from gridspec import GridSpec

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


class ArakawaC(GridSpec):
    """
    Description
    -----------

    This is the base-class object for all Arakawa-C type reduced grid
    projections generated from a specified supergrid structure; it is
    a sub-class of GridSpec.

    Parameters
    ----------

    options_obj: object

        A Python object containing the command line argument
        attributes.

    """

    def __init__(self: GridSpec, options_obj: object):
        """
        Description
        -----------

        Creates a new ArakawaC object.

        """

        # Define the base-class attributes.
        super().__init__(options_obj=options_obj)

        # Define the reduced grid variable attributes.
        self.reduce_grid_dict = {
            "qlat": ["nyp", "nxp"],
            "qlon": ["nyp", "nxp"],
            "tlat": ["ny", "nx"],
            "tlon": ["ny", "nx"],
            "ulat": ["ny", "nxp"],
            "ulon": ["ny", "nxp"],
            "vlat": ["nyp", "nx"],
            "vlon": ["nyp", "nx"],
        }

    def compute_grid(self: GridSpec) -> None:
        """
        Description
        -----------

        This method defines the Arakawa-C type grid projection from
        the supergrid structure collected from the specified
        netCDF-formatted file path; the Arakawa-C grid projection is
        encapsulated within the base-class object reduce_grid_obj upon
        exit; the base-class netCDF-formatted output file dimension
        attributes are defined within the ncdim_obj object upon exit.

        """

        # Define the netCDF grid-coordinate dimensions.
        nx = (self.grids_obj.mask.ncvalues.shape)[1]
        self.ncdim_obj = parser_interface.object_setattr(
            object_in=self.ncdim_obj, key="nx", value=nx
        )
        self.ncdim_obj = parser_interface.object_setattr(
            object_in=self.ncdim_obj, key="nxp", value=(nx + 1)
        )

        ny = (self.grids_obj.mask.ncvalues.shape)[0]
        self.ncdim_obj = parser_interface.object_setattr(
            object_in=self.ncdim_obj, key="ny", value=ny
        )
        self.ncdim_obj = parser_interface.object_setattr(
            object_in=self.ncdim_obj, key="nyp", value=(ny + 1)
        )

        # Define the grid-cell corner point geographical coordinate
        # values using the supergrid.
        qlat = self.grids_obj.latitude.ncvalues[::2, ::2]
        qlon = self.grids_obj.longitude.ncvalues[::2, ::2]

        # Define the mass variable geographical grid coordinate
        # values using the supergrid.
        tlat = self.grids_obj.latitude.ncvalues[1::2, 1::2]
        tlon = self.grids_obj.longitude.ncvalues[1::2, 1::2]

        # Define the topography grid values; these are defined at the
        # mass variable geographical grid coordinate values of the
        # reduced .
        topog = self.grids_obj.topography.ncvalues[:, :]

        # Define the zonal-velocity variable geographical coordinate
        # values using the supergrid.
        ulat = self.grids_obj.latitude.ncvalues[1::2, ::2]
        ulon = self.grids_obj.longitude.ncvalues[1::2, ::2]

        # Define the meridional-velocity variable geographical
        # coordinate values using the supergrid.
        vlat = self.grids_obj.latitude.ncvalues[::2, 1::2]
        vlon = self.grids_obj.longitude.ncvalues[::2, 1::2]

        # Update the base-class attribute containing the reduced grid
        # attributes.
        for reduce_grid in self.reduce_grid_list:

            self.reduce_grid_obj = parser_interface.object_setattr(
                object_in=self.reduce_grid_obj, key=reduce_grid, value=eval(reduce_grid)
            )

    def prepare_ncfile(self: GridSpec) -> None:
        """
        Description
        -----------

        This method prepares the output netCDF-formatted output file
        containing the reduced grid attributes; the base-class
        attributes ncvar_obj object, containing the reduced grid
        attributes, is defined upon exit.

        """

        # Build the object containing the reduced grid variables to be
        # written to the netCDF-formatted file path.
        for grid_var in list(vars(self.reduce_grid_obj).keys()):

            # Define the netCDF variable attributes.
            dims = parser_interface.dict_key_value(
                dict_in=self.reduce_grid_dict, key=grid_var, force=True, no_split=True
            )
            if dims is None:
                msg = (
                    "The grid dimension variable names could not be determined "
                    f"for grid coordinate variable {grid_var}. Aborting!!!"
                )
                raise GridSpecError(msg=msg)

            ncvar_dict = {
                grid_var: {
                    "varname": grid_var,
                    "dims": dims,
                    "type": "float64",
                    "values": parser_interface.object_getattr(
                        self.reduce_grid_obj, key=grid_var
                    ),
                }
            }

            # Update the base-class netCDF variable object.
            self.ncvar_obj = self.update_ncvar(
                ncvar_dict=ncvar_dict, ncvar_obj=self.ncvar_obj
            )

    def run(self: GridSpec) -> None:
        """
        Description
        -----------

        This method performs the following tasks:

        (1) Collects the netCDF attributes describing the supergrid
            structure from the specified file paths.

        (2) Defines/computes the reduced grid as defined by the
            supergrid structure specified upon entry.

        (3) Prepares the netCDF-formatted output file strutures.

        (4) Writes the netCDF-formatted output file path in accordance
            with the experiment configuaration.

        """

        # Collect the necessary variables from the netCDF-formatted
        # file path.
        self.read_ncfile()

        # Compute the values describing the respective grid-type.
        self.compute_grid()

        # Prepare the base-class attributes for the netCDF-formatted
        # file path creation.
        self.prepare_ncfile()

        # Write the netCDF-formatted file path containing the reduced
        # grid attributes.
        self.write_ncfile()
