# =========================================================================

# Module: ush/innov_stats/compute/soca_ocean.py

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

    soca_ocean.py

Description
-----------

    This module contains the base-class module for all Sea-ice and
    Ocean Analysis (SOCA) ocean innovation statistics diagnostics.

Classes
-------

    SOCAOcean(options_obj, basedir)

        This is the base-class object for all Sea-ice and Ocean
        Analysis (SOCA) ocean innovation statistics diagnostics; it is
        a sub-class of InnovStats.

Author(s)
---------

    Henry R. Winterbottom; 28 January 2023

History
-------

    2023-01-28: Henry Winterbottom -- Initial implementation.

"""

# ----

# pylint: disable=eval-used
# pylint: disable=too-many-locals

# ----

from dataclasses import dataclass

import numpy
from ioapps import netcdf4_interface
from tools import parser_interface

from innov_stats.compute import InnovStats

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


@dataclass
class SOCAOcean(InnovStats):
    """
    Description
    -----------

    This is the base-class object for all Sea-ice and Ocean Analysis
    (SOCA) ocean innovation statistics diagnostics; it is a sub-class
    of InnovStats.

    Parameters
    ----------

    options_obj: object

        A Python object containing the command line argument
        attributes.

    basedir: str

        A Python string specifying the top-level of the working
        directory tree.

    """

    def __init__(self, options_obj: object, basedir: str):
        """
        Description
        -----------

        Creates a new SOCAOcean object.

        """

        # Define the base-class attributes.
        super().__init__(options_obj=options_obj, basedir=basedir)

        # Define the SOCA ocean application innovation statistic
        # attributes.
        self.column_frmt = "meters_%03d"
        self.column_scale = 1000.0
        self.get_depths()

        # Initialize the netCDF-formatted output file attributes.
        (self.ncdim_obj, self.ncvar_obj) = [
            parser_interface.object_define() for i in range(2)
        ]

        # Define the total number of vertical levels.
        self.nlevs = len(self.levels_obj.layer_mean)

    def regional_stats(self, ncfilename: str, innovinfo_obj: object) -> None:
        """
        Description
        -----------

        This method computes the innovation statistics from the input
        netCDF formatted file and the respective innovation statistics
        information for the respective variables and defined within
        the user specified experiment configuration; the innovation
        statistics are computed as a function of the regions of
        interest specified in the user experiment configuration; the
        base-class attribute ncvar_obj is updated accordingly.

        Parameters
        ----------

        ncfilename: str

            A Python string specifying the netCDF formatted filename
            path for the respective variable.

        innovinfo_obj: object

            A Python object containing the innovation statistics
            information for the requested variable.

        """

        # Loop through each region of interest to compute the
        # respective innovation statistics; proceed accordingly.
        for region in vars(self.regions_obj):

            # Determine the observation locations that are valid for
            # the specified region of interest.
            msg = f"Computing innovation statistics for region {region}."
            self.logger.info(msg=msg)

            obslocs = self.get_obslocs(
                ncfilename=ncfilename, innovinfo_obj=innovinfo_obj, region=region
            )

            # Quality check the respective observation locations;
            # proceed accordingly.
            qcval = netcdf4_interface.ncreadvar(
                ncfile=ncfilename, ncvarname=innovinfo_obj.ncqc
            )
            qcchk = numpy.logical_and(qcval[obslocs] == 0, qcval[obslocs] <= 0)

            depth = netcdf4_interface.ncreadvar(
                ncfile=ncfilename, ncvarname=innovinfo_obj.ncdepth
            )
            depth = depth[obslocs]
            depth = depth[qcchk]

            # Collect the diagnostic variable "first-guess" values and
            # quality check accordingly.
            omf = netcdf4_interface.ncreadvar(
                ncfile=ncfilename, ncvarname=innovinfo_obj.ncomf
            )
            omf = omf[obslocs]
            omf = omf[qcchk]

            # Compute the innovation statistics for the respective
            # region as a function of depth.
            (bias, count, rmsd) = [numpy.empty(self.nlevs) for i in range(3)]

            for level in range(0, self.nlevs):
                depth_obs = numpy.logical_and(
                    depth >= self.levels_obj.layer_top[level],
                    depth <= self.levels_obj.layer_bottom[level],
                )
                bias[level] = numpy.mean(omf[depth_obs])
                count[level] = numpy.ma.count(omf[depth_obs])
                rmsd[level] = numpy.sqrt(numpy.mean(
                    omf[depth_obs] * omf[depth_obs]))

            # Define that geographical attributes describing the
            # respective region.
            region_info_dict = parser_interface.object_getattr(
                object_in=self.regions_obj, key=region
            )
            (lat_min, lat_max, lon_min, lon_max) = [
                region_info_dict[region_thresh]
                for region_thresh in ["lat_min", "lat_max", "lon_min", "lon_max"]
            ]
            attrs_dict = {
                "lat_min": lat_min,
                "lat_max": lat_max,
                "lon_min": lon_min,
                "lon_max": lon_max,
            }

            # Loop through all innovation statistic attributes and
            # update the base-class netCDF object accordingly.
            for stat in self.stats_type_list:
                ncvar_dict = {
                    f"{stat}_{region}": {
                        "varname": f"{stat}_{region}",
                        "dims": "depth",
                        "type": "float64",
                        "values": eval(stat),
                        "attrs": attrs_dict,
                    }
                }

                self.ncvar_obj = self.update_ncvar(
                    ncvar_dict=ncvar_dict, ncvar_obj=self.ncvar_obj
                )

    def run(self) -> None:
        """
        Description
        -----------

        This method performs the following tasks.

        (1) Initializes the SQLite3 database and corresponding tables
            format.

        (2) Loops through each diagnostics variable specified in the
            user experiment configuration and computes the regional
            innovation statistics as specified within the user
            configuration.

        (3) The innovation statistics for the respective diagnostics
            variable, as a function of region, are written to a
            specified output netCDF formatted file.

        (4) The innovation statistics for the respective diagnostics
            variable are written to a specified SQLite3 database table
            within the user experiment configuration specified SQLite3
            database file path.

        """

        # Initialize the SQLite3 database and tables; proceed
        # accordingly.
        self.build_database(
            column_frmt=self.column_frmt,
            levels_list=self.levels_obj.layer_mean,
            column_scale=self.column_scale,
        )

        # Compute the innovation statistics for each specified
        # diagnostics variable in the experiment configuration.
        for diagsvar in self.diagsvars:

            # Collect the attributes for the respective diagnostic
            # variable.
            vardict = parser_interface.object_getattr(
                object_in=self.diagsinfo_obj, key=diagsvar
            )

            innovinfo_obj = self.get_innovinfo(
                vardict=vardict, addinfo_list=["ncdepth", "ncomf", "ncqc"]
            )

            # Update the base-class object containing the respective
            # netCDF variable attributes.
            dimsdict = {
                "depth": {
                    "size": len(self.levels_obj.layer_mean),
                    "type": "float32",
                    "values": self.levels_obj.layer_mean,
                }
            }

            self.ncvar_obj = self.build_ncdims(
                dimsdict=dimsdict, ncvar_obj=self.ncvar_obj
            )

            # Define netCDF-formatted file path containing the
            # respective diagnostics variable; proceed accordingly.
            ncfilename = self.get_ncfilename(vardict=vardict, ncdim="nlocs")

            if ncfilename is not None:

                # Update the base-class object containing the
                # respective netCDF variable dimension attributes.
                self.ncdim_obj = parser_interface.object_setattr(
                    object_in=self.ncdim_obj, key="depth", value=self.nlevs
                )

                # Compute the respective diagnostic variables
                # innovation statistics as a function of the region(s)
                # specified in the experiment configuration.
                self.regional_stats(ncfilename=ncfilename,
                                    innovinfo_obj=innovinfo_obj)

                # Write the netCDF-formatted file for the respective
                # diagnostics variable and write/update the SQLite3
                # database table(s) for the respective diagnostics
                # variable.
                self.write_ncout(
                    vardict=vardict,
                    variable=diagsvar,
                    ncdim_obj=self.ncdim_obj,
                    ncvar_obj=self.ncvar_obj,
                )
                self.write_database(
                    vardict=vardict,
                    variable=diagsvar,
                    column_frmt=self.column_frmt,
                    column_scale=self.column_scale,
                )
