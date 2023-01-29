# =========================================================================

# Module: ush/innov_stats/gsi_atmos.py

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

    gsi_atmos.py

Description
-----------

    This module contains the base-class module for Gridpoint
    Statistical Interpolation (GSI) atmosphere innovation statistics
    diagnostics.

Classes
-------

    GSIAtmos(options_obj, basedir)

        This is the base-class object for all Gridpoint Statistical
        Interpolation (GSI) atmosphere innovation statistics
        diagnostics; it is a sub-class of InnovStats.


Author(s)
---------

    Henry R. Winterbottom; 28 January 2023

History
-------

    2023-01-28: Henry Winterbottom -- Initial implementation.

"""

# ----

# pylint: disable=broad-except
# pylint: disable=eval-used
# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
# pylint: disable=unused-variable

# ----

from dataclasses import dataclass

import numpy
from ioapps import netcdf4_interface
from tools import parser_interface

from innov_stats import InnovStats, error

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


@dataclass
class GSIAtmos(InnovStats):
    """
    Description
    -----------

    This is the base-class object for all Gridpoint Statistical
    Interpolation (GSI) atmosphere innovation statistics diagnostics;
    it is a sub-class of InnovStats.

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

        Creates a new GSIAtmos object.

        """

        # Define the base-class attributes.
        super().__init__(options_obj=options_obj, basedir=basedir)

        # Define the GSI atmosphere application innovation statistic
        # attributes.
        self.column_frmt = "hPa_%04d"
        self.column_scale = 1.0
        self.get_pressures()
        self.innov_var_list = ["is_spechumid", "is_temperature", "is_wind"]

        # Initialize the netCDF-formatted output file attributes.
        (self.ncdim_obj, self.ncvar_obj) = [
            parser_interface.object_define() for i in range(2)
        ]

    def compute_innovstat(
        self,
        innovinfo_obj: object,
        region: str,
        pres: numpy.array,
        omf: numpy.array,
        ncfilename: str,
        vardict: dict,
        is_wind: bool = False,
    ) -> None:
        """
        Description
        -----------

        This method computes the innovation statistics, as function of
        specified region, using the "first-guess" value (contained in
        the omf numpy array) and updates the base-class netCDF object
        in preparation for the final output file creation.

        Parameters
        ----------

        innovinfo_obj: object

            A Python object containing the innovation statistics
            information for the requested variable.

        region: str

            A Python string specifying the region of interest for
            which to compute the respective innovation statistics.

        pres: numpy.array

            A Python numpy array containing the pressure levels at
            which to bin the respective innovation statistics.

        omf: numpy.array

            A Python numpy array containing the "first-guess" values
            from which to compute the respective innovation
            statistics.

        ncfilename: str

            A Python string specifying the netCDF formatted filename
            path for the respective variable.

        vardict: dict

            A Python dictionary containing the respective variable
            attributes and diagnostics information.

        Keywords
        --------

        is_wind: bool, optional

            A Python boolean valued variable specifying whether the
            "first-guess" variable values are derived from the total
            wind.

        """

        # Define the observation locations within the respective
        # region.
        region_info_dict = parser_interface.object_getattr(
            object_in=self.regions_obj, key=region
        )
        obslocs = self.regional_obs(
            ncfilename=ncfilename,
            innovinfo_obj=innovinfo_obj,
            region=region,
        )
        obslocs = self.usage_obs(
            ncfilename=ncfilename,
            innovinfo_obj=innovinfo_obj,
            vardict=vardict,
            obslocs=obslocs,
            is_wind=is_wind,
        )

        # Compute the innovation statistics for the respective
        # region as a function of pressure level.
        (bias, count, rmsd) = [numpy.empty(self.nlevs) for i in range(3)]

        for level in range(0, self.nlevs):
            preslocs = numpy.logical_and(
                obslocs,
                numpy.logical_and(
                    pres >= self.levels_obj.layer_bottom[level],
                    pres <= self.levels_obj.layer_top[level],
                ),
            )
            bias[level] = numpy.mean(omf[preslocs])
            count[level] = numpy.ma.count(omf[preslocs])
            rmsd[level] = numpy.sqrt(numpy.mean(omf[preslocs] * omf[preslocs]))

        # Define that geographical attributes describing the
        # respective region.
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
                    "dims": "plevs",
                    "type": "float64",
                    "values": eval(stat),
                    "attrs": attrs_dict,
                }
            }
            self.ncvar_obj = self.update_ncvar(
                ncvar_dict=ncvar_dict, ncvar_obj=self.ncvar_obj
            )

    def get_innov_var_bool(self, vardict: dict) -> object:
        """
        Description
        -----------

        This method defines a boolean value specifying the type of
        variable for which the GSI innovation statistics are to be
        computed.

        Parameters
        ----------

        vardict: dict

            A Python dictionary containing the respective variable
            attributes and diagnostics information.

        Returns
        -------

        innov_var_obj: object

            A Python object containing the boolean values for the
            supported variable type innovation statistics.

        Raises
        ------

        InnovStatsError:

            * raised if an error is encountered while constructing the
              Python object containing the boolean values for the
              supported variable type innovation statistics
              (innov_var_obj).

        """

        # Collect the innovation statistic attributes from the
        # experiment configuration and define the local object;
        # proceed accordingly.
        try:
            innov_var_obj = parser_interface.object_define()
            for innov_var in self.innov_var_list:
                value = parser_interface.dict_key_value(
                    dict_in=vardict, key=innov_var, force=True
                )
                if value is None:
                    value = False

                innov_var_obj = parser_interface.object_setattr(
                    object_in=innov_var_obj, key=innov_var, value=value
                )

        except Exception as errmsg:
            msg = (
                "The innovation variable type declaration failed with "
                "error {errmsg}. Aborting!!!"
            )
            error(msg=msg)

        return innov_var_obj

    def innov_spechumid(self, ncfilename: str, vardict: dict) -> None:
        """
        Description
        -----------

        This method computes the innovation statistics metrics for
        specific humidity and updates the base-class attribute
        ncvar_obj as a function of the regions specified in the user
        experiment configuration.

        Parameters
        ----------

        ncfilename: str

            A Python string specifying the netCDF formatted filename
            path for the respective variable.

        vardict: dict

            A Python dictionary containing the respective variable
            attributes and diagnostics information.

        """

        # Read the necessary values from the netCDF-formatted file.
        innovinfo_obj = self.get_innovinfo(
            vardict=vardict,
            addinfo_list=["ncobstype", "ncpres", "ncomf", "ncsatsh", "ncuse"],
        )

        satsh = netcdf4_interface.ncreadvar(
            ncfile=ncfilename, ncvarname=innovinfo_obj.ncsatsh
        )
        omf = netcdf4_interface.ncreadvar(
            ncfile=ncfilename, ncvarname=innovinfo_obj.ncomf
        )
        pres = netcdf4_interface.ncreadvar(
            ncfile=ncfilename, ncvarname=innovinfo_obj.ncpres
        )

        # Scale the specific humidity variable value.
        omf = 100.0 * omf * satsh

        # Loop through each region of interest and proceed
        # accordingly.
        for region in vars(self.regions_obj):

            # Compute the specific humidity innovation statistics for
            # the respective region; update the base-class netCDF
            # object.
            self.compute_innovstat(
                innovinfo_obj=innovinfo_obj,
                region=region,
                pres=pres,
                omf=omf,
                ncfilename=ncfilename,
                vardict=vardict,
            )

    def innov_temperature(self, ncfilename: str, vardict: dict) -> None:
        """
        Description
        -----------

        This method computes the innovation statistics metrics for
        temperature and updates the base-class attribute ncvar_obj as
        a function of the regions specified in the user experiment
        configuration.

        Parameters
        ----------

        ncfilename: str

            A Python string specifying the netCDF formatted filename
            path for the respective variable.

        vardict: dict

            A Python dictionary containing the respective variable
            attributes and diagnostics information.

        """

        # Compute the innovation statistics for the temperature
        # variable; proceed accordingly.

        # Read the necessary values from the netCDF-formatted file.
        innovinfo_obj = self.get_innovinfo(
            vardict=vardict, addinfo_list=["ncobstype", "ncpres", "ncomf", "ncuse"]
        )

        omf = netcdf4_interface.ncreadvar(
            ncfile=ncfilename, ncvarname=innovinfo_obj.ncomf
        )
        pres = netcdf4_interface.ncreadvar(
            ncfile=ncfilename, ncvarname=innovinfo_obj.ncpres
        )

        # Loop through each region of interest and proceed
        # accordingly.
        for region in vars(self.regions_obj):

            # Compute the specific humidity innovation statistics for
            # the respective region; update the base-class netCDF
            # object.
            self.compute_innovstat(
                innovinfo_obj=innovinfo_obj,
                region=region,
                pres=pres,
                omf=omf,
                ncfilename=ncfilename,
                vardict=vardict,
            )

    def innov_wind(self, ncfilename: str, vardict: dict) -> None:
        """
        Description
        -----------

        This method computes the innovation statistics metrics for
        temperature and updates the base-class attribute ncvar_obj as
        a function of the regions specified in the user experiment
        configuration.

        Parameters
        ----------

        ncfilename: str

            A Python string specifying the netCDF formatted filename
            path for the respective variable.

        vardict: dict

            A Python dictionary containing the respective variable
            attributes and diagnostics information.

        """

        # Read the necessary values from the netCDF-formatted file.
        innovinfo_obj = self.get_innovinfo(
            vardict=vardict,
            addinfo_list=["ncobstype", "ncpres", "ncomfu", "ncomfv", "ncuse"],
        )

        omfu = netcdf4_interface.ncreadvar(
            ncfile=ncfilename, ncvarname=innovinfo_obj.ncomfu
        )
        omfv = netcdf4_interface.ncreadvar(
            ncfile=ncfilename, ncvarname=innovinfo_obj.ncomfv
        )
        pres = netcdf4_interface.ncreadvar(
            ncfile=ncfilename, ncvarname=innovinfo_obj.ncpres
        )

        # Compute the total wind speed variable.
        omf = numpy.sqrt(omfu * omfu + omfv * omfv)

        # Loop through each region of interest and proceed
        # accordingly.
        for region in vars(self.regions_obj):

            # Compute the specific humidity innovation statistics for
            # the respective region; update the base-class netCDF
            # object.
            self.compute_innovstat(
                innovinfo_obj=innovinfo_obj,
                region=region,
                pres=pres,
                omf=omf,
                ncfilename=ncfilename,
                vardict=vardict,
                is_wind=True,
            )

    def regional_obs(self, ncfilename: str, innovinfo_obj: object, region: str):
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

        innovinfo_obj: obj

            A Python object containing the innovation statistics
            information for the requested variable.

        region: str

            A Python string specifying the region of interest for
            which to compute the respective innovation statistics.

        """

        msg = f"Computing innovation statistics for region {region}."
        self.logger.info(msg=msg)
        nclat = netcdf4_interface.ncreadvar(
            ncfile=ncfilename, ncvarname=innovinfo_obj.nclat
        )
        nclon = netcdf4_interface.ncreadvar(
            ncfile=ncfilename, ncvarname=innovinfo_obj.nclon
        )

        # Determine the observation locations that are valid for the
        # specified region of interest.
        region_info_obj = parser_interface.object_define()
        region_info_dict = parser_interface.object_getattr(
            object_in=self.regions_obj, key=region
        )

        region_info_list = ["lat_min", "lat_max", "lon_min", "lon_max"]

        for item in region_info_list:
            value = parser_interface.dict_key_value(
                dict_in=region_info_dict, key=item, no_split=True
            )
            region_info_obj = parser_interface.object_setattr(
                object_in=region_info_obj, key=item, value=value
            )

        # Define all valid observation locations and proceed
        # accordingly.
        obslocs = self.get_obslocs(
            ncfilename=ncfilename, innovinfo_obj=innovinfo_obj, region=region
        )
        region_cond = numpy.logical_and(
            numpy.logical_and(
                nclat >= region_info_obj.lat_min, nclat <= region_info_obj.lat_max
            ),
            nclon >= region_info_obj.lon_min,
            nclon <= region_info_obj.lon_max,
        )
        obslocs = numpy.logical_and(obslocs, region_cond)

        return obslocs

    def usage_obs(self, ncfilename, innovinfo_obj, vardict, obslocs, is_wind=False):
        """
        Description
        -----------

        This method updates an array of boolean-values variables
        specifying the valid observations for values to be used in the
        computation of the respective innovation statistics.

        Parameters
        ----------

        ncfilename: str

            A Python string specifying the netCDF formatted filename
            path for the respective variable.

        innovinfo_obj: obj

            A Python object containing the innovation statistics
            information for the requested variable.

        vardict: dict

            A Python dictionary containing the respective variable
            attributes and diagnostics information.

        obslocs: array-type

            A Python array-type variable containing boolean values
            specifying the valid observations to be used for the
            computation of the respective variable innovation
            statistics.

        Keywords
        --------

        is_wind: bool, optional

            A Python boolean valued variable specifying whether the
            observations are wind observations.

        Returns
        -------

        obslocs: array-type

            A Python array-type variable containing boolean values
            specifying the valid observations to be used for the
            computation of the respective variable innovation
            statistics based on the usage flags collected within this
            method.

        Raises
        ------

        InnovStatsError:

            * raised if an error is encountered while determining the
              valid observations to be used for the respective
              variable innovation statistics computations.

        """
        try:
            # Define the observation usage; this step will exclude
            # specific observations based on the GSI observation type.
            varinfo_obj = parser_interface.object_define()
            for item in self.innov_var_list:
                value = parser_interface.dict_key_value(
                    dict_in=vardict, key=item, force=True
                )
                if value is None:
                    value = False

                varinfo_obj = parser_interface.object_setattr(
                    object_in=varinfo_obj, key=item, value=value
                )

            # Determine what observations are to be used for the
            # respective diagnostic variable innovation statistics
            # computations.
            pres = netcdf4_interface.ncreadvar(
                ncfile=ncfilename, ncvarname=innovinfo_obj.ncpres
            )
            ncobstype = netcdf4_interface.ncreadvar(
                ncfile=ncfilename, ncvarname=innovinfo_obj.ncobstype
            )

            ncuse = netcdf4_interface.ncreadvar(
                ncfile=ncfilename, ncvarname=innovinfo_obj.ncuse
            )
            usage = ncuse == 1
            usage = numpy.logical_and(usage, numpy.isfinite(pres))

            # Exclude observation types for the wind innovation
            # statistics computations.
            if is_wind:
                obslocs = numpy.logical_and(ncobstype >= 280, ncobstype <= 282)
                obslocs = numpy.logical_or(
                    obslocs, numpy.logical_or(ncobstype == 220, ncobstype == 221)
                )
                obslocs = numpy.logical_or(
                    obslocs, numpy.logical_and(ncobstype >= 230, ncobstype <= 235)
                )

            obslocs = numpy.logical_and(obslocs, usage)

        except Exception as errmsg:
            msg = (
                "The determination of valid observations to be used for "
                "the computation of the innovation statistics metrics "
                f"failed with error {errmsg}. Aborting!!!"
            )
            error(msg=msg)

        return obslocs

    def run(self):
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

            innov_var_obj = self.get_innov_var_bool(vardict=vardict)

            # Update the base-class object containing the
            # respective netCDF variable attributes.
            dimsdict = {
                "plevs": {
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

                # Compute the specific humidity innovation statistics.
                if innov_var_obj.is_spechumid:
                    self.innov_spechumid(ncfilename=ncfilename, vardict=vardict)

                # Compute the temperature innovation statisitics.
                if innov_var_obj.is_temperature:
                    self.innov_temperature(ncfilename=ncfilename, vardict=vardict)

                # Compute the wind innovation statistics.
                if innov_var_obj.is_wind:
                    self.innov_wind(ncfilename=ncfilename, vardict=vardict)

                # Update the base-class object containing the
                # respective netCDF variable dimension attributes.
                self.ncdim_obj = parser_interface.object_setattr(
                    object_in=self.ncdim_obj, key="plevs", value=self.nlevs
                )

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
