# =========================================================================

# $$$ MODULE DOCUMENTATION BLOCK

# UFS-RNR-containers :: innov_stats/py/innov_stats/gsi_atmos.py

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

   GSIAtmos(yaml_file, basedir)
    
       This is the base-class object for all Gridpoint Statistical
       Interpolation (GSI) atmosphere innovation statistics
       diagnostics; it is a sub-class of InnovStats.


Author(s)
---------

   Henry R. Winterbottom; 04 August 2022

History
-------

   2022-08-04: Henry Winterbottom -- Initial implementation.

"""

# ----

import inspect
import numpy
import statistics

from ioapps import netcdf4_interface
from innov_stats import InnovStats
from innov_stats import InnovStatsError
from tools import parser_interface

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


class GSIAtmos(InnovStats):
    """
    Description
    -----------

    This is the base-class object for all Gridpoint Statistical
    Interpolation (GSI) atmosphere innovation statistics diagnostics;
    it is a sub-class of InnovStats.

    Parameters
    ----------

    yaml_file: str

        A Python string specifying the path to the user experiment
        configuration file.

    basedir: str

        A Python string specifying the top-level of the working
        directory tree.

    """

    def __init__(self, yaml_file, basedir):
        """
        Description
        -----------

        Creates a new GSIAtmos object.

        """

        # Define the base-class attributes.
        super(GSIAtmos, self).__init__(
            yaml_file=yaml_file, basedir=basedir)

        self.cycle = self.get_cycle()
        self.exptname = self.get_exptname()
        self.regions_obj = self.get_regions()
        self.diagsvars = self.get_diagsinfo()

        self.column_frmt = 'hPa_%04d'
        self.column_scale = 1.0
        self.get_pressures()
        self.innov_var_list = ['is_spechumid', 'is_temperature', 'is_wind']

    def get_innov_var_bool(self, vardict, varname):
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

        innov_var_obj: obj

            A Python object containing the boolean values for the
            supported variable type innovation statistics.

        Raises
        ------

        InnovStatsError:

            * raised if an error is encountered while constructing the
              Python object containing the boolean values for the
              supported variable type innovation statistics
              (innov_var_obj).

            * raised if multiple instances of True are defined for the
              respective diagnostics variable type declaration within
              the user experiment configuration.

        """
        try:
            innov_var_obj = parser_interface.object_define()
            for innov_var in self.innov_var_list:
                value = parser_interface.dict_key_value(
                    dict_in=vardict, key=innov_var, force=True)
                if value is None:
                    value = False
                innov_var_obj = parser_interface.object_setattr(
                    object_in=innov_var_obj, key=innov_var, value=value)
        except Exception as error:
            msg = ('The innovation variable type declaration failed with '
                   'error {0}. Aborting!!!'.format(error))
            raise InnovStatsError(msg=msg)
        vartype_list = list()
        for innov_var in self.innov_var_list:
            vartype_list.append(parser_interface.object_getattr(
                object_in=innov_var_obj, key=innov_var))
        if sum(vartype_list) > 1:
            method_name = inspect.currentframe().f_code.co_name
            msg = ('The user experiment configuration contains multiple '
                   'declarations for variable {0}. Upon entry the method {1} '
                   'received the following:\n\n'.format(varname, method_name))
            for innov_var in vars(innov_var_obj):
                value = parser_interface.object_getattr(
                    object_in=innov_var_obj, key=innov_var, force=True)
                if value is not None:
                    msg = msg + '{0}: {1}\n'.format(innov_var, value)
            msg = msg + '\nAborting!!!'
            raise InnovStatsError(msg=msg)
        return innov_var_obj

    def get_pressures(self):
        """
        Description
        -----------

        This method collects the layer information from the user
        experiment configuration; the base-class attributes nlevs and
        levels_obj is defined and contains the following:

        layer_bottom: the bottom interface of a given layer; this is
                      collected from the user experiment
                      configuration.

        layer_top: the top interface of a given layer; this is
                   collected from the user experiment configuration.

        layer_mean: this is the mean value of the layer_bottom and
                    layer_top value for a given layer; it is computed
                    within this method.

        Raises
        ------

        InnovStatsError:

            * raised if the pressure levels information cannot be
              determined from the user experiment configuration.

        Notes
        -----

        * The orientation of the layers for the atmosphere column is
          flipped (i.e., the layer_top value is the bottom of the
          layer and layer_bottom is the top of the layer) due to the
          orientation of the atmosphere pressure (isobaric) profile.

        """

        # Define the pressure-level information from the user
        # experiment configuration and proceed accordingly.
        pressure_levels_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key='pressure_levels',
            force=True)
        if pressure_levels_dict is None:
            msg = ('The pressure levels information could not be '
                   'determined from the user experiment configuration. '
                   'Aborting!!!')
            raise InnovStatsError(msg=msg)

        # Define the base-class attributes levels_obj; this namespace
        # contains both the layer bounding regions (i.e., top and
        # bottom layer interfaces for the respective level) as well as
        # the mean value for the respective interval (i.e., the layer
        # mean -- the middle of the layer).
        levels_attr_dict = {'plevs_bottom': 'layer_top', 'plevs_top':
                            'layer_bottom'}
        for levels_attr in levels_attr_dict.keys():
            value = parser_interface.dict_key_value(
                dict_in=pressure_levels_dict, key=levels_attr,
                force=True)
            if value is None:
                msg = ('The attribute {0} could not be determined from '
                       'the user experiment configuration. Aborting!!!'.
                       format(levels_attr))
                raise InnovStatsError(msg=msg)

            levels = [float(level) for level in value]

            self.levels_obj = parser_interface.object_setattr(
                object_in=self.levels_obj, key=levels_attr_dict[levels_attr],
                value=levels)

        value = [statistics.mean(k) for k in zip(self.levels_obj.layer_bottom,
                                                 self.levels_obj.layer_top)]
        self.levels_obj = parser_interface.object_setattr(
            object_in=self.levels_obj, key='layer_mean', value=value)

        # Define the total number of vertical levels.
        self.nlevs = len(self.levels_obj.layer_mean)

    def innov_spechumid(self, ncfilename, vardict):
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

        # Compute the innovation statistics for the specific humidity
        # variable; proceed accordingly.
        innovinfo_obj = self.get_innovinfo(
            vardict=vardict, addinfo_list=['ncobstype', 'ncpres', 'ncomf',
                                           'ncsatsh', 'ncuse'])
        satsh = netcdf4_interface.ncreadvar(ncfile=ncfilename,
                                            ncvarname=innovinfo_obj.ncsatsh)
        omf = netcdf4_interface.ncreadvar(ncfile=ncfilename,
                                          ncvarname=innovinfo_obj.ncomf)
        pres = netcdf4_interface.ncreadvar(ncfile=ncfilename,
                                           ncvarname=innovinfo_obj.ncpres)
        omf = 100.0*omf*satsh

        # Loop through each region of interest and proceed
        # accordingly.
        for region in vars(self.regions_obj):
            (bias, count, rmsd) = [numpy.empty(self.nlevs) for i in range(3)]
            region_info_dict = parser_interface.object_getattr(
                object_in=self.regions_obj, key=region)
            obslocs = self.regional_obs(
                ncfilename=ncfilename, innovinfo_obj=innovinfo_obj, vardict=vardict, region=region)
            obslocs = self.usage_obs(
                ncfilename=ncfilename, innovinfo_obj=innovinfo_obj, vardict=vardict, obslocs=obslocs)
            for i in range(0, self.nlevs):
                preslocs = numpy.logical_and(obslocs, numpy.logical_and(pres >=
                                                                        self.levels_obj.layer_bottom[i],
                                                                        pres <= self.levels_obj.layer_top[i]))
                bias[i] = numpy.mean(omf[preslocs])
                count[i] = numpy.ma.count(omf[preslocs])
                rmsd[i] = numpy.sqrt(numpy.mean(omf[preslocs]*omf[preslocs]))
            (lat_min, lat_max, lon_min, lon_max) = [region_info_dict[region_thresh] for
                                                    region_thresh in ['lat_min', 'lat_max',
                                                                      'lon_min', 'lon_max']]
            attrs_dict = {'lat_min': lat_min, 'lat_max': lat_max, 'lon_min': lon_min,
                          'lon_max': lon_max, }

            # Loop through all innovation statistic attributes and
            # update the base-class netCDF object accordingly.
            for stat in self.stats_type_list:
                ncvar_dict = {'{0}_{1}'.format(stat, region):
                              {'varname': '{0}_{1}'.format(stat, region),
                               'dims': 'plevs', 'type': 'float64', 'values': eval(stat),
                               'attrs': attrs_dict}}
                self.ncvar_obj = self.update_ncvar(
                    ncvar_dict=ncvar_dict, ncvar_obj=self.ncvar_obj)

    def innov_temperature(self, ncfilename, vardict):
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
        innovinfo_obj = self.get_innovinfo(
            vardict=vardict, addinfo_list=['ncobstype', 'ncpres', 'ncomf',
                                           'ncuse'])
        omf = netcdf4_interface.ncreadvar(ncfile=ncfilename,
                                          ncvarname=innovinfo_obj.ncomf)
        pres = netcdf4_interface.ncreadvar(ncfile=ncfilename,
                                           ncvarname=innovinfo_obj.ncpres)

        # Loop through each region of interest and proceed
        # accordingly.
        for region in vars(self.regions_obj):
            (bias, count, rmsd) = [numpy.empty(self.nlevs) for i in range(3)]
            region_info_dict = parser_interface.object_getattr(
                object_in=self.regions_obj, key=region)
            obslocs = self.regional_obs(
                ncfilename=ncfilename, innovinfo_obj=innovinfo_obj, vardict=vardict, region=region)
            obslocs = self.usage_obs(
                ncfilename=ncfilename, innovinfo_obj=innovinfo_obj, vardict=vardict, obslocs=obslocs)
            for i in range(0, self.nlevs):
                preslocs = numpy.logical_and(obslocs, numpy.logical_and(pres >=
                                                                        self.levels_obj.layer_bottom[i],
                                                                        pres <= self.levels_obj.layer_top[i]))
                bias[i] = numpy.mean(omf[preslocs])
                count[i] = numpy.ma.count(omf[preslocs])
                rmsd[i] = numpy.sqrt(numpy.mean(omf[preslocs]*omf[preslocs]))
            (lat_min, lat_max, lon_min, lon_max) = [region_info_dict[region_thresh] for
                                                    region_thresh in ['lat_min', 'lat_max',
                                                                      'lon_min', 'lon_max']]

            # Loop through all innovation statistic attributes and
            # update the base-class netCDF object accordingly.
            attrs_dict = {'lat_min': lat_min, 'lat_max': lat_max, 'lon_min': lon_min,
                          'lon_max': lon_max, }
            for stat in self.stats_type_list:
                ncvar_dict = {'{0}_{1}'.format(stat, region):
                              {'varname': '{0}_{1}'.format(stat, region),
                               'dims': 'plevs', 'type': 'float64', 'values': eval(stat),
                               'attrs': attrs_dict}}
                self.ncvar_obj = self.update_ncvar(
                    ncvar_dict=ncvar_dict, ncvar_obj=self.ncvar_obj)

    def innov_wind(self, ncfilename, vardict):
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
        innovinfo_obj = self.get_innovinfo(
            vardict=vardict, addinfo_list=['ncobstype', 'ncpres', 'ncomfu',
                                           'ncomfv', 'ncuse'])
        omfu = netcdf4_interface.ncreadvar(ncfile=ncfilename,
                                           ncvarname=innovinfo_obj.ncomfu)
        omfv = netcdf4_interface.ncreadvar(ncfile=ncfilename,
                                           ncvarname=innovinfo_obj.ncomfv)
        pres = netcdf4_interface.ncreadvar(ncfile=ncfilename,
                                           ncvarname=innovinfo_obj.ncpres)
        omf = numpy.sqrt(omfu*omfu + omfv*omfv)

        # Loop through each region of interest and proceed
        # accordingly.
        for region in vars(self.regions_obj):
            (bias, count, rmsd) = [numpy.empty(self.nlevs) for i in range(3)]
            region_info_dict = parser_interface.object_getattr(
                object_in=self.regions_obj, key=region)
            obslocs = self.regional_obs(
                ncfilename=ncfilename, innovinfo_obj=innovinfo_obj, vardict=vardict, region=region)
            obslocs = self.usage_obs(
                ncfilename=ncfilename, innovinfo_obj=innovinfo_obj, vardict=vardict, obslocs=obslocs, is_wind=True)
            for i in range(0, self.nlevs):
                preslocs = numpy.logical_and(obslocs, numpy.logical_and(pres >=
                                                                        self.levels_obj.layer_bottom[i],
                                                                        pres <= self.levels_obj.layer_top[i]))
                bias[i] = numpy.mean(omf[preslocs])
                count[i] = numpy.ma.count(omf[preslocs])
                rmsd[i] = numpy.sqrt(numpy.mean(omf[preslocs]*omf[preslocs]))
            (lat_min, lat_max, lon_min, lon_max) = [region_info_dict[region_thresh] for
                                                    region_thresh in ['lat_min', 'lat_max',
                                                                      'lon_min', 'lon_max']]

            # Loop through all innovation statistic attributes and
            # update the base-class netCDF object accordingly.
            attrs_dict = {'lat_min': lat_min, 'lat_max': lat_max, 'lon_min': lon_min,
                          'lon_max': lon_max}
            for stat in self.stats_type_list:
                ncvar_dict = {'{0}_{1}'.format(stat, region):
                              {'varname': '{0}_{1}'.format(stat, region),
                               'dims': 'plevs', 'type': 'float64', 'values': eval(stat),
                               'attrs': attrs_dict}}
                self.ncvar_obj = self.update_ncvar(ncvar_dict=ncvar_dict,
                                                   ncvar_obj=self.ncvar_obj)

    def regional_obs(self, ncfilename, innovinfo_obj, vardict, region):
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

        """

        # Determine the observation locations that are valid for the
        # specified region of interest.
        region_info_list = ['lat_min', 'lat_max', 'lon_min', 'lon_max']
        nclat = netcdf4_interface.ncreadvar(
            ncfile=ncfilename, ncvarname=innovinfo_obj.nclat)
        nclon = netcdf4_interface.ncreadvar(
            ncfile=ncfilename, ncvarname=innovinfo_obj.nclon)
        region_info_obj = parser_interface.object_define()
        region_info_dict = parser_interface.object_getattr(
            object_in=self.regions_obj, key=region)
        for item in region_info_list:
            value = parser_interface.dict_key_value(
                dict_in=region_info_dict, key=item, no_split=True)
            region_info_obj = parser_interface.object_setattr(
                object_in=region_info_obj, key=item, value=value)
        msg = ('Computing innovation statistics for region {0}.'.
               format(region))
        self.logger.info(msg=msg)

        # Define all valid observation locations and proceed
        # accordingly.
        obslocs = self.get_obslocs(ncfilename=ncfilename,
                                   innovinfo_obj=innovinfo_obj,
                                   region=region)
        region_cond = numpy.logical_and(numpy.logical_and(
            nclat >= region_info_obj.lat_min, nclat <= region_info_obj.lat_max),
            nclon >= region_info_obj.lon_min,
            nclon <= region_info_obj.lon_max)
        obslocs = numpy.logical_and(obslocs, region_cond)
        return obslocs

    def usage_obs(self, ncfilename, innovinfo_obj, vardict, obslocs,
                  is_wind=False):
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
                value = parser_interface.dict_key_value(dict_in=vardict,
                                                        key=item, force=True)
                if value is None:
                    value = False
                varinfo_obj = parser_interface.object_setattr(
                    object_in=varinfo_obj, key=item, value=value)
            ncobstype = netcdf4_interface.ncreadvar(
                ncfile=ncfilename, ncvarname=innovinfo_obj.ncobstype)
            ncuse = netcdf4_interface.ncreadvar(ncfile=ncfilename,
                                                ncvarname=innovinfo_obj.ncuse)
            pres = netcdf4_interface.ncreadvar(ncfile=ncfilename,
                                               ncvarname=innovinfo_obj.ncpres)
            usage = ncuse == 1
            usage = numpy.logical_and(usage, numpy.isfinite(pres))

            # Exclude specific variables for the wind innovation
            # statistics computations.
            if is_wind:
                obslocs = numpy.logical_and(ncobstype >= 280, ncobstype <= 282)
                obslocs = numpy.logical_or(obslocs, numpy.logical_or(
                    ncobstype == 220, ncobstype == 221))
                obslocs = numpy.logical_or(obslocs, numpy.logical_and(
                    ncobstype >= 230, ncobstype <= 235))
                obslocs = numpy.logical_and(obslocs, usage)
            if not is_wind:
                obslocs = numpy.logical_and(obslocs, usage)
        except Exception as error:
            msg = ('The determination of valid observations to be used for '
                   'the computation of the innovation statistics metrics '
                   'failed with error {0}. Aborting!!!'.format(error))
            raise InnovStatsError(msg=msg)
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
        self.build_database(column_frmt=self.column_frmt,
                            levels_list=self.levels_obj.layer_mean, column_scale=self.column_scale)
        for diagsvar in self.diagsvars:
            (self.ncdim_obj, self.ncvar_obj) = \
                [parser_interface.object_define() for i in range(2)]
            vardict = parser_interface.object_getattr(
                object_in=self.diagsinfo_obj, key=diagsvar)
            ncfilename = self.get_ncfilename(vardict=vardict, ncdim='nlocs')
            if ncfilename is None:
                pass
            if ncfilename is not None:
                innov_var_obj = self.get_innov_var_bool(
                    vardict=vardict, varname=diagsvar)
                if innov_var_obj.is_spechumid:
                    self.innov_spechumid(
                        ncfilename=ncfilename, vardict=vardict)
                if innov_var_obj.is_temperature:
                    self.innov_temperature(
                        ncfilename=ncfilename, vardict=vardict)
                if innov_var_obj.is_wind:
                    self.innov_wind(ncfilename=ncfilename, vardict=vardict)
                self.ncdim_obj = parser_interface.object_setattr(
                    object_in=self.ncdim_obj, key='plevs', value=self.nlevs)
                dimsdict = {'plevs': {'size': len(self.levels_obj.layer_mean),
                                      'type': 'float32',
                                      'values': self.levels_obj.layer_mean}}
                self.ncvar_obj = self.build_ncdims(
                    dimsdict=dimsdict, ncvar_obj=self.ncvar_obj)
                self.write_ncout(vardict=vardict, variable=diagsvar,
                                 ncdim_obj=self.ncdim_obj, ncvar_obj=self.ncvar_obj)
                self.write_database(vardict=vardict, variable=diagsvar,
                                    column_frmt=self.column_frmt, column_scale=self.column_scale)
