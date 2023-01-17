# =========================================================================

# $$$ MODULE DOCUMENTATION BLOCK

# UFS-RNR-containers :: innov_stats/py/innov_stats/soca_ocean.py

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

   SOCAOcean(yaml_file, basedir)

       This is the base-class object for all Sea-ice and Ocean
       Analysis (SOCA) ocean innovation statistics diagnostics; it is
       a sub-class of InnovStats.

Author(s)
---------

   Henry R. Winterbottom; 02 August 2022

History
-------

   2022-08-02: Henry Winterbottom -- Initial implementation.

"""

# ----

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


class SOCAOcean(InnovStats):
    """
    Description
    -----------

    This is the base-class object for all Sea-ice and Ocean Analysis
    (SOCA) ocean innovation statistics diagnostics; it is a sub-class
    of InnovStats.

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

        Creates a new SOCAOcean object.

        """

        # Define the base-class attributes.
        super(SOCAOcean, self).__init__(
            yaml_file=yaml_file, basedir=basedir)
        self.column_frmt = 'meters_%03d'
        self.column_scale = 1000.0
        self.get_depths()

    def get_depths(self):
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

            * raised if the depth levels information cannot be
              determined from the user experiment configuration.

        """

        # Define the depth profile information from the user
        # experiment configuration and proceed accordingly.
        depth_levels_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key='depth_levels',
            force=True)
        if depth_levels_dict is None:
            msg = ('The depth levels information could not be '
                   'determined from the user experiment configuration. '
                   'Aborting!!!')
            raise InnovStatsError(msg=msg)

        # Define the base-class attributes levels_obj; this namespace
        # contains both the layer bounding regions (i.e., top and
        # bottom layer interfaces for the respective level) as well as
        # the mean value for the respective interval (i.e., the layer
        # mean -- the middle of the layer).
        levels_attr_dict = {'depth_bottom': 'layer_bottom', 'depth_top':
                            'layer_top'}
        for levels_attr in levels_attr_dict.keys():
            value = parser_interface.dict_key_value(
                dict_in=depth_levels_dict, key=levels_attr,
                force=True)
            if value is None:
                msg = ('The attribute {0} could not be determined from '
                       'the user experiment configuration. Aborting!!!'.
                       format(levels_attr))
                raise InnovStatsError(msg=msg)
            self.levels_obj = parser_interface.object_setattr(
                object_in=self.levels_obj, key=levels_attr_dict[levels_attr],
                value=value)
        value = [statistics.mean(k) for k in zip(self.levels_obj.layer_bottom,
                                                 self.levels_obj.layer_top)]
        self.levels_obj = parser_interface.object_setattr(
            object_in=self.levels_obj, key='layer_mean', value=value)

        # Define the total number of vertical levels.
        self.nlevs = len(self.levels_obj.layer_mean)

    def regional_stats(self, ncfilename, innovinfo_obj):
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

        # Loop through each region of interest to compute the
        # respective innovation statistics; proceed accordingly.

        for region in vars(self.regions_obj):
            msg = ('Computing innovation statistics for region {0}.'.
                   format(region))
            self.logger.info(msg=msg)
            obslocs = self.get_obslocs(ncfilename=ncfilename,
                                       innovinfo_obj=innovinfo_obj,
                                       region=region)
            qc = netcdf4_interface.ncreadvar(ncfile=ncfilename,
                                             ncvarname=innovinfo_obj.ncqc)
            qcchk = numpy.logical_and(qc[obslocs] == 0, qc[obslocs] <= 0)
            depth = netcdf4_interface.ncreadvar(ncfile=ncfilename,
                                                ncvarname=innovinfo_obj.ncdepth)
            depth = depth[obslocs]
            depth = depth[qcchk]
            omf = netcdf4_interface.ncreadvar(ncfile=ncfilename,
                                              ncvarname=innovinfo_obj.ncomf)
            omf = omf[obslocs]
            omf = omf[qcchk]
            (bias, count, rmsd) = [numpy.empty(self.nlevs) for i in range(3)]
            for i in range(0, self.nlevs):
                depth_obs = numpy.logical_and(depth >= self.levels_obj.layer_top[i],
                                              depth <= self.levels_obj.layer_bottom[i])
                bias[i] = numpy.mean(omf[depth_obs])
                count[i] = numpy.ma.count(omf[depth_obs])
                rmsd[i] = numpy.sqrt(numpy.mean(omf[depth_obs]*omf[depth_obs]))
            region_info_dict = parser_interface.object_getattr(
                object_in=self.regions_obj, key=region)
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
                               'dims': 'depth', 'type': 'float64', 'values': eval(stat),
                               'attrs': attrs_dict}}
                self.update_ncvar(ncvar_dict=ncvar_dict)

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
            innovinfo_obj = self.get_innovinfo(
                vardict=vardict, addinfo_list=['ncdepth', 'ncomf', 'ncqc'])
            ncfilename = self.get_ncfilename(vardict=vardict, ncdim='nlocs')
            if ncfilename is None:
                pass
            if ncfilename is not None:
                self.ncdim_obj = parser_interface.object_setattr(
                    object_in=self.ncdim_obj, key='depth', value=self.nlevs)
                dimsdict = {'depth': {'size': len(self.levels_obj.layer_mean),
                                      'type': 'float32',
                                      'values': self.levels_obj.layer_mean}}
                self.build_ncdims(dimsdict=dimsdict)
                self.regional_stats(ncfilename=ncfilename,
                                    innovinfo_obj=innovinfo_obj)
                self.write_ncout(vardict=vardict, variable=diagsvar,
                                 ncdim_obj=self.ncdim_obj, ncvar_obj=self.ncvar_obj)
                self.write_database(vardict=vardict, variable=diagsvar)
