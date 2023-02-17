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

from dataclasses import dataclass

import numpy
import warnings
import xarray

from typing import Dict, List, Tuple

from ioapps import netcdf4_interface
from tools import fileio_interface, parser_interface

from remapper.remap import esmf, maskupdate
from remapper.remapio import xarray_interface

from exceptions import RemapperError
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

        self.gridfill_kwargs = {'eps': 1.e-2, 'relax': 0.6, 'itermax':
                                1500, 'initzonal': True, 'cyclic': True}

    def build_cice_cfmetadata(self, grid_obj: object) -> object:
        """
        Description
        -----------

        This method defines the dimension variable attributes for the
        output netCDF formatted file assuming the CF metadata
        convention.

        Parameters
        ----------

        grid_obj: object

            A Python object containing the contents of a xarray
            formatted dataset.

        Returns
        -------

        dims_obj: object

            A Python object containing the dimension variable
            attributes for the output netCDF formatted file (assuming
            the CF metadata convention).

        """
        dims_obj = parser_interface.object_define()
        kwargs = {'object_in': self.dstgrid_obj, 'key': 'mass', 'force':
                  True}
        mass_grid_dict = parser_interface.object_getattr(**kwargs)
        if mass_grid_dict is None:
            msg = ('The destination mass grid geographical location '
                   'attributes could not be determined from the user '
                   'experiment configuration. Aborting!!!')
            raise RemapperError(msg=msg)
        kwargs = {'dict_in': mass_grid_dict, 'key': 'nclat', 'force':
                  True, 'no_split': True}
        nclat = parser_interface.dict_key_value(**kwargs)
        if nclat is None:
            msg = ('The latitude coordinate variable %s could not be '
                   'found within the user-specified output grid file '
                   '%s. Aborting!!!' % (nclat, filename))
            raise RemapperError(msg=msg)
        kwargs = {'object_in': grid_obj, 'key': nclat}
        latm = parser_interface.object_getattr(**kwargs)
        latm = latm.mean(axis=1).values
        kwargs = {'dict_in': mass_grid_dict, 'key': 'nclon', 'force':
                  True, 'no_split': True}
        nclon = parser_interface.dict_key_value(**kwargs)
        if nclon is None:
            msg = ('The longitude coordinate variable %s could not be '
                   'found within the user-specified output grid file '
                   '%s. Aborting!!!' % (nclon, filename))
            raise RemapperError(msg=msg)
        kwargs = {'object_in': grid_obj, 'key': nclon}
        lonm = parser_interface.object_getattr(**kwargs)
        lonm = lonm.mean(axis=0).values
        (ncat, time) = (numpy.zeros([self.nlevs]), numpy.zeros([1]))
        ncat = list(numpy.arange(0, self.nlevs, 1))
        dims_list = ['latm', 'lonm', 'ncat', 'time']
        for item in dims_list:
            kwargs = {'object_in': dims_obj, 'key': item, 'value':
                      eval(item)}
            dims_obj = parser_interface.object_setattr(**kwargs)
        return dims_obj

    def build_cice_landmask(self, srcgrid_obj: object, dstgrid_obj: object,
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
        msg = ('Interpolating the CICE land/sea mask.')
        self.logger.info(msg=msg)
        landmask_obj = parser_interface.object_define()
        kwargs = {'object_in': srcgrid_obj, 'key': 'ncfile', 'force':
                  True}
        ncfile = parser_interface.object_getattr(**kwargs)
        kwargs = {'ncfile': ncfile, 'ncvarname': 'wet'}
        ncvar = netcdf4_interface.ncreadvar(**kwargs)
        invar = numpy.where(ncvar > 0.0, 1.0, 0.0)
        kwargs = {'object_in': dstgrid_obj, 'key': 'ncfile', 'force':
                  True}
        ncfile = parser_interface.object_getattr(**kwargs)
        kwargs = {'ncfile': ncfile, 'ncvarname': 'wet'}
        dstgrid_mask = netcdf4_interface.ncreadvar(**kwargs)
        dstgrid_mask = numpy.where(dstgrid_mask > 0.0, 1.0, 0.0)
        kwargs = {'invar': invar, 'interp_type': 'bilinear',
                  'srcgrid_stagger': 'mass', 'dstgrid_stagger': 'mass',
                  'srcgrid_obj': srcgrid_obj, 'dstgrid_obj':
                  dstgrid_obj, 'remap_obj': remap_obj}
        interp_mask = self.remap_app(**kwargs)
        interp_mask = numpy.where(interp_mask < 0.99, 0.0, 1.0)
        (ni, nj, nlevs) = (len(interp_mask[0, :]), len(interp_mask[:, 0]),
                           nlevs)
        landmask_list = ['dstgrid_mask', 'interp_mask']
        for item in landmask_list:
            landmask = numpy.zeros([nlevs, nj, ni])
            kwargs = {'object_in': landmask_obj, 'key': item,
                      'value': landmask}
            landmask_obj = parser_interface.object_setattr(**kwargs)
        landmask_obj.interp_mask[:, ...] = interp_mask[:, :]
        landmask_obj.dstgrid_mask[:, ...] = dstgrid_mask[:, :]
        grid_dims_list = ['ni', 'nj', 'nlevs']
        for grid_dim in grid_dims_list:
            kwargs = {'object_in': landmask_obj, 'key': grid_dim,
                      'value': eval(grid_dim)}
            landmask_obj = parser_interface.object_setattr(**kwargs)
        return landmask_obj

    def build_cice_ncoutput(self, dstgrid_obj: object, landmask_obj: object,
                            varinfo_obj: object, nlevs: int, variable_list: List,
                            output_netcdf: str) -> None:
        """
        Description
        -----------

        This method builds the external netCDF formatted file to
        contain the interpolated input variables.

        Parameters
        ----------

        dstgrid_obj: object

            A Python object containing the destination grid attributes
            collected from the user experiment configuration.

        landmask_obj: object

            A Python object containing the interpolated source grid
            landmask, the destination grid landmask, and the
            destination grid dimensions.

        varinfo_obj: object

            A Python object containing the variable attributes for the
            variables to be interpolated (regridded) from a source
            grid projection to a destination grid projection.

        nlevs: int

            A Python integer value specifying the total number of
            unstaggered vertical levels for the CICE destination grid.

        variable_list: list

            A Python list containing the list of user-specified
            variables to be interpolated (from the user experiment
            configuration).

        output_netcdf: str

            A Python string specifying the path to the external netCDF
            formatted file to contain the CICE initial conditions.

        """
        msg = ('Preparing output file %s for interpolated CICE variables.' %
               output_netcdf)
        self.logger.info(msg=msg)
        kwargs = {'object_in': dstgrid_obj, 'key': 'ncfile', 'force':
                  True}
        ncfile = parser_interface.object_getattr(**kwargs)
        if ncfile is None:
            msg = ('The destination grid netCDF formatted file path '
                   'could not be determined from the user experiment '
                   'configuration. Aborting!!!')
            raise RemapperError(msg=msg)
        kwargs = {'ncfile': ncfile}
        grid_obj = xarray_interface.open(**kwargs)
        kwargs = {'grid_obj': grid_obj}
        dims_obj = self.build_cice_cfmetadata(**kwargs)
        grid_obj.close()
        varobj_list = list()
        for variable in variable_list:
            (nx, ny) = (landmask_obj.ni, landmask_obj.nj)
            kwargs = {'object_in': varinfo_obj, 'key': variable,
                      'force': True}
            varinfo_dict = parser_interface.object_getattr(**kwargs)
            if varinfo_dict is None:
                msg = ('The user experiment configuration does not specify '
                       'the attributes for variable %s and/or could not be '
                       'determined from the user experiment configuration. '
                       'Aborting!!!' % variable)
                raise RemapperError(msg=msg)
            if varinfo_dict['grid_stagger'].lower() == 'mass':
                coords = {'nj': (['nj'], dims_obj.latm),
                          'ni': (['ni'], dims_obj.lonm), 'Time': (['Time'], dims_obj.time)}
                try:
                    coords.update({varinfo_dict['zdim_name']:
                                   (varinfo_dict['zdim_name'], dims_obj.ncat)})
                    ncat = dims_obj.ncat
                    dims = ['Time', 'ncat', 'nj', 'ni']
                    msg = ('Build netCDF variable %s of (x,y,z) dimension (%s,%s,%s).' %
                           (variable, nx, ny, nlevs))
                    varval = numpy.zeros([1, nlevs, ny, nx])
                except KeyError:
                    varinfo_dict['zdim_name'] = None
                    dims = ['Time', 'nj', 'ni']
                    msg = ('Build netCDF variable %s of (x,y) dimension (%s,%s).' %
                           (variable, nx, ny))
                    varval = numpy.zeros([1, ny, nx])
                self.logger.info(msg=msg)
                kwargs = {'dict_in': varinfo_dict,
                          'key': 'dst_ncvarname', 'no_split': True}
                ncvarname = parser_interface.dict_key_value(**kwargs)
                kwargs = {'varval': varval, 'coords': coords, 'dims': dims, 'ncvarname':
                          ncvarname}
                varobj = xarray_interface.varobj(**kwargs)
                varobj_list.append(varobj)
        kwargs = {'ncfile': output_netcdf, 'varobj_list': varobj_list,
                  'unlimitdim': 'Time'}
        xarray_interface.dataset(**kwargs)

    def build_glorys_landmask(self, srcgrid_obj: object, dstgrid_obj: object,
                              remap_obj: object, nlevs: int, fillvalue: float,
                              ncvarname: str) -> object:
        """
        Description
        -----------

        This method identifies the GLORYS landmask and interpolates
        the respective (source) grid landmask to the destination grid
        projection; the destination grid landmask is also defined;
        this information is used to define a Python object containing
        the respective interpolated source grid landmask, the
        destination grid landmask, and the respective dimensions for
        the destination grid.

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
            unstaggered vertical levels for the MOM6 destination grid.

        fillvalue: float

            A Python float valued variable specifying the missing data
            value to be used to construct the GLORYS landmask.

        ncvarname: str

            A Python string specifying the netCDF variable name for
            which the fillvalue was determined/computed.

        Returns
        -------

        landmask_obj: object

            A Python object containing the interpolated source grid
            landmask, the destination grid landmask, and the
            destination grid dimensions.

        """
        landmask_obj = parser_interface.object_define()
        kwargs = {'object_in': srcgrid_obj,
                  'key': 'ncfile', 'force': True}
        ncfile = parser_interface.object_getattr(**kwargs)
        if ncfile is None:
            msg = ('The grid containing the source grid attributes could not be '
                   'determined from the user experiment configuration. Aborting!!!')
            raise RemapperError(msg=msg)
        kwargs = {'ncfile': ncfile, 'ncvarname': ncvarname}
        srcgrid_mask_obj = xarray_interface.read(**kwargs)
        msg = ('Resetting/updating values within source grid land mask.')
        self.logger.info(msg=msg)
        srcgrid_mask = xarray.where(srcgrid_mask_obj.values[0, :, :] >=
                                    fillvalue, 0.0, 1.0)
        srcgrid_mask = numpy.where(srcgrid_mask > 0.0, 1.0, 0.0)
        kwargs = {'object_in': dstgrid_obj, 'key': 'ncfile', 'force': True}
        ncfile = parser_interface.object_getattr(**kwargs)
        if ncfile is None:
            msg = ('The grid containing the destination grid attributes could not be '
                   'determined from the user experiment configuration. Aborting!!!')
            raise RemapperError(msg=msg)
        kwargs = {'ncfile': ncfile, 'ncvarname': 'wet'}
        dstgrid_mask_obj = xarray_interface.read(**kwargs)
        msg = ('Resetting/updating values within destination grid land mask.')
        self.logger.info(msg=msg)
        dstgrid_mask = dstgrid_mask_obj.values
        dstgrid_mask = numpy.where(dstgrid_mask > 0.0, 1.0, 0.0)
        msg = ('Interpolating the source grid land mask to the destination grid '
               'projection.')
        self.logger.info(msg=msg)
        kwargs = {'invar': srcgrid_mask, 'interp_type': 'bilinear',
                  'srcgrid_stagger': 'mass', 'dstgrid_stagger': 'mass',
                  'srcgrid_obj': srcgrid_obj, 'dstgrid_obj':
                  dstgrid_obj, 'remap_obj': remap_obj}
        interp_mask = self.remap_app(**kwargs)
        interp_mask = numpy.where(interp_mask < 0.99, 0.0, 1.0)
        (ni, nj) = (len(interp_mask[0, :]), len(interp_mask[:, 0]))
        landmask_list = ['dstgrid_mask', 'interp_mask']
        for item in landmask_list:
            landmask = numpy.zeros([nj, ni])
            kwargs = {'object_in': landmask_obj,
                      'key': item, 'value': landmask}
            landmask_obj = parser_interface.object_setattr(**kwargs)
        landmask_obj.dstgrid_mask[:, :] = dstgrid_mask[:, :]
        landmask_obj.interp_mask[:, :] = interp_mask[:, :]
        grid_dim_list = ['ni', 'nj']
        for grid_dim in grid_dim_list:
            kwargs = {'object_in': landmask_obj, 'key': grid_dim, 'value':
                      eval(grid_dim)}
            landmask_obj = parser_interface.object_setattr(**kwargs)
        return landmask_obj

    def build_socacice_ncoutput(self, dstgrid_obj: object, landmask_obj: object,
                                varinfo_obj: object, variable_list: List, output_netcdf: str) -> None:
        """
        Description
        -----------

        This method builds the external netCDF formatted file to
        contain the interpolated input variables for a respective SOCA
        sea-ice analysis.

        Parameters
        ----------

        dstgrid_obj: object

            A Python object containing the destination grid attributes
            collected from the user experiment configuration.

        landmask_obj: object

            A Python object containing the interpolated source grid
            landmask, the destination grid landmask, and the
            destination grid dimensions.

        varinfo_obj: object

            A Python object containing the variable attributes for the
            variables to be interpolated (regridded) from a source
            grid projection to a destination grid projection.

        variable_list: list

            A Python list containing the list of user-specified
            variables to be interpolated (from the user experiment
            configuration).

        output_netcdf: str

            A Python string specifying the path to the external netCDF
            formatted file to contain the MOM6 initial conditions.

        """
        msg = ('Preparing output file %s for interpolated SOCA CICE '
               'variables.' % output_netcdf)
        self.logger.info(msg=msg)
        kwargs = {'object_in': dstgrid_obj, 'key': 'ncfile', 'force':
                  True}
        ncfile = parser_interface.object_getattr(**kwargs)
        if ncfile is None:
            msg = ('The destination grid netCDF formatted file path '
                   'could not be determined from the user experiment '
                   'configuration. Aborting!!!')
            raise RemapperError(msg=msg)
        kwargs = {'ncfile': ncfile}
        grid_obj = xarray_interface.open(**kwargs)
        kwargs = {'grid_obj': grid_obj}
        dims_obj = self.build_cice_cfmetadata(**kwargs)
        grid_obj.close()
        varobj_list = list()
        for variable in variable_list:
            (nx, ny) = (landmask_obj.ni, landmask_obj.nj)
            kwargs = {'object_in': varinfo_obj, 'key': variable,
                      'force': True}
            varinfo_dict = parser_interface.object_getattr(**kwargs)
            if varinfo_dict is None:
                msg = ('The user experiment configuration does not specify '
                       'the attributes for variable %s and/or could not be '
                       'determined from the user experiment configuration. '
                       'Aborting!!!' % variable)
                raise RemapperError(msg=msg)
            if varinfo_dict['grid_stagger'].lower() == 'mass':
                coords = {'y_axis': (['y_axis'], dims_obj.latm),
                          'x_axis': (['x_axis'], dims_obj.lonm), 'Time': (['Time'],
                                                                          dims_obj.time)}
                dims = ['Time', 'y_axis', 'x_axis']
                msg = ('Build netCDF variable %s of (x,y) dimension (%s,%s).' %
                       (variable, nx, ny))
                varval = numpy.zeros([1, ny, nx])
                self.logger.info(msg=msg)
                kwargs = {'dict_in': varinfo_dict,
                          'key': 'dst_ncvarname', 'no_split': True}
                ncvarname = parser_interface.dict_key_value(**kwargs)
                kwargs = {'varval': varval, 'coords': coords, 'dims': dims, 'ncvarname':
                          ncvarname}
                varobj = xarray_interface.varobj(**kwargs)
                varobj_list.append(varobj)
        kwargs = {'ncfile': output_netcdf, 'varobj_list': varobj_list,
                  'unlimitdim': 'Time'}
        xarray_interface.dataset(**kwargs)

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
        var_obj = parser_interface.object_define()
        kwargs = {'object_in': var_obj, 'key': 'ncvarname', 'value':
                  ncvarname}
        var_obj = parser_interface.object_setattr(**kwargs)
        return var_obj

    def esmf_remap(self, invar: numpy.array, interp_type: str,
                   srcgrid_stagger: str,
                   dstgrid_stagger: str, srcgrid_obj: object, dstgrid_obj: object,
                   remap_obj: object, reuse_weights: bool = True, src2src: bool = False,
                   dst2dst: bool = False) -> numpy.array:
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
        kwargs = {'srcgrid_obj': srcgrid_obj, 'dstgrid_obj':
                  dstgrid_obj, 'remap_obj': remap_obj}
        esmf_remapper = esmf.ESMF(**kwargs)
        kwargs = {'interp_type': interp_type, 'srcgrid_stagger':
                  srcgrid_stagger, 'dstgrid_stagger': dstgrid_stagger,
                  'invar': invar, 'reuse_weights': reuse_weights,
                  'src2src': src2src, 'dst2dst': dst2dst}
        outvar = esmf_remapper.run(**kwargs)
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

        """
        if interp_scheme == 'esmf':
            remap_app = self.esmf_remap
        if interp_scheme == 'slint':
            msg = ('The SLINT interpolation scheme is not yet supported. '
                   'Aborting!!!')
            raise RemapperError(msg=msg)
        return remap_app

    def maskland(self, landmask_obj: object, output_netcdf: str):
        """
        Description
        -----------

        This method masks interpolated values which occur overland and
        subsequently updates the respective variable arrays within the
        output netCDF formatted file.

        Parameters
        ----------

        landmask_obj: object

            A Python object containing the interpolated source grid
            landmask, the destination grid landmask, and the
            destination grid dimensions.

        output_netcdf: str

            A Python string specifying the path to the external netCDF
            formatted file to contain the MOM6 initial conditions.

        """
        ncvarname_list = ['h', 'Salt', 'Temp', 'u', 'v']
        for ncvarname in ncvarname_list:
            kwargs = {'ncfile': output_netcdf, 'ncvarname': ncvarname,
                      'squeeze': True, 'axis': 0}
            invar = netcdf4_interface.ncreadvar(**kwargs)
            nlevs = invar.shape[0] - 1
            outvar = invar
            outvar = numpy.where(landmask_obj.dstgrid_mask == 0.0,
                                 0.0, invar)
            outvar[0, -1] = outvar[0, -2]
            kwargs = {'ncvarname': ncvarname}
            var_obj = self.build_varobj(**kwargs)
            kwargs = {'ncfile': output_netcdf, 'var_obj': var_obj,
                      'var_arr': outvar}
            xarray_interface.write(**kwargs)
        kwargs = {'ncfile': output_netcdf, 'ncvarname': 'h', 'squeeze':
                  True, 'axis': 0}
        invar = netcdf4_interface.ncreadvar(**kwargs)
        outvar = invar
        outvar = numpy.where(landmask_obj.dstgrid_mask == 0.0,
                             0.0, invar)
        kwargs = {'ncvarname': 'h'}
        var_obj = self.build_varobj(**kwargs)
        kwargs = {'ncfile': output_netcdf, 'var_obj': var_obj,
                  'var_arr': outvar}
        xarray_interface.write(**kwargs)
