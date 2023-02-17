# =========================================================================

# Module: ush/remapper/models/ocean/__init__.py

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


"""

# ----

import numpy

from typing import Dict, List, Tuple

from exceptions import RemapperError

from remapper.models import Models

from ioapps import netcdf4_interface
from tools import fileio_interface, parser_interface

from remapper.remapio import xarray_interface

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


class Ocean(Models):
    """


    """

    def __init__(self: Models):
        """
        Description
        -----------

        Creates a new Ocean object.

        """

        # Define the base-class attributes.
        super().__init__()

        self.bathy_adjust_default_dict = {'angstrom_z': 1.0e-10, 'htolerance':
                                          1.0e-10}

    def bathy_adjust(self: Models, bathy_obj: object, landmask_obj: object,
                     output_netcdf: str) -> None:
        """
        Description
        -----------

        This method adjusts the destination grid in accordance with
        the destination grid bathymetry and writes both the updated
        sea-surface height (sfc) and thickness profile (h) to the
        output netCDF file (base-class attribute output_netcdf); this
        method follows from the MOM6 MOM_state_initialization module
        adjust adjustetatofitbathymetry subroutine.

        Parameters
        ----------

        bathy_obj: object

            A Python object containing the MOM6 destination grid
            bathymetry attributes.

        landmask_obj: object

            A Python object containing the interpolated source grid
            landmask, the destination grid landmask, and the
            destination grid dimensions.

        output_netcdf: str

            A Python string specifying the path to the external netCDF
            formatted file to contain the MOM6 initial conditions.

        """

        # Update the sea-surface height variable using the land/sea
        # mask and update the netCDF-formatted output file.
        msg = ('Applying the land/sea mask to the sea-surface height variable.')
        self.logger.info(msg=msg)

        sfc = netcdf4_interface.ncreadvar(ncfile=output_netcdf,
                                          ncvarname="sfc", squeeze=True, axis=0)
        sfc = numpy.where(landmask_obj.dstgrid_mask[0, :, :] > 0, sfc, 0.0)
        var_obj = self.build_varobj(ncvarname="sfc")

        xarray_interface.write(ncfile=output_netcdf,
                               var_obj=var_obj, var_arr=sfc)

        # Adjust the interpolated variable thickness profile with
        # respect to the destination grid bathymetry.
        msg = ('Adjusting the thickness profile based on the destination grid '
               'bathymetry.')
        self.logger.info(msg=msg)

        h = netcdf4_interface.ncreadvar(ncfile=output_netcdf, ncvarname="h",
                                        squeeze=True, axis=0)
        eta = numpy.zeros([(landmask_obj.nlevs + 1), landmask_obj.nj,
                           landmask_obj.ni])
        eta[0, :, :] = sfc[:, :]

        # The following code block is adopted from the MOM6
        # MOM_state_initialization module adjust
        # adjustetatofitbathymetry subroutine
        for level in range(0, landmask_obj.nlevs):
            eta[(level+1), :, :] = eta[level, :, :] - h[level, :, :]
        dz = eta[-1, ...] + bathy_obj.depth[...]

        for level in range(0, (landmask_obj.nlevs + 1)):
            eta[level, :, :] = numpy.where(-1.0*eta[level, :, :] >
                                           (bathy_obj.depth[:, :] +
                                            bathy_obj.htolerance),
                                           -1.0*bathy_obj.depth[:, :],
                                           eta[level, :, :])

        for level in range((landmask_obj.nlevs - 1), 0, -1):
            h[level, :, :] = numpy.where(eta[level, :, :] < (eta[(level + 1), :, :] +
                                                             bathy_obj.angstrom_z),
                                         bathy_obj.angstrom_z, (eta[level, :, :] -
                                                                eta[(level + 1), :, :]))

        for xidx in range(landmask_obj.ni):
            for yidx in range(landmask_obj.nj):
                if (-eta[landmask_obj.nlevs, yidx, xidx] <
                        (bathy_obj.depth[yidx, xidx] - bathy_obj.htolerance)):
                    if eta[0, yidx, xidx] <= eta[landmask_obj.nlevs, yidx, xidx]:
                        h[:, yidx, xidx] = (eta[1, yidx, xidx] + bathy_obj.depth[yidx, xidx]) \
                            / (float(landmask_obj.nlevs))
                    else:
                        dilate = (eta[0, yidx, xidx] + bathy_obj.depth[yidx, xidx]) \
                            / (eta[0, yidx, xidx] - eta[landmask_obj.nlevs, yidx, xidx])
                        h[:, yidx, xidx] = h[:, yidx, xidx]*dilate

        # Apply the land mask to the updated/adjusted thickness
        # profile and update the netCDF-formatted output file.
        msg = ('Applying land mask to updated thickness profile.')
        self.logger.info(msg=msg)
        h = numpy.where(landmask_obj.dstgrid_mask > 0, h, 0.0)

        msg = ('Writing updated thickness profile to netCDF file %s.' %
               self.output_netcdf)
        self.logger.info(msg=msg)
        var_obj = self.build_varobj(ncvarname="h")

        xarray_interface.write(ncfile=output_netcdf,
                               var_obj=var_obj, var_arr=h)

    def build_bathy(self: Models, varinfo_obj: object) -> object:
        """
        Description
        -----------

        This method builds the base-class attribute bathy object
        (bathy); this object contains the MOM6 contains the bathymetry
        for the output file which has been editted/updated if a
        bathymetry edits file (bathy_edits_file) is available to the
        MOM6 remapper application.

        Parameters
        ----------

        varinfo_obj: object

            A Python object containing the variable attributes for the
            variables to be interpolated (regridded) from a source
            grid projection to a destination grid projection.

        Returns
        -------

        bathy_obj: object

            A Python object containing the MOM6 destination grid
            bathymetry attributes.

        """
        bathy_obj = parser_interface.object_define()
        kwargs = {'object_in': varinfo_obj, 'key': 'bathy_file', 'force':
                  True}
        bathy_file = parser_interface.object_getattr(**kwargs)
        if bathy_file is None:
            msg = ('The netCDF formatted file containing the bathymetry grid '
                   '(bathy_file) could not be determined from the user experiment '
                   'configuration. Aborting!!!')
            raise RemapperError(msg=msg)
        msg = ('Reading bathymetry from file %s.' % bathy_file)
        self.logger.info(msg=msg)
        kwargs = {'ncfile': bathy_file, 'ncvarname': 'depth'}
        depth_obj = xarray_interface.read(**kwargs)
        depth = depth_obj.values
        kwargs = {'object_in': varinfo_obj, 'key': 'bathy_edits_file', 'force':
                  True}
        bathy_edits_file = parser_interface.object_getattr(**kwargs)
        if bathy_edits_file is None:
            msg = ('The user experiment configuration file does not specify a '
                   'netCDF file path containing the bathymetry grid updates (i.e., '
                   'edits; bathy_edits_file); no bathymetry edits will be made; this may'
                   'cause the application to not perform as expected.')
            self.logger.warn(msg=msg)
        if bathy_edits_file is not None:
            msg = ('Reading bathymetry edits from file %s.' % bathy_edits_file)
            self.logger.info(msg=msg)
            bathy_edits_obj = parser_interface.object_define()
            bathy_edits_attr_list = ['iEdit', 'jEdit', 'nEdits', 'zEdit']
            for bathy_edits_attr in bathy_edits_attr_list:
                kwargs = {'ncfile': bathy_edits_file,
                          'ncvarname': bathy_edits_attr}
                ncvar_obj = xarray_interface.read(**kwargs)
                kwargs = {'object_in': bathy_edits_obj, 'key': bathy_edits_attr,
                          'value': ncvar_obj}
                bathy_edits_obj = parser_interface.object_setattr(
                    **kwargs)
            for idx in bathy_edits_obj.nEdits.values:
                xcoord = bathy_edits_obj.iEdit[idx].values
                ycoord = bathy_edits_obj.jEdit[idx].values
                input_bathy = depth[ycoord, xcoord]
                depth[ycoord, xcoord] = bathy_edits_obj.zEdit[idx].values
                output_bathy = depth[ycoord, xcoord]
                msg = ('Resetting input bathymetry coordinate (%s,%s) value of %s '
                       'to a value of %s.' % (ycoord, xcoord, input_bathy,
                                              output_bathy))
                self.logger.warn(msg=msg)
        depth[:, :] = numpy.where(depth < 9.5, 9.5, depth)
        depth[:, :] = numpy.where(depth > 6500.0, 6500., depth)
        kwargs = {'object_in': bathy_obj, 'key': 'depth', 'value': depth}
        bathy_obj = parser_interface.object_setattr(**kwargs)
        for key in self.bathy_adjust_default_dict.keys():
            kwargs = {'dict_in': self.bathy_adjust_default_dict, 'key': key, 'no_split':
                      True}
            value = parser_interface.dict_key_value(**kwargs)
            kwargs = {'object_in': bathy_obj, 'key': key, 'value': value}
            bathy_obj = parser_interface.object_setattr(**kwargs)
        return bathy_obj

    def build_mom6_cfmetadata(self, grid_obj: object, dstgrid_obj: object,
                              nlevs: int) -> object:
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

        dstgrid_obj: object

            A Python object containing the destination grid attributes
            collected from the user experiment configuration.

        nlevs: int

            A Python integer value specifying the total number of
            unstaggered vertical levels for the MOM6 destination grid.

        Returns
        -------

        dims_obj: object

            A Python object containing the dimension variable
            attributes for the output netCDF formatted file (assuming
            the CF metadata convention).

        """
        dims_obj = parser_interface.object_define()
        kwargs = {'object_in': dstgrid_obj, 'key': 'mass', 'force':
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
        kwargs = {'object_in': self.dstgrid_obj, 'key': 'uvel', 'force':
                  True}
        uvel_grid_dict = parser_interface.object_getattr(**kwargs)
        if uvel_grid_dict is None:
            msg = ('The destination zonal velocity grid geographical location '
                   'attributes could not be determined from the user '
                   'experiment configuration. Aborting!!!')
            raise RemapperError(msg=msg)
        kwargs = {'dict_in': uvel_grid_dict, 'key': 'nclat', 'force':
                  True, 'no_split': True}
        nclat = parser_interface.dict_key_value(**kwargs)
        if nclat is None:
            msg = ('The latitude coordinate variable %s could not be '
                   'found within the user-specified output grid file '
                   '%s. Aborting!!!' % (nclat, filename))
            raise RemapperError(msg=msg)
        kwargs = {'object_in': grid_obj, 'key': nclat}
        latu = parser_interface.object_getattr(**kwargs)
        latu = latu.mean(axis=1).values
        kwargs = {'dict_in': uvel_grid_dict, 'key': 'nclon', 'force':
                  True, 'no_split': True}
        nclon = parser_interface.dict_key_value(**kwargs)
        if nclon is None:
            msg = ('The longitude coordinate variable %s could not be '
                   'found within the user-specified output grid file '
                   '%s. Aborting!!!' % (nclon, filename))
            raise RemapperError(msg=msg)
        kwargs = {'object_in': grid_obj, 'key': nclon}
        lonu = parser_interface.object_getattr(**kwargs)
        lonu = lonu.mean(axis=0).values
        kwargs = {'object_in': dstgrid_obj, 'key': 'vvel', 'force':
                  True}
        vvel_grid_dict = parser_interface.object_getattr(**kwargs)
        if vvel_grid_dict is None:
            msg = ('The destination meridional velocity grid geographical '
                   'location attributes could not be determined from the user '
                   'experiment configuration. Aborting!!!')
            raise RemapperError(msg=msg)
        kwargs = {'dict_in': vvel_grid_dict, 'key': 'nclat', 'force':
                  True, 'no_split': True}
        nclat = parser_interface.dict_key_value(**kwargs)
        if nclat is None:
            msg = ('The latitude coordinate variable %s could not be '
                   'found within the user-specified output grid file '
                   '%s. Aborting!!!' % (nclat, filename))
            raise RemapperError(msg=msg)
        kwargs = {'object_in': grid_obj, 'key': nclat}
        latv = parser_interface.object_getattr(**kwargs)
        latv = latv.mean(axis=1).values
        kwargs = {'dict_in': vvel_grid_dict, 'key': 'nclon', 'force':
                  True, 'no_split': True}
        nclon = parser_interface.dict_key_value(**kwargs)
        if nclon is None:
            msg = ('The longitude coordinate variable %s could not be '
                   'found within the user-specified output grid file '
                   '%s. Aborting!!!' % (nclon, filename))
            raise RemapperError(msg=msg)
        kwargs = {'object_in': grid_obj, 'key': nclon}
        lonv = parser_interface.object_getattr(**kwargs)
        lonv = lonv.mean(axis=0).values
        (layers, time) = (numpy.zeros([nlevs]), numpy.zeros([1]))
        dims_list = ['latm', 'latu', 'latv', 'layers', 'lonm', 'lonu',
                     'lonv', 'time']
        for item in dims_list:
            kwargs = {'object_in': dims_obj, 'key': item, 'value':
                      eval(item)}
            dims_obj = parser_interface.object_setattr(**kwargs)
        return dims_obj

    def build_mom6_landmask(self, srcgrid_obj: object, dstgrid_obj: object,
                            remap_obj: object, nlevs: int) -> object:
        """
        Description
        -----------

        This method interpolates the source grid landmask to the
        destination grid projection; the destination grid landmask is
        also defined; this information is used to define a Python
        object containing the respective interpolated source grid
        landmask, the destination grid landmask, and the respective
        dimensions for the destination grid.

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

        Returns
        -------

        landmask_obj: object

            A Python object containing the interpolated source grid
            landmask, the destination grid landmask, and the
            destination grid dimensions.

        """
        landmask_obj = parser_interface.object_define()
#        kwargs = {'object_in': srcgrid_obj,
#                  'key': 'ncfile', 'force': True}
#        ncfile = parser_interface.object_getattr(**kwargs)
#        if ncfile is None:
#            msg = ('The grid containing the source grid attributes could not be '
#                   'determined from the user experiment configuration. Aborting!!!')
#            raise RemapperError(msg=msg)
        # kwargs = {'ncfile': ncfile, 'ncvarname': 'wet'}

        # kwargs = {'object_in': srcgrid_obj,
        #          'key': 'mask_ncfile', 'force': True}
        #mask_ncfile = parser_interface.object_getattr(**kwargs)

        # if mask_ncfile is None:
        #    msg = ('The grid containing the source grid attributes could not be '
        #           'determined from the user experiment configuration. Aborting!!!')
        #    raise RemapperError(msg=msg)

        #kwargs = {'ncfile': mask_ncfile, 'ncvarname': 'mask'}
        #srcgrid_mask_obj = xarray_interface.read(**kwargs)
        msg = ('Resetting/updating values within source grid land mask.')
        self.logger.info(msg=msg)
        srcgrid_mask = srcgrid_mask_obj.values
        srcgrid_mask = numpy.where(srcgrid_mask > 0.0, 1.0, 0.0)
        kwargs = {'object_in': dstgrid_obj, 'key': 'ncfile', 'force': True}
        ncfile = parser_interface.object_getattr(**kwargs)
        if ncfile is None:
            msg = ('The grid containing the destination grid attributes could not be '
                   'determined from the user experiment configuration. Aborting!!!')
            raise RemapperError(msg=msg)

        kwargs = {'object_in': dstgrid_obj,
                  'key': 'mask_ncfile', 'force': True}
        mask_ncfile = parser_interface.object_getattr(**kwargs)

        if mask_ncfile is None:
            msg = ('The grid containing the source grid attributes could not be '
                   'determined from the user experiment configuration. Aborting!!!')
            raise RemapperError(msg=msg)

        kwargs = {'ncfile': mask_ncfile, 'ncvarname': 'mask'}
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
        (ni, nj, nlevs) = (len(interp_mask[0, :]), len(
            interp_mask[:, 0]), nlevs)
        landmask_list = ['dstgrid_mask', 'interp_mask']
        for item in landmask_list:
            landmask = numpy.zeros([nlevs, nj, ni])
            kwargs = {'object_in': landmask_obj,
                      'key': item, 'value': landmask}
            landmask_obj = parser_interface.object_setattr(**kwargs)
        landmask_obj.dstgrid_mask[:, ...] = dstgrid_mask[:, :]
        landmask_obj.interp_mask[:, ...] = interp_mask[:, :]
        grid_dim_list = ['ni', 'nj', 'nlevs']
        for grid_dim in grid_dim_list:
            kwargs = {'object_in': landmask_obj, 'key': grid_dim, 'value':
                      eval(grid_dim)}
            landmask_obj = parser_interface.object_setattr(**kwargs)
        return landmask_obj

    def build_mom6_ncoutput(self, dstgrid_obj: object, varinfo_obj: object,
                            nlevs: int, variable_list: List, output_netcdf: str,
                            update_hcoord: bool = False, update_zcoord: bool = False) -> None:
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

        varinfo_obj: object

            A Python object containing the variable attributes for the
            variables to be interpolated (regridded) from a source
            grid projection to a destination grid projection.

        nlevs: int

            A Python integer value specifying the total number of
            unstaggered vertical levels for the MOM6 destination grid.

        variable_list: list

            A Python list containing the list of user-specified
            variables to be interpolated (from the user experiment
            configuration).

        output_netcdf: str

            A Python string specifying the path to the external netCDF
            formatted file to contain the MOM6 initial conditions.

        Keywords
        --------

        update_zcoord: bool, optional

            A Python boolean-type variable specifying whether the
            respective input model vertical coordinate needs to be
            updated to comply with the expectations of MOM6.

        """
        msg = ('Preparing output file %s for interpolated MOM6 variables.' %
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
        kwargs = {'grid_obj': grid_obj, 'dstgrid_obj': dstgrid_obj,
                  'nlevs': nlevs}
        dims_obj = self.build_mom6_cfmetadata(**kwargs)
        grid_obj.close()
        varobj_list = list()
        for variable in variable_list:
            coords = dict()
            kwargs = {'object_in': varinfo_obj, 'key': variable, 'force':
                      True}
            varinfo_dict = parser_interface.object_getattr(**kwargs)
            if varinfo_dict is None:
                msg = ('The user experiment configuration does not specify '
                       'the attributes for variable %s and/or could not be '
                       'determined from the user experiment configuration. '
                       'Aborting!!!' % variable)
                raise RemapperError(msg=msg)
            if varinfo_dict['grid_stagger'].lower() == 'mass':
                if update_hcoord:
                    msg = ('Updating horizontal mass coordinates for variable '
                           '%s.' % variable)
                    self.logger.warn(msg=msg)
                    varinfo_dict['xdim_name'] = 'lonh'
                    varinfo_dict['ydim_name'] = 'lath'
                coords = {varinfo_dict['ydim_name']: ([varinfo_dict['ydim_name']],
                                                      dims_obj.latm),
                          varinfo_dict['xdim_name']: ([varinfo_dict['xdim_name']],
                                                      dims_obj.lonm),
                          'Time': (['Time'], dims_obj.time)}
                (nx, ny) = (len(dims_obj.lonm), len(dims_obj.latm))
            if varinfo_dict['grid_stagger'].lower() == 'uvel':
                if update_hcoord:
                    msg = ('Updating horizontal u-velocity coordinates for variable '
                           '%s.' % variable)
                    self.logger.warn(msg=msg)
                    varinfo_dict['xdim_name'] = 'lonq'
                    varinfo_dict['ydim_name'] = 'lath'
                coords = {varinfo_dict['ydim_name']: ([varinfo_dict['ydim_name']],
                                                      dims_obj.latm),
                          varinfo_dict['xdim_name']: ([varinfo_dict['xdim_name']],
                                                      dims_obj.lonm),
                          'Time': (['Time'], dims_obj.time)}
                (nx, ny) = (len(dims_obj.lonu), len(dims_obj.latu))
            if varinfo_dict['grid_stagger'].lower() == 'vvel':
                if update_hcoord:
                    msg = ('Updating horizontal v-velocity coordinates for variable '
                           '%s.' % variable)
                    self.logger.warn(msg=msg)
                    varinfo_dict['xdim_name'] = 'lonh'
                    varinfo_dict['ydim_name'] = 'latq'
                coords = {varinfo_dict['ydim_name']: ([varinfo_dict['ydim_name']],
                                                      dims_obj.latm),
                          varinfo_dict['xdim_name']: ([varinfo_dict['xdim_name']],
                                                      dims_obj.lonm),
                          'Time': (['Time'], dims_obj.time)}
                (nx, ny) = (len(dims_obj.lonv), len(dims_obj.latv))
            if update_zcoord:
                kwargs = {'dict_in': varinfo_dict, 'key': 'zdim_name', 'force':
                          True, 'no_split': True}
                zcoord = parser_interface.dict_key_value(**kwargs)
                if zcoord is None:
                    msg = ('The variable %s does not have a z-coordinate; '
                           'skipping.' % variable)
                    self.logger.warn(msg=msg)
                if zcoord is not None:
                    msg = ('Updating z-coordinate for variable %s from %s to '
                           'Layer.' % (variable, zcoord))
                    self.logger.warn(msg=msg)
                    varinfo_dict['zdim_name'] = 'Layer'
            try:
                coords.update({varinfo_dict['zdim_name']: ([varinfo_dict['zdim_name']],
                                                           dims_obj.layers)})
                dims = ['Time', varinfo_dict['zdim_name'], varinfo_dict['ydim_name'],
                        varinfo_dict['xdim_name']]
                msg = ('Building netCDF array for variable %s of (x,y,z) dimension '
                       '(%s,%s,%s).' % (variable, nx, ny, nlevs))
            except KeyError:
                varinfo_dict['zdim_name'] = None
                dims = ['Time', varinfo_dict['ydim_name'],
                        varinfo_dict['xdim_name']]
                msg = ('Build netCDF variable %s of (x,y) dimension (%s,%s).' %
                       (variable, nx, ny))
            self.logger.info(msg=msg)
            if varinfo_dict['zdim_name'] is None:
                varval = numpy.zeros([1, ny, nx])
            if varinfo_dict['zdim_name'] is not None:
                varval = numpy.zeros([1, nlevs, ny, nx])
            kwargs = {'dict_in': varinfo_dict,
                      'key': 'dst_ncvarname', 'no_split': True}
            ncvarname = parser_interface.dict_key_value(**kwargs)
            kwargs = {'varval': varval, 'coords': coords, 'dims': dims, 'ncvarname':
                      ncvarname}
            varobj = xarray_interface.varobj(**kwargs)
            varobj_list.append(varobj)
        kwargs = {'ncfile': output_netcdf, 'varobj_list': varobj_list,
                  'unlimitdim': 'Time'}  # , 'variable_list': variable_list}
        xarray_interface.dataset(**kwargs)
        kwargs = {'ncfile': output_netcdf, 'ncvarname': 'lonq',
                  'ncvar': dims_obj.lonu}
        netcdf4_interface.ncwritevar(**kwargs)
        kwargs = {'ncfile': output_netcdf, 'ncvarname': 'latq',
                  'ncvar': dims_obj.latv}
        netcdf4_interface.ncwritevar(**kwargs)

    def build_mom6_varattrobjs(self, varinfo_dict: Dict) -> Tuple[object, object]:
        """
        Description
        -----------

        This method collects the attributes for the MOM6 variables to
        be interpolated (regridded) from the user experiment
        configuration file.

        Parameters
        ----------

        varinfo_dict: dict

            A Python dictionary containing the user experiment
            configuration MOM6 variable attributes.

        Returns
        -------

        ncattr_obj: object

            A Python object containing the netCDF file attributes
            collected from the user experiment configuration.

        varattr_obj: object

            A Python object containing the MOM6 variable attributes
            collected from the user experiment configuration.

        """
        (ncattr_obj, varattr_obj) = [parser_interface.object_define()
                                     for i in range(2)]
        ncattrs_list = ['ncfilename', 'src_ncvarname', 'dst_ncvarname']
        for ncattr in ncattrs_list:
            kwargs = {'dict_in': varinfo_dict, 'key': ncattr, 'force': True,
                      'no_split': True}
            value = parser_interface.dict_key_value(**kwargs)
            if value is None:
                msg = ('The netCDF attribute %s could not be determined from '
                       'the user experiment configuration for variable %s. '
                       'Aborting!!!' % (ncattr, variable))
                raise RemapperError(msg=msg)
            kwargs = {'object_in': ncattr_obj,
                      'key': ncattr, 'value': value}
            ncattr_obj = parser_interface.object_setattr(**kwargs)
            kwargs = {'path': ncattr_obj.ncfilename}
            exist = fileio_interface.fileexist(**kwargs)
            if not exist:
                msg = ('The netCDF file %s does not exist. Aborting!!!' %
                       ncattr_obj.ncfilename)
                raise RemapperError(msg=msg)
        varattr_list = ['grid_stagger',
                        'interp_type', 'xdim_name', 'ydim_name']
        for varattr in varattr_list:
            kwargs = {'dict_in': varinfo_dict, 'key': varattr, 'force': True,
                      'no_split': True}
            value = parser_interface.dict_key_value(**kwargs)
            if value is None:
                msg = ('The variable attribute %s could not be determined for '
                       'variable %s from the user experiment configuration. '
                       'Aborting!!!' % (varattr, variable))
                raise RemapperError(msg=msg)
            kwargs = {'object_in': varattr_obj,
                      'key': varattr, 'value': value}
            varattr_obj = parser_interface.object_setattr(**kwargs)
            kwargs = {'dict_in': varinfo_dict, 'key': 'zdim_name', 'force': True,
                      'no_split': True}
            value = parser_interface.dict_key_value(**kwargs)
            kwargs = {'object_in': varattr_obj,
                      'key': 'zdim_name', 'value': value}
            varattr_obj = parser_interface.object_setattr(**kwargs)
        return (ncattr_obj, varattr_obj)

    def build_oras5_landmask(self, srcgrid_obj: object, dstgrid_obj: object,
                             remap_obj: object, nlevs: int) -> object:
        """
        Description
        -----------

        This method interpolates the source grid landmask to the
        destination grid projection; the destination grid landmask is
        also defined; this information is used to define a Python
        object containing the respective interpolated source grid
        landmask, the destination grid landmask, and the respective
        dimensions for the destination grid.

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
        kwargs = {'ncfile': ncfile, 'ncvarname': 'thetao_oras'}
        srcgrid_mask_obj = xarray_interface.read(**kwargs)
        msg = ('Resetting/updating values within source grid land mask.')
        self.logger.info(msg=msg)
        srcgrid_mask = xarray.where(xarray.ufuncs.isnan(
            srcgrid_mask_obj.values[0, 0, :, :]), 0.0, 1.0)
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
        (ni, nj, nlevs) = (len(interp_mask[0, :]), len(
            interp_mask[:, 0]), nlevs)
        landmask_list = ['dstgrid_mask', 'interp_mask']
        for item in landmask_list:
            landmask = numpy.zeros([nlevs, nj, ni])
            kwargs = {'object_in': landmask_obj,
                      'key': item, 'value': landmask}
            landmask_obj = parser_interface.object_setattr(**kwargs)
        for i in range(nlevs):
            landmask_obj.dstgrid_mask[i, :, :] = dstgrid_mask[:, :]
            landmask_obj.interp_mask[i, :, :] = interp_mask[:, :]
        grid_dim_list = ['ni', 'nj', 'nlevs']
        for grid_dim in grid_dim_list:
            kwargs = {'object_in': landmask_obj, 'key': grid_dim, 'value':
                      eval(grid_dim)}
            landmask_obj = parser_interface.object_setattr(**kwargs)
        return landmask_obj
