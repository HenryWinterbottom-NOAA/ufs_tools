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
