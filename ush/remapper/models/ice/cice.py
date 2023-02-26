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
