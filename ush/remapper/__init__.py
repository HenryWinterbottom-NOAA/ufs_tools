from dataclasses import dataclass

import os

from remapper.models.ice.cice import CICE
from remapper.models.ocean.mom6 import MOM6

from confs.yaml_interface import YAML

from typing import Tuple, Union

from exceptions import RemapperError
from utils.logger_interface import Logger

from tools import fileio_interface, parser_interface

# ----


@dataclass
class Remapper:
    """

    """

    def __init__(self, options_obj: object):
        """
        Description
        -----------

        Creates a new Remapper object.

        """

        # Define the base-class attributes.
        self.options_obj = options_obj
        self.logger = Logger()
        self.yaml_file = self.options_obj.yaml_file
        self.yaml_dict = YAML().read_yaml(yaml_file=self.yaml_file)

        # Define the forecast model attributes.
        self.forecast_model = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key="forecast_model", force=True,
            no_split=True)

        if self.forecast_model is None:
            msg = ('The forecast model for which to perform the variable '
                   'remapping could not be determined from the user experiment '
                   'configuration. Aborting!!!')
            raise RemapperError(msg=msg)

        msg = f"Remapping variables for forecast model {self.forecast_model}."
        self.logger.info(msg=msg)

        # Define the available remapping attributes.
        self.interp_schemes = ["esmf"]
        self.grid_types_list = ['mass', 'uvel', 'vvel']
        self.grid_items_list = ['arakawa', 'grid_ncfile', 'topo_ncfile']
        self.coord_types_list = ['nclat', 'nclon']

        # Define the supported forecast model remapping applications.
        self.models_dict = {"cice": CICE,
                            "mom6": MOM6
                            }

        # Define the supported remapping types; these are the YAML
        # keys (and value pairs) corresponding to the types of
        # remapping coefficients contained within the
        # respective/provided files.
        self.remap_types_list = ['dstuvel2dstmass_bilinear',
                                 'dstuvel2dstmass_nrstnghbr',
                                 'dstvvel2dstmass_bilinear',
                                 'dstvvel2dstmass_nrstnghbr',
                                 'dstmass2dstuvel_bilinear',
                                 'dstmass2dstuvel_nrstnghbr',
                                 'dstmass2dstvvel_bilinear',
                                 'dstmass2dstvvel_nrstnghbr',
                                 'srcmass2dstmass_bilinear',
                                 'srcmass2dstmass_nrstnghbr',
                                 'srcmass2dstuvel_bilinear',
                                 'srcmass2dstuvel_nrstnghbr',
                                 'srcmass2dstvvel_bilinear',
                                 'srcmass2dstvvel_nrstnghbr',
                                 'srcuvel2dstmass_bilinear',
                                 'srcuvel2dstmass_nrstnghbr',
                                 'srcvvel2dstmass_bilinear',
                                 'srcvvel2dstmass_nrstnghbr',
                                 'srcuvel2srcmass_bilinear',
                                 'srcuvel2srcmass_nrstnghbr',
                                 'srcvvel2srcmass_bilinear',
                                 'srcvvel2srcmass_nrstnghbr'
                                 ]

    def getgrid_info(self, grid_name: str) -> object:
        """
        Description
        -----------

        This method parses the user experiment configuration and
        collects the input and output grid configuration attributes
        and defines Python objects for the respective source and
        destination grids; the returned object attributes specify the
        netCDF variable names for the respective geographical
        coordinate variables (e.g., latitude and longitude) in
        accordance with the respective Arakawa grid staggerings as
        well as the netCDF formatted file paths for the respective
        source and desination grids.

        Parameters
        ----------

        grid_name: str

            A Python string specifying the grid (e.g., source or
            destination) to be collected from the user experiment
            configuration.

        Returns
        -------

        grid_obj: object

            A Python object containing the respective grid (e.g.,
            source or destination) attributes.

        """

        # Collect the source and destination grid attributes.
        matching = [s for s in ['source', 'destination'] if grid_name in s]
        if not matching:
            msg = (f'The specified grid type {grid_name} is not supported or cannot '
                   'be determined from the user experiment configuration. '
                   'Aborting!!!')
            raise RemapperError(msg=msg)

        # Collect the respective grid attributes and define the Python
        # object accordingly.
        grid_obj = parser_interface.object_define()
        kwargs = {'dict_in': self.yaml_dict, 'key': 'grids', 'force': True}
        grids_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key="grids", force=True)
        if grids_dict is None:
            msg = ('The grids attribute could not be determined from the '
                   f'user experiment configuration in YAML file {self.yaml_file}. '
                   'Aborting!!!')
            raise RemapperError(msg=msg)

        grid_dict = parser_interface.dict_key_value(
            dict_in=grids_dict, key=grid_name, force=True)

        if grid_dict is None:
            msg = (f'The grid type {grid_name} could not be determined from the user '
                   'experiment grid attributes. Aborting!!!')
            raise RemapperError(msg=msg)

        for grid_type in self.grid_types_list:
            grid_type_dict = parser_interface.dict_key_value(
                dict_in=grid_dict, key=grid_type, force=True)
            if grid_type_dict is None:
                msg = (f'The attributes for grid type {grid_type} could not be '
                       'determined from the user experiment configuration. '
                       'Aborting!!!')
                raise RemapperError(msg=msg)

            # Define the coordinate attributes.
            coord_dict = {}

            for coord_type in self.coord_types_list:
                coord_name = parser_interface.dict_key_value(
                    dict_in=grid_type_dict, key=coord_type, force=True, no_split=True)
                if coord_name is None:
                    msg = (f'The netCDF variable name for coordinate type {coord_type} for '
                           f'grid type {grid_type} could not be determined from the user '
                           'experiment configuration. Aborting!!!')
                    raise RemapperError(msg=msg)

                coord_dict[coord_type] = coord_name
            grid_obj = parser_interface.object_setattr(
                object_in=grid_obj, key=grid_type, value=coord_dict)

        for grid_item in self.grid_items_list:
            value = parser_interface.dict_key_value(
                dict_in=grid_dict, key=grid_item, force=True, no_split=True)
            if value is None:
                msg = (f'The {grid_name} grid attribute {grid_item} could not be '
                       'determined from the user experiment configuration. Aborting!!!')
                raise RemapperError(msg=msg)

            grid_obj = parser_interface.object_setattr(
                object_in=grid_obj, key=grid_item, value=value)

            if grid_name.lower() == "destination":
                output_netcdf = parser_interface.dict_key_value(
                    dict_in=self.yaml_dict, key="output_netCDF", force=True,
                    no_split=True)

                if output_netcdf is None:
                    msg = (f'The experiment configuration file {self.yaml_file} does not '
                           'contain the attribute output_netCDF. Aborting!!!')
                    raise RemapperError(msg=msg)

                grid_obj = parser_interface.object_setattr(
                    object_in=grid_obj, key="output_netcdf", value=output_netcdf)

                fileio_interface.dirpath_tree(
                    path=os.path.dirname(output_netcdf))

        return grid_obj

    def getremap_info(self):
        """
        Description
        -----------

        This method collects the remapping attributes within the
        experiment configuration; the returned remap_obj contains the
        experiment configuration specified remapping type and
        associated source to destination grid remapping attributes
        netCDF file paths as well as the interpolation scheme; if the
        user experiment configuration does not contain a netCDF file
        for the respective remapping type, the remap_obj attribute is
        defined as NoneType.

        Returns
        -------

        remap_obj: obj

            A Python object containing the user experiment
            configuration netCDF file paths for the respective
            supported remapping types.

        """

        # Collect the remapping attributes from the experiment
        # configuration file.
        remap_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key="remap", force=True)
        if remap_dict is None:
            msg = ('The remapping information and/or the remap attribute '
                   'could not be determined from the user experiment '
                   'configuration. Aborting!!!')
            raise RemapperError(msg=msg)

        # Define a Python object containing each of the remapping
        # types (i.e., netCDF file paths) specified within the
        # experiment configuration.
        remap_obj = parser_interface.object_define()

        for remap_type in self.remap_types_list:
            remap_file = parser_interface.dict_key_value(
                dict_in=remap_dict, key=remap_type, force=True, no_split=True)
            if remap_file is None:
                msg = ('The user experiment configuration does not specify a '
                       f'netCDF file path for remap type {remap_type}; the respective '
                       'remapping type will not be supported.')
                self.logger.warn(msg=msg)

            if remap_file is not None:
                exist = fileio_interface.fileexist(path=remap_file)
                if not exist:
                    msg = f'The file path {remap_file} does not exist. Aborting!!!'
                    raise RemapperError(msg=msg)

                msg = (f'The remapping type {remap_type} will use the attributes '
                       f'within netCDF file path {remap_file}.')
                self.logger.info(msg=msg)

            remap_obj = parser_interface.object_setattr(
                object_in=remap_obj, key=remap_type, value=remap_file)

        # Define rotation of current options.
        rotate_currents = parser_interface.dict_key_value(
            dict_in=remap_dict, key="rotate_currents", force=True)
        if rotate_currents is None:
            rotate_currents = False

        remap_obj = parser_interface.object_setattr(
            object_in=remap_obj, key="rotate_currents", value=rotate_currents)

        # Determine and define the interpolation scheme attributes.
        interp_scheme = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key="interp_scheme", force=True, no_split=True)
        if interp_scheme is None:
            msg = ("The interpolation scheme as not been define within configuration "
                   f"file {self.yaml_file}; the supported scheme(s) is(are) "
                   f"{self.interp_schemes}. Aborting!!!")

        if not [s for s in self.interp_schemes if interp_scheme in s]:
            msg = (f'The specified interpolation scheme {interp_scheme} is currently not '
                   'supported. Aborting!!!')
            raise RemapperError(msg=msg)

        remap_obj = parser_interface.object_setattr(
            object_in=remap_obj, key="interp_scheme", value=interp_scheme)

        # Define the vertical level attributes.
        nlevs = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key="nlevs", force=True)
        if nlevs is None:
            msg = ('The total number of levels for the remapped/interpolated '
                   'variable(s) could not be determined from experiment '
                   'configuration file {self.yaml_file}. Aborting!!!')
            raise RemapperError(msg=msg)

        remap_obj = parser_interface.object_setattr(
            object_in=remap_obj, key="nlevs", value=nlevs)

        return remap_obj

    def getvar_info(self) -> Tuple[object, Union[object, None]]:
        """
        Description
        -----------

        This method parses the user experiment configuration and
        defines a Python object containing the respective variable(s)
        to be remapped attributes and the respective forecast model
        from the respective variables originate.

        Returns
        -------

        varinfo_obj: object

            A Python object containing the respective variable(s) to
            be remapped attributes and the respective forecast model
            from which the respective variables originate.

        soca_varinfo_obj: object or NoneType

            A Python object containing the respective SOCA variable(s)
            to be remapped attributes and the respective forecast
            model from which the respective variables originate;
            currently this return applies to only the GLORYS forecast
            model; for other forecast model variables the returned
            value is NoneType.

        """

        # Define the variable and respective forecast model attributes
        # accordingly.
        (varinfo_obj, soca_varinfo_obj) = \
            [parser_interface.object_define(), None]
        varinfo_obj = parser_interface.object_setattr(
            object_in=varinfo_obj, key="forecast_model", value=self.forecast_model)

        # Collect the bathymetry arguments from the experiment
        # configuration file.
        bathy_attrs_list = ['bathy_edits_file', 'bathy_file']

        for bathy_attr in bathy_attrs_list:
            value = parser_interface.dict_key_value(
                dict_in=self.yaml_dict, key=bathy_attr, force=True, no_split=True)
            varinfo_obj = parser_interface.object_setattr(
                object_in=varinfo_obj, key=bathy_attr, value=value)

        # Build the Python objects containing the variable attributes.
        variables_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key="variables", force=True)
        if variables_dict is None:
            msg = ('The variables to be remapped and their respective attributes '
                   'could not be determined from the user experiment configuration. '
                   'Aborting!!!')
            raise RemapperError(msg=msg)

        variable_list = []

        for variable in variables_dict:
            vardict = parser_interface.dict_key_value(
                dict_in=variables_dict, key=variable, no_split=True)
            varinfo_obj = parser_interface.object_setattr(
                object_in=varinfo_obj, key=variable, value=vardict)
            variable_list.append(variable)

        varinfo_obj = parser_interface.object_setattr(
            object_in=varinfo_obj, key="variable_list", value=variable_list)

        return (varinfo_obj, soca_varinfo_obj)

    def remap(self) -> None:
        """ """
        forecast_model = self.varinfo_obj.forecast_model.lower()
        if self.forecast_model in self.models_dict.keys():
            model_obj = parser_interface.dict_key_value(
                dict_in=self.models_dict, key=self.forecast_model, no_split=True)

            forecast_model_obj = model_obj(srcgrid_obj=self.srcgrid_obj,
                                           dstgrid_obj=self.dstgrid_obj,
                                           remap_obj=self.remap_obj,
                                           varinfo_obj=self.varinfo_obj)
            forecast_model_obj.run()

        else:
            msg = (f"The forecast model {self.forecast_model} is not supported. "
                   "Aborting!!!")
            raise RemapperError(msg=msg)

    def run(self) -> None:
        """ """
        self.srcgrid_obj = self.getgrid_info(grid_name="source")
        self.dstgrid_obj = self.getgrid_info(grid_name="destination")
        self.remap_obj = self.getremap_info()
        (self.varinfo_obj, self.soca_varinfo_obj) = self.getvar_info()

        self.remap()
