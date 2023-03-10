# =========================================================================

# Module: ush/tcdiags/io/__init__.py

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

    This module contains the base-class object for all
    netCDF-formatted file and respective variable(s) reading and
    writing.

Classes
-------

    TCDiagsIO(yaml_dict)

        This is the base-class object for all netCDF-formatted file
        and respective variable(s) reading and writing.

Requirements
------------

- ufs_pytils; https://github.com/HenryWinterbottom-NOAA/ufs_pyutils

Author(s)
---------

    Henry R. Winterbottom; 08 March 2023

History
-------

    2023-03-08: Henry Winterbottom -- Initial implementation.

"""

# ----

from typing import Dict

from tcdiags.atmos.heights import height_from_pressure
from tcdiags.atmos.pressures import pressure_from_thickness
from metpy.units import units
import numpy
from exceptions import TCDiagsError
from ioapps import netcdf4_interface
from tools import parser_interface
from utils.logger_interface import Logger

# ----

# Check that the input file contains all of the mandatory
# variables.
INPUTS_DICT = {
    "latitude": {"name": "lat", "units": "degree"},
    "longitude": {"name": "lon", "units": "degree"},
    "pressure": {"name": "pres", "units": "pascals"},
    "specific_humidity": {"name": "spfh", "units": "kg/kg"},
    "surface_height": {"name": "zsfc", "units": "gpm"},
    "surface_pressure": {"name": "psfc", "units": "Pa"},
    "temperature": {"name": "temp", "units": "K"},
    "uwind": {"name": "uwnd", "units": "m/s"},
    "vwind": {"name": "vwnd", "units": "m/s"},
}

# Define the variable attributes.
VARIN_ATTRS_DICT = {
    "flip_lat": False,
    "flip_z": False,
    "method": -99,
    "ncfile": None,
    "ncvarname": None,
    "scale_add": 0.0,
    "scale_mult": 1.0,
    "squeeze": False,
    "squeeze_axis": 0,
}

# Pressure profile computation methods.
PRES_PROF_COMP_METHODS_DICT = {
    1: pressure_from_thickness
}

# ----


class TCDiagsIO:
    """
    Description
    -----------

    This is the base-class object for all netCDF-formatted file and
    respective variable(s) reading and writing.

    Parameters
    ----------

    inputs_dict: dict

        A Python dictionary containing the input attributes as defined
        by the `inputs_yaml` attribute within the experiment
        configuration.

    Raises
    ------

    TCDiagsError:

        * raised if the attribute `inputs` can not be determined from
          the experiment configuration file.

    """

    def __init__(self, yaml_dict: Dict):
        """
        Description
        -----------

        Creates a new TCDiagsIO object.

        """

        # Define the base-class attributes.
        self.logger = Logger()
        self.yaml_dict = yaml_dict
        self.inputs_dict = parser_interface.dict_key_value(
            dict_in=self.yaml_dict, key="inputs_yaml", force=True
        )

        if self.inputs_dict is None:
            msg = (
                "The attribute `inputs` could not be determined from the experiment "
                "configuration file. Aborting!!!"
            )
            raise TCDiagsError(msg=msg)

        self.variable_range_msg = "Variable %s range values: (%s, %s) %s."

    def _get_pressure(self, inputs_obj: object) -> object:
        """
        Description
        -----------

        This method computes the pressure profile in accordance with
        the experiment configuration attributes.

        Parameters
        ----------

        inputs_obj: object

            A Python object containing the mandatory input variables;
            this includes, at minimum, the pressure and surface
            pressure values; see below for additional information.

        Returns
        -------

        inputs_obj: object

            A Python object updated to contain the pressure profile
            array in accordance with the method defined in the
            experiment configuration.

        Raises
        ------

        TCDiagsError:

            * raised if the pressure variable attributes cannot be
              determined from the experiment configuration or are not
              defined within this module (see
              `PRES_PROF_COMP_METHODS_DICT` above).

        Notes
        -----

        For the respective pressure profile computations, the
        respective methods assume the pressure array upon entry
        contains the following:

        1: pressure_from_thickness :: pres = the layer thickness; this
           is used to derive the pressure profile by integrating from
           the top layer thickness to the surface.

        """

        # Define the method to be used for the pressure profile
        # computation; proceed accordingly.
        var_dict = parser_interface.dict_key_value(
            dict_in=self.inputs_dict, key="pressure", force=True, no_split=True)

        if var_dict is None:
            msg = ("The pressure variable attributes could not be determined from "
                   "the experiment configuration. Aborting!!!"
                   )
            raise TCDiagsError(msg=msg)

        method = parser_interface.dict_key_value(
            dict_in=var_dict, key="method", force=True)
        if method is None:
            msg = ("The attribute `method` could not be determined for variable "
                   "`pressure` in the experiment configuration. Aborting!!!"
                   )
            raise TCDiagsError(msg=msg)

        # Compute the pressure profile accordingly.
        app = parser_interface.dict_key_value(
            dict_in=PRES_PROF_COMP_METHODS_DICT, key=method, force=True)
        if app is None:
            msg = f"The pressure profile method {method} is not defined. Aborting!!!"
            raise TCDiagsError(msg=msg)

        inputs_obj = app(inputs_obj=inputs_obj)

        msg = (self.variable_range_msg % ("pressure", numpy.array(inputs_obj.pres).min(),
                                          numpy.array(inputs_obj.pres).max(),
                                          inputs_obj.pres.units))
        self.logger.debug(msg=msg)

        return inputs_obj

    def read_inputs(self) -> object:
        """
        Description
        -----------

        This method collects the mandatory variables from the
        experiment configuration specified YAML-formatted file
        containing the input variable attributes.

        Returns
        -------

        inputs_obj: object

            A Python object containing the mandatory input variables.

        Raises
        ------

        TCDiagsError:

            * raised if a mandatory input variable has not been
              specified in the YAML-formatted file containing the
              input variable attributes.

            * raised if a required input variable attribute is
              NoneType after parsing the YAML-formatted file
              containing the input variable attributes.

        """

        # Check that all required variables are defined; proceed
        # accordingly.
        mand_inputs = set(list(sorted(INPUTS_DICT.keys())))
        yaml_inputs = set(list(sorted(self.inputs_dict.keys())))
        missing_vars = list(sorted(mand_inputs - yaml_inputs))

        if len(missing_vars) > 0:
            msg = (
                "The following mandatory input variables could not be found in the "
                f"experiment configuration: {missing_vars}. Aborting!!!"
            )
            raise TCDiagsError(msg=msg)

        # Build the input variables object; proceed accordingly.
        inputs_obj = parser_interface.object_define()

        for yaml_key in INPUTS_DICT:
            varin_obj = parser_interface.object_define()

            var_dict = parser_interface.dict_key_value(
                dict_in=self.inputs_dict, key=yaml_key, force=True, no_split=True
            )

            for varin_attr in VARIN_ATTRS_DICT:
                value = parser_interface.dict_key_value(
                    dict_in=var_dict, key=varin_attr, force=True, no_split=True
                )

                if value is None:
                    value = parser_interface.dict_key_value(
                        dict_in=VARIN_ATTRS_DICT,
                        key=varin_attr,
                        force=True,
                        no_split=True,
                    )

                if value is None:
                    msg = (
                        f"The mandatory attribute {varin_attr} for variable {yaml_key} "
                        "could not be determined from the experiment configuration. "
                        "Aborting!!!"
                    )
                    raise TCDiagsError(msg=msg)

                varin_obj = parser_interface.object_setattr(
                    object_in=varin_obj, key=varin_attr, value=value
                )

            # Collect the respective variable and scale as necessary.
            msg = f"Reading variable {yaml_key} from netCDF-formatted file path {varin_obj.ncfile}."
            self.logger.info(msg=msg)

            values = netcdf4_interface.ncreadvar(
                ncfile=varin_obj.ncfile,
                ncvarname=varin_obj.ncvarname,
                squeeze=varin_obj.squeeze,
                axis=varin_obj.squeeze_axis,
            )

            # Manipulate the input variable values/grid projection
            # accordingly.
            if varin_obj.flip_z:
                if len(values.shape) > 2:
                    values = numpy.flip(values[:, :, :], axis=0)
                    if varin_obj.flip_lat:
                        values = numpy.flip(values[:, :, :], axis=1)

                else:
                    if varin_obj.flip_lat:
                        try:
                            values = numpy.flip(values[:, :], axis=0)
                        except IndexError:
                            pass

            values = varin_obj.scale_mult * (values) + varin_obj.scale_add

            # Define the variable object accordingly.
            (var_name, var_units) = [parser_interface.dict_key_value(
                dict_in=INPUTS_DICT[yaml_key], key=key, no_split=True) for key in ["name", "units"]]

            values = units.Quantity(values, var_units)
            inputs_obj = parser_interface.object_setattr(
                object_in=inputs_obj, key=var_name, value=values)

            msg = (self.variable_range_msg % (yaml_key, numpy.array(values.min()),
                                              numpy.array(values.max()), values.units))
            self.logger.debug(msg=msg)

        # Compute/define the remaining diagnostic variables.
        inputs_obj = self._get_pressure(inputs_obj=inputs_obj)
        inputs_obj = height_from_pressure(inputs_obj=inputs_obj)

        values = inputs_obj.hght
        msg = (self.variable_range_msg % ("height", numpy.array(values.min()),
                                          numpy.array(values.max()), values.units))
        self.logger.debug(msg=msg)

        return inputs_obj
