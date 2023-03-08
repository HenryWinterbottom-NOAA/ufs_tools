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

- ufs_pytils; https://github.com/HenryWinterbottom-NOAA/ufs_pyutils]

Author(s)
---------

    Henry R. Winterbottom; 08 March 2023

History
-------

    2023-03-08: Henry Winterbottom -- Initial implementation.

"""

# ----

from dataclasses import dataclass
from typing import Dict

import numpy
from exceptions import TCDiagsError
from ioapps import netcdf4_interface
from tools import parser_interface
from utils.logger_interface import Logger

# ----


@dataclass
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

        # Check that the input file contains all of the mandatory
        # variables.
        read_inputs_dict = {
            "latitude": "lat_degrees",
            "longitude": "lon_degrees",
            "temperature": "tmp_K",
        }

        mand_inputs = set(list(sorted(read_inputs_dict.keys())))
        yaml_inputs = set(list(sorted(self.inputs_dict.keys())))
        missing_vars = list(sorted(mand_inputs - yaml_inputs))

        if len(missing_vars) > 0:
            msg = (
                "The following mandatory input variables could not be found in the "
                f"experiment configuration: {missing_vars}. Aborting!!!"
            )
            raise TCDiagsError(msg=msg)

        # Define the mandatory variables object; proceed accordingly.
        inputs_obj = parser_interface.object_define()

        # Define the variable attributes.
        varin_attrs_dict = {
            "flip_lat": False,
            "flip_z": False,
            "ncfile": None,
            "ncvarname": None,
            "scale_add": 0.0,
            "scale_mult": 1.0,
            "squeeze": False,
            "squeeze_axis": 0,
        }

        # Build the input variables object; proceed accordingly.
        for yaml_key in yaml_inputs:
            varin_obj = parser_interface.object_define()

            var_dict = parser_interface.dict_key_value(
                dict_in=self.inputs_dict, key=yaml_key, force=True, no_split=True
            )

            for varin_attr in varin_attrs_dict:
                value = parser_interface.dict_key_value(
                    dict_in=var_dict, key=varin_attr, force=True, no_split=True
                )

                if value is None:
                    value = parser_interface.dict_key_value(
                        dict_in=varin_attrs_dict,
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
            msg = f"Variable {yaml_key} range values: ({values.min()}, {values.max()})."
            self.logger.debug(msg=msg)

            # Update the input variables object.
            inputs_obj = parser_interface.object_setattr(
                object_in=inputs_obj, key=yaml_key, value=values
            )

        return inputs_obj
