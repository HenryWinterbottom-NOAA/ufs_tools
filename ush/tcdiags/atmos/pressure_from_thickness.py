from typing import Dict
import numpy

from metpy.units import units

from tools import parser_interface
from ioapps import netcdf4_interface


def pressure_from_thickness(inputs_obj: object) -> numpy.array:
    """ """

    # Initialize the pressure profile.
    dpres = numpy.array(inputs_obj.pres)
    pres = numpy.array(inputs_obj.pres)
    pres[0, :, :] = numpy.array(inputs_obj.psfc)[:, :]

    # Compute the pressure profile using the surface pressure and
    # layer thickness; proceed accordingly.
    for zlev in range(pres.shape[0] - 2, 0, -1):

        # Compute the pressure profile.
        pres[zlev, :, :] = pres[zlev + 1, :, :] + dpres[zlev, :, :]

    # Correct units and update the input variable object.
    pres = units.Quantity(pres, "Pa")

    inputs_obj = parser_interface.object_setattr(
        object_in=inputs_obj, key="pres", value=pres)

    return inputs_obj
