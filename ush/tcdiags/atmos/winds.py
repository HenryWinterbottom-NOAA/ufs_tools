# =========================================================================

# Module: ush/tcdiags/atmos/winds.py

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

    winds.py

Description
-----------

    This module contains functions to compute total wind field
    diagnostics.

Functions
---------

    _cleanup(xspharm)

        This function attempts to destroy the initialized spherical
        harmonic transform object and subsequently collect any
        remaining/residual object attributes.

    _get_lev_uv(inputs_obj, lev)

        This function returns the total wind vector component
        quantities at the specified vertical level.

    _init_spharm(array_in)

        This method initializes the spherical harmonic transform
        object.

    global_divg(inputs_obj)

        This function computes the global divergence field using
        spherical harmonic transforms.

    global_vort(inputs_obj)

        This function computes the global vorticity field using
        spherical harmonic transforms.

    global_wind_part(inputs_obj)

        This function computes the components of the total wind field
        using spherical harmonic transforms; the methodology follows
        from Lynch [1988].

Requirements
------------

- metpy; https://unidata.github.io/MetPy/latest/index.html

- pyspharm; https://github.com/jswhit/pyspharm

- ufs_pytils; https://github.com/HenryWinterbottom-NOAA/ufs_pyutils

Author(s)
---------

    Henry R. Winterbottom; 09 March 2023

History
-------

    2023-03-09: Henry Winterbottom -- Initial implementation.

"""

# ----

# pylint: disable=invalid-name
# pylint: disable=too-many-locals

# ----

import gc
from typing import Tuple, Union

import numpy
import spharm
from exceptions import TCDiagsError
from metpy.units import units
from pint import Quantity
from tools import parser_interface
from utils.logger_interface import Logger

# ----

# Define all available functions.
__all__ = ["global_divg", "global_vort", "global_wind_part"]

# ----

logger = Logger()

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


def _cleanup(xspharm: spharm.Spharmt) -> None:
    """
    Description
    -----------

    This function attempts to destroy the initialized spherical
    harmonic transform object and subsequently collect any
    remaining/residual object attributes.

    Parameters
    ----------

    xspharm: object

        A Python object containing the initialized spherical harmonic
        transform object.

    """

    # Destroy the initialized spherical harmonic transform object.
    del xspharm
    gc.collect()


# ----


def _get_lev_uv(inputs_obj: object, lev: int) -> Tuple[Quantity, Quantity]:
    """
    Description
    -----------

    This function returns the total wind vector component quantities
    at the specified vertical level.

    Parameters
    ----------

    inputs_obj: object

        A Python object containing, at minimum, the zonal and
        meridional total wind components; units are meters per second.

    lev: int

        A Python integer specifying the vertical level from which to
        collect the zonal and meridional total wind components.

    Returns
    -------

    xuwnd: array-type or pint.Quantity

        A Python array-type variable containing the zonal wind
        component collected for the specified vertical level.

    xvwnd: array-type or pint.Quantity

        A Python array-type variable containing the meridional wind
        component collected for the specified vertical level.

    """

    # Define the zonal and meridional total wind components.
    (uwnd, vwnd) = [
        parser_interface.object_getattr(object_in=inputs_obj, key=wndvar)
        for wndvar in ["uwnd", "vwnd"]
    ]

    # Define the zonal and meridional total wind components for the
    # respective/specified vertical level; proceed accordingly.
    if uwnd.ndim == 3 and vwnd.ndim == 3:
        (xuwnd, xvwnd) = [uwnd[lev, :, :], vwnd[lev, :, :]]

    elif uwnd.ndim == 2 and vwnd.ndim == 2:
        (xuwnd, xvwnd) = [uwnd[:, :], vwnd[:, :]]

    else:
        msg = (
            "The wind vector components could not be parsed and/or "
            "the level attribute is invalid. Aborting!!!"
        )
        raise TCDiagsError(msg=msg)

    return (xuwnd, xvwnd)


# ----


def _init_spharm(array_in: Union[numpy.array, Quantity]) -> spharm.Spharmt:
    """
    Description
    -----------

    This method initializes the spherical harmonic transform object.

    Parameters
    ----------

    array_in: array-type or pint.Quantity

        A Python array-type variable; the shape of the respective
        array will be used to initialize the spherical harmonic
        transform object.

    Returns
    -------

    xspharm: object

        A Python object containing the initialized spherical harmonic
        transform object.

    """

    # Initialize the spherical harmonic transform object accordingly.
    if array_in.ndim == 3:
        (nx, ny) = (array_in.shape[2], array_in.shape[1])

    else:
        (nx, ny) = (array_in.shape[1], array_in.shape[2])

    xspharm = spharm.Spharmt(nx, ny)

    return xspharm


# ----


def global_divg(inputs_obj: object) -> object:
    """
    Description
    -----------

    This function computes the global divergence field using spherical
    harmonic transforms.

    Parameters
    ----------

    inputs_obj: object

        A Python object containing, at minimum, the zonal and
        meridional wind components; units are meters per second.

    Returns
    -------

    inputs_obj: object

        A Python object updated to now contain the global divergence
        field `divg`.

    """

    # Initialize the local variable and objects.
    xdivg = numpy.zeros(inputs_obj.uwnd.shape)
    xspharm = _init_spharm(array_in=xdivg)
    nlevs = xdivg.shape[0]

    # Compute the divergence field; proceed accordingly.
    msg = f"Computing global diverenge array of dimension {xdivg.shape}."
    logger.info(msg=msg)

    for lev in range(nlevs):
        (u, v) = _get_lev_uv(inputs_obj=inputs_obj, lev=lev)

        (_, dataspec) = xspharm.getvrtdivspec(ugrid=u, vgrid=v)
        xdivg[lev, :, :] = xspharm.spectogrd(dataspec=dataspec)

    # Define the correct units with respect to the input variable.
    xdivg = units.Quantity(xdivg, "1 / second")

    inputs_obj = parser_interface.object_setattr(
        object_in=inputs_obj, key="divg", value=xdivg
    )

    msg = (
        f"Global divergence values range({numpy.array(xdivg).min()}, "
        f"{numpy.array(xdivg).max()}) {xdivg.units}."
    )
    logger.debug(msg=msg)

    # Deallocate memory for the spherical harmonic transform object.
    _cleanup(xspharm=xspharm)

    return inputs_obj


# ----


def global_vort(inputs_obj: object) -> object:
    """
    Description
    -----------

    This function computes the global vorticity field using spherical
    harmonic transforms.

    Parameters
    ----------

    inputs_obj: object

        A Python object containing, at minimum, the zonal and
        meridional wind components; units are meters per second.

    Returns
    -------

    inputs_obj: object

        A Python object updated to now contain the global vorticity
        field `vort`.

    """

    # Initialize the local variable and objects.
    xvort = numpy.zeros(inputs_obj.uwnd.shape)
    xspharm = _init_spharm(array_in=xvort)
    nlevs = xvort.shape[0]

    msg = f"Computing global vorticity array of dimension {xvort.shape}."
    logger.info(msg=msg)

    # Compute the vorticity field; proceed accordingly.
    for lev in range(nlevs):
        (u, v) = _get_lev_uv(inputs_obj=inputs_obj, lev=lev)

        (dataspec, _) = xspharm.getvrtdivspec(ugrid=u, vgrid=v)
        xvort[lev, :, :] = xspharm.spectogrd(dataspec=dataspec)

    # Define the correct units with respect to the input variable.
    xvort = units.Quantity(xvort, "1 / second")

    inputs_obj = parser_interface.object_setattr(
        object_in=inputs_obj, key="vort", value=xvort
    )

    msg = (
        f"Global vorticity values range({numpy.array(xvort).min()}, "
        f"{numpy.array(xvort).max()}) {xvort.units}."
    )
    logger.debug(msg=msg)

    # Deallocate memory for the spherical harmonic transform object.
    _cleanup(xspharm=xspharm)

    return inputs_obj


# ----


def global_wind_part(inputs_obj: object) -> object:
    """
    Description
    -----------

    This function computes the components of the total wind field
    using spherical harmonic transforms; the methodology follows from
    Lynch [1988].

    Parameters
    ----------

    inputs_obj: object

        A Python object containing, at minimum, the zonal and
        meridional wind components; units are meters per second.

    Returns
    -------

    inputs_obj: object

        A Python object updated to now contain the global total wind
        field paritioning; the respective quantities are as follows.

        - (udiv, vdiv); the components of the divergent component of
          the global total wind field.

        - (uhrm, vhrm); the components of the harmonic (i.e.,
          residual) component of the global total wind field.

        - (uvor, vvor); the components of the rotational component of
          the global total wind field.

    References
    ----------

    Lynch, P., 1994: Paritioning the Wind in a Limited
    Domain. Mon. Wea. Rev., 116, 86-93.

    https://doi.org/10.1175/1520-0493(1988)116<0086:DTWFVA>2.0.CO;2

    """

    # Initialize the local variable and objects.
    (xudiv, xuhrm, xuvor, xvdiv, xvhrm, xvvor) = [
        numpy.zeros(inputs_obj.uwnd.shape) for idx in range(6)
    ]

    xspharm = _init_spharm(array_in=xuvor)
    nlevs = xuvor.shape[0]

    msg = f"Computing global partitioned total wind arrays of dimension {xuvor.shape}."
    logger.info(msg=msg)

    # Compute the total wind components; proceed accordingly.
    for lev in range(nlevs):
        (u, v) = _get_lev_uv(inputs_obj=inputs_obj, lev=lev)
        (zspec_save, dspec_save) = xspharm.getvrtdivspec(ugrid=u, vgrid=v)

        # Compute the divergent component of the total wind field.
        (zspec, dspec) = [numpy.zeros(zspec_save.shape), dspec_save]
        (xudiv[lev, :, :], xvdiv[lev, :, :]) = xspharm.getuv(
            vrtspec=zspec, divspec=dspec
        )

        # Compute the rotational component of the total wind field.
        (zspec, dspec) = [zspec_save, numpy.zeros(dspec_save.shape)]
        (xuvor[lev, :, :], xvvor[lev, :, :]) = xspharm.getuv(
            vrtspec=zspec, divspec=dspec
        )

        # Compute the residual (i.e., harmonic) component of the total
        # wind field.
        xuhrm[lev, :, :] = numpy.array(
            u[:, :]) - (xuvor[lev, :, :] + xudiv[lev, :, :])
        xvhrm[lev, :, :] = numpy.array(
            v[:, :]) - (xvvor[lev, :, :] + xvdiv[lev, :, :])

    # Define the correct units with respect to the input variable.
    (xudiv, xuhrm, xuvor, xvdiv, xvhrm, xvvor) = [
        inputs_obj.uwnd.units for wcmpn in [xudiv, xuhrm, xuvor, xvdiv, xvhrm, xvvor]
    ]

    inputs_obj = [
        parser_interface.object_setattr(
            object_in=inputs_obj, key=key, value=value)
        for key in ["udiv", "uhrm", "uvor", "vdiv", "vhrm", "vvor"]
        for value in [xudiv, xuhrm, xuvor, xvdiv, xvhrm, xvvor]
    ]

    # Deallocate memory for the spherical harmonic transform object.
    _cleanup(xspharm=xspharm)

    return inputs_obj
