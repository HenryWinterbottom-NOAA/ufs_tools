# =========================================================================

# Module: ush/tcdiags/metrics/tcmpi.py

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

from metpy.calc import mixing_ratio_from_specific_humidity as mrfsh
from metpy.units import units

from tcpyPI import pi

# ----

#logger = Logger()

# ----

__author__ = "Henry R. Winterbottom"
__maintainer__ = "Henry R. Winterbottom"
__email__ = "henry.winterbottom@noaa.gov"

# ----


def tcmpi(spfh: numpy.array) -> object:
    """ """

    mix_ratio = mrfsh(spfh)

    print(mix_ratio)
