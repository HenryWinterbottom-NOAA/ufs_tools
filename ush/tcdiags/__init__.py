# =========================================================================

# Module: ush/tcdiags/__init__.py

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


# ----

from dataclasses import dataclass

from confs.yaml_interface import YAML

from tcdiags.io import TCDiagsIO

from exceptions import TCDiagsError
from utils.logger_interface import Logger

# ----


@dataclass
class TCDiags:
    """

    """

    def __init__(self, options_obj: object):
        """
        Description
        -----------

        Creates a new TCDiags object.

        """

        # Define the base-class attributes.
        self.options_obj = options_obj
        self.logger = Logger()
        self.yaml_file = self.options_obj.yaml_file
        self.yaml_dict = YAML().read_yaml(yaml_file=self.yaml_file)
        self.tcdiags_io = TCDiagsIO(yaml_dict=self.yaml_dict)

    def run(self) -> None:
        """
        """

        inputs_obj = self.tcdiags_io.read_inputs()

        # print(inputs_obj)
