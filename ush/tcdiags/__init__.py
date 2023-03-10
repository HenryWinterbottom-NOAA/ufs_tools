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

from tcdiags.tc import FilterVortex

from tools import parser_interface

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

        # Define the available application options.
        self.apps_dict = {
            "tcfilt": FilterVortex}

    def run(self) -> None:
        """
        """

        inputs_obj = self.tcdiags_io.read_inputs()

        # Execute each of the specified applications.
        for app in self.apps_dict:

            # Check whether the application is to be executed; proceed
            # accordingly.
            opt_attr = parser_interface.object_getattr(
                object_in=self.options_obj, key=app, force=True)

            if opt_attr is not None:

                if parser_interface.str_to_bool(opt_attr):

                    # Launch the respective application.
                    app_class = parser_interface.dict_key_value(
                        dict_in=self.apps_dict, key=app, no_split=True)

                    app_class(inputs_obj=inputs_obj).run()
