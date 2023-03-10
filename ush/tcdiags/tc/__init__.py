from dataclasses import dataclass

from tcdiags.atmos import winds

from exceptions import TCDiagsError
from tools import parser_interface


@dataclass
class FilterVortex:
    """

    """

    def __init__(self, inputs_obj: object):
        """ 
        Description
        -----------

        Creates a new FilterVortex object.

        """

        # Define the base-class attributes.
        self.inputs_obj = inputs_obj
        (self.thermo_obj, self.winds_obj) = [
            parser_interface.object_define() for idx in range(2)]

        # Initialize the local Python object containing the wind field
        # attributes; proceed accordingly.
        wind_comps_list = ["uwnd", "vwnd"]

        for wind_comp in wind_comps_list:
            value = parser_interface.object_getattr(
                object_in=self.inputs_obj, key=wind_comp, force=True)

            if value is None:
                msg = (f"The mandatory wind component variable {wind_comp} "
                       "could not be determined from the input variable "
                       "attributes object `inputs_obj`. Aborting!!!"
                       )
                raise TCDiagsError(msg=msg)

            self.winds_obj = parser_interface.object_setattr(
                object_in=self.winds_obj, key=wind_comp, value=value)

    def compute_winds(self):
        """ """

        # Compute the wind field components.
        self.winds_obj = winds.global_divg(inputs_obj=self.winds_obj)
        self.winds_obj = winds.global_vort(inputs_obj=self.winds_obj)

    def run(self) -> None:
        """

        """

        self.compute_winds()
