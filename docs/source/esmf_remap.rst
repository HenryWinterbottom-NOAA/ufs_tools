#########################################
Building ESMF Remapping Coefficient Files
#########################################

The Earth System Modeling Framework (ESMF) applications provide a
means to compute various types (e.g., bilinear, nearest-neighbor,
etc.,) remapping coefficient values. The ``esmf_remap`` application
provides a script-level interface around the `xESMF
<https://xesmf.readthedocs.io/en/latest/>`_ package.

^^^^^^^^^^^^^^^^^^^^^^^^^
Application Configuration
^^^^^^^^^^^^^^^^^^^^^^^^^

The ``esmf_remap`` application for a remapping of :math:`0.25^{\circ}`
grid-cell centers :math:`5.0^{\circ}` grid cell centers is as follows.

::

     interp_type: bilinear
     grids:

          destination:
               is_grib: False
	       is_netcdf: True
               gribfile: Null
               ncfile: /run/gridspec.tripolar.0p25.nc
	       nclat: tlat
               nclon: tlon
               ncxdim_name: nx
               ncydim_name: ny          

          source:
               is_grib: False
	       is_netcdf: True
               gribfile: Null
               ncfile: /run/gridspec.tripolar.5p0.nc
	       nclat: tlat
               nclon: tlon
               ncxdim_name: nx
               ncydim_name: ny

      output_netCDF: /run/tripolar.0p25_5p0.bilinear.t2t_remap.nc


+------------------------------+---------------------------------------------------------------------------+
| Configuration Variable       | Description                                                               |
+==============================+===========================================================================+
| ``interp_type``              | The interpolation type for which to compute the ESMF attributes; both     |
|                              | ``bilinear`` and ``nearest_s2d`` are supported.                           |
+------------------------------+---------------------------------------------------------------------------+
| ``grids``                    | The respective ``source`` and ``destination`` grid informations; see the  |
|                              | following table for further information.                                  | 
+------------------------------+---------------------------------------------------------------------------+
| ``output_netcdf``            | The netCDF-formatted file path to contain the computed ESMF attributes.   |
+------------------------------+---------------------------------------------------------------------------+

The allowable attributes for the respective ``source`` and
``destination`` grid-types are as follows.

+------------------------------+---------------------------------------------------------------------------+
| Grid-type Variable           | Description                                                               |
+==============================+===========================================================================+
| ``is_grib``                  | Specifies whether the respective grid is contained/defined within a       |
|                              | `WMO GRIB <https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/>`_     |
|                              | formatted file path; this option is not yet supported.                    |
+------------------------------+---------------------------------------------------------------------------+
| ``is_netcdf``                | Specifies whether the respective grid is contained/defined within a       |
|                              | `netCDF-formatted <https://www.unidata.ucar.edu/software/netcdf/>`_ file  |
|                              | path.                                                                     |
+------------------------------+---------------------------------------------------------------------------+
| ``gribfile``                 | Defines the WMO-GRIB formatted file path containing the respective grid;  |
|                              | used only if ``is_grib`` is ``True``.                                     |
+------------------------------+---------------------------------------------------------------------------+
| ``ncfile``                   | Defines the netCDF-formatted file path containing the respective grid     |
|                              | attributes; used only if ``is_netcdf`` is ``True``.                       |
+------------------------------+---------------------------------------------------------------------------+
| ``nclat``                    | Defines the netCDF-formatted file path variable name for the latitude     |
|                              | coordinate values.                                                        |
+------------------------------+---------------------------------------------------------------------------+
| ``nclon``                    | Defines the netCDF-formatted file path variable name for the longitude    |
|                              | coordinate values.                                                        |
+------------------------------+---------------------------------------------------------------------------+
| ``ncxdim_name``              | Defines the netCDF-formatted file path `x`-dimension variable name.       |
+------------------------------+---------------------------------------------------------------------------+
| ``ncydim_name``              | Defines the netCDF-formatted file path `y`-dimension variable name.       |
+------------------------------+---------------------------------------------------------------------------+

^^^^^^^^^^^^^^^^^^^^^^^^^
Launching the Application
^^^^^^^^^^^^^^^^^^^^^^^^^

The ``esmf_remap`` application may be launched from within the UFS
Tools package for the example described in the previous section as
follows.

::

   user@host:$ cd scripts
   user@host:$ python compute_esmf_remap.py --yaml_file ../parm/esmf_remap.yaml

If successful, a netCDF-formatted file containing the ESMF remapping
attributes will be created containing the following metadata
attributes.

::

   user@host:$ ncdump -h /run/tripolar.0p25_5p0.bilinear.t2t_remap.nc

   netcdf tripolar.0p25_5p0.bilinear.t2t_remap {
   dimensions:
           n_s = 113472 ;
   variables:
   	  double S(n_s) ;
		  S:_FillValue = NaN ;
	  int col(n_s) ;
	  int row(n_s) ;
   }

A description of the respective ESMF remapping variables/attributes
can be found `here
<https://earthsystemmodeling.org/docs/release/ESMF_8_0_1/ESMF_refdoc/node3.html#regridoutput>`_.

In addition to the stand-alone script example above, a wrapper
application around the above-mentioned script is included with the UFS
Tools package and allows the computation of various remapping
attributes to run in succession. This

::

   user@host:$ cd wrappers
   user@host:$ python wrapper_esmf_remap.py --yaml_file ../parm/esmf_remap/wrapper_esmf_remap.yaml --yaml_template ../parm/esmf_remap/esmf_remap_tmpl.yaml --script_path ../scripts/compute_esmf_remap.py
