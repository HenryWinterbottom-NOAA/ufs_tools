.. role:: red
  :class: red

#########################################
Remapping Analyses and Initial Conditions
#########################################

Analysis remapping, supported by the available remapping types (i.e.,
ESMF, SLINT, etc.,), can be accomplished by using the ``remapper``
application. For optimal results using the ``remapper`` application,
the user should first generate the tripolar grid projections using the
``gridspec`` for the respective grid-scale resolution. Next the ESMF
remapping coefficients should then be computed using the
``esmf_remap`` application.

^^^^^^^^^^^^^^^^^^^^^^^^^
Application Configuration
^^^^^^^^^^^^^^^^^^^^^^^^^

The ``remapper`` application configuration for a :math:`0.25^{\circ}`
`MOM6 <https://www.gfdl.noaa.gov/mom-ocean-model/>`_ analysis to a
:math:`1.0^{\circ}` MOM6 tripolar grid remapping for sea-surface
height and the current velocity vector components is as follows.

::

   bathy_edits_file: /path/to/topo/ocean_topog_edits.1p0.nc
   bathy_file: /path/to/topo/ocean_topog.1p0.nc
   forecast_model: mom6

   grids:	  
	destination:
	     arakawa: C  
             grid_ncfile: /path/to/gridspec/gridspec.tripolar.1p0.nc      
             mass:
                  nclat: tlat
		  nclon: tlon
             topo_ncfile: /path/to/topo/ocean_topog.1p0.nc  
	     uvel:
	          nclat: ulat
		  nclon: ulon
	     vvel:
                  nclat: vlat
		  nclon: vlon   
        source:
	     arakawa: C
             grid_ncfile: /path/to/gridspec/gridspec.tripolar.0p25.nc
             mass:
                  nclat: tlat
		  nclon: tlon
             topo_ncfile: /path/to/topo/ocean_topog.0p25.nc
	     uvel:
                  nclat: ulat
		  nclon: ulon
	     vvel:
	          nclat: vlat
		  nclon: vlon

   interp_scheme: esmf
   nlevs: 75
   output_netCDF: /path/to/remapper/outputs/MOM6.ic.1p0.nc

   remap:
        dstmass2dstuvel_bilinear: /path/to/esmf_remap/tripolar.1p0_1p0.bilinear.t2u_remap.nc
	dstmass2dstvvel_bilinear: /path/to/esmf_remap/tripolar.1p0_1p0.bilinear.t2v_remap.nc
        rotate_currents: True
	srcmass2dstmass_bilinear: /path/to/esmf_remap/tripolar.0p25_1p0.bilinear.t2t_remap.nc
	srcuvel2dstmass_bilinear: /path/to/esmf_remap/tripolar.0p25_1p0.bilinear.u2t_remap.nc
	srcvvel2dstmass_bilinear: /path/to/esmf_remap/tripolar.0p25_1p0.bilinear.v2t_remap.nc
	srcuvel2srcmass_bilinear: /path/to/esmf_remap/tripolar.0p25_0p25.bilinear.u2t_remap.nc
	srcvvel2srcmass_bilinear: /path/to/esmf_remap/tripolar.0p25_0p25.bilinear.v2t_remap.nc

   variables:

        sea_surface_height:
             src_ncvarname: sfc
             dst_ncvarname: sfc          
             ncfilename: /path/to/remapper/inputs/MOM.res_1.nc
             grid_stagger: mass
             interp_type: bilinear
             xdim_name: lonh
             ydim_name: lath

	ucurr:
	     src_ncvarname: u
             dst_ncvarname: u
             ncfilename: /path/to/remapper/inputs/MOM.res.nc
             grid_stagger: uvel
             interp_type: bilinear
             xdim_name: lonq
             ydim_name: lath
             zdim_name: Layer

	vcurr:
             src_ncvarname: v
             dst_ncvarname: v
             ncfilename: /path/to/remapper/inputs/MOM.res_1.nc
             grid_stagger: vvel
             interp_type: bilinear
             xdim_name: lonh
             ydim_name: latq
             zdim_name: Layer



+------------------------------+---------------------------------------------------------------------------+
| Configuration Variable       | Description                                                               |
+==============================+===========================================================================+
| ``bathy_edits_file``         | This is an optional variable specifying the file path to a                |
|                              | netCDF-formatted file containing bathymetry value edits; if not being     |
|                              | used this should be set to ``Null`` or commented out.                     |
+------------------------------+---------------------------------------------------------------------------+
| ``bathy_file``               | The file path to the netCDF formatted file containing the destination     |
|                              | grid projection bathymetry.                                               |
+------------------------------+---------------------------------------------------------------------------+
| ``forecast_model``           | The forecast model from which the analysis was produced; currently        |
|                              | supported are MOM6, `CICE <https://tinyurl.com/cice-model>`_,             |
|                              | `ORAS5 <https://tinyurl.com/oras5-sst>`_, and                             |
|                              | `GLORYS <https://tinyurl.com/glorys-seaice>`_ sea-ice.                    |
+------------------------------+---------------------------------------------------------------------------+
| ``grids``                    | The respective ``source`` and ``destination`` grid informations; see the  |
|                              | the corresponding table below for further information.                    | 
+------------------------------+---------------------------------------------------------------------------+
| ``interp_scheme``            | The interpolation scheme to be used for the remapping; currently          |
|                              | supported schemes are ESMF and SLINT; **Note that the files provided      |
|                              | beneath the YAML** ``remap`` **key must be generated from the same        |
|                              | interpolation scheme.**                                                   |
+------------------------------+---------------------------------------------------------------------------+
| ``nlevs``                    | The total number of vertical levels to be remapped to the destination     |
|                              | grid projection.                                                          |
+------------------------------+---------------------------------------------------------------------------+
| ``output_netcdf``            | The netCDF-formatted file path to the remapped (MOM6) analysis variables. |
+------------------------------+---------------------------------------------------------------------------+
| ``remap``                    | The remapping attributes; see the table below for further information.    |
+------------------------------+---------------------------------------------------------------------------+
| ``variables``                | The analysis variable remapping attributes; see the table below for more  |
|                              | information.                                                              |
+------------------------------+---------------------------------------------------------------------------+

+----------------------------------+---------------------------------------------------------------------------+
| ``grids`` Configuration Variable | Description                                                               |
+==================================+===========================================================================+
| ``arakawa``                      | The Arakawa type for the respective grid (e.g., destination or source).   |
+----------------------------------+---------------------------------------------------------------------------+
| ``grid_ncfile``                  | The path to the netCDF-formatted file containing the respective grid      |
|                                  | geographical locations.                                                   |
+----------------------------------+---------------------------------------------------------------------------+
| ``mass``, ``uvel``, ``vvel``     | The respective grid coordinate types for which variables may be defined   |
|                                  | on an Arakawa type grid; the ``nclat`` and ``nclon`` variables are the    |
|                                  | netCDF-formatted file variable names for the respective latitude and      |
|                                  | longitude geographical locations.                                         |
+----------------------------------+---------------------------------------------------------------------------+
| ``topo_ncfile``                  | The path to the netCDF-formatted file containing the respective grid      |
|                                  | topography and corresponding landmask.                                    |
+----------------------------------+---------------------------------------------------------------------------+

+----------------------------------+---------------------------------------------------------------------------+
| ``remap`` Configuration Variable | Description                                                               |
+==================================+===========================================================================+
| ``dstuvel2dstmass_bilinear``     | Bilinear remapping coefficients from the destination ``uvel`` grid        |
|                                  | locations to the destination grid ``mass`` grid locations.                |
+----------------------------------+---------------------------------------------------------------------------+
| ``dstuvel2dstmass_nrstnghbr``    | Nearest-neighbor remapping coefficients from the destination ``uvel``     |
|                                  | grid locations to the destination grid ``mass`` grid locations.           |
+----------------------------------+---------------------------------------------------------------------------+
| ``dstvvel2dstmass_bilinear``     | Bilinear remapping coefficients from the destination ``vvel`` grid        |
|                                  | locations to the destination grid ``mass`` grid locations.                |
+----------------------------------+---------------------------------------------------------------------------+
| ``dstvvel2dstmass_nrstnghbr``    | Nearest-neighbor remapping coefficients from the destination ``vvel``     |
|                                  | grid locations to the destination grid ``mass`` grid locations.           |
+----------------------------------+---------------------------------------------------------------------------+
| ``dstmass2dstuvel_bilinear``     | Bilinear remapping coefficients from the destination ``mass`` grid        |
|                                  | locations to the destination grid ``uvel`` grid locations.                |
+----------------------------------+---------------------------------------------------------------------------+
| ``dstmass2dstuvel_nrstnghbr``    | Nearest-neighbor remapping coefficients from the destination ``mass``     |
|                                  | grid locations to the destination grid ``uvel`` grid locations.           |
+----------------------------------+---------------------------------------------------------------------------+
| ``dstmass2dstvvel_bilinear``     | Bilinear remapping coefficients from the destination ``mass`` grid        |
|                                  | locations to the destination grid ``vvel`` grid locations.                |
+----------------------------------+---------------------------------------------------------------------------+
| ``dstmass2dstvvel_nrstnghbr``    | Nearest-neighbor remapping coefficients from the destination ``mass``     |
|                                  | grid locations to the destination grid ``vvel`` grid locations.           |
+----------------------------------+---------------------------------------------------------------------------+
| ``rotate_currents``              | A boolean variable specifying whether the current vector components must  |
|                                  | be rotated following remapping; this is ``True`` for all staggered        |
|                                  | Arakawa grids.                                                            |
+----------------------------------+---------------------------------------------------------------------------+
| ``srcmass2dstmass_bilinear``     | Bilinear remapping coefficients from the source ``mass`` grid locations   |
|                                  | to the destination grid ``mass`` grid locations.                          |
+----------------------------------+---------------------------------------------------------------------------+
| ``srcmass2dstmass_nrstnghbr``    | Nearest-neighbor remapping coefficients from the source ``mass`` grid     |
|                                  | locations to the destination grid ``mass`` grid locations.                |
+----------------------------------+---------------------------------------------------------------------------+
| ``srcmass2dstuvel_bilinear``     | Bilinear remapping coefficients from the source ``mass`` grid locations   |
|                                  | to the destination grid ``uvel`` grid locations.                          |
+----------------------------------+---------------------------------------------------------------------------+
| ``srcmass2dstuvel_nrstnghbr``    | Nearest-neighbor remapping coefficients from the source ``mass`` grid     |
|                                  | locations to the destination grid ``uvel`` grid locations.                |
+----------------------------------+---------------------------------------------------------------------------+
| ``srcmass2dstvvel_bilinear``     | Bilinear remapping coefficients from the source ``mass`` grid locations   |
|                                  | to the destination grid ``vvel`` grid locations.                          |
+----------------------------------+---------------------------------------------------------------------------+
| ``srcmass2dstvvel_nrstnghbr``    | Nearest-neighbor remapping coefficeitns from the source ``mass`` grid     |
|                                  | locations to the destination grid ``vvel`` grid locations.                |
+----------------------------------+---------------------------------------------------------------------------+
| ``srcuvel2dstmass_bilinear``     | Bilinear remapping coefficeitns from the source ``mass`` grid locations   |
|                                  | to the destination grid ``mass`` grid locations.                          |
+----------------------------------+---------------------------------------------------------------------------+
| ``srcuvel2dstmass_nrstnghbr``    | Nearest-neighbor remapping coefficients from the source ``mass`` grid     |
|                                  | locations to the destination grid ``mass`` grid locations.                |
+----------------------------------+---------------------------------------------------------------------------+
| ``srcuvel2srcmass_bilinear``     | Bilinear remapping coefficients from the source ``uvel`` grid             |
|                                  | locations to the source grid ``mass`` grid locations.                     |
+----------------------------------+---------------------------------------------------------------------------+
| ``srcuvel2srcmass_nrstnghbr``    | Nearest-neighbor remapping coefficients from the source ``uvel``          |
|                                  | grid locations to the source grid ``mass`` grid locations.                |
+----------------------------------+---------------------------------------------------------------------------+
| ``srcvvel2srcmass_bilinear``     | Bilinear remapping coefficients from the source ``vvel`` grid             |   
|                                  | locations to the source grid ``mass`` grid locations.                     |
+----------------------------------+---------------------------------------------------------------------------+
| ``srcvvel2srcmass_nrstnghbr``    | Nearest-neighbor remapping coefficients from the source ``vvel``          |  
|                                  | grid locations to the source grid ``mass`` grid locations.                |  
+----------------------------------+---------------------------------------------------------------------------+

The configuration keys for the respective ``variables`` are as follows.

+------------------------------+---------------------------------------------------------------------------+
| Configuration Variable       | Description                                                               |
+==============================+===========================================================================+
| ``dst_ncvarname``            | The destination grid netCDF variable name.                                |
+------------------------------+---------------------------------------------------------------------------+
| ``grid_stagger``             | The destination grid grid-staggered locations (e.g., ``mass``, ``uvel``,  |
|                              | or ``vvel``).                                                             |
+------------------------------+---------------------------------------------------------------------------+
| ``interp_type``              | The interpolation type to be used for remapping; may be either            |
|                              | ``bilinear`` or ``nrstnghbr``.                                            |
+------------------------------+---------------------------------------------------------------------------+
| ``ncfilename``               | The netCDF-formatted file path containing the respective source grid      |
|                              | variable values.                                                          |
+------------------------------+---------------------------------------------------------------------------+
| ``src_ncvarname``            | The source grid netCDF variable name.                                     |
+------------------------------+---------------------------------------------------------------------------+
| ``xdim_name``                | The netCDF x-dimension coordinate variable name for the destination grid  |
|                              | (i.e., remapped) varible.                                                 |
+------------------------------+---------------------------------------------------------------------------+
| ``ydim_name``                | The netCDF y-dimension coordinate variable name for the destination grid  |
|                              | (i.e., remapped) varible.                                                 |
+------------------------------+---------------------------------------------------------------------------+
| ``zdim_name``                | The netCDF z-dimension coordinate variable name for the destination grid  |
|                              | (i.e., remapped) varible; this is only required for 3-dimensional         |
|                              | variables.                                                                |
+------------------------------+---------------------------------------------------------------------------+

Example configuration files for both MOM6 and CICE analysis remappings
can be found `here
<https://github.com/HenryWinterbottom-NOAA/ufs_tools/tree/develop/parm/remapper>`_.

^^^^^^^^^^^^^^^^^^^^^^^^^
Launching the Application
^^^^^^^^^^^^^^^^^^^^^^^^^

The ``remapper`` application may be launched within the UFS Tools
package for the MOM6 example described in the previous section as
follows.

::

   user@host:$ cd scripts
   user@host:$ python compute_remapper.py --yaml_file ../parm/remapper/remapper.yaml

If successful the netCDF-formatted file path defined by
``output_netcdf`` will be created. The following are selected images
illustrating the remapping of MOM6 and CICE :math:`0.25^{\circ}`
(left) resolution analyses to :math:`1.0^{\circ}` (center) and
:math:`5.0^{\circ}` (right) tripolar grid projections.

.. image:: _images/MOM6.sst.0p25.png
   :width: 33%

.. image:: _images/MOM6.sst.1p0.png
   :width: 33%

.. image:: _images/MOM6.sst.5p0.png
   :width: 33%

.. image:: _images/CICE.icefrac.0p25.npstere.png
   :width: 33%

.. image:: _images/CICE.icefrac.1p0.npstere.png
   :width: 33%

.. image:: _images/CICE.icefrac.5p0.npstere.png
   :width: 33%

.. image:: _images/CICE.icefrac.0p25.spstere.png
   :width: 33%

.. image:: _images/CICE.icefrac.1p0.spstere.png
   :width: 33%

.. image:: _images/CICE.icefrac.5p0.spstere.png
   :width: 33%
