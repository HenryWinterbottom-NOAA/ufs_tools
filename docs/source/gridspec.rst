#######################################
Building Grid Projection Files
#######################################

^^^^^^^^^^^^^^^^^^^^^^^^^
Application Configuration
^^^^^^^^^^^^^^^^^^^^^^^^^

The ``gridspec`` application is configured using YAML-formatted
files as follows.

::

   grid_type: C
   is_tripolar: True
   is_wrap_lons: True
     
   # Define the grid-projection angle orientation information attributes.
   angle:

       ncfile: /work/gridspec/mom6_cice_grids.gridspec/1p0/ocean_hgrid.nc
       ncvarname: angle_dx
       ncxdim: nxp
       ncydim: nyp

   # Define the latitude information attributes.
   latitude:

       ncfile: /work/gridspec/mom6_cice_grids.gridspec/1p0/ocean_hgrid.nc
       ncvarname: y
       ncxdim: nx
       ncydim: ny
     
   # Define the longitude information attributes.
   longitude:

       ncfile: /work/gridspec/mom6_cice_grids.gridspec/1p0/ocean_hgrid.nc
       ncvarname: x
       ncxdim: nx
       ncydim: ny

    output_netcdf: /run/gridspec.tripolar.1p0.nc

+------------------------------+---------------------------------------------------------------------------+
| Configuration Variable       | Description                                                               |
+==============================+===========================================================================+
| ``angle``                    | The input grid-projection angle orientation; this is used for instances   |
|                              | of rotating velocity vectors for different grid projections -- namely     |
|                              | tripolar type.                                                            |
+------------------------------+---------------------------------------------------------------------------+
| ``gridtype``                 | The Arakawa-type grid; a complete list of Arakawa-type grids can be found |
|                              | `here <https://tinyurl.com/arakawa-type-grids>`_.                         |
+------------------------------+---------------------------------------------------------------------------+
| ``is_tripolar``              | Indicates whether the respective input grid projection is that of a       |
|                              | tripolar grid.                                                            |
+------------------------------+---------------------------------------------------------------------------+
| ``is_wraplons``              | Indicates whether to wrap the respective longitude array values; if True, |
|                              | this will also shift the grid-projection angle by a factor of :math:`\pi`.|
+------------------------------+---------------------------------------------------------------------------+
| ``latitude``                 | The input grid-projection latitude coordinate attributes.                 |
+------------------------------+---------------------------------------------------------------------------+
| ``longitude``                | The input grid-projection longitude coordinate attributes.                |
+------------------------------+---------------------------------------------------------------------------+
| ``output_netcdf``            | The netCDF-formatted file path to contain the computed grid-projection.   |
+------------------------------+---------------------------------------------------------------------------+

For the respective grid-projection keys (i.e., ``angle``,
``latitude``, and ``longitude``), the corresponding attributes are as
follow.

+------------------------------+---------------------------------------------------------------------------+
| Configuration Variable       | Description                                                               |
+==============================+===========================================================================+
| ``ncfile``                   | The path to the netCDF-formatted file containing the                      |
|                              | `supergrid                                                                |
|                              | <https://mom6.readthedocs.io/en/main/api/generated/pages/                 |
|                              | Discrete_Grids.html>`_ for the respective grid-projection.                |
+------------------------------+---------------------------------------------------------------------------+
| ``ncvarname``                | The netCDF-formatted file variable name for the respective geographical   |
|                              | coordinate/attribute.                                                     |
+------------------------------+---------------------------------------------------------------------------+
| ``ncxdim``                   | The zonal-coordinate (i.e., :math:`x`) netCDF dimension name.             |
+------------------------------+---------------------------------------------------------------------------+
| ``ncydim``                   | The meridional-coordinate (i.e., :math:`y`) netCDF dimension name.        |
+------------------------------+---------------------------------------------------------------------------+

Using the a :math:`1.0^{\circ}` `MOM6
<https://mom6.readthedocs.io/en/main/index.html>`_ tripolar-projection
supergrid, the netCDF attributes listed in the table above may be
found as follows.

::
   
   user@host:$ ncdump -h /work/gridspec/mom6_cice_grids.gridspec/1p0/ocean_hgrid.nc

   netcdf ocean_hgrid {
   dimensions:
           nyp = 641 ;
	   nxp = 721 ;
	   ny = 640 ;
	   nx = 720 ;
	   string = 255 ;
   variables:
	   char tile(string) ;
	   double y(nyp, nxp) ;
		   y:units = "degrees" ;
	   double x(nyp, nxp) ;
	           x:units = "degrees" ;
      	   double dy(ny, nxp) ;
		   dy:units = "meters" ;
	   double dx(nyp, nx) ;
		   dx:units = "meters" ;
	   double area(ny, nx) ;
	           area:units = "m2" ;
	   double angle_dx(nyp, nxp) ;
		   angle_dx:units = "degrees" ;
    }

Note that in the above example that the ``latitude`` and ``longitude``
grid-dimension attributes ``ncxdim`` and ``ncydim`` are set to ``nx``
and ``ny`` despite the netCDF-formatted file metadata attributes. This
is due to how the tripolar supergrid is reduced to compute the nominal
resolution ``gridspec`` grid.

^^^^^^^^^^^^^^^^^^^^^^^^^
Launching the Application
^^^^^^^^^^^^^^^^^^^^^^^^^

The ``gridspec`` application may be launched from within the UFS Tools
package for the MOM6 :math:`1.0^{\circ}` example above as follows.

::

   user@host:$ cd scripts
   user@host:$ python compute_gridspec.py --yaml_file ../parm/gridspec/gridspec.tripolar.1p0.yaml

If successful, a netCDF-formatted file containing the reduced-grid
``gridspec`` attributes will be created containing the following
metadata attributes.

::

   user@host:$ ncdump -h /run/gridspec.tripolar.1p0.nc

   netcdf gridspec.tripolar.1p0 {
   dimensions:
	   nx = 360 ;
	   ny = 320 ;
   variables:
	   double angle(ny, nx) ;
		   angle:description = "Grid projection rotation angle; radians." ;
	   double qlat(ny, nx) ;
		   qlat:description = "Array of q-grid latitudes; degrees north." ;
	   double qlon(ny, nx) ;
		   qlon:description = "Array of q-grid longitudes; degrees east." ;
	   double tlat(ny, nx) ;
		   tlat:description = "Array of t-grid latitudes; degrees north." ;
	   double tlon(ny, nx) ;
		   tlon:description = "Array of t-grid longtudes; degrees east." ;
	   double ulat(ny, nx) ;
		   ulat:description = "Array of u-grid latitudes; degrees north." ;
	   double ulon(ny, nx) ;
		   ulon:description = "Array of u-grid longitudes; degrees east." ;
	   double vlat(ny, nx) ;
		   vlat:description = "Array of v-grid latitudes; degrees north." ;
	   double vlon(ny, nx) ;
		   vlon:description = "Array of v-grid longitudes; degrees east." ;

   // global attributes:
		   :_FillValue = NaN ;
		   :date = "00:00:00 UTC 01 January 2000" ;
	           :created = "Ham E. Spam" ;
   }

The following table describes each of the respective metadata attributes.

+------------------------------+---------------------------------------------------------------------------+
| MetaData Variable(s)         | Description                                                               |
+==============================+===========================================================================+
| ``angle``                    | The grid projection rotation angle; this is used for non-Arakawa A grids  |
|                              | to rotate velocity vectors relative to the center of the respective grid  |
|                              | cell; units are degrees.                                                  |
+------------------------------+---------------------------------------------------------------------------+
| ``qlat``, ``qlon``           | The geographical coordinates for the respective grid cell vertices; units |
|                              | are degrees.                                                              |
+------------------------------+---------------------------------------------------------------------------+
| ``tlat``, ``tlon``           | The geographical coordinates for the respective grid cell centers; these  |
|                              | are often referred to as the mass-variable locations; units are degrees.  |
+------------------------------+---------------------------------------------------------------------------+
| ``ulat``, ``ulon``           | The geographical coordinates for the zonal velocity components for the    |
|                              | respective Arakawa-type grid; units are degrees.                          | 
+------------------------------+---------------------------------------------------------------------------+
| ``vlat``, ``vlon``           | The geographical coordinates for the meridional velocity components for   |
|                              | the respective Arakawa-type grid; units are degrees.                      | 
+------------------------------+---------------------------------------------------------------------------+

Example configuration files can be found beneath ``parm/gridspec``. 
