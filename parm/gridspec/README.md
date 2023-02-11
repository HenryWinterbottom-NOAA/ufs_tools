# Configuring the Gridspec Application

This document describes the configuration options for the `gridspec`
application. The following is a snippet of a YAML-formatted
configuration file for the `gridspec` application.

~~~
output_netcdf: /run/gridspec.tripolar.0p25.nc
grid_type: C

latitude:

     ncfile: /run/mom6_cice_grids/0p25/ocean_hgrid.nc
     ncxdim: nx
     ncydim: ny
     ncvarname: y

longitude:
     ncfile: /run/mom6_cice_grids/0p25/ocean_hgrid.nc
     ncxdim: nx
     ncydim: ny
     ncvarname: x
~~~