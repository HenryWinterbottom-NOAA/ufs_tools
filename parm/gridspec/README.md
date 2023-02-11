# Configuring the Gridspec Application

This document describes the configuration options for the `gridspec`
application. The following is a snippet of a YAML-formatted
configuration file for the `gridspec` application.

~~~
grid_type: C
is_tripolar: True
is_wrap_lons: True

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

output_netcdf: /run/gridspec.tripolar.0p25.nc

~~~

The following table describes each of the respective variables above.

<div align="center">

| Variable | Description |
| :-------------: | :-----------: |
| `grid_type` | <div align="left">This variable value specifies the type of Arakawa-type grids to define/compute; a comprehensive list of Arakawa-type grids can be found [here](https://en.wikipedia.org/wiki/Arakawa_grids); currently only Arakawa-C type grids are supported.</div> |
| `is_tripolar` | <div align="left">This variable value specifies whether the output gridspec-formatted file is a tri-polar grid projection, a discussion of which can be found [here](https://github.com/dtcenter/MET/issues/1231); for Arakawa-C type tripolar grid projections the dimensions of the zonal- and meridional-velocity grid geographical locations must be the same as the mass grids geographical location array dimensions; this requires that the arrays are define accordingly.</div> |
| `is_wrap_lons` | <div align="left">If `True`, the longitude coordinate values will be reset within the range $[0, 2\pi)$.</div>|
| `output_netcdf` | <div align="left">The path to the netCDF-formatted gridspec-formatted output file; note that if this application is run within the provided Docker (or Singularity) container, the directory-tree component of the path must be a path bound to the container when it is launched.</div> | 

</div>

The `latitude` and `longitude` YAML keys within the example
configuration define the attributes for the respective geographical
coordinates. Seperate YAML blocks are implemented to support instances
when the respective geographical coordinate values are collected from
different file path sources. The following section describes the
YAML-keys within the YAML bloeks.

<div align="center">

| Variable | Description |
| :-------------: | :-----------: |
| `ncfile` | <div align="left">The path to the [netCDF](https://www.unidata.ucar.edu/software/netcdf/)-formatted filepath containing the respective geographical coordinate variable values; note that if this application is run within the provided Docker (or Singularity) container, the directory-tree component of the path must be a path bound to the container when it is launched.</div> |

</div>