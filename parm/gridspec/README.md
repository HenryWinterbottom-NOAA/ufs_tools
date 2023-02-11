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


</div>
