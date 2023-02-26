# Configuring the EMSF Remapping Application

This document describes the configuration options for the supported
`esmf_remap` applications. The following is a snippet for a
YAML-formatted `esmf_remap` configuration for the
`scripts/compute_esmf_remap.py` application.

~~~
grids:

     destination:
          is_grib: False
          gribfile: Null
          is_netcdf: True
          ncfile: /run/gridspec.tripolar.0p25.nc        
          nclat: tlat
          nclon: tlon
          ncxdim_name: nx
          ncydim_name: ny

     source:
          is_grib: False
          gribfile: Null
          is_netcdf: True
          ncfile: /run/gridspec.tripolar.5p0.nc      
          nclat: tlat
          nclon: tlon
          ncxdim_name: nx
          ncydim_name: ny    

interp_type: bilinear
output_netCDF: /run/tripolar.0p25_5p0.bilinear.t2t_remap.nc
~~~

The attributes `destination` and `source` denote the respective grid
types[^1]. The following table describes the YAML-block attributes
common to both grid types. 

[^1]: Source grid implies that the variables to be remapped will be
defined on the source grid. Destination grid implies that the source
variables will be remapped to the respective grid using the remapping
coefficients defined within the corresponding netCDF-formatted output
file (see `output_netCDF`).

<div align="center">

| Variable | Description |
| :-------------: | :-----------: |
| `is_grib` | <div align="left">A boolean valued variable defining whether the respective grid's projection is defined (and to be collected from) a [WMO-GRIB](https://community.wmo.int/en/activity-areas/wis/grib3) formatted file (not yet supported).</div> | 
| `gribfile` | <div align="left">The path to the WMO-GRIB formatted file containing the respective grid's projection.</div> | 
| `is_netcdf` | <div align="left">A boolean valued variable defining whether the respective grid's projection is defined (and to be collected from) a [netCDF](https://www.unidata.ucar.edu/software/netcdf/)-formatted file.</div> |
| `ncfile` | <div align="left">The path to the netCDF-formatted file containing the respective grid's projection.</div> | 
| `nclat` | <div align="left">The netCDF-formatted file variable name for the latitude geographical location variable array.</div> | 
| `nclon` | <div align="left">The netCDF-formatted file variable name for the longitude geographical location variable array.</div> | 
| `ncxdim_name` | <div align="left">The x-dimension for the respective geographical location variable arrays (see `nclat` and `nclon`).</div> |
| `ncydim_name` | <div align="left">The y-dimension for the respective geographical location variable arrays (see `nclat` and `nclon`).</div> |

</div>

The `interp_type` specifies the type of
remapping (i.e., how the ESMF coeffients are defined/computed). The
currently supported options are `bilinear`, `nearest_d2s`, and
`nearest_s2d`.

The above-mentioned wrapper application,
`wrappers/wrapper_esmf_remap.py` allows the user to execute several
sequential instances of the ESMF remapping application. The following
code snippet describes the required YAML-formatted configuration file.

~~~
tripolar.0p25_5p0.bilinear.t2t_remap:

     destination:
          is_grib: False
          gribfile: Null
          is_netcdf: True
          ncfile: /run/gridspec.tripolar.5p0.nc
          nclat: tlat
          nclon: tlon
          ncxdim_name: nx
          ncydim_name: ny     
     
     source:

          is_grib: False
          gribfile: Null
          is_netcdf: True
          ncfile: /run/gridspec.tripolar.0p25.nc
          nclat: tlat
          nclon: tlon
          ncxdim_name: nx
          ncydim_name: ny

     interp_type: bilinear
     output_netCDF: /run/tripolar.0p25_5p0.bilinear.t2t_remap.nc
          
tripolar.0p25_5p0.nearest_s2d.t2t_remap:

     destination:

          is_grib: False
          gribfile: Null
          is_netcdf: True
          ncfile: /run/gridspec.tripolar.5p0.nc
          nclat: tlat
          nclon: tlon
          ncxdim_name: nx
          ncydim_name: ny     

     source:

          is_grib: False
          gribfile: Null
          is_netcdf: True
          ncfile: /run/gridspec.tripolar.0p25.nc
          nclat: tlat
          nclon: tlon
          ncxdim_name: nx
          ncydim_name: ny

     interp_type: nearest_s2d
     output_netCDF: /run/tripolar.0p25_5p0.nearest_s2d.t2t_remap.nc
~~~

In the above example, YAML keys define a unique name for the remapping
instance. The remainder of the variables within the respective block
are the same as in the table above. Finally, a `yaml_template`
attribute example is provided in `esmf_remap_tmpl.yaml`. Note that
this is a file identical to that used for the `compute_esmf_remap.py`
application described above. It should not be modified by the user.

#

Please direct questions to [Henry
R. Winterbottom](mailto:henry.winterbottom@noaa.gov?subject=[ufs_tools])