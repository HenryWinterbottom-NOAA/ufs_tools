# Applications

The following sections describe the supported applications. Users
should note that when running in a shell (i.e., not in the provided
Docker container), the following should be done before invoking any of
the provided applications.

~~~
user@host:$ export PYTHONPATH=/path/to/ufs_pyutils/:/path/to/ufs_tools/ush:$PYTHONPATH
~~~

When invoking any of the provided applications within the provided
Docker container all paths, including those within the respective
YAML-formatted configurtaion files, must be defined relative to the
Docker container directory tree. See the examples throughout for
additional detail.

## GridSpec

The `gridspec` application provides an interface for generating
[gridspec](https://arxiv.org/pdf/1911.08638.pdf) formatted
grids. Grids of this format are utilized by the various remapping
applications.

To launch the application in the user shell, do as follows.

~~~
user@host:$ cd scripts/
user@host:$ python compute_gridspec.py --yaml_file /path/to/gridspec/yaml
~~~

In the above example, the `yaml_file` attribute is a YAML-formatted
file containing the directives to define the gridspec format for the
specified grid-spacing resolution and projection. An example
YAML-formatted configuration file may be found
[here](./parm/gridspec/gridspec.yaml).

If choosing to launch the application within the available Docker
container, do as follows.

~~~
user@host:$ docker container run -v /path/to/ufs_tools:/ufs_tools -v /path/to/run:/run run -it ubuntu20.04-miniconda_ufs_pyutils.ufs_tools:latest
user@host:$ export PYTHONPATH=/ufs_tools/ush:$PYTHONPATH
user@host:$ cd /ufs_tools/scripts
user@host:$ python compute_gridspec.py --yaml_file=/path/to/gridspec/yaml
~~~

Note that the `yaml_file` attribute path is relative to the Docker
container directory tree.

#

Please direct questions to [Henry
R. Winterbottom](mailto:henry.winterbottom@noaa.gov?subject=[ufs_tools])