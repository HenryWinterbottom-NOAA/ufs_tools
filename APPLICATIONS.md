# Applications

The following sections describe the supported applications.

### Gridspec

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