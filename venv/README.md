# UFS Tools Python Virtual Environment

The following describes how to configure a Python virtual environment
to allow for the use of selected `ufs_tools` applications. Please note
that the remapper applications are not supported by this virtual
environment/stack. The Docker and/or Singularity containers should be
used instead for those applications (see [here](../Docker)).

### Building the Virtual Environment

The Python virtual environment may be configured and install as
follows.

~~~
user@host:$ /path/to/python -m venv /path/to/virtual_environment
user@host:$ cd /path/to/virtual_environment
user@host:$ ./bin/python -m pip install --upgrade pip
user@host:$ ./bin/pip install -r ../stack-requirements.txt
~~~

If successful this will yield the path `/path/to/virtual_environment`
containing the Python virtual environment. The Python virtual
environment can be installed as follows.

### Activating the Virtual Environment

The action below assumes a POSIX bash environment (e.g., shell). For
activation support for alternative shells, see
[here](https://tinyurl.com/activate-venv). 

~~~
user@host:$ . ./bin/activate
~~~

Supported `ufs_tools` applications may now be executed within the
Python virtual environment.

### Deactivating the Virtual Environment

When finished, exit the Python virtual environment by doing the
following.

~~~
user@host:$ deactivate
~~~

#

Please direct questions to [Henry
R. Winterbottom](mailto:henry.winterbottom@noaa.gov?subject=[UFS-Tools)