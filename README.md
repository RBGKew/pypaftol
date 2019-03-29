# pypaftol: Home of `paftools`

This repository contains the `paftol` Python module which provides
functionality for data analysis in the PAFTOL project. Some of this
functionality will be applicable and hopefully found to be useful
beyond PAFTOL, e.g. for processing HybSeq data in general.
Functionality can be accessed either via the `paftools` script or the
Python API.

The documentation provided here is intended to help beta-testers with
getting started. It is not complete and any feedback is welcome (and
will be rewarded with biscuits).


## Prerequisites

Building and installing the `paftol` module requires

* Python 2.7.x (Python 3.x is currently not supported)
* Python 2.7.x development package
* BioPython 1.66 (newer versions are likely to work as well)
* setuptools
* GNU C compiler (`gcc`) and associated tools
* GNU make

The following bioinformatics applications and suites are required for
full functionality of the module and the `paftools` script:

* samtools
* bwa
* exonerate
* mafft
* clustalo (aka clustal-omega)
* emboss
* embassy-phylip
* spades
* fastqc (currently exactly version 0.11.5 is required)

Additional prerequisites for PAFTOL internal use include:

* Python `mysqlc.connector`

These prerequisites should generally be provided on the cluster.


## Quick Installation Guide

1. Clone the repository
```
git clone https://github.com/RBGKew/pypaftol.git
```

2. Install by running the command
```
make hinstall
```
This will install the package in `$HOME/lib/python`, which is the
standard directory for installing Python modules for use in your
account only. You'll need to ensure that your `PYTHONPATH` environment
variable includes this directory, see Tips section below.

3. Check that the installation was successful by running
```
paftools -h
```
This should give you a help message listing the `paftools` subcommands
currently available.

At the time of writing this `README`, this installation process works
on the cluster. Sharing any feedback is very welcome, of course.


## Docker Image

The `paftol-base` Docker image has pypaftol and prerequisites
installed, with `paftools` ready to be used. If you have Docker
installed on your computer, you should be able to pull this image and
then use it to run a container using the commands:
```
docker pull jttkim/paftol-base
docker run -t -i jttkim/paftol-base
```

This image is primarily for use by the PAFTOL Data Analysis team.
Therefore, the default version (`latest`) may change without prior
announcements. You can append a version number to the image name in
the `docker pull` and `docker run` commands if you prefer to work with
the same image regardless of any updates, please see
https://hub.docker.com/r/jttkim/paftol-base for a list of available
versions.

For using the Docker image for developing your own code using
`paftools` or the `paftol` Python module I suggest you create a
directory for your work on the system you use as your Docker host
(e.g. your laptop) and bind mount that directory on the container, see
https://docs.docker.com/engine/reference/commandline/run/ .


## Tips

### Setting up `PYTHONPATH`

`PYTHONPATH` is an environment variable which the Python interpreter
uses to obtain a list of directories to search for modules when
executing an `import` statement. By default, this variable won't
include any directories in your login directory, so if you want to
install any modules in your personal space, you'll need to add the
directory where you install modules for your personal use. This can be
done by the following snippet of bash code:
```
if test -z "$PYTHONPATH" ; then
  PYTHONPATH=${HOME}/lib/python
else
  PYTHONPATH="${HOME}/lib/python:${PYTHONPATH}"
fi
export PYTHONPATH
```
Please review / scrutinise this code by consulting the bash and Python
documentation and adapt it as necessary. Once you're satisfied, you
may consider adapting your login script (e.g. `~/.profile`)
accordingly, so `PYTHONPATH` is set up automatically each time you log
in.
