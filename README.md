# pypaftol: Home of `paftools`

This repository contains the `paftools` Python module which provides functionality for recovery and data analysis of target capture data in the [PAFTOL Project][https://www.kew.org/science/our-science/projects/plant-and-fungal-trees-of-life]. Some of this functionality will be applicable and hopefully found to be useful beyond PAFTOL, e.g. for processing HybSeq data in general. Functionality can be accessed either via the `paftools` script or the Python API.

The documentation provided here is intended to help beta-testers with getting started.

Please also check out the [tutorial](paftools_tutorial.md) here.


## Prerequisites

Building and installing the `paftol` module requires

* Python 2.7.x (Python 3.x is currently not supported)
* Python 2.7.x development package (`libpython-all-dev`)
* BioPython 1.66 (newer versions are likely to work as well)  (`python-biopython`, and also `python-biopython-sql`)
* setuptools (`python-setuptools`)
* epydoc (`python-epydoc`)
* GNU C compiler (`gcc`) and associated tools
* GNU make

The following bioinformatics applications and suites are required for full functionality of the module and the `paftools` script:

* Trimmomatic - to use Trimmomatic via Paftools a little shell script is required called ```trimmomatic``` that needs to be available from the command line
```
#! /bin/bash
args=$@
java -jar <FULL_PATH_TO>/trimmomatic-0.39.jar ${args[@]}
```
* spades
* samtools
* bwa
* exonerate
* mafft
* clustalo (aka clustal-omega)
* emboss
* embassy-phylip
* fastqc (currently exactly version 0.11.5 is required)

Additional prerequisites for PAFTOL internal use include:

* Python `mysql.connector`  <!--Paul B. - changed from mysqlc.connector -->

These prerequisites should generally be provided on the cluster.


## Installation Guide

1. Clone the repository
```shell
git clone https://github.com/RBGKew/pypaftol
```

2. Install by running the command
```
make hinstall
```
This will install the package in `$HOME/lib/python`, which is the standard directory for installing Python modules for use in your
account only. You'll need to ensure that your `PYTHONPATH` environment variable includes this directory, see Tips section below.

3. Check that the installation was successful by running
```
paftools -h
```
This should give you a help message listing the `paftools` subcommands
currently available.

4. If you like a HTML version of the APIs provided by the `paftools` package
and its subpackages, run the command
```
make doc
```

At the time of writing this `README`, this installation process works on the cluster. Sharing any feedback is very welcome, of course.


## Tools Using the PAFTOL Databases

There are two PAFTOL databases:

* The **production** database is used to keep track of the data
  production process, from sourcing specimens to the generation of
  `fastq` files. It is developed and operated by the PAFTOL Data
  Production Team, in collaboration with IT's Science Application
  Development Team.

* The **analysis** database is designed to record data analysis,
  currently focusing on recording recovery of target coding sequences.
  Additional uses also include recording of external accession
  numbers, e.g. from the ENA / SRA, both upon submitting PAFTOL data
  to external archives and in the context of adding data (including
  NGS data, assembled transcriptomes, and others), to PAFTOL's
  phylogeny computations.

There are `paftools` subcommands to interact with these databases
programmatically, which requires database access details. All
`paftool` tools expect these details in the `~/.paftol` directory. The
details for accessing the production database and the analysis
database are expected in the files `productiondb.cfg` and
`analysisdb.cfg`, respectively. These files must have the following
structure:

```
username: <username>
password: <password>
host: <hostname>
dbname: <databasename>
```

The elements in angled braces must be replaced with the actual details
which can be obtained from the Data Production and Data Analysis
Teams, respectively.

**N.B.:** As these files contain database passwords they should
/never/ be attached to emails or included in git repositories.


## Docker Image

The `paftol-base` Docker image has pypaftol and prerequisites
installed, with `paftools` ready to be used. If you have Docker
installed on your computer, you should be able to pull this image and
then use it to run a container using the commands:
```
docker pull jttkim/paftol-base:latest
docker run -t -i jttkim/paftol-base:latest
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


## Developer Notes

### Adding a `paftools` command

To add a command, `foo` in this example:

* identify the (mandatory) parameters and options required by `foo`
* write a function `addFooParser` that takes an argparse parser as an argument and adds the relevant parameters and options to that
* write a `runFoo` function that takes a single `argNamespace` parameter, containing the argument namespace generated by the parser, and uses the attributes in that namespace to execute the command
* finally, wire everything up by calling `addFooParser` in `paftoolsMain` and by calling `p.add_default(func=runFoo)` on the subparser
