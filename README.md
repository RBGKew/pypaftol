# pypaftol: Home of `paftools`

This repository contains the `paftol` Python module which provides
functionality for data analysis in the PAFTOL project. Some of this
functionality will be applicable and hopefully found to be useful
beyond PAFTOL, e.g. for processing HybSeq data in general.
Functionality can be accessed either via the `paftools` script or the
Python API.


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

2. Ensure that your `PYTHONPATH` environment variable includes the
standard module path for the `--home` installation scheme of setuptools
(see Tips section below).

3. Install by running the command
```
make hinstall
```

4. Check that the installation was successful by running
```
paftools -h
```
This should give you a help message listing the `paftools` subcommands
currently available.


## Tips

### Setting up `PYTHONPATH`

This can be done by the following snippet of bash code:
```
if test -z "$PYTHONPATH" ; then
  PYTHONPATH=${HOME}/lib/python
else
  PYTHONPATH="${HOME}/lib/python:${PYTHONPATH}"
fi
export PYTHONPATH
```

Please review / scrutinise this code and consult the Python
documentation as necessary. Once you're informed and satisfied, you
may consider adapting your login script (e.g. `~/.profile`)
accordingly, so `PYTHONPATH` is set up automatically each time you log
in.
