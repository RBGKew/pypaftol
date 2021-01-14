# pypaftol: Home of `paftools`

This repository contains the `paftools` Python module which provides functionality for recovery and data analysis of target capture data in the [PAFTOL Project]([https://www.kew.org/science/our-science/projects/plant-and-fungal-trees-of-life). Some of this functionality will be applicable and hopefully found to be useful beyond PAFTOL, e.g. for processing HybSeq data in general. Functionality can be accessed either via the `paftools` script or the Python API.

The documentation provided here is intended to help beta-testers with getting started.

Description and examples of usage are provided in our [usage](paftools_usage.md) and [tutorial](paftools_tutorial.md) readme.


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

* trimmomatic
* spades
* samtools
* bwa
* exonerate
* mafft
* clustalo (aka clustal-omega)
* emboss
* embassy-phylip
* fastqc (currently exactly version 0.11.5 is required)


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

Additional information about installation are available [here][Advanced_Install.md]

## Tips

### Setting up `PYTHONPATH`

`PYTHONPATH` is an environment variable which the Python interpreter uses to obtain a list of directories to search for modules when
executing an `import` statement. By default, this variable won't include any directories in your login directory, so if you want to
install any modules in your personal space, you'll need to add the directory where you install modules for your personal use. This can be
done by the following snippet of bash code:

```
if test -z "$PYTHONPATH" ; then
  PYTHONPATH=${HOME}/lib/python
else
  PYTHONPATH="${HOME}/lib/python:${PYTHONPATH}"
fi
export PYTHONPATH
```

Please review / scrutinise this code by consulting the bash and Python documentation and adapt it as necessary. Once you're satisfied, you may consider adapting your login script (e.g. `~/.profile`) accordingly, so `PYTHONPATH` is set up automatically each time you login.