#!/usr/bin/env python


import setuptools

paftolScripts = ['paftools']
paftolDescription = 'packages and tools for the Plant and Fungal Trees of Life (PAFTOL) project'
paftolUrl = 'https://www.paftol.org/'
paftolEntryPoints = {'console_scripts': ['paftools = paftol.cli:paftoolsMain']}

setuptools.setup(name='paftol', version='0.0.1', packages=setuptools.find_packages(), author='PAFTOL Bioinformatics Team', author_email='info@paftol.org', description=paftolDescription, keywords='bioinformatics phylogeny', url=paftolUrl, entry_points = paftolEntryPoints)
# could also include long_description, download_url, classifiers, etc.


