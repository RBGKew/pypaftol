#!/usr/bin/env python


import setuptools



clib = setuptools.Extension('paftol.clib',
                            sources = ['src/clib.c'],
                            include_dirs = [],
                            library_dirs = [],
                            libraries = [],
                            extra_compile_args = ['-Wall', '-pedantic', '-Wno-long-long', '-fPIC'],
			    extra_link_args = ['-fPIC'])
paftolScripts = ['paftools']
paftolDescription = 'packages and tools for the Plant and Fungal Trees of Life (PAFTOL) project'
paftolUrl = 'https://www.paftol.org/'
paftolEntryPoints = {'console_scripts': ['paftools = paftol.cli:paftoolsMain']}

setuptools.setup(name='paftol',
                 version='0.0.1',
                 packages=setuptools.find_packages(),
                 author='PAFTOL Bioinformatics Team',
                 author_email='info@paftol.org',
                 description=paftolDescription,
                 keywords='bioinformatics phylogeny',
                 url=paftolUrl,
                 entry_points = paftolEntryPoints,
                 ext_modules = [clib])
# could also include long_description, download_url, classifiers, etc.
