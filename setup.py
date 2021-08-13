#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: valengo
"""

# always prefer setuptools over distutils
from setuptools import setup, find_packages
# to use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))
name = 'amplicnv'

# get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name=name,

    # versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version='1.0.0b6',

    description='A Python3.x package for CNV detection',
    long_description=long_description,

    # the project's main homepage.
    url='https://github.com/ip4-team/amplicnv',

    # author details
    author='IP4 Team',
    author_email='ip4.developers@gmail.com',

    # is this ok?
    license='MIT',

    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    # What does your project relate to?
    keywords='CNV analysis development',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=['img', 'test', 'poster', 'venv']),

    # Alternatively, if you want to distribute just a my_module.py, uncomment
    # this:
    # py_modules=["nrrhandler", "vcfhandler", "amplicnv", "bedloader"],

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=['pysam', 'pandas', 'plotly', 'numpy', 'scikit-learn',
                      'scipy', 'numexpr', 'matplotlib', 'bedhandler'],

    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
    # extras_require={
    #    'dev': ['check-manifest'],
    #    'test': ['coverage'],
    # },

    python_requires='>=3',

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    package_data={
       name: ['data/*.txt'],
    },

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    # data_files=[('my_data', ['data/data_file'])],

    # To provide executable wrapper, use entry points in preference to the
    # "wrapper" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
       'console_scripts': [
           'amplicnv=amplicnv.wrapper.command_line:main'
       ]
    }
)
