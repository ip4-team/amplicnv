#!/usr/bin/env python3
"""

@author: valengo
"""

import argparse
from enum import Enum


class Modules(Enum):
    BEDLOADER = 'bedloader'
    CNVDETECTOR = 'cnvdetector'
    NRRHANDLER = 'nrrhandler'
    VCFHANDLER = 'vcfhandler'


def parsecmd():
    parser = argparse.ArgumentParser(description='CNVfinder is a Python 3.x package for copy number (CNV) variation '
                                                 'detection on whole exome sequencing (WES) data from amplicon-based '
                                                 'enrichment technologies.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # module group
    module_group = parser.add_argument_group('module', 'Available modules')
    module_group.add_argument('-m', '--module', type=str, nargs='?', default=Modules.CNVDETECTOR,
                              choices={Modules.BEDLOADER, Modules.CNVDETECTOR, Modules.NRRHANDLER, Modules.VCFHANDLER},
                              help='Module to run')

    # Modules.BEDLOADER
    bedloader_group = parser.add_argument_group(Modules.BEDLOADER, '{} parameters'.format(Modules.BEDLOADER))

    # Modules.CNVDETECTOR
    cnvdetector_group = parser.add_argument_group(Modules.CNVDETECTOR, '{} parameters'.format(Modules.CNVDETECTOR))

    # Modules.NRRHANDLER
    nrrhandler_group = parser.add_argument_group(Modules.NRRHANDLER, '{} parameters'.format(Modules.NRRHANDLER))

    # Modules.VCFHANDLER
    vcfhandler_group = parser.add_argument_group(Modules.VCFHANDLER, '{} parameters'.format(Modules.VCFHANDLER))


def main():
    pass


if __name__ == '__main__':
    main()
