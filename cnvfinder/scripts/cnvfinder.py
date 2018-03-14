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
    module_group.add_argument('--module', type=str, nargs='?', default=Modules.CNVDETECTOR.value,
                              choices={Modules.BEDLOADER.value, Modules.CNVDETECTOR.value,
                                       Modules.NRRHANDLER.value, Modules.VCFHANDLER.value},
                              help='Module to run')

    # Modules.BEDLOADER
    bedloader_group = parser.add_argument_group(Modules.BEDLOADER.value,
                                                'Module for preprocessing amplicons defined in a BED file.')
    bedloader_group.add_argument('--target', type=str, required=True,
                                 help='Path to file in which sequencing amplicons are listed.')
    bedloader_group.add_argument('--region', type=str,
                                 help='Limit target definition to a given region. It should be in the form: '
                                      'chr1:10000-90000')
    bedloader_group.add_argument('--spacing', type=int, default=20,
                                 help='Number of nucleotides to ignore at amplicon start and end, to avoid '
                                      'overlapping reads.')
    bedloader_group.add_argument('--min-data', type=int, default=50,
                                 help='Minimum number of nucleotides for a valid target.')
    bedloader_group.add_argument('--max-pool', type=int, default=None,
                                 help='Maximum number of origin pools allowed for a valid target.')

    # Modules.CNVDETECTOR
    cnvdetector_group = parser.add_argument_group(Modules.CNVDETECTOR.value,
                                                  '{} parameters'.format(Modules.CNVDETECTOR.value))

    # Modules.NRRHANDLER
    nrrhandler_group = parser.add_argument_group(Modules.NRRHANDLER.value,
                                                 '{} parameters'.format(Modules.NRRHANDLER.value))

    # Modules.VCFHANDLER
    vcfhandler_group = parser.add_argument_group(Modules.VCFHANDLER.value,
                                                 '{} parameters'.format(Modules.VCFHANDLER.value))

    opts = parser.print_help()


def main():
    parsecmd()


if __name__ == '__main__':
    main()
