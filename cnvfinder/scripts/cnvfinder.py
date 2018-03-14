#!/usr/bin/env python3
"""

@author: valengo
"""

import argparse
import sys
from enum import Enum


class Modules(Enum):
    BEDLOADER = 'bedloader'
    CNVDETECTOR = 'cnvdetector'
    NRRHANDLER = 'nrrhandler'
    VCFHANDLER = 'vcfhandler'


def create_parser(description: str, usage: str = None) -> argparse.ArgumentParser:
    return argparse.ArgumentParser(description=description, usage=usage,
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)


class CNVFinder(object):
    def __init__(self):
        self.args = None
        parser = create_parser('CNVfinder is a Python 3.x package for copy number (CNV) '
                               'variation detection on whole exome sequencing (WES) data from '
                               'amplicon-based enrichment technologies',
                               usage='''cnvfinder <command> [<args>]
                               
Available commands used in various situations:

    {}
    \tPreprocess amplicons defined in a BED file
    
    {}
    \tDetect copy number variation on whole exome sequencing data
    \tCreate plots
    
    {}
    \tCompute read depth for each sequencing target
    \tCompare baseline samples x test sample (ratio)
    \tCreate read depth and ratio related plots
    
    {}
    \tFilter variants according to a baseline
    \tCompute B-allele frequency (BAF)
    \tCreate BAF related plots
    
For getting help of a specific command use: cnvfinder <command> --help
    '''.format(Modules.BEDLOADER.value, Modules.CNVDETECTOR.value, Modules.NRRHANDLER.value, Modules.VCFHANDLER.value))

        parser.add_argument('command', help='Module to run')
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            sys.exit('Unrecognized {} command'.format(args.command))

        # use dispatch pattern to invoke method with same name
        getattr(self, args.command)()

    def bedloader(self):
        parser = create_parser('Preprocess amplicons defined in a BED file',
                               usage='''cnvfinder bedloader [<args>]''')
        parser.add_argument('--target', type=str, required=True,
                            help='Path to file in which sequencing amplicons are listed')
        parser.add_argument('--region', type=str,
                            help='Limit target definition to a given region. It should be in the form: '
                                 'chr1:10000-90000')
        parser.add_argument('--spacing', type=int, default=20,
                            help='Number of nucleotides to ignore at amplicon start and end, to avoid '
                                 'overlapping reads')
        parser.add_argument('--min-data', type=int, default=50,
                            help='Minimum number of nucleotides for a valid target')
        parser.add_argument('--max-pool', type=int, default=None,
                            help='Maximum number of origin pools allowed for a valid target')

        self.args = parser.parse_args(sys.argv[2:])

    def cnvdetector(self):
        pass
        # Modules.CNVDETECTOR
        # cnvdetector_group = parser.add_argument_group(Modules.CNVDETECTOR.value,
        #
        #                                               '{} parameters'.format(Modules.CNVDETECTOR.value))

    def nrrhandler(self):
        pass
        # Modules.NRRHANDLER
        # nrrhandler_group = parser.add_argument_group(Modules.NRRHANDLER.value,
        #                                              '{} parameters'.format(Modules.NRRHANDLER.value))

    def vcfhandler(self):
        pass
        # Modules.VCFHANDLER
        # vcfhandler_group = parser.add_argument_group(Modules.VCFHANDLER.value,
        #                                              '{} parameters'.format(Modules.VCFHANDLER.value))


def main():
    CNVFinder()


if __name__ == '__main__':
    main()
