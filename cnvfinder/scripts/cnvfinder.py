#!/usr/bin/env python3
"""

@author: valengo
"""

import argparse
import sys
from enum import Enum


class SubCommands(Enum):
    BEDLOADER = 'bedloader'
    CNVDETECTOR = 'cnvdetector'
    NRRHANDLER = 'nrrhandler'
    VCFHANDLER = 'vcfhandler'
    COUNT = 'count'
    COMPARE = 'compare'


def create_parser(description: str, usage: str = None) -> argparse.ArgumentParser:
    return argparse.ArgumentParser(description=description, usage=usage,
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)


def parse_sub_command(parser: argparse.ArgumentParser) -> argparse.Namespace:
    return parser.parse_args(sys.argv[2:])


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
    \tCount the number of reads aligned to each target
    
    {}
    \tCompare a test sample with a baseline of samples considering read depth
    
    {}
    \tDetect copy number variation on whole exome sequencing data
    \tCreate plots
    
    {}
    \tCreate read depth and ratio related plots
    
    {}
    \tFilter variants according to a baseline
    \tCompute B-allele frequency (BAF)
    \tCreate BAF related plots
    
For getting help of a specific command use: cnvfinder <command> --help
    '''.format(SubCommands.BEDLOADER.value, SubCommands.COUNT.value, SubCommands.COMPARE.value, SubCommands.CNVDETECTOR.value, SubCommands.NRRHANDLER.value, SubCommands.VCFHANDLER.value))

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
        self.args = parse_sub_command(parser)

    def count(self):
        parser = create_parser('Count the number of reads aligned to each target',
                               usage='''cnvfinder count [<args>]''')
        parser.add_argument('--target', type=str, required=True,
                            help='Path to file in which sequencing amplicons are listed')
        parser.add_argument('--bamfile', type=str, required=True,
                            help='alignment filename (bam format)')
        parser.add_argument('--region', type=str,
                            help='Limit target definition to a given region. It should be in the form: '
                                 'chr1:10000-90000')
        parser.add_argument('--parallel', dest='parallel', action='store_true', default=True,
                            help='Count target read depth in parallel')
        self.args = parse_sub_command(parser)

    def compare(self):
        parser = create_parser('Compare a test sample with a baseline of samples considering read depth',
                               usage='''cnvfinder count [<args>]''')
        self.args = parse_sub_command(parser)

    def cnvdetector(self):
        pass
        # Modules.CNVDETECTOR
        # cnvdetector_group = parser.add_argument_group(Modules.CNVDETECTOR.value,
        #
        #                                               '{} parameters'.format(Modules.CNVDETECTOR.value))

    def nrrhandler(self):
        parser = create_parser('')

    def vcfhandler(self):
        pass
        # Modules.VCFHANDLER
        # vcfhandler_group = parser.add_argument_group(Modules.VCFHANDLER.value,
        #                                              '{} parameters'.format(Modules.VCFHANDLER.value))


def main():
    CNVFinder()


if __name__ == '__main__':
    main()
