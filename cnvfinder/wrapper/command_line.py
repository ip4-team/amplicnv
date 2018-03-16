#!/usr/bin/env python3
"""

@author: valengo
"""

import argparse
import sys
from abc import ABCMeta
from abc import abstractmethod
from enum import Enum

from cnvfinder.nrrhandler import NRR
from ..utils import overrides


def create_parser(description: str, command: str = 'command', usage: str = None) -> argparse.ArgumentParser:
    if usage is None:
        usage = 'cnvfinder {} [<args>]'.format(command)
    return argparse.ArgumentParser(description=description, usage=usage,
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)


def parse_sub_command(parser: argparse.ArgumentParser) -> argparse.Namespace:
    return parser.parse_args(sys.argv[2:])


class Command(metaclass=ABCMeta):
    def __init__(self, name: str, description: str):
        self.name = name
        self.description = description


class Count(Command):
    def __init__(self):
        super().__init__(self.__class__.__name__.lower(), 'Count the number of reads aligned to each target')

    @staticmethod
    def run(bedfile, bamfile, region, parallel, output):
        nrr = NRR(bedfile=bedfile,
                  bamfile=bamfile,
                  region=region,
                  parallel=parallel)
        nrr.save(output)


class SubCommands(Enum):
    # bedloader module
    BEDLOADER = 'bedloader'
    # cnvdetector module
    DETECT = 'detect'
    # vcfhandler module
    BAFCOMPUTE = 'bafcompute'
    VCFCOMPARE = 'vcfcompare'
    VCFPLOT = 'vcfplot'  # TODO: add this functionality
    # nrrhandler module
    COMPARE = 'compare'


class ArgDesc(Enum):
    target = 'Path to file in which sequencing amplicons are listed'
    region = 'Limit target definition to a given region. It should be in the form: chr1:10000-90000'
    spacing = 'Number of nucleotides to ignore at amplicon start and end, to avoid overlapping reads'
    min_data = 'Minimum number of nucleotides for a valid target'
    max_pool = 'Maximum number of origin pools allowed for a valid target'
    bamfile = 'alignment filename (bam format)'
    parallel = 'Count target read depth in parallel'
    output = 'Output filename'


class ICNVfinder(metaclass=ABCMeta):
    @abstractmethod
    def detect(self):
        pass

    @abstractmethod
    def count(self):
        pass

    @abstractmethod
    def compare(self):
        pass

    @abstractmethod
    def bafcompute(self):
        pass

    @abstractmethod
    def vcfcompare(self):
        pass

    @abstractmethod
    def bedloader(self):
        pass


class Wrapper(ICNVfinder):
    def __init__(self):
        self.args = None
        parser = create_parser('CNVfinder is a Python 3.x package for copy number (CNV) '
                               'variation detection on whole exome sequencing (WES) data from '
                               'amplicon-based enrichment technologies',
                               usage='''cnvfinder <command> [<args>]
                               
Available commands used in various situations:

    {}
    \tDetect copy number variation in a test sample applying read depth and variant data (optional)
    
    {}
    \t{}
    
    {}
    \tCompare a test sample with a baseline of samples considering read depth
            
    {}
    \tCompute B-allele frequency (BAF)
    
    {}
    \tCompare a test sample with a baseline of samples considering B-allele frequency and other variant data

    {}
    \tPreprocess amplicons defined in a BED file
    
For getting help of a specific command use: cnvfinder <command> --help
    '''.format(SubCommands.DETECT.value, Count().name, Count().description, SubCommands.COMPARE.value,
               SubCommands.BAFCOMPUTE.value, SubCommands.VCFCOMPARE.value, SubCommands.BEDLOADER.value))

        parser.add_argument('command', help='Module to run')
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            sys.exit('Unrecognized {} command'.format(args.command))

        getattr(self, args.command)()

    @overrides(ICNVfinder)
    def bedloader(self):
        parser = create_parser('Preprocess amplicons defined in a BED file',
                               command=SubCommands.BEDLOADER.value)
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

    @overrides(ICNVfinder)
    def count(self):
        parser = create_parser(Count().description,
                               command=Count().name)
        parser.add_argument('--target', type=str, required=True,
                            help='Path to file in which sequencing amplicons are listed')
        parser.add_argument('--bamfile', type=str, required=True,
                            help='alignment filename (bam format)')
        parser.add_argument('--region', type=str,
                            help='Limit target definition to a given region. It should be in the form: '
                                 'chr1:10000-90000')
        parser.add_argument('--parallel', dest='parallel', action='store_true', default=True,
                            help='Count target read depth in parallel')
        parser.add_argument('--output', type=str, required=True,
                            help='Output filename')
        self.args = parse_sub_command(parser)

    @overrides(ICNVfinder)
    def compare(self):
        parser = create_parser('Compare a test sample with a baseline of samples considering read depth',
                               command=SubCommands.COMPARE.value)
        parser.add_argument('--baseline', required=True, action='append', nargs='?',
                            help='Path to baseline sample bamfile. This parameter can be passed multiple times: '
                                 '--baseline file1.bam --baseline file2.bam')
        parser.add_argument('--test', type=str, required=True,
                            help='Path to test sample bamfile')
        parser.add_argument('--target', type=str, required=True,
                            help='Path to file in which sequencing amplicons are listed')
        parser.add_argument('--size', type=int, default=200,
                            help='Block size when sliding window')
        parser.add_argument('--step', type=int, default=10,
                            help='Step size when sliding window')
        parser.add_argument('--metric', type=str, default='std', choices={'std', 'IQR'},
                            help='param used to define which metric should be used when comparing')
        parser.add_argument('--interval-range', type=float, default=3,
                            help='Value to multiply metric by')
        parser.add_argument('--min-read', type=int, default=30,
                            help='Minimum number of reads expected for valid targets')
        parser.add_argument('--below-cutoff', type=float, default=0.7,
                            help='Filter out data (ratios) below this cutoff')
        parser.add_argument('--above-cutoff', type=float, default=1.3,
                            help='Filter out data (ratios) above this cutoff')
        parser.add_argument('--max-dist', type=int, default=15000000,
                            help='Maximum distance allowed of a cnv-like block, to its closest cnv block, for it be a '
                                 'cnv as well')
        parser.add_argument('--cnv-like-range', type=float, default=0.7,
                            help='Value to multiply interval_range by in order to detect cnv-like (CNVs when applying '
                                 'looser calculations)')
        parser.add_argument('--bins', type=int, default=500,
                            help='Number of bins to use when plotting ratio data')
        parser.add_argument('--method', type=str, default='chr_group', choices={'chr_group'})
        parser.add_argument('--output', type=str, default='results',
                            help='Output directory name')
        self.args = parse_sub_command(parser)

    @overrides(ICNVfinder)
    def detect(self):
        parser = create_parser('Detect copy number variation in a test sample applying read depth and variant data ('
                               'optional)',
                               command=SubCommands.DETECT.value)
        self.args = parse_sub_command(parser)

    @overrides(ICNVfinder)
    def bafcompute(self):
        parser = create_parser('Compute B-allele frequency (BAF)',
                               command=SubCommands.BAFCOMPUTE.value)
        parser.add_argument('--vcf', type=str, required=True,
                            help='Path to variant file (VCF)')
        parser.add_argument('--output', type=str, required=True,
                            help='Output filename')
        self.args = parse_sub_command(parser)

    @overrides(ICNVfinder)
    def vcfcompare(self):
        parser = create_parser('Compare a test sample with a baseline of samples considering B-allele frequency',
                               command=SubCommands.VCFCOMPARE.value)
        parser.add_argument('--baseline', required=True, action='append', nargs='?',
                            help='Path to baseline sample vcf file. This parameter can be passed multiple times: '
                                 '--baseline file1.vcf.gz --baseline file2.vcf.gz')
        parser.add_argument('--test', type=str, required=True,
                            help='Path to test sample vcf file')
        parser.add_argument('--metric', type=str, default='std', choices={'std', 'IQR'},
                            help='param used to define which metric should be used when comparing')
        parser.add_argument('--interval-range', type=float, default=3,
                            help='Value to multiply metric by')
        parser.add_argument('--size', type=int, default=400,
                            help='Block size when sliding window')
        parser.add_argument('--step', type=int, default=40,
                            help='Step size when sliding window')
        parser.add_argument('--cnv-like-range', type=float, default=0.7,
                            help='Value to multiply interval_range by in order to detect cnv-like (CNVs when applying '
                                 'looser calculations)')
        parser.add_argument('--max-dist', type=int, default=15000000,
                            help='Maximum distance allowed of a cnv-like block, to its closest cnv block, for it be a '
                                 'cnv as well')
        parser.add_argument('--output', type=str, default='results',
                            help='Output directory name')
        self.args = parse_sub_command(parser)

    def show_args(self):
        opts = [opt for opt in dir(self.args) if not opt.startswith('_')]
        for opt in opts:
            print('{}: {}'.format(opt, getattr(self.args, opt)))


def main():
    opts = Wrapper()
