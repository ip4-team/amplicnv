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


def get_arg_name_from_enum(arg: Enum):
    if type(arg.value) == dict:
        return '--' + list(arg.value.keys())[-1]
    return '--' + arg.name.replace('_', '-')


def get_arg_help_from_enum(arg: Enum):
    if type(arg.value) == dict:
        return list(arg.value.items())[-1][-1]
    return arg.value


class Command(object):
    def __init__(self, name: str, description: str):
        self.name = name
        self.description = description


class Detect(Command):
    def __init__(self):
        super().__init__(self.__class__.__name__.lower(), 'Detect copy number variation in a test sample applying '
                                                          'read depth and variant data (optional)')

    @staticmethod
    def run():
        pass


class Bedloader(Command):
    def __init__(self):
        super().__init__(self.__class__.__name__.lower(), 'Preprocess amplicons defined in a BED file')

    @staticmethod
    def run():
        pass


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


class Compare(Command):
    def __init__(self):
        super().__init__(self.__class__.__name__.lower(), 'Compare a test sample with a baseline of samples '
                                                          'considering read depth')

    @staticmethod
    def run():
        pass


class BafCompute(Command):
    def __init__(self):
        super().__init__(self.__class__.__name__.lower(), 'Compute B-allele frequency (BAF)')

    @staticmethod
    def run():
        pass


class VCFCompare(Command):
    def __init__(self):
        super().__init__(self.__class__.__name__.lower(), 'Compare a test sample with a baseline of samples '
                                                          'considering B-allele frequency and other variant data')

    @staticmethod
    def run():
        pass


class ArgDesc(Enum):
    target = 'Path to file in which sequencing amplicons are listed'
    region = 'Limit target definition to a given region. It should be in the form: chr1:10000-90000'
    spacing = 'Number of nucleotides to ignore at amplicon start and end, to avoid overlapping reads'
    min_data = 'Minimum number of nucleotides for a valid target'
    max_pool = 'Maximum number of origin pools allowed for a valid target'
    bamfile = 'alignment filename (bam format)'
    parallel = 'Count target read depth in parallel'
    output = 'Output filename'
    b_baseline = {'baseline': 'Path to baseline sample bamfile. This parameter can be passed multiple times: '
                              '--baseline file1.bam --baseline file2.bam'}
    b_test = {'test': 'Path to test sample bamfile'}
    v_baseline = {'baseline': 'Path to baseline sample vcf file. This parameter can be passed multiple times: '
                              '--baseline file1.vcf.gz --baseline file2.vcf.gz'}
    size = 'Block size when sliding window'
    step = 'Step size when sliding window'
    metric = 'Define which metric should be used when comparing'
    interval_range = 'Value to multiply metric by'
    min_read = 'Minimum number of reads expected for valid targets'
    below_cutoff = 'Filter out data (ratios) below this cutoff'
    above_cutoff = 'Filter out data (ratios) above this cutoff'
    max_dist = 'Maximum distance allowed of a cnv-like block, to its closest cnv block, for it be a cnv as well'
    cnv_like_range = 'Value to multiply interval_range by in order to detect cnv-like (CNVs when applying looser ' \
                     'calculations) '
    bins = 'Number of bins to use when plotting ratio data'
    method = ''
    outdir = 'Output directory name'
    vcf = 'Path to variant file (VCF)'


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
    \t{}
    
    {}
    \t{}
    
    {}
    \t{}
            
    {}
    \t{}
    
    {}
    \t{}

    {}
    \t{}
    
For getting help of a specific command use: cnvfinder <command> --help'''.format(Detect().name, Detect().description,
                                                                                 Count().name, Count().description,
                                                                                 Compare().name, Compare().description,
                                                                                 BafCompute().name,
                                                                                 BafCompute().description,
                                                                                 VCFCompare().name,
                                                                                 VCFCompare().description,
                                                                                 Bedloader().name,
                                                                                 Bedloader().description))

        parser.add_argument('command', help='Module to run')
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            sys.exit('Unrecognized {} command'.format(args.command))

        getattr(self, args.command)()

    @overrides(ICNVfinder)
    def bedloader(self):
        parser = create_parser(Bedloader().description,
                               command=Bedloader().name)
        parser.add_argument(get_arg_name_from_enum(ArgDesc.target), type=str, required=True,
                            help=get_arg_help_from_enum(ArgDesc.target))
        parser.add_argument(get_arg_name_from_enum(ArgDesc.region), type=str,
                            help=get_arg_help_from_enum(ArgDesc.region))
        parser.add_argument(get_arg_name_from_enum(ArgDesc.spacing), type=int, default=20,
                            help=get_arg_help_from_enum(ArgDesc.spacing))
        parser.add_argument(get_arg_name_from_enum(ArgDesc.min_data), type=int, default=50,
                            help=get_arg_help_from_enum(ArgDesc.min_data))
        parser.add_argument(get_arg_name_from_enum(ArgDesc.max_pool), type=int, default=None,
                            help=get_arg_help_from_enum(ArgDesc.max_pool))
        self.args = parse_sub_command(parser)

    @overrides(ICNVfinder)
    def count(self):
        parser = create_parser(Count().description,
                               command=Count().name)
        parser.add_argument(get_arg_name_from_enum(ArgDesc.target), type=str, required=True,
                            help=get_arg_help_from_enum(ArgDesc.target))
        parser.add_argument(get_arg_name_from_enum(ArgDesc.bamfile), type=str, required=True,
                            help=get_arg_help_from_enum(ArgDesc.bamfile))
        parser.add_argument(get_arg_name_from_enum(ArgDesc.region), type=str,
                            help=get_arg_help_from_enum(ArgDesc.region))
        parser.add_argument(get_arg_name_from_enum(ArgDesc.parallel), dest=ArgDesc.parallel, action='store_true',
                            default=True, help=get_arg_help_from_enum(ArgDesc.parallel))
        parser.add_argument(get_arg_name_from_enum(ArgDesc.output), type=str, required=True,
                            help=get_arg_help_from_enum(ArgDesc.output))
        self.args = parse_sub_command(parser)

    @overrides(ICNVfinder)
    def compare(self):
        parser = create_parser(Compare().description,
                               command=Compare().name)
        parser.add_argument(get_arg_name_from_enum(ArgDesc.b_baseline), required=True, action='append', nargs='?',
                            help=get_arg_help_from_enum(ArgDesc.b_baseline))
        parser.add_argument(get_arg_name_from_enum(ArgDesc.b_test), type=str, required=True,
                            help=get_arg_help_from_enum(ArgDesc.b_test))
        parser.add_argument(get_arg_name_from_enum(ArgDesc.target), type=str, required=True,
                            help=get_arg_help_from_enum(ArgDesc.target))
        parser.add_argument(get_arg_name_from_enum(ArgDesc.size), type=int, default=200,
                            help=get_arg_help_from_enum(ArgDesc.size))
        parser.add_argument(get_arg_name_from_enum(ArgDesc.step), type=int, default=10,
                            help=get_arg_help_from_enum(ArgDesc.step))
        parser.add_argument(get_arg_name_from_enum(ArgDesc.metric), type=str, default='std', choices={'std', 'IQR'},
                            help=get_arg_help_from_enum(ArgDesc.metric))
        parser.add_argument(get_arg_name_from_enum(ArgDesc.interval_range), type=float, default=3,
                            help=get_arg_help_from_enum(ArgDesc.interval_range))
        parser.add_argument(get_arg_name_from_enum(ArgDesc.min_read), type=int, default=30,
                            help=get_arg_help_from_enum(ArgDesc.min_read))
        parser.add_argument(get_arg_name_from_enum(ArgDesc.below_cutoff), type=float, default=0.7,
                            help=get_arg_help_from_enum(ArgDesc.below_cutoff))
        parser.add_argument(get_arg_name_from_enum(ArgDesc.above_cutoff), type=float, default=1.3,
                            help=get_arg_help_from_enum(ArgDesc.above_cutoff))
        parser.add_argument(get_arg_name_from_enum(ArgDesc.max_dist), type=int, default=15000000,
                            help=get_arg_help_from_enum(ArgDesc.max_dist))
        parser.add_argument(get_arg_name_from_enum(ArgDesc.cnv_like_range), type=float, default=0.7,
                            help=get_arg_help_from_enum(ArgDesc.cnv_like_range))
        parser.add_argument(get_arg_name_from_enum(ArgDesc.bins), type=int, default=500,
                            help=get_arg_help_from_enum(ArgDesc.bins))
        parser.add_argument(get_arg_name_from_enum(ArgDesc.method), type=str,
                            default='chr_group', choices={'chr_group'})
        parser.add_argument(get_arg_name_from_enum(ArgDesc.outdir), type=str, default='results',
                            help=get_arg_help_from_enum(ArgDesc.outdir))
        print(parser)
        self.args = parse_sub_command(parser)

    @overrides(ICNVfinder)
    def detect(self):
        parser = create_parser(Detect().description,
                               command=Detect().name)
        self.args = parse_sub_command(parser)

    @overrides(ICNVfinder)
    def bafcompute(self):
        parser = create_parser(BafCompute().description,
                               command=BafCompute().name)
        parser.add_argument('--vcf', type=str, required=True,
                            help='Path to variant file (VCF)')
        parser.add_argument('--output', type=str, required=True,
                            help='Output filename')
        self.args = parse_sub_command(parser)

    @overrides(ICNVfinder)
    def vcfcompare(self):
        parser = create_parser(VCFCompare().description,
                               command=VCFCompare().name)
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
