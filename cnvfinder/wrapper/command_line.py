#!/usr/bin/env python3
"""

@author: valengo
"""

import argparse
import sys
from abc import ABCMeta
from abc import abstractmethod
from enum import Enum

from cnvfinder.bedloader import ROI, bedwrite
from cnvfinder.nrrhandler import NRR
from cnvfinder.vcfhandler import VCF


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
    v_test = 'Path to test sample vcf file'
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


class CNVFinderWrapper(object):
    def __init__(self):

        commands = [command.lower() for command in dir(self) if not command.startswith('_')]
        available_commands = ''
        for command in commands:
            try:
                available_commands += '''
    {}
    \t{}
            '''.format(getattr(self, command.capitalize())().name,
                       getattr(self, command.capitalize())().description)
            except AttributeError:
                sys.exit('Class \"{}\" is not eligible for a command. Does it extend {} class?'.format(
                    command, self._Command.__name__))

        parser = create_parser('CNVfinder is a Python 3.x package for copy number (CNV) '
                               'variation detection on whole exome sequencing (WES) data from '
                               'amplicon-based enrichment technologies',
                               usage='''cnvfinder <command> [<args>]

Available commands:
  {}

For getting help of a specific command use: cnvfinder <command> --help'''.format(available_commands))

        parser.add_argument('command', help='Module to run')
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command.capitalize()):
            sys.exit('Unrecognized {} command'.format(args.command))

        getattr(self, args.command.capitalize())().run()

    class _Command(metaclass=ABCMeta):
        def __init__(self, name: str, description: str):
            self.name = name
            self.description = description

        @abstractmethod
        def run(self):
            pass

    class Detect(_Command):
        def __init__(self):
            super().__init__(self.__class__.__name__.lower(), 'Detect copy number variation in a test sample applying '
                                                              'read depth and variant data (optional)')

        def run(self):
            parser = create_parser(self.description,
                                   command=self.name)

            args = parse_sub_command(parser)

    class Count(_Command):
        def __init__(self):
            super().__init__(self.__class__.__name__.lower(), 'Count the number of reads aligned to each target')

        def run(self):
            parser = create_parser(self.description,
                                   self.name)
            parser.add_argument(get_arg_name_from_enum(ArgDesc.target), type=str, required=True,
                                help=get_arg_help_from_enum(ArgDesc.target))
            parser.add_argument(get_arg_name_from_enum(ArgDesc.bamfile), type=str, required=True,
                                help=get_arg_help_from_enum(ArgDesc.bamfile))
            parser.add_argument(get_arg_name_from_enum(ArgDesc.region), type=str,
                                help=get_arg_help_from_enum(ArgDesc.region))
            parser.add_argument(get_arg_name_from_enum(ArgDesc.parallel), dest=ArgDesc.parallel.name,
                                action='store_true', default=True, help=get_arg_help_from_enum(ArgDesc.parallel))
            parser.add_argument(get_arg_name_from_enum(ArgDesc.output), type=str, required=True,
                                help=get_arg_help_from_enum(ArgDesc.output))
            args = parse_sub_command(parser)

            nrr = NRR(bedfile=args.target,
                      bamfile=args.bamfile,
                      region=args.region,
                      parallel=args.parallel)
            nrr.save(args.output)

    class Compare(_Command):
        def __init__(self):
            super().__init__(self.__class__.__name__.lower(), 'Compare a test sample with a baseline of samples '

                                                              'considering read depth')

        def run(self):
            parser = create_parser(self.description,
                                   command=self.name)
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
            args = parse_sub_command(parser)

    class Bafcompute(_Command):
        def __init__(self):
            super().__init__(self.__class__.__name__.lower(), 'Compute B-allele frequency (BAF)')

        def run(self):
            parser = create_parser(self.description,
                                   command=self.name)
            parser.add_argument(get_arg_name_from_enum(ArgDesc.vcf), type=str, required=True,
                                help=get_arg_help_from_enum(ArgDesc.vcf))
            parser.add_argument(get_arg_name_from_enum(ArgDesc.output), type=str, required=True,
                                help=get_arg_help_from_enum(ArgDesc.output))
            args = parse_sub_command(parser)

            vcf = VCF(args.vcf)
            df = vcf.variants.drop(columns=['info'])
            bedwrite(args.output, df)

    class Vcfcompare(_Command):
        def __init__(self):
            super().__init__(self.__class__.__name__.lower(), 'Compare a test sample with a baseline of samples '
                                                              'considering B-allele frequency and other variant data')

        def run(self):
            parser = create_parser(self.description,
                                   self.name)
            parser.add_argument(get_arg_name_from_enum(ArgDesc.v_baseline), required=True, action='append', nargs='?',
                                help=get_arg_help_from_enum(ArgDesc.v_baseline))
            parser.add_argument(get_arg_name_from_enum(ArgDesc.v_test), type=str, required=True,
                                help=get_arg_help_from_enum(ArgDesc.v_test))
            parser.add_argument(get_arg_name_from_enum(ArgDesc.metric), type=str, default='IQR', choices={'std', 'IQR'},
                                help=get_arg_help_from_enum(ArgDesc.metric))
            parser.add_argument(get_arg_name_from_enum(ArgDesc.interval_range), type=float, default=1.5,
                                help=get_arg_help_from_enum(ArgDesc.interval_range))
            parser.add_argument(get_arg_name_from_enum(ArgDesc.size), type=int, default=400,
                                help=get_arg_help_from_enum(ArgDesc.size))
            parser.add_argument(get_arg_name_from_enum(ArgDesc.step), type=int, default=40,
                                help=get_arg_help_from_enum(ArgDesc.step))
            parser.add_argument(get_arg_name_from_enum(ArgDesc.cnv_like_range), type=float, default=0.7,
                                help=get_arg_help_from_enum(ArgDesc.cnv_like_range))
            parser.add_argument(get_arg_name_from_enum(ArgDesc.max_dist), type=int, default=15000000,
                                help=get_arg_help_from_enum(ArgDesc.max_dist))
            parser.add_argument(get_arg_name_from_enum(ArgDesc.outdir), type=str, default='results',
                                help=get_arg_help_from_enum(ArgDesc.outdir))
            args = parse_sub_command(parser)

    class Bedloader(_Command):
        def __init__(self):
            super().__init__(self.__class__.__name__.lower(), 'Preprocess amplicons defined in a BED file')

        def run(self):
            parser = create_parser(self.description,
                                   command=self.name)
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
            parser.add_argument(get_arg_name_from_enum(ArgDesc.output), type=str, required=True,
                                help=get_arg_help_from_enum(ArgDesc.output))
            args = parse_sub_command(parser)

            roi = ROI(args.target, region=args.region, spacing=args.spacing,
                      mindata=args.min_data, maxpool=args.max_pool)

            bedwrite(args.output, roi.targets)


def main():
    CNVFinderWrapper()
