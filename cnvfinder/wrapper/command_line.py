#!/usr/bin/env python3
"""

@author: valengo
"""

import sys
from abc import ABCMeta
from abc import abstractmethod
from collections import defaultdict

from cnvfinder import CNVTest
from cnvfinder.bedloader import ROI, bedwrite
from cnvfinder.nrrhandler import NRR, NRRList, NRRTest
from cnvfinder.vcfhandler import VCF, VCFList, VCFTest
from cnvfinder.wrapper.argdesc import ArgDesc, Strings
from cnvfinder.wrapper.utils import get_arg_name_from_enum, get_arg_help_from_enum, getattr_by, parse_sub_command, \
    create_parser


class CNVFinderWrapper(object):
    def __init__(self):

        commands = [command.lower() for command in dir(self) if command[0].isupper()]
        self.instances = defaultdict(lambda: None)

        available_commands = ''
        for command in commands:
            try:
                # create instance of "command" classes
                self.instances[command] = getattr(self, command.capitalize())(self)
                available_commands += '''
    {}
    \t{}
            '''.format(self.instances[command].name,
                       self.instances[command].description)
            except AttributeError:
                sys.exit(Strings.wrong_command.value.format(
                    command, self._Command.__name__))

        parser = create_parser(Strings.description.value,
                               usage='''{}

{}
  {}

{}'''.format(Strings.usage.value.format(Strings.command.value), Strings.available_commands.value, available_commands,
             Strings.getting_help.value))

        parser.add_argument(Strings.command.value, help=Strings.command_help.value)
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command.capitalize()):
            sys.exit(Strings.unrecognized_command.value.format(args.command))

        self.instances[args.command].run()

    def get_instance_by_name(self, name: str):
        ins = self.instances.get(name)
        if ins is None:
            sys.exit(Strings.get_command_error.value.format(name, Strings.issue_url.value))
        return ins

    class _Command(metaclass=ABCMeta):
        name = ''
        description = ''

        def __init__(self, name: str, description: str, outer_instance):
            type(self).name = name
            type(self).description = description
            self.outer_instance = outer_instance

        @abstractmethod
        def run(self):
            pass

    class Detect(_Command):
        def __init__(self, outer_instance):
            super().__init__(self.__class__.__name__.lower(), Strings.detect_description.value, outer_instance)

        def run(self):
            parser = create_parser(self.description,
                                   command=self.name)
            parser.add_argument(get_arg_name_from_enum(ArgDesc.no_rd), dest=ArgDesc.no_rd.name,
                                action='store_true', default=False, help=get_arg_help_from_enum(ArgDesc.no_rd))
            parser.add_argument(get_arg_name_from_enum(ArgDesc.no_vcf), dest=ArgDesc.no_vcf.name,
                                action='store_true', default=False, help=get_arg_help_from_enum(ArgDesc.no_vcf))
            parser.add_argument(get_arg_name_from_enum(ArgDesc.outdir), type=str, default=Strings.default_outdir.value,
                                help=get_arg_help_from_enum(ArgDesc.outdir))

            args = parse_sub_command(parser)
            ratio = not getattr_by(ArgDesc.no_rd, args)
            variant = not getattr_by(ArgDesc.no_vcf, args)

            nrrtest = self.outer_instance.get_instance_by_name('compare').run() if ratio else None
            vcftest = self.outer_instance.get_instance_by_name('vcfcompare').run() if ratio else None

            cnvtest = CNVTest(nrrtest, vcftest,
                              path=getattr_by(ArgDesc.outdir, args))
            cnvtest.analyze_plot()
            rfname = '{}/cnvfrombam.csv'.format(cnvtest.path2table)
            vfname = '{}/cnvfromvcf.csv'.format(cnvtest.path2table)
            if ratio:
                cnvtest.blocks[0].to_csv(rfname, sep='\t', index=False)
            if variant:
                cnvtest.blocks[-1].to_csv(vfname, sep='\t', index=False)
            cnvtest.save()
            cnvtest.summary()
            print('Done!')

    class Count(_Command):
        def __init__(self, outer_instance):
            super().__init__(self.__class__.__name__.lower(), Strings.count_description.value, outer_instance)

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

            nrr = NRR(bedfile=getattr_by(ArgDesc.target, args),
                      bamfile=getattr_by(ArgDesc.bamfile, args),
                      region=getattr_by(ArgDesc.region, args),
                      parallel=getattr_by(ArgDesc.parallel, args))
            nrr.save(getattr_by(ArgDesc.output, args))

    class Compare(_Command):
        def __init__(self, outer_instance):
            super().__init__(self.__class__.__name__.lower(), Strings.compare_description.value, outer_instance)

        def run(self) -> NRRTest:
            parser = create_parser(self.description,
                                   command=self.name)
            parser.add_argument(get_arg_name_from_enum(ArgDesc.baseline), required=True, action='append', nargs='?',
                                help=get_arg_help_from_enum(ArgDesc.baseline))
            parser.add_argument(get_arg_name_from_enum(ArgDesc.test), type=str, required=True,
                                help=get_arg_help_from_enum(ArgDesc.test))
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
            parser.add_argument(get_arg_name_from_enum(ArgDesc.outdir), type=str, default=Strings.default_outdir.value,
                                help=get_arg_help_from_enum(ArgDesc.outdir))
            args = parse_sub_command(parser)

            # load amplicons and define targets
            roi = ROI(getattr_by(ArgDesc.target, args))
            # load test sample and count number of reads in regions of targets
            sample = NRR(bamfile=getattr_by(ArgDesc.test, args),
                         bed=roi)
            # load baseline samples and count their number of reads in region of targets
            baseline = NRRList(bamfiles=getattr_by(ArgDesc.baseline, args), bed=roi)
            # make test
            nrrtest = NRRTest(baseline, sample,
                              path=getattr_by(ArgDesc.outdir, args),
                              size=getattr_by(ArgDesc.size, args),
                              step=getattr_by(ArgDesc.step, args),
                              metric=getattr_by(ArgDesc.metric, args),
                              interval_range=getattr_by(ArgDesc.interval_range, args),
                              minread=getattr_by(ArgDesc.min_read, args),
                              below_cutoff=getattr_by(ArgDesc.below_cutoff, args),
                              above_cutoff=getattr_by(ArgDesc.above_cutoff, args),
                              maxdist=getattr_by(ArgDesc.max_dist, args),
                              cnv_like_range=getattr_by(ArgDesc.cnv_like_range, args),
                              bins=getattr_by(ArgDesc.bins, args),
                              method=getattr_by(ArgDesc.method, args))
            nrrtest.makeratio()
            if nrrtest.ratios:
                print('Creating plots at {}'.format(nrrtest.path2plot))
                nrrtest.plot()
                filename = '{}/nrrtest.csv'.format(nrrtest.path2table)
                print('Writing table at {}'.format(filename))
                nrrtest.df.to_csv(filename, sep='\t', index=False)
                print('Done!')

            return nrrtest

    class Bafcompute(_Command):
        def __init__(self, outer_instance):
            super().__init__(self.__class__.__name__.lower(), Strings.bafcompute_description.value, outer_instance)

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
        def __init__(self, outer_instance):
            super().__init__(self.__class__.__name__.lower(), Strings.vcfcompare_description.value, outer_instance)

        def run(self) -> VCFTest:
            parser = create_parser(self.description,
                                   self.name)
            parser.add_argument(get_arg_name_from_enum(ArgDesc.vcf_baseline), required=True, action='append', nargs='?',
                                help=get_arg_help_from_enum(ArgDesc.vcf_baseline))
            parser.add_argument(get_arg_name_from_enum(ArgDesc.vcf_test), type=str, required=True,
                                help=get_arg_help_from_enum(ArgDesc.vcf_test))
            parser.add_argument(get_arg_name_from_enum(ArgDesc.vcf_metric), type=str, default='IQR',
                                choices={'std', 'IQR'}, help=get_arg_help_from_enum(ArgDesc.vcf_metric))
            parser.add_argument(get_arg_name_from_enum(ArgDesc.vcf_interval_range), type=float, default=1.5,
                                help=get_arg_help_from_enum(ArgDesc.vcf_interval_range))
            parser.add_argument(get_arg_name_from_enum(ArgDesc.vcf_size), type=int, default=400,
                                help=get_arg_help_from_enum(ArgDesc.vcf_size))
            parser.add_argument(get_arg_name_from_enum(ArgDesc.vcf_step), type=int, default=40,
                                help=get_arg_help_from_enum(ArgDesc.vcf_step))
            parser.add_argument(get_arg_name_from_enum(ArgDesc.vcf_cnv_like_range), type=float, default=0.7,
                                help=get_arg_help_from_enum(ArgDesc.vcf_cnv_like_range))
            parser.add_argument(get_arg_name_from_enum(ArgDesc.vcf_max_dist), type=int, default=15000000,
                                help=get_arg_help_from_enum(ArgDesc.vcf_max_dist))
            parser.add_argument(get_arg_name_from_enum(ArgDesc.to_filter), dest=ArgDesc.to_filter.name,
                                action='store_true', default=True, help=get_arg_help_from_enum(ArgDesc.to_filter))
            parser.add_argument(get_arg_name_from_enum(ArgDesc.outdir), type=str, default=Strings.default_outdir.value,
                                help=get_arg_help_from_enum(ArgDesc.outdir))
            args = parse_sub_command(parser)

            sample = VCF(getattr_by(ArgDesc.vcf_test, args))
            baseline = VCFList(getattr_by(ArgDesc.vcf_baseline, args))
            vcftest = VCFTest(baseline, sample,
                              metric=getattr_by(ArgDesc.vcf_metric, args),
                              interval_range=getattr_by(ArgDesc.vcf_interval_range, args),
                              size=getattr_by(ArgDesc.vcf_size, args),
                              step=getattr_by(ArgDesc.vcf_step, args),
                              cnv_like_range=getattr_by(ArgDesc.vcf_cnv_like_range, args),
                              maxdist=getattr_by(ArgDesc.vcf_max_dist, args),
                              path=getattr_by(ArgDesc.outdir, args))

            if getattr_by(ArgDesc.to_filter, args):
                vcftest.filter()
                vcftest.load_filtered_out()
                if vcftest.filtered_out_pos is None:
                    vcftest.compare()
                    vcftest.save_filtered_out()
                vcftest.eliminate_vars()
            vcftest.split()
            print('Creating plots at {}'.format(vcftest.path2plot))
            vcftest.vcfplot()
            filename = '{}/vcftest.csv'.format(vcftest.path2table)
            print('Writing table at {}'.format(filename))
            vcftest.df.to_csv(filename, sep='\t', index=False)
            print('Done!')

            return vcftest

    class Bedloader(_Command):
        def __init__(self, outer_instance):
            super().__init__(self.__class__.__name__.lower(), Strings.bedloader_description.value, outer_instance)

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

            roi = ROI(getattr_by(ArgDesc.target, args),
                      region=getattr_by(ArgDesc.region, args),
                      spacing=getattr_by(ArgDesc.spacing, args),
                      mindata=getattr_by(ArgDesc.min_data, args),
                      maxpool=getattr_by(ArgDesc.max_pool, args))

            bedwrite(getattr_by(ArgDesc.output, args), roi.targets)


def main():
    CNVFinderWrapper()
