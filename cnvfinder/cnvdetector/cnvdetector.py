#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: valengo
"""
from cnvfinder.ideogram import Ideogram
from cnvfinder.nrrhandler import NRRTest
from cnvfinder.utils.utils import bedwrite, sort_chroms
from cnvfinder.vcfhandler import VCFTest
from ..utils import ConfigfileParser
from ..nrrhandler import NRRConfig
from ..vcfhandler import VCFConfig
from plotly.offline import plot
from pandas import DataFrame
from ..utils import validstr
from ..graphics import create_subplots
from ..utils import createdir
from ..utils import appenddir


@validstr('path', empty_allowed=False)
class CNVTest(object):
    """
    **Detect and call CNVs based on read depth and/or variant info**

    :param NRRTest nrrtest: read depth data. When it's None, analysis employing "read depth" won't be performed
    :param vcftest: variant data. When it's None, analysis employing "read depth" won't be performed
    :param path: path to output directory
    """
    def __init__(self, nrrtest: NRRTest = None, vcftest: VCFTest = None, path: str = 'results'):
        self.nrrtest = nrrtest
        self.vcftest = vcftest
        self.cnvdf = None
        self.blocks = None
        self.path = path
        self.path2plotcnv = appenddir(self.path, 'plots/bam-vcf')
        self.path2bed = appenddir(self.path, 'bed')
        self.path2table = appenddir(self.path, 'tables/bam-vcf')
        self.rootname = 'analyzed'
        createdir(self.path2table)
        createdir(self.path2plotcnv)
        createdir(self.path2bed)

    def analyze_plot(self):
        """
        Detect CNVs applying ratio and/or BAF data. Output (list and plots) will be saved at self.path
        """
        blocks = []
        rtitles, rlayout, rtraces = [], [], []
        vtitles, vlayout, vtraces = [], [], []

        # compute and analyze
        if self.nrrtest:
            rdf, rtraces, rtitles, rlayout = self.nrrtest.compute_plot(mode='analyzed',
                                                                       filename='all-chrom-ratios.html')
            blocks.append(rdf)
        if self.vcftest:
            vdf, vtraces, vtitles, vlayout = self.vcftest.compute_plot(mode='analyzed', filename='all-chrom-bafs.html')
            blocks.append(vdf)

        # merge cnv blocks
        if self.nrrtest:
            self.cnvdf = DataFrame(self.nrrtest.merge(blocks))

        self.blocks = blocks
        if isinstance(self.cnvdf, DataFrame) and not self.cnvdf.empty:
            print('{} gene size potential CNVs detected!'.format(len(self.cnvdf[3].unique())))
        elif self.nrrtest:
            print('{} gene size potential CNVs detected!'.format(0))

        # plot results
        print('Plotting results')
        if self.vcftest and self.nrrtest:
            # TODO: solve for when there is no data for a chr that is not Y
            titles = []
            # plot for chrom and create titles for whole subplot
            for i, title in enumerate(rtitles):
                titles.extend([title, ''])
                try:
                    rtrace = [rtraces[i]]
                except IndexError:
                    rtrace = [[]]
                try:
                    vtrace = [vtraces[i]]
                except IndexError:
                    vtrace = [[]]
                fig = create_subplots([rtrace, vtrace],
                                      [title, ''],
                                      [rlayout, vlayout],
                                      height=600)
                filename = '{}/chr{}-{}.html'.format(self.path2plotcnv,
                                                     title.split(' ')[-1],
                                                     self.rootname)
                plot(fig, filename=filename, auto_open=False)
            # plot all together
            fig = create_subplots([rtraces, vtraces],
                                  titles,
                                  [rlayout, vlayout])
            plot(fig,
                 filename='{}/{}'.format(self.path2plotcnv,
                                         'all-chrom-ratios-bafs.html'),
                 auto_open=False)

    def save(self):
        """
        Save list of detected CNVs
        """
        try:
            filename = '{}/{}{}'.format(self.path2bed,
                                        self.nrrtest.sample.bamfile.split('/')[-1],
                                        '.bed')
            bedwrite(filename, self.cnvdf, header=False)
        except AttributeError:
            if self.nrrtest is not None:
                print('There is no CNV yet!')
                print('You should run cnvtest.detect() ' +
                      'before saving anything!')
            else:
                pass

    def summary(self):
        """
        Write a summary of the tests at self.path/summary.txt
        """
        no_cnv = 'No potential CNVs were found under the specified parameters\n'
        cnv = 'Potential CNVs were found in chromosome(s): {}\n'
        filename = '{}/summary.txt'.format(self.path)
        chroms = []

        with open(filename, 'w') as file:
            if self.blocks[0].empty and self.blocks[-1].empty:
                file.write(no_cnv)
            else:
                for block in self.blocks:
                    if not block.empty:
                        chroms.extend(block['chrom'].unique())
                file.write(cnv.format(', '.join(set(chroms))))

    def create_ideogram(self):
        """
        Create ideograms for the current test
        """
        print('Creating chromosome ideograms...')
        chroms = sort_chroms(self.cnvdf[0].unique())
        ideo = Ideogram(chroms=chroms)

        gain = '#0ab26c'
        loss = '#cc1231'

        # cnvs found using ratios
        if self.nrrtest is not None:
            df = self.blocks[0]
            df_loss = df[df.call == 'loss']
            df_gain = df[df.call == 'gain']
            if not df_gain.empty:
                ideo.add_data_above(df_gain, color=gain)
            if not df_loss.empty:
                ideo.add_data_below(df_loss, color=loss)

        # cnvs found using variant
        if self.vcftest is not None and not self.blocks[-1].empty:
            ideo.add_data_above(self.blocks[-1], color=gain)

        # legend
        legend = [dict(color=gain, label='gain'),
                  dict(color=loss, label='loss')]
        ideo.add_legend(legend)

        ideo.save('{}/{}'.format(self.path, 'ideogram.png'))


@validstr('filename', empty_allowed=False)
class CNVConfig(object):
    """
    **A wrapper for CNVTest. This class loads test parameters from a configuration file**

    :param str filename: path to configuration file
    :param bool ratio: specify the usage of read depth data
    :param bool variant: specify the usage of variant data
    """

    def __init__(self, filename: str, ratio: bool = True, variant: bool = True):
        self.filename = filename
        self.sections_params = {
            'output': 'm'
        }
        # load configfile
        self.config = ConfigfileParser(self.filename,
                                       self.sections_params)
        nrrtest = NRRConfig(filename).nrrtest if ratio else None
        vcftest = VCFConfig(filename).vcftest if variant else None
        self.cnvtest = CNVTest(nrrtest, vcftest,
                               path=self.config.sections['output']['path'])
        self.cnvtest.analyze_plot()
        rfname = '{}/cnvfrombam.csv'.format(self.cnvtest.path2table)
        vfname = '{}/cnvfromvcf.csv'.format(self.cnvtest.path2table)
        if ratio:
            self.cnvtest.blocks[0].to_csv(rfname, sep='\t', index=False)
        if variant:
            self.cnvtest.blocks[-1].to_csv(vfname, sep='\t', index=False)
        self.cnvtest.save()
        self.cnvtest.summary()
        self.cnvtest.create_ideogram()
        print('Done!')
