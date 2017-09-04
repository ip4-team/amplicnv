#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: valengo
"""

from pysam import VariantFile
from collections import defaultdict
from ..utils import Region
from pandas import DataFrame
from plotly.offline import plot
from ..graphics import y_scatter
from ..graphics import histogram
from ..graphics import create_subplots
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from ..commons import ChromDF as cdf
from ..utils import createdir
from ..utils import validstr
from ..utils import ConfigfileParser
from ..utils import ExtensionManager as em
from ..utils import overrides
from ..stats import compute_metric
from ..stats import above_range
from ..stats import isbimodal
from ..utils import appenddir


@validstr('vcffile', empty_allowed=False)
class VCF(object):
    '''
    Hold information on variants loaded from a VCF file
    using pysam.VariantFile
    '''
    def __init__(self, vcffile):
        '''
        Parameters:
             vcffile (str): vcf filename
        '''
        self.vcffile = vcffile
        if self.vcffile.endswith('bam'):
            self.vcffile = em(self.vcffile, 'bam', 'vcf.gz').new_filename
        self.variants = self.__loadVars()

    def __loadVars(self):
        print('Loading {} variants'.format(self.vcffile))

        columns = ['chrom',
                   'start',
                   'stop',
                   'alleles',
                   'alts',
                   'contig',
                   'id',
                   'info',
                   'pos',
                   'qual',
                   'ref',
                   'rid',
                   'rlen',
                   'dp',
                   'baf']

        variants = []

        try:
            vcf_in = VariantFile(self.vcffile)
            lines = list(vcf_in.fetch())
        except OSError as error:
            print('File {} not found!'.format(self.vcffile))
            return None
        except AttributeError as error:
            print('Failed getting data from {} file!'.format(self.vcffile))
            return None

        for line in lines:
            row = []
            for i in range(len(columns) - 2):
                row.append(getattr(line, columns[i]))
            row.append(row[7]['DP'])
            row.append(self.calc_baf(row[7]))
            variants.append(tuple(row))

        # update column names
        columns[1] = 'chromStart'
        columns[2] = 'chromEnd'

        print('{} variants were loaded!'.format(len(variants)))
        return DataFrame(variants, columns=columns)

    def calc_baf(self, info):
        '''
        Compute the B allele frequency

        Parameters:
             info (object): variant info object from pysam

        Returns:
             baf (number): B allele frequency
        '''

        baf = info['AO'][0]/(info['AO'][0] + info['RO'])

        return baf


class VCFList(object):
    def __init__(self, vcffiles):
        '''
        Parameters:
             vcffiles (list): vcf filenames
        '''
        self.vcffiles = vcffiles
        self.list = []

        for i, vcffile in enumerate(vcffiles):
            if vcffile.endswith('bam'):
                self.vcffiles[i] = em(vcffile, 'bam', 'vcf.gz').new_filename
            self.list.append(VCF(self.vcffiles[i]))


@validstr('configfile', empty_allowed=True)
class VCFTest(cdf):
    '''
    Hold information about tests between a VCF baseline
    and a VCF test sample
    '''
    def __init__(self, baseline, sample, metric='IQR',
                 interval_range=1.5, size=400, step=40,
                 cnv_like_range=0.7, maxdist=15000000,
                 path='results'):
        '''
        Parameters:
             baseline (VCFList): obj representing variants of the baseline
             sample (VCF): obj representing variants of the sample test

             metric (str): param used to define which metric should be used
             when detecting CNVs:
                  'cutoff' = predefined cutoffs

             interval_range (number): value to multiply metric by
             size (number): block size when sliding window
             step (number): step size when sliding window

             cnv_like_range (number): value to multiply interval_range by in
             order to detect cnv-like (cnvs when applying looser calculations)

             maxdist (number): maximum distance allowed of a cnv-like block, to
             its closest cnv block, for it be a cnv as well.

             path (str): where to save (path/filename) analysis results
        '''

        self.baseline = baseline
        self.sample = sample
        self.size = size
        self.step = step
        self.metric = metric
        self.interval_range = interval_range
        self.cnv_like_range = cnv_like_range
        self.maxdist = maxdist
        self.path = path
        self.rootname = 'bafs'
        self.path2plot = appenddir(self.path, 'plots/vcf/baf')
        self.path2table = appenddir(self.path, 'tables/vcf')
        self.path2plotcnv = appenddir(self.path, 'plots/vcf/cnv')
        createdir(self.path2plot)
        createdir(self.path2table)
        createdir(self.path2plotcnv)
        self._createdf()

    def compare(self, nmin=None, maxdiff=0.00001):
        '''
        Apply filters on sample variants considering variants
        from baseline files

        Current filters:
             similar_ratios: given a variant, filter it out if all AO/DP
             are similar considering sample and baseline files

        Parameters:
             nmin (number): min number samples (in baseline) to one variant to
             be considered error when appearing with the same frequency
             and position
        '''
        if nmin is None:
            nmin = len(self.baseline.list)

        print('Analyzing {} file'.format(self.sample.vcffile))
        eliminated = 0
        filtered_out = []
        filtered_out_pos = []
        for variant in self.sample.variants.itertuples():
            similar_ratios = 0

            for i in range(len(self.baseline.list)):
                # add filter here
                similar_ratios += self.__similar_ratios(variant,
                                                        self.baseline.list[i].variants,
                                                        maxdiff)
            if similar_ratios >= nmin:
                eliminated += 1
                filtered_out.append([variant.pos,
                                     variant.chrom,
                                     variant.rlen,
                                     variant.alleles,
                                     variant.alts,
                                     variant.ref,
                                     variant.info['AO'],
                                     variant.info['RO'],
                                     variant.info['DP']])
                filtered_out_pos.append((variant.chrom, variant.pos))
        self.eliminated = eliminated
        self.filtered_out = filtered_out
        self.filtered_out_pos = filtered_out_pos
        self.sample.eliminate_vars = self.filtered_out_pos
        self._createdf()
        self.split()

        print('Call VCFTest.save_filtered_out() if you want ' +
              'to save filtered_out variants on a text file!')

    @overrides(cdf)
    def _createdf(self):
        self.df = self.sample.variants.copy()

    def __similar_ratios(self, row, copy2compare, maxdiff):
        similar_vars = copy2compare.loc[(copy2compare['pos'] == row.pos) &
                                        (copy2compare['chrom'] == row.chrom)]
        if not similar_vars.empty:
            for variant in similar_vars.itertuples():
                try:
                    diff = abs(variant.baf - row.baf)
                    if diff < maxdiff:
                        return 1
                except KeyError as error:
                    print(variant.info.keys(), variant.pos)
                    print(row.info.keys(), row.pos)
        return 0

    def save_filtered_out(self):
        '''
        For each sample, this method saves filtered out variants on file
        '''
        try:
            print('Saving filtered out variants of {}'.format(self.sample.vcffile))
            with open(self.sample.vcffile + '.txt', 'w') as file:
                for variant in self.filtered_out:
                    for i in range(len(variant) - 1):
                        file.write('{}\t'.format(variant[i]))
                    file.write('{}\n'.format(variant[-1]))
                print('Done!')

        except FileNotFoundError as error:
            print(error)

    def load_filtered_out(self):
        '''
        Load filtered out variants from file
        '''
        variants = []
        variants_pos = []

        print('Loading list of filtered out variants' +
              ' from {} file'.format(self.sample.vcffile))

        try:
            with open(self.sample.vcffile + '.txt', 'r') as file:
                lines = file.readlines()
                for line in lines:
                    variants.append(line.strip('\n').split('\t'))
                    variants_pos.append((line.split('\t')[1],
                                         int(line.split('\t')[0])))

                self.filtered_out = variants
                self.filtered_out_pos = variants_pos
                print('Done!')
        except FileNotFoundError as error:
            print('"{}.txt" file not found!'.format(self.sample.vcffile))
            print('In this case, you could run VCFTest.compare()')
            self.filtered_out = None
            self.filtered_out_pos = None

    def __splitVars(self, homozygous_freq=0.90):
        '''
        From self.df, this method creates two new dataframes:
        one for homozygous variants and another for heterozygous
        variants

        Parameters:
           homozygous_freq (number): relative homozygous frequency
        '''
        df = self.df
        return df[df.baf >= homozygous_freq], df[df.baf < homozygous_freq]

    @overrides(cdf)
    def getin(self, region=None, column=None, mode=None):
        df = self.getview(mode=mode)
        rows = self.getrows(df, region)
        if rows is not None and column:
            return list(rows.loc[:, column])
        else:
            return rows

    def getview(self, mode=None):
        '''
        Return a view of self.df depending on the mode

        Parameters:
             mode (str): define which view shoul be returned:
                  'homozygous' = only homozygous variants
                  'heterozygous' = only heterozygous variants
                  default is None = all variants

        Returns:
             df (Pandas.DataFrame): dataframe view

        '''
        if mode is 'homozygous':
            return self.hom_vars
        if mode is 'heterozygous':
            return self.het_vars
        else:
            return self.df

    def cluster_baf(self, region=None, mode=None, n_clusters=2):
        '''
        cluster BAFs in region applying KMeans.
        return labels, cluster centers, and evaluation score
        for clustering

        Parameters:
             region (string): must stick to the format:
                  chrom:chromStart-chromEnd

             mode (None/str): None, 'heterozygous',  or 'heterozygous'
             n_clusters (number): number of clusters for clustering

        Returns:
             dic_labels (defaultdict): labels dictionary
             centers (list): cluster centers
             score: silhouette score for the resulting clustering
        '''
        dic_labels = defaultdict(list)
        bafs = self.getin(column='baf', region=region, mode=mode)
        x = [[baf] for baf in bafs]

        predictor = KMeans(n_clusters=n_clusters)
        try:
            labels = predictor.fit_predict(x)
        except ValueError:
            print('n_samples={} should be >= n_clusters={}'.format(len(x),
                                                                   n_clusters))
            labels = []
        try:
            score = silhouette_score(x, labels, metric='euclidean')
        except ValueError:
            print('Found  {} clusters'.format(len(labels)))
            print('Cannot compute silhouette_score for {} block!'.format(region))
            score = -1

        for i in range(len(labels)):
            if dic_labels[labels[i]] is None:
                dic_labels[labels[i]] = []
            dic_labels[labels[i]].append(bafs[i])

        try:
            centers = [center[0] for center in predictor.cluster_centers_]
        except AttributeError:
            centers = []

        return dic_labels, centers, score

    def eliminate_vars(self):
        '''
        filter out given variants

        Parameters:
             to_eliminate (list): list of variants to filter out
        '''
        if self.filtered_out_pos is not None:
            initial_len = len(self.df)
            print('Filtering out variants')

            for i, item in enumerate(self.filtered_out_pos):
                self.df = self.df[~((self.df.chrom == item[0]) &
                                    (self.df.pos == item[1]))]
            end_len = len(self.df)
            print('{} sample variants were filtered out!'.format(initial_len -
                                                                 end_len))

    def filter(self, mindp=60):
        '''
        apply filters on variants

        Parameters:
             mindp (number): minimum # reads
        '''
        print('Filtering by: {}'.format('mindp'))
        initial_len = len(self.df)
        self.df = self.df[self.df.dp >= mindp]
        self.split()
        end_len = len(self.df)
        print('{} sample variants were filtered out!'.format(initial_len -
                                                             end_len))

    def split(self):
        '''
        split self.df in two dataframes:
        "self.het_vars" for heterozygous variants and
        "self.hom_vars" for homozygous variants
        '''
        if self.df is not None:
            self.hom_vars, self.het_vars = self.__splitVars()

    @overrides(cdf)
    def iterchroms(self, mode=None, df=None):
        '''
        create groups of variants by chromosome

        Parameters:
             mode (None/str): None, 'heterozygous',  or 'heterozygous'
             df (None/pandas.DataFrame): dataframe to split

        Returns:
             list of dictionaries:
                  {'id': 'chrom:chromStart-chromEnd', 'df': dataframe}
        '''
        if df is None:
            return self.creategroups(self.getview(mode=mode))
        else:
            return self.creategroups(df)

    def vcfplot(self, region=None, mode=None, filename=None,
                auto_open=False):
        '''
        plot variants within region

        Parameters:
             region (string): it must stick to the format:
                  chrom:chromStart-chromEnd
        '''

        # define layout for x and y axis for histogram
        hist_layout = defaultdict(lambda: None)
        hist_layout['y'] = dict(title='# of variants', zeroline=False,
                                range=[0, 1000])
        hist_layout['x'] = dict(zeroline=False)

        # define layout for x and y axis for scatter
        scat_layout = defaultdict(lambda: None)
        scat_layout['y'] = dict(title='BAFs', zeroline=False)
        scat_layout['x'] = dict(zeroline=False)

        if region is None:
            scatter_traces, hist_traces, titles = self.__noregion_plot(mode,
                                                                       scat_layout,
                                                                       hist_layout,
                                                                       auto_open,
                                                                       filename)
        else:
            scatter_traces, hist_traces, titles = self.__region_plot(region,
                                                                     mode,
                                                                     scat_layout,
                                                                     hist_layout,
                                                                     auto_open,
                                                                     filename)

    def __region_plot(self, region, mode, scat_layout, hist_layout,
                      auto_open, filename=None):
        '''
        plot variants within region
        region -- must stick to the format:
        chrom:chromStart-chromEnd
        '''
        hist_traces = []
        scatter_traces = []
        titles = []

        bafs = self.getin(column='baf', region=region,
                          mode=mode)
        pos = self.getin(column='pos', region=region,
                         mode=mode)
        scatter_traces.append([
            y_scatter(bafs, x=pos, size=3, toplot=False,
                      color='rgb(153, 204, 255)')
        ])
        hist_traces.append([
            histogram(bafs, toplot=False, autobinx=False,
                      color='rgb(153, 204, 255)')
        ])
        titles.append(region)
        titles.append('')

        if filename is None:
            filename = '{}/{}-{}.html'.format(self.path2plot,
                                             region.split(':')[0],
                                             self.rootname)
        plot(create_subplots([scatter_traces, hist_traces],
                             titles, layouts=[scat_layout, hist_layout],
                             height=600),
             auto_open=auto_open, filename=filename)
        return scatter_traces, hist_traces, titles

    def __noregion_plot(self, mode, scat_layout, hist_layout,
                        auto_open, filename=None):
        '''
        plot variants when region is None
        '''
        hist_traces = []
        scatter_traces = []
        titles = []

        # define filename
        if filename is None:
            filename = '{}/{}-{}.html'.format(self.path2plot,
                                             'all-chrom',
                                             self.rootname)

        for group in self.iterchroms():
            st, ht, t = self.__region_plot(group['id'], mode, scat_layout,
                                           hist_layout, auto_open)
            scatter_traces.extend(st)
            hist_traces.extend(ht)
            titles.extend(t)
        plot(create_subplots([scatter_traces, hist_traces],
                             titles, layouts=[scat_layout,
                             hist_layout]),
             auto_open=auto_open, filename=filename)
        return scatter_traces, hist_traces, titles

    @overrides(cdf)
    def _compute(self):

        print('Computing stats on variant data')
        columns = ['chrom', 'chromStart', 'chromEnd', 'region',
                   'isbimodal', 'score', 'diff']
        all_blocks = []
        for group in self.iterchroms():
            print('Working on blocks of {}'.format(group['id']))
            for block in self.iterblocks(group, size=self.size,
                                         step=self.step):
                # normality test
                het_bafs = self.getin(column='baf', region=block['id'],
                                      mode='heterozygous')
                bimodal = isbimodal(het_bafs)

                # clustering
                labels, centers, score = self.cluster_baf(region=block['id'],
                                                          mode='heterozygous')
                if centers:
                    centers_diff = abs(centers[0] - centers[-1])
                else:
                    centers_diff = -1
                chrom, chromStart, chromEnd = Region(block['id']).as_tuple
                all_blocks.append([chrom, chromStart, chromEnd, block['id'],
                                   bimodal, score, centers_diff])
        return self.__IQR(DataFrame(all_blocks, columns=columns))

    def __IQR(self, blocks):
        '''
        compute IQR on silhouette score and centers
        difference for each block

        Parameters:
             blocks (list): potential CNV bocks
        '''
        print('Working on clustering (BAF) data')
        all_blocks = []
        columns = list(blocks.columns)
        columns.extend(['isoutlier_1st', 'isoutlier_2nd'])
        # compute IQR metric based on silhouette score and centers difference
        cb_score, ca_score, m_score = compute_metric(blocks, -2, self.metric)
        cb_diff, ca_diff, m_diff = compute_metric(blocks, -1, self.metric)
        ranges = [self.interval_range,
                  self.interval_range * self.cnv_like_range]
        for row in blocks.itertuples():
            block_data = list(row[1:len(row)])
            # apply IQR for each range
            for i, interval_range in enumerate(ranges):
                block_data.append(
                    above_range(row.score, m_score, ca_score, interval_range) and
                    above_range(row.diff, m_diff, ca_diff, interval_range) and
                    row.isbimodal)
            all_blocks.append(block_data)
        return DataFrame(all_blocks, columns=columns)

    @overrides(cdf)
    def _make_subplots(self, cnv, value_column='baf',
                       pos_column='pos', cnvlike=None):
        # define layout for x and y axis
        layout = defaultdict(lambda: None)
        layout['y'] = dict(title='BAFs', zeroline=False,
                           range=[0, 1])
        layout['x'] = dict(zeroline=False)

        traces, titles = self._create_traces(cnv, value_column,
                                             pos_column, cnvlike=cnvlike,
                                             plotting=True, layout=layout)
        return traces, titles, layout


@validstr('filename', empty_allowed=False)
class VCFConfig(object):
    '''
    Detect CNVs based on read B allele frequency (BAF) data
    '''
    def __init__(self, filename, tofilter=True):
        '''
        Parameters:
             filename (str): config file filename
             tofilter (boolean): whether to apply filters on variants
        '''
        self.filename = filename
        self.sections_params = {
            'baseline': 'm',
            'sample': 'm',
            'output': 'm',
            'vartest': 'o'
        }
        # load configfile
        self.config = ConfigfileParser(self.filename,
                                       self.sections_params)
        # load sample test
        sample = VCF(self.config.sections['sample']['bamfile'])
        # load baseline test
        baseline = VCFList(self.config.sections['baseline']['bamfiles'])
        # make test
        if self.config.sections['vartest']:
            self.vcftest = VCFTest(baseline, sample,
                                   **self.config.sections['vartest'],
                                   path=self.config.sections['output']['path'])
        else:
            self.vcftest = VCFTest(baseline, sample,
                                   path=self.config.sections['output']['path'])

        if tofilter:
            self.vcftest.filter()
            self.vcftest.load_filtered_out()
            if self.vcftest.filtered_out_pos is None:
                self.vcftest.compare()
                self.vcftest.save_filtered_out()
            self.vcftest.eliminate_vars()
        self.vcftest.split()
        print('Creating plots at {}'.format(self.vcftest.path2plot))
        self.vcftest.vcfplot()
        filename = '{}/vcftest.csv'.format(self.vcftest.path2table)
        print('Writing table at {}'.format(filename))
        self.vcftest.df.to_csv(filename, sep='\t',  index=False)
        print('Done!')
