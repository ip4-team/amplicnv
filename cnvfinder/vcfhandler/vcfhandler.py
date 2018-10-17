#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: valengo
"""
from typing import Union

from pysam import VariantFile
from collections import defaultdict
from ..utils import Region
from pandas import DataFrame, Series
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
    """
    **Hold information on variants loaded from a VCF file using pysam.VariantFile**

    :param str vcffile: path to vcf file

    """

    def __init__(self, vcffile: str):
        self.vcffile = vcffile
        if self.vcffile.endswith('bam'):
            self.vcffile = em(self.vcffile, 'bam', 'vcf.gz').new_filename
        self.variants = self.__load_vars()

    def __load_vars(self):
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
        except OSError:
            print('File {} not found!'.format(self.vcffile))
            return None
        except AttributeError:
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

    @staticmethod
    def calc_baf(info: Series) -> float:
        """
        Compute the frequency of the first alternative allele

        :param Series info: single variant data
        :return: computed B-allele frequency
        """
        baf = info['AO'][0] / (info['AO'][0] + info['RO'])

        return baf


class VCFList(object):
    """
    **VCF object list management**

    :param list vcffiles: list of paths to variant (VCF) files
    """

    def __init__(self, vcffiles: list):
        self.vcffiles = vcffiles
        self.list = []

        for i, vcffile in enumerate(vcffiles):
            if vcffile.endswith('bam'):
                self.vcffiles[i] = em(vcffile, 'bam', 'vcf.gz').new_filename
            self.list.append(VCF(self.vcffiles[i]))


@validstr('configfile', empty_allowed=True)
class VCFTest(cdf):
    """
    **Hold information about tests between a VCF baseline (VCFList) and a VCF test sample**

    :param VCFList baseline: represents baseline's VCF files
    :param VCF sample: represents sample's VCF file
    :param str metric: param used to define which metric should be used. For instance, only 'IQR' is available for BAF
    :param float interval_range: value to multiply metric by
    :param int size: block's size when sliding window
    :param int step: block's size when sliding window
    :param float cnv_like_range: value to multiply interval_range by in order to detect cnv-like
    :param int maxdist: maximum distance allowed of a cnv-like block, to its closest cnv block, for it be a cnv as well
    :param str path: output directory path
    """

    def __init__(self, baseline: VCFList, sample: VCF, metric: str = 'IQR',
                 interval_range: float = 1.5, size: int = 400, step: int = 40,
                 cnv_like_range: float = 0.7, maxdist: int = 15000000,
                 path: str = 'results'):
        super().__init__(None)
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
        self.eliminated = None
        self.filtered_out = None
        self.filtered_out_pos = None
        self.hom_vars, self.het_vars = None, None
        createdir(self.path2plot)
        createdir(self.path2table)
        createdir(self.path2plotcnv)
        self._createdf()

    def compare(self, nmin: int = None, maxdiff: float = 0.00001):
        """
        Apply filters on sample variants considering variants from baseline files. Current filters:

            - **similar ratios**: given a variant, filter it out if BAF values are identical considering all or 'nmin' samples

        :param int nmin: minimum number of baseline samples to one variant to be considered a false positive
        :param float maxdiff: max difference between numbers for them to be equal
        """
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

    @staticmethod
    def __similar_ratios(row: Series, copy2compare: DataFrame, maxdiff: float) -> int:
        similar_vars = copy2compare.loc[(copy2compare['pos'] == row.pos) &
                                        (copy2compare['chrom'] == row.chrom)]
        if not similar_vars.empty:
            for variant in similar_vars.itertuples():
                try:
                    diff = abs(variant.baf - row.baf)
                    if diff < maxdiff:
                        return 1
                except KeyError:
                    print(variant.info.keys(), variant.pos)
                    print(row.info.keys(), row.pos)
        return 0

    def save_filtered_out(self):
        """
        Save filtered out variants of the test sample. Output file is defined as: self.vcffile + '.txt'
        """
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
        """
        Load filtered out variants from file. Path to file is defined as: self.vcffile + '.txt'
        """
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
        except FileNotFoundError:
            print('"{}.txt" file not found!'.format(self.sample.vcffile))
            print('In this case, you could run VCFTest.compare()')
            self.filtered_out = None
            self.filtered_out_pos = None

    def __split_vars(self, homozygous_freq: float = 0.90) -> tuple:
        """
        From self.df, this method creates two new dataframes grouped by:
            - homozygous variants
            - heterozygous variants

        :param float homozygous_freq: relative homozygous frequency
        :return: split dataframe
        """
        df = self.df
        return df[df.baf >= homozygous_freq], df[df.baf < homozygous_freq]

    @overrides(cdf)
    def getin(self, region: str = None, column: str = None, mode: str = None) -> Union[list, None, DataFrame]:
        """
        Get data in 'column' for variants that are located in region considering 'mode'. Available columns:

        - chrom
        - start
        - stop
        - alleles
        - alts
        - contig
        - id
        - info
        - pos
        - qual
        - ref
        - rid
        - rlen
        - dp
        - baf

        :param str region: it should be in the form: chr1:10000-90000
        :param str column: column's name
        :param mode: define which view should be returned: homozygous or heterozygous
        :return: dataframe containing the view, or a list if 'column' is passed
        """
        df = self.getview(mode=mode)
        rows = self.getrows(df, region)
        if rows is not None and column:
            return list(rows.loc[:, column])
        else:
            return rows

    def getview(self, mode: str = None) -> DataFrame:
        """
        Create a view of self.df depending on the mode. Available modes:

        - homozygous: create view with homozygous variants
        - heterozygous: create a view with heterozygous variants

        :param str mode: view mode
        :return: dataframe's view or dataframe, if mode = None
        """
        if mode is 'homozygous':
            return self.hom_vars
        if mode is 'heterozygous':
            return self.het_vars
        else:
            return self.df

    def cluster_baf(self, region: str = None, mode: str = None, n_clusters: int = 2) -> tuple:
        """
        Cluster BAF values in 'region', considering 'mode', and applying the KMeans method.

        :param str region: it should be in the form: chr1:10000-90000
        :param mode: homozygous or heterozygous
        :param int n_clusters: number of expected clusters
        :return: labels (defaultdict), cluster centers (list), and silhouette score (float)
        """
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
        """
        Filter out the given variants from self.df
        """
        if self.filtered_out_pos is not None:
            initial_len = len(self.df)
            print('Filtering out variants')

            for i, item in enumerate(self.filtered_out_pos):
                self.df = self.df[~((self.df.chrom == item[0]) &
                                    (self.df.pos == item[1]))]
            end_len = len(self.df)
            print('{} sample variants were filtered out!'.format(initial_len -
                                                                 end_len))

    def filter(self, mindp: int = 60):
        """
        Apply filters on variants located at self.df

        :param int mindp: minimum number of reads to pass the filter
        """
        print('Filtering by: {}'.format('mindp'))
        initial_len = len(self.df)
        self.df = self.df[self.df.dp >= mindp]
        self.split()
        end_len = len(self.df)
        print('{} sample variants were filtered out!'.format(initial_len -
                                                             end_len))

    def split(self):
        """
        Split self.df in two dataframes

        - **self.het_vars** for heterozygous variants and
        - **self.hom_vars** for homozygous variants
        """
        if self.df is not None:
            self.hom_vars, self.het_vars = self.__split_vars()

    @overrides(cdf)
    def iterchroms(self, mode: str = None, df: DataFrame = None) -> list:
        """
        Group variants by chromosome

        :param str mode: homozygous or heterozygous
        :param Dataframe df: dataframe to split
        :return: list of dictionary groups {'id': 'chrom:chromStart-chromEnd', 'df': dataframe}
        """
        if df is None:
            return self.create_groups(self.getview(mode=mode))
        else:
            return self.create_groups(df)

    def vcfplot(self, region: str = None, mode: str = None, filename: str = None,
                auto_open: bool = False):
        """
        Plot variants within region

        :param region: it should be in the form: chr1:10000-90000
        :param str mode: homozygous or heterozygous
        :param str filename: path to output file
        :param bool auto_open: whether to automatically open the resulting plot
        """
        # define layout for x and y axis for histogram
        hist_layout = defaultdict(lambda: dict)
        hist_layout['y'] = dict(title='# of variants', zeroline=False,
                                range=[0, 1000])
        hist_layout['x'] = dict(zeroline=False)

        # define layout for x and y axis for scatter
        scat_layout = defaultdict(lambda: dict)
        scat_layout['y'] = dict(title='BAFs', zeroline=False)
        scat_layout['x'] = dict(zeroline=False)

        if region is None:
            self.__noregion_plot(mode,
                                 scat_layout,
                                 hist_layout,
                                 auto_open,
                                 filename)
        else:
            self.__region_plot(region,
                               mode,
                               scat_layout,
                               hist_layout,
                               auto_open,
                               filename)

    def __region_plot(self, region: str, mode: str, scat_layout: defaultdict, hist_layout: defaultdict,
                      auto_open: bool, filename: str = None) -> tuple:
        """
        Plot variants within region

        :param str region:  it should be in the form: chr1:10000-90000
        :param str mode: heterozygous or heterozygous
        :param defaultdict scat_layout: scatter layout
        :param defaultdict hist_layout: histogram layout
        :param bool auto_open: whether to auto open plot on web browser
        :param str filename: path to output file
        :return: scatter traces, histogram traces, and their titles
        """
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

    def __noregion_plot(self, mode: str, scat_layout: defaultdict, hist_layout: defaultdict,
                        auto_open: bool, filename: str = None) -> tuple:
        """
        Plot variants when region is None

        :param str mode: heterozygous or heterozygous
        :param defaultdict scat_layout: scatter layout
        :param defaultdict hist_layout: histogram layout
        :param bool auto_open: whether to auto open plot on web browser
        :param str filename: path to output file
        :return: scatter traces, histogram traces, and their titles
        """
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
    def _compute(self) -> DataFrame:
        """
        Compute stats on B-allele frequency data. Firstly, apply "isbimodal" method in order to detect whether
        BAF values within blocks have evidence of a bimodal distribution. Secondly, BAF values are clustered generating
        generating labels, cluster centers, and a corresponding silhouette score. The absolute difference among
        the cluster's centers and the score are evaluated using IQR analysis.

        :return: resulting dataframe
        """

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
                chrom, chrom_start, chrom_end = Region(block['id']).as_tuple
                all_blocks.append([chrom, chrom_start, chrom_end, block['id'],
                                   bimodal, score, centers_diff])
        return self.__iqr(DataFrame(all_blocks, columns=columns))

    def __iqr(self, blocks: DataFrame) -> DataFrame:
        """
        For each block, compute IQR on silhouette score and the absolute difference among the cluster's centers

        :param Dataframe blocks: data
        :return: dataframe containing IQR analysis results
        """
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
        return self._call(DataFrame(all_blocks, columns=columns))

    @overrides(cdf)
    def _make_subplots(self, cnv: DataFrame, value_column: str = 'baf',
                       pos_column: str = 'pos', cnvlike: DataFrame = None) -> tuple:
        """
        Create subplots for showing analysis results for each chromosome

        :param DataFrame cnv: potential CNV
        :param str value_column: column to lookup for Y axis data
        :param str pos_column: column to lookup for X axis data
        :param DataFrame cnvlike: potential CNV-like
        :return:
        """
        # define layout for x and y axis
        layout = defaultdict(lambda: dict)
        layout['y'] = dict(title='BAFs', zeroline=False,
                           range=[0, 1])
        layout['x'] = dict(zeroline=False)

        traces, titles = self._create_traces(cnv, value_column,
                                             pos_column, cnvlike=cnvlike,
                                             toplot=True, layout=layout)
        return traces, titles, layout

    @overrides(cdf)
    def _call(self, df: DataFrame) -> DataFrame:
        df['call'] = df['isbimodal'].apply(lambda x: 'gain' if True else 'loss')
        return df


@validstr('filename', empty_allowed=False)
class VCFConfig(object):
    """
    **Detect CNVs based on read B-allele frequency (BAF) data**

    :param str filename: path to configuration file
    :param bool tofilter: whether to apply filters on variants
    """

    def __init__(self, filename: str, tofilter: bool = True):
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
        if self.config.sections['sample']['vcffile']:
            sample = VCF(self.config.sections['sample']['vcffile'])
        else:
            print('\'vcffile\' option path not defined for sample\'s VCF file. Using bamfile rootname instead.')
            sample = VCF(self.config.sections['sample']['bamfile'])

        # load baseline test
        if self.config.sections['baseline']['vcffiles']:
            baseline = VCFList(self.config.sections['baseline']['vcffiles'])
        else:
            print('\'vcffile\' option path not defined for baseline\'s VCF files. Using bamfile rootname instead.')
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
        self.vcftest.df.to_csv(filename, sep='\t', index=False)
        print('Done!')
