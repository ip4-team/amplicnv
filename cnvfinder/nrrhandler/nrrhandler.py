#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: valengo
"""
import sys
from collections import defaultdict
from math import sqrt
from multiprocessing import cpu_count
from typing import Union

import pysam
from numpy import log
from pandas import DataFrame

from cnvfinder.tsvparser import CoverageFileParser
from ..bedloader import ROI
from ..commons import ChromDF as cdf
from ..graphics import scatter
from ..mphandler import MPPoolHandler
from ..stats import above_range
from ..stats import below_range
from ..stats import classify_by_count
from ..stats import compute_metric
from ..stats import filter_by_cutoff
from ..stats import iqr
from ..utils import ConfigfileParser
from ..utils import NumberProperty
from ..utils import Region
from ..utils import appenddir
from ..utils import createdir
from ..utils import overrides
from ..utils import validstr


def readcount(region_list: list, filename: str) -> list:
    """
    Count the number of reads in regions using pysam.AlignmentFile.count method

    :param list region_list: regions in the format: 'chr1:1000-10000'
    :param str filename: path to bamfile
    :return counters: list of number of reads in regions
    """
    try:
        with pysam.AlignmentFile(filename, 'rb') as bamfile:
            counters = []
            for region in region_list:
                try:
                    counters.append(bamfile.count(region=region))
                except ValueError:
                    print("Failed counting the number of reads in {} from {}".format(region, filename))
                    counters.append(None)
            return counters
    except OSError:
        sys.exit("Failed counting the number of reads in region from {}".format(filename))


def attr_from_header(line):
    if line.startswith('#'):
        return line.split(':', 1)[1].strip()
    else:
        print('{0} is not a header line'.format(line))
        return None


def metric2label(counters, q1, q3, metric, interval_range):
    labels = []
    for counter in counters:
        if below_range(counter, metric, q1, interval_range):
            labels.append('-')
        elif above_range(counter, metric, q3, interval_range):
            labels.append('+')
        else:
            labels.append('o')
    return labels


@validstr('bedfile', empty_allowed=True)
@validstr('bamfile', empty_allowed=True)
@validstr('region', empty_allowed=True)
class NRR(object):
    """

    NRR stands for "Number of Reads in Region" loaded from a BAM file

    :param str bedfile: path to bedfile where amplicons are listed in
    :param str bamfile: path to alignment file (bam format)
    :param str region: limit target definition to a given region. It should be in the form: chr1:10000-90000
    :param list counters: list of read depth counters
    :param ROI bed: amplicons already loaded in memory
    :param bool parallel: whether to count target read depth in parallel
    :param bool to_label: whether to label targets regarding each target read depth in comparison to the mean
    :param str covfile: path to amplicon.cov file

    """

    def __init__(self, bedfile: str = None, bamfile: str = None, region: str = None,
                 counters: list = [], bed: Union[ROI, CoverageFileParser] = None, parallel: bool = True,
                 to_label: bool = False, covfile: str = None):
        self.covfile = covfile
        self.bedfile = bedfile
        self.bamfile = bamfile
        self.region = region
        self._counters = self.counters = counters
        self._bed = self.bed = bed
        self.reads_by_pool = defaultdict(int)
        self._nreads = self.nreads = 0
        self.normalized_counters = []
        self.labels = None
        self.labels_by_pool = None

        # load or count rd
        if self.covfile is not None:
            self.bed = CoverageFileParser(self.covfile)
            self.counters = self.bed.counters
        elif self.load(self.bamfile + '.txt') is None:
            self.count(parallel=parallel)
            self.save()

        if self.counters:
            self.reads_by_pool = self.__count_pools()
            self.normalized_counters = self.__norm()

        if len(self.counters) > 0 and to_label:
            print('Labeling targets')
            self.labels = self.__label_targets(mode='log')
            self.labels_by_pool = self.count_label_by_pool()

    @property
    def bed(self):
        return self._bed

    @bed.setter
    def bed(self, value):
        if (value is None and
                self.bedfile is not None and self.covfile is None):
            self._bed = ROI(self.bedfile)
        else:
            self._bed = value

    @property
    def counters(self):
        return self._counters

    @counters.setter
    def counters(self, value):
        self._counters = value

    @property
    def nreads(self):
        return self._nreads

    @nreads.setter
    def nreads(self, value):
        self._nreads = value

    def count(self, cores: int = None, parallel: bool = False):
        """

        For each defined target, count read depth

        :param int cores: number of cores to be used when counting in parallel mode. Default: all available
        :param bool parallel: whether to count read depth in parallel

        """
        if parallel:
            if cores is None:
                cores = cpu_count()
            self.counters = self.__parallel_count(cores)
        else:
            self.counters = self.__count()

    def load(self, filename: str) -> Union[int, None]:
        """

        Load a single NRR from a text file

        :param str filename: path to count file. Normally something like: bamfile.bam.txt
        :return: 1 when data loading is successful, None otherwise
        """
        print('Loading {0} read counters'.format(filename))

        try:
            with open(filename, 'r') as file:
                lines = file.readlines()
        except IOError:
            print('There is no counter file for {0} file'.format(self.bamfile))
            return None

        try:
            self.bamfile = attr_from_header(lines[0])
            self.bedfile = attr_from_header(lines[1])
            self.region = attr_from_header(lines[2])
            if self.region == 'None':
                self.region = None
        except IndexError:
            print('Not enough information on {} header'.format(filename))
            print('Aborting!')
            return None

        print('Extracting read counters from file')
        counters = []
        for line in lines:
            if line.startswith('chr'):
                counters.append(int(line.split('\t')[4].strip()))

        self.counters = counters
        if self.bed is None:
            self.bed = ROI(self.bedfile)

        return 1

    def save(self, filename: str = None):
        """
        Save a single NRR on a text file

        :param str filename: path to output
        """
        print('Saving {0} read counters'.format(self.bamfile))
        if filename is None:
            filename = self.bamfile + '.txt'

        try:
            with open(filename, 'w') as file:
                file.write('# BAM source: {0}\n'.format(self.bamfile))
                file.write('# bed source: {0}\n'.format(self.bedfile))
                file.write('# Regions: {0}\n'.format(self.region))

                t = list(self.bed.targets.itertuples())
                for i in range(len(t)):
                    file.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(t[i][1],
                                                                  t[i][2],
                                                                  t[i][3],
                                                                  t[i][5],
                                                                  self.counters
                                                                  [i]))
            print('Done!')
        except IOError as error:
            print(error)

    def __count(self):
        try:
            bamfile = pysam.AlignmentFile(self.bamfile, 'rb')
        except OSError as error:
            print(error)
            return None

        print('Counting number of reads of {}'.format(
            self.bamfile))
        read_counters = []
        for row in self.bed.targets.itertuples():
            read_counters.append(bamfile.count(row[1], row[2], row[3]))
        return read_counters

    def __parallel_count(self, cores: int):
        print('Counting number of reads of {}, using {} cores'.format(
            self.bamfile, cores))

        # get target regions
        targets = self.bed.targets
        # row[2] + 1 -> faking that we use one based index
        region_list = ['{}:{}-{}'.format(row[1], row[2] + 1, row[3])
                       for row in targets.itertuples()]

        # define chunksize and create chunks
        chunksize, extra = divmod(len(region_list), cores * 4)
        if extra:
            chunksize += 1

        chunks = []
        for i in range(0, len(region_list), chunksize):
            chunks.append(region_list[i:i + chunksize])

        # define tasks
        tasks = [(readcount, chunk, self.bamfile) for chunk in chunks]

        with pysam.AlignmentFile(self.bamfile, 'rb'):
            counters = MPPoolHandler(tasks, cores).run()

        if counters:
            counters = [c for counter in counters for c in counter]
        return counters

    def __count_pools(self):
        """
        Count number of reads per pool.
        This method is used when loading a NRR from a file
        """
        print('Counting reads by pools')
        targets = list(self.bed.targets.itertuples())
        if len(targets) != len(self.counters):
            print('Number of targets and their read counters differ.')
            print('Aborting!')
            return None

        reads_by_pool = defaultdict(int)
        for i, target in enumerate(targets):
            pools = target.pools.unique_flattened()
            for pool in pools:
                reads_by_pool[pool] += self.counters[i] / len(pools)
        self.nreads = sum(self.counters)

        return reads_by_pool

    def __norm(self, mag: int = 1000000):
        print('Normalizing counters')
        normalized = []
        targets = list(self.bed.targets.itertuples())
        if len(targets) != len(self.counters):
            print('Number of targets and their read counters differ.')
            print('Aborting!')
            return None

        for i, target in enumerate(targets):
            current_pools_counter = []
            pools = target.pools.unique_flattened()
            for pool in pools:
                current_pools_counter.append((self.counters[i] /
                                              self.reads_by_pool[pool]))
            normalized.append(mag * (sum(current_pools_counter) /
                                     len(pools)))

        return normalized

    def __label_targets(self, iqr_range: float = 1.5, std_range: float = 1.5, mode: str = 'normalized') -> list:
        """
        Label targets considering IQR on counters (normalized or not - user's choice) for discovering whether a target
        is in the lower (-) quartile, upper (+) quartile, or middle (o)

        :param float iqr_range: value to multiply IQR by
        :param float std_range: value to multiply std by
        :param str mode: 'normalized' for employing self.normalized_counters; 'log' for log(self.normalized_counters)
        :return: labels (list) a list of size (len.counters) representing the labels for the targets

        """
        # compute metric (IQR)
        if mode == 'normalized':
            counters = self.normalized_counters
        elif mode == 'log':
            counters = log(self.normalized_counters)
        else:
            counters = self.counters
        df = DataFrame(counters)
        q1, q3, metric = iqr(df, 0)
        # label targets considering interquartile range (IQR)
        labels = metric2label(counters, q1, q3, metric, iqr_range)
        # filter out counters + or - labeled
        df = DataFrame(labels)
        df.loc[:, len(df.columns)] = counters
        filtered = df[df[0] == 'o']
        # label targets considering std
        mean = filtered.iloc[:, 1].mean()
        q1, q3, metric = compute_metric(filtered, 1, 'std', center=mean)
        labels = metric2label(counters, q1, q3, metric, std_range)
        return labels

    def count_label_by_pool(self) -> DataFrame:
        """
        Count the number of targets arranged by label in each pool

        :return: ldf number of targets considering pools x labels
        """
        df = DataFrame(self.bed.targets.iloc[:, 4])
        df[len(df.columns)] = self.labels

        counters = defaultdict(defaultdict)
        for label in df[1].unique():
            counters[label] = defaultdict(lambda: 0)

        for row in df.itertuples():
            for pool in row[1].pools:
                counters[row[2]][pool] += 1

        ldf = DataFrame.from_dict(counters)
        ldf.reset_index(inplace=True)
        ldf.rename(columns={'index': 'pool'}, inplace=True)
        return ldf


@validstr('bedfile', empty_allowed=True)
@validstr('region', empty_allowed=True)
class NRRList(object):
    """
    NRR object list management

    :param str bedfile: path to file (BED) where the amplicons are listed in
    :param list bamfiles: list of paths to alignment (BAM) files
    :param str region: limit target definition to a given region. It should be in the form: chr1:10000-90000
    :param ROI bed: amplicons already loaded into memory
    :param bool parallel: whether to count defined targets read depth in parallel
    :param bool to_classify: whether to classify defined targets regarding their read depth
    :param list covfiles: list of paths to amplicon.cov files
    """

    def __init__(self, bedfile: str = None, bamfiles: list = None, region: str = None,
                 bed: ROI = None, parallel: bool = True, to_classify: bool = False, covfiles: list = None):
        self.bedfile = bedfile
        self._bamfiles = self.bamfiles = bamfiles
        self._covfiles = self.covfiles = covfiles
        self.region = region
        self.bed = bed
        self.list = []
        self.mean = []
        self.sd = []
        self.median = []
        self.normalized_mean = []
        # caches counters
        self.__counters = []
        self.__normalized_counters = []
        self.labels = None

        if self.covfiles:
            for i, coverage_filename in enumerate(self.covfiles):
                self.list.append(NRR(region=region,
                                     covfile=coverage_filename,
                                     to_label=to_classify))
                if self.list[i].counters:
                    self.__counters.append(self.list[i].counters)
                    self.__normalized_counters.append(self.list[i].
                                                      normalized_counters)
        elif self.bamfiles:
            for i, bamfile in enumerate(self.bamfiles):
                self.list.append(NRR(bedfile=bedfile,
                                     bamfile=bamfile,
                                     bed=bed,
                                     region=region,
                                     parallel=parallel,
                                     to_label=to_classify))

                bed = self.list[i].bed
                if self.list[i].counters:
                    self.__counters.append(self.list[i].counters)
                    self.__normalized_counters.append(self.list[i].
                                                      normalized_counters)
        if to_classify:
            print('Classifying targets')
            self.df, self.labels = self.__label_targets()

    @property
    def bamfiles(self):
        return self._bamfiles

    @bamfiles.setter
    def bamfiles(self, value):
        if value is not None and not isinstance(value, list):
            self._bamfiles = [value]
        else:
            self._bamfiles = value

    def count(self):
        for nrr in self.list:
            nrr.count()
            if nrr.counters:
                self.__counters.append(nrr.counters)
                self.__normalized_counters.append(nrr.normalized_counters)

    def save(self):
        for nrr in self.list:
            nrr.save(filename=nrr.bamfile + '.txt')

    def make_mean(self):
        """
        Compute baseline read depth mean for each defined target
        """
        self.normalized_mean = []
        self.sd = []
        if self.__counters and len(self.__counters) == len(self.list):
            try:
                self.mean = [sum(c) / len(c) for c in zip(*self.__counters)]

                for counters in zip(*self.__normalized_counters):
                    mean = sum(counters) / len(counters)
                    deviations = [(c - mean) ** 2 for c in counters]
                    variance = sum(deviations) / len(counters)
                    sd = sqrt(variance)
                    self.normalized_mean.append(mean)
                    self.sd.append(sd)
            except TypeError:
                print('Amount of read counters (baseline and sample) differ')
                print('Are their region targets different? Aborting!')

        else:
            for l in self.list:
                if not l.counters:
                    if l.bamfile:
                        print('Not possible to calculate mean when ' +
                              '{0} counter list is missing!'.format(l.bamfile))
                    else:
                        print('Not possible to calculate mean when ' +
                              '{0} counter list is missing!'.format(l.bamfile + '.txt'))

    @staticmethod
    def _add_pools(targets):
        """
        Add pools info to the df in order to verify whether +, -, or o
        amplicons are pool related
        """
        pools = [target[5].pools for target in targets.itertuples()]
        return pools

    def __label_targets(self):
        # targets
        targets = self.list[0].bed.targets
        # make a copy of targets df base (chrom:chromStart-chromEnd)
        df = targets.iloc[:, :4].copy()
        # add pools into table
        df.loc[:, len(df.columns)] = self._add_pools(targets)
        for sample in self.list:
            df.loc[:, len(df.columns)] = sample.labels
        filtered_df = df.iloc[:, 5:len(df.columns)]
        return df, filtered_df.apply(classify_by_count, axis=1)


@validstr('configfile', empty_allowed=True)
@validstr('path', empty_allowed=False)
class NRRTest(cdf):
    """
    Hold information about tests between a NRR baseline (NRRList) and a NRR test sample

    :param NRRList baseline: represents bamfiles of the baseline
    :param NRR sample: represents the bamfile of a test sample
    :param str path: output directory path
    :param int size: block's size when sliding window
    :param int step: step's size when sliding window
    :param str metric: param used to define which metric should be used 'std' or 'IQR'
    :param float interval_range: value to multiply metric by
    :param int minread: minimum number of reads used to filter targets
    :param float below_cutoff: filter out data (ratios) below this cutoff
    :param float above_cutoff: filter out data (ratios) above this cutoff
    :param int maxdist: maximum distance allowed of a cnv-like block, to its closest cnv block, for it be a cnv as well
    :param float cnv_like_range: value to multiply interval_range by in order to detect cnv-like
    :param int bins: number of bins to use when plotting ratio data
    :param str method: method used in order to group rations when plotting
    """

    def __init__(self, baseline: NRRList, sample: NRR, path: str = 'results',
                 size: int = 200, step: int = 10, metric: str = 'std', interval_range: float = 3,
                 minread: int = 30, below_cutoff: float = 0.7, above_cutoff: float = 1.3,
                 maxdist: int = 15000000, cnv_like_range: float = 0.7,
                 bins=500, method='chr_group'):
        super().__init__(None)
        self._baseline = self.baseline = baseline
        self._sample = self.sample = sample
        self.ratios = []
        self.path = path
        self.size = size
        self.step = step
        self.metric = metric
        self.interval_range = interval_range
        self.minread = minread
        self.below_cutoff = below_cutoff
        self.above_cutoff = above_cutoff
        self.maxdist = maxdist
        self.cnv_like_range = cnv_like_range
        self.method = method
        self.bins = bins
        self.columns = ['chrom',
                        'chromStart',
                        'chromEnd',
                        'gene',
                        'amplicon',
                        'counter',
                        'mean',
                        'norm_counter',
                        'norm_mean',
                        'sd',
                        'ratio']
        self.rootname = 'ratios'
        self.path2plot = appenddir(self.path, 'plots/bam/ratios')
        self.path2table = appenddir(self.path, 'tables/bam')
        self.path2plotcnv = appenddir(self.path, 'plots/bam/cnv')
        self.normal_ratio = 1
        createdir(self.path2plot)
        createdir(self.path2plotcnv)
        createdir(self.path2table)

    # properties
    bins = NumberProperty('_bins')

    minread = NumberProperty('_minread')

    below_cutoff = NumberProperty('_below_cutoff')

    above_cutoff = NumberProperty('_above_cutoff')

    @property
    def baseline(self):
        return self._baseline

    @baseline.setter
    def baseline(self, value):
        if isinstance(value, NRRList):
            self._baseline = value
        else:
            print('Baseline type must be NRRList and not {0}'.format(
                type(value)))

    @property
    def sample(self):
        return self._sample

    @sample.setter
    def sample(self, value):
        if isinstance(value, NRR):
            self._sample = value
        else:
            print('Sample type must be NRR and not {0}'.format(
                type(value)))

    @overrides(cdf)
    def _createdf(self):
        """
        create NRRTest dataframe with sample, baseline, and
        test data
        """
        df = self.sample.bed.targets.copy()
        df.loc[:, len(df.columns)] = self.sample.counters
        df.loc[:, len(df.columns)] = self.baseline.mean
        df.loc[:, len(df.columns)] = self.sample.normalized_counters
        df.loc[:, len(df.columns)] = self.baseline.normalized_mean
        df.loc[:, len(df.columns)] = self.baseline.sd
        df.loc[:, len(df.columns)] = self.ratios
        df.columns = self.columns
        self.df = df

    def make_ratio(self):
        """
        Compute ratio between the standardized read depth values from a test sample
        and the mean baseline read depth
        """
        print('Calculating reference bam files # reads mean')
        self.baseline.make_mean()

        print('Calculating # reads test sample/ baseline')
        ratios = []

        try:
            for i in range(len(self.sample.counters)):
                if self.baseline.mean[i] > self.minread:
                    ratios.append(self.sample.normalized_counters[i] /
                                  self.baseline.normalized_mean[i])
                else:
                    if ratios:
                        ratios.append(ratios[-1])
                    else:
                        ratios.append(1)

            self.ratios = ratios
            if self.ratios:
                self._createdf()
        except IndexError:
            print('Amount of read counters (baseline and sample) differ')
            print('Are their region targets different? Aborting!')
        except TypeError:
            print('Amount of read counters (baseline and sample) differ')
            print('Are their region targets different? Aborting!')

    def plot(self):
        """
        Plot all defined target ratios
        """
        if self.ratios:
            if self.method == 'chr_group':
                xs, ys, offsets, names = self.__regroup(self.bins)
                scatter(xs, ys, yax=[0, 5], offset=offsets,
                        name=names,
                        filename='{}/{}-{}.html'.format(self.path2plot,
                                                        'all-chrom',
                                                        self.rootname))
            else:
                x, y = self.__remean(self.ratios, self.bins, self.method)
                scatter(x, y, yax=[0, 5],
                        filename='{}/{}-{}.html'.format(self.path2plot,
                                                        'all-chrom',
                                                        self.rootname))
        else:
            print('There is no data to plot!')
            print('You should call NRRTest.make_ratio() ' +
                  'before plotting something!')

    def filter(self, df: DataFrame) -> DataFrame:
        """
        Filter dataframe by mean >= minread

        :param DataFrame df: dataframe to be filtered
        :return: filtered dataframe
        """
        return df[df.loc[:, 'mean'] >= self.minread]

    @staticmethod
    def __remean(ratios, bins, method):

        if isinstance(bins, str):
            bins = int(bins)

        if len(ratios) < bins:
            bins = len(ratios)

        if method == 'blocks':
            xs = [i * bins for i in range(len(ratios) // bins)]
        else:
            xs = [i for i in range(len(ratios) - bins + 1)]
        ys = [sum(ratios[x:x + bins]) / bins for x in xs]
        return xs, ys

    def __regroup(self, bins):

        if isinstance(bins, str):
            bins = int(bins)

        current_offset = 0
        offsets = []
        xs = []
        ys = []
        groups = self.iterchroms()
        unique_chroms = self.df.chrom.unique()

        for group in groups:
            offsets.append(current_offset)
            current_offset += len(group['df'])
            serie = [row[-1] for row in group['df'].itertuples()]
            x, y = self.__remean(serie, bins, 'blocks')
            xs.append(x)
            ys.append(y)

        return xs, ys, offsets, unique_chroms

    def summary(self, npools: int = 12) -> DataFrame:
        """
        Create a summary of NRRTest samples (baseline and test)

        :param int npools: number of pools
        :return: summary as a dataframe
        """
        # define column names
        columns = ['ntotal']
        poolnames = ['pool_{}'.format(i) for i in range(1, npools + 1)]
        columns.extend(poolnames)
        index = []

        # select data from dataframe
        n_totals = []
        for sample in self.baseline.list:
            n_sample = [sample.nreads]
            index.append(sample.bamfile.split('/')[-1])
            for i in range(1, npools + 1):
                n_sample.append(sample.reads_by_pool[i])
            n_totals.append(n_sample)
        df = DataFrame(n_totals, columns=columns, index=index)
        print(df.describe())
        return df

    def merge(self, potential_cnvs: list) -> list:
        """
        Merge CNV blocks

        :param list potential_cnvs: list of potential CNV blocks
        :return: list of CNV blocks merged
        """
        cnv_dict = defaultdict(list)
        for p_cnv in potential_cnvs:
            for cnv in p_cnv.itertuples():
                region = '{}:{}-{}'.format(cnv[1], cnv[2], cnv[3])
                targets = self.getin(region)
                for target in targets.itertuples():
                    cnv_dict['{}:{}-{}'.format(target.chrom,
                                               target.chromStart,
                                               target.chromEnd)] = [
                        target.chrom,
                        target.chromStart,
                        target.chromEnd,
                        target.gene, target[-1]]
        cnvs = [value for k, value in cnv_dict.items()]
        cnvs.sort(key=lambda x: (x[0], x[1], x[2]))
        return cnvs

    @overrides(cdf)
    def _make_subplots(self, cnv, value_column='ratio',
                       pos_column='chromStart', cnvlike=None):
        # define layout for x and y axis
        layout = defaultdict(lambda: dict)
        layout['y'] = dict(title='Ratios', zeroline=False,
                           range=[0, 3])
        layout['x'] = dict(zeroline=False)

        traces, titles = self._create_traces(cnv, value_column,
                                             pos_column, cnvlike=cnvlike,
                                             toplot=True, layout=layout)
        return traces, titles, layout

    @overrides(cdf)
    def _compute(self):
        print('Computing stats on ratio data')

        if not self.ratios:
            self.make_ratio()

        all_blocks = []
        columns = ['chrom', 'chromStart', 'chromEnd', 'region', 'median',
                   'isoutlier_1st', 'isoutlier_2nd']

        # filter targets by minread and ratios (below and above cutoff)
        filtered_targets = filter_by_cutoff(self.filter(self.get_autosomal()),
                                            -1, self.below_cutoff,
                                            self.above_cutoff)
        # compute and use metric to find blocks of outliers for each range
        cb, ca, m = compute_metric(filtered_targets, -1, self.metric)
        ranges = [self.interval_range, self.interval_range * self.cnv_like_range]
        for group in self.iterchroms():
            print('Working on blocks of {}'.format(group['id']))

            for block in self.iterblocks(group, size=self.size, step=self.step):
                med_ratio = self.filter(block['df']).loc[:, 'ratio'].median()
                chrom, chrom_start, chrom_end = Region(block['id']).as_tuple
                block_data = [chrom, chrom_start, chrom_end, block['id'], med_ratio]
                for i, interval_range in enumerate(ranges):
                    block_data.append(above_range(med_ratio, m, ca, interval_range) or
                                      below_range(med_ratio, m, cb, interval_range))
                all_blocks.append(block_data)
        return self._call(DataFrame(all_blocks, columns=columns))

    @overrides(cdf)
    def _call(self, df):
        df['call'] = df['median'].apply(lambda x: 'loss' if x < self.normal_ratio else 'gain')
        return df


@validstr('filename', empty_allowed=False)
class NRRConfig(object):
    """
    Detect CNVs based on read depth data

    :param str filename: path to configuration file

    """

    def __init__(self, filename: str):
        self.filename = filename
        self.sections_params = {
            'baseline': 'm',
            'bed': 'm',
            'sample': 'm',
            'output': 'm',
            'targtest': 'o'
        }
        # load configfile
        self.config = ConfigfileParser(self.filename,
                                       self.sections_params)
        # load sample test
        sample = NRR(bamfile=self.config.sections['sample']['bamfile'],
                     covfile=self.config.sections['sample']['covfile'],
                     bedfile=self.config.sections['bed']['bedfile'])
        # load baseline test
        baseline = NRRList(bamfiles=self.config.sections['baseline']['bamfiles'],
                           covfiles=self.config.sections['baseline']['covfiles'],
                           bedfile=self.config.sections['bed']['bedfile'])
        # make test
        if self.config.sections['targtest']:
            self.nrrtest = NRRTest(baseline, sample,
                                   **self.config.sections['targtest'],
                                   path=self.config.sections['output']['path'])
        else:
            self.nrrtest = NRRTest(baseline, sample,
                                   path=self.config.sections['output']['path'])

        self.nrrtest.make_ratio()
        if self.nrrtest.ratios:
            print('Creating plots at {}'.format(self.nrrtest.path2plot))
            self.nrrtest.plot()
            filename = '{}/nrrtest.csv'.format(self.nrrtest.path2table)
            print('Writing table at {}'.format(filename))
            self.nrrtest.df.to_csv(filename, sep='\t', index=False)
            print('Done!')
