#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: valengo
"""

import pysam
from ..graphics import scatter
from pandas import DataFrame
from ..bedloader import ROI
from ..utils import validstr
from collections import defaultdict
from math import sqrt
from ..utils import Region
from ..utils import createdir
from ..utils import NumberProperty
from ..utils import ConfigfileParser
from ..utils import overrides
from ..commons import ChromDF as cdf
from ..stats import compute_metric
from ..stats import above_range
from ..stats import below_range
from ..stats import filter_by_cutoff
from ..stats import IQR
from ..stats import classify_by_count
from multiprocessing import cpu_count
from ..mphandler import MPPoolHandler
from numpy import log
from ..utils import appenddir


def readcount(region_list, filename):
    '''
    count the number of reads in region using
    pysam.AlignmentFile.count method

    Parameters:
         region_list (list): list of regions in format -> 'chr1:1000-10000'
         filename (str): bamfile from where to count

    Returns:
         counter (number/None): counter number or None in case it fails
    '''
    try:
        with pysam.AlignmentFile(filename, 'rb') as bamfile:
            counters = []
            for region in region_list:
                try:
                    counters.append(bamfile.count(region=region))
                except ValueError as error:
                    print(error)
                    counters.append(None)
            return counters
    except OSError as error:
        print(error)
        return None


@validstr('bedfile', empty_allowed=True)
@validstr('bamfile', empty_allowed=True)
@validstr('region', empty_allowed=True)
class NRR(object):
    '''
    NRR stands for "Number of Reads in Region" loaded from a BAM file.
    '''
    def __init__(self, bedfile=None, bamfile=None, region=None,
                 counters=[], bed=None, parallel=True):
        '''
        Parameters:
             bedfile (str): filename (BED) where the amplicons are listed in
             bamfile (str): aligment file name (bam format)

             region (str): region limit for analysis (view). It must stick to
             the format: chrom:chromStart-chromEnd

             counters (list): list of read depth counters
             bed (ROI): list of amplicons already loaded into memory
             parallel: whether count target read depth in parallel

        '''
        self.bedfile = bedfile
        self.bamfile = bamfile
        self.region = region
        self.counters = counters
        self.bed = bed
        self.reads_by_pool = defaultdict(int)
        self.nreads = 0
        self.normalized_counters = []
        self.labels = None
        self.labels_by_pool = None

        if self.load(self.bamfile + '.txt') is None:
            self.count(parallel=parallel)
            self.save()

        if len(self.counters) > 0:
            print('Labeling targets')
            self.labels = self.__label_targets(mode='log')
            self.labels_by_pool = self.count_label_by_pool()

    @property
    def bed(self):
        return self._bed

    @bed.setter
    def bed(self, value):
        if (value is None and
                self.bedfile is not None):
            self._bed = ROI(self.bedfile, self.region)
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

    def count(self, cores=None, parallel=False):
        '''
        Count read depth, in the alignment file (bamfile), for each target
        that was loaded from the BED file

        Parameters:
             cores (int): number of cores to be used in case parallel is True.
             If it's no specified, multiprocessing.cpu_count() is used in order
             to verify how many processors are in the machine.

             parallel (boolean): whether to count in parallel
        '''
        if parallel:
            if cores is None:
                cores = cpu_count()
            self.counters = self.__parallel_count(cores)
        else:
            self.counters = self.__count()
        if self.counters and self.bed:
            self.reads_by_pool = self.__count_pools()
            self.normalized_counters = self.__norm()

    def load(self, filename):
        '''
        Load a single NRR from a text file

        Parameters:
             filename (str): count file filename. Normally something like:
                  bamfile.bam.txt

        Returns:
             1: when data loading is successful
             None: when it fails

        '''
        print('Loading {0} read counters'.format(filename))

        try:
            with open(filename, 'r') as file:
                lines = file.readlines()
        except IOError as error:
            print('There is no counter file for {0} file'.format(self.bamfile))
            return None

        try:
            self.bamfile = self.__attr_from_header(lines[0])
            self.bedfile = self.__attr_from_header(lines[1])
            self.region = self.__attr_from_header(lines[2])
            if self.region == 'None':
                self.region = None
        except IndexError as error:
            print('Not enough information on file header.')
            print('Aborting!')
            return None

        print('Extracting read counters from file')
        counters = []
        for line in lines:
            if line.startswith('chr'):
                counters.append(int(line.split('\t')[4].strip()))

        self.counters = counters
        if self.bed is None:
            self.bed = ROI(self.bedfile, self.region)
        if self.counters and self.bed:
            self.reads_by_pool = self.__count_pools()
            self.normalized_counters = self.__norm()

        return 1

    def save(self, filename=None):
        '''
        Save a single NRR on a text file

        Parameters:
             filename (str): the name of the file
        '''
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

    def __attr_from_header(self, line):
        if line.startswith('#'):
            return line.split(':', 1)[1].strip()
        else:
            print('{0} is not a header line'.format(line))
            return None

    def __count(self):
        try:
            bamfile = pysam.AlignmentFile(self.bamfile, 'rb')
        except OSError as error:
            print(error)
            return None

        print('Couting number of reads of {}'.format(
            self.bamfile))
        read_counters = []
        for row in self.bed.targets.itertuples():
            read_counters.append(bamfile.count(row[1], row[2], row[3]))
        return read_counters

    def __parallel_count(self, cores):
        print('Couting number of reads of {}, using {} cores'.format(
            self.bamfile, cores))
        counters = []

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

        with pysam.AlignmentFile(self.bamfile, 'rb') as bamfile:
            counters = MPPoolHandler(tasks, cores).run()

        if counters:
            counters = [c for counter in counters for c in counter]
        return counters

    def __count_pools(self):
        '''
        Count number of reads per pool.
        This method is used when loading a NRR from a file
        '''
        print('Counting reads by pools')
        targets = list(self.bed.targets.itertuples())
        if (len(targets) != len(self.counters)):
            print('Number of targets and their read counters differ.')
            print('Aborting!')
            return None

        reads_by_pool = defaultdict(int)
        for i, t in enumerate(targets):
            for pool in t[5].pools:
                reads_by_pool[pool] += self.counters[i]/len(t[5].pools)
        self.nreads = sum(self.counters)

        return reads_by_pool

    def __norm(self, mag=1000000):
        print('Normalizing counters')
        normalized = []
        targets = list(self.bed.targets.itertuples())
        if (len(targets) != len(self.counters)):
            print('Number of targets and their read counters differ.')
            print('Aborting!')
            return None

        for i, t in enumerate(targets):
            current_pools_counter = []
            for pool in t[5].pools:
                current_pools_counter.append((self.counters[i] /
                                             self.reads_by_pool[pool]))
            normalized.append(mag * (sum(current_pools_counter) /
                                     len(t[5].pools)))

        return normalized

    def __metric2label(self, counters, Q1, Q3, metric, interval_range):
        labels = []
        for counter in counters:
            if below_range(counter, metric, Q1, interval_range):
                labels.append('-')
            elif above_range(counter, metric, Q3, interval_range):
                labels.append('+')
            else:
                labels.append('o')
        return labels

    def __label_targets(self, iqr_range=1.5, std_range=1.5, mode='normalized'):
        '''
        Label targets considering IQR on counters
        (normalized or not - user's choice) for
        discoverying whether a target is in
        the lower (-) quartile, upper (+) quartile, or middle (o)

        Parameters:
            iqr_range (number): value to multiply IQR by
            std_range (number): value to multiply std by
            mode (str): 'normalized' for employing self.normalized_counters,
            anything else for applying self.counters

        Returns:
            labels (list): a list of size len(self.counters) representing the
            labels for the targets
        '''
        # compute metric (IQR)
        if mode == 'normalized':
            counters = self.normalized_counters
        elif mode == 'log':
            counters = log(self.normalized_counters)
        else:
            counters = self.counters
        df = DataFrame(counters)
        Q1, Q3, metric = IQR(df, 0)
        # label targets considering interquartile range (IQR)
        labels = self.__metric2label(counters, Q1, Q3, metric, iqr_range)
        # filter out counters + or - labeled
        df = DataFrame(labels)
        df.loc[:, len(df.columns)] = counters
        filtered = df[df[0] == 'o']
        # label targets considering std
        mean = filtered.iloc[:, 1].mean()
        Q1, Q3, metric = compute_metric(filtered, 1, 'std', center=mean)
        labels = self.__metric2label(counters, Q1, Q3, metric, std_range)
        return labels

    def count_label_by_pool(self):
        '''
        Count the number of targets arranged by label in each pool

        Returns:
        ldf (DataFrame): number of targets considering pools x labels
        '''

        df = DataFrame(self.bed.targets.iloc[:, 4])
        df[len(df.columns)] = self.labels

        counters = defaultdict(str)
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
    '''
    NRR list management
    '''
    def __init__(self, bedfile=None, bamfiles=None, region=None,
                 bed=None, parallel=True, to_classify=False):
        '''
        Parameters:
             bedfile (str): filename (BED) where the amplicons are listed in
             bamfiles (list): list of aligment file names (bam format)

             region (str): region limit for analysis (view). It must stick to
             the format: chrom:chromStart-chromEnd

             bed (ROI): list of amplicons already loaded into memory
             parallel: whether count target read depth in parallel
        '''
        self.bedfile = bedfile
        self.bamfiles = bamfiles
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

        for i, bamfile in enumerate(self.bamfiles):
            self.list.append(NRR(bedfile=bedfile,
                                 bamfile=bamfile,
                                 bed=bed,
                                 region=region,
                                 parallel=parallel))

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

    def makemean(self):
        '''
        Compute baseline counter mean
        '''
        self.normalized_mean = []
        self.sd = []
        if self.__counters and len(self.__counters) == len(self.list):
            try:
                self.mean = [sum(c)/len(c) for c in zip(*self.__counters)]

                for counters in zip(*self.__normalized_counters):
                    mean = sum(counters)/len(counters)
                    deviations = [(c - mean) ** 2 for c in counters]
                    variance = sum(deviations)/len(counters)
                    sd = sqrt(variance)
                    self.normalized_mean.append(mean)
                    self.sd.append(sd)
            except TypeError as error:
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

    def _addpools(self, targets):
        '''
        add pools info to the df in order to verify whether +, -, or o
        amplicons are pool related
        '''
        pools = [target[5].pools for target in targets.itertuples()]
        return pools

    def __label_targets(self):
        # targets
        targets = self.list[0].bed.targets
        # make a copy of targets df base (chrom:chromStart-chromEnd)
        df = targets.iloc[:, :4].copy()
        # add pools into table
        df.loc[:, len(df.columns)] = self._addpools(targets)
        for sample in self.list:
            df.loc[:, len(df.columns)] = sample.labels
        filtered_df = df.iloc[:, 5:len(df.columns)]
        return df, filtered_df.apply(classify_by_count, axis=1)


@validstr('configfile', empty_allowed=True)
@validstr('path', empty_allowed=False)
class NRRTest(cdf):
    '''
    Hold information about tests between a NRR baseline
    and a NRR test sample
    '''
    def __init__(self, baseline, sample, path='results',
                 size=200, step=10, metric='std', interval_range=3,
                 minread=30, below_cutoff=0.7, above_cutoff=1.3,
                 maxdist=15000000, cnv_like_range=0.7,
                 bins=500, method='chr_group'):
        '''
        Construct and initialize a NRRTest object.
        Parameters:
             baseline (NRRList): obj representing exomes of the baseline
             sample (NRR): obj representing the exome of a test sample
             path (str): where to save (path/filename) analysis results
             size (int): block size when sliding window
             step (int): step size when sliding window

             metric (str): param used to define which metric should be used
             when detecting CNVs:
                  'std' = standard deviation
                  'IQR' = Interquartile Range

             interval_range (number): value to multiply metric by
             minread (number): minimum number of reads used to filter targets
             below_cutoff (number): filter out data (ratios) below this cutoff
             above_cutoff (number): filter out data (ratios) above this cutoff

             maxdist (number): maximum distance allowed of a cnv-like block, to
             its closest cnv block, for it be a cnv as well.

             cnv_like_range (number): value to multiply interval_range by in
             order to detect cnv-like (cnvs when applying looser calculations)

             bins (int): number of bins to use when plotting ratio data
             method (str): method used in order to group rations when plotting
        '''
        self.baseline = baseline
        self.sample = sample
        self.ratios = []
        self.df = None
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
        '''
        create NRRTest dataframe with sample, baseline, and
        test data
        '''
        df = self.sample.bed.targets.copy()
        df.loc[:, len(df.columns)] = self.sample.counters
        df.loc[:, len(df.columns)] = self.baseline.mean
        df.loc[:, len(df.columns)] = self.sample.normalized_counters
        df.loc[:, len(df.columns)] = self.baseline.normalized_mean
        df.loc[:, len(df.columns)] = self.baseline.sd
        df.loc[:, len(df.columns)] = self.ratios
        df.columns = self.columns
        self.df = df

    def makeratio(self):
        '''
        Compute between the standardized read depth values from a test sample
        and the mean baseline read depth
        '''
        print('Calculating reference bam files # reads mean')
        self.baseline.makemean()

        print('Calculating # reads test sample/ baseline')
        ratios = []

        try:
            for i in range(len(self.sample.counters)):
                if (self.baseline.mean[i] > self.minread):
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
        except IndexError as error:
            print('Amount of read counters (baseline and sample) differ')
            print('Are their region targets different? Aborting!')
        except TypeError as error:
            print('Amount of read counters (baseline and sample) differ')
            print('Are their region targets different? Aborting!')

    def plot(self):

        if self.ratios:
            if self.method == 'chr_group':
                xs, ys, offsets, names = self.__regroup(self.bins, self.method)
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
            print('You should call NRRTest.makeratio() ' +
                  'before plotting something!')

    def filter(self, df):
        return df[df.loc[:, 'mean'] >= self.minread]

    def __remean(self, ratios, bins, method):

        if isinstance(bins, str):
            bins = int(bins)

        if len(ratios) < bins:
            bins = len(ratios)

        if method == 'blocks':
            xs = [i*bins for i in range(len(ratios)//bins)]
        else:
            xs = [i for i in range(len(ratios)-bins+1)]
        ys = [sum(ratios[x:x+bins])/bins for x in xs]
        return xs, ys

    def __regroup(self, bins, method):

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

    def summary(self, npools=12):
        '''
        Create a summary of NRRTest samples (baseline and test)
        '''
        # define column names
        columns = ['ntotal']
        poolnames = ['pool_{}'.format(i) for i in range(1, npools + 1)]
        columns.extend(poolnames)
        index = []

        # select data from dataframe
        ntotals = []
        for sample in self.baseline.list:
            nsample = []
            nsample.append(sample.nreads)
            index.append(sample.bamfile.split('/')[-1])
            for i in range(1, npools + 1):
                nsample.append(sample.reads_by_pool[i])
            ntotals.append(nsample)
        df = DataFrame(ntotals, columns=columns, index=index)
        print(df.describe())
        return df

    def merge(self, potential_cnvs):
        '''
        Merge CNV blocks

        Parameters:
             potential_cnvs(list): list of potential CNVs (CNV blocks)

        Returns
             cnvs (list): list of CNVs merged from blocks
        '''
        cnv_dict = defaultdict(str)
        for pcnv in potential_cnvs:
            for cnv in pcnv.itertuples():
                region = '{}:{}-{}'.format(cnv[1], cnv[2], cnv[3])
                targets = self.getin(region)
                for row in targets.itertuples():
                    cnv_dict['{}:{}-{}'.format(row[5].chrom,
                                               row[5].chromStart,
                                               row[5].chromEnd)] = [
                        row[5].chrom,
                        row[5].chromStart,
                        row[5].chromEnd,
                        row[4], row[-1]]
        cnvs = [value for k, value in cnv_dict.items()]
        cnvs.sort(key=lambda x: (x[0], x[1], x[2]))
        return cnvs

    @overrides(cdf)
    def _make_subplots(self, cnv, value_column='ratio',
                       pos_column='chromStart', cnvlike=None):
        # define layout for x and y axis
        layout = defaultdict(lambda: None)
        layout['y'] = dict(title='Ratios', zeroline=False,
                           range=[0, 3])
        layout['x'] = dict(zeroline=False)

        traces, titles = self._create_traces(cnv, value_column,
                                             pos_column, cnvlike=cnvlike,
                                             plotting=True, layout=layout)
        return traces, titles, layout

    @overrides(cdf)
    def _compute(self):
        print('Computing stats on ratio data')

        if not self.ratios:
            self.makeratio()

        all_blocks = []
        columns = ['chrom', 'chromStart', 'chromEnd', 'region', 'median',
                   'isoutlier_1st', 'isoutlier_2nd']

        # filter targets by minread and ratios (below and above cutoff)
        filtered_targets = filter_by_cutoff(self.filter(self.get_autossomic()),
                                            -1, self.below_cutoff,
                                            self.above_cutoff)
        # compute and use metric to find blocks of outliers for each range
        cb, ca, m = compute_metric(filtered_targets, -1, self.metric)
        ranges = [self.interval_range, self.interval_range * self.cnv_like_range]
        for group in self.iterchroms():
            print('Working on blocks of {}'.format(group['id']))

            for block in self.iterblocks(group, size=self.size, step=self.step):
                med_ratio = self.filter(block['df']).loc[:, 'ratio'].median()
                chrom, chromStart, chromEnd = Region(block['id']).as_tuple
                block_data = [chrom, chromStart, chromEnd, block['id'], med_ratio]
                for i, interval_range in enumerate(ranges):
                    block_data.append(above_range(med_ratio, m, ca, interval_range) or
                                      below_range(med_ratio, m, cb, interval_range))
                all_blocks.append(block_data)
        return DataFrame(all_blocks, columns=columns)


@validstr('filename', empty_allowed=False)
class NRRConfig(object):
    '''
    Detect CNVs based on read depth data
    '''
    def __init__(self, filename):
        '''
        Parameters:
             filename (str): the configfile's name
        '''
        self.filename = filename
        self.sections_params = {
            'baseline': 'm',
            'bed': 'm',
            'sample': 'm',
            'output':  'm',
            'targtest': 'o'
        }
        # load configfile
        self.config = ConfigfileParser(self.filename,
                                       self.sections_params)
        # load sample test
        sample = NRR(**self.config.sections['sample'],
                     bedfile=self.config.sections['bed']['bedfile'])
        # load baseline test
        baseline = NRRList(**self.config.sections['baseline'],
                           bedfile=self.config.sections['bed']['bedfile'])
        # make test
        if self.config.sections['targtest']:
            self.nrrtest = NRRTest(baseline, sample,
                                   **self.config.sections['targtest'],
                                   path=self.config.sections['output']['path'])
        else:
            self.nrrtest = NRRTest(baseline, sample,
                                   path=self.config.sections['output']['path'])

        self.nrrtest.makeratio()
        if self.nrrtest.ratios:
            print('Creating plots at {}'.format(self.nrrtest.path2plot))
            self.nrrtest.plot()
            filename = '{}/nrrtest.csv'.format(self.nrrtest.path2table)
            print('Writing table at {}'.format(filename))
            self.nrrtest.df.to_csv(filename, sep='\t',  index=False)
            print('Done!')
