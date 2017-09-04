#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: valengo
"""

from collections import defaultdict
from ..utils import Region
from ..utils import NumberProperty
from pandas import DataFrame
from plotly.offline import plot
from ..graphics import create_subplots
from scipy.signal import medfilt
from abc import ABCMeta
from abc import abstractmethod
from ..graphics import y_scatter


class ChromDF(metaclass=ABCMeta):
    '''
    This class manages a pandas DataFrame that must have its three
    first columns named as:
    "chrom", "chromStart", and "chromEnd"
    '''
    def __init__(self, df):
        '''
        Parameters:
            df (pandas.DataFrame): the dataframe itself
        '''
        self.df = df

    # properties
    size = NumberProperty('_size')

    step = NumberProperty('_step')

    interval_range = NumberProperty('_interval_range')

    minvars = NumberProperty('_minvars')

    alpha = NumberProperty('_alpha')

    diff_cutoff = NumberProperty('_diff_cutoff')

    silhouette_cutoff = NumberProperty('_silhouette_cutoff')

    cnv_like_range = NumberProperty('_cnv_like_range')

    maxdist = NumberProperty('_maxdist')

    size = NumberProperty('_size')

    step = NumberProperty('_step')

    interval_range = NumberProperty('_interval_range')

    minread = NumberProperty('_minread')

    below_cutoff = NumberProperty('_below_cutoff')

    above_cutoff = NumberProperty('_above_cutoff')

    maxdist = NumberProperty('_maxdist')

    cnv_like_range = NumberProperty('_cnv_like_range')

    @property
    def mlabels(self):
        return ['chrom',
                'chromStart',
                'chromEnd']

    @mlabels.setter
    def mlabels(self, value):
        self._mlabels = ['chrom',
                         'chromStart',
                         'chromEnd']

    @property
    def df(self):
        return self._df

    @df.setter
    def df(self, value):
        '''
        Assure that df stick to the rules regarding
        columns len and names.

        Parameters:
            value (pandas.DataFrame): a pandas DataFrame
        '''
        if value is not None:
            columns = value.columns
            if len(columns) < len(self.mlabels):
                raise ValueError('ChromDF must have at least' +
                                 ' {} columns'.format(len(self.mlabels)))
            for i in range(len(self.mlabels)):
                if columns[i] != self.mlabels[i]:
                    raise ValueError('ChromDF column at {}'.format(i) +
                                     ' must be named as "{}"'.format(self.mlabels[i]) +
                                     ', but "{}" was given.'.format(columns[i]))
        self._df = value

    def iterchroms(self, df=None):
        '''
        Create a list to iterate through self.df
        based on its 'chrom' column values

        Parameters:
            df (pandas.DataFrame): a pandas DataFrame
        '''
        if df is None:
            return self.creategroups(self.df)
        else:
            return self.creategroups(df)

    def iterblocks(self, group, size=400, step=10):
        '''
        Given a group (of variants), this method creates
        blocks from it by creating (sliding) windows of a
        given size at a given step

        Parameters:
            group (pandas.DataFrame): a pandas DataFrame
            size (int): block size
            step (int): distance between two consecutive blocks
        '''
        return self.__createblocks(group, size, step)

    def getrows(self, df, region):
        '''
        Get df rows in region

        Parameters:
            df (pandas.DataFrame): pandas DataFrame
            region (str): the region itself. It must stick to the format:
             chrom:chromStart-chromEnd

        Returns:
            rows (pandas.Dataframe)
        '''
        ras_tuple = Region(region).as_tuple
        if ras_tuple:
            if len(ras_tuple) >= 3:
                rows = df[(df.chrom == ras_tuple[0]) &
                          (df.chromStart >= ras_tuple[1]) &
                          (df.chromEnd <= ras_tuple[2])]
            elif len(ras_tuple) >= 2:
                rows = df[(df.chrom == ras_tuple[0]) &
                          (df.chromStart >= ras_tuple[1])]
            else:
                rows = df[df.chrom == ras_tuple[0]]
            return rows
        return df

    def getin(self, region=None, column=None):
        '''
        Get self.df values in region and/or column

        Paramters:
            region (str): the region itself. It must stick to the format:
             chrom:chromStart-chromEnd
            column (str): column name

        Returns:
            rows (pandas.DataFrame): in case column is None
            list(rows) (list): in case column is not None

        '''
        rows = self.getrows(self.df, region)
        if rows is not None and column:
            return list(rows.loc[:, column])
        else:
            return rows

    def medfilter(self, column, kernel_size=9):
        '''
        Apply median filter on self.df column

        Parameters:
            column (str): column name (where to apply median filter)
            kernel_size (int): kernel size for median filter.
        '''
        med_ratios = []
        for group in self.iterchroms():
            med_ratios.extend(medfilt(list(group['df'].loc[:, column]),
                              kernel_size=kernel_size))
        self.ratios = med_ratios
        self.df.loc[:, column] = med_ratios

    def get_autossomic(self):
        '''
        Filter self.df by chromosomes, so rows where df.chrom == chrX or
        df.chrom == chrY are filtered out.

        Returns:
            df (pandas.DataFrame) that contains data of autossomic chromosomes
        '''
        return self.df[(self.df.chrom != 'chrX') &
                       (self.df.chrom != 'chrY')]

    def createid(self, block):
        '''
        Create block id for a given block. The id is based on the block
        region (where it starts and ends)

        Parameters:
            block (pandas.DataFrame): a pandas DataFrame

        Returns:
            id (str): a block id in the format: chrom:chromStart-chromEnd
            None: in case block id creation fails
        '''
        try:
            head = block.head(1)
            tail = block.tail(1)

            chrom = head.iloc[0, 0]
            chrom_start = head.iloc[0, 1]
            chrom_end = tail.iloc[0, 2]

            return '{}:{}-{}'.format(chrom, chrom_start, chrom_end)
        except AttributeError as error:
            print('Failed on getting group/block region!')
        except IndexError as error:
            print('Failed on getting group/block region!')
        return None

    def creategroups(self, df):
        '''
        create a list to iterate through self.df
        based on its 'chrom' column values. In other words, each group
        will contain data of a unique chromosome.

        Parameters:
            df (pandas.DataFrame): a pandas DataFrame

        Returns:
            groups (defaultdict): a dictionary whose keys are regions
            (ids) created by self.createid and items are pandas.DataFrame obj
        '''
        unique_chroms = df.chrom.unique()
        groups = []
        for i, chrom in enumerate(unique_chroms):
            group = defaultdict(lambda: None)
            group['df'] = df[df.chrom == chrom]
            group['id'] = self.createid(group['df'])
            groups.append(group)
        return groups

    def compute(self):
        '''
        compute stats on self.df data. What is computed depends on the class
        that implements self._compute()
        '''
        return self._compute()

    def compute_analyze(self):
        '''
        detect stats on self.df data and analyze them. What is computed
        depends on the class that implements: self._compute()

        '''
        return self._unify_blocks(self._analyze_blocks(self._compute(),
                                                       self.maxdist))

    def compute_plot(self, mode='filtered',
                     filename='cnv-subplot.html'):
        '''
        detect and plot CNVs
        '''
        df = self._compute()
        if mode == 'analyzed':
            cnvs = self._unify_blocks(self._analyze_blocks(df, self.maxdist),
                                      todf=True)
            traces, titles, layout = self._make_subplots(cnvs)
        else:
            cnvs, potential_cnvs = self._filter_blocks(df, todf=True)
            traces, titles, layout = self._make_subplots(self._unify_blocks(cnvs, todf=True),
                                       cnvlike=self._unify_blocks(potential_cnvs,
                                               todf=True))
        fig = create_subplots([traces], titles,
                              layouts=[layout],
                              height=5000)
        plot(fig, filename='{}/{}'.format(self.path2plotcnv,
                                          filename),
             auto_open=False)
        if mode == 'analyzed':
            return cnvs, traces, titles, layout
        else:
            return df, traces, titles, layout

    def _unify_blocks(self, blocks, todf=False):
        '''
        Unify blocks that have overlaps (in regions/id)

        Parameters:
             blocks -- blocks DataFrame
             todf -- whether to return bocks as a pandas.DataFrame or list

        Returns:
             unified_blocks (pandas.DataFrame or list)
        '''
        unified_blocks = []
        for group in self.iterchroms(df=blocks):
            print('Unifying blocks of {}'.format(group['id']))
            rows = list(group['df'].itertuples())
            i = 0
            while i < len(rows):
                chrom = rows[i].chrom
                chromStart = rows[i].chromStart
                chromEnd = rows[i].chromEnd
                j = i + 1
                while j < len(rows) and rows[j].chromStart < chromEnd:
                    chromEnd = rows[j].chromEnd
                    j += 1
                if todf:
                    unified_blocks.append([chrom, chromStart, chromEnd])
                else:
                    unified_blocks.append(['{}:{}-{}'.format(chrom,
                                                             chromStart,
                                                             chromEnd)])
                i = j
        if todf:
            return DataFrame(unified_blocks,
                             columns=['chrom', 'chromStart', 'chromEnd'])
        else:
            return unified_blocks

    def _analyze_blocks(self, blocks, maxdist):
        '''
        Analyze blocks in order to define which CNV-like blocks
        are actual CNVs

        Parameters:
             blocks (pandas.DataFrame): a pandas DataFrame of cnv and
             cnv-like data
             maxdist (number): maximum distance allowed of a cnv-like block, to
             its closest cnv block, for it be a cnv as well.

        Returns:
             cnv_blocks (pandas.DataFrame): cnv + cnv-like that were detected
             as cnvs
        '''
        cnvs, potential_cnvs = self._filter_blocks(blocks, todf=True)
        cnv_blocks = []

        print('Analyzing blocks...')
        for group in self.iterchroms(df=cnvs):
            for row in group['df'].itertuples():
                cnv_blocks.append([row[1], row[2], row[3]])
                pot = potential_cnvs[(potential_cnvs.chrom == row.chrom) &
                                     (potential_cnvs.chromEnd >= row.chromStart - maxdist) &
                                     (potential_cnvs.chromStart <= row.chromEnd + maxdist)]
                for row in pot.itertuples():
                    cnv_blocks.append([row[1], row[2], row[3]])
        cnv_blocks.sort(key=lambda x: (x[0], x[1], x[2]))
        return DataFrame(cnv_blocks, columns=blocks.columns[0:3])

    def _filter_blocks(self, blocks, todf=False):
        '''
        Split blocks dataframe in two dataframes or list: one for CNV blocks
        (1 * interval_range) and another for potential CNV blocks
        (interval_range * self.cnv_like_range)

        Parameters:
             blocks (pandas.DataFrame): a pandas DataFrame of cnv and
             cnv-like data
             todf (boolean): whether return "cnvs" and "potential cnvs"
             as pandas.DataFrame

        Returns:
             cnvs, potential_cnvs (pandas.DataFrame or list)
        '''
        cnvs = self._unify_blocks(blocks[blocks.isoutlier_1st], todf=todf)
        potential_cnvs = self._unify_blocks(blocks[(blocks.isoutlier_1st == False) &
                                                   (blocks.isoutlier_2nd)], todf=todf)
        return cnvs, potential_cnvs

    def __createblocks(self, group, size, step):
        '''
        Given a group (pandas.DataFrame by chromosome), this method
        creates blocks (sliding windows) of it

        Parameters:
             group (pandas.DataFrame): a pandas DataFrame
             size (int): block's size
             step (int): step to take between two consecutives blocks

        Returns:
             blocks (defaultdict):
        '''
        if size > len(group['df']):
            size = len(group['df'])

        blocks = []
        for i in range(0, len(group['df'])-size+step, step):
            block = defaultdict(lambda: None)
            block['df'] = group['df'].iloc[i:i+size, :]
            block['id'] = self.createid(block['df'])
            blocks.append(block)
        return blocks

    def _create_traces(self, cnv, value_column, pos_column, cnvlike=None,
                       plotting=False, layout=None):
        '''
        Create traces with for making subplots

        Parameters:
             cnv (pandas.DataFrame): detected CNVs
             value_column (str): name of the column where to find target
             ratio/baf pos_column (str): name of the column where to find
             target position cnvlike (pandas.DataFrame): detected cnv-like

        Returns:
             traces (list): list of traces
             titles (list): list of titles
        '''

        traces = []
        titles = []
        for group in self.iterchroms():
            # get all values and pos
            chrom, chromStart, chromEnd = Region(group['id']).as_tuple
            titles.append('Chromosome {}'.format(chrom.split('chr')[-1]))
            values = self.getin(column=value_column, region=group['id'])
            values_pos = self.getin(column=pos_column, region=group['id'])

            # get only cnv and/or cnv-like values and pos
            cnvdf_list = [cnv, cnvlike]
            cnv_values, cnv_pos = [], []
            for cnvdf in cnvdf_list:
                current_values, current_pos = [], []
                if cnvdf is not None:
                    chromdf = cnvdf[cnvdf.chrom == chrom]
                    for row in chromdf.itertuples():
                        region = '{}:{}-{}'.format(row[1], row[2], row[3])
                        current_values.extend(self.getin(column=value_column,
                                                         region=region))
                        current_pos.extend(self.getin(column=pos_column,
                                                      region=region))
                cnv_values.append(current_values)
                cnv_pos.append(current_pos)

            # create traces
            traces.append([
                y_scatter(values, x=values_pos, toplot=False, size=3,
                          color='rgb(153, 204, 255)'),
                y_scatter(cnv_values[-1], x=cnv_pos[-1], toplot=False, size=3,
                          color='rgb(255, 102, 0)', name='Potential-CNVs'),
                y_scatter(cnv_values[0], x=cnv_pos[0], toplot=False,
                          size=3, name='CNVs',
                          color='rgb(153, 0, 0)')
                ])
            chrom_trace = [
                y_scatter(values, x=values_pos, toplot=False, size=3,
                          color='rgb(153, 204, 255)'),
                y_scatter(cnv_values[-1], x=cnv_pos[-1], toplot=False, size=3,
                          color='rgb(255, 102, 0)', name='Potential-CNVs'),
                y_scatter(cnv_values[0], x=cnv_pos[0], toplot=False,
                          size=3, name='CNVs',
                          color='rgb(153, 0, 0)')
                ]
            chrom_title = 'Chromosome {}'.format(chrom.split('chr')[-1])
            filename = '{}/{}-{}-analyzed.html'.format(self.path2plotcnv,
                                                       group['id'].split(':')[0],
                                                       self.rootname)
            fig = create_subplots([[chrom_trace]], [chrom_title],
                                  layouts=[layout],
                                  height=400)
            plot(fig, filename=filename, auto_open=False)

        return traces, titles

    @abstractmethod
    def _make_subplots(self):
        pass

    @abstractmethod
    def _compute(self):
        pass

    @abstractmethod
    def _createdf(self):
        pass
