#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: valengo
"""

from collections import defaultdict
from typing import Union

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
    """
    **This class manages a DataFrame that must have its three first columns named as:
    "chrom", "chromStart", and "chromEnd"**

    :param DataFrame df: the dataframe itself
    """

    def __init__(self, df: Union[DataFrame, None]):
        self.rootname = None
        self.path2plotcnv = None
        self._df = self.df = df
        self.ratios = []

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

    minread = NumberProperty('_minread')

    below_cutoff = NumberProperty('_below_cutoff')

    above_cutoff = NumberProperty('_above_cutoff')

    @property
    def mlabels(self):
        return ['chrom',
                'chromStart',
                'chromEnd']

    @property
    def df(self):
        return self._df

    @df.setter
    def df(self, value: Union[DataFrame, None]):
        """
        Assure that df stick to the rules regarding columns len and names.

        :param DataFrame value: the dataframe itself
        """
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

    def iterchroms(self, df: DataFrame = None) -> list:
        """
        Create a list to iterate through self.df

        :param DataFrame df: dataframe itself
        :return: dataframe grouped by chromosomes
        """
        if df is None:
            return self.create_groups(self.df)
        else:
            return self.create_groups(df)

    def iterblocks(self, group: DataFrame, size: int = 400, step: int = 10) -> list:
        """
        Given a chromosome group, this method creates blocks from it by creating (sliding) windows of a
        given size at a given step

        :param Dataframe group: dataframe containing data
        :param int size: block's size
        :param int step: distance between two consecutive blocks
        :return: dictionary of blocks {}
        """
        return self.__createblocks(group, size, step)

    @staticmethod
    def getrows(df: DataFrame, region: str) -> DataFrame:
        """
        Get dataframe rows in region

        :param DataFrame df: dataframe itself
        :param str region: it should be in the form: chr1:10000-90000
        :return: a new dataframe containing only rows within region
        """
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

    def getin(self, region: str = None, column: str = None) -> Union[DataFrame, list, None]:
        """
        Get self.df values in region and/or column

        :param str region: it should be in the form: chr1:10000-90000
        :param str column: column's name
        :return: None if no data is found; a dataframe in case column is None; a list if column is not None
        """
        rows = self.getrows(self.df, region)
        if rows is not None and column:
            return list(rows.loc[:, column])
        else:
            return rows

    def medfilter(self, column: str, kernel_size: int = 9):
        """
        Apply median filter on self.df column. Values are stored at this very same column

        :param str column: column to lookup for values that will be median filtered
        :param int kernel_size: median filter kernel's size
        """
        med_ratios = []
        for group in self.iterchroms():
            med_ratios.extend(medfilt(list(group['df'].loc[:, column]),
                                      kernel_size=kernel_size))
        self.ratios = med_ratios
        self.df.loc[:, column] = med_ratios

    def get_autosomal(self):
        """
        Filter self.df by chromosomes, so rows where df.chrom == chrX or
        df.chrom == chrY are filtered out

        :return: dataframe that contains data only for autosomal chromosomes
        """
        return self.df[(self.df.chrom != 'chrX') &
                       (self.df.chrom != 'chrY')]

    @staticmethod
    def create_id(block: DataFrame) -> Union[str, None]:
        """
        Create block id for a given block. The id is based on the block
        region (where it starts and ends)

        :param DataFrame block: data
        :return: a block id in the format: chrom:chromStart-chromEnd or None in case block id creation fails
        """
        try:
            head = block.head(1)
            tail = block.tail(1)

            chrom = head.iloc[0, 0]
            chrom_start = head.iloc[0, 1]
            chrom_end = tail.iloc[0, 2]

            return '{}:{}-{}'.format(chrom, chrom_start, chrom_end)
        except AttributeError:
            print('Failed on getting group/block region!')
        except IndexError:
            print('Failed on getting group/block region!')
        return None

    def create_groups(self, df: DataFrame) -> list:
        """
        Create a list to iterate through self.df based on its 'chrom' column values. In other words, each group
        will contain data of a unique chromosome

        :param DataFrame df: data
        :return: list of groups represented as: {'id': 'chrom:chromStart-chromEnd', 'df': DataFrame}
        """
        unique_chroms = df.chrom.unique()
        groups = []
        for i, chrom in enumerate(unique_chroms):
            group = defaultdict(lambda: None)  # type: Union[DataFrame, str]
            group['df'] = df[df.chrom == chrom]
            group['id'] = self.create_id(group['df'])
            groups.append(group)
        return groups

    def compute(self) -> DataFrame:
        """
        compute stats on self.df data. What is computed depends on the class
        that implements self._compute()
        """
        return self._compute()

    def compute_analyze(self):
        """
        detect stats on self.df data and analyze them. What is computed
        depends on the class that implements: self._compute()
        """
        return self._unify_blocks(self._analyze_blocks(self._compute(),
                                                       self.maxdist))

    def compute_plot(self, mode: str = 'filtered',
                     filename: str = 'cnv-subplot.html') -> tuple:
        """
        Detect and plot CNVs

        :param str mode: 'analyzed' -> CNV and CNV-like blocks are merged; 'filtered' -> non CNV blocks are filtered out
        :param str filename: path to output file
        :return: analyzed or filtered CNVs as dataframe, traces for plotting, plot titles, and layout
        """
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

    def _unify_blocks(self, blocks: DataFrame, todf: bool = False) -> Union[DataFrame, list]:
        """
        Unify blocks that have overlaps in regions/id

        :param DataFrame blocks: data
        :param bool todf: whether to return blocks as pandas.DataFrame or list
        :return: unified blocks as list or DataFrame
        """
        unified_blocks = []
        for group in self.iterchroms(df=blocks):
            print('Unifying blocks of {}'.format(group['id']))
            rows = list(group['df'].itertuples())
            i = 0
            while i < len(rows):
                chrom = rows[i].chrom
                chrom_start = rows[i].chromStart
                chrom_end = rows[i].chromEnd
                call = rows[i].call
                j = i + 1
                while j < len(rows) and rows[j].chromStart < chrom_end and rows[j].call == call:
                    chrom_end = rows[j].chromEnd
                    j += 1
                if todf:
                    unified_blocks.append([chrom, chrom_start, chrom_end, call])
                else:
                    unified_blocks.append(['{}:{}-{}'.format(chrom,
                                                             chrom_start,
                                                             chrom_end)])
                i = j
        if todf:
            return DataFrame(unified_blocks,
                             columns=['chrom', 'chromStart', 'chromEnd', 'call'])
        else:
            return unified_blocks

    def _analyze_blocks(self, blocks: DataFrame, maxdist: int):
        """
        Analyze blocks in order to define which CNV-like blocks are actual CNVs

        :param DataFrame blocks: a pandas DataFrame of cnv and cnv-like data
        :param int maxdist: maximum distance allowed of a CNV-like block to its closest CNV block, it's also a CNV block
        :return:
        """
        cnvs, potential_cnvs = self._filter_blocks(blocks, todf=True)
        cnv_blocks = []
        print('Analyzing blocks...')
        for group in self.iterchroms(df=cnvs):
            for row in group['df'].itertuples():
                cnv_blocks.append([row[1], row[2], row[3], row[-1]])
                pot = potential_cnvs[(potential_cnvs.chrom == row.chrom) &
                                     (potential_cnvs.chromEnd >= row.chromStart - maxdist) &
                                     (potential_cnvs.chromStart <= row.chromEnd + maxdist)]
                for r in pot.itertuples():
                    cnv_blocks.append([r[1], r[2], r[3], r[-1]])
        cnv_blocks.sort(key=lambda x: (x[0], x[1], x[2]))
        return DataFrame(cnv_blocks, columns=['chrom', 'chromStart', 'chromEnd', 'call'])

    def _filter_blocks(self, blocks: DataFrame, todf: bool = False) -> tuple:
        """
        Split blocks dataframe in two dataframes or list: one for CNV blocks (1 * interval_range) and another for
        CNV-like blocks (interval_range * self.cnv_like_range)

        :param DataFrame blocks: data
        :param bool todf: whether to return as DataFrame
        :return:
        """
        cnvs = self._unify_blocks(blocks[blocks.isoutlier_1st], todf=todf)
        potential_cnvs = self._unify_blocks(blocks[(blocks.isoutlier_1st == False) &
                                                   blocks.isoutlier_2nd], todf=todf)
        return cnvs, potential_cnvs

    def __createblocks(self, group: DataFrame, size: int, step: int) -> list:
        """
        Given a group, this method creates blocks using a sliding window approach

        :param DataFrame group: data
        :param int size: block's size
        :param int step: step to take between two consecutive blocks
        :return: a list of blocks
        """
        if size > len(group['df']):
            size = len(group['df'])

        blocks = []
        for i in range(0, len(group['df']) - size + step, step):
            block = defaultdict(lambda: None)  # type: Union[DataFrame, str]
            block['df'] = group['df'].iloc[i:i + size, :]
            block['id'] = self.create_id(block['df'])
            blocks.append(block)
        return blocks

    def _create_traces(self, cnv: DataFrame, value_column: str, pos_column: str, cnvlike: DataFrame = None,
                       toplot: bool = False, layout: defaultdict = None, auto_open: bool = False) -> tuple:
        """
        Create traces for making subplots

        :param DataFrame cnv: detected CNVs
        :param str value_column: column's name to lookup for Y axis values
        :param str pos_column: column's name to lookup for X axis values
        :param DataFrame cnvlike: detected CNV-like
        :param list bool toplot: whether to plot intermediate plots
        :param defaultdict layout: plot's layout
        :param bool auto_open: whether to auto open plots in the browser
        :return traces and titles:
        """
        traces = []
        titles = []
        for group in self.iterchroms():
            # get all values and pos
            chrom, chrom_start, chrom_end = Region(group['id']).as_tuple
            titles.append('Chromosome {}'.format(chrom.split('chr')[-1]))
            values = self.getin(column=value_column, region=group['id'])
            values_pos = self.getin(column=pos_column, region=group['id'])

            # get only cnv and/or cnv-like values and pos
            cnvdf_list = [cnv, cnvlike]
            cnv_values, cnv_pos = [], []
            for cnvdf in cnvdf_list:
                current_values, current_pos = [], []
                if cnvdf is not None:
                    chrom_df = cnvdf[cnvdf.chrom == chrom]
                    for row in chrom_df.itertuples():
                        region = '{}:{}-{}'.format(row[1], row[2], row[3])
                        current_values.extend(self.getin(column=value_column,
                                                         region=region))
                        current_pos.extend(self.getin(column=pos_column,
                                                      region=region))
                cnv_values.append(current_values)
                cnv_pos.append(current_pos)

            # create traces
            traces.append([
                y_scatter(values, x=values_pos, toplot=toplot, size=3,
                          color='rgb(153, 204, 255)'),
                y_scatter(cnv_values[-1], x=cnv_pos[-1], toplot=toplot, size=3,
                          color='rgb(255, 102, 0)', name='Potential-CNVs'),
                y_scatter(cnv_values[0], x=cnv_pos[0], toplot=toplot,
                          size=3, name='CNVs',
                          color='rgb(153, 0, 0)')
            ])
            chrom_trace = [
                y_scatter(values, x=values_pos, toplot=toplot, size=3,
                          color='rgb(153, 204, 255)'),
                y_scatter(cnv_values[-1], x=cnv_pos[-1], toplot=toplot, size=3,
                          color='rgb(255, 102, 0)', name='Potential-CNVs'),
                y_scatter(cnv_values[0], x=cnv_pos[0], toplot=toplot,
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
            plot(fig, filename=filename, auto_open=auto_open)

        return traces, titles

    @abstractmethod
    def _make_subplots(self, cnv: DataFrame, value_column: str = '',
                       pos_column: str = '', cnvlike: DataFrame = None) -> tuple:
        pass

    @abstractmethod
    def _compute(self) -> DataFrame:
        pass

    @abstractmethod
    def _createdf(self):
        pass

    @abstractmethod
    def _call(self, df: DataFrame) -> DataFrame:
        pass
