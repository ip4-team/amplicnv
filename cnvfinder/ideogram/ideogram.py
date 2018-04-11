#!/usr/bin/env python3

import pandas
from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
from matplotlib.patches import Patch

from cnvfinder.utils import resource_path_or_exit

"""
 attribution: this file (ideogram.py) is based on https://gist.github.com/daler/c98fc410282d7570efc3
"""


def chromosome_collections(df, y_positions, height, to_log=False, **kwargs):
    """
    yields BrokenBarHCollection of features that can be added to an Axes object
    :param to_log: whether to log info
    :param df:
    :param y_positions:
    :param height:
    :param kwargs: are passed to BrokenBarHCollection
    :return:
    """
    del_width = False
    if 'width' not in df.columns:
        del_width = True
        df['width'] = df['chromEnd'] - df['chromStart']
    for chrom, group in df.groupby('chrom'):
        if to_log:
            print(chrom)
        yrange = (y_positions[chrom], height)
        xranges = group[['chromStart', 'width']].values
        yield BrokenBarHCollection(
            xranges, yrange, facecolors=group['colors'], **kwargs)
    if del_width:
        del df['width']


class Ideogram(object):
    def __init__(self, file: str = None, chroms: list = None, chrom_height: float = 1,
                 chrom_spacing: float = 1.5, fig_size: tuple = (6, 8), colors: dict = None,
                 to_log=False):
        """
        Keyword arguments
        :param file: file to load chromosome bands data. BED format with 'chrom', 'chromStart', 'chromEnd', 'name',
        and 'gieStain' as columns. First line is skipped.
        Default: 'cytoBand' table from https://genome.ucsc.edu/cgi-bin/hgTables.
        :param chroms: plot only chromosomes that are in this list.
        Default: ['chr%s' % i for i in list(range(1, 23)) + ['M', 'X', 'Y']]
        :param chrom_height: height of each ideogram
        :param chrom_spacing: spacing between consecutive ideogram
        :param fig_size: width and height in inches
        :param colors: colors for different chromosome stains
        :param to_log: whether to print log info
        """

        self.to_log = to_log
        self.chrom_height = chrom_height
        self.chrom_spacing = chrom_spacing
        self.fig_size = fig_size
        self._colors = self.colors = colors
        self._chroms = self.chroms = chroms
        self.file = file
        self.__pgk_bands_file = 'data/cytoBand.txt'
        self.df = self.load_bands()

        # chromosomes
        self.ybase = 0
        self.chrom_ybase = {}
        self.chrom_centers = {}
        self.fig, self.ax = self.add_chromosomes()

    @property
    def colors(self):
        return self._colors

    @colors.setter
    def colors(self, value):
        default_colors = {
            'gneg': (1., 1., 1.),
            'gpos25': (.6, .6, .6),
            'gpos50': (.4, .4, .4),
            'gpos75': (.2, .2, .2),
            'gpos100': (0., 0., 0.),
            'acen': (.8, .4, .4),
            'gvar': (.8, .8, .8),
            'stalk': (.9, .9, .9)}
        if value is None:
            self._colors = default_colors
        else:
            self._colors = value

    @property
    def chroms(self):
        return self._chroms

    @chroms.setter
    def chroms(self, value):
        if value is None:
            self._chroms = ['chr%s' % i for i in list(range(1, 23)) + ['M', 'X', 'Y']]
        else:
            self._chroms = value

    def add_chromosomes(self):

        for chrom in self.chroms[::-1]:
            self.chrom_ybase[chrom] = self.ybase
            self.chrom_centers[chrom] = self.ybase + self.chrom_height / 2.
            self.ybase += self.chrom_height + self.chrom_spacing

        fig = plt.figure(figsize=self.fig_size)
        ax = fig.add_subplot(111)

        if self.to_log:
            print("adding ideograms...")
        for collection in chromosome_collections(self.df, self.chrom_ybase, self.chrom_height, edgecolors=(0, 0, 0)):
            ax.add_collection(collection)

        # axes tweaking
        ax.set_yticks([self.chrom_centers[i] for i in self.chroms])
        ax.set_yticklabels(self.chroms)
        ax.xaxis.set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)

        return fig, ax

    def save(self, filename: str, **kwargs):
        """
        Save ideograms in a file
        :param filename: filename
        :param kwargs: are passed to pyplot.savefig
        """
        self.ax.axis('tight')
        self.fig.savefig(filename, **kwargs)

    def show(self):
        self.ax.axis('tight')
        self.fig.show()

    def load_bands(self):
        file = self.file if self.file is not None else resource_path_or_exit(self.__pgk_bands_file)
        df = pandas.read_table(file, skiprows=1, names=['chrom', 'chromStart', 'chromEnd', 'name', 'gieStain'])
        df = self.filter_by_chroms(df)
        df['width'] = df.chromEnd - df.chromStart
        df['colors'] = df.gieStain.apply(lambda x: self.colors[x])
        return df

    def filter_by_chroms(self, df: pandas.DataFrame) -> pandas.DataFrame:
        return df[df.chrom.apply(lambda x: x in self.chroms)].copy()

    def add_data(self, df: pandas.DataFrame, height: float = 0.5, padding: float = 0.1,
                 color: str = '#2243a8', alpha: float = 0.5, linewidths: float = 0, **kwargs):
        """
        Add (genomic) data in the plot
        :param df: data
        :param height: height of genomic track. Should be smaller than 'chrom_spacing'
        :param padding: padding between the top of a genomic track and its corresponding ideogram
        :param color: track's color. It will be used in case 'colors' not in df.columns
        :param alpha: alpha value used for blending
        :param linewidths: line widths
        :param kwargs: are passed to BrokenBarHCollection
        """
        df = self.filter_by_chroms(df)

        if 'colors' not in df.columns:
            df['colors'] = color

        data_ybase = {}

        for chrom in self.chroms:
            data_ybase[chrom] = self.chrom_ybase[chrom] + (height + padding)

        for collection in chromosome_collections(df, data_ybase, abs(height),
                                                 alpha=alpha, linewidths=linewidths, **kwargs):
            self.ax.add_collection(collection)

    def add_data_above(self, df: pandas.DataFrame, color: str = None):
        """
        Wrapper for adding data above ideograms
        :param color: bars color
        :param df: data
        """
        self.add_data(df, height=0.5, padding=0.6, color=color)

    def add_data_below(self, df: pandas.DataFrame, color: str = None):
        """
        Wrapper for adding data below ideograms
        :param color: bars color
        :param df: data
        """
        self.add_data(df, height=-0.5, padding=-0.1, color=color)

    def add_legend(self, to_patches: list, loc='lower right', **kwargs):
        """
        Create a legend base on to_patches list
        :param to_patches: list of dict -> {color: color, label: label}
        :param loc: legend location
        :param kwargs: are passed to pyplot.legend
        """
        patches = [Patch(color=p['color'], label=p['label']) for p in to_patches]
        self.fig.legend(handles=patches, loc=loc, **kwargs)
