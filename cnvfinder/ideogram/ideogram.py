#!/usr/bin/env python3
import sys

import pandas
import pkg_resources
from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection

from cnvfinder.utils.utils import get_package_name


def chromosome_collections(df, y_positions, height, **kwargs):
    """
    Yields BrokenBarHCollection of features that can be added to an Axes
    object.
    Parameters
    ----------
    df : pandas.DataFrame
        Must at least have columns ['chrom', 'start', 'end', 'color']. If no
        column 'width', it will be calculated from start/end.
    y_positions : dict
        Keys are chromosomes, values are y-value at which to anchor the
        BrokenBarHCollection
    height : float
        Height of each BrokenBarHCollection
    Additional kwargs are passed to BrokenBarHCollection
    """
    del_width = False
    if 'width' not in df.columns:
        del_width = True
        df['width'] = df['end'] - df['start']
    for chrom, group in df.groupby('chrom'):
        print(chrom)
        yrange = (y_positions[chrom], height)
        xranges = group[['start', 'width']].values
        yield BrokenBarHCollection(
            xranges, yrange, facecolors=group['colors'], **kwargs)
    if del_width:
        del df['width']


class Ideogram(object):
    def __init__(self, file: str = None, chroms: list = None, chrom_height: int = 1,
                 chrom_spacing: int = 1, fig_size: tuple = (6, 8),
                 colors: dict = None):
        self.chrom_height = chrom_height
        self.chrom_spacing = chrom_spacing
        self.fig_size = fig_size
        self._colors = self.colors = colors
        self._chroms = self.chroms = chroms
        self.file = file
        self.df = self.load_bands()

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

    def m_plot(self):
        ybase = 0
        chrom_ybase = {}
        chrom_centers = {}

        for chrom in self.chroms[::-1]:
            chrom_ybase[chrom] = ybase
            chrom_centers[chrom] = ybase + self.chrom_height / 2.
            ybase += self.chrom_height + self.chrom_spacing

        fig = plt.figure(figsize=self.fig_size)
        ax = fig.add_subplot(111)

        print("adding ideograms...")
        for collection in chromosome_collections(self.df, chrom_ybase, self.chrom_height, edgecolors=(0, 0, 0)):
            ax.add_collection(collection)

        # Axes tweaking
        ax.set_yticks([chrom_centers[i] for i in self.chroms])
        ax.set_yticklabels(self.chroms)
        ax.axis('tight')

        plt.show()

    @staticmethod
    def get_filename_from_pkg():
        file = 'data/cytoBandIdeo.txt'
        pkg_name = get_package_name()

        if not pkg_resources.resource_exists(pkg_name, file):
            sys.exit('{} not found! Was it added in setup.py?')

        return pkg_resources.resource_filename(pkg_name, file)

    def load_bands(self):
        file = self.file if self.file is not None else self.get_filename_from_pkg()
        df = pandas.read_table(file, skiprows=1, names=['chrom', 'start', 'end', 'name', 'gieStain'])
        df = df[df.chrom.apply(lambda x: x in self.chroms)]
        df['width'] = df.end - df.start
        df['colors'] = df.gieStain.apply(lambda x: self.colors[x])
        return df
