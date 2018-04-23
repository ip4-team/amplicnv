#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: valengo
"""
from typing import Union

from plotly.offline import plot
import plotly.graph_objs as go
from collections import defaultdict
from plotly import tools
from ..utils import ismultlist


def scatter(x: list, y: list, yax: Union[str, list] = 'auto', offset: list = None, name: list = None,
            mode: str = 'markers', filename: str = 'temp-plot.html', mult: int = 2):
    """
    Plot custom scatter traces

    :param list x: x axis data
    :param list y: y axis data
    :param str yax: y axis range
    :param list offset: offsets for x axis
    :param list name: scatter's names
    :param str mode: markers, lines, or lines+markers
    :param str filename: path to output file
    :param int mult: value to multiply each y axis value
    """
    if (offset is not None and
            name is not None):
        data = create_traces(x, y, offset, name, mode, mult)

    else:
        trace = go.Scatter(
            y=[dy * mult for dy in y],
            x=x,
            mode=mode,
            marker=dict(size=3)
        )
        data = [trace]

    if yax == 'auto':
        layout = go.Layout(yaxis=dict(autorange=True),
                           font=dict(size=10))
    else:
        layout = go.Layout(yaxis=dict(range=yax),
                           font=dict(size=10))
    plot({'data': data, 'layout': layout}, filename=filename, auto_open=False)


def create_traces(xs: list, ys: list, offsets: list, names: list, mode: str, mult: int, size: float = 3) -> list:
    """
    Create custom scatter traces

    :param list xs: x axis data
    :param list ys: y axis data
    :param list offsets: offsets for x axis
    :param list names: trace's names
    :param str mode: markers, lines, or lines+markers
    :param int mult: value to multiply each y axis value
    :param float size: marker's size
    :return: list of traces
    """
    traces = []
    for i in range(len(names)):
        trace = go.Scatter(
            y=[dy * mult for dy in ys[i]],
            x=[dx + offsets[i] for dx in xs[i]],
            mode=mode,
            marker=dict(size=size),
            name=names[i]
        )
        traces.append(trace)
    return traces


class HeatMap(object):
    """
    This class creates a heatmap using reads counting values from a BAM file

    :param list counters: it's a list of read counters. Each value was generated using pysam.count() see nrrhandler.NRR
    :param list pools: it's a list of pools, where counters[i] match a region that is actually from pools[i]
    :param int nbins: number of bins
    """

    def __init__(self, counters, pools, nbins=10):
        self.counters = counters
        self.pools = pools
        self.nbins = nbins
        self.__create()

    def __create(self):
        """
        Creates the heatmap itself
        """
        z = [[] for z in range(max(self.pools))]
        x = ['Pool ' + str(x) for x in range(1, max(self.pools) + 1)]

        binned_counts = [defaultdict(int) for i in range(max(self.pools))]
        for i in range(len(self.counters)):
            binned_counts[self.pools[i] - 1][1 + self.counters[i] // self.nbins] += 1

        for i in range(len(binned_counts)):
            for area in range(max(binned_counts[i])):
                z[i].append(binned_counts[i].get(area, 0))

        data = [go.Heatmap(x=x, z=z, transpose=True,
                           colorscale='Reds', dy=self.nbins)]
        plot({"data": data})


def y_scatter(y: list, x: list = None, filename: str = 'temp-plot.html', toplot: bool = True,
              mode: str = 'markers', size: float = 1, name: str = '', color: str = None,
              auto_open: bool = False) -> go.Scatter:
    """
    Create a scatter plot with or without x axis data

    :param list y: y axis data
    :param list x: x axis data
    :param str filename: path to output file
    :param bool toplot: whether to write the plot
    :param str mode: markers, lines, or lines+markers
    :param float size: marker's size
    :param str name: scatter's name
    :param str color: scatter's color
    :param bool auto_open: whether to auto open the plot
    :return: trace
    """
    if x is None:
        data = [
            go.Scatter(
                y=y,
                mode=mode,
                marker=dict(size=size,
                            color=color),
                name=name
            )
        ]
    else:
        data = [
            go.Scatter(
                x=x,
                y=y,
                mode=mode,
                marker=dict(size=size,
                            color=color),
                name=name
            )
        ]

    if toplot is True:
        plot(data, filename=filename, auto_open=auto_open)
    return data[0]


def histogram(x: list, filename: str = 'temp-plot.html', toplot: bool = True,
              autobinx: bool = True, color: str = None, auto_open: bool = False) -> go.Histogram:
    """
    Create histogram

    :param list x: x axis data
    :param str filename: path to output file
    :param bool toplot: whether to write plot to file
    :param bool autobinx: whether to auto bin histogram
    :param str color: bar's color
    :param bool auto_open: whether to auto open plot in the browser
    :return: histogram
    """
    data = [
        go.Histogram(
            x=x,
            autobinx=autobinx,
            xbins=dict(
                start=0.2,
                end=0.8,
                size=0.2),
            marker=dict(color=color),
        )
    ]
    if toplot is True:
        plot(data, filename=filename, auto_open=auto_open)
    return data[0]


def create_subplots(mtraces: list, titles: list, layouts: list = None,
                    height: int = 10000, width: int = 1250) -> go.Figure:
    """
    Make subplots using plotly.tools

    :param list mtraces: multidimensional list of plotly scatter/hist
    :param list titles: one dimensional list of titles
    :param list layouts: one dimensional list of (dicts) layouts
    :param int height: plot's height
    :param int width: plot's width
    :return: Figure
    """
    if not ismultlist(mtraces):
        raise ValueError('mtraces must be multidimensional!')

    rows = len(mtraces) * len(mtraces[0])
    fig = tools.make_subplots(rows=rows, cols=1,
                              subplot_titles=titles,
                              print_grid=False)
    for i in range(len(mtraces)):
        for j in range(1 + i, len(mtraces[i]) + 1 + i):
            m = len(mtraces)
            s = len(mtraces) - 1
            for trace in mtraces[i][j - (1 + i)]:
                fig.append_trace(trace, j * m - s - i, 1)
            if layouts is not None:
                fig['layout']['xaxis{}'.format(j * m - s - i)].update(layouts[i]['x'])
                fig['layout']['yaxis{}'.format(j * m - s - i)].update(layouts[i]['y'])
    fig['layout'].update(height=height, width=width, showlegend=False)
    return fig
