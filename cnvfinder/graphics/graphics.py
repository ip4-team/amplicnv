#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: valengo
"""

from plotly.offline import plot
import plotly.graph_objs as go
from collections import defaultdict
from plotly import tools
from ..utils import ismultlist


def scatter(x, y, yax='auto',
            offset=None, name=None,
            mode='markers',
            filename='temp-plot.html',
            mult=2):

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


def create_traces(xs, ys, offsets, names, mode, mult, size=3):
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
    '''
    This class creates a heatmap using reads counting
    values from a BAM file
    '''
    def __init__(self, counters, pools, nbins=10):
        '''
        Calls functions in order to create the heatmap

        Parameters:
             counters (list): it is a list of reads counters. Each value
             was generated using pysam.count() - see baseline.py

             pools (list): it is a list of pools that are related to
             the counters:
                  counters[i] is the amount of reads that "match" a region
                  that is acctually from pools[i]

             nbins (int): number of bins
        '''
        self.__create(counters, pools, nbins)

    def __create(self, counters, pools, nbins):
        '''
        Creates the heatmap itself
        '''
        z = [[] for z in range(max(pools))]
        x = ['Pool ' + str(x) for x in range(1, max(pools) + 1)]

        binned_counts = [defaultdict(int) for i in range(max(pools))]
        for i in range(len(counters)):
            binned_counts[pools[i]-1][1+counters[i]//nbins] += 1

        for i in range(len(binned_counts)):
            for area in range(max(binned_counts[i])):
                z[i].append(binned_counts[i].get(area, 0))

        data = [go.Heatmap(x=x, z=z, transpose=True,
                colorscale='Reds', dy=nbins)]
        plot({"data": data})


def y_scatter(y, x=None, filename='temp-plot.html', toplot=True,
              mode='markers', size=1, name='', color=None,
              auto_open=False):

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


def histogram(x, filename='temp-plot.html', toplot=True,
              autobinx=True, color=None, auto_open=False):

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


def create_subplots(mtraces, titles, layouts=None,
                    height=10000, width=1250):
    '''
    Make subplots using plotly.tools

    Parameters:
         mtraces (list of lists): multidimensional array of plotly scatter/hist
         titles (list): one dimensional array of titles
         layouts (list): one dimensional array of (dicts) layouts
         height (number): plot height
         width (number): plot width
    '''
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
