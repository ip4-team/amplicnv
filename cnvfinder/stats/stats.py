#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: valengo
"""
from numpy import histogram
from statistics import mode
from statistics import StatisticsError
from collections import defaultdict


def classify_by_count(row, default='u', prop=0, attrs=['-', '+']):
    '''
    Classify a given list (register) according to the
    appearance proportion of unique attributes (column values).
    The main ideia is that some attributes defined by 'attrs' parameter
    cannot display simultaneously a proportion higher than 'prop' considering
    the total number of attributes (prop rule).
    Wheen a given register doesn't brake this rule, it is classified according
    to its own mode.

    Parameters:
        row (list): data
        default (str): name of the class to be returned when register (list)
        brakes the proportion rule
        prop (float): max proportion allowed for each attr (simultaneously)
        attrs (list): list of attr that are used in the prop rule

    Returns:
        (typeof(row[0])): register's class
    '''
    counters = defaultdict(lambda: 0)
    # count unique label appearances
    for c in row:
        counters[c] += 1
    # count how many labels have counters above len(row) * prop
    above_prop = 0
    for attr in attrs:
        if counters[attr] > len(row) * prop:
            above_prop += 1
    # uncertain register
    if (above_prop) == (len(attrs)):
        return default
    # otherwise, classify by mode
    return classify_by_mode(row, default)


def classify_by_mode(row, default='u'):
    '''
    Classify a given list according to the most frequent element.
    In case there is no most frequent, the class is defined
    by default str.

    Parameters:
        row (list): data
        default (str): name of the class when there is no most frequent element

    Returns:
        (typeof(row[0])): register's class
    '''
    try:
        row_mode = mode(row)
    except StatisticsError:
        row_mode = default
    return row_mode


def filter_by_cutoff(df, col, below_cutoff, above_cutoff):
    '''
    Filter a given df depending on values
    (below_cutoff and above_cutoff) that are in a given col

    Parameters:
         df (pandas.Dataframe): dataframe with target data (to be filtered)
         col (int): column where to look for data
         below_cutoff (number): filter out data below this cutoff
         above_cutoff (number): filter out data above this cutoff

    Returns:
         filtered dataframe (pandas.DataFrame)
    '''
    return df[(df.iloc[:, col] < above_cutoff) &
              (df.iloc[:, col] > below_cutoff)]


def IQR(df, col):
    '''
    Compute Interquartile Range (IQR)

    Parameters:
         df (pandas.DataFrame): dataframe with target data
         col (number): column where to look for data

    Returns:
         Q1, Q3, and IQR
    '''
    Q1 = df.iloc[:, col].quantile(0.25)
    Q3 = df.iloc[:, col].quantile(0.75)
    return Q1, Q3, Q3-Q1


def std(df, col):
    '''
    Compute Standard Deviation (std)

    Parameters:
         df (pandas.DataFrame): dataframe with target data
         col (number): column where to look for data

     Returns:
         std
    '''
    return df.iloc[:, col].std()


def below_range(test_value, metric, center, interval_range):
    '''
    Check whether test_value is smaller (<) than a certain value
    that is computed as: center - metric * interval_range

    Parameters:
         test_value (number): value to be compared
         metric (number): value of metric used to comparison (IQR or std)
         center (number): center of comparison (for example, Q1)
         interval_range (number): value to multiply metric by
    Returns:
         True or False
    '''
    return(test_value < (center - metric * interval_range))


def above_range(test_value, metric, center, interval_range):
    '''
    Check whether test_value is larger (>) than a certain value
    that is computed as: center + metric * interval_range

    Parameters:
         test_value (number): value to be compared
         metric (number): value of metric used to comparison (IQR or std)
         center (number): center of comparison (for example, Q1)
         interval_range (number): value to multiply metric by

    Returns:
         True or False
    '''
    return(test_value > (center + metric * interval_range))


def isbimodal(bafs, bins=[0.2, 0.4, 0.6, 0.8], interval_range=1.5):
    '''
    Verify if data (bafs) have bimodal distribution

    Parameters:
         bafs (list): list of values
         bins (list): list of bins for histogram
         interval_range (number): range for comparison among bins

     Returns:
         True or False
    '''
    hist = histogram(bafs, bins)
    # len(bin[0.4:0.6]) * 1.5 < len(bin[0.2:0.4] + len(bin(0.6:0.8)))
    return(hist[0][1] * interval_range < hist[0][0] + hist[0][2])


def compute_metric(df, col, metric, center=1):
    '''
    Compute IQR or STD depending on "metric" param

    Parameters:
         df (pandas.DataFrame): dataframe that contains target data
         col (int): column where to look for data
         metric (str): 'IQR' or 'std'
         center (number): center of values

    Returns:
         IQR() or std()
    '''

    if metric == 'IQR':
        return IQR(df, col)
    else:
        return center, center, std(df, col)
