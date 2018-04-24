#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: valengo
"""
from typing import Union

from numpy import histogram, float64
from statistics import mode
from statistics import StatisticsError
from collections import defaultdict

from pandas import DataFrame


def classify_by_count(row: list, default: str = 'u', prop: int = 0, attrs: Union[list, set] = {'-', '+'}):
    """
    Classify a given list (register) according to the
    appearance proportion of unique attributes (column values).
    The main idea is that some attributes defined by 'attrs' parameter
    cannot display simultaneously a proportion higher than 'prop' considering
    the total number of attributes (prop rule).
    When a given register doesn't brake this rule, it is classified according
    to its own mode.

    :param list row: data
    :param str default: name of the class to be returned when register brakes the proportion rule
    :param float prop: max proportion allowed for each attr (simultaneously)
    :param list attrs: list of attr that are used in the prop rule
    :return: register's class (typeof(row[0]))
    """
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
    if above_prop == (len(attrs)):
        return default
    # otherwise, classify by mode
    return classify_by_mode(row, default)


def classify_by_mode(row: list, default: str = 'u'):
    """
    Classify a given list according to the most frequent element.
    In case there is no most frequent, the class is defined
    by default str

    :param list row: data
    :param str default: name of the class when there is no most frequent element
    :return: register's class (typeof(row[0]))
    """
    try:
        row_mode = mode(row)
    except StatisticsError:
        row_mode = default
    return row_mode


def filter_by_cutoff(df: DataFrame, col: int, below_cutoff: float, above_cutoff: float) -> DataFrame:
    """
    Filter a given df depending on values (below_cutoff and above_cutoff) that are in a given col

    :param DataFrame df: dataframe with target data
    :param int col: column where to look for data
    :param float below_cutoff: filter out rows that present value below this cutoff at column col
    :param float above_cutoff: filter out rows that present value above this cutoff at column col
    :return: filtered dataframe
    """
    return df[(df.iloc[:, col] < above_cutoff) &
              (df.iloc[:, col] > below_cutoff)]


def iqr(df: DataFrame, col: int) -> tuple:
    """
    Compute Interquartile Range (IQR)

    :param DataFrame df: dataframe with target data
    :param int col: column where to look for data
    :return: Q1, Q3, and IQR
    """
    q1 = df.iloc[:, col].quantile(0.25)
    q3 = df.iloc[:, col].quantile(0.75)
    return q1, q3, q3 - q1


def std(df: DataFrame, col: int) -> float64:
    """
    Compute Standard Deviation (std)

    :param DataFrame df: dataframe with target data
    :param int col: column where to look for data
    :return: std
    """
    return df.iloc[:, col].std()


def below_range(test_value: float, metric: float, center: float, interval_range: float) -> bool:
    """
    Check whether test_value is smaller (<) than a certain value
    that is computed as: center - metric * interval_range

    :param float test_value: value to be compared
    :param float metric: value of metric used to comparison (IQR or std)
    :param float center: center of comparison (for example, Q1)
    :param float interval_range: value to multiply metric by
    :return: whether test value is below range
    """
    return test_value < (center - metric * interval_range)


def above_range(test_value: float, metric: float, center: float, interval_range: float) -> bool:
    """
    Check whether test_value is larger (>) than a certain value
    that is computed as: center + metric * interval_range

    :param float test_value: value to be compared
    :param float metric: value of metric used to comparison (IQR or std)
    :param float center: center of comparison (for example, Q1)
    :param float interval_range: value to multiply metric by
    :return: whether test value is above range
    """
    return test_value > (center + metric * interval_range)


def isbimodal(bafs: list, bins: Union[list, set] = {0.2, 0.4, 0.6, 0.8}, interval_range: float = 1.5) -> bool:
    """
    Verify if data (bafs) have bimodal distribution

    :param list bafs: values
    :param list bins: bins for histogram
    :param float interval_range: range for comparison among bins
    :return: whether data has evidence of bimodal distribution
    """
    as_list = list(bins)
    as_list.sort()
    hist = histogram(bafs, as_list)
    return hist[0][1] * interval_range < hist[0][0] + hist[0][2]


def compute_metric(df: DataFrame, col: int, metric: str, center: float = 1) -> tuple:
    """
    Compute IQR or STD depending on "metric" param

    :param DataFrame df: dataframe that contains target data
    :param int col: column where to look for data
    :param str metric: IQR or std
    :param number center: center of values
    :return: iqr() or std()
    """

    if metric == 'IQR':
        return iqr(df, col)
    else:
        return center, center, std(df, col)
