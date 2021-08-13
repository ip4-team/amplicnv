#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: valengo
"""
import contextlib
import sys
from collections import defaultdict
import configparser
import os
import errno
from typing import Union, Sized

import pkg_resources
from pandas import DataFrame


class GenericDescriptor(object):
    """
    Get and set attribute data in a given instance by using its getter and setter
    """
    def __init__(self, getter, setter):
        self.getter = getter
        self.setter = setter

    def __get__(self, instance, owner=None):
        if instance is None:
            return self
        return self.getter(instance)

    def __set__(self, instance, value):
        return self.setter(instance, value)


def validstr(attr_name: str, empty_allowed: str = True):
    """
    Check if a given attribute is a valid str

    :param str attr_name: name of the attribute
    :param bool empty_allowed: whether att_name is allowed to be None/empty
    :return: decorator
    """

    def decorator(cls):
        """
        Take a class definition, process it, and return a modified class
        with the same name as the one it was passed
        """
        name = '_' + attr_name

        def getter(self):
            return getattr(self, name)

        def setter(self, value):
            if value is not None:
                assert isinstance(value, str), (
                        'Attribute ' + attr_name + ' must be str')
            if not empty_allowed and not value:
                raise ValueError('Attribute {0} cannot be empty'.format(
                    attr_name))
            setattr(self, name, value)

        setattr(cls, attr_name, GenericDescriptor(getter, setter))
        return cls

    return decorator


def number(s: str) -> Union[int, float]:
    """
    Try to cast a string (if any) to a number (int or float)

    :param str s: string to be converted
    :return: converted string
    :raise ValueError: if casting fails
    """
    if isinstance(s, str):
        try:
            return int(s)
        except ValueError:
            return float(s)
    elif isinstance(s, float) or isinstance(s, int):
        return s
    else:
        raise ValueError('Failed casting {} to a number'.format(s))


class NumberProperty(object):
    """
    Define number properties

    :param str name: property name
    """

    def __init__(self, name):
        self.name = name

    def __get__(self, obj, objtype):
        return getattr(obj, self.name)

    def __set__(self, obj, value):
        try:
            setattr(obj, self.name, number(value))
        except ValueError:
            name = self.name.split('_')[-1]
            print('Warning: "{}" must be a number or str(number), '.format(name) +
                  'but "{}" was passed!'.format(value))
            print('Things will probably fail in a while...')


@validstr('region', empty_allowed=True)
class Region(object):
    """
    Handle region string (chr:chrStart-chrEnd)

    :param str region: should be in the form: chr:start-end, for example: chr1:10000-90000" or "chr1:10000" or "chr1
    """

    def __init__(self, region):
        self.as_string = region
        self._as_tuple = self.as_tuple = self.__to_tuple(region)

    @property
    def as_tuple(self):
        return self._as_tuple

    @as_tuple.setter
    def as_tuple(self, value):
        if value is None or (isinstance(value, tuple) and 3 >= len(value) >= 1):
            self._as_tuple = value
        else:
            print('Invalid value {0} for Region.as_tuple'.format(value))

    def __repr__(self):
        if self.as_string is None:
            return None
        elif len(self.as_tuple) == 3:
            return '{0}:{1}-{2}'.format(*self.as_tuple)
        elif len(self.as_tuple) == 2:
            return '{0}:{1}'.format(*self.as_tuple)
        elif len(self.as_tuple) == 1:
            return '{0}'.format(*self.as_tuple)

    def __to_tuple(self, region) -> Union[tuple, None]:
        """
        Write string "chr:chrStart-chrEnd" as tuple: (chr, chrStart, chrEnd)

        :param str region: the region itself
        :return: location (tuple): if region exists, otherwise: None
        """
        if self.as_string is None:
            return None

        location = []

        split4chrom = region.split(':')
        location.append(split4chrom[0])

        if len(split4chrom) > 1:
            split4bounds = split4chrom[1].split('-')
            try:
                location.append(int(split4bounds[0]))
            except ValueError:
                print('Region start {0} '.format(split4bounds[0]) +
                      'is not a number. Skipping!')

            if len(split4bounds) > 1 and len(location) > 1:
                try:
                    location.append(int(split4bounds[1]))
                except ValueError:
                    print('Region end {0} '.format(split4bounds[1]) +
                          'is not a number. Skipping!')
        return tuple(location)


class ConfigfileParser(object):
    """
    Handle config.cfg or config.ini file parsing. An example for the 'params' parameter is:

    params = {
                'baseline': 'm',
                'bed': 'm',
                'sample': 'm',
                'output':  'm',
                'targtest': 'o'
            }

    :param str filename: path to configuration file
    :param params: dictionary describing whether params are mandatory
    """

    def __init__(self, filename, params):
        self.filename = filename
        self.sections_params = params
        self.sections = self.__parseconfig()

    def __getopts(self, section: str) -> Union[defaultdict, None]:
        """
        Get options from section

        :param str section: label for section
        :return: if 'section' exists, return a dictionary containing 'section' options. Otherwise return None
        """
        config = configparser.ConfigParser()
        config.read(self.filename)

        try:
            options = config.options(section)
        except KeyError:
            print('Section [{0}] not provided. Skipping!'.format(section))
            return None
        except configparser.NoSectionError:
            print('Section [{0}] not provided. Skipping!'.format(section))
            return None

        optdict = defaultdict(lambda: None)  # type: Union[Sized, list, defaultdict]
        for opt in options:
            try:
                optdict[opt] = config.get(section, opt).strip().split('\n')
                if len(optdict[opt]) == 1:
                    optdict[opt] = optdict[opt][0]
            except KeyError:
                print('Error getting {} '.format(opt) +
                      'option from section [{}]'.format(section))
                optdict[opt] = None
        return optdict

    def __parseconfig(self):
        """
        Parse config file
        """
        sections = defaultdict(lambda: None)
        for key, value in self.sections_params.items():
            sections[key] = self.__getopts(key)
            if sections[key] is None and self.sections_params[key] is 'm':
                raise ValueError('You must provide a [{0}] '.format(key) +
                                 'section on {0} file'.format(self.filename))
        return sections


@validstr('filename', empty_allowed=False)
@validstr('orig_ext', empty_allowed=False)
@validstr('prom_ext', empty_allowed=True)
class ExtensionManager(object):
    """
    Manages two different file extensions (rename files)

    :param str filename: path to file
    :param str orig_ext: original extension
    :param str prom_ext: new extension
    """
    def __init__(self, filename: str, orig_ext: str, prom_ext: str):
        self.orig_ext = orig_ext
        self.prom_ext = prom_ext
        self._filename = self.filename = filename
        self.new_filename = self.filename.replace(orig_ext, prom_ext)

    @property
    def filename(self):
        return self._filename

    @filename.setter
    def filename(self, value):
        if value.endswith(self.orig_ext):
            self._filename = value
        else:
            print('{0} does not end with {1}'.format(value,
                                                     self.orig_ext))

    def __repr__(self):
        if self.filename is None:
            return None
        return '{0} -> {1}'.format(self.filename,
                                   self.new_filename)


def ismultlist(x: list):
    """
    Verify if a list is multidimensional

    :param list x: list itself
    :return: whether list is multidimensional
    """
    if x:
        return ((isinstance(x, list)) and
                isinstance(x[0], list))
    else:
        return False


def createdirs(dirs: list):
    """
    Try to create every dir in dirs meanwhile tell the user
    what is going on

    :param list dirs: list of dirs
    """
    print('In case they don\'t exist, creating directories:')
    for d in dirs:
        createdir(d)
        print('{}: success'.format(d))


def createdir(path: str):
    """
    Try to create a dir if it does not exist already

    :param path to directory
    :raise OsError: when an error other than errno.EEXIST happens
    """
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def overrides(interface_class: str):
    """
    Define "overrides" property

    :param str interface_class: name of the super class
    :return: the property itself
    """

    def overrider(method):
        assert (method.__name__ in dir(interface_class))
        return method

    return overrider


def appenddir(path: str, dirname: str) -> str:
    """
    Append a dirname to a path string

    :param str path: path where to append
    :param str dirname: path to be appended
    :return: path/dirname
    """
    if path.endswith('/'):
        return '{}{}'.format(path, dirname)
    else:
        return '{}/{}'.format(path, dirname)


def bedwrite(filename: str, df: DataFrame, index: bool = False, header: bool = True, sep: str = '\t'):
    """
    Write df as bedfile

    :param str filename: where to write
    :param DataFrame df: what to write
    :param bool index: whether to write row index
    :param bool header: whether to write column names
    :param str sep: column separator
    """
    print('Writing targets to file: "{0}"'.format(filename))
    with contextlib.suppress(FileNotFoundError):
        os.remove(filename)
        os.remove(filename + '.idx')
    df.to_csv(filename, sep=sep, index=index, header=header)


def get_package_name():
    """
    Get package name

    :return: name
    """
    return __name__.split('.')[0]


def resource_exists(filename: str):
    """
    :param filename:
    :return: whether a file exists in package
    """
    return pkg_resources.resource_exists(get_package_name(), filename)


def resource_path_or_exit(filename: str):
    """
    :param filename:
    :return: filename path in package
    """
    if not resource_exists:
        sys.exit('\"{}\" not found! Was it added in setup.py?'.format(filename))
    return pkg_resources.resource_filename(get_package_name(), filename)


def sortable_chromosome(chromosome: str) -> str:
    """
    Format a chromosome, so its sortable

    :param chromosome: chromosome itself
    :return: formatted chromosome
    """
    try:
        return '{:0>2}'.format(int(chromosome))
    except ValueError:
        return chromosome


def sort_chroms(chrom_list: list) -> list:
    """
    Sort chromosomes in a list

    :param chrom_list: list itself
    :return: sorted list
    """
    to_sort = [[sortable_chromosome(chrom.strip('chr')), chrom] for chrom in chrom_list]
    to_sort.sort(key=lambda x: x[0])
    return [chrom[-1] for chrom in to_sort]
