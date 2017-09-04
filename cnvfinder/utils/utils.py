#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: valengo
"""

from collections import defaultdict
import configparser
import os
import errno


class GenericDescriptor(object):
    '''
    Get and set attribute data in a given instance
    by using its getter and setter
    '''
    def __init__(self, getter, setter):
        self.getter = getter
        self.setter = setter

    def __get__(self, instance, owner=None):
        if instance is None:
            return self
        return self.getter(instance)

    def __set__(self, instance, value):
        return self.setter(instance, value)


def validstr(attr_name, empty_allowed=True):
    '''
    Check if a given attribute is a valid str

    Parameters:
    attr_name (str): name of the attribute
    empty_allowed (boolean): whether attr_name is allowed to be None/empty

    Returns:
         decorator
    '''
    def decorator(cls):
        '''
        Take a class definition, process it, and return a modified class
        with the same name as the one it was passed
        '''
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


def number(s):
    '''
    try to cast a string (if any) to a number (int or float)

    Parameters:
         s (str): string to be converted

    Returns:
         s (int/float) if it's already a float or int
         n (int/float) s converted to int or float
    '''
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
    '''
    Define number properties
    '''
    def __init__(self, name):
        '''
        Parameters:
             name (str): property name
        '''
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
    '''
    Handle region string (chr:chrStart-chrEnd)
    '''
    def __init__(self, region):
        '''
        Parameters:
             region (str): the region itself
        '''
        self.as_string = region
        self.as_tuple = self.__totuple(region)

    @property
    def as_tuple(self):
        return self._as_tuple

    @as_tuple.setter
    def as_tuple(self, value):
        if (value is None or (isinstance(value, tuple) and len(value) <= 3 and
                              len(value) >= 1)):
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

    def __totuple(self, region):
        '''
        Write string "chr:chrStart-chrEnd" as tuple: (chr, chrStart, chrEnd)

        Returns:
             location (tuple): if regions exists
             None: if there is n region
        '''
        if self.as_string is None:
            return None

        location = []

        split4chrom = region.split(':')
        location.append(split4chrom[0])

        if len(split4chrom) > 1:
            split4bounds = split4chrom[1].split('-')
            try:
                location.append(int(split4bounds[0]))
            except ValueError as error:
                print('Region start {0} '.format(split4bounds[0]) +
                      'is not a number. Skipping!')

            if len(split4bounds) > 1 and len(location) > 1:
                try:
                    location.append(int(split4bounds[1]))
                except ValueError as error:
                    print('Region end {0} '.format(split4bounds[1]) +
                          'is not a number. Skipping!')

        return tuple(location)


class ConfigfileParser(object):
    '''
    Handle config.ini file parsing
    '''
    def __init__(self, filename, params):
        '''
        Parameters:
             filename (str): filename
             params (dict): dictionary discribing whether params
             are mandatory. For example:
                  self.sections_params = {
                      'baseline': 'm',
                      'bed': 'm',
                      'sample': 'm',
                      'output':  'm',
                      'targtest': 'o'
                  }
        '''
        self.filename = filename
        self.sections_params = params
        self.sections = self.__parseconfig()

    def __getopts(self, section):
        '''
        Get options from section

        Parameters:
             section (str): label for a section

        Returns:
             optdict (defaultdict): dictionary containing "section" options
        '''
        config = configparser.ConfigParser()
        config.read(self.filename)

        try:
            options = config.options(section)
        except KeyError as error:
            print('Section [{0}] not provided. Skipping!'.format(section))
            return None
        except configparser.NoSectionError as error:
            print('Section [{0}] not provided. Skipping!'.format(section))
            return None

        optdict = defaultdict(lambda: None)
        for opt in options:
            try:
                optdict[opt] = config.get(section, opt).strip().split('\n')
                if len(optdict[opt]) == 1:
                    optdict[opt] = optdict[opt][0]
            except:
                print('Error getting {0} '.format(opt) +
                      'option from section [{1}]'.format(section))
                optdict[opt] = None
        return optdict

    def __parseconfig(self):
        '''
        Parse config file
        '''
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
    '''
    Manages two different file extensions (rename files)
    '''
    def __init__(self, filename, orig_ext, prom_ext):
        '''
        Parameters:
             filename (str): filename
             orig_ext (str): original extension
             prom_ext (str): new extension
        '''
        self.orig_ext = orig_ext
        self.prom_ext = prom_ext
        self.filename = filename
        self.new_filename = self.filename.replace(orig_ext, prom_ext)

    @property
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


def ismultlist(x):
    '''
    Verify if a list is multidimensional

    Parameters:
         x (list): the list itself

     Returns:
          True/False: whether it's a multimensional list
    '''
    if x:
        return ((isinstance(x, list)) and
                isinstance(x[0], list))
    else:
        return False


def createdirs(dirs):
    '''
    Try to create every dir in dirs meanwhile tell the user
    what is going on

    Parameters:
        dirs (list): list of dirs
    '''
    print('In case they don\'t exist, creating directories:')
    for d in dirs:
        createdir(d)
        print('{}: success'.format(d))


def createdir(path):
    '''
    Try to create a dir if it does not exist already

    Parameters:
         path (str): dir path
    '''
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def overrides(interface_class):
    '''
    Define "overrides" property

    Parameters:
         interface_class (str): overrided class name

    Returns:
         overrider (method/obj): the property
    '''
    def overrider(method):
        assert(method.__name__ in dir(interface_class))
        return method
    return overrider


def appenddir(path, dirname):
    '''
    Append a dirname to a path string.

    Parameters:
        path (str): path itself
        dirname (str): dirname itself

    Returns:
        (str): path/dirname
    '''
    if path.endswith('/'):
        return '{}{}'.format(path, dirname)
    else:
        return '{}/{}'.format(path, dirname)
