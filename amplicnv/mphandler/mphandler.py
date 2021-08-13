#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: valengo
"""

from multiprocessing import Pool
from multiprocessing import freeze_support
from multiprocessing import cpu_count


class MPPoolHandler(object):
    """
    Handle tasks multiprocessing using multiprocessing.Pool

    :param list tasks: a list of tasks. Each task is defined as a tuple
    :param int cores: number of cores to be used when multiprocessing
    :param int chunksize: chunksize of tasks to be send each time for a given process

    """
    def __init__(self, tasks: list, cores: int = None, chunksize: int = None):
        self.tasks = tasks
        self._cores = self.cores = cores
        self._chunksize = self.chunksize = chunksize
        self.pool = None

    @property
    def cores(self):
        return self._cores

    @cores.setter
    def cores(self, value):
        if value is None:
            value = cpu_count()
        self._cores = value

    @property
    def chunksize(self):
        return self._chunksize

    @chunksize.setter
    def chunksize(self, value):
        if value is None:
            value, extra = divmod(len(self.tasks), self.cores * 4)
            if extra:
                value += 1
        if value >= len(self.tasks):
            value = 1

        self._chunksize = value

    def __getstate__(self):
        self_dict = self.__dict__.copy()
        del self_dict['pool']
        return self_dict

    def __setstate__(self, state):
        self.__dict__.update(state)

    @staticmethod
    def _runstar(args):
        return args[0](*args[1:len(args)])

    def run(self):
        """
        Run tasks

        :return: a list containing the resulting stuff of tasks
        """
        freeze_support()
        self.pool = Pool(self.cores)
        res = list(self.pool.imap(self._runstar, self.tasks))
        self.pool.close()
        return res
