#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: valengo
"""

from multiprocessing import Pool
from multiprocessing import freeze_support
from multiprocessing import cpu_count


class MPPoolHandler(object):
    '''
    handle tasks multiprocessing using multiprocessing.Pool
    '''
    def __init__(self, tasks, cores=None, chunksize=None):
        '''
         Parameters:
              tasks (list): a list of tuples (tasks)
              cores (int): number of cores to be used when multiprocessing
              chunksize (int): chunksize of tasks to be send each time for a
              given process
         '''
        self.tasks = tasks
        self.cores = cores
        self.chunksize = chunksize

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

    def _runstar(self, args):
        return args[0](*args[1:len(args)])

    def run(self):
        '''
        Run tasks

        Returns:
             res (list): resulting stuff of tasks
        '''
        freeze_support()
        self.pool = Pool(self.cores)
        res = list(self.pool.imap(self._runstar, self.tasks))
        self.pool.close()
        return res
