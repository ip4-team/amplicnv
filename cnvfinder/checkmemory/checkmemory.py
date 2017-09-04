#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: valengo
"""


class Disk(object):
    '''
    This class holds information related to
    cat /sys/block/sda/stat
    '''
    def __init__(self):
        self.start = None
        self.end = None
        self.attr = ['read_ios', 'read_merges',
                     'read_sectors', 'read_ticks', 'write_ios',
                     'write_merges', 'write_sectors', 'write_ticks'
                     'in_flight', 'io_ticks', 'time_in_queue']

    def mark_start(self, filename):
        '''
        This method reads a file. It's supposed to be used
        before mark_end method

        Parameters:
            filename (str): filename
        '''
        self.start = self.__readfile(filename)

    def mark_end(self, filename):
        '''
        This method reas a file. It's supposed to be used
        before mark_end

        Parameters:
            filename (str): filename
        '''
        self.end = self.__readfile(filename)

    def show_diff(self):
        '''
        This method shows the difference of the attrs of
        /sys/block/sda/stat file
        '''
        if self.start is None:
            print('Cannot diff when "start" is None')
            print('Call Disk.mark_start(filename) before diff')
        elif self.end is None:
            print('Cannot diff when "end" is None')
            print('Call Disk.mark_end(filename) before diff')
        else:
            split_start = self.start[0].split()
            split_end = self.end[0].split()

            for i in range(len(self.attr)):
                print(int(split_end[i]) - int(split_start[i])),

    def __readfile(self, filename):
        '''
        This method reads a file

        Parameters:
            filename (str): filename
        '''
        try:
            with open(filename, 'r') as file:
                lines = file.readlines()
        except FileNotFoundError as error:
            print(error)
            return None
        else:
            return lines
