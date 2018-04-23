#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: valengo
"""
from typing import Union


class Disk(object):
    """
    This class holds information related to 'cat /sys/block/sda/stat'
    """

    def __init__(self):
        self.start = None
        self.end = None
        self.attr = ['read_ios', 'read_merges',
                     'read_sectors', 'read_ticks', 'write_ios',
                     'write_merges', 'write_sectors', 'write_ticks'
                                                      'in_flight', 'io_ticks', 'time_in_queue']

    def mark_start(self, filename: str):
        """
        This method reads a file. It's supposed to be used before mark_end method

        :param str filename: path to file
        """
        self.start = self.__readfile(filename)

    def mark_end(self, filename: str):
        """
        This method reads a file. It's supposed to be used before mark_end

        :param str filename: path to file
        """
        self.end = self.__readfile(filename)

    def show_diff(self):
        """
        This method shows the difference of the attrs of /sys/block/sda/stat file
        """
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

    @staticmethod
    def __readfile(filename: str) -> Union[list, None]:
        """
        This method reads a file

        :param str filename: path to file
        :return file's lines
        """
        try:
            with open(filename, 'r') as file:
                lines = file.readlines()
        except FileNotFoundError:
            print("{} not found".format(filename))
            return None
        else:
            return lines
