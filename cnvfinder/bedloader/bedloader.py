#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 19:35:14 2016

@author: roberto
"""
from typing import Union

from bedhandler.domain import RegionId, GeneId, Pool
from pandas import DataFrame
from bedhandler.handler import BedFileLoader


def ntuple(string: str) -> tuple:
    """
    Return a tuple of integers from a string of comma-separated integers. When casting a substring to a int fails,
    it uses the substring itself

    :param str string: string of comma-separated integers
    :return: tuple of int
    """
    return tuple([maybeint(number) for number in string.strip().split(',')])


def maybeint(string: str) -> Union[str, int]:
    """
    Try to cast a string to a int

    :param str string: to cast
    :return: int(string) or string in case cast fails
    """
    try:
        return int(string)
    except ValueError:
        return string


class ROI(object):
    """
    This class stores "regions of interest" in the .bed file

    - self.amplicons is a list of Amplicon objects
    - self.targets is a DataFrame of valid targets for CNV analysis

    :param str bedfile: path to file in which sequencing amplicons are listed (BED)
    :param str region: should be in the form: chr:start-end, for example: chr1:10000-90000" or "chr1:10000" or "chr1
    :param int spacing: is the number of nucleotides to ignore at amplicon start and end, to avoid overlapping reads
    :param int mindata: is the minimum number of nucleotides for a valid target
    :param int maxpool: is the maximum number of origin pools allowed for a valid target
    """

    def __init__(self, bedfile: str, region: str = None, spacing: int = 20,
                 mindata: int = 50, maxpool: int = None):
        # spacing should be >= 0
        if not type(spacing) == int or spacing < 0:
            print("ROI: Invalid spacing: {0}; defaulting to 0".format(spacing))
            spacing = 0

        # mindata should be > 0
        if not type(mindata) == int or mindata <= 0:
            print("ROI: Invalid mindata: {0}; defaulting to 1".format(mindata))
            mindata = 1

        print("Loading .bed data.")
        if maxpool:
            print("ROI limited to max number of pools = {0}".format(maxpool))
        self.source = bedfile

        bed_file = BedFileLoader(bedfile)
        lines = bed_file.expand_columns()

        self.amplicons = self.load_amplicons(lines, region, maxpool, 8)
        if self.amplicons:
            self.targets = self.define_targets(spacing, mindata)
        else:
            return

    def __repr__(self):
        return 'ROI("{0}"; {1} amplicons, {2} targets)'.format(
            self.source, len(self.amplicons), len(self.targets))

    @staticmethod
    def load_amplicons(lines: list, region: str, maxpool: int, pool_loc: int) -> list:
        """
        Return a sorted list of amplicons

        :param list lines: list of lines loaded from bed file
        :param str region: it should be in the form: chr1:10000-90000'
        :param int maxpool: maximum number of origin pools allowed for a valid target
        :return: sorted list of amplicons
        """
        print("Loading amplicons...")
        amplicons = []
        total = 0
        skipped = 0
        for line in lines:
            newamplicon = None
            total += 1
            if pool_loc is None:
                try:
                    pool_loc = [i for i in range(len(line)) if "Pool=" in line[i]][0]
                except IndexError:
                    print("Could not find Pool data.")
                    pool_loc = False
                    break
            if pool_loc:
                newamplicon = Amplicon(line, pool_loc, region, maxpool)
            if newamplicon is not None:
                amplicons.append(newamplicon)
            else:
                skipped += 1

        print("{0} amplicons read from the .bed file.".format(total))
        print("{0} amplicons loaded.".format(len(amplicons)))
        print("{0} amplicons out of range.".format(skipped))
        print('Sorting by chromosome and range...')
        amplicons.sort(key=lambda x: (x.chromosome, x.chromStart, x.chromEnd))
        return amplicons

    def define_targets(self, spacing, mindata) -> DataFrame:
        """
        Return valid targets for further analysis. Each target is in the form: (chromosome, start, end, Amplicon)

        :param int spacing: number of nucleotides to ignore at amplicon start and end, to avoid overlapping reads
        :param int mindata: minimum number of nucleotides for a valid target
        :return: a dataframe of targets
        """
        print('Defining targets...')
        targets = []
        this_chrom = ''
        last_chrom_end = 0
        amp_size = len(self.amplicons)

        for i, amplicon in enumerate(self.amplicons):
            if amplicon.chrom != this_chrom:
                print("Analyzing chromosome: " + amplicon.chrom)
                this_chrom = amplicon.chrom
                last_chrom_end = 0

            # Skip amplicons that are contained within the last amplicon
            if amplicon.chromEnd - spacing <= last_chrom_end:
                continue

            # Avoid regions of overlap with the previous amplicon
            this_start = max(last_chrom_end + spacing,
                             amplicon.chromStart + spacing)
            this_end = amplicon.chromEnd - spacing

            # Avoid regions of overlap with the next amplicon, if applicable
            if i < (amp_size - 1) and self.amplicons[i + 1].chrom == amplicon.chrom:
                this_end = min(this_end,
                               self.amplicons[i + 1].chromStart - spacing)

            # Valid amplicons must be at least `mindata` nt long:
            if this_end - this_start >= mindata:
                targets.append((
                    this_chrom,
                    this_start,
                    this_end,
                    amplicon.genename,
                    amplicon
                ))

            last_chrom_end = amplicon.chromEnd

        print("{0} targets acquired.".format(len(targets)))
        return DataFrame(targets)


class Amplicon(object):
    """
    Hold data for each amplicon listed in the BED file

    :param list beddata: is the line in the .bed file corresponding to the amplicon
    :param int pool_loc: location of pool data in beddata
    :param str region: should be in the form: chr1:10000-90000'
    :param int maxpool: maximum number of origin pools allowed for a valid target
    """

    fields = [('chrom', str),
              ('chromStart', int),
              ('chromEnd', int),
              ('regionId', RegionId),
              ('score', maybeint),
              ('strand', str),
              ('frame', str),
              ('gene', GeneId),
              ('pools', Pool),
              ('submittedRegion', maybeint)]

    fieldnames = {field[0] for field in fields}

    def __new__(cls, beddata: list, pool_loc: int, region: str = None, maxpool: int = None):
        """
        Verify the need to create an Amplicon object

        :param list beddata: is the line in the .bed file corresponding to the amplicon
        :param int pool_loc: location of pool data in beddata
        :param str region: should be in the form: chr1:10000-90000'
        :param int maxpool: maximum number of origin pools allowed for a valid target
        :return: None or cls
        """
        location = [beddata[0], int(beddata[1]), int(beddata[2])]
        pools = beddata[pool_loc]

        # No region constraint and no # maxpool defined
        if not region and not maxpool:
            return object.__new__(cls)

        if maxpool:
            # Number of pools > maxpool -> False
            if len(pools) > maxpool:
                return None

        if region:
            # Wrong chromosome -> False
            if location[0] != region[0] or len(region) > 3:
                return None

            # Chromosome and start are defined; location is before start -> False
            elif len(region) >= 2 and max(location[1:]) < region[1]:
                return None

            # End is defined; location is after end -> False
            elif len(region) == 3 and min(location[1:]) > region[2]:
                return None

        return object.__new__(cls)

    def __init__(self, beddata: list, pool_loc: int, region: str = None, maxpool: int = None):
        # Load attribute data
        for i in range(len(beddata)):
            setattr(self, self.fields[i][0], self.fields[i][1](beddata[i]))
        self.pools = beddata[pool_loc]
        self.genename = self.gene
        chromosome = self.chrom.split('chr')[-1]
        try:
            self.chromosome = '{:0>2}'.format(int(chromosome))
        except ValueError:
            self.chromosome = chromosome

    def __getattr__(self, attribute):
        """
        Override getattr() to handle unset valid fields
        """
        if attribute in self.fieldnames:
            return None
        else:
            raise AttributeError

    def __repr__(self):
        return "Amplicon({0};{1}:{2}-{3})".format(self.name,
                                                  self.chrom,
                                                  self.chromStart,
                                                  self.chromEnd)
