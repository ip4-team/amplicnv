#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 19:35:14 2016

@author: roberto
"""

import pandas
import os
import contextlib
from ..utils import Region


# Helper functions
def ntuple(string):
    '''
    Return a tuple of integers from a string of comma-separated integers.
    Defaults to string(s) if this fails.
    '''
    return tuple([maybeint(number) for number in string.strip().split(',')])


def maybeint(string):
    '''
    Try to cast a string to int. Return a string if this fails.
    '''

    try:
        return int(string)
    except ValueError:
        return string


def bedwrite(filename, df):
    print('Writing targets to file: "{0}"'.format(filename))
    with contextlib.suppress(FileNotFoundError):
        os.remove(filename)
        os.remove(filename + '.idx')
    df.to_csv(filename, sep='\t', index=False, header=False)


# Main classes
class ROI(object):
    '''
    This class stores "regions of interest" in the .bed file.
    self.amplicons is a list of Amplicon objects.
    self.targets is a pandas.DataFrame of valid targets for CNV analysis.
    '''

    def __init__(self, bedfile, region=None, spacing=20,
                 mindata=50, maxpool=None):
        '''
        Load amplicons and targets for CNV analysis.
        "bedfile" is the path to the .bed file.
        "region" should be in the form: chr:start-end, for example:
        chr1:10000-90000" or "chr1:10000" or "chr1"
        "spacing" is the number of nucleotides to ignore at amplicon start and
            end, to avoid overlapping reads.
        "mindata" is the minimum number of nucleotides for a valid target.
        '''

        # spacing should be >= 0
        if not type(spacing) == int or spacing < 0:
            print("ROI: Invalid spacing: {0}; defaulting to 0".format(spacing))
            spacing = 0

        # mindata should be > 0
        if not type(mindata) == int or mindata <= 0:
            print("ROI: Invalid mindata: {0}; defaulting to 1".format(mindata))
            mindata = 1

        print("Loading .bed data.")

        region = Region(region).as_tuple
        if region:
            regstrings = ['dummy', 'pter', 'qter']
            for i, item in enumerate(region):
                regstrings[i] = item
            print("ROI limited to {0}:{1}-{2}".format(*regstrings))
        if maxpool:
            print("ROI limited to max number of pools = {0}".format(maxpool))
        self.source = bedfile

        try:
            with open(bedfile, 'r') as file:
                lines = file.readlines()
        except FileNotFoundError as error:
            print(error)
            return None
        self.amplicons = self.load_amplicons(lines, region, maxpool)
        if self.amplicons:
            self.targets = self.define_targets(spacing, mindata)
        else:
            return None

    def __repr__(self):
        return 'ROI("{0}"; {1} amplicons, {2} targets)'.format(
                self.source, len(self.amplicons), len(self.targets))

    def load_amplicons(self, lines, region, maxpool):
        '''
        Return a sorted list of amplicons.
        '''
        print("Loading amplicons...")
        amplicons = []
        total = 0
        skipped = 0
        pool_loc = None
        for line in lines:
            newamplicon = None
            if line.startswith('chr'):
                total += 1
                beddata = line.strip().split('\t')
                if pool_loc is None:
                    try:
                        pool_loc = [i for i in range(len(beddata)) if "Pool=" in beddata[i]][0]
                    except IndexError:
                        print("Could not find Pool data.")
                        pool_loc = False
                        break
                if pool_loc:
                    newamplicon = Amplicon(beddata, pool_loc, region, maxpool)
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

    def define_targets(self, spacing, mindata):
        '''
        Return a pandas DataFrame array of valid targets for CNV analysis.
        Each target is in the form:
            (chromosome, start, end, Amplicon)
        '''
        print('Defining targets...')
        targets = []
        this_chrom = ''
        last_chromEnd = 0
        ampsize = len(self.amplicons)

        for i, amplicon in enumerate(self.amplicons):
            if amplicon.chrom != this_chrom:
                print("Analyzing chromosome: " + amplicon.chrom)
                this_chrom = amplicon.chrom
                last_chromEnd = 0

            # Skip amplicons that are contained within the last amplicon
            if amplicon.chromEnd - spacing <= last_chromEnd:
                continue

            # Avoid regions of overlap with the previous amplicon
            this_start = max(last_chromEnd + spacing,
                             amplicon.chromStart + spacing)
            this_end = amplicon.chromEnd - spacing

            # Avoid regions of overlap with the next amplicon, if applicable
            if i < (ampsize-1) and self.amplicons[i+1].chrom == amplicon.chrom:
                this_end = min(this_end,
                               self.amplicons[i+1].chromStart - spacing)

            # Valid amplicons must be at least `mindata` nt long:
            if this_end - this_start >= mindata:
                targets.append((
                                this_chrom,
                                this_start,
                                this_end,
                                amplicon.genename,
                                amplicon
                                ))

            last_chromEnd = amplicon.chromEnd

        print("{0} targets acquired.".format(len(targets)))
        return pandas.DataFrame(targets)


class Amplicon(object):
    '''
    Hold data for each amplicon listed in the BED file
    '''

    # column names from: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
    fields = [('chrom', str),
              ('chromStart', int),
              ('chromEnd', int),
              ('name', str),
              ('score', maybeint),
              ('strand', str),
              ('thickStart', maybeint),
              ('thickEnd', maybeint),
              ('itemRgb', ntuple),
              ('blockCount', maybeint),
              ('blockSizes', ntuple),
              ('blockStarts', ntuple)]

    fieldnames = {field[0] for field in fields}

    def __new__(cls, beddata, pool_loc, region=None, maxpool=None):
        '''
        Verify the need to create an Amplicon object.
        '''
        location = [beddata[0], int(beddata[1]), int(beddata[2])]
        pools = [int(n) for n in beddata[pool_loc].split('Pool=')[1].split(',')]

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

    def __init__(self, beddata, pool_loc, region=None, maxpool=None):
        '''
        `bedline` is the line in the .bed file corresponding to the amplicon
        `region` is a region constraint (if any)
        '''
        # Load attribute data
        for i in range(len(beddata)):
            setattr(self, self.fields[i][0], self.fields[i][1](beddata[i]))
        self.pools = list([int(n) for n in beddata[pool_loc].split('Pool=')[1].split(',')])
        self.genename = self.name.split('_')[0]
        chromosome = self.chrom.split('chr')[-1]
        try:
            self.chromosome = '{:0>2}'.format(int(chromosome))
        except ValueError:
            self.chromosome = chromosome

    def __getattr__(self, attribute):
        '''
        Override getattr() to handle unset valid fields
        '''
        if attribute in self.fieldnames:
            return None
        else:
            raise AttributeError

    def __repr__(self):
        return "Amplicon({0};{1}:{2}-{3})".format(self.name,
                                                  self.chrom,
                                                  self.chromStart,
                                                  self.chromEnd)

if __name__ == '__main__':
    test = ROI('AmpliSeqExome.20131001.designed.bed',
               region=None,
               spacing=20,
               mindata=50)
    bedwrite('test.bed', test.targets)
    print("Done.")
