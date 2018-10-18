import sys
from collections import defaultdict

from bedhandler.handler import BedFileLoader
from pandas import DataFrame


class CoverageFileParser:
    """
    Parses a amplicon coverage file loaded by BedFileLoader

    :param str filename: path to amplicon.cov file
    """

    def __init__(self, filename: str):
        bed_file = BedFileLoader(filename)
        print('Loading coverage data from {}'.format(filename))
        # noinspection PyProtectedMember
        if bed_file.file_type != bed_file._BedFileLoader__amplicon_cov:
            print('{} is not a valid amplicon_coverage file'.format(filename))
            sys.exit(1)
        self.targets, self.counters = self.define_targets(bed_file.expand_columns(), bed_file.columns)

    def define_targets(self, lines, columns) -> tuple:
        """
        Extract columns from lines based on columns of interest and
        split them in two entities: one DataFrame representing actual targets [chrom, chromStart, chromEnd...] and
        a list representing the number of reads for each target.

        :param list lines: actual data
        :param list columns: list of columns of interest

        :return: a DataFrame describing the targets and a list of counters
        """
        targets = []
        counters = []
        columns_map = self.create_column_map(columns)
        for line in lines:
            targets.append([line[columns_map['chrom']],
                            int(line[columns_map['chrom_start']]),
                            int(line[columns_map['chrom_end']]),
                            str(line[columns_map['gene']]),
                            line[columns_map['pools']]])
            counters.append(int(line[columns_map['total_reads']]))
        return DataFrame(targets, columns=['chrom', 'chromStart', 'chromEnd', 'gene', 'pools']), counters

    @staticmethod
    def create_column_map(columns) -> defaultdict:
        """
        Create a dict based on columns

        :param list columns: list of columns
        :return: a dictionary mapping 'columns' values and indexes
        """
        columns_map = defaultdict(int)
        for i, column in enumerate(columns):
            columns_map[column] = i
        return columns_map
