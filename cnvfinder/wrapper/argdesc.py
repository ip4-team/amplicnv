#!/usr/bin/env python3
"""

@author: valengo
"""
from enum import Enum, unique


@unique
class Strings(Enum):
    target = 'Path to file in which sequencing amplicons are listed'
    region = 'Limit target definition to a given region. It should be in the form: chr1:10000-90000'
    outdir = 'Output directory name'
    output = 'Output filename'
    spacing = 'Number of nucleotides to ignore at amplicon start and end, to avoid overlapping reads'
    min_data = 'Minimum number of nucleotides for a valid target'
    max_pool = 'Maximum number of origin pools allowed for a valid target'
    bamfile = 'alignment filename (bam format)'
    no_parallel = 'Disable counting target read depth in parallel'
    baseline = 'Path to baseline sample bamfile. This parameter can be passed multiple times: --baseline file1.bam ' \
               '--baseline file2.bam '
    test = 'Path to test sample bamfile'
    size = 'Block size when sliding window'
    step = 'Step size when sliding window'
    metric = 'Define which metric should be used when comparing ratios'
    interval_range = 'Value to multiply metric by'
    min_read = 'Minimum number of reads expected for valid targets'
    below_cutoff = 'Filter out data (ratios) below this cutoff'
    above_cutoff = 'Filter out data (ratios) above this cutoff'
    max_dist = 'Maximum distance allowed of a cnv-like block, to its closest cnv block, for it be a cnv as well'
    cnv_like_range = 'Value to multiply interval_range by in order to detect cnv-like (CNVs when applying looser ' \
                     'calculations) '
    bins = 'Number of bins to use when plotting ratio data'
    method = ''
    vcf = 'Path to variant file (VCF)'
    vcf_baseline = 'Path to baseline sample vcf file. This parameter can be passed multiple times: --vcf-baseline ' \
                   'file1.vcf.gz --vcf-baseline file2.vcf.gz '
    vcf_test = 'Path to test sample vcf file'
    to_filter = 'whether to apply filters on variants'
    ideo_file = 'Path to (BED) file where chromosome bands are listed. If none is passed, we use \"CytoBand\" table ' \
                'from https://genome-euro.ucsc.edu/cgi-bin/hgTables '
    ideo_chroms = 'Create ideograms only for chromosomes that are in this list'
    ideo_chrom_height = 'Height of each ideogram'
    ideo_chrom_spacing = 'Spacing between consecutive ideograms'
    ideo_fig_size = 'Resulting figure (ideogram plot) size. Default will be defined as: 12 for width and 1 * number ' \
                    'of chromosomes for height. For example: 12 1 '
    ideo_track_file = 'Path to (BED) file where genomic data is listed'
    ideo_above = 'Each track is plotted above its respective ideogram'
    ideo_below = 'Each track is plotted below its respective ideogram'
    ideo_track_height = 'Height of genomic track. Its absolute value should be smaller than the space between ' \
                        'consecutive ideograms '
    ideo_track_padding = 'Padding between the top of a genomic track and its corresponding ideogram'
    ideo_track_color = 'Track\'s color. It will be used in case \'colors\' is not defined in the corresponding ' \
                       'track\'s file '
    ideo_track_alpha = 'Alpha value used for blending'
    ideo_track_linewidths = 'Track\'s line widths'
    ideo_legend = 'Add a legend entry to the figure {label color}. This parameter can be passed multiple times: ' \
                  '--legend gain 0ab26c --legend loss cc1231 '

    # errors
    wrong_command = 'Class \"{}\" is not eligible for a command. Does it extend {} class?'
    unrecognized_command = 'Unrecognized {} command'
    parsing_error = 'Cannot parse \"{}\". This is likely a bug. Please, report it at: {}'
    get_command_error = 'Cannot find {} \"command\". This is likely a bug. Please, report it at {}'

    # usage
    usage = 'cnvfinder {} [<args>]'
    issue_url = 'https://github.com/ip4-team/cnvfinder/issues'
    command = 'command'
    command_help = 'Module to run'
    description = 'CNVfinder is a Python 3.x package for copy number (CNV) variation detection on whole exome ' \
                  'sequencing (WES) data from amplicon-based enrichment technologies '
    getting_help = 'For getting help of a specific command use: cnvfinder <command> --help'
    available_commands = 'Available commands:'
    detect_description = 'Detect copy number variation in a test sample applying read depth and variant data'
    count_description = 'Count the number of reads aligned to each target'
    compare_description = 'Compare a test sample with a baseline of samples considering read depth'
    bafcompute_description = 'Compute B-allele frequency (BAF)'
    vcfcompare_description = 'Compare a test sample with a baseline of samples considering B-allele frequency and ' \
                             'other variant data'
    bedloader_description = 'Define targets from amplicons located in a BED file'
    ideogram_description = 'Create ideograms for visualization of chromosomes and other genomic tracks'
    cfg = '\n\tpassing arguments using a configuration file.\n\tSee example: ' \
          'http://cnvfinder.readthedocs.io/en/latest/#configuration-file '
    no_rd = 'Do not apply read depth data'
    no_vcf = 'Do not apply variant data'
    cfg_file = 'Path to configuration file'

    # defaults
    default_outdir = 'cnvfinder_results'


class ArgDesc(Enum):
    # common arguments
    target = Strings.target.value
    region = Strings.region.value
    outdir = Strings.outdir.value
    output = Strings.output.value

    # bedloader arguments
    spacing = Strings.spacing.value
    min_data = Strings.min_data.value
    max_pool = Strings.max_pool.value

    # count arguments
    bamfile = Strings.bamfile.value
    no_parallel = Strings.no_parallel.value

    # compare arguments
    baseline = Strings.baseline.value
    test = Strings.test.value
    size = Strings.size.value
    step = Strings.step.value
    metric = Strings.metric.value
    interval_range = Strings.interval_range.value
    min_read = Strings.min_read.value
    below_cutoff = Strings.below_cutoff.value
    above_cutoff = Strings.above_cutoff.value
    max_dist = Strings.max_dist.value
    cnv_like_range = Strings.cnv_like_range.value
    bins = Strings.bins.value
    method = ''

    # bafcompute arguments
    vcf = Strings.vcf.value

    # vcfcompare arguments
    vcf_baseline = Strings.vcf_baseline.value
    vcf_test = Strings.vcf_test.value
    vcf_size = Strings.size.value
    vcf_step = Strings.step.value
    vcf_metric = Strings.metric.value
    vcf_interval_range = Strings.interval_range.value
    vcf_max_dist = Strings.max_dist.value
    vcf_cnv_like_range = Strings.cnv_like_range.value
    to_filter = Strings.to_filter.value

    # detect arguments
    no_rd = Strings.no_rd.value
    no_vcf = Strings.no_vcf.value

    # cfg
    cfg_file = Strings.cfg_file.value

    # ideogram
    file = Strings.ideo_file.value
    chroms = Strings.ideo_chroms.value
    chrom_height = Strings.ideo_chrom_height.value
    chrom_spacing = Strings.ideo_chrom_spacing.value
    fig_size = Strings.ideo_fig_size.value
    legend = Strings.ideo_legend.value
    above_track_file = Strings.ideo_above.value
    above_track_height = Strings.ideo_track_height.value
    above_track_padding = Strings.ideo_track_padding.value
    above_track_color = Strings.ideo_track_color.value
    above_track_alpha = Strings.ideo_track_alpha.value
    above_track_linewidths = Strings.ideo_track_linewidths.value
    below_track_file = Strings.ideo_below.value + ' '  # TODO change our solution for "strings" and so on
    below_track_height = Strings.ideo_track_height.value + ' '
    below_track_padding = Strings.ideo_track_padding.value + ' '
    below_track_color = Strings.ideo_track_color.value + ' '
    below_track_alpha = Strings.ideo_track_alpha.value + ' '
    below_track_linewidths = Strings.ideo_track_linewidths.value + ' '
