.. CNVfinder documentation master file, created by
   sphinx-quickstart on Mon Sep  4 09:41:10 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

CNVfinder's documentation
--------------------------

.. toctree::
   :maxdepth: 1

   code



What is **CNVfinder**?
======================
CNVfinder is a **Python 3.x** package for copy number variation (CNV) detection on whole exome sequencing (WES) data from amplicon-based enrichment technologies.
CNVfinder was designed to please three kind of people: those familiar with coding, those more comfortable with user friendly programs, and those who fit both styles.

CNVfinder is split in several modules and the main modules are: **nrrhandler**; **vcfhandler**; and **cnvdetector**. Their name speak for themselves: nrrhandler handles number of reads in regions loaded from BAM files using pysam.AlignmentFile, vcfhandler handles VCF files using pysam.VariantFile, and finally, cnvdetector puts all modules together for CNV detection.

Installation
============
CNVfinder and its dependencies can be easily installed using the `Python Package Index <https://pypi.python.org/pypi/pip>`_ (PIP)::

	pip3 install -U cnvfinder

Dependencies
============
All dependencies are installed when using the command above, however, if you want to install/use CNVfinder in other ways, you should install the lastest version of the distributions:

* **scipy**: ``pip3 install -U scipy``
* **scikit-learn**: ``pip3 install -U scikit-learn``
* **pysam**: ``pip3 install -U pysam``
* **plotly**: ``pip3 install -U plotly``
* **pandas**: ``pip3 install -U pandas``
* **numpy**: ``pip3 install -U numpy``
* **numexpr**: ``pip3 install -U numexpr``

Before we start
===============
When running the cnvdetector module, keep in mind the following assumptions:
 
* your user has permission to write in the directory where BAM files are located.

Running tests
=============
After installing CNVfinder, all you have to do is to create a config file (`here it's how to do it <#configuration-file>`_), do some preprocessing (`check it here <#preprocessing>`_), and type the following on the Python command line::

	from cnvfinder import CNVConfig
	CNVConfig('/path/to/file/filename.ini')

This will provide you with some beautiful graphics and a BED file for visualization of the detected CNVs on `IGV <http://software.broadinstitute.org/software/igv/>`_. After all, one thing to keep in mind is that BAM files are normally very large, what means that the above script will take a while to finish. However, you can use our nrrhandler module in order to speed the things up for future CNV detection tests that will be using the same BAM files (`see here <#preprocessing>`_).

Configuration file
===================
In the following snippet you can find an example of a `configuration file <https://docs.python.org/3/library/configparser.html>`_. There are four mandatory sections: bed, sample, baseline, and output. Basically, the arguments of each mandatory section represent the names of the files that can be used in the CNV detection process::

	[bed]
		bedfile = /path/to/AmpliSeqExome.20131001.designed.bed

	[sample]
		bamfile = /path/to/test/file.bam
		covfile = /path/to/test/file.tsv
		vcffile = /path/to/test/file.vcf.gz

	[baseline]
		bamfiles = 
			/path/to/baseline/file_0.bam
			/path/to/baseline/file_1.bam
			/path/to/baseline/file_2.bam
			/path/to/baseline/file_3.bam
			/path/to/baseline/file_4.bam
		covfiles =
		         /path/to/baseline/file_0.tsv
		         /path/to/baseline/file_1.tsv
		         /path/to/baseline/file_2.tsv
		         /path/to/baseline/file_3.tsv
		         /path/to/baseline/file_4.tsv
		vcffiles =
		         /path/to/baseline/file_0.vcf.gz
		         /path/to/baseline/file_1.vcf.gz
		         /path/to/baseline/file_2.vcf.gz
		         /path/to/baseline/file_3.vcf.gz
		         /path/to/baseline/file_4.vcf.gz

	[targtest]
		metric = std
		interval_range = 3
		below_cutoff = 0.7
		above_cutoff = 1.3
		minread = 30
		bins = 500
		method = chr_group

	[vartest]
		size = 400
		step = 40
		metric = IQR
		interval_range = 1.5

	[output]
		path = test/results

The sample section arguments "bamfile" and "covfile" are equivalent and just one of them must be provided. The difference is that
"covfile" already contains the number of reads found in each target of the given sample. When working with "bamfile", the "bedfile" of the "bed"
section is used to define targets from the amplicons listed in this file. Furthermore, the number of reads for each defined target will be counted
from the "bamfile". This is also valid for "bamfiles" and "covfiles" of the "baseline" section. Keep in mind that you should not mix the usage of "covfiles" and
"bamfiles". An explanation about amplicon coverage files ("covfiles") can be found `here <#amplicon-coverage-files>`_.

There are other sections and parameters that can be specified:

* **targtest**
	* **minread**: minimum number of reads for a target to be considered when calculating the sample/baseline # reads ratio
	* **below_cutoff**:
	* **above_cutoff**: these two parameters are used in order to filter out ratios that are out of range (< below_cutoff or > above_cutoff) when computing std/IQR
	* **bins**: number of bins to create when plotting the graphic for target x ratio
	* **method**: method used for bins creation (block, window, or chr_group are available)
	* **size**: when using mode 'block', the block size can be specified. Default: 200
	* **step**: when using mode 'block', the step between two consecutive blocks can be specified. Default: 10
	* **metric**: measure of variability to be used when analyzing sample ratios. Standard deviation (std) and interquartile range (IQR) are available. Default: std/IQR.
	* **interval_range**: multiplier used when applying IQR or std. Default: 3 for std and 1.5 for IQR.
	* **cnv_like_range**: multiplier used when applying IQR or std. Default: 0.7 * interval_range.
	* **maxdist**: maximum distance allowed between a cnv-like block to a cnv block in order to the cnvlike block be a valid cnv. Default: 15000000.
* **vartest**
	* **mindp**: minimum number of read depth for a variant to be considered
	* **mode**: method used to group up targets when analyzing (by target, by gene, and block are available)
	* **size**: when using mode 'block', the block size can be specified. Default: 400
	* **step**: when using mode 'block', the step between two consecutive blocks can be specified. Default: 40
	* **metric**: measure of variability to be used when analyzing variants. Interquartile range (IQR) is available.
	* **interval_range**: multiplier used when applying IQR. Default: 1.5.
	* **cnv_like_range**: multiplier used when applying IQR. Default: 0.7 * interval_range.
	* **maxdist**: maximum distance allowed between a cnv-like block to a cnv block in order to the cnv-like block be a valid cnv. Default: 15000000.


Preprocessing
=============
The idea of preprocessing when applying CNVfinder is that you can speed up the CNV detection if you can
provide pre-calculated amplicon read counts for a sample (unknown or baseline) instead of counting them
from the BAM file.   

A file containing forward and reverse counts for each amplicon ("Amplicon Coverage Summary File") is
automatically generated for each sample by the Torrent Suite's CoverageAnalysis plugin, and can be
downloaded from the Coverage Analysis Report page. This file can be used in lieu of preprocessed data,
as CNVfinder is able to recognize the file structure and will extract target data accordingly.  

Alternatively, **BAM preprocessing** can be achieved using the nrrhandler module::

	from cnvfinder.nrrhandler import NRR
	nrr = NRR(bedfile='path/to/file.bed', bamfile='path/to/file.bam')
	nrr.save()

The above snippet will count the number of reads in regions (defined by the BED file) found 
in the BAM file and store this information on a file named as '{bamfile name}.bam.txt'. The next time 
the {bamfile name}.bam is used (file.bam in our example), the nrrhandler will get the number of reads 
from the text file, what is way faster than extract this information from the BAM file. This is very useful 
for the creation of an exome baseline and run multiple tests using it.  

Please note that if no preprocessed coverage file is listed as input in the configuration file, this 
extraction task will be performed from BAM files during CNVfinder execution. Due to the length of the process,
we strongly advise to use preprocessed coverage data files instead of BAM files whenever a sample is expected to
be used multiple times in CNVfinder analyses.

Amplicon coverage files
=======================
Amplicon coverage files, arguments "covfile" and "covfiles", are provided as output of Ion Torrent :sup:`TM` sequencers from Thermo Fisher Scientific.
They can be obtained from the Torrent Server by going to "Completed Runs & Results", choosing a report, clicking on "coverageAnalysis.html", clicking on the link of the desired barcode,
sliding to the bottom of the page, clicking on "Download the amplicon coverage summary file". You'll notice this file is in *xls* format and you should
convert it to a tab-separated values (TSV) file. This can be done in any spreadsheet editor.

Having issues using CNVfinder?
You can contact us at our `Github project page <https://github.com/ip4-team/cnvfinder>`_ and we'll be very happy in helping you. 
