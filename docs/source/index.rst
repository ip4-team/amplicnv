.. CNVfinder documentation master file, created by
   sphinx-quickstart on Mon Sep  4 09:41:10 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to **CNVfinder's** documentation!
-----------------------------------------

.. toctree::
   :maxdepth: 2



What is **CNVfinder**?
======================
CNVfinder is a **Python 3.x** package for copy number (CNV) variation detection on whole exome sequencing (WES) data from amplicon-based enrichment technologies. Besides the `BED <https://genome.ucsc.edu/FAQ/FAQformat#format1>`_ and `BAM <http://software.broadinstitute.org/software/igv/bam>`_ (binary SAM) files, the user may specify a `VCF <http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/>`_ (Variant Call Format) file for each BAM file, allowing a more accurate analysis.
CNVfinder was designed to please three kind of people: those familiar with coding, those more comfortable with user friendly programs, and those who fit both styles.

CNVfinder is split in four main modules: **bedloader**, **nrrhandler**; **vcfhandler**; and **cnvdetector**. Their name speak for themselves: bedloader loads/writes BED files, nrrhandler handles number of reads in regions loaded from BAM files using pysam.AlignmentFile, vcfhandler handles VCF files using pysam.VariantFile, and finally, cnvdetector puts all these modules together for CNV detection.

Installation
============
CNVfinder and its dependencies can be easily installed using the `Python Package Index <https://pypi.python.org/pypi/pip>`_ (PIP)::

	pip3 install -U cnvfinder

Dependencies
============
All dependencies are installed when using pip3, however if you want to install/use CNVfinder in other ways, you should install the lastest version of the distributions:

* **scipy**: ``pip3 install -U scipy``
* **scikit-learn**: ``pip3 install -U scikit-learn``
* **pysam**: ``pip3 install -U pysam``
* **plotly**: ``pip3 install -U plotly``
* **pandas**: ``pip3 install -U pandas``
* **numpy**: ``pip3 install -U numpy``
* **numexpr**: ``pip3 install -U numexpr``

Running tests
=============
After installing CNVfinder, all you have to do is to create a config file (`here it's how to do it <#>`_), do some preprocessing (`check it here <#>`_), and type the following on the Python command line::

	from cnvfinder import CNVConfig
	CNVConfig('/path/to/file/filename.ini')

This will provide you with some beautiful graphics and a BED file for visualization of the detected CNVs on `IGV <http://software.broadinstitute.org/software/igv/>`_. After all, one thing to keep in mind is that BAM files are normally very large, what means that the above script will take a while to finish. However, you can use our nrrhandler module in order to speed the things up for future CNV detection tests that will be using the same BAM files (`see here <#>`_).

Configuration file
===================
In the following snippet you can find a example of a `configuration file <https://docs.python.org/3/library/configparser.html>`_. There are four mandatory sections: bed, sample, baseline, and output. Each of these sections has its own mandatory keys and values: bedfile, bamfile, bamfiles, and path, respectively. Basically, these aforementioned values represent the names of the files that will be used in the CNV detection process::

	[bed]
		bedfile = /path/to/AmpliSeqExome.20131001.designed.bed

	[sample]
		bamfile = 
			/path/to/test/file.bam

	[baseline]
		bamfiles = 
			/path/to/baseline/file_0.bam
			/path/to/baseline/file_1.bam
			/path/to/baseline/file_2.bam
			/path/to/baseline/file_3.bam
			/path/to/baseline/file_4.bam

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

There are other sections and parameters that can be specified:

* **targtest**
	* **minread**: minimum number of reads for a target to be considered when calculating the sample/baseline # reads ratio
	* **below_cutoff**:
	* **above_cutoff**: these two parameters are used in order to filter out rations that are out of range (< below_cutoff or > above_cutoff) when computing std/IQR
	* **bins**: number of bins to create when plotting the graphic for target x ratio
	* **method**: method used for bins creation (block, window, or chr_group are available)
* **vartest**
	* **mindp**: minimum number of read depth for a variant to be considered
* **targtest and vartest**
	* **mode**: method used to group up targets when analyzing (by target, by gene, and block are available)
	* **size**: when using mode 'block', the block size can be specified. Default: 200/400
	* **step**: when using mode 'block', the step between two consecutive blocks can be specified. Default: 10/400
	* **metric**: measure of variability to be used when analyzing sample ratios. Standard deviation (std) and interquartile range (IQR) are available. Default: std/IQR.
	* **targtest/vartest**: interval_range: multiplier used when applying IQR or std. Default: 3 for std and 1.5 for IQR.
	* **cnv_like_range**: multiplier used when applying IQR or std. Default: 0.7 * interval_range.
	* **maxdist**: maximum distance allowed between a cnvlike block to a cnv block in order to the cnvlike block be a valid cnv. Default: 15000000.

Preprocessing
=============
The idea of preprocessing when applying CNVfinder is that you can speed up the CNV detection 
if you store some information extracted from BAM and/or VCF files in the first time when using them. 
Actually, this extraction task will happen every time if you don't take this 
proprocessing step seriously! You really DO NOT want this! **BAM preprocessing** can be achieved using the nrrhandler module::

	from cnvfinder.nrrhandler import NRR
	nrr = NRR(bedfile='path/to/file.bed', bamfile='path/to/file.bam')
	nrr.save()

The above snippet will count the number of reads in regions (defined by the BED file) found 
in the BAM file and store this information on a file named as '{bamfile name}.bam.txt'. The next time 
the {bamfile name}.bam is used (file.bam in our example), the nrrhandler will get the number of reads 
from the text file, what is way faster than extract this information from the BAM file. This is very useful 
for the creation of a exome baseline and run multiple tests using it.  