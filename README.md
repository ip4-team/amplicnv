**Here you'll find which steps are necessary to get CNVfinder up and running.**

# CNVfinder #

* CNVfinder is a Python 3.x package for copy number (CNV) variation detection on whole exome sequencing (WES) data from amplicon-based enrichment technologies. Besides the [BED](https://genome.ucsc.edu/FAQ/FAQformat#format1) and [BAM](http://software.broadinstitute.org/software/igv/bam) (binary SAM) files, the user may specify a [VCF](http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/) (Variant Call Format) file for each BAM file, allowing a more accurate analysis.   
* Version 1.0.0b3

## How do I get set up? ##
### Install with pip ###

```
#!shell

pip3 install cnvfinder
```

### How to run tests ###
CNVfinder was designed to please three kind of people: those familiar with coding, those more comfortable with user friendly programs, and those who fit both styles. That being said, you just need to know who you are and read the sections that best suit your needs.

CNVfinder is split in four main modules: [bedloader](#beloader-module); [nrrhandler](#nrrhandler-module); [vcfhandler](#vcfhandler-module); and [cnvfinder](#cnvfinder-module). Their name speak for themselves:
**bedloader** loads/writes BED files, **nrrhandler** handles **n**umber of **r**eads in **r**egions loaded from BAM files using [pysam.AlignmentFile](http://pysam.readthedocs.io/en/latest/usage.html#opening-a-file), **vcfhandler** handles VCF files using [pysam.VariantFile](http://pysam.readthedocs.io/en/latest/usage.html#working-with-vcf-bcf-formatted-files), and finally, cnvfinder puts all these modules together for CNV detection. 
