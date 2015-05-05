# Requirements #

Corset requires [samtools](http://samtools.sourceforge.net/) to be installed. It will use the library and headers files, so these must be in a place that can be accessed by Corset.

We also suggest using gcc version 4.3 or greater because some newer container classes
are available in c++0x. This is not a requirement but code compiled with older
versions of gcc will run slower.

# Installation #

To install, please down load the tar ball from [here](Download.md). Unzip and untar the source code.

Run configure:
```
./configure
```

Note that you may have to specify the directory location of sam.h and libbam with
```
--with-bam_inc=<directory_containing_sam.h>
```
and
```
--with-bam_lib=<directory_containing_bam_library>
```
if they are not in the usual paths. **Note that you need to give the absolute path
and not a relative path**.

Then run:
```
make
make install
```


# Usage #

In the simplest case, Corset can be run in the directory containing your bam files simply by typing:
```
corset *.bam
```

The inputs should be one bam files for each sample. The bam files should have been produced by multi-mapping the reads to the transcriptome. For example with bowtie/bowtie2 you should use the parameter --all (or -k with a large number).

The usage information provided by corset is:
```

Usage: corset [options] <input bam files>

input bam files:
	 The input files should be multi-mapped bam files. They can be single, paired-end or mixed
	 and do not need to be indexed. A space separated list should be given.
	 e.g. corset sample1.bam sample2.bam sample3.bam
	 or just: corset sample*.bam

	 If you want to combine the results from different transcriptomes. i.e. the same reads have
	 been mapped twice or more, you can used a comma separated list like below:
	 corset sample1_Trinity.bam,sample1_Oases.bam sample2_Trinity.bam,sample2_Oases.bam ...

options are:

	 -d <double list> A comma separated list of distance thresholds. The range must be
	                  between 0 and 1. e.g -d 0.4,0.5. If more than one distance threshold
	                  is supplied, the output filenames will be of the form:
	                  counts-<threshold>.txt and clusters-<threshold>.txt
	                  Default: 0.3

	 -D <double>      The value used for thresholding the log likelihood ratio. The default
	                  value will depend on the number of degrees of freedom (which is 1 - the
	                  number of groups). By default D = 17.5 + 2.5 * ndf, which corresponds
	                  approximately to a p-value threshold of 10^-5, when there are fewer than
	                  10 groups.

	 -m <int>         Filter out any transcripts with fewer than this many reads aligning.
	                  Default: 10

	 -g <list>        Specifies the grouping. i.e. which samples belong to which experimental
	                  groups. The parameter must be a comma separated list (no spaces), with the
	                  groupings given in the same order as the bam filename. For example:
	                  -g Group1,Group1,Group2,Group2 etc. If this option is not used, each sample
	                  is treated as an independent experimental group.

	 -p <string>      Prefix for the output filenames. The output files will be of the form
	                  <prefix>-counts.txt and <prefix>-clusters.txt. Default filenames are:
	                  counts.txt and clusters.txt

	 -f <true/false>  Specifies whether the outputfiles should be overwritten if they already exist.
	                  Default: false

	 -n <string list> Specifies the sample names to be used in the header of the output count file.
	                  This should be a comma separated list without spaces.
	                  e.g. -n Group1-ReplicateA,Group1-ReplicateB,Group2-ReplicateA etc.
	                  Default: the input filenames will be used.
```




# Output #

By default corset will output two files: `clusters.txt` and  `counts.txt`. If you have specified multiple
distance thresholds, then the output will be of the form `clusters-<threshold>.txt` and  `counts-<threshold>.txt`.

`clusters.txt` is a tab delimited table with one line for each transcript. The first column contains the transcript ids and the second column is the cluster id it has been assigned to.

`counts.txt` is also a tab delimited table. It lists the number of reads assigned to each cluster, one
per row. There is one columns for each sample.

The cluster naming is of the form `Clusters-X.Y`. The `X` is the _super-cluster_ ID. Any transcript which shares even a single read with another transcript will have the same _super-cluster_ ID. The `Y` indicates the cluster number within the _super-cluster_ (ie. those which resulted from the hierarchical clustering and expression testing. If you run Corset with the option `-m 0` (no filtering on the number of reads), you might also find clusters with the prefix "NoReadsCluster". These clusters have no reads which map to them, and are therefore excluded from the counts file.