Corset is a command-line software program to go from a de novo transcriptome assembly to gene-level counts. Our software takes a set of reads that have been multi-mapped to the transcriptome (where multiple alignments per read were reported) and hierarchically clusters the transcripts based on the proportion of shared reads and expression patterns. It will report the clusters and gene-level counts for each sample, which are easily tested for differential expression with count based tools such as edgeR and DESeq.

  * Downloaded the source code for [Corset beta version](Download.md) - this is a very early version of the software, so feedback is most welcome

  * Instruction for installation and usage can be found [here](InstallingRunningUsage.md).

  * A full example pipeline to go from RNA-Seq reads to differential gene expression results on non-model organisms is also provided [here](Example.md).

  * Questions should be posted to our new [google group](https://groups.google.com/forum/#!forum/corset-project)

Authors: Nadia M. Davidson and Alicia Oshlack

## News ##

**July 11th 2014**: We are finally ready to release [corset version 1.00!](Download.md) This version has some major improvements to the amount of memory used and some minor improvements for speed. Please keep giving feedback on the google group.

**July 26th 2014**: [Our paper](http://genomebiology.com/2014/15/7/410/abstract) is now out at Genome Biology