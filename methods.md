### Methods

The data was parsed as follows.  First, DNASTAR® SeqManNGEN® (©2017 DNASTAR) was used to align Illumina® sequenced reads to a reference file of *M. leprae*.  SeqManPro® was used to filter and call SNP's using the described criteria for two levels, one at a low coverage, and one at a normal coverage.
Then, [samtools](http://samtools.sourceforge.net) was used to call SNPs using the fasta file and the aligned bam file, thus identifying the following conditions:
* Called SNP: sequences called as SNP by SeqManPro
* low Coverage possible SNP: sequences called as SNP by SeqManPro at the low coverage filter, but not at the normal coverage filter
* low Coverage possible WT: sequences not called as SNP by SeqManPro at any filter, but were identified with sufficient coverage by SamTools
* No Info: no SNP info; a determination could not be made on SNP status

The python script then output the data in the format:

       sample1,sample2,...sampleN,SNP,SNP_info1,SNP_info2,...SNPinfoX

