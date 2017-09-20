# STAB
Snp TAble Builder: Python tool to build a table of multiple SNP files, either as an intermediate step or as a final product

#### Usage:

python stab.py -c concatenatedList -l lowCovList -u highOrNormalCovList

outputs to STDOUT

### Requirements:

concatenatedList is a concatenated list of SNP locations of all samples to be scanned using samtools and concatenated together

lowCovList is a tab-delimited text file of low coverage SNP output, usually from DNASTAR SeqManPro, in the format:

IgnoredSymbolOrField\tSampleName\tReferenceName\tRContigPosition\tReferencePosition\tSNP_VNTR_orOther\RefBase\tCalledBase\tGenotype\tImpact\tHomopolymer\tSNP%\tP_not_ref\tQ_call\tFeature_type\tFeature_name\tOthers

highOrNormalCovList is a tab-delimited text file of normal or upper bound SNP output, usually from SeqManPro, in the same format as above

