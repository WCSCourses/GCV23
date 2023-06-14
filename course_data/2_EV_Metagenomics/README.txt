This folder contains the files necessary to run the Enterovirus metagenomics clinical scenario for the WTAC in Clinical Virology and Genomics.
There is a set of paired-end FASTQ files (Illumina MiSeq, v2, 151nt) from which the human reads have been largely removed. Suriviving reads were trimmed and cleaned with Prinseq-lite (entropy 70, derep all)
Additionally, there are three zipped protein databases for use with DIAMOND (ev.prot.diamondDB.zip & vh.prot.diamondDB.zip) and PALADIN (ev.prot.paladinDB.zip). The first two of these expand to single files with .dmnd extensions, whereas the latter expands to a set of files with the same prefix, analagous to a BWA index.
subsets.local is a small tabular file taking those members of the NCBI taxonomy dump file prot.accession2taxid.gz that will be relevant following DIAMOND or PALADIN analysis of the FASTQ files.

13/02/2020
Dr. David Bibby
UKHSA
