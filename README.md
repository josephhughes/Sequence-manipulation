# Sequence Manipulation tools
###### by Joseph Hughes, University of Glasgow

#### Consensus.pl 

Consensus.pl is a perl script which creates consensus sequences from a fasta alignement.
This script uses use Bioperl modules.

Usage: Consensus.pl <list of arguments>
-in <txt> - the input alignment in fasta format
-out <txt> - the directory path for the output consensus in fasta format
-t <txt> - an optional threshold parameter as a percentage (0-100), so that positions in the alignment
           with lower percent-identity than the threshold are marked by ? in the consensus
-iupac <txt> - an optional parameter to make a consensus using IUPAC ambiguity codes from DNA and RNA.
-help        - Get this help

**Example:**
`perl ../../Script/Consensus.pl -in alignment.fasta -out ../Consensus/consensus.fasta -iupac`

It is best to run this from within the directory where you have your alignment so that the 
identifier of the consensus sequence in the output does not have a path direcotry name.


#### fastaQual2fastq.pl
A perl script to convert FASTA + QUAL files to FASTQ

**Usage:** `fastaQual2fastq.pl path/to/input/file`
The input file must have the extension .fasta or .fna and the associated quality file
must have the same name with extension .qual

#### Fasta2Phylip.pl
A perl script to convert a fasta file to phylip format for running PHYML or RaxML for example.

**Usage:** `Fasta2Phylip.pl path/to/input/fastafile` 

#### SplitFasta.pl
A perl script to split a multiple fasta file into a user specified number of files. 
Useful for splitting a job to run on multiple threads or server nodes.

**Usage:** `SplitFasta.pl -in path/to/multifastafile -s int`
int is the number of files you want it split into.

#### SelectSeq.pl
A perl script to pull out individual fasta sequences from a multi fasta file based of id matching,
description matching, exact id match and randomly selecting a sequence

**Usage:** SelectSeq.pl -in path/to/multifastafile -out output/file [options]

#### gb2fa.pl
A perl script that uses bioperl and converts a genbank file containing 1 or more sequences
into a multi fasta file. 

**Usage:** `gb2fa.pl -in infile.gb -out outfile.fa`

#### fa2oneline.pl
A perl script for converting a fasta file wrapped over multiple lines into a fasta file 
with the sequence on a single line.

**Usage:** `fa2oneline.pl multiline.fa > oneline.fa`

