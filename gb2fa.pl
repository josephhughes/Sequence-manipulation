#!/usr/bin/perl -w
#
# Convert fasta to genbank

use strict;
use Getopt::Long; 
use Bio::SeqIO;

my ($infile,$outfile);
&GetOptions(
	    'in:s'      => \$infile,#genbank input 
	    'out:s'    =>\$outfile,#in fasta
           );

my $seq_in  = Bio::SeqIO->new(-file => "$infile" , '-format' => 'genbank');

my $seq_out = Bio::SeqIO->new( -file   => ">$outfile",
                               -format => 'fasta',
                             );

# write each entry in the input file to the output file
while (my $inseq = $seq_in->next_seq) {
   $seq_out->write_seq($inseq);
}
