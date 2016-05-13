#!/usr/bin/perl -w
#
# Convert fasta to genbank
# Part of Sequence-manipulation by J. Hughes (joseph(dot)hughes(at)glasgow(dot)ac(dot)uk

use strict;
use Getopt::Long; 
use Bio::SeqIO;

my ($infile,$outfile,$gi);
&GetOptions(
	    'in:s'      => \$infile,#genbank input 
	    'out:s'    =>\$outfile,#in fasta
	    'gi' =>\$gi,#specify that you want the gi|number|gb|accession identifier
           );

my $seq_in  = Bio::SeqIO->new(-file => "$infile" , '-format' => 'genbank');

my $seq_out = Bio::SeqIO->new( -file   => ">$outfile",
                               -format => 'fasta',
                             );

# write each entry in the input file to the output file
if ($gi){
  while (my $inseq = $seq_in->next_seq) {
    my $gi = $inseq->primary_id();
    my $descr=$inseq->desc();
    my $acc = $inseq->accession_number();
    my $new_id="gi|$gi|gb|$acc";
    $inseq->id($new_id);
    $seq_out->write_seq($inseq);
  }
}else{
  while (my $inseq = $seq_in->next_seq) {
    $seq_out->write_seq($inseq);
  }
}
