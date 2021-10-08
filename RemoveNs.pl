#!/usr/bin/perl

#  this is a perl script for removing sequences that contain Ns from a fastq file

use strict;
use Bio::SeqIO;
use Getopt::Long; 
my ($infile, $total);
################



&GetOptions(
	    'in:s'     => \$infile,#multi-fasta file
           );

my $outNs=$infile."_Ns.fa";
my $out=$infile."_clean.fa";
my $seq_in  = Bio::SeqIO->new(-format => 'fasta',-file   => $infile);
my $seq_Ns  = Bio::SeqIO->new(-format => 'fasta',-file   => ">$outNs");
my $seq_out  = Bio::SeqIO->new(-format => 'fasta',-file   => ">$out");

my $cntN=0;
my $cntClean=0;
my $total=0;
while( my $seq = $seq_in->next_seq() ) {
  my $id=$seq->display_id();
  my $seq_str=$seq->seq();
  $total++;
  if ($seq_str=~/N/i){
    $seq_Ns->write_seq($seq);
    $cntN++;
  }else{
    $seq_out->write_seq($seq);
    $cntClean++;
  }
}
print "Number of seqs $total\nNumber with Ns $cntN\nNumber clean $cntClean\n";


