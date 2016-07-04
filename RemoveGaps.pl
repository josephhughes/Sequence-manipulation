#!/usr/bin/perl -w
# A perlscript written by Joseph Hughes, University of Glasgow
# use this perl script to remove all gaps in a fasta sequence


use strict;
use Getopt::Long; 
use Bio::SeqIO;

my ($in,$out,$help,$r);
&GetOptions(
	    'in:s'      => \$in,#a multifasta sequence file 
	    'out:s'      => \$out,#output file 
	    'r'     =>  \$r,
	    "help"  => \$help,  # provides help with usage
           );

if (($help)&&!($help)||!($in)){
 print "Usage : RemoveGaps.pl <list of arguments>\n";
 print " -in <txt> - fasta files\n";
 print " -out <txt> - output fasta files\n";
 print "-r  - replace gaps wirh Ns\n";
 print " -help        - Get this help\n";
 exit();
 }

my $seqio_obj = Bio::SeqIO->new(-file => "$in", -format => 'fasta' );
my $seq_out;
if ($out){
 $seq_out = Bio::SeqIO->new(-file => ">$out", -format => 'fasta' );
}

while( my $seq = $seqio_obj->next_seq ) {
  my $seq_str=$seq->seq;
  if ($r){
    $seq_str=~s/\-/N/g;
    $seq->seq($seq_str);
    $seq_out->write_seq($seq);
  }else{
    $seq_str=~s/\-//g;
    $seq->seq($seq_str);
    $seq_out->write_seq($seq);
  }
}
