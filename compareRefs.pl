#!/usr/bin/perl -w
# A perlscript written by Joseph Hughes, University of Glasgow
# use this perl script to compare two reference sequences in a fasta file

use strict;
use Getopt::Long; 
use Bio::SeqIO;

my ($in,$out,$help,$r);
&GetOptions(
	    'in:s'      => \$in,#a multifasta sequence file 
	    'out:s'      => \$out,#output file 
	    "help"  => \$help,  # provides help with usage
           );

if (($help)&&!($help)||!($in)){
 print "Usage : compareRefs.pl <list of arguments>\n";
 print " -in <txt> - fasta files\n";
 print " -out <txt> - output fasta files\n";
 print " -help        - Get this help\n";
 exit();
 }

my $seqio_obj = Bio::SeqIO->new(-file => "$in", -format => 'fasta' );
my $seq_out;
if ($out){
 $seq_out = Bio::SeqIO->new(-file => ">$out", -format => 'fasta' );
}

my %seq;
my $seq1 = $seqio_obj->next_seq;
my $seq2 =  $seqio_obj->next_seq;

my $seq_str1=$seq1->seq;
my $seq_str2=$seq2->seq;
my $cnt=0;
if ($r){
  my @nucs1=split(//,$seq_str1);
  my @nucs2=split(//,$seq_str2);
  for (my $i=0; $i<scalar(@nucs1); $i++){
      if ($nucs1[$i] ne $nucs2[$i]){
        $cnt++;
      }
  }
  $seq1->seq($seq_str1);
  $seq2->seq($seq_str2);
  $seq_out->write_seq($seq1);
  $seq_out->write_seq($seq2);
}

print "Total number of differences $cnt\n";