#!/usr/bin/perl -w
# A perlscript written by Joseph Hughes, University of Glasgow
# use this perl script to get the average sequence length for
# a multi-fasta file, get the range (min - max) and the total bases in the file
# total As, Cs, Ts, Gs, Ns 

use strict;
use Getopt::Long; 
use Bio::SeqIO;
use List::Util qw( min max );

my ($in,$out,$help,$r);
&GetOptions(
	    'in:s'      => \$in,#a multifasta sequence file 
	    'out:s'      => \$out,#output file text tab-delimited
	    "help"  => \$help,  # provides help with usage
           );

if (($help)&&!($help)||!($in)||!($out)){
 print "Usage : FastaStats.pl <list of arguments>\n";
 print " -in <txt> - fasta files\n";
 print " -out <txt> - output text-tab delimited\n";
 print " -help        - Get this help\n";
 exit();
 }

my $seqio_obj = Bio::SeqIO->new(-file => "$in", -format => 'fasta' );
my (%basecnt,@length);
my $totalbases=0;
my $nbseq=0;
while( my $seq = $seqio_obj->next_seq ) {
  my $seq_str=$seq->seq;
  $totalbases=$totalbases+length($seq_str);
  push(@length,length($seq_str));
  $nbseq++;
  my @bases=split(//,$seq_str);
  foreach my $base(@bases){
    $basecnt{$base}++;
  }
}

print "# seqs\tAverage length\tMin length\tMax length\tBase cnt\tTotal bases\n";
print "$nbseq\t".($totalbases/$nbseq);
print "\t".(min @length)."\t".(max @length)."\t";
foreach my $base (keys %basecnt){
  print "$base:$basecnt{$base};";
}
print "\t$totalbases\n";