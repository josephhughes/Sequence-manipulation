#!/usr/bin/perl

use strict;
use Bio::SeqIO;
use Getopt::Long;


my ($in,$out,$unknown,$help);
my $unknown="X";
GetOptions(
       'h|help|?'  => \$help,#het the help information
	   'in:s'  => \$in,
	   'out:s'  => \$out,
	   'u|unknown:s' =>\$unknown,
	   );

if (($help)||!$in||!$out){
  print "usage: translate_seq [-h] [-in text] [-out text] [-unknown Character]\n";
  print "Part of Sequence-manipulation by J. Hughes (joseph(dot)hughes(at)glasgow(dot)ac(dot)uk\n\n";

  print "   [-h]            = This helpful help screen.\n";
  print "   [-in text]      = Multi-fasta file name input\n";
  print "   [-out text]     = File name for the protein output.\n";
  print "   [-unknown char] = Single character to use for unknown amino acid translations, default is X.\n";
  exit();
}


my $seqin = Bio::SeqIO->new( -format => 'fasta', -file => $in); 
my $seqout = Bio::SeqIO->new( -format => 'fasta', -file => ">$out" );

while( (my $seq = $seqin->next_seq()) ) {
  #print $seq->seq."\n";
  my $pseq = $seq->translate(-unknown => $unknown);
  #print "$pseq\n";
  $seqout->write_seq($pseq);
}

