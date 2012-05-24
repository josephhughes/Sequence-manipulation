#!/usr/bin/perl -w
# A perlscript written by Joseph Hughes, University of Glasgow
# use this perl script to format differnet fastq files

use strict;
use Getopt::Long; 
use Bio::SeqIO;   
use Bio::Seq::Quality;

my ($infile,$outfile,$iformat,$oformat);

my $result = GetOptions ("input=s" => \$infile,
	"output=s" => \$outfile,
	"if:s" => \$iformat,  #fastq-illumina fastq 'sanger' 'fastq-solexa'
	"of:s" => \$oformat  #fastq-illumina fastq 'sanger' 'fastq-solexa'
);

print "converting $iformat to $oformat\n";

my $in = Bio::SeqIO->new(-format    => "$iformat",
                           -file      => $infile);
                           
# grabs the FASTQ parser, specifies the Illumina variant
my $out = Bio::SeqIO->new(-format    => "$oformat",
                           -file      => ">$outfile");

# $seq is a Bio::Seq::Quality object
  while (my $seq = $in->next_seq) {
      $out->write_seq($seq); 
  }
