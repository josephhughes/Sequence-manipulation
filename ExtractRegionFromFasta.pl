#!/usr/bin/perl

##################################################################### 
# A perlscript written by Joseph Hughes, University of Glasgow
# Extract a given region from a sequence in a multi-fasta file
# Provide the identifier of the sequence and start and stop coordinates
#####################################################################

use strict;
use Getopt::Long; 
use warnings;
use Bio::DB::Fasta;

my ($help,$in,$start,$stop,$ref,$out); 
&GetOptions(
	    'in:s'  => \$in, #multifasta file
	    'out:s'  => \$out, #multifasta file
	    "ref:s" => \$ref,#the identifier of the sequence to extract from
	    'start:s'  => \$start, #the start site
	    'stop:s'  => \$stop, #the stop site
        "help:s" => \$help,
           );

if (($help)&&!($help)||!($in)||!($ref)||!($start)||!($stop)){
 print "Usage : perl ExtractRegionFromFasta.pl -in 43346_ref_Bison_UMD1.0_chrUn.fa -ref \"gi|735089004|ref|NW_011494785.1|\" -start 12075790 -stop 12090981 -out ISG96region_Bison.fa \n";
 print " -in <txt>  - multifasta file\n";
 print " -ref <txt>  - the identifier of the sequence to extract from\n";
 print " -start <txt>  - the start site\n";
 print " -stop <txt>  - the stop site\n";
 print " -out <txt>  - fasta output for the region of interest\n";
 print " -help        - Get this help\n";
 exit();
}

my $db = Bio::DB::Fasta->new($in);
my $seq = $db->seq($ref, $start => $stop);
if ($out){
  open(OUT,">$out")|| die "Can't open $out\n";
  print OUT ">$ref\:$start-$stop\n$seq\n";
}else{
  print ">$ref\:$start-$stop\n$seq\n";
}

