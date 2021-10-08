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
use Bio::SeqIO;

my $ext=0;
my ($help,$in,$start,$stop,$ref,$out,$all,$join); 
&GetOptions(
	    'in:s'  => \$in, #multifasta file
	    'out:s'  => \$out, #multifasta file
	    "ref:s" => \$ref,#the identifier of the sequence to extract from
	    'start:s'  => \$start, #the start site
	    'stop:s'  => \$stop, #the stop site
	    'ext:i'  => \$ext, #number of bases to extend by at either end
	    'all'  => \$all,
	    'join:s' => \$join,
        "help:s" => \$help,
           );

if (($help)&&!($help)||!($in)){
 print "Usage : perl ExtractRegionFromFasta.pl -in 43346_ref_Bison_UMD1.0_chrUn.fa -ref \"gi|735089004|ref|NW_011494785.1|\" -start 12075790 -stop 12090981 -out ISG96region_Bison.fa \n";
 print " -in <txt>  - multifasta file\n";
 print " -ref <txt>  - the identifier of the sequence to extract from\n";
 print " -start <txt>  - the start site\n";
 print " -stop <txt>  - the stop site\n";
 print " -ext <txt>  - the number of bases to extend by at either end\n";
 print " -out <txt>  - fasta output for the region of interest\n";
 print " -all - output all the sequences for the region of interest rather than default which is just for the ref identifer\n";
 print " -join <txt> - multi-region that needs to be joined, e.g.: 266..13468,13468..21555";
 print " -help        - Get this help\n";
 exit();
}


if ($all){
  my $infile  = Bio::SeqIO->new(-file => "$in" ,
                         -format => 'fasta');
  if ($out){
    open(OUT,">$out")|| die "Can't open $out\n";
	while (my $seq = $infile->next_seq()){
	 my $id=$seq->id;
     my $sub_seq = $seq->subseq($start,$stop);
     print OUT ">$id\n$sub_seq\n";
    }
  }else{
	while (my $seq = $infile->next_seq()){
	  my $id=$seq->id;
#     $seq_str=$seq->seq();
#     my $sub_seq  = substr $seq_str, ($start-$ext), ($stop=$start+1); 
     my $sub_seq = $seq->subseq($start,$stop);
     print ">$id\n$sub_seq\n"; 
    }
  }
}elsif ($ref){
  my $db = Bio::DB::Fasta->new($in);
  my $seq = $db->seq($ref, $start-$ext => $stop+$ext);
  if ($out){
    open(OUT,">$out")|| die "Can't open $out\n";
    print OUT ">$ref\:".($start-$ext)."-".($stop+$ext)."\n$seq\n";
  }else{
    print ">$ref\:".($start-$ext)."-".($stop+$ext)."\n$seq\n";
  }
}elsif ($join){
   my @regions=split(/,/,$join);
   my $infile  = Bio::SeqIO->new(-file => "$in" ,
                         -format => 'fasta');
  if ($out){
    open(OUT,">$out")|| die "Can't open $out\n";
	while (my $seq = $infile->next_seq()){
	 my $id=$seq->id;
	 my $out_str=">$id\n";
	 foreach my $region (@regions){
	   #print "$join $region\n";
	   my ($start,$stop)=split(/\.\./,$region);
	   my $sub_seq = $seq->subseq($start,$stop);
	   $out_str=$out_str.$sub_seq;
	 }
     print OUT "$out_str\n";
    }
  }
}

