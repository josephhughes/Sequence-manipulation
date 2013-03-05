#!/usr/bin/perl

# use this perl script to trim down a fasta alignment file

use Bio::AlignIO;
use strict;
use Getopt::Long; 

my ($help,$first,$last,$trim,$in,$out,$cut);
&GetOptions(
	    'h|help|?'  => \$help,#het the help information
	    'f:s'       => \$first, #first base to keep
	    'l:s'         => \$last, #last base to keep
	    't:s'        => \$trim, #trim bases/aa from end of read
	    'c:s'       =>\$cut, #region to cut 102-106
	    '-i:s'     => \$in, #input fasta file
	    '-o:s'   => \$out, #output fasta file
           );

if (($help)||!$in||!$out){
  #print "usage: TrimAlignment [-h] [-f N] [-l N] [-t N] [-m MINLEN] [-z] [-v] [-i INFILE] [-o OUTFILE]\n";
  print "usage: TrimAlignment [-h] [-f N] [-l N] [-t N] [-c N:N] [-i INFILE] [-o OUTFILE]\n";
  print "Part of Sequence-manipulation by J. Hughes (joseph(dot)hughes(at)glasgow(dot)ac(dot)uk\n\n";

  print "   [-h]         = This helpful help screen.\n";
  print "   [-f N]       = First base to keep. Default is 1 (=first base).\n";
  print "   [-l N]       = Last base to keep. Default is entire read.\n";
  print "   [-t N]       = Trim N nucleotides from the end of the read.\n";
  print "                  '-t'  can not be used with '-l' and '-f'.\n";
  #print "   [-m MINLEN]  = With [-t], discard sequences shorter than MINLEN.\n";
  print "   [-c N:N]       = Cut-out a region e.g. bases 102:106\n";
  print "   [-i INFILE]  = FASTA input file.\n";
  print "   [-o OUTFILE] = FASTA output file.\n";
  exit();
}else{
  my $inaln  = Bio::AlignIO->new(-file => "$in" ,
                         -format => 'fasta');
  my $outaln = Bio::AlignIO->new(-file => ">$out",
                         -format => 'fasta');
 while ( my $aln = $inaln->next_aln ) { 
  # Manipulate
  #$aln->remove_seq($seq);
  my $mini_aln;
  if ($first && $last){
    $mini_aln = $aln->slice($first,$last);  # get a block of columns
  #$mini_aln = $aln->select_noncont(1,3,5,7,11); # select certain sequences
  }
  if ($cut){
    my @coord=split(/:/,$cut);
    $mini_aln = $aln->remove_columns([$coord[0],$coord[1]]); # remove by position
  }
  #$new_aln = $aln->remove_columns(['mismatch']); # remove by property
   $outaln->write_aln($mini_aln); 
 }
}