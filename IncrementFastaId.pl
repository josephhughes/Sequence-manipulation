#!/usr/bin/perl -w
# this perlscript adds a digit to the fasat id incrementally for each sequence in a file

use strict;
use Getopt::Long; 


my ($infile,$outfile,$help);
&GetOptions(
	    'in:s'      => \$infile,#input fatsa file
	    'out:s'   => \$outfile,
	    "help"  => \$help,  # provides help with usage
           );

if (($help)&&!($help)||!($infile)||!($outfile)){
 print "Usage : CountSynNonSyn.pl <list of arguments>\n";
 print " -in <txt>  - input fastta file\n";
 print " -out <txt> - output fasta file\n";
 print " -help        - Get this help\n";
 exit();
 }
open (IN,"<$infile")||die "can't open input $infile\n";
open (OUT,">$outfile")||die "can't open output $outfile\n";
my $cnt=0;
while (<IN>){ 
  my $line=$_;
  if ($line=~/>\w+/){
    $line=~s/(>\w+)/$1$cnt/;
    $cnt++;
  }
  print OUT "$line";
}