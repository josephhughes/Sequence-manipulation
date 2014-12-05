#!/usr/bin/perl

#rename the same file to remove : if the identifier 

use strict;
use Getopt::Long; 

my ($infile,$outfile);
&GetOptions(
	    'in:s'      => \$infile,#sam file
	    'out:s'         => \$outfile, # sam file  
           );
open (IN,"<$infile")|| die "Can't open $infile\n";

while(<IN>) {
  if ($_=~/^(\w.{5,30)\t.+/){
    print $1;
  }
}