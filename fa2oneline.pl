#!/usr/bin/perl -w

# perl script to convert a fasta file wrapped over multiple lines into a single line fasta
# Part of Sequence-manipulation by J. Hughes (joseph(dot)hughes(at)glasgow(dot)ac(dot)uk\

use strict;

my $input_fasta=$ARGV[0];
open(IN,"<$input_fasta") || die ("Error opening $input_fasta $!");

my $line = <IN>; 
print $line;

while ($line = <IN>){
  chomp $line;
  if ($line=~m/^>.+/) { 
    print "\n",$line,"\n";
  }else{ 
    print $line;
  }
}
print "\n";