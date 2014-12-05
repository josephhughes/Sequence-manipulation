#!/usr/bin/perl -w
#
# use this to rename the fasta ids in a multi fasta file 
# use a table to rename (1 column is the old name, second is the newname)
#tab delimited table

use strict;
use Getopt::Long; 

my ($infile,$outfile,$table, %names);
&GetOptions(
	    'in:s'      => \$infile,#input fastafile
	    'table:s'     => \$table,#input table with the names 
	    'out:s'   => \$outfile,#output with renamed IDs
           );

my (%newname,%info);
open(TABLE,"$table")||die "Can't open input table!\n";
while (<TABLE>){
  chomp($_);
  my @line=split(/\t/,$_);
  $newname{$line[0]}=$line[1];
  $info{$line[0]}=$line[2];
}

open (FASTA,"<$infile")||die "Can't open $infile\n";
open (OUT,">$outfile")||die "Can't open $outfile\n";
while (<FASTA>){
  my $line=$_;
  for my $old_id (keys %newname){
    if ($line=~/$old_id/){
      $line=~s/$old_id/$newname{$old_id}/;
    }
  }
  print OUT "$line";
}

