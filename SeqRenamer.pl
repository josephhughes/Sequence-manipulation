#!/usr/bin/perl

# use this perl script to rename a fasta/fastq file

use strict;
use Getopt::Long; 

my ($help,$in,$out,$fmt,$pre,$k);
&GetOptions(
	    'h|help|?'  => \$help,#het the help information
	    'i:s'       => \$in, #input in fasta or fastq format
	    'o:s'         => \$out, #last base to keep
	    'fmt:s'        => \$fmt, #the format
	    'pre:s'       =>\$pre, #the prefix to use for the sequence id
           );

if (($help)||!$in||!$out||!$fmt){
  print "usage: SeqRenamer [-h] [-i INFILE] [-o OUTFILE] [-fmt <txt>] [-pre <txt>]\n";
  print "Part of Sequence-manipulation by J. Hughes (joseph(dot)hughes(at)glasgow(dot)ac(dot)uk\n\n";

  print "   [-h]         = This helpful help screen.\n";
  print "   [-fmt <txt>]       = The format (fastq or fasta)\n";
  print "   [-pre <txt>]       = The prefix to use for the name\n";
  print "   [-i INFILE]  = FASTA input file.\n";
  print "   [-o OUTFILE] = FASTA output file.\n";
  exit();
}
open (IN, "<$in") || die "Can't open $in input file\n";
open (OUT, ">$out")|| die "can't open output file\n";

my @lines;
my $cnt=0;
while (<IN>){
  chomp($_);
  push (@lines, $_);
}
my $cnt=0;
my $seq_name;
if ($fmt=~/fastq/){
	for ($k=0; $k<scalar(@lines); $k++){
		if ($lines[$k]=~/^\@(\w+.+\d+)/){
			$seq_name=$1;
			$lines[$k]=~s/^\@$seq_name/\@$pre$cnt $seq_name/;
			$cnt++;
		}
		if ($lines[$k]=~/^\+$seq_name/){
			$lines[$k]=~s/^\+$seq_name/\+$pre$cnt $seq_name/;
		}
		print OUT $lines[$k]."\n";
	}
}
if ($fmt=~/fasta/){
	for ($k=0; $k<scalar(@lines); $k++){
		if ($lines[$k]=~/^>(.+)/){
			my $seq_name=$1;
			$lines[$k]=~s/^\>$seq_name/\>$pre$cnt $seq_name/;
			$cnt++;
		}
		print OUT $lines[$k]."\n";
	}
}


close(IN);
close(OUT);
