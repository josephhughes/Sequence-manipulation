#! /usr/bin/perl -w

################################################################################
# This script converts the alignment sequences in phylip sequential format to 
# fasta format 
# Author: Wenjie Deng
# usage: Phylip2Fasta.pl infile outfile
################################################################################
use strict;

my $usage = "usage: Phylip2Fasta.pl infile outfile\n";
my $infile = shift or die $usage;	# input phylip file name
my $outfile = shift or die $usage;	# output fasta file
my $unixfile = $infile."_unix";
my $count = 0;

system ("tr '\r' '\n' <$infile >$unixfile");
open(IN, $unixfile) || die "Can't open $unixfile: $!\n";
open(OUT, ">$outfile") || die "Can't open out $outfile: $!\n";

while(my $line = <IN>) {
	chomp $line;
	next if (!$line || $line =~ /\d+\s+\d+/);
	if($line =~ /^(.*)\s+(\S+)\s*$/) {
		my $name = $1;
		my $sequence = $2;
		$name =~ s/\s+$//;
		print OUT ">$name\n$sequence\n";
		$count++;
	}
}
close IN;
close OUT;
unlink $unixfile;
print "There are total $count sequences\n";