#!/usr/bin/perl

use Getopt::Long;

$result = GetOptions ("input=s" => \$infile,
	"output=s" => \$outfile,
);

$outfile = "$infile\_singleline.fasta" unless $outfile;

open INFILE, "<$infile" or die $!;
open OUTFILE, ">$outfile" or die $!;
while (<INFILE>) {
  chomp($_);
  if (/^\>/) {
	print OUTFILE "\n$_\n";  
  }else{
    print OUTFILE "$_";
  }
}
