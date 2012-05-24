#!/usr/bin/perl

use Getopt::Long;
$ascii = 64;

$result = GetOptions ("input=s" => \$infile,
	"output=s" => \$outfile,
	"ascii:i" => \$ascii
);

$outfile = "$infile.fasta" unless $outfile;

open INFILE, "<$infile" or die $!;
open FASTA, ">$outfile" or die $!;
open QUAL, ">$outfile.qual" or die $!;

while (<INFILE>) {
	if (/^@/) {
		tr'@:'>_';
		my $fasta_header = $_;
		print FASTA $fasta_header;
		print QUAL  $fasta_header;

		$seq = <INFILE>; # get sequence in next line
		print FASTA $seq;

		<INFILE>; # skip qual header
		$qual = <INFILE>; # get qual string;
		chomp $qual;
		
		@quals = split '', $qual; # split qual string into characters
		foreach (@quals) {
			$qval = eval (ord ($_) - $ascii);
			$qval = int (10 * log(1 + 10 ** ($qval / 10.0)) / log(10)) if $ascii == 64; # convert solexa quality to phred quality
			print QUAL $qval . " ";
		}
		print QUAL "\n";
	}
}
