#!/usr/bin/perl -w
# obtained from Yu-Wei's Bioinformatics playground 
# http://yuweibioinfo.blogspot.com/2009/01/fasta-to-phylip-converter.html

use strict;

MFAtoPHYLIP($ARGV[0]);

sub MFAtoPHYLIP
{
	my $inline;
	my $outfile = "$_[0]\.phy";
	my $count = 0;
	my $len;
	my $substate = 0;
	my @subheader;
	my @subcontent;
	my $m;
	my $n;

	open (FILE, "<$_[0]");
	while (defined($inline = <FILE>))
	{
		chomp($inline);
		if ($inline =~ /^>([A-Za-z0-9.\-_:]+)/)
		{
			$subheader[$count] = $1;
			$subcontent[$count] = "";
			$count++;
		}
		else
		{
			$subcontent[$count - 1] = $subcontent[$count - 1] . " $inline";
		}
	}
	close (FILE);

	# Calculate the content length
	$n = length($subcontent[0]);
	$len = $n;
	for ($m = 0; $m < $n; $m++)
	{
		if (substr($subcontent[0], $m, 1) eq " ")
		{
			$len--;
		}
	}

	open (FILE, ">$outfile");
	print FILE "   $count    $len\n";
	for ($m = 0; $m < $count; $m++)
	{
		$len = 10 - length($subheader[$m]);
		print FILE "$subheader[$m]";
		for ($n = 0; $n < $len; $n++)
		{
			print FILE " ";
		}
		print FILE " $subcontent[$m]\n";
	}
	close (FILE);
}
