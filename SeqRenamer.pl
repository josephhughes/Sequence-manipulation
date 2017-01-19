#!/usr/bin/perl

# use this perl script to rename a fasta/fastq file

use strict;
use Getopt::Long; 

my ($help,$in,$out,$fmt,$pre,$k,$table,$col,$newcol,%newnames);
&GetOptions(
	    'h|help|?'  => \$help,#het the help information
	    'i:s'       => \$in, #input in fasta or fastq format
	    'o:s'         => \$out, #last base to keep
	    'fmt:s'        => \$fmt, #the format
	    'pre:s'       =>\$pre, #the prefix to use for the sequence id
	    't:s'   => \$table, #the table to rename from
	    'col:i'  => \$col, #the column containing the old name
	    'newcol:i'  => \$newcol, #the column containing the new names

           );

if (($help)||!$in||!$out||!$fmt){
  print "usage: SeqRenamer [-h] [-i INFILE] [-o OUTFILE] [-fmt <txt>] [-pre <txt>] [-t <txt>]\n";
  print "Part of Sequence-manipulation by J. Hughes (joseph(dot)hughes(at)glasgow(dot)ac(dot)uk\n\n";

  print "   [-h]         = This helpful help screen.\n";
  print "   [-fmt <txt>]       = The format (fastq or fasta)\n";
  print "   [-pre <txt>]       = The prefix to use for the name\n";
  print "   [-t <txt>]       = An input table in tab-delimited format\n";
  print "   [-col <txt>]       = The column with the old name in the table\n";
  print "   [-newcol <txt>]       = The column with the new name in the table\n";
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
if ($table){
	$col=$col-1;
	$newcol=$newcol-1;
	open(TABLE,"$table")||die "Can't open input $table!\n";
	while (<TABLE>){
	  chomp($_);
	  my @line=split(/\t/,$_);
	  $newnames{$line[$col]}=$line[$newcol];
	  print "$line[$col]\t$line[$newcol]\n";
	}
	my $seq_name;
   if ($fmt=~/fastq/){
		for ($k=0; $k<scalar(@lines); $k++){
			if ($lines[$k]=~/^\@(\w+.+(\:|\/)\d+)$/){
				$seq_name=$1;
				my $newname=$newnames{$seq_name};
				$lines[$k]=~s/^\@$seq_name/\@$newname $seq_name/;
			}
			if ($lines[$k]=~/^\+$seq_name/){
			    my $newname=$newnames{$seq_name};
				$lines[$k]=~s/^\+$seq_name/\+$newname $seq_name/;
			}
			print OUT $lines[$k]."\n";
		}
	}
	if ($fmt=~/fasta/){
		for ($k=0; $k<scalar(@lines); $k++){
			if ($lines[$k]=~/^>(.+)/){
				my $seq_name=$1;
				my $newname=$newnames{$seq_name};
				$lines[$k]=~s/^\>\Q$seq_name\E/\>$newname $seq_name/;
				print "Rename $seq_name $newname\n";
			}
			print OUT $lines[$k]."\n";
		}
	}
}else{
    open (TABLE, ">$pre\_nametable.txt")|| die "can't open output file $pre\_nametable.txt\n";
	my $cnt=0;
	my $seq_name;
	if ($fmt=~/fastq/){
		for ($k=0; $k<scalar(@lines); $k++){
			if ($lines[$k]=~/^\@(\w+.+(\:|\/)\d+)/){
				$seq_name=$1;
				$lines[$k]=~s/^\@$seq_name/\@$pre$cnt $seq_name/;
				$cnt++;
				print TABLE "$seq_name\t$pre$cnt\n";
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
				print TABLE "$seq_name\t$pre$cnt\n";
			}
			print OUT $lines[$k]."\n";
		}
	}
    close(TABLE);
}
close(IN);
close(OUT);
