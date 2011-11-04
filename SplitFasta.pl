#!/usr/bin/perl -w
# A perlscript written by Joseph Hughes, University of Glasgow
# use this perl script to split a multifasta file into user defined number of files


use strict;
use Getopt::Long; 
use Bio::SeqIO;


my ($fasta,$split,$help);
&GetOptions(
	    'in:s'      => \$fasta,#multifastafile
	    's:s'   => \$split,#number of files to split it into
	    "help"  => \$help,  # provides help with usage
           );

my $in = Bio::SeqIO->new(-file => "$fasta" , '-format' => 'fasta');
my $grep="grep ".'">"'." $fasta |wc -l";
#print "$grep\n";
my $count=`$grep`;
#print "count $count\n";
my $seqsplit=$count/$split ;
#print "nb of sequences $seqsplit\n";
my $cntseq=0;
my $i=1;
my $filename=$fasta;
$filename=~s/(.+)\.(fasta|fa|fna)$/$1\_$i\.fasta/;
#print "$filename\n";
my $out = Bio::SeqIO->new(-file => ">$filename" , '-format' => 'fasta');
while ( my $seq = $in->next_seq() ) {
 $cntseq++;
 my $batch=$seqsplit*$i;
 #print "$cntseq $batch\n";
  if ($cntseq<=($seqsplit*$i)){
    $out->write_seq($seq);
  }else{
    $i++;
    $filename=$fasta;
    $filename=~s/(.+)\.(fasta|fa|fna)$/$1\_$i\.fasta/;
    $out = Bio::SeqIO->new(-file => ">$filename" , '-format' => 'fasta');
    $out->write_seq($seq);
  }
  
}
