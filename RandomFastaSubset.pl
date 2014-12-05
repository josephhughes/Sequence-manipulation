#!/usr/bin/perl -w
# A perlscript written by Joseph Hughes, University of Glasgow
# use this perl script to select a random subset of sequences from a multi fasta file

use strict;
use Getopt::Long; 
use Bio::SeqIO;                    
 
my ($infile, $outfile, $nbsub, $help);
&GetOptions(
	    'in:s'      => \$infile,#the multifasta file to subsample
	    'out:s'   => \$outfile,#the subset of random sequences
	    "sub:s"  => \$nbsub, #number of sequences to sample
	    "help"  => \$help,  # provides help with usage
           );

if (($help)&&!($help)||!($infile)||!($outfile)){
 print "Usage : RandomFastaSubset.pl <list of arguments>\n";
 print " -in <txt> - the input multifasta\n";
 print " -out <txt> - the name of your output fasta file\n";
 print " -sub <txt> - the number of sequences to subsamplings\n";
 print " -help        - Get this help\n";
 exit();
 }
my $in  = Bio::SeqIO->new(-file => "$infile" ,-format => 'fasta');  
my $out = Bio::SeqIO->new(-file => ">$outfile" , '-format' => 'fasta');

my (%sequences);
my $cnt=0;
while (my $seq = $in->next_seq()){
  my $id = $seq->id();
  $sequences{$cnt}=$seq;
  print "$sequences{$cnt}\n";
  $cnt++;
}

my @array = 0 .. ($cnt-1);
my (@num);
for (my $i=0; $i<$nbsub; $i++){
  push @num, splice(@array, int rand @array, 1);
}
print scalar(@num);  
foreach my $num (@num){
  my $rdmseq=$sequences{$num};
  $out->write_seq($rdmseq);
}
