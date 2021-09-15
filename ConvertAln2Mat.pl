#!/usr/bin/perl -w
# A perlscript written by Joseph Hughes, University of Glasgow
# use this perl script to convert a multiple sequence alignment into a tex-tab delimited matrix with positions as column names


use strict;
use Getopt::Long; 
use Bio::SeqIO;                    
 
my ($in,$out,$help,$gene);
&GetOptions(
	    'in:s'      => \$in,#fastafile
	    'out:s'   => \$out,#output text-tabe matrix
	    'gene:s'  => \$gene,#gene label
	    "help"  => \$help,  # provides help with usage
           );
if (($help)&&!($help)||!($in)||!($out)){
 print "Usage : ConvertAln2Mat.pl <list of arguments>\n";
 print " -in <txt> - the input fasta file\n";
 print " -out <txt> - the name of your output text-tab delimited matrix\n";
 print " -gene <txt> - gene label\n";
 print " -help        - Get this help\n";
 exit();
 }

my $inseq  = Bio::SeqIO->new(-file => "$in" ,
                         -format => 'fasta');

open(OUT,">$out")||die "Can't open $out\n";
my $total_sites=0;
my %seqs;
while (my $seq = $inseq->next_seq()){
  my $id=$seq->id;
  my $seq_str=$seq->seq();
  $seqs{$id}=$seq_str;
  $total_sites=length($seq_str);
}
print OUT "sequence_name";
for (my $i=0; $i<$total_sites; $i++){
  if ($gene){ 
    print OUT "\t".$gene."position".($i+1);
  }else{
    print OUT "\tposition".($i+1);
  }
}
print OUT "\n";

for my $id (keys %seqs){
  my @bases=split(//,$seqs{$id}); #the input can also be a protein alignment so these could equally be amino acids
  print OUT $id."\t".join("\t",@bases)."\n";
}