#!/usr/bin/perl -w
# A perlscript written by Joseph Hughes, University of Glasgow
# use this perl script to split and combine fasta files into the different genes
# typically used after a reference assembly

use strict;
use Getopt::Long; 
use Bio::SeqIO;


my ($fasta,$help,$sub);
&GetOptions(
	    'in:s'      => \$fasta,#a comma-separated list of multifasta sequences to rename and split by gene
	    'sub:s'   => \$sub, #string to substitute from the fasta file name to get the sample name
	    "help"  => \$help,  # provides help with usage
           );

my @files=split(/,/,$fasta);
my %sequences;
for my $file (@files){
  my $in = Bio::SeqIO->new(-file => "$file" , '-format' => 'fasta');
  while ( my $seq = $in->next_seq() ) {
    my $id = $seq->display_id();
    my $seq = $seq->seq();
    (my $newid = $file)=~s/_bwa_E_cons.fa//;
    $sequences{$id}{$newid."_".$id}=$seq;
  }
}

for my $gene (keys %sequences){
  my $out = Bio::SeqIO->new(-file => ">$gene\.fa" , '-format' => 'fasta');
  for my $newid (keys %{$sequences{$gene}}){
     my $seq = Bio::Seq->new(-seq => $sequences{$gene}{$newid}, -display_id => $newid );
     $out->write_seq($seq);
  }
}