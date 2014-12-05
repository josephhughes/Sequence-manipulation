#!/usr/bin/perl -w

# use this script to count the number of mutations relative to a reference 
# sequences must be aligned


use strict;
use Bio::SeqIO;
use Getopt::Long; 

#Quickly calculate hamming distance between two strings
#
#NOTE: Strings must be same length.
#      returns number of different characters.
#see  http://www.perlmonks.org/?node_id=500235
sub mismatch_count($$) { length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }

my ($infile,$ref,$rfile);
&GetOptions(
	    'in:s'      => \$infile,#fasta file aligned
	    'ref:s'   => \$ref,#name of the reference
	    'rfile:s' =>\$rfile, # the fasta file of the reference assuming it is aligned to in multifasta
           );

print "$infile\n";
my $in = Bio::SeqIO->new(-file=>"$infile", -format=>"fasta");
my (%seqs);
while (my $obj = $in->next_seq) { 
  my $seq=$obj->seq;
  my $id=$obj->id;
  $seqs{$id}=$seq;
}
if ($rfile){
  my $refin = Bio::SeqIO->new(-file=>"$rfile", -format=>"fasta");
  while (my $obj = $refin->next_seq) { 
	my $refseq=$obj->seq;
	$ref=$obj->id;
	$seqs{$ref}=$refseq;
  }
}

my $totalbases=0;
my $totalseqs=0;
my $totalmismatches=0;
my $seqwithmismatch=0;
for my $id (keys %seqs){
  $totalbases=$totalbases+length($seqs{$id});
  my $mm=mismatch_count($seqs{$ref},$seqs{$id});
  $totalmismatches=$totalmismatches+$mm;
  if ($mm>0){
	$seqwithmismatch++;
  }
  $totalseqs++;

}
print "Total number of sequences = $totalseqs\n";
print "Total number of bases = $totalbases\n";
print "Total number of mismatches = $totalmismatches\n";
print "Total sequences with mismatches = $seqwithmismatch\n";

