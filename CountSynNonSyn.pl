#!/usr/bin/perl -w
# A perlscript written by Joseph Hughes, University of Glasgow
# use this perl script to count the number of nonsynonynous and synonymous mutations
# per sequence in an aligned multi-fasta file and the number of substitutions per sequence
# compared to a reference sequence
# Plots the distribution counts to the terminal and the data to a file
# specify whether it is DNA or RNA sequences

use strict;
use Getopt::Long; 
use Bio::SeqIO;


my ($infasta,$out,$code,$help,$inref);
&GetOptions(
	    'in:s'      => \$infasta,#aligned multi-fastafile or a single sequence
	    'out:s'   => \$out,#printing output data
	    'ref:s'   => \$inref, #reference sequence
	    "code:s"  => \$code,  # genetic code
	    "help"  => \$help,  # provides help with usage
           );

if (($help)&&!($help)||!($infasta)||!($out)||!($code)||!($inref)){
 print "Usage : CountSynNonSyn.pl <list of arguments>\n";
 print " -in <txt>  - aligned according to reference input file in fasta format, sequences must start with the first codon position\n";
 print " -out <txt> - the name of the output file for the data\n";
 print " -code <dna/rna> - a parameter to specify the genetic code to use dna or rna.\n";
 print " -help        - Get this help\n";
 exit();
 }
my (%genetic_code,@refcodon_list,%nonsyn,%syn,%nuc);
print "You have selected the $code genetic code\n";
%genetic_code=getCode($code);

my $ref  = Bio::SeqIO->new(-file => "$inref" ,
                         -format => 'fasta');

while (my $refseq = $ref->next_seq()){
  my $refseq_str=uc($refseq->seq);  
  #@refcodon_list = split (/(\S{3})/,$refseq_str);
  @refcodon_list = $refseq_str =~ /(\S{3})/g;
}
my $in  = Bio::SeqIO->new(-file => "$infasta" ,
                         -format => 'fasta');
open (OUT, ">$out")||die "Can't open output $out\n";

while (my $seq = $in->next_seq()){
  my $seq_str=uc($seq->seq); 
  #print "$seq_str\n";
  my $nonsyncount=0;
  my $syncount=0;
  my $nuccount=0;
  #my @codon_list = split (/(\S{3})/,$seq_str);
  my @codon_list = $seq_str =~ /(\S{3})/g;
  for (my $i=0; $i<scalar(@codon_list); $i++){
    #print "$codon_list[$i] $refcodon_list[$i]\n";
   # print "test ".$genetic_code{$codon_list[$i]}." test ".$genetic_code{$refcodon_list[$i]}."\n";
    my @refsite=split(//,$refcodon_list[$i]);
    my @site=split(//,$codon_list[$i]);
    my $codoncount=0;
    for (my $j=0; $j<scalar(@refsite); $j++){
      if ($site[$j] ne $refsite[$j]){
        $nuccount++;
        $codoncount++;
      }
    }
    if ($genetic_code{$codon_list[$i]} ne $genetic_code{$refcodon_list[$i]}){
      $nonsyncount=$nonsyncount+$codoncount;
    }elsif ($genetic_code{$codon_list[$i]} eq $genetic_code{$refcodon_list[$i]}){
      $syncount=$syncount+$codoncount;
    }
  }
  $nonsyn{$nonsyncount}++;
  if ($nonsyncount>4){
    my $seq_id=$seq->id();
    print "$seq_id has $nonsyncount nonsynonymous mutations\n";
  }
  $syn{$syncount}++;
  $nuc{$nuccount}++;
}
print "Number of sequences with x synonymous mutations\n";
for my $nbsyn (keys %syn){
  print "$nbsyn $syn{$nbsyn}\n";
}
print "Number of sequences with x nonsynonymous mutations\n";
for my $nbnonsyn (keys %nonsyn){
  print "$nbnonsyn $nonsyn{$nbnonsyn}\n";
}
print "Number of sequences with x nucleotide mutations\n";
for my $nbnuc (keys %nuc){
  print "$nbnuc $nuc{$nbnuc}\n";
}

sub getCode {
my $code_type=$_[0];
print "$code_type in subroutine\n";
my %code;
  if($code_type=~/rna/){
(%code) = (
'UCA' => 'S', # Serine
'UCC' => 'S', # Serine
'UCG' => 'S', # Serine
'UCU' => 'S', # Serine
'UUC' => 'F', # Phenylalanine
'UUU' => 'F', # Phenylalanine
'UUA' => 'L', # Leucine
'UUG' => 'L', # Leucine
'UAC' => 'Y', # Tyrosine
'UAU' => 'Y', # Tyrosine
'UAA' => '*', # Stop
'UAG' => '*', # Stop
'UGC' => 'C', # Cysteine
'UGU' => 'C', # Cysteine
'UGA' => '*', # Stop
'UGG' => 'W', # Tryptophan
'CUA' => 'L', # Leucine
'CUC' => 'L', # Leucine
'CUG' => 'L', # Leucine
'CUU' => 'L', # Leucine
'CCA' => 'P', # Proline
'CAU' => 'H', # Histidine
'CAA' => 'Q', # Glutamine
'CAG' => 'Q', # Glutamine
'CGA' => 'R', # Arginine
'CGC' => 'R', # Arginine
'CGG' => 'R', # Arginine
'CGU' => 'R', # Arginine
'AUA' => 'I', # Isoleucine
'AUC' => 'I', # Isoleucine
'AUU' => 'I', # Isoleucine
'AUG' => 'M', # Methionine
'ACA' => 'T', # Threonine
'ACC' => 'T', # Threonine
'ACG' => 'T', # Threonine
'ACU' => 'T', # Threonine
'AAC' => 'N', # Asparagine
'AAU' => 'N', # Asparagine
'AAA' => 'K', # Lysine
'AAG' => 'K', # Lysine
'AGC' => 'S', # Serine
'AGU' => 'S', # Serine
'AGA' => 'R', # Arginine
'AGG' => 'R', # Arginine
'CCC' => 'P', # Proline
'CCG' => 'P', # Proline
'CCU' => 'P', # Proline
'CAC' => 'H', # Histidine
'GUA' => 'V', # Valine
'GUC' => 'V', # Valine
'GUG' => 'V', # Valine
'GUU' => 'V', # Valine
'GCA' => 'A', # Alanine
'GCC' => 'A', # Alanine
'GCG' => 'A', # Alanine
'GCU' => 'A', # Alanine
'GAC' => 'D', # Aspartic Acid
'GAU' => 'D', # Aspartic Acid
'GAA' => 'E', # Glutamic Acid
'GAG' => 'E', # Glutamic Acid
'GGA' => 'G', # Glycine
'GGC' => 'G', # Glycine
'GGG' => 'G', # Glycine
'GGU' => 'G'  # Glycine
);
}elsif ($code_type=~/dna/){
(%code) = (
'TCA' => 'S', # Serine
'TCC' => 'S', # Serine
'TCG' => 'S', # Serine
'TCT' => 'S', # Serine
'TTC' => 'F', # Phenylalanine
'TTT' => 'F', # Phenylalanine
'TTA' => 'L', # LeTcine
'TTG' => 'L', # LeTcine
'TAC' => 'Y', # Tyrosine
'TAT' => 'Y', # Tyrosine
'TAA' => '*', # Stop
'TAG' => '*', # Stop
'TGC' => 'C', # Cysteine
'TGT' => 'C', # Cysteine
'TGA' => '*', # Stop
'TGG' => 'W', # Tryptophan
'CTA' => 'L', # LeTcine
'CTC' => 'L', # LeTcine
'CTG' => 'L', # LeTcine
'CTT' => 'L', # LeTcine
'CCA' => 'P', # Proline
'CAT' => 'H', # Histidine
'CAA' => 'Q', # GlTtamine
'CAG' => 'Q', # GlTtamine
'CGA' => 'R', # Arginine
'CGC' => 'R', # Arginine
'CGG' => 'R', # Arginine
'CGT' => 'R', # Arginine
'ATA' => 'I', # IsoleTcine
'ATC' => 'I', # IsoleTcine
'ATT' => 'I', # IsoleTcine
'ATG' => 'M', # Methionine
'ACA' => 'T', # Threonine
'ACC' => 'T', # Threonine
'ACG' => 'T', # Threonine
'ACT' => 'T', # Threonine
'AAC' => 'N', # Asparagine
'AAT' => 'N', # Asparagine
'AAA' => 'K', # Lysine
'AAG' => 'K', # Lysine
'AGC' => 'S', # Serine
'AGT' => 'S', # Serine
'AGA' => 'R', # Arginine
'AGG' => 'R', # Arginine
'CCC' => 'P', # Proline
'CCG' => 'P', # Proline
'CCT' => 'P', # Proline
'CAC' => 'H', # Histidine
'GTA' => 'V', # Valine
'GTC' => 'V', # Valine
'GTG' => 'V', # Valine
'GTT' => 'V', # Valine
'GCA' => 'A', # Alanine
'GCC' => 'A', # Alanine
'GCG' => 'A', # Alanine
'GCT' => 'A', # Alanine
'GAC' => 'D', # Aspartic Acid
'GAT' => 'D', # Aspartic Acid
'GAA' => 'E', # GlTtamic Acid
'GAG' => 'E', # GlTtamic Acid
'GGA' => 'G', # Glycine
'GGC' => 'G', # Glycine
'GGG' => 'G', # Glycine
'GGT' => 'G'  # Glycine
);
}
return %code;
}