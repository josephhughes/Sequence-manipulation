#!/usr/bin/perl -w

# Remove the sequences that are shorter than a certain length and plot distribution of lengths etc...
# perl ~/Documents/Repo/Sequence-manipulation/SeqLengthDistribution.pl -inseq ~/Downloads/matrix_nogaps.fasta -out matrix_length -lower 397 -upper 395

use strict;
use Getopt::Long; 
use Bio::SeqIO;

my $lower=0;
my $sum_len=0;
my $new_sum_len=0;
my ($inseq,%freq,%newfreq,$outfile,$upper,$seq_out,$size,$help);
&GetOptions(
	    'inseq:s'   => \$inseq, #the sequences
	    'out:s'   => \$outfile,#output fasta file
	    'size:i' =>\$size,#length of sequences to select
	    'upper:i'   => \$upper,#cutoff size for sequences
	    'lower:i'   => \$lower,#lower cutoff size for sequences
	    'help'  => \$help,  # provides help with usage
           );
if (($help)&&!($help)||!($inseq)){
 print "Usage : SelectSeq.pl <list of arguments>\n";
 print " -inseq <txt> - the input fasta file\n";
 print " -out <txt> - the name of your output fasta file\n";
 print " -size <int> - length of sequences to select\n";
 print " -upper <int> - upper cut-off size for sequences\n";
 print " -lower <txt> - lower cut-off size for sequences\n";
 print " -help        - Get this help\n";
 exit();
 }


my $seqio_obj = Bio::SeqIO->new(-file => "$inseq", -format => 'fasta' );

if ($outfile){
 $seq_out = Bio::SeqIO->new(-file => ">$outfile", -format => 'fasta' );
}
while( my $seq = $seqio_obj->next_seq ) {
  my $length=$seq->length;
  $freq{$length}++;
  $sum_len=$sum_len+$length;
  if ($outfile && $upper && !$lower){
    if ($length>$upper){
      $seq_out->write_seq($seq);
      $newfreq{$length}++;
      $new_sum_len=$new_sum_len+$length;
    }
  }
  if ($outfile && !$upper && $lower){
    if ($length<$lower){
      $seq_out->write_seq($seq);
      $newfreq{$length}++;
      $new_sum_len=$new_sum_len+$length;
    }
  }
  if ($outfile && $upper && $lower){
    if ($length<$lower && $length>$upper){
      $seq_out->write_seq($seq);
      $newfreq{$length}++;
      $new_sum_len=$new_sum_len+$length;
    }
  }
  if ($outfile && $size){
    if ($length==$size){
      $seq_out->write_seq($seq);
      $newfreq{$length}++;
      $new_sum_len=$new_sum_len+$length;
    }
  }
}

print "x sequences with y length\n";
foreach my $seqnb (sort {$a<=>$b} keys %freq){
  print "$seqnb\t$freq{$seqnb}\t";
#   for (my $i=0; $i<$freq{$seqnb}; $i++){
#     print ":";
#   }
  print "\n";
}

print "Sum of lengths is $sum_len\n";
print "New distribution:\nx sequences with y length\n";
foreach my $seqnb (sort {$a<=>$b} keys %newfreq){
  print "$seqnb\t$freq{$seqnb}\t";
#   for (my $i=0; $i<$freq{$seqnb}; $i++){
#     print ":";
#   }
  print "\n";
}

print "Sum of lengths in the newfile is $new_sum_len\n";
