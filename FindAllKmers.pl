#!/usr/bin/perl -w
# A perlscript written by Joseph Hughes, University of Glasgow
# use this perl script to find all the k-mers from a multifasta file 

use strict;
use Bio::SeqIO;
my $verbose="1";

if (!$ARGV[0]){
  print "USAGE: FindAllKmers.pl your_fasta_file k output_file_name\n";
}else{
  my $fasta_in = $ARGV[0];
  my $k = $ARGV[1]; 
  my $in_seq_num = 0;
  my $inseq = Bio::SeqIO->new( -file => "<$fasta_in",
					 -format => 'fasta');
  open (FASTAOUT, ">$ARGV[2]") || die "Can not open temp fasta file:\n $ARGV[2]\n";

  while (my $seq = $inseq->next_seq) {
	# Calculate base cooridate data
	my $id=$seq->id();
	my $seq_len = $seq->length();
	my $max_start = $seq->length() - $k;
	# Print some summary data
	print STDERR "\n==============================\n" if $verbose;
	print STDERR "SEQ LEN: $seq_len\n" if $verbose;
	print STDERR "MAX START: $max_start\n" if $verbose;
	print STDERR "==============================\n" if $verbose;
	# CREATE FASTA FILE OF ALL K LENGTH OLIGOS
	# IN THE INPUT SEQUENCE
	print STDERR "Creating oligo fasta file\n" if $verbose;
	for my $i (0..$max_start) {
	  my $start_pos = $i + 1;
	  my $end_pos = $start_pos + $k - 1;
	  my $oligo = $seq->subseq($start_pos, $end_pos);
	  # Set counts array to zero
	  print FASTAOUT ">$id\_$start_pos\n";
	  print FASTAOUT "$oligo\n";
	}
  }
  close (FASTAOUT);
}