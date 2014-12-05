#!/usr/bin/perl

use strict;
use Bio::SeqIO;
use Getopt::Long;


my ($in,$out);
GetOptions(
	   'in:s'  => \$in,
	   'out:s'  => \$out,
	   );


my $seqin = Bio::SeqIO->new( -format => 'fasta', -file => $in); 
my $seqout = Bio::SeqIO->new( -format => 'fasta', -file => ">$out" );


while( (my $seq = $seqin->next_seq()) ) {
    print $seq->seq."\n";
	my $pseq = $seq->translate(-unknown => '-');
	print "$pseq\n";
	$seqout->write_seq($pseq);
}
