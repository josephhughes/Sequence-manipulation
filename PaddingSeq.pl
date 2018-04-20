#!/usr/bin/perl -w
# A perlscript written by Joseph Hughes, University of Glasgow
# use this perl script to add Ns on either end of a sequence
# specify start end and length

use strict;
use Getopt::Long; 
use Bio::SeqIO;

my ($in,$help,$out,$where,$size);
&GetOptions(
	    'in:s'   => \$in,#fasta file - multi-fasta
	    'out:s'  => \$out, #fasta output file with the padding
	    'where:s'  => \$where, # can be start, end or both
	    'size:i'    => \$size, #specify the size in nucleotides
	    "help"   => \$help,  # provides help with usage
           );

if (($help||!$in||!$out||!$where||!$size)){
 print "Usage : PaddingSeq.pl <list of arguments>\n";
 print " -in <txt>    - a multi-fasta file\n";
 print " -out <txt>   - the padded multi-fasta file\n";
 print " -where <txt> - this can be start end or both\n";
 print " -size <txt>  - the length of Ns to pad with\n";
 print " -help        - Get this help\n";
 exit();
 }

my $inseq  = Bio::SeqIO->new(-file => "$in" ,
                             -format => 'fasta');
my $outseq  = Bio::SeqIO->new(-file => ">$out" ,
                             -format => 'fasta');
while ( my $seq_obj = $inseq->next_seq() ) {  
  my $id=$seq_obj->display_id;
  my $seq=$seq_obj->seq;
  if ($where eq "start"){
    $seq="N" x $size . $seq;
    my $out_seq_obj = Bio::Seq->new(-seq        => $seq,
                         -display_id => $id);
    $outseq->write_seq($out_seq_obj);
  }
  if ($where eq "end"){
    $seq=$seq . "N" x $size;
    my $out_seq_obj = Bio::Seq->new(-seq        => $seq,
                         -display_id => $id);

    $outseq->write_seq($out_seq_obj);
  }
  if ($where eq "both"){
    $seq="N" x $size . $seq . "N" x $size;
    my $out_seq_obj = Bio::Seq->new(-seq        => $seq,
                         -display_id => $id);

    $outseq->write_seq($out_seq_obj);
  }
 
  
}

