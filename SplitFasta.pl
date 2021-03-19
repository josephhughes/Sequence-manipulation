#!/usr/bin/perl

# script to split up a fasta file based on a dataframe with sequence in first row 
# and group in another row which can be specified

use strict;
use Bio::SeqIO;
use Getopt::Long;


my ($in,$fasta,$colname,$help);
GetOptions(
       'h|help|?'  => \$help,#het the help information
	   'in:s'  => \$in,
	   'fasta:s'  => \$fasta,
	   'colname:s' =>\$colname,
	   );

if (($help)||!$in){
  print "usage: SplitFasta [-h] [-in text] [-out text] [-unknown Character]\n";
  print "Part of Sequence-manipulation by J. Hughes (joseph(dot)hughes(at)glasgow(dot)ac(dot)uk\n\n";

  print "   [-h]            = This helpful help screen.\n";
  print "   [-in text]      = Text-tab delimited file with column names\n";
  print "   [-fasta text]   = Big fasta file that needs to be split.\n";
  print "   [-colname text] = The column name to use for splitting by.\n";
  exit();
}

open(METADATA,"<$in")||die "Can't open $in\n";
my $header=<METADATA>;
chomp($header);
my @colnames=split(/\t/,$header);
my $col_int;
for (my $i=0; $i<scalar(@colnames);$i++){
  if ($colname==$colnames[$i]){
    $col_int=$i;
  }
}

my %grouping;
while(<METADATA>){
  chomp($_);
  my @values=split(/\t/,$_);
  $grouping{$values[0]}=$values[$col_int];
#  print "$values[0]\t$values[$col_int]\n";
}


my $seqin = Bio::SeqIO->new( -format => 'fasta', -file => $fasta); 
#my $seqout = Bio::SeqIO->new( -format => 'fasta', -file => ">$out" );

my %sequences;
while( (my $seq = $seqin->next_seq()) ) {
  my $seq_id=$seq->display_id();
  #print "$seq_id ".$grouping{$seq_id}."\n";
  if ($grouping{$seq_id}=~/./){
    $sequences{$grouping{$seq_id}}{$seq_id}=$seq->seq();#The back-slash in-front of the hash provides us with a reference to the hash
    #print $seq."\n";
  }
  #my $seqout = Bio::SeqIO->new( -format => 'fasta', -file => ">$out" );

#  $seqout->write_seq($pseq);
}

for my $group (keys %sequences){
  my $seqout = Bio::SeqIO->new( -format => 'fasta', -file => ">$group" );
  for my $seq_id (keys %{$sequences{$group}}){
    my $seq_obj = Bio::Seq->new(-seq        => $sequences{$group}{$seq_id},
                         -display_id => $seq_id );
    $seqout->write_seq($seq_obj);
  }
}

