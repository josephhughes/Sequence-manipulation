#!/usr/bin/perl -w
#

# use this to rename the identifiers in a fasta file according to the names in a table
# you can specify the column for the current name and the column for the new name
# and specify the output file name

use strict;
use Getopt::Long; 
use Bio::SeqIO;


my ($in,$out,$table, $col, $newcol,%newnames,%seen);
&GetOptions(
	    'in:s'      => \$in,#the input fasta file
	    'table:s'     => \$table,#input table with the names (tabl delimited)
	    'out:s'   => \$out,#the output file name
	    'col:s'  => \$col, #the column containing the old name
	    'newcol:s'  => \$newcol, #the column containing the new names
           );

$col=$col-1;
$newcol=$newcol-1;
open(TABLE,"$table")||die "Can't open input $table!\n";
while (<TABLE>){
  chomp($_);
  my @line=split(/\t/,$_);
  $newnames{$line[$col]}=$line[$newcol];
}


my $inseq  = Bio::SeqIO->new(-file => "$in" ,
                         -format => 'fasta');
my $outseq = Bio::SeqIO->new(-file => ">$out" , '-format' => 'fasta');
#           my $newseq = Bio::Seq->new(-seq => "$uniqseq",  
#                          -desc => "$final_hash{$uniqseq}",
#                          -display_id => $new_id);
#                          # -desc => "$final_hash{$uniqseq}",
#          $outseq->write_seq($newseq);



while ( my $seq_obj = $inseq->next_seq() ) {
    my $id=$seq_obj->display_id;
    $id=~s/cds://g;
    $id=~s/\/\d+-\d+//g;
    
    if ($newnames{$id}){
      $id=$newnames{$id};
      if ($seen{$id}){
        print "SEEN $id already\n";
      }
    }else{
      print "No name for $id\n";
    }
    $seq_obj->display_id($id);
    $seen{$id}++;
    $outseq->write_seq($seq_obj);
}


