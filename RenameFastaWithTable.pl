#!/usr/bin/perl -w
#

# use this to rename the identifiers in a fasta file according to the names in a table
# you can specify the column for the current name and the column for the new name
# and specify the output file name

use strict;
use Getopt::Long; 
use Bio::SeqIO;


my ($in,$out,$table, $col, $newcol,%newnames,%seen,$help);
&GetOptions(
	    'in:s'      => \$in,#the input fasta file
	    'table:s'     => \$table,#input table with the names (tabl delimited)
	    'out:s'   => \$out,#the output file name
	    'col:s'  => \$col, #the column containing the old name
	    'newcol:s'  => \$newcol, #the column containing the new names
	    'help' => \$help,
           );
if ($help||!$in||!$table||!$out||!$col||!$newcol){
  print "Usage: perl RenameFastaWithTable.pl -in in.fa -table metadata.txt -out renamed.fa -col 1 -newcol 3\n";
  print "-in <txt> - input fasta sequence\n";
  print "-table <txt> - metadata table with names of sequences\n";
  print "-out <fasta> - renamed output fasta sequence\n";
  print "-col <txt> - original column number with the original name\n";
  print "-newcol <txt> - column number for the new name\n";
  print "-help        - Get this help\n";
  exit();
}


$col=$col-1;
$newcol=$newcol-1;
open(TABLE,"$table")||die "Can't open input $table!\n";
while (<TABLE>){
  chomp($_);
  my @line=split(/\t/,$_);
  (my $newname=$line[$newcol])=~s/\s/_/g;
  $newnames{$line[$col]}=$newname;
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


