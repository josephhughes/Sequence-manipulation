#!/usr/bin/perl -w
# A perlscript written by Joseph Hughes, University of Glasgow
# use this perl script to remove a set of sequences based on id matching

use strict;
use Getopt::Long; 
use Bio::SeqIO;                    
 

my ($infasta,$outfile,$help,$idfile);
&GetOptions(
	    'in:s'      => \$infasta,#fastafile
	    'out:s'   => \$outfile,#output fasta file
	    'idfile:s' => \$idfile, #text file with mutliple ids to pull out
	    "help"  => \$help,  # provides help with usage
           );
if (($help)&&!($help)||!($infasta)||!($outfile)){
 print "Usage : SelectSeq.pl <list of arguments>\n";
 print " -in <txt> - the input fasta file\n";
 print " -out <txt> - the name of your output fasta file\n";
 print " -idfile <txt> - a file with multiple ids you want to match pull out\n";
 print " -help        - Get this help\n";
 exit();
 }
if ($idfile){
  print "Looking for multiple ids from $idfile in fasta file $infasta...\n";
  open(IDS,"<$idfile")||die "Can't open $idfile\n";
  my (%list);
  my $idcnt=0;
  my $hits=0;
  while (<IDS>){
    chomp($_);
    my @elements=split(/\t/,$_);
    my $id=$elements[0];
    $id=~s/\.\d+$//;
    #print "$id\n";
    $list{$id}++;
    $idcnt++;
  }


	my $in  = Bio::SeqIO->new(-file => "$infasta" ,
							 -format => 'fasta');
	my $outid = Bio::SeqIO->new(-file => ">In$outfile" , '-format' => 'fasta');
	my $out = Bio::SeqIO->new(-file => ">Not$outfile" , '-format' => 'fasta');
	while (my $seq = $in->next_seq()){
	  my $id=$seq->id;
  
	  # remove the accession version
	  $id=~s/\.\d+$//;
	 ## foreach my $listid (keys %list){
		if ($list{$id}){
		 $outid->write_seq($seq);
		 $hits++;
		 delete $list{$id};
		}else{
		  $out->write_seq($seq);
		}
	}
}