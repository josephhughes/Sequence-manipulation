#!/usr/bin/perl -w
# A perlscript written by Joseph Hughes, University of Glasgow
# use this perl script to remove description from fasta file


use strict;
use Getopt::Long; 
use Bio::SeqIO;                    
 

my ($infasta,$outfile,$help);
&GetOptions(
	    'in:s'      => \$infasta,#fastafile
	    'out:s'   => \$outfile,#output fasta file
	    "help"  => \$help,  # provides help with usage
           );

if (($help)&&!($help)||!($infasta)||!($outfile)){
 print "Usage : Consensus.pl <list of arguments>\n";
 print " -in <txt> - the input fasta file\n";
 print " -out <txt> - the name of your output fasta file\n";
 print " -help        - Get this help\n";
 exit();
 }

print "Removing descr in file $infasta...\n";
my $in  = Bio::SeqIO->new(-file => "$infasta" ,
                         -format => 'fasta');
my $out = Bio::SeqIO->new(-file => ">$outfile" , '-format' => 'fasta');

while (my $seq = $in->next_seq()){
    my $id=$seq->id;
    my $seq_str=$seq->seq;  
    $seq = Bio::Seq->new(-seq => "$seq_str",  
                        -display_id => $id);

    $out->write_seq($seq);
}
