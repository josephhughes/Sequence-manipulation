#!/usr/bin/perl -w 

use Bio::AlignIO;
#Bioperl for format conversions
$infile=$ARGV[0];
$output=$ARGV[1];
print "$infile\n$output\n";

#open (MYFILE, "$infile" ) || die;
#open (DATA, ">$output" )  ||die;

use Bio::AlignIO;
$in  = Bio::AlignIO->new(-file => "$infile" , '-format' => 'fasta');
$out = Bio::AlignIO->new(-file => ">$output" , '-format' => 'nexus');
    # note: we quote -format to keep older perls from complaining.

while ( my $aln = $in->next_aln() ) {
    $aln->uppercase();
    $out->write_aln($aln);
}

