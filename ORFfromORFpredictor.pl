#!/usr/bin/perl -w
# A perlscript written by Joseph Hughes, University of Glasgow
# use this perl script to get the DNA ORF that corresponds to the translation of ORF predictor

use strict;
use Getopt::Long; 
use Bio::SeqIO;
use Bio::AlignIO;

my ($infasta,$ORFpredictor,$output,$help);
&GetOptions(
	    'in:s'      => \$infasta,#fastafile
	    'orfinfo:s'   => \$ORFpredictor,#the output from ORFpredictor that has position of reading frame
	    "out:s"  => \$output,  # output file of cds in fasta format
	    "help"  => \$help,  # provides help with usage
           );


my $in = Bio::SeqIO->new(-file => "$infasta" , '-format' => 'fasta');
my $orf = Bio::SeqIO->new(-file => "$ORFpredictor" , '-format' => 'fasta');
my $outfasta = Bio::SeqIO->new(-file => ">$output" , '-format' => 'fasta');
my (%orfsite);
while ( my $orfseq = $orf->next_seq() ) {
    	my $descr=$orfseq->desc();
    	my $id=$orfseq->id();
        $orfsite{$id}=$descr;
}
foreach my $id (keys %orfsite){
  print "$id\t$orfsite{$id}\n";
}

while (my $dnaseq = $in->next_seq()){
  my $id=$dnaseq->id();
  #$outfasta->write_seq($dnaseq);
  my @cdsinfo=split(/\s+/,$orfsite{$id});
  if ($cdsinfo=/\-/){
    my $reversed = $dnaseq->revcom;
    #$outfasta->write_seq($reversed);
    print "reversed $cdsinfo[1],$cdsinfo[2]\n";
    my $subseq = $reversed->subseq($cdsinfo[1],$cdsinfo[2]);
    my $cds = Bio::Seq->new( -seq => $subseq,
                        -id  => $id,
				        -desc => $orfsite{$id},
    );
    $outfasta->write_seq($cds);
  }elsif ($cdsinfo=/\+/){
    print "forward $cdsinfo[1],$cdsinfo[2]\n";
    my $subseq = $reversed->subseq($cdsinfo[1],$cdsinfo[2]);
    my $cds = Bio::Seq->new( -seq => $subseq,
                        -id  => $id,
				        -desc => $orfsite{$id},
    );
    $outfasta->write_seq($cds);

  }else{
    print "$id missing $cdsinfo[1],$cdsinfo[2]\n";
  }
}
  