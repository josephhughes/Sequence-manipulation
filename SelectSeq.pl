#!/usr/bin/perl -w
# A perlscript written by Joseph Hughes, University of Glasgow
# use this perl script to select individual fasta file based on id or description matching
# or randomly select a file


use strict;
use Getopt::Long; 
use Bio::SeqIO;                    
 

my ($infasta,$outfile,$fid,$fdescr,$rdm,$help,$mid,$idfile,$singletons,$frequency,$longest);
&GetOptions(
	    'in:s'      => \$infasta,#fastafile
	    'out:s'   => \$outfile,#output fasta file
	    'fid:s'   => \$fid,#find id
	    'mid:s'  => \$mid,#match id
	    'idfile:s' => \$idfile, #text file with mutliple ids to pull out
	    "fdescr:s"  => \$fdescr,  # find descr
	    "rdm"  => \$rdm,  # randomly select a fasta sequence
	    "singletons" => \$singletons, # split file into singletons and non-singletons
	    "longest" => \$longest,
	    "freq:s"  => \$frequency, #frequency a sequence needs to be found at
	    "help"  => \$help,  # provides help with usage
           );
if (($help)&&!($help)||!($infasta)||!($outfile)){
 print "Usage : SelectSeq.pl <list of arguments>\n";
 print " -in <txt> - the input fasta file\n";
 print " -out <txt> - the name of your output fasta file\n";
 print " -fid <txt> - the id of the fasta file you want to match - exact match\n";
 print " -mid <txt> - the id of the fasta file you want to match - partial match\n";
 print " -idfile <txt> - a file with multiple ids you want to match pull out\n";
 print " -fdescr <txt> - the description of a fasta sequence you want to match - partial match\n";
 print " -rdm - randomly select a sequence from a fasta file\n";
 print " -singletons - split file into singletons and no singletons fasta\n";
 print " -longest - output the longest sequence from a amulitfasta file\n";
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
  my $id=$1 if $_=~/^(.+)\s*/;
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
 ## foreach my $listid (keys %list){
    if ($list{$id}){
     $outid->write_seq($seq);
     $hits++;
     delete $list{$id};
    }else{
      $out->write_seq($seq);
    }
 ## }
}
my @nomatches=keys %list;
print "There are $hits matches from the $idcnt ids\n";
print "No matches for @nomatches \n";

}elsif ($fid){  
print "Looking for $fid (exact id) in your fasta file $infasta...\n";
my $in  = Bio::SeqIO->new(-file => "$infasta" ,
                         -format => 'fasta');
my $out = Bio::SeqIO->new(-file => ">$outfile" , '-format' => 'fasta');
my $hits=0;
  while (my $seq = $in->next_seq()){
    my $id=$seq->id;
    if ($id eq $fid){
     $out->write_seq($seq);
     $hits++;
     }
  }
print "There are $hits matches\n";
}elsif ($mid){  
print "Looking for $mid match in the id of your fasta file $infasta...\n";
my $in  = Bio::SeqIO->new(-file => "$infasta" ,
                         -format => 'fasta');
my $out = Bio::SeqIO->new(-file => ">$outfile" , '-format' => 'fasta');
my $hits=0;
  while (my $seq = $in->next_seq()){
    my $id=$seq->id;
    if ($id =~/$mid/){
     $out->write_seq($seq);
     $hits++;
     }
  }
print "There are $hits matches\n";
}elsif($fdescr){
print "Looking for $fdescr in your fasta file $infasta...\n";
my $in  = Bio::SeqIO->new(-file => "$infasta" ,
                         -format => 'fasta');
my $out = Bio::SeqIO->new(-file => ">$outfile" , '-format' => 'fasta');
my $hits=0;
  while (my $seq = $in->next_seq()){
    my $desc=$seq->desc;
    if ($desc =~/$fdescr/){
     $out->write_seq($seq);
     $hits++;
    }
  }
print "there are $hits matches\n";
}elsif($rdm){
print "Randomly selecting a sequence...\n";
my (@seqs);
my $in  = Bio::SeqIO->new(-file => "$infasta" ,
                         -format => 'fasta');
my $out = Bio::SeqIO->new(-file => ">$outfile" , '-format' => 'fasta');
  while (my $seq = $in->next_seq()){
    push(@seqs,$seq);
  }
  my $index   = rand @seqs;
  my $rdmseq = $seqs[$index];
  print "$rdmseq\n";
  $out->write_seq($rdmseq);
}elsif($singletons){
my %seen;
my $in  = Bio::SeqIO->new(-file => "$infasta" ,
                         -format => 'fasta');
my $nosingletons = Bio::SeqIO->new(-file => ">$outfile.nosingletons.fa" , '-format' => 'fasta');
my $singletons = Bio::SeqIO->new(-file => ">$outfile.singletons.fa" , '-format' => 'fasta');

  while (my $seq = $in->next_seq()){
    my $id=$seq->id;
    my $seqstr=$seq->seq;
    $seen{$seqstr}{$id}++;
  }
  foreach my $seqstr (keys %seen){
    my $nbseqs=keys %{$seen{$seqstr}};
    print "$nbseqs $seqstr\n";
    if ($nbseqs>1){
      foreach my $id (keys %{$seen{$seqstr}}){
         my $seq_obj = Bio::Seq->new(-seq => $seqstr, 
                          -display_id => $id);#-desc => $nbseqs, 
        $nosingletons->write_seq($seq_obj);
      }
    }else{
      foreach my $id (keys %{$seen{$seqstr}}){
         my $seq_obj = Bio::Seq->new(-seq => $seqstr,
                          -display_id => $id);# -desc => $nbseqs, 
        $singletons->write_seq($seq_obj);
      }
    }
  }

}elsif($frequency){
my %seen;
my $in  = Bio::SeqIO->new(-file => "$infasta" ,
                         -format => 'fasta');
my $above = Bio::SeqIO->new(-file => ">$outfile.above$frequency.fa" , '-format' => 'fasta');
my $below = Bio::SeqIO->new(-file => ">$outfile.below$frequency.fa" , '-format' => 'fasta');
my $totalseq=0;
  while (my $seq = $in->next_seq()){
    my $id=$seq->id;
    my $seqstr=$seq->seq;
    $seen{$seqstr}{$id}++;
    $totalseq++;
  }
  foreach my $seqstr (keys %seen){
    my $nbseqs=keys %{$seen{$seqstr}};
    my $seqfreq=$nbseqs/$totalseq;
    print "$nbseqs $seqfreq $seqstr\n";
    if ($seqfreq>$frequency){
      foreach my $id (keys %{$seen{$seqstr}}){
         my $seq_obj = Bio::Seq->new(-seq => $seqstr, 
                          -display_id => $id);#-desc => $nbseqs, 
        $above->write_seq($seq_obj);
      }
    }elsif ($seqfreq<=$frequency){
      foreach my $id (keys %{$seen{$seqstr}}){
         my $seq_obj = Bio::Seq->new(-seq => $seqstr,
                          -display_id => $id);# -desc => $nbseqs, 
        $below->write_seq($seq_obj);
      }
    }
  }

}elsif($longest){
my $in  = Bio::SeqIO->new(-file => "$infasta" ,
                         -format => 'fasta');
my $out = Bio::SeqIO->new(-file => ">$outfile\_longest.fa" , '-format' => 'fasta');
my $length=0;
my ($longest_id,$longest_seq);
my $totalseq=0;
  while (my $seq = $in->next_seq()){
    my $id=$seq->id;
    my $seqstr=$seq->seq;
    if (length($seqstr)>$length){
      $longest_id=$id;
      $longest_seq=$seqstr;
      $length=length($seqstr);
    }
  }
  print "The longest sequence is $longest_id which is $length bp in length\n";
  my $seq_obj = Bio::Seq->new(-seq => $longest_seq, 
                              -display_id => $longest_id);#-desc => $nbseqs, 
  $out->write_seq($seq_obj);

}

