#!/usr/bin/perl

# script that provides the nucleotide pairwise distance,
# as a sliding window with repsect to one provided sequence
# input is:
# 1) nucleotide sequence alignment in FASTA
# 2) the name of the sequence to compare to
# 3) the size of the sliding window

# load module
use strict;
use Bio::SearchIO;
use Bio::SeqIO;
use Getopt::Long;

# get the information from the command line
#global vars
my ($aln, $window, $refname, $help, $out, $start, $end, $rm);

GetOptions( "aln=s"  => \$aln,     # alignment in fasta format
            "window=i" => \$window, # window size in amino acids    
            "start=i" => \$start, # start 350
            "end=i" => \$end, # end 9802 9500
            "refname=s" => \$refname, # name of the reference sequence     
		    "out=s" => \$out, # tsv file with pairwise distance per window
		    "remove" => \$rm, # remove gaps and Ns
		    "help" => \$help, );   # to get help
# get command-line argument, or die with a usage statement
if(!($help)&&!($aln)||!($window)||!($refname)||!($out))
 {
 print "Usage : SlidingWindowNuc.pl <list of arguments>, all arguments are necessary\n";
 print "Example: perl ~/Documents/Repo/Sequence-manipulation/SlidingWindowNuc.pl -aln Alignment.fa -window 30 -start 350 -end 9500 -refname 288_diversi -out test -remove\n";
 print " -aln <txt> - nucleotide alignment in fasta format\n";
 print " -window <int> - window size in aa\n";
 print " -start <int> - start of region of interest\n";
 print " -end <int> - end of region of interest\n";
 print " -refname <txt> -name of the reference sequence\n";
 print " -remove - if specified, removes any Ns, Xs and gaps\n";
 print " -out <txt> - name of output in tsv\n";
 print " -help          - Get more detailed help\n";
 exit();
 }
 
if($help)
 {
 print "Usage : SlidingWindow.pl <list of arguments>, all arguments are necessary\n";

print "---------------------------------------------------------------------\n";
 print " -aln <txt> - nucleotide alignment in fasta format\n";
 print " -window <int> - window size in aa\n";
 print " -start <int> - start of region of interest\n";
 print " -end <int> - end of region of interest\n";
 print " -refname <txt> -name of the reference sequence\n";
 print " -remove - if specified, removes any Ns, Xs and gaps\n";
 print " -out <txt> - name of output in tsv\n";
 print " -help          - Get more detailed help\n";
 exit();
 }

# create output file
open(OUT,">$out")||die "Can't open $out\n";

my $nucwindow=3*$window;
# read in the codon alignment
my $inseq = Bio::SeqIO->new('-file' => "<$aln",
                            '-format' => 'fasta');

if (!$start){
  $start=1;
}
if (!$end){
  $end=length($inseq->next_seq->seq());
}
print "Start $start\nEnd $end\n";

# translate the codon alignment
my (%translation,%nuc);
while (my $seq = $inseq->next_seq) {
  my $orf_str = $seq->subseq($start,$end);
  my $id=$seq->id;
  my $orf = Bio::Seq->new(-display_id => $id, -seq => $orf_str);
  my $prot=$orf->translate;
  print $id."\n".$prot->seq."\n";
  $nuc{$id}=$seq->seq;
  $translation{$id}=$prot->seq;
}

#print "$refname $translation{$refname}\n";
#my $tmp = Bio::Seq->new(-display_id => "Test", -seq => "tcagctccctctCTA");
#print "\n".$tmp->translate()->seq."\n";

# calculating the total mismatches
my (%dist);
for my $id1 (keys %translation){
  for my $id2 (keys %translation){
    if ($dist{"nuc"}{$id1}{$id2} || $dist{"nuc"}{$id2}{$id1} || $id1 eq $id2){
      #do nothing
    }else{ 
      my ($aa1,$aa2,$nuc1,$nuc2);
      if ($rm){
        ($aa1,$aa2)=removeX($translation{$id1},$translation{$id2});
        ($aa1,$aa2)=removeGap($aa1,$aa2);
      }else{
        ($aa1,$aa2)=($translation{$id1},$translation{$id2});
      }
      my $total_mm=mismatch_count($aa1,$aa2);
      my $aa_dist=$total_mm/length($translation{$id1});
      #print "Amino acid\t$id1\t$id2\t$aa_dist\n";
      if ($rm){
        ($nuc1,$nuc2)=removeX($nuc{$id1},$nuc{$id2});
        ($nuc1,$nuc2)=removeGap($nuc1,$nuc2);
      }else{
        ($nuc1,$nuc2)=($nuc{$id1},$nuc{$id2});
      }
      my $total_nuc_mm=mismatch_count($nuc1,$nuc2);
      my $nuc_dist=$total_nuc_mm/length($nuc{$id1});
      #print "Nucleotide\t$id1\t$id2\t$nuc_dist\n";
      $dist{"nuc"}{$id1}{$id2}=$nuc_dist;
      $dist{"aa"}{$id1}{$id2}=$aa_dist;
    }
  }
}

for my $id1 (keys %{$dist{"nuc"}}){
  for my $id2 (keys %{$dist{"nuc"}{$id1}}){
    print "Nuc\t$id1\t$id2\t". $dist{"nuc"}{$id1}{$id2}."\n";
  }
}
for my $id1 (keys %{$dist{"aa"}}){
  for my $id2 (keys %{$dist{"aa"}{$id1}}){
    print "AA\t$id1\t$id2\t". $dist{"aa"}{$id1}{$id2}."\n";
  }
}

# calculating the mismatches for each window
for my $id (keys %translation){
  if ($id ne $refname){
    my @aa=split(//,$translation{$id});
    for (my $i=0; $i<(scalar(@aa)-$window);$i++){
      my $windseq=substr($translation{$id},$i,$window);
      my $refseq=substr($translation{$refname},$i,$window);
      my ($string1,$string2);
      if ($rm){
        ($string1,$string2)=removeX($windseq,$refseq);
        ($string1,$string2)=removeGap($string1,$string2);
      }else{
        ($string1,$string2)=($windseq,$refseq);
      }
      my $mm=mismatch_count($string1,$string2);
      
      print OUT "AA\t$id\t$refname\t$i\t".($i*3+$start)."\t$mm\n";
    }

    my @nucs=split(//,$nuc{$id});
    for (my $j=0; $j<(scalar(@nucs)-$nucwindow); $j++){
      my $nucseq=substr($nuc{$id},$j,$nucwindow);
      my $refseq=substr($nuc{$refname},$j,$nucwindow);
      my ($string1,$string2);
      if ($rm){
        ($string1,$string2)=removeN($nucseq,$refseq);
        ($string1,$string2)=removeGap($string1,$string2);
      }else{
        ($string1,$string2)=($nucseq,$refseq);
      }
      my $mm=mismatch_count($string1,$string2);
      print OUT "Nuc\t$id\t$refname\t$j\t$j\t$mm\n";
    }
  }
}



#Quickly calculate hamming distance between two strings
#
#NOTE: Strings must be same length.
#      returns number of different characters.
#see  http://www.perlmonks.org/?node_id=500235
sub mismatch_count($$) { length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }

sub removeX {
  my ($result1,$result2);
  my ($string1,$string2)=@_;
  #print "1:$string1\n2:$string2\n";
  for(0 .. length($string1)) {
    my $char2 = substr($string2, $_, 1);
    my $char1 = substr($string1, $_, 1);
    if($char2=~/X/i || $char1=~/X/i) {
      #print "X found\n";
    } else {
        $result1 .= $char1;
        $result2 .= $char2;
    }
  }
  return($result1,$result2);
}


sub removeN {
  my ($result1,$result2);
  my ($string1,$string2)=@_;
  #print "1:$string1\n2:$string2\n";
  for(0 .. length($string1)) {
    my $char2 = substr($string2, $_, 1);
    my $char1 = substr($string1, $_, 1);
    if($char2=~/N/i || $char1=~/N/i) {
      #print "N found\n";
    } else {
        $result1 .= $char1;
        $result2 .= $char2;
    }
  }
  return($result1,$result2);
}

sub removeGap {
  my ($result1,$result2);
  my ($string1,$string2)=@_;
  #print "1:$string1\n2:$string2\n";
  for(0 .. length($string1)) {
    my $char2 = substr($string2, $_, 1);
    my $char1 = substr($string1, $_, 1);
    if($char2=~/-/i || $char1=~/-/i) {
      #print "N found\n";
    } else {
        $result1 .= $char1;
        $result2 .= $char2;
    }
  }
  return($result1,$result2);
}