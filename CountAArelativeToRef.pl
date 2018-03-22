#!/usr/bin/perl -w
#
# use this to count the number of AA at each site relative to a given reference in the protein alignment

use strict;
use Getopt::Long; 
use Bio::SeqIO;

# Output: the different characters
#A	C	D	E	F	G	H	I	K	L	M	N	P	Q	R	S	T	V	W	Y	GAP
# a residue position per row
# additional column with total sequences, consensus and frequency

# additional input could be a list of sequences either from a particular lineage or from a particular host

#   perl CountAArelativeToRef.pl -in PB1.prot.reduced.fa -out PB1_residue_count.txt -r A_Wilson-Smith_1933 -p


my @aminoacids=qw/- A B C D E F G H I J K L M N P Q R S T V W X Y Z/;

my ($help,$infile,$out,$list,$ref,%list,%base,%gap,%ref,$site);
&GetOptions(
	    'in:s'      => \$infile,#multifasta file of the proteins
	    'out:s'    =>\$out,# a text-tab delimited 
	    'l:s'  =>\$list, #list of identifiers
	    'r:s' => \$ref, #the id of the reference you want to use. All columns with gaps in the reference will be stripped A_Puerto_Rico_8_1934
        's:i' => \$site, #site of interest
        'h:s' => \$help,
           );

if(!($help)&&!($ref))
 {
 print "Usage : CountAArelativeToRef.pl <list of arguments>, all arguments are necessary\n";
 print "Part of Sequence-manipulation by J. Hughes (joseph(dot)hughes(at)glasgow(dot)ac(dot)uk\n\n";
 print "The columns/sites with gaps in the reference are removed and the position provided by the user\n";
 print "is from the ungapped reference\n";
 print " -r <txt> - name of the reference sequence to use\n";
 print " -l <txt> list of identifiers to filter by, from a particular lineage or host\n";
 print " -in <txt> the protein alignment\n";
 print " -out <txt> a text tab-delimited output file\n";
 print " -s <int> site of interest\n";
 print " -help          - Get more detailed help\n";
 exit();
 }
 
if($help)
 {
 print "Usage : CountAArelativeToRef.pl <list of arguments>, all arguments are necessary\n";
 print "Part of Sequence-manipulation by J. Hughes (joseph(dot)hughes(at)glasgow(dot)ac(dot)uk\n\n";
 print "The columns/sites with gaps in the reference are removed and the position provided by the user\n";
 print "is from the ungapped reference\n";
 print " -r <txt> - name of the reference sequence to use\n";
 print " -l <txt> list of identifiers to filter by\n";
 print " -in <txt> the protein alignment\n";
 print " -out <txt> a text tab-delimited output file\n";
 print " -s <int> site of interest\n";
 print " -help          - Get more detailed help\n";

exit();
 }


           
my $in  = Bio::SeqIO->new(-file => "$infile" , '-format' => 'fasta');

open(OUT,">$out")||die "Can't open $out\n";
if ($list){
  print "A list of identifiers has been provided so the alignment will be filtered to keep these.\n";
  open(LIST,"<$list")||die "can't open $list\n";
  my $cnt_list=0;
  my $cnt_match=0;
  my $no_match=0;
  my $total_seq=0;
  while(<LIST>){
    chomp($_);
    if ($_=~/.+/){
      my $str=$_;
      $str =~ s/^\s+|\s+$//g;
      my @elements=split(/\t/,$str);
      $list{$elements[0]}++;
      $cnt_list++;
    }
  }
  if (!$list{$ref}){
    print "The reference id is not in the list of identifiers to filter\n";
    exit;
  }else{  
      
	  while ( my $seq_obj = $in->next_seq() ) {
		my $id=$seq_obj->display_id;
		my @splitseq=split(//,uc($seq_obj->seq));
		$total_seq++;
		if ($list{$id} && $id ne $ref){
		  print "$id is matched\n";
		  $cnt_match++;
		  my $position=1;
		  foreach my $aa (@splitseq){
			$base{$position}{$aa}++;
			$position++;
		  }
		}
		if ($id eq $ref){ #locating the gaps in the reference alignment
		  my $position=1;
		    foreach my $aa (@splitseq){
			$ref{$position}=$aa;
			  if ($aa eq "-"){
			  #print "Gap at $position in $ref\n";
			  $gap{$position}++;
			  }
			  $position++;
		    }
	    }elsif (!$list{$id} && $id ne $ref){
	      $no_match++;
	    }
	  }
  # now printing out and changing the coordinates according to the gaps in the reference
  print OUT join("\t","Position relative to $ref","reference AA",@aminoacids,"total sequences","Consensus AA","Frequency",),"\n";
  
  my $newposition=1;
  foreach my $position (sort { $a <=> $b} keys %base){
    my $consensus_aa;
    my $consensus_cnt=0;
    my $sum=0;
    if ($gap{$position}){
      $newposition=$newposition;
    }else{
      print OUT "$newposition\t$ref{$position}"; 
      foreach my $aa (@aminoacids){
        if ($base{$position}{$aa}){
          print OUT "\t$base{$position}{$aa}";
          $sum=$sum+$base{$position}{$aa};
          if ($base{$position}{$aa}>$consensus_cnt){
            $consensus_cnt=$base{$position}{$aa};
            $consensus_aa=$aa;
          }
        }else{
          print OUT "\t0";
        }
      }
      $newposition++;
      print OUT "\t$sum\t$consensus_aa\t".(sprintf "%.2f",$consensus_cnt/$sum)."\n";
    }
  }
 }
  print "$cnt_list values in the list provided\n";
  print "$cnt_match values matched in the fasta file\n";
  print "No matches for $no_match\n";
  print "Total sequences $total_seq\n";
  
}if (!$list){
  print "No list of identifiers has been provided\n";
  while ( my $seq_obj = $in->next_seq() ) {
    my $id=$seq_obj->display_id;
    my @splitseq=split(//,uc($seq_obj->seq)); 
    if ($id ne $ref){
      my $position=1;
      foreach my $aa (@splitseq){
        $base{$position}{$aa}++;
        $position++;
      }
    }
    if ($id eq $ref){ #locating the gaps in the reference alignment
      my $position=1;
      foreach my $aa (@splitseq){
        $ref{$position}=$aa;
        if ($aa eq "-"){
        #print "Gap at $position in $ref\n";
        $gap{$position}++;
        }
        $position++;
      }
    }
       
  }
  # now printing out and changing the coordinates according to the gaps in the reference
  print OUT join("\t","Position relative to $ref","reference AA",@aminoacids,"total sequences","Consensus AA","Frequency",),"\n";
  
  my $newposition=1;
  foreach my $position (sort { $a <=> $b} keys %base){
    my $consensus_aa;
    my $consensus_cnt=0;
    my $sum=0;
    if ($gap{$position}){
      $newposition=$newposition;
    }else{
      print OUT "$newposition\t$ref{$position}"; 
      foreach my $aa (@aminoacids){
        if ($base{$position}{$aa}){
          print OUT "\t$base{$position}{$aa}";
          $sum=$sum+$base{$position}{$aa};
          if ($base{$position}{$aa}>$consensus_cnt){
            $consensus_cnt=$base{$position}{$aa};
            $consensus_aa=$aa;
          }
        }else{
          print OUT "\t0";
        }
      }
      $newposition++;
      print OUT "\t$sum\t$consensus_aa\t".(sprintf "%.2f",$consensus_cnt/$sum)."\n";
    }
  }  
}

my $i;
open(MAT,"<$out")||die "Can't open gapless matrix\n";
my $head=<MAT>;
chomp($head);
my @labels=split(/\t/,$head);
while(<MAT>){
  chomp($_);
  if ($_=~/^$site\t/){
    my @vals=split(/\t/,$_);
    print "$labels[0]\n$labels[1]\t$vals[1]\n";
    print "$labels[28]\t$vals[28]\t$labels[29]\t$vals[29]\n";
    for ($i=2; $i<27; $i++){
      if ($vals[$i]>0){
        print "$labels[$i]\t".(sprintf "%.2f",$vals[$i]/$vals[27])."\n";
      }
      
    }
  }
}
