#!/usr/bin/perl -w

# A perlscript written by Joseph Hughes, University of Glasgow
##################################################################### 
# A perlscript written by Joseph Hughes, University of Glasgow
# Compile the constellation of mutations relative to a reference from an alignment
# input is an alignment and the identifier of the reference
# output is a comma, semicolon or tab delimiter dependent on the -sep value provided
# option to keep all sites or only keep the sites that have at least one mutation relative to the reference
# each output sequence has the same number of output values
# works for proteins or nucleotides
# provide the reference identifier as an input
# if gene label is provided then output mutation format is S:M23Y
#####################################################################
use strict;
use Getopt::Long; 
use Bio::SeqIO;                    
 
my ($in,$out,$help,$gene,$refid,$unknown);
my $sep="tab";
my $dropinvar = ''; 
&GetOptions(
	    'in:s'      => \$in,#fastafile
	    'out:s'   => \$out,#output text-tabe matrix
	    'gene:s'  => \$gene,#gene label
	    'sep:s'   => \$sep,
	    'ref:s'   => \$refid,
	    'dropinvar' => \$dropinvar,
	    'unknown:s' => \$unknown,
	    "help"  => \$help,  # provides help with usage
           );
if (($help)&&!($help)||!($in)||!($out)||!($unknown)){
 print "Usage : ConvertAln2Mat.pl <list of arguments>\n";
 print " -in <txt> - the input fasta file\n";
 print " -out <txt> - the name of your output text-tab delimited matrix\n";
 print " -gene <txt> - gene label\n";
 print " -sep <txt> - separator to use in the output can be one of comma OR semicolon or tab [default tab]\n";
 print " -dropinvar - use this argument if you want to drop invariable sites\n";
 print " -refid <txt> - the identifier of the reference\n";
 print " -unknown <txt> - unknown character usually X for protein and N for nucleotide alignments\n";
 print " -help        - Get this help\n";
 exit();
 }

my $inseq  = Bio::SeqIO->new(-file => "$in" ,
                         -format => 'fasta');
my $delim="\t";
if ($sep=~/comma/){
  $delim=",";
}
if ($sep=~/semicolon/){
  $delim=";";
}

open(OUT,">$out")||die "Can't open $out\n";
my $total_sites=0;
my %seqs;
while (my $seq = $inseq->next_seq()){
  my $id=$seq->id;
  my $seq_str=$seq->seq();
  $seqs{$id}=uc($seq_str);
  $total_sites=length($seq_str);
}
if (!$dropinvar){
	# just printing out a matrix
	print OUT "sequence_name\t";
	for (my $i=0; $i<$total_sites; $i++){
	  if ($gene){ 
	    if ($i==0){
	       print OUT $gene."position".($i+1);
	    }else{
		  print OUT $delim.$gene."position".($i+1);
		}
	  }else{
	    if ($i==0){
	      print OUT "position".($i+1);
	    }else{
	      print OUT $delim."position".($i+1);
	    }
	  }
	}
	print OUT "\n";
	if (!$refid){
		for my $id (keys %seqs){
		  my @bases=split(//,$seqs{$id}); #the input can also be a protein alignment so these could equally be amino acids
		  print OUT $id."\t".join($delim,@bases)."\n";
		}
	}elsif($refid){
	    my @refbases=split(//,$seqs{$refid});
		for my $id (keys %seqs){
		  if ($id!~$refid){
		    my @bases=split(//,$seqs{$id}); #the input can also be a protein alignment so these could equally be amino acids
		    print OUT $id."\t";
		    for(my $i=0; $i<@bases; $i++){
		      if ($gene){
		        if ($i==0){
		          print OUT $gene.":".$refbases[$i].($i+1).$bases[$i];
		        }else{
		          print OUT $delim.$gene.":".$refbases[$i].($i+1).$bases[$i];
		        }
		      }else{
		        #print $i."\n";
		        if ($i==0){
		          print OUT $refbases[$i].($i+1).$bases[$i];
		        }else{
		          print OUT $delim.$refbases[$i].($i+1).$bases[$i];
		        }
		      }
		    }
		    print OUT "\n";
		  }
		}
	}
}elsif($dropinvar){
  my %polymorphic;
  my %matrix;
  if ($refid){
    #keep track of all sites where there is a difference relative to the reference
    my @refbases=split(//,$seqs{$refid});
    for my $id (keys %seqs){
      if ($id!~$refid){
        my @bases=split(//,$seqs{$id});
        for(my $i=0; $i<scalar(@bases);$i++){
          if ($bases[$i] ne $refbases[$i] && $bases[$i] ne $unknown){
          print "$in $id not equal to $refid at $i $bases[$i] <=> $refbases[$i] unknown is >$unknown<\n";
            $polymorphic{$i}++;
          }
          if ($gene){
            $matrix{$id}{$i}=$gene.":".$refbases[$i].($i+1).$bases[$i];
          }elsif(!$gene){
            $matrix{$id}{$i}=$refbases[$i].($i+1).$bases[$i];
          }
        }
      }
    }
    print OUT "sequence_name\t";
    my @positions=@_;
    foreach my $position (sort {$a <=> $b} keys %polymorphic){
      push(@positions,($position+1));
    }  
    print OUT join($delim,@positions)."\n";
    foreach my $id (keys %matrix){
      my @mutations=@_;
      foreach my $position (sort {$a <=> $b} keys %polymorphic){
        push(@mutations,$matrix{$id}{$position});
      }
      print OUT "$id\t".join($delim,@mutations)."\n";
    }
    
    
  }else{
    # no reference identifier provided, need to keep a track of the different types 
    # of nucleotides or residues at each site and exclude later those where there is no variation
    for my $id (keys %seqs){
      my @bases=split(//,$seqs{$id});
      for (my $i=0; $i<scalar(@bases); $i++){
        if ($bases[$i] ne $unknown){
          $polymorphic{$i}{$bases[$i]}++;
        }
        if ($gene){
          $matrix{$id}{$i}=$gene.":".$bases[$i];
        }elsif(!$gene){
          $matrix{$id}{$i}=$bases[$i];
        }      
      }
    }
    # printing out the header
    print OUT "sequence_name\t";
    my @positions=@_;
    for my $position (sort {$a <=> $b} keys %polymorphic){
      my $nb_mutations=keys %{$polymorphic{$position}};
      if ($nb_mutations>1){
        push(@positions,$position+1);
      }
    }
    print OUT join($delim,@positions)."\n";
    foreach my $id (keys %matrix){
      my @mutations=@_;
      foreach my $site (@positions){
        push(@mutations,$matrix{$id}{$site-1});
      }
      print OUT "$id\t".join($delim,@mutations)."\n";
    }
  }
}

