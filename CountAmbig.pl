#!/usr/bin/perl -w
# A perlscript written by Joseph Hughes, University of Glasgow
# use this perl script to count ambiguities in a consensus sequence 


use strict;
use Getopt::Long; 
use Bio::SeqIO;                    
 
my ($infasta,$help,$detail,$ambig);
&GetOptions(
	    'in:s'      => \$infasta,#fastafile
	    'detail'  => \$detail,
	    'ambigonly'  => \$ambig,
	    "help"  => \$help,  # provides help with usage
           );
if (($help)&&!($help)||!($infasta)){
 print "Usage : CountAmbig.pll <list of arguments>\n";
 print " -in <txt> - the input fasta file\n";
 print " -detail <txt> - [optionally] provide of summary of each ambiguous IUPAC\n";
 print " -ambigonly <txt> - [optionally] provide of summary of just the ambiguous sites (not Ns)\n";
 print " -help        - Get this help\n";
 exit();
 }

my $in  = Bio::SeqIO->new(-file => "$infasta" ,
                         -format => 'fasta');

my @char =qw/R Y S W K M B D H V N/; 
my @onlyambigchar =qw/R Y S W K M B D H V/;                         
# sequence summary count Ns, R etc... 
my $sum=0;  
my %details;                     
while (my $seq = $in->next_seq()){
	my $id=$seq->id;
	my $str=$seq->seq;
  if ($ambig){
    for my $c (@onlyambigchar){
      my $count = () = $str =~ /\Q$c/gi;
      #print "$id\t count $count of $c\n";
      $sum=$sum+$count;
    }
    print "$id\t$sum\n";
  }
  else{
    for my $c (@char){
      my $count = () = $str =~ /\Q$c/gi;
      #print "$id\t count $count of $c\n";
      $sum=$sum+$count;
      $details{$c}=$count;
    }  
    if ($detail){
      foreach my $char (keys %details){
        print "$char\t$details{$char}\n";
      }
    }else{
      print "$id\t$sum\n";
    }
  }
}