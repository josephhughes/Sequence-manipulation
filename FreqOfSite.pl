#!/usr/bin/perl

# load module
use strict;
use Getopt::Long;
use Bio::SeqIO;

# Use this script to get the frequency of a particular nucleotide at a position
# from an epidemiological set for example.
# Users provides an alignment that contains a reference
# The alignment may be gapped, so you provide the name of your reference you are comparing to.
# The user needs to know how many gaps are in the reference

#global vars
my ( $help, $ref, $set, $site);

GetOptions( "ref=s" => \$ref, # name of the reference to use
		    "set=s" => \$set, #the name of the aligned fasta file which contains the alignment 
		    "site=i" => \$site, #
		    "help" => \$help, );   # to get help

if(!($help)&&!($ref))
 {
 print "Usage : FreqofSite.pl <list of arguments>, all arguments are necessary\n";
 print "Part of Sequence-manipulation by J. Hughes (joseph(dot)hughes(at)glasgow(dot)ac(dot)uk\n\n";
 print "The columns/sites with gaps in the reference are removed and the position provided by the user\n";
 print "is from the ungapped reference\n";
 print " -ref <txt> - name of the reference sequence to use\n";
 print " -set <txt> an alignment of fasta sequences\n";
 print " -site <txt> the site wehere there is a mutation of interest\n";
 print " -help          - Get more detailed help\n";
 exit();
 }
 
if($help)
 {
print "To help with getting the frequency of a mutations at a site\n\n";

print "---------------------------------------------------------------------\n";
print " -ref <txt> - name of the reference sequence to use\n";
print " -set <txt> an alignment of fasta sequences\n";
print " -site <txt> the site wehere there is a mutation of interest\n";
print " -help    - Gets this help information\n";

exit();
 }
 
# matrix used to store the proportion of ATCG at a site 
my (%matrix,%refungapped,@ref);
my $inseq  = Bio::SeqIO->new(-file => "$set" ,
                         -format => 'fasta');
print "$ref\n";
my $gapcnt=0;
while ( my $seq_obj = $inseq->next_seq() ) {
  my $id=$seq_obj->display_id;
  #print "$id\n";
  if ($id=~/$ref/){
    #print "$id $ref\n";
    @ref=split(//,uc($seq_obj->seq));
    for (my $i=0;$i<scalar(@ref); $i++){
      #print "$ref[$i]\n";
      if ($ref[$i]=~/\-/){
        #print "gap\n";
        $gapcnt++;
        my $ungappedsite=$i+1-$gapcnt;
        my $alignsite=$i+1;
        #print "$ungappedsite $alignsite\n";
        $refungapped{$ungappedsite}=$i+1;
      }else{
        my $ungappedsite=$i+1-$gapcnt;
        $refungapped{$ungappedsite}=$i+1;
      }
    }
  }
  my @nucs=split(//,uc($seq_obj->seq));
  for (my $j=0;$j<scalar(@nucs); $j++){
    $matrix{$j+1}{$nucs[$j]}++;
  }
}
print "Site $site Reference $ref[$refungapped{$site}]\n";
print "Site with gaps $refungapped{$site}\n";
print "As: $matrix{$refungapped{$site}}{A}\n";
print "Cs: $matrix{$refungapped{$site}}{C}\n";
print "Ts: $matrix{$refungapped{$site}}{T}\n";
print "Gs: $matrix{$refungapped{$site}}{G}\n";
   

