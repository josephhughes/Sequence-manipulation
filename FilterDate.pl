# Perl script to filter a column based on a date
# Modified 2021-01-08 to allow for filtering before and after a certain date

use strict;
use Getopt::Long; 
use warnings;
use POSIX qw(strftime);
#use DateTime qw();
#my $date =DateTime->now->strftime('%Y%m%d');
use Time::ParseDate;
use Time::Piece;

my ($help,$in,$out,$before,$after,$colname);
&GetOptions(
	    'in:s'  => \$in, # text-tab file formatted with anonymisation
	    'out:s' => \$out, # filtered output file
	    'before:s' => \$before,
	    'after:s' => \$after,
	    'colname:s' => \$colname, #the column name that contains the date to fiter
        "help:s" => \$help,
           );

if ($help||!$in||!$out||!$colname){
  print "Usage: perl \n";
  print " -in <txt> - text-tab delimited file anonymised\n";
  print " -out <txt> -text-tab output\n";
  print " -after <txt> - keep everything after this date cut-off\n";
  print " -before <txt> - keep everything before this date cut-off\n";
  print " -colname <txt> - the column name that contains the dates\n";
  print " -help      - Get this help\n";
  exit();
}

open(OUT,">$out")||die "Can't open $out\n";
my ($after_cutoff,$before_cutoff);
if ($after){
  $after_cutoff=parsedate($after);
  print "Keeping everything after the following date: $after\n";
}
if ($before){
  $before_cutoff=parsedate($before);
  print "Keeping everything before the following date: $before\n";
}

open(IN,"<$in")||die "Can't open $in\n";
my $header=<IN>;
chomp($header);
print OUT "$header\n";
my @colnames=split(/\t/,$header);
my $colnumber;
for (my $i=0; $i<scalar(@colnames);$i++){
  if ($colname eq $colnames[$i]){
    $colnumber=$i;
  }
}
print "Using column $colnumber => $colname\n"; 

my $exc=0;
my $inc=0;
while(<IN>){
  chomp($_);
  my @col=split(/\t/,$_);
  my $tp;
  if ($col[$colnumber]=~/\d{4}-\d{2}-\d{2}/){# best format
    #print "Best format >$col[$colnumber]<\n";
    #2020-09-24
    #2020-09-22
    $tp = parsedate($col[$colnumber]);   
  }elsif ($col[$colnumber]=~/^\d{4}\-\d{2}$/){
    $tp = parsedate("$col[$colnumber]"."-01");
    #print "Assuming $col[$colnumber] is "."$col[$colnumber]"."-01"."\n";
  }elsif ($col[$colnumber]=~/^(\d{4}\-\d{2})-XX$/){
    $tp = parsedate("$1"."-01");
    #print "Assuming $col[$colnumber] is "."$1"."-01"."\n";
  }else{# try your best to parse
    #print "Trying my best to parse >$col[$colnumber]<\n";
    $tp = parsedate($col[$colnumber]);
  }
  #print "$col[0]\t$col[$colnumber]\t$tp\t$before\t$before_cutoff\n";
  if ($after && $tp>=$after_cutoff && $before && $tp<=$before_cutoff){
    my $date = localtime($tp)->strftime("%m/%d/%Y");
    #print $col[$colnumber],"\t",$tp,"\t",$date,"\n";
    $inc++;
    print OUT "$_\n";
  }elsif ($after && $tp>=$after_cutoff && !$before){
    $inc++;
    print OUT "$_\n";   
  }elsif (!$after && $before && $tp<=$before_cutoff){ 
    $inc++;
    print OUT "$_\n";       
  }else{
    #print "$col[$colnumber] is excluded\n";
    $exc++;
  }

}
print "$exc records are excluded\n";
print "$inc records are included\n";