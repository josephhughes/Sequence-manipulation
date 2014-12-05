#usr/bin/perl

#
#  Converst SCARF to numeric
#  from seqanswer http://seqanswers.com/forums/showthread.php?t=2343
$start = time();
open (IN, "@ARGV[0]") || die ("could not open file\n");
open (OUT,">@ARGV[1]")|| die ("could not close the file\n");;
	while (<IN>) {	{	{	{
					@in=split(/:/);
					@same = "@in[0..5]\n";
		    			 } 	
				$diff =  "@in[6]\n";
				} 			
			
			chomp($diff);	
			chomp ($diff);		
			@ASCII_converted= unpack("C*",$diff);
			foreach $i(@ASCII_converted)	{
							$i-=64;	
							}
			}
                
		$same = join('',@same);
		$same =~ s/ /:/g;
		@Same = $same;
		chomp (@Same);
		@final=  "@Same:@ASCII_converted\n";  
		$final = join ('',@final);
		print OUT $final;
		   }	
$end = time();		
print "Hey!  It only took ", ($end-$start), " seconds!";	
close IN;
close OUT;
