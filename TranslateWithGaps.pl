#!/usr/bin/perl
#
# TranslateWithGaps.pl  
# Script written by Joseph Hughes - University of Glasgow
#
#    Usage:  TranslateWithGaps -in aligned_nuc.fasta  -out aligned_prot.fasta -gc flu 
#
#        pep.aln:    protein alignment either in CLUSTAL or FASTA format
#
#        nuc.fasta:  DNA sequences (single multi-fasta or separated files)
#
#        Options:  -gc genetic code (flu|universal(default)|degen|vmitocondria)
#

use strict;
use Getopt::Long; 
use Bio::SeqIO;


my ($in,$out,$gencode,$help,$r);
&GetOptions(
	    'in:s'      => \$in,#a multifasta aaligned nucletoide sequence file 
	    'out:s'      => \$out,#output file aligned fasta protein
	    'gc:s'      => \$gencode,#genetic code
	    "help"  => \$help,  # provides help with usage
           );

if (($help)&&!($help)||!($in)){
 print "Usage : TranslateWithGaps.pl   <list of arguments>\n";
 print " -in <txt> - aligned nucleotide fasta files\n";
 print " -out <txt> - output aligned protein fasta files\n";
 print " -gc <txt> - gencode to use one of flu|universal(default)|degen|vmitocondria\n";
 print " -help        - Get this help\n";
 exit();
 }


my %c2p;

if ($gencode eq "universal") {
    #-----------#
    # universal
    #-----------#
    $c2p{$_} = "L" for qw(CTT CTC CTA CTG TTA TTG);
    $c2p{$_} = "R" for qw(CGT CGC CGA CGG AGA AGG);
    $c2p{$_} = "S" for qw(TCT TCC TCA TCG AGT AGC);
    $c2p{$_} = "A" for qw(GCT GCC GCA GCG);
    $c2p{$_} = "G" for qw(GGT GGC GGA GGG);
    $c2p{$_} = "P" for qw(CCT CCC CCA CCG);
    $c2p{$_} = "T" for qw(ACT ACC ACA ACG);
    $c2p{$_} = "V" for qw(GTT GTC GTA GTG);
    $c2p{$_} = "I" for qw(ATT ATC ATA);
    $c2p{$_} = "-" for qw(---);
    $c2p{$_} = "*" for qw(TAA TGA TAG);
    $c2p{$_} = "C" for qw(TGT TGC);
    $c2p{$_} = "D" for qw(GAT GAC);
    $c2p{$_} = "E" for qw(GAA GAG);
    $c2p{$_} = "F" for qw(TTT TTC);
    $c2p{$_} = "H" for qw(CAT CAC);
    $c2p{$_} = "K" for qw(AAA AAG);
    $c2p{$_} = "N" for qw(AAT AAC);
    $c2p{$_} = "Q" for qw(CAA CAG);
    $c2p{$_} = "Y" for qw(TAT TAC);
    $c2p{$_} = "M" for qw(ATG);
    $c2p{$_} = "W" for qw(TGG);
    $c2p{$_} = "X" for qw(NNN);

} elsif ($gencode eq "vmitochondria") {
    #--------------------------#
    # vertebrate mitochondrial
    #--------------------------#
    print "Codon table not ready yet\n";
#     %p2c = (
#         "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
#         "R" => "(CG.)",
#         "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
#         "A" => "(GC.)",
#         "G" => "(GG.)",
#         "P" => "(CC.)",
#         "T" => "(AC.)",
#         "V" => "(G(U|T).)",
#         "I" => "(A(U|T)(U|T|C|Y))",
#         "_" => "(((U|T)A(A|G|R))|(AG(A|G|R)))",
#         "*" => "(((U|T)A(A|G|R))|(AG(A|G|R)))",
#         "C" => "((U|T)G(U|T|C|Y))",
#         "D" => "(GA(U|T|C|Y))",
#         "E" => "(GA(A|G|R))",
#         "F" => "((U|T)(U|T)(U|T|C|Y))",
#         "H" => "(CA(U|T|C|Y))",
#         "K" => "(AA(A|G|R))",
#         "N" => "(AA(U|T|C|Y))",
#         "Q" => "(CA(A|G|R))",
#         "Y" => "((U|T)A(U|T|C|Y))",
#         "M" => "(A(U|T)(A|G|R))",
#         "W" => "((U|T)G(A|G|R))",
#         "X" => "...",
#     );
} elsif ($gencode eq "flu") {
    #--------------------------#
    # flu
# Ambiguous Amino Acids	3-Letter	1-Letter
# Asparagine or aspartic acid	Asx	B
# Glutamine or glutamic acid	Glx	Z
# Leucine or Isoleucine	Xle	J
# Unspecified or unknown amino acid	Xaa	X
    #--------------------------#
   print "Codon table not ready yet\n";

#    %p2c = (
#         "Z" => "(SAA)",
#         "J" => "(MT(T|A))",
#         "B" => "((AA(C|T))|(GA(C|T))|(RAT)|(RAC))",
#         "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R))|(YT(A|G)))",
#         "R" => "((CG.)|(AG(A|G|R))|(MG(A|G)))",
#         "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
#         "A" => "(GC.)",
#         "G" => "(GG.)",
#         "P" => "(CC.)",
#         "T" => "(AC.)",
#         "V" => "(G(U|T).)",
#         "I" => "(A(U|T)(U|T|C|Y|A|W))",
#         "_" => "(((U|T)A(A|G|R))|((T|U)GA))",
#         "*" => "(((U|T)A(A|G|R))|((T|U)GA))",
#         "C" => "((U|T)G(U|T|C|Y))",
#         "D" => "(GA(U|T|C|Y))",
#         "E" => "(GA(A|G|R))",
#         "F" => "((U|T)(U|T)(U|T|C|Y))",
#         "H" => "(CA(U|T|C|Y))",
#         "K" => "(AA(A|G|R))",
#         "N" => "(AA(U|T|C|Y))",
#         "Q" => "(CA(A|G|R))",
#         "Y" => "((U|T)A(U|T|C|Y))",
#         "M" => "(A(U|T)G)",
#         "W" => "((U|T)GG)",
#         "X" => "...",
#    );
}elsif ($gencode eq "degen") {
    #--------------------------#
    # degenerate table
# Ambiguous Amino Acids	3-Letter	1-Letter
# Asparagine or aspartic acid	Asx	B
# Glutamine or glutamic acid	Glx	Z
# Leucine or Isoleucine	Xle	J
# Unspecified or unknown amino acid	Xaa	X
    #--------------------------#

    $c2p{$_} = "L" for qw(CTT CTC CTA CTG TTA TTG TTR YTR CTN CTY CTR CTK CTM CTS CTW CTH CTD);
    $c2p{$_} = "R" for qw(CGT CGC CGA CGG AGA AGG CGN AGR CGY CGR MGR CGK CGM CGS CGW CGH CGD);
    $c2p{$_} = "S" for qw(TCT TCC TCA TCG AGT AGC TCN AGY TCY TCR TCK TCM TCS TCW TCH TCD);
    $c2p{$_} = "A" for qw(GCT GCC GCA GCG GCN GCY GCR GCK GCM GCS GCW GCH GCD);
    $c2p{$_} = "G" for qw(GGT GGC GGA GGG GGN GGY GGR GGK GGM GGS GGW GGH GGD);
    $c2p{$_} = "P" for qw(CCT CCC CCA CCG CCN CCY CCR CCK CCM CCS CCW CCH CCD);
    $c2p{$_} = "T" for qw(ACT ACC ACA ACG ACN ACY ACR ACK ACM ACS ACW ACH ACD);
    $c2p{$_} = "V" for qw(GTT GTC GTA GTG GTN GTY GTR GTK GTM GTS GTW GTH GTC);
    $c2p{$_} = "I" for qw(ATT ATC ATA ATH ATY ATW ATM);
    $c2p{$_} = "-" for qw(---);
    $c2p{$_} = "*" for qw(TAA TGA TAG TRA TAR);
    $c2p{$_} = "C" for qw(TGT TGC TGY);
    $c2p{$_} = "D" for qw(GAT GAC GAY);
    $c2p{$_} = "E" for qw(GAA GAG GAR);
    $c2p{$_} = "F" for qw(TTT TTC TTY);
    $c2p{$_} = "H" for qw(CAT CAC CAY);
    $c2p{$_} = "K" for qw(AAA AAG AAR);
    $c2p{$_} = "N" for qw(AAT AAC AAY);
    $c2p{$_} = "Q" for qw(CAA CAG CAR);
    $c2p{$_} = "Y" for qw(TAT TAC TAY);
    $c2p{$_} = "M" for qw(ATG);
    $c2p{$_} = "W" for qw(TGG);
    $c2p{$_} = "X" for qw(NNN NTA AT- NTN TNN ANN NNT CNA -TT -TC A-- CNT TAN AAN T-- GNT GNG NNA NGA TTN T-A ANT NNA NAT);
    $c2p{$_} = "X" for qw(GKT TNT NCT TT- NAA NYC TA- GT- NCN AKA NCA TGN NNC NCC CNN TNC NGT TNA -AT -AA NTC);
    $c2p{$_} = "X" for qw(NNG GAN GCK KAC NGC ATN GNN CAN NTT GNA G-- NAT NTT GYG ANA ANG YAG --T NTG -TA);
    $c2p{$_} = "X" for qw(NAN NAG AWC YTT YTA TYA NGG AYA GYT AGN --A ANC -CN -CT NGN TYT GGS GNC AC-);
    $c2p{$_} = "X" for qw(TNG TYT NAC KCA ATK ARG KGK AAK KTT TKT ARG -TG -AC -GA CNG CYA KGT KTT);
    $c2p{$_} = "X" for qw(KCC NCG NCG YCA CWT CA- AMA RAA KTC C-- GYA AYT CAK RCA YAT KTA RGA CYT RTT);
    $c2p{$_} = "X" for qw(KGC YGT KAT CKG CRA AYC TYC TTK AKT YAC SCT ARA -GT);
    $c2p{$_} = "X" for qw(KCT STG RCT WGC SAA YTK YTC YCT RAR AST YCG GRA -MR AAW MAA YYT WRG WTT KKK WAR RAW WWT KRW TWW WAA);
    $c2p{$_} = "X" for qw(AWW WWW TYY SMM MMS WAW TTW RKK SST GAW YYY AAM GRT TKA AMC WCT RTG GAS GKW TTM AWT TWT AMT WGG KGG RTA RAT GRC THT KTG KCG);
    $c2p{$_} = "X" for qw(ASA CRT KNN GMT SAT GAM CMT CAW RAC TYG G-T AKC GAK YAA MTT);
    $c2p{$_} = "X" for qw(GG- NTS RTC TRT WAT -GG GWG RGT CKT -NN C-N WCA WCN GKA CWA AGW SCA KGA TST);
    $c2p{$_} = "X" for qw(MGT TWA ATD T-T CMA CKA GMA STC MAC CKK MTA WGT CAM WGA MCA);
    $c2p{$_} = "X" for qw(MCT GKG AGK GWA TTS GKC CWG TMG WCG GRY RCC ART);
    $c2p{$_} = "X" for qw(SAG GWT GSG NVT RRT YAW TAK WCC WTG ASG ASC WTC MKT CC- MGA MAG AMG -NT CNC GC- -AG AKG SCC -AN TMT SCG GYC GRG GWC SNN GRN AG- CWC CYC CRG AWG ARC AYG ATS TDT TGR KAG CVT RKT RAG GRK --N GRW TYN YTG YCC TWC AA- MTC NMT -CC NKT SGA AWA SGT SAC GDT G-G WMA RNN KAA RGG WAC YCW RGC A-N WAG NNY NST GST GSA TSG CWR TGK AMY ATR TGW CSA CST WWA ASM WTA MAT TKN TGS TKG YGC YGG TKC AAS STA STT AGM TC- CAS);
    $c2p{$_} = "X" for qw(--C -CA GMC TRN VAT CKC MCG --G TWS CYG GMG SGC NCY NNM G-N CMG GA- NGW MTG BAT NNK MTN TAW HAT KMA NCW CSG T-G G-A A-T NAY WCW CNY TRC C-G TWG);
    $c2p{$_} = "X" for qw(MGC GSC NGY YGA);
 }   
#---------------------#
#  Get nuc sequences
#---------------------#

my $seqio_obj = Bio::SeqIO->new(-file => "$in", -format => 'fasta' );
open(OUT,">$out")||die "Can't open $out\n";
my %dodgy;
while (my $seq = $seqio_obj->next_seq()){
  my $id=$seq->id;
  my $seq_str=$seq->seq();
  my $nuclen = length($seq_str);
  my $mod = $nuclen % 3;
  if ($mod == 0){ 
    my @codons = ( $seq_str =~ m/.../g );
    print OUT ">$id\n";
    my $aapos=1;
    foreach my $codon (@codons){
      if ($c2p{uc($codon)}){
        print OUT $c2p{uc($codon)};
        $aapos++;
      }else{
        print OUT "X";
        $aapos++;

        #print "$codon at aa position $aapos is not in codon table $gencode\n";
        #exit;
        $dodgy{$codon}++;
      }
    }
    print OUT "\n";
  }else{
    print "The sequence $id is not a multiple of 3. Not a codon alignment\n";
    exit;
  }

}

print "Dodgy codons\n";
for my $c (keys %dodgy){
  print "$c ";
}           
print "\n";           
