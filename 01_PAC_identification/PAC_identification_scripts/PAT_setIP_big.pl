use strict;
use Getopt::Long;
use Bio::Perl;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Seq;
use lib '/share/pub/xingsl/shilai/pipeline/tools/star-fusion/lib';
#use lib '/media/bmi/8E42CB1742CB0347/APAdb/STAR'; #U can set you library path
require ('/share/pub/xingsl/shilai/project/APA/PAC_3seq_map/PAC_identification_scripts/funclib.pl');

##########################################################################################
#  parameters & usage
##########################################################################################
my %opts=();
GetOptions(\%opts,"h","in=s","suf:s","flds:s","skip:i","chr=s");
my $USAGE=<< "_USAGE_";
Usage: 

PAT_setIP_big.pl -in F:/script_out_2/MCPA/1965_A_run321_TAGCTTGT_S69_L004_R1_001_f.reg1.A.PATs -skip 0 -suf "" -flds 0:1:2 -chr "F:/sys_rawdata/cecad_mm/mm10.fa"

-h=help
-in=input file (HAS columns: chr,strand,coord) or flds
-suf=suffix for output file (default is infile and XX.IP) if suf then <infile.real.suf and infile.IP.suf>
-flds=0:1:2(default)
-skip=0 default to skip N lines
-chr=chromosome file (can be large)
_USAGE_


#############################################################################
#  invoke the program                             
#############################################################################
die $USAGE if $opts{'h'}||!$opts{'in'}||!$opts{'chr'};
my $infile=$opts{'in'};
my $flds=$opts{'flds'};
my $nskip=0;
$nskip=$opts{'skip'} if $opts{'skip'} and $opts{'skip'}>0;

$flds='0:1:2' if (!$flds);
my @cols=split(/:/,$flds);
die "error flds=$flds" if scalar(@cols)!=3;

my $chrfile=$opts{'chr'};
die "chr=$chrfile not exists!" if !defined($chrfile) or isFileEmptyOrNotExist($chrfile);

my $oreal=$infile.'.real';
my $oIP=$infile.'.IP';
#if ($oreal=~/\.PA|PAT|PAs|PAS|PATs$/) {

my $suf=$opts{'suf'};
if ($suf) {
  $oreal=$oreal.".$suf";
  $oIP=$oIP.".$suf";
}
my ($cntIP,$cntReal)=(0,0);


# index (only the first time) and connect to the fasta file
my $db = Bio::DB::Fasta->new( $chrfile);

print "infile(in)=$infile\nchrfile(chr)=$chrfile\nflds=$flds\n";

open(IN,"<$infile") or die "cannot read in=$infile\n";

open(REAL ,">$oreal" ) or die "cannot write real=$oreal\n";
open( IP, ">$oIP" ) or die "cannot write IP=$oIP\n";
#open(OUT,">tets_PA2.txt");

while (my $line=<IN>) {
  chomp($line);
  if ($nskip) {
	$line=<IN>;
	$nskip--;
	next;
  }
  my @items=split(/\t/,$line);
  #print "@items\n";
  #my @items=split(/\s/,$line);
  my $chrName=$items[$cols[0]];
  my $strand=$items[$cols[1]];
  my $pos=$items[$cols[2]];
  #print OUT "$chrName\t$strand\t$pos\n";
  my ($start,$end)=($pos-9,$pos+10); #默认位点左右各有9,10nt (-10~-1[PA],1~10)
    if ($pos-10<0) { 
	  $start=1;
      $end=10+$pos;
	}
  my $subseq = uc($db->seq($chrName, $start => $end));
  if ($strand eq '-') {
	$subseq=reverseAndComplement($subseq,1,1);
  }
  if (isThisIP($subseq)) {
	$cntIP++;
	print IP "$line\n";
  } else {
	$cntReal++;
	print REAL "$line\n";
  }
}
close(REAL);
close(IP);
close(IN);
#close(OUT);
#if (!$suf) {
#  rename($oreal,$infile);
#}
print "Real\t$cntReal\nIP\t$cntIP\n";


###############################################
sub isThisIP {
	my($subseq)=shift;
	my $s=0;
	my @wseq=split(//,$subseq);
	return 1 if $subseq=~/A{6,}/i;
	my $e=$s+9;	
	my $cnt=0;
	for my $i ($s..$e) {
		$cnt++ if($wseq[$i] eq 'A');
	}
	return 1 if($cnt>=7);
	my $ee=$s+19;
	$e++;	
    while ($e<=$ee) {
		$cnt++ if($wseq[$e] eq 'A');
		$cnt-- if($wseq[$s] eq 'A');
		$s++;
		$e++;
        return 1 if($cnt==7);
    }
	return 0;
}


