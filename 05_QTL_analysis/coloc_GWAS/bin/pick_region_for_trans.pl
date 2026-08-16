use strict;
use warnings;

die "perl $0 <in> <out>" unless @ARGV==2;

my $in=shift;
my $out=shift;

open IN,"<$in";
open OUT,">$out";
<IN>;
while(<IN>){
	chomp;
	#Chr1_228779089_TA_T	Zm00001d015115_T001|Zm00001d015115_T001|Chr5|-	Inf	2.2250738585072e-308	6.25925405898699e-300
	my @argv=split /\t/;
	my @atm=split /_/,$argv[0];
	my $chr=$atm[0];
	my $pos=$atm[1];
	my $start=$pos-500000;
	if($start < 0){$start=0;}
	my $end=$pos+500000;
	print OUT "$chr\t$start\t$end\n";
}
close IN;
close OUT;
