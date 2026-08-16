use strict;
use warnings;

die "perl $0 <in> <out>" unless @ARGV==2;

my $in=shift;
my $out=shift;

open IN,"<$in";
open OUT,">$out";
while(<IN>){
	chomp;
	my @atm=split /\t/;
	my $pos=$atm[0];
	my $strand=$atm[4];
	my $ratio=0;
	my $length = abs($atm[3]-$atm[2])+1;
	if($strand eq "+"){
		$ratio = int( (($pos - $atm[2]) / $length) * 100);
	}elsif($strand eq "-"){
		$ratio = int( (($atm[3] - $pos) / $length) * 100);	
	}
	if($ratio < 0){$ratio = 0;}
	if($ratio > 100){$ratio = 100;}
	print OUT "$ratio\n";
}
close IN;
close OUT;
