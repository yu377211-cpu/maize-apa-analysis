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
	my $pos=$atm[3];
	my $strand=$atm[11];
	my $ratio=0;
	my $length = abs($atm[8]-$atm[7])+1;
	if($atm[10] eq "three_prime_UTR"){
		if($strand eq "+"){
			$ratio = int( (($pos - $atm[7]) / $length) * 100);
		}elsif($strand eq "-"){
			$ratio = int( (($atm[8] - $pos) / $length) * 100);	
		}
		if($ratio < 0){$ratio = 0;}
		if($ratio > 100){$ratio = 100;}
		print OUT "$ratio\n";
	}	
}
close IN;
close OUT;
