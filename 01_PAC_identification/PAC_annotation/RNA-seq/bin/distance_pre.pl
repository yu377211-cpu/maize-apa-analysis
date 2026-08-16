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
	if($atm[10] eq "three_prime_UTR"){
		if($strand eq "+"){
			print OUT "$pos\t$atm[7]\t+\n";
		}elsif($strand eq "-"){
			print OUT "$pos\t$atm[8]\t-\n";	
		}
	}	
}
close IN;
close OUT;
