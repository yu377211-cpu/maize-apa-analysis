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
	my $len=0;
	my $strand=$atm[2];
	if($strand eq "+"){
		$len= $atm[0] - $atm[1];
	}elsif($strand eq "-"){
		$len= $atm[1] - $atm[0];		
		}
	if($len <0){
		$len = 0;
	}
	print OUT "$len\n";	
}
close IN;
close OUT;
