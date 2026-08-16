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
	if($atm[6] eq "."){next;}
	
	my $pos=$atm[3];
	my $strand=$atm[5];
	my $ratio=0;
	my $length = abs($atm[8]-$atm[7])+1;
	if($strand eq "+"){
		$ratio = int( (($pos - $atm[7]) / $length) * 100);
	}elsif($strand eq "-"){
		$ratio = int( (($atm[8] - $pos) / $length) * 100);	
	}
	my $flag=0;
	if($ratio <=25){
		$flag="S";
	}elsif($ratio >25 && $ratio <75){
		$flag="M";
	}elsif($ratio >=75){
		$flag="L";
	}
	print OUT "$atm[4]\t$atm[9]\t$ratio\t$flag\n";
}
close IN;
close OUT;
