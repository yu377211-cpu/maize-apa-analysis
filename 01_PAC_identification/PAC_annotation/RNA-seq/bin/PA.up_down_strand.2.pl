use strict;
use warnings;

die "perl $0 <in.bed> <out>" unless @ARGV==2;

my $in=shift;
my $out=shift;

open IN,"<$in";
open OUT,">$out";

while(<IN>){
	chomp;
	my @atm=split /\t/;
	my $pos=$atm[3];
	my $chr=$atm[0];
	my $strand=$atm[5];
	my ($pos_up,$pos_down,$pos2);
	if($strand eq "+"){
		$pos_up=$pos-35;
		$pos_down=$pos-10;
		print OUT "$chr\t$pos_up\t$pos_down\tupstream\t+\n";
	}elsif($strand eq "-"){
		$pos_up=$pos+10;
		$pos_down=$pos+35;
		print OUT "$chr\t$pos_up\t$pos_down\tdownstream\t-\n";
	}
}
close IN;
close OUT;
