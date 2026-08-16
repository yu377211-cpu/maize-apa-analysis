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
	$pos2=$pos+1;
	if($strand eq "+"){
		$pos_up=$pos-100;
		$pos_down=$pos+50;
		print OUT "$chr\t$pos_up\t$pos\tupstream\t+\n$chr\t$pos2\t$pos_down\tdownstream\t+\n";
	}elsif($strand eq "-"){
		$pos_up=$pos+100;
		$pos_down=$pos-50;
		print OUT "$chr\t$pos_down\t$pos\tdownstream\t-\n$chr\t$pos2\t$pos_up\tupstream\t-\n";
	}
}
close IN;
close OUT;
