use strict;
use warnings;

die "perl $0 <in1.key> <in2.infor> <out>" unless @ARGV==3;

my $lst=shift;
my $total=shift;
my $out=shift;

my %hash;

open IN1,"<$lst";
while(<IN1>){
	chomp;
	my @atm=split/\t/;
	$hash{$atm[0]}=$_;
}
close IN1;

open IN2, "<$total";
open OUT, ">$out";
while(<IN2>){
	chomp;
	my @tmp=split /\t/;
	if(exists $hash{$tmp[0]}){
		print OUT "$tmp[1]|$tmp[0]\t$tmp[2]\t$tmp[3]\t$hash{$tmp[0]}\n";
	}
}
close IN2;
close OUT;
