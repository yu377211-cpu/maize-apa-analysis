use strict;
use warnings;

die "perl $0 <ERP011069_pair.APA.PAC.info_DP2> <ERP011069_DP2.intersect.3UTR.pos> <number of sample> <out>" unless @ARGV==4;

my $lst=shift;
my $total=shift;
my $num=shift;
my $out=shift;

my %hash;

open IN1,"<$lst";
while(<IN1>){
	chomp;
	my @atm=split /\t/;
	my $chr= $atm[0];
	$chr =~ s/chr/Chr/g;
	my $str="$chr\t$atm[1]\t$atm[2]";
	$hash{$str}= $atm[7]/$num;
}
close IN1;

open IN2, "<$total";
open OUT, ">$out";
while(<IN2>){
	chomp;
	my @tmp=split /\t/;
	my $string="$tmp[0]\t$tmp[1]\t$tmp[2]";
	if(exists $hash{$string}){
		print OUT "$tmp[0]:$tmp[1]_$tmp[2]\t$hash{$string}\n";
	}
}
close IN2;
close OUT;
