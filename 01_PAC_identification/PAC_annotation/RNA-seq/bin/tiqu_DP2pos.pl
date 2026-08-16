use strict;
use warnings;

die "perl $0 <../ERP011069_pair.APA.PAC.info_DP2> <../ERP011069.PAC.DP2.3UTR.pos_strand.bed> <out>" unless @ARGV==3;

my $lst=shift;
my $total=shift;
my $out=shift;

my %hash;

open IN1,"<$lst";
while(<IN1>){
	chomp;
	my @atm=split /\t/;
	$atm[0] =~ s/chr/Chr/g;
	my $string= "$atm[0]\t$atm[1]\t$atm[2]\t$atm[3]";
	$hash{$string}=$atm[6];
}
close IN1;

open IN2, "<$total";
open OUT, ">$out";
while(<IN2>){
	chomp;
	my @tmp=split /\t/;
	my $str= "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[5]";
	if(exists $hash{$str}){
		print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$hash{$str}\t$tmp[4]\t$tmp[5]\t$tmp[6]\n";
	}
}
close IN2;
close OUT;
