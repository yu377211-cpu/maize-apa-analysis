use strict;
use warnings;

die "perl $0 <ERP011069.all.PAC.PATcount.txt> <in2.pos_strand> <out>" unless @ARGV==3;

my $lst=shift;
my $total=shift;
my $out=shift;

my %hash;

open IN1,"<$lst";
while(<IN1>){
	chomp;
	my @atm=split /\t/;
	my $string= "$atm[1]\t$atm[2]\t$atm[3]\t$atm[7]\t$atm[4]";
	$hash{$string}=$atm[0];
}
close IN1;

open IN2, "<$total";
open OUT, ">$out";
#my $title=<IN2>;
#print OUT "$title";
while(<IN2>){
	chomp;
	my @tmp=split /\t/;
	my $str= "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[5]";
	if(exists $hash{$str}){
		print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$hash{$str}\t$tmp[5]\t$tmp[6]\n";
	}
}
close IN2;
close OUT;
