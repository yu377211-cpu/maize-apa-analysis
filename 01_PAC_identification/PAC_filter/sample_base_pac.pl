use strict;
use warnings;

die "perl $0 <in1.sampleSta> <in2.QC> <out>" unless @ARGV==3;

my $lst=shift;
my $total=shift;
my $out=shift;

my %hash;

open IN1,"<$lst";
while(<IN1>){
	chomp;
	my @atm=split /\t/;
	$hash{$atm[0]}=$atm[1];
}
close IN1;

open IN2, "<$total";
open OUT, ">$out";
print OUT "SampleName\tCleanBases\tPAC\n";
while(<IN2>){
	chomp;
	my @tmp=split /\t/;
	if(exists $hash{$tmp[0]}){
		print OUT "$tmp[0]\t$tmp[2]\t$hash{$tmp[0]}\n";
	}
}
close IN2;
close OUT;
