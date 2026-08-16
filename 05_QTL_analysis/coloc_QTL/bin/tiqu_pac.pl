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
	$hash{$_}=1;
}
close IN1;

open IN2, "<$total";
open OUT, ">$out";
#my $title=<IN2>;
#print OUT "$title";
while(<IN2>){
	chomp;
	my @tmp=split /\t/;
	if(exists $hash{$tmp[0]}){
		print OUT "$tmp[1]\n";
	}
}
close IN2;
close OUT;
