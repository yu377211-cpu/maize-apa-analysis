use strict;
use warnings;
use Scalar::Util qw(looks_like_number);

die "perl $0 <in> <pairList> <out>" unless @ARGV==3;

my $in=shift;
my $lst=shift;
my $out=shift;

open IN,"<$in";
open OUT1,">$out.changeID.txt";

open IN2,"<$lst";
my %hash;
while(<IN2>){
	chomp;
	my @tmp=split /\t/;
	my $value=$tmp[1];
	$value =~ s/282set_//g;
	$hash{$tmp[0]}=$value;
}
close IN2;

while(<IN>){
	chomp;
	my @atm=split /\t/;
	print OUT1 "$atm[0]\t";
	foreach(1 .. scalar @atm-1){
		if(exists $hash{$atm[$_]}){
			$atm[$_] = $hash{$atm[$_]};
		}
		print OUT1 "$atm[$_]\t";
		
	}
	print OUT1 "\n";
}
close IN;
close OUT1;

