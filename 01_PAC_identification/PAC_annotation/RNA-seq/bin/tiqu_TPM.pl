use strict;
use warnings;

die "perl $0 <ERP011069.PAC.CS.3UTR.pos_strand.bed> <ERP011069.all.PAC.hq.TPM> <out>" unless @ARGV==3;

my $lst=shift;
my $total=shift;
my $out=shift;

my %hash;

open IN1,"<$lst";
while(<IN1>){
	chomp;
	my @atm=split /\t/;
	$hash{$atm[4]}=1;
}
close IN1;

open IN2, "<$total";
open OUT, ">$out";
my $title=<IN2>;
print OUT "$title";
while(<IN2>){
	chomp;
	my @tmp=split /\t/;
	if(exists $hash{$tmp[0]}){
		print OUT "$_\n";
	}
}
close IN2;
close OUT;
