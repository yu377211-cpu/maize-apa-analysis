use strict;
use warnings;

die "perl $0 <in> <out>" unless @ARGV==2;

my $in=shift;
my $out=shift;

open IN,"<$in";
open OUT,">$out";
<IN>;
while(<IN>){
	chomp;
	my @argv=split /\t/;
	my @atm=split /_/,$argv[1];
	my $chr=$atm[0];
	my $pos=$atm[1];
	my $start=$pos-500000;
	if($start < 0){$start=0;}
	my $end=$pos+500000;
	print OUT "$chr\t$start\t$end\n";
}
close IN;
close OUT;
