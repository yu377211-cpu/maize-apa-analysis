use strict;
use warnings;

die "perl $0 <in1> <in2> <out>" unless @ARGV==3;

my $target=shift;
my $all=shift;
my $out=shift;

open INS,"<$target";
my %hash;
<INS>;
while(<INS>){
        chomp;
	my $line=$_;
        $hash{$line}=1;
}
close INS;

open ALL, "<$all";
my $title=<ALL>;
open OUT,">$out";
print OUT "$title";
while(<ALL>){
	chomp;
	my @atm=split /\t/,$_;
	if(exists $hash{$atm[0]}){
		print OUT "$_\n";
	}
}
close ALL;
close OUT;
