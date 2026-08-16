use strict;
use warnings;

die "perl $0 <in.id> <in.info> <out>" unless @ARGV==3;

my $id=shift;
my $in=shift;
my $out=shift;

my %hash;
open ID,"<$id";
while(<ID>){
	chomp;
	$hash{$_}=1;
}
close ID;

open IN,"<$in";
my $title=<IN>;
open OUT, ">$out";
print OUT "$title";
while(<IN>){
	chomp;
	my @atm=split /\t/;
	if(exists $hash{$atm[0]}){
		print OUT "$_\n";
	}
}
close IN;
close OUT;
