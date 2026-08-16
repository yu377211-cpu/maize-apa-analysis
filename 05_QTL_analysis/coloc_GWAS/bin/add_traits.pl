use strict;
use warnings;

die "perl $0 <in.traits> <in.nearest_match> <out>" unless @ARGV==3;

my $trait=shift;
my $match=shift;
my $out=shift;

open IN,"<$trait";
my %hash;
<IN>;
while(<IN>){
        chomp;
	#6	chr6_31730503	31730503	T	C	5.32917521912815e-06	NA	100grainweight
	my @argv=split /\t/;
	my $pvalue = sprintf("%.10e", $argv[5]);  # 强制科学计数法
	my $trait=$argv[7];
	my $string=$argv[1];
	$hash{$string}{$pvalue}=$trait;
}
close IN;

open IN2,"<$match";
open OUT,">$out";
my $title=<IN2>;
chomp $title;
print OUT "$title\ttrait.y\n";
while(<IN2>){
	chomp;
	my @atm=split /\t/;
	my $gid=$atm[9];
	my $pvalue = sprintf("%.10e", $atm[13]);
	if(exists $hash{$gid} && exists $hash{$gid}{$pvalue}){
		print OUT "$_\t$hash{$gid}{$pvalue}\n";
	}else{
		print OUT "$_\tNA\n";
	}
}
close IN2;
close OUT;
