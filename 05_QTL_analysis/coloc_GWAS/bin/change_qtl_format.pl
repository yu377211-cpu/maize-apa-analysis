use strict;
use warnings;

die "perl $0 <in.qtl> <out.format>" unless @ARGV==2;

my $qtl=shift;
my $out=shift;

open IN,"<$qtl";
open OUT,">$out";
<IN>;
print OUT "variantID\tchr\tpos\tref\talt\tpvalue\ttargetGene\n";
while(<IN>){
	chomp;
	#Chr1_228779089_TA_T	Zm00001d015115_T001|Zm00001d015115_T001|Chr5|-	Inf	2.2250738585072e-308	6.25925405898699e-300
	my @argv=split /\t/;
	my @atm=split /_/,$argv[0];
	my $chr=$atm[0];
	my $pos=$atm[1];
	my $ref=$atm[2];
	my $alt=$atm[3];
	my $id=$atm[0]."_".$atm[1];
	my @tmp=split /\|/,$argv[1];
	my $gene=$tmp[0];
	my $pvalue=$argv[3];
	#variantID	chr	pos	ref	alt	pvalue	gene
	print OUT "$id\t$chr\t$pos\t$ref\t$alt\t$pvalue\t$gene\n";
}
close IN;
close OUT;
