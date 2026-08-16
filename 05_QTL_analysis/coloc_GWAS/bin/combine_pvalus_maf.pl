use strict;
use warnings;

die "perl $0 <in.freq> <in.pvalue> <in.site> <out>" unless @ARGV==4;

my $freq=shift;
my $pval=shift;
my $site=shift;
my $out=shift;

open IN,"<$freq";
open OUT,">$out";
my %hash;
while(<IN>){
	chomp;
	#1      chr1_1741    T    C     0.001595     1254
	my @argv=split /\t/;
	my $id=$argv[1];
	my $ref=$argv[2];
	my $alt=$argv[3];
	my $maf=$argv[4];
	my $string=$id."_".$ref."_".$alt;
	$hash{$string}=$maf;
}
close IN;

open PV,"<$pval";
my %hashp;
#Chr9_13110784_A_T	Zm00001d045111_T001|Zm00001d045111_T001|Chr9|+	546.924436364295	8.3236611534628e-158	9.12632594871518e-151
while(<PV>){
	chomp;
	my @array=split/\t/;
	my $pvalue=$array[3];
	$array[0]=~ s/Chr/chr/g;
	$hashp{$array[0]}=$pvalue;
}
close PV;

open IN2,"<$site";
print OUT "chr\tvariantID\tpos\tref\talt\tpvalue\tmaf\n";
while(<IN2>){
	chomp;
	my @tmp=split /\t/;
#Zm00001d021979_T007|Zm00001d021979_T007|Chr7|+	Chr7_166829063_C_G	0.987634625762289	L2	1	1
	$tmp[1]=~ s/Chr/chr/g;
	my @tmp2=split/_/,$tmp[1];
	my $chr=$tmp2[0];
	my $pos=$tmp2[1];
	my $ref=$tmp2[2];
	my $alt=$tmp2[3];
	my $id=$chr."_".$pos;
	my ($maf,$pvalue)=("NA","0");
	if(exists $hash{$tmp[1]}){
		$maf=$hash{$tmp[1]};
	}
	if(exists $hashp{$tmp[1]}){
		$pvalue=$hashp{$tmp[1]};
	}
	print OUT "$chr\t$id\t$pos\t$ref\t$alt\t$pvalue\t$maf\n";
}
close IN2;
close OUT;
