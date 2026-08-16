use strict;
use warnings;

die "perl $0 <in.freq> <in.site> <out>" unless @ARGV==3;

my $freq=shift;
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
	my $string=$id."\t".$ref."\t".$alt;
	$hash{$string}=$maf;
}
close IN;

open IN2,"<$site";
print OUT "chr\tvariantID\tpos\tref\talt\tpvalue\tmaf\n";
while(<IN2>){
	chomp;
	my @tmp=split /\t/;
	my $chr=$tmp[2];
	my $pos=$tmp[3];
	my $ref=$tmp[5];
	my $alt=$tmp[6];
	my $pvalue=0;
	if ($tmp[4] =~ /^\d+(\.\d+)?$/) {
		$pvalue= 10**(-$tmp[4]);
	}
	$tmp[1] =~ s/\.s_/_/g;
	my $id=$tmp[1];
	my $string2=$id."\t".$ref."\t".$alt;
	if(exists $hash{$string2}){
		print OUT "$chr\t$id\t$pos\t$ref\t$alt\t$pvalue\t$hash{$string2}\n";
	}else{
		print OUT "$chr\t$id\t$pos\t$ref\t$alt\t$pvalue\tNA\n";
	}
}
close IN2;
close OUT;
