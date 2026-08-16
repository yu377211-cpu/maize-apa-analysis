#!/usr/bin/perl
use strict;
use warnings;

die "perl $0 <in.target> <in.pv> <out>" unless @ARGV==3;

my($file1,$file2,$out)=@ARGV;

#首先读取并存储第二个文件的所有位点信息
my%pos2pvalue;
my@positions;

open(my$fh2,'<',$file2)or die "无法打开文件 $file2: $!";
while(<$fh2>){
	chomp;
	my@fields=split/\s+/;
	if($fields[0]=~/^Chr(\d+)_(\d+)_[ACGT]_[ACGT]$/){
	my($chr,$pos)=($1,$2);
	$pos2pvalue{$chr}{$pos}=$fields[3];
	push @positions,[$chr,$pos];
	}
}
close$fh2;

open OUT, ">$out";
#处理第一个文件
open(my$fh1,'<',$file1)or die "无法打开文件 $file1: $!";
while(<$fh1>){
	chomp;
	my@fields=split/\s+/;
	my($chr,$pos)=($fields[0]=~/^chr(\d+)$/i?$1:(),$fields[2]);
	if($fields[5] eq "0"){
		my $nearest_pvalue= find_nearest_pvalue($chr,$pos,\%pos2pvalue,\@positions);
		print OUT "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\t$nearest_pvalue\t$fields[6]\n";
	}else{
		print OUT "$_\n";
	}
}
close$fh1;
close OUT;

#寻找最近位点的子程序
sub find_nearest_pvalue{
	my($target_chr,$target_pos,$pos2pvalue_ref,$positions_ref)=@_;
	my$min_dist=1e12;#初始化为一个大数
	my$best_pvalue="NA";

	foreach my $pos_ref(@$positions_ref){
		my($chr,$pos)=@$pos_ref;
		next unless $chr eq $target_chr;

		my $dist=abs($pos-$target_pos);
		if($dist<$min_dist){
			$min_dist=$dist;
			$best_pvalue=$pos2pvalue_ref -> {$chr}{$pos};
		}
	}
	return $best_pvalue;
}
