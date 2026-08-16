use strict;
use warnings;

#----------------------------------------------------------------------
die "perl $0 <gene_leng_file>  <reads_count_files> <Ouput_file>" unless @ARGV==3;
#
##其中 gene_length_file没有表头，格式如下：
#
##Gene1  2312
##Gene2  3212
#
##其中 reads_count_files 有表头，格式如下：
#
##GeneID sample1 sample2 sample3
##Gene1    22    2123    122
##Gene2    21    112     1233
##Gene3    211   2112    31233
##-----------------------------code--------------------------------------


open FA,"$ARGV[0]"; #打开基因组length文件
my %len;
while(<FA>){
	chomp;
	my($id,$len)=split /\t/,$_,2;
	#print "$id\t$len\n";
	$len{$id}=$len/1000;
}
close FA;

open FA1,"$ARGV[1]"; #打开reads count文件
my $head=<FA1>;#去除行名,并获得其内容，后面打印用
my %hash;
while(<FA1>){
	chomp;
	my ($ID,$count)=split /\t/,$_,2;
	@{$hash{$ID}}=split /\t/,$count;
}
close FA1;

my $n=0; #用于计算$s
my $m=0; #用于计算$s
my $leng;
foreach $leng(keys %len){
	foreach (@{$hash{$leng}}){
		$_=($_)/($len{$leng}); #得到RPK
	#print "$leng\t$_\n";
		$n+=1;
	}
	$m+=1;
}

my $s=$n/$m; #两次循环使得列数变为$n,所以要除以$m;n是列，m是行。

#print "$s\n";

my %ALL=();

my $i;
for($i=0; $i<= $s-1;$i+=1){
	my $total=0;
	foreach $leng(keys %len){
		$total+=${$hash{$leng}}[$i]; #计算每一列的总reads count数（均一化基因长度以后的）
	}
	#print "$total\n";
	$ALL{$i}=$total; #每一列对应一个总的reads count，存到%ALL中
}

foreach $leng(keys %len){
	#print "$leng\n";
	my $index=0;
	#print "$ALL{$index}\n";
	foreach (@{$hash{$leng}}){
		#print "$_\n";
		$_=($_*1000000)/($ALL{$index}); #计算TPM值；
		#print "$leng\t$_\n";
		$index+=1;
	}
}

#foreach $leng(sort keys %hash){
	#$data=@{hash{$leng}};
#
#print "$leng\t",join ("\t",@{$data}),"\n"; #打印数据
#}


open OU ,">$ARGV[2]";
print OU "$head"; #打印行名，不需要加换行，因为读的时候以及读入换行符号了

foreach $leng(sort keys %len){ #打印数据，这里不能直接调用%hash,需要调用#%len；很奇怪
	print OU "$leng\t"; #打印基因名
	foreach (@{$hash{$leng}}){ #循环打印每个基因对应的几列数据
		printf OU "%.3f", "$_"; #设置小数位3位
		print OU "\t";
	}
	print OU "\n" #打印完成后，换行！
}
close OU;
