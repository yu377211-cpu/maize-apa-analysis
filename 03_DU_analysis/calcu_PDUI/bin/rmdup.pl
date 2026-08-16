#!/usr/bin/perl
use strict;
use warnings;

# 检查命令行参数
die "Usage: $0 <input_file>\n" unless @ARGV == 1;
my $input_file = $ARGV[0];

# 打开文件
open my $fh, '<', $input_file or die "Cannot open file '$input_file': $!";

# 用于存储行的哈希表
my %seen;

while (my $line = <$fh>) {
	chomp $line;
	my @fields = split "\t", $line;

    # 检查是否有足够的列
	die "Line does not have enough fields: $line" unless @fields > 1;

     #取出第一列并分割
	my @first_column_parts = split "_", $fields[0];
	my $tmp0 = $first_column_parts[0];
    
        # 构建一个用于比较的键
	my $key = join "\t", $tmp0, @fields[1..$#fields];
    
        #检查是否已经见过这个键
	if (exists $seen{$key}) {
        	next;  #如果已经见过，跳过这一行
	} else {
        	$seen{$key} = 1;  #标记这个键为已见
		print "$line\n";  #打印这一行
        }
}

close $fh;
