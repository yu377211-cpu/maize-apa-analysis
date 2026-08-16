#!/usr/bin/env perl
use strict;
use warnings;

die "Usage: perl $0 <UTR.bed> <intersect.bed> <out>\n" unless @ARGV == 3;

my ($utr_bed, $intersect_file, $out_file) = @ARGV;

# ------------------------------------------------------------
# 1. 读取 3'UTR BED，构建每个转录本的串联外显子信息
# ------------------------------------------------------------
my %transcripts;
open(my $fh_utr, '<', $utr_bed) or die "Cannot open $utr_bed: $!";
while (<$fh_utr>) {
    chomp;
    my @f = split /\t/;
    next if @f < 6;
    my ($chr, $start, $end, $tid, $type, $strand) = @f[0..5];
    push @{$transcripts{$tid}->{raw}}, { start => $start, end => $end };
    $transcripts{$tid}->{strand} = $strand;
}
close $fh_utr;

# 对每个转录本：按 5'→3' 排序外显子，计算每个外显子的长度和前缀累加长度
foreach my $tid (keys %transcripts) {
    my $strand = $transcripts{$tid}->{strand};
    my @exons = @{$transcripts{$tid}->{raw}};
    
    # 排序：正链按 start 升序，负链按 start 降序（因为 5' 端在右侧）
    if ($strand eq '+') {
        @exons = sort { $a->{start} <=> $b->{start} } @exons;
    } else {
        @exons = sort { $b->{start} <=> $a->{start} } @exons;
    }
    
    my @exon_info;
    my $total_len = 0;
    foreach my $exon (@exons) {
        my $len = abs($exon->{end} - $exon->{start}) + 1;
        push @exon_info, {
            start  => $exon->{start},
            end    => $exon->{end},
            len    => $len,
            prefix => $total_len   # 当前外显子之前的累积长度
        };
        $total_len += $len;
    }
    $transcripts{$tid}->{exons}     = \@exon_info;
    $transcripts{$tid}->{total_len} = $total_len;
}

# ------------------------------------------------------------
# 2. 处理 intersect 文件，直接使用其提供的 UTR 外显子坐标
# ------------------------------------------------------------
open(my $fh_in, '<', $intersect_file) or die "Cannot open $intersect_file: $!";
open(my $fh_out, '>', $out_file) or die "Cannot write $out_file: $!";
print $fh_out "PACid\tTranscript\tRatio\tFlag\n"; 

while (<$fh_in>) {
    chomp;
    my @f = split /\t/;
    next if @f < 6;
    
    # 跳过没有 UTR 重叠的行（第 7 列为 "."）
    next if $f[6] eq '.';
    
    # 提取关键字段（索引从 0 开始）
    my $pac_id     = $f[4];          # PAC 名称
    my $pos        = $f[3];          # PAC 位点（单碱基）
    my $tid        = $f[9];          # 转录本 ID
    my $exon_start = $f[7];          # UTR 外显子起始（来自 intersect 的 -wb 输出）
    my $exon_end   = $f[8];          # UTR 外显子终止
    my $utr_strand = $f[11];         # UTR 链（来自 intersect 的第 12 列）
    
    # 检查转录本是否存在于 UTR BED 中
    unless (exists $transcripts{$tid}) {
        warn "Warning: Transcript '$tid' not found in UTR BED, skipping\n";
        next;
    }
    
    my $info = $transcripts{$tid};
    my @exons = @{$info->{exons}};
    my $total_len = $info->{total_len};
    
    # 根据 intersect 给出的外显子坐标，在排序列表中找到对应的外显子，获取其前缀
    my ($prefix, $found) = (0, 0);
    foreach my $exon (@exons) {
        if ($exon->{start} == $exon_start && $exon->{end} == $exon_end) {
            $prefix = $exon->{prefix};
            $found = 1;
            last;
        }
    }
    
    # 如果未找到（理论上不会，因为坐标来自同一文件），则尝试用位置扫描作为后备（但不跳过）
    unless ($found) {
        warn "Warning: Exon $exon_start-$exon_end not found in sorted list for $tid, using position scan fallback\n";
        # 后备：通过位置确定外显子
        foreach my $exon (@exons) {
            if ($pos >= $exon->{start} && $pos <= $exon->{end}) {
                $prefix = $exon->{prefix};
                $exon_start = $exon->{start};
                $exon_end   = $exon->{end};
                $found = 1;
                last;
            }
        }
        # 如果还是找不到，则输出比例 0 并标记 S（或跳过，但按用户要求必须输出）
        unless ($found) {
            warn "Warning: PAC $pac_id at $pos cannot be assigned to any exon of $tid, setting ratio=0\n";
            print $fh_out "$pac_id\t$tid\t0\tS\n";
            next;
        }
    }
    
    # 计算 PAC 在当前外显子内的偏移（5' 端为基准）
    my $offset;
    if ($utr_strand eq '+') {
        $offset = $pos - $exon_start;
    } else {   # 负链
        $offset = $exon_end - $pos;
    }
    
    # 不再检查偏移量是否越界，直接计算相对位置
    my $relative_pos = $prefix + $offset;
    
    # 计算比例并限定在 0~100
    my $ratio = ($relative_pos / $total_len) * 100;
    #$ratio = 0   if $ratio < 0;
    #$ratio = 100 if $ratio > 100;
    
    # 分类
    my $flag;
    if ($ratio <= 25) {
        $flag = 'S';
    } elsif ($ratio > 25 && $ratio < 75) {
        $flag = 'M';
    } else {
        $flag = 'L';
    }
    
    # 输出四列：PAC_ID, Transcript_ID, Ratio, Flag
    print $fh_out "$pac_id\t$tid\t".sprintf("%.1f", $ratio)."\t$flag\n";
}

close $fh_in;
close $fh_out;

print "Done. Results saved to $out_file\n";
