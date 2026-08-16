#!/usr/bin/perl
use strict;
use warnings;

die "perl <in> <out>" unless @ARGV==2;

my $in=shift;
my $out=shift;

open IN, "<$in";
open OUT, ">$out";
while(<IN>){
	chomp;
	my @fields = split /\t/;
	next unless @fields >= 6 && $fields[0] eq $fields[3];
	my @numbers = ($fields[1], $fields[2], $fields[4], $fields[5]);
	@numbers = grep { defined $_ && $_ =~ /^-?\d+\.?\d*$/ } @numbers;
	next unless @numbers >= 3;
	@numbers = sort { $a <=> $b } @numbers;
	print OUT "$fields[0]\t$numbers[1]\t$numbers[2]\n";
}
close IN;
close OUT;
