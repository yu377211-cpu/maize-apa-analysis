use strict;
use warnings;

die "perl $0 <in.fa> <out.fa>" unless @ARGV==2;

my $in=shift;
my $out=shift;
open IN,"<$in";
open OUT,">$out";

my %complement = ('A' => 'T', 'T' => 'A', 'G' => 'C', 'C' => 'G', 'a' => 't', 't' => 'a', 'g' => 'c', 'c' => 'g');

while(<IN>){
	chomp;
	if(/^>/){
		print OUT "$_\n";
	}else{
		my $reversed_seq = reverse $_;
		my $complement_seq = '';
		foreach my $base (split //, $reversed_seq) {
			if (exists $complement{$base}) {
				$complement_seq .= $complement{$base};
			}else{
				$complement_seq .= $base;
			}
		}
		print OUT "$complement_seq\n";
	}
}
close IN;
close OUT;
