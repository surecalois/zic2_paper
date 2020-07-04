#!/usr/bin/perl
use strict;
use warnings;

my $in=<>;
chomp $in;
my @old = split /\t/,$in;
#my @outs;
my $n = scalar @old;
while(<>){
	chomp;
	my @iterms = split /\t/;
	if($iterms[0] eq $old[0]){
		for my $ii (1..($n-1)){
			$old[$ii]+=$iterms[$ii];
		}
	}else{
		#push @outs,join("\t",@old);
		print join("\t",@old),"\n";
		@old=@iterms;
	}
}
