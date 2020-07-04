#!/usr/bin/perl
use strict;
use warnings;
my $insam = shift @ARGV;

open IN,"$insam";
open HF,"align-$insam"
open HF2,"noalign-$insam"

while(<>){
	next if /^@/;
	chomp;
	my @iterms = split /\t/;
	
}
close IN;
close HF;
close HF2;