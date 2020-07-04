#!/usr/bin/perl
use strict;
use warnings;

my $in;
my $total=0;
my $good=0;
my $name = $ARGV[0];
my @bins=(0,0,0);
my $line = 0;


open IN,"$name";
$name =~ s/\.fastq//i;
open HF2,">${name}_nosgRNA.txt";

while($in=<IN>){
	next unless $line++ % 4 == 1;
	$line = $line % 4;
	$total++;
	if($in =~ /CACC(\w{19,21})GTTT/){
		$good++;
		print $1,"\n";
		$bins[length($1)-19]++;
	}else{
		print HF2 $in;
	}
}

close IN;
close HF2;
#print "#\t",$name,"\t";
#print $good,"\t",$total,"\t",$good/$total,"\t";
#print join("\t",@bins),"\n";