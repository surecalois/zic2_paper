#!/usr/bin/perl
#sgRNA_target
use strict;
use warnings;

my $genome=$ENV{genome};#"/Users/jiejiaxu/genome/bowtie_index/hg19/hg19";
my $bed=$ENV{bed};#"/Users/jiejiaxu/genome/annotation/hg19/refgene19_jj.bed";

while(my $sgRNA=<>){
	chomp $sgRNA;
	my $cmd = "bowtie2 -x $genome -c $sgRNA --no-head";
	my $line =`$cmd 2>/dev/null`;
	chomp $line;	
	my ($flag,$chr,$loc,$seq) = (split /\t/,$line)[1,2,3,9];
	print "$sgRNA\ttarget not found!\n" and next if $chr eq "*";
	my $sig = $flag^0x0010?'+':'-';
	my $len = length $seq;	
	my $loc1 = $loc;#faidx no need to -1
	my $loc2 = $loc+$len-1;
	if($flag^0x0010){
		#+strand
		$loc2 = $loc2+3;	
	}else{
		#-strand
		$loc1 = $loc1-3;	
	}
	my $mdz = ($line =~ /\tNM:i:([\d\w]+)\t/)[0];	
	$cmd="samtools faidx $genome.fa $chr:$loc1-$loc2\n";
	my $x = qx/$cmd/;
	my @fa = split /\n/,$x;
	my $fa_note = shift @fa;
	my $fa_dna = join '',@fa;
	
	if($sig eq '-'){
		#$seq=$fa_dna;
		$fa_dna = reverse $fa_dna;
		$fa_dna =~ tr/ATCGatcg/TAGCtagc/;
	}
	
	print "$sgRNA\tNM:i:$mdz\t";
	print "$chr:$loc1-$loc2:$sig\t$fa_dna\t";
	
	my $fh;
	open $fh,"|bedtools intersect -a stdin -b $bed -loj|head -1|cut -f10";
	print $fh "$chr\t$loc1\t$loc2\tnoname\t$mdz\t$sig";
	close $fh;
	#print "\n";
}