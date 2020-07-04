#!/usr/bin/perl
use warnings;
use strict;
#take a fa file
#my $bed = "/Users/jiejiaxu/genome/annotation/hg19/refGene19_jj.bed";
my $fa = shift @ARGV;
my @sgRNA_bed = loc_sgRNAs($fa);

my @un_target;
my @target;

foreach (@sgRNA_bed){
	if(/^\*/){
		push @un_target,$_;
	}else{
		push @target,$_;
	}
}

#@un_target = map{(split /\t/,$_)[3]} @un_target;

my $fh;
$fa =~ s/\.fa/_ctrl.txt/;
open $fh,"|cut -f4,7- >$fa" or die;
print $fh join("\n",@un_target);
close $fh;

open $fh,"|sortBed" or die;
print $fh join("\n",@target);
close $fh;

sub loc_sgRNAs{
	my $fa = shift;
    my $genome=$ENV{genome};#"/Users/jiejiaxu/genome/bowtie_index/hg19/hg19";
	my $bowtie_para = " -N 1 -L 15";#" -D 15 -R 2";
	my $cmd = "bowtie2 -x $genome -f $fa -p 32 --rdg 10,10 --rfg 10,10".$bowtie_para;
	$fa =~ s/\.fa$/\.sam/;	
	#my $x =`$cmd 2>/dev/null`; #capture STDOUT discard STDERR	
	system("$cmd -S $fa");
	my @res;
	my $fh;
	my $line;
	open $fh,"$fa" or die "no sam files.\n";		
	while($line=<$fh>){
	next if $line =~ /^@/;
	my ($name,$flag,$chr,$loc,$seq) = (split /\t/,$line)[0,1,2,3,9];
	my $sig = $flag^0x0010?'+':'-';
	my ($loc1,$loc2);
	my $len = length $seq;
# 	my $mdz = ($line =~ /\tMD:Z:([\d\w]+)\t/)[0];
	
	$loc1 = $loc-1;
	$loc2 = $loc+$len-1;

	if($seq !~ /^CCN|NGG$/){	
		if($flag^0x0010){
			#+strand
			$loc2 = $loc2+3;	
		}else{
			#-strand
			$loc1 = $loc1-3;	
		}
	}
    my $sgRNAname= $sig eq "+" ? "$chr:$loc2:+" : "$chr:$loc1:-";
	my $out = "$chr\t$loc1\t$loc2\t$sgRNAname\t.\t$sig";
       $out.="\t".join("\t",split /:/,$name);
       $out.="\t".$seq;
       push @res,$out;
	}
	close $fh;
	return @res;
}