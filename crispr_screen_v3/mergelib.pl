#!/usr/bin/perl
use strict;
use warnings;
my %table;
my @files = @ARGV;
my $file_id=0;

foreach my $file (@files){
    my %subtable;
    my @item;
    
    open HF, $file;
    while(<HF>){
        chomp;
        @item = split /\s/;
        $subtable{$item[0]} += $item[1];
    }
    close HF;
    
    if($file_id == 0){
        %table = %subtable;
    }else{
        merge(\%table,\%subtable,$file_id);
    }
    $file_id++;
}

my $total=0;
foreach my $gene (keys %table){
	my @counts = split /\t/,$table{$gene};
	$total = 0;
	map{$total += $_} @counts;
	push @counts,$total;
	$table{$gene} = \@counts;
}

foreach (sort {${$table{$b}}[-1] <=> ${$table{$a}}[-1]} keys %table){
    print $_,"\t",join("\t",@{$table{$_}}),"\n";
}

#map{print $_,"\t",join("\t",@{$table{$_}}),"\n"} keys %table;

#sub functions
sub merge{
	my $a = shift;
	my $b = shift;
	my $n = shift;
	
		#if the key exists in both a and b then append b.key's value
		#if the key exists in b but not in a then insert 0
		foreach(keys %$b){
			if(exists $$a{$_}){
				$$a{$_} = "$$a{$_}\t$$b{$_}";
			}else{
				$$a{$_} = "0\t"x$n."$$b{$_}";
			}
		}
		
		#if the key exists in a but not in b then append 0
		foreach(keys %$a){
			$$a{$_} = "$$a{$_}\t0" unless exists $$b{$_};
		}

}