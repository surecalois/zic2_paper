#!/usr/bin/perl
use strict;
use warnings;

my $old="chr0:0:*";
my $loc = 0;
my @counts;
my %table;
my %pam;
my $flag;
#the input file is sorted, do a counting first;
while(<>){
	chomp;
	my @iterms = split /\t/;
	my $sgRNA = $iterms[3];
	my $seq = $iterms[-1];
	if(!same_sgRNA($sgRNA,$old)){
		if($loc > 0){
			($old,$flag) = find_sgRNA(\%table,\%pam,$old);
			my ($chr,$pp,$sig) = split /:/,$old;
			my ($st,$ed) = $sig eq "+" ? ($pp-23,$pp) : ($pp,$pp+23);
			print "$chr\t$st\t$ed\t$old\t.\t$sig\t",join("\t",@counts),"\t$flag\n";
		};
		$loc++;
		undef %table;
		undef %pam;
		
		@counts = @iterms[6..((scalar @iterms)-3)];
		$old = $sgRNA;	
		$table{$sgRNA}++;		
		$pam{$sgRNA} = ($seq =~ /gg$/i) ? 1 : 0;
		
	}else{
		my @temp_con = @iterms[6..((scalar @iterms)-3)];
		map{$counts[$_]+=$temp_con[$_]} (0..((scalar @temp_con)-1));
		$table{$sgRNA}++;
		$pam{$sgRNA} = ($seq =~ /gg$/i) ? 1 : 0;
		#$count++;
	}
	
}

sub find_sgRNA{
	my $table = shift;
	my $pam = shift;
	my $id = shift;
	my $p=$pam->{$id};
	my $c=$table->{$id};
	foreach(keys %{$table}){
		next unless $pam->{$_};
		if($table->{$id} > $c){
			$id = $_;
			$p=$pam->{$id};
	 		$c=$table->{$id};
		}
	}
	return ($id,$p);
}

sub same_sgRNA{
	$_=shift;
	my @a = split /:/;
	$_=shift;
	my @b = split /:/;
	return 0 unless $a[0] eq $b[0];
	return 0 unless $a[2] eq $b[2];
	return 0 unless abs($a[1] - $b[1]) < 7;
	return 1;
}