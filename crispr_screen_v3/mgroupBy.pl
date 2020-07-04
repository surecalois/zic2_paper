#!/usr/bin/perl
use strict;
use warnings;

my $flag = 1;
while(<>){
    if($flag){
        my @iterms = split /\t/;
        my $fn = scalar @iterms;
        my $c = "-c ".join(",",(5..$fn));
        my $param = $c." -o distinct".",distinct"x($fn-5);
        #print $param,"\n";
        open HF,"|groupBy -g 1,2,3,4 $param" or die;
        $flag = 0;
    }
    print HF;
}
close HF;
