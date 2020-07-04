#!/bin/bash
alignbed=$1
get_seqs.sh $alignbed|sort -k4|std_sgRNA.pl >std-$alignbed

get_genes.sh std-$alignbed |mgroupBy.pl >gene-$alignbed 

get_seqs.sh gene-$alignbed >res-$alignbed
rm std-$alignbed
rm gene-$alignbed
