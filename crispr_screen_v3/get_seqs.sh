#!/bin/bash
inputbed=$1
#genome=/Users/jiejiaxu/genome/bowtie_index/hg19/hg19
bedtools getfasta -fi $genome.fa -bed $inputbed -tab -s -fo ${inputbed%.bed}_tpam.txt
cut -f2 ${inputbed%.bed}_tpam.txt|paste -d"\t" $inputbed - #> ${inputbed%.bed}_pam.bed
rm ${inputbed%.bed}_tpam.txt
