#!/bin/bash
inputbed=$1
#bed=/Users/jiejiaxu/genome/annotation/hg19/refgene19_jjgene.bed
bedtools intersect -a $inputbed -b $bed -loj \
|perl -lane '$F[4]=$F[-3];splice(@F,-6);print join("\t",@F)'
#|perl -pe 'my @iterms = split /\t/;$iterms[4]=$iterms[-3];splice(@iterms,-6);$_ = join("\t",@iterms)."\n"' #>${inputbed%.bed}_gene.bed