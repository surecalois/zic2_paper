#!/bin/bash
insam=$1
awk '$3 != "*"' $insam |cut -f 1,3 \
|perl -lane '$F[0] =~ s/:/\t/g;print "$F[1]\t$F[0]"' > ${insam%.sam}-align.txt
sort ${insam%.sam}-align.txt >z.txt
../groupTable.pl z.txt >zg.txt
sort -k5 -rn zg.txt >szg.txt
grep NonTargetingControlGuideForHuman szg.txt >ctrl_szg.txt
grep hsa-mir szg.txt >mir_szg.txt
grep -v hsa-mir szg.txt | grep -v NonTargetingControlGuideForHuman > gene_szg.txt
