#!/bin/bash
ko_lab=$1
tmp=''
tmp2=''
for ii in *.fastq
do
	#echo $ii
	extrac_sgRNA.pl $ii > ${ii%.fastq}.seq
	sortcount.sh ${ii%.fastq}.seq > ${ii%.fastq}.txt
    sortcount.sh ${ii%.fastq}_nosgRNA.txt > nosgRNA-${ii%.fastq}.txt
    rm ${ii%.fastq}.seq
    rm ${ii%.fastq}_nosgRNA.txt
	tmp="${tmp} ${ii%.fastq}.txt"
    tmp2="${tmp2} nosgRNA-${ii%.fastq}.txt"
done
mergelib.pl $tmp >${ko_lab}-table.txt
rm $tmp
mergelib.pl $tmp2 >${ko_lab}-nosgRNA-table.txt
rm $tmp2
perl -ne 'chomp;@a=split /^\w+\K\t/;$a[1] =~ s/\t/:/g;print ">$a[1]\n$a[0]\n"' \
${ko_lab}-table.txt > ${ko_lab}.fa

#|perl -pe 's/\t/:/g' \
#|awk -F'NGG:' '{print ">"$2"\n"$1"NGG"}'

#join -a1 -a2 -1 1 -2 1 -o 0,1.2,2.2 -e "0" $tmp >join.txt
#awk -F' ' '$2+$3 > 100' join.txt >fed_join.txt
