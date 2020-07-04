#!/bin/bash
screen=$1
echo "maketable.sh {screen}"
maketable.sh ${screen}
echo "loc_sgRNAs.pl ${screen}.fa >${screen}.bed"
loc_sgRNAs.pl ${screen}.fa >${screen}.bed
echo "makeresult.sh ${screen}.bed"
makeresult.sh ${screen}.bed
key=`head -1 res-${screen}.bed|awk '{print NF-2}'`
sort -k $key,$key -rn res-${screen}.bed > res-${screen}.txt
