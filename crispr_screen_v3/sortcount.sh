#!/bin/bash
inputlist=$1
sort $inputlist|uniq -c\
|sed "s/^\s+//" \
|awk -F' ' '{print $2"\t"$1}'