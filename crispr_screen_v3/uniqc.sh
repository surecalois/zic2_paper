#!/bin/bash
uniq -c|perl -lne 's/^\s*(\d*)\s/$1\t/;@F= split /\t/;$a=shift @F;push @F,$a; print join("\t",@F)'