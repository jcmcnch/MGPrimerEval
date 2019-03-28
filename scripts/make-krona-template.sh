#!/bin/bash

if [[ ${#1} -eq 0 ]] ; then
    echo 'Please provide a .tax input file from RDP classifier (output from qiime1).'
    exit 0
fi

filestem=`basename $1 .tax`

python ../scripts/filter-RDP.py $1 > $filestem.tmp

cat $filestem.tmp | sort | uniq -c | sed -e 's/^\s\+\([[:digit:]]\+\)\s\(.\+\)/\1\t\2/' | sed 's/;/\t/g' | sed 's/D_[[:digit:]]__//g' > $filestem.krona.input

source activate kronatools-env

ktImportText -c -n $filestem -o $filestem.html $filestem.krona.input

conda deactivate

rm $filestem.tmp
