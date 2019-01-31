#!/bin/bash

if [[ ${#1} -eq 0 ]] ; then
    echo 'Please enter an input fasta file that has repeat regions masked (converted to Ns).'
    exit 0
fi

if [[ ${#2} -eq 0 ]] ; then
    echo 'Please enter a file name for an output fasta file that will contain the non-N-containing sequences.'
    exit 0
fi

if [[ ${#3} -eq 0 ]] ; then
    echo 'Please enter a file name for an output fast that contain the N-containing sequences.'
    exit 0
fi

prinseq-lite.pl -ns_max_n 0 -fasta $1 -out_good $2 -out_bad $3
