#!/bin/bash
source activate repeatmasker-env

if [[ ${#1} -eq 0 ]] ; then
    echo 'Please enter an input fasta file as the first argument after the script name.'
    exit 0
fi

RepeatMasker -pa 5 -s -noisy -html -xm -source -poly -species arabidopsis $1

source deactivate
