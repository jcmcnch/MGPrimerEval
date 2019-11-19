#!/bin/bash

mkdir paired-kronas

PS3='Please choose a group of interest: '
options=("Archaea" "Bacteria" "Cyanobacteria (including Chloroplasts)" "Eukaryotes")
select opt in "${options[@]}"
do
    case $opt in
        "Archaea")
            group="ARCH" && echo "you chose $opt"; break
            ;;
        "Bacteria")
            group="BACT-NON-CYANO" && echo "you chose $opt"; break
            ;;
        "Cyanobacteria (including Chloroplasts)")
            group="BACT-CYANO" && echo "you chose $opt"; break
            ;;
        "Eukaryotes")
            group="EUK" && echo "you chose $opt"; break
            ;;
        *) echo "invalid option $REPLY";;
    esac
done

PS3='Please choose a primer of interest: '
options=("515Y" "806RB" "926R")
select opt in "${options[@]}"
do
    case $opt in
        "515Y")
            primer=$opt && echo "you chose $opt"; break
            ;;
        "806RB")
            primer=$opt && echo "you chose $opt"; break
            ;;
        "926R")
            primer=$opt && echo "you chose $opt"; break
            ;;
        *) echo "invalid option $REPLY";;
    esac
done

PS3='Please choose the number of mismatches: '
options=("0-mismatch" "1-mismatch" "2-mismatch")
select opt in "${options[@]}"
do
    case $opt in
        "0-mismatch")
            mismatches=$opt && echo "you chose $opt"; break
            ;;
        "1-mismatch")
            mismatches=$opt && echo "you chose $opt"; break
            ;;
        "2-mismatch")
            mismatches=$opt && echo "you chose $opt"; break
            ;;
        *) echo "invalid option $REPLY";;
    esac
done

for item in `ls 13-classified/individual/*${group}*${primer}*${mismatches}*`; do

  cat $item >> paired-kronas/$group.$primer.$mismatches.mismatches-to-primer.tax

done

for item in `ls 15-matches-classified/individual/*${group}*${primer}*${mismatches}*`; do

  cat $item >> paired-kronas/$group.$primer.$mismatches.matches-to-primer.subsample.tax

done

cd paired-kronas

../scripts/make-krona-template.sh $group.$primer.$mismatches.mismatches-to-primer.tax
../scripts/make-krona-template.sh $group.$primer.$mismatches.matches-to-primer.subsample.tax

cd ..
