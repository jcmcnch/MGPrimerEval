#!/bin/bash
mkdir mismatch-locations-test

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
options=("515Y" "806RB" "926R" "V4F" "V4RB" "341F" "785R" "27F" "338R")
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
        "V4F")
            primer=$opt && echo "you chose $opt"; break
            ;;
        "V4RB")
            primer=$opt && echo "you chose $opt"; break
            ;;
        "341F")
            primer=$opt && echo "you chose $opt"; break
            ;;
        "785R")
            primer=$opt && echo "you chose $opt"; break
            ;;
        "27F")
            primer=$opt && echo "you chose $opt"; break
            ;;
        "338R")
            primer=$opt && echo "you chose $opt"; break
            ;;
        *) echo "invalid option $REPLY";;
    esac
done

PS3='Please choose the number of mismatches: '
options=("1-mismatch" "2-mismatch" "30pc-mismatch")
select opt in "${options[@]}"
do
    case $opt in
        "1-mismatch")
            mismatches=$opt && echo "you chose $opt"; break
            ;;
        "2-mismatch")
            mismatches=$opt && echo "you chose $opt"; break
            ;;
				"30pc-mismatch")
            mismatches=$opt && echo "you chose $opt"; break
            ;;
        *) echo "invalid option $REPLY";;
    esac
done



for item in `ls archive/10-checked/926R/0-mismatch/*nohit.fastq`; do

	cutadapt --overlap=20 --no-indels --no-trim -b AAACTYAAAKRAATTGRCGG --error-rate=0.3 --info-file=mismatch-locations-test/`basename $item .fastq`.info.tsv $item

done
