#!/bin/bash

for item in `ls 03-low-complexity-filtered/ | cut -f1 -d\. | sort | uniq`; do

	 SSUnum=`cat 03-low-complexity-filtered/$item* | grep -c "^>"`

	 totalSeqNum=`zcat /data/bkp1/00_clean_MG_reads/$item* | grep -c "^@"`

         SSUFrac=`bc <<< "scale=8; $SSUnum/$totalSeqNum"` 
        
         printf "$item\t$SSUFrac\n"

done
