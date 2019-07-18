#!/bin/bash

mkdir subsetted-info-files

for item in tax-intermediate/*ids; do

	outfile=`basename $item .ids`.info

	group=`echo $item | cut -d\. -f2`

	primer=`echo $item | cut -d\. -f3`

	infofile=`ls catted-mismatch-info-files/$primer.$group*`
	
	./scripts/subsetInfoFile.py --inputids $item --inputinfo $infofile --output subsetted-info-files/$outfile

done
