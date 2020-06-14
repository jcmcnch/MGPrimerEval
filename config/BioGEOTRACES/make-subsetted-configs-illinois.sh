#!/bin/bash

mkdir -p subsetted-illinois

for item in {01..16} ; do 

	cp config-BioGEOTRACES-illinois-EUK-primers.yaml subsetted-illinois/config-BioGEOTRACES-$item.yaml

done

sed -i 's/#926R /926R /' subsetted-illinois/*01.yaml
sed -i 's/#806RB /806RB /' subsetted-illinois/*02.yaml
sed -i 's/#515Y /515Y /' subsetted-illinois/*03.yaml
sed -i 's/#V4F /V4F /' subsetted-illinois/*04.yaml
sed -i 's/#V4R /V4R /' subsetted-illinois/*05.yaml
sed -i 's/#V4RB /V4RB /' subsetted-illinois/*06.yaml
sed -i 's/#338R /338R /' subsetted-illinois/*07.yaml
sed -i 's/#341F /341F /' subsetted-illinois/*08.yaml
sed -i 's/#785R /785R /' subsetted-illinois/*09.yaml
sed -i 's/#27F /27F /' subsetted-illinois/*10.yaml
sed -i 's/#1389F /1389F /' subsetted-illinois/*11.yaml
sed -i 's/#1510R /1510R /' subsetted-illinois/*12.yaml
sed -i 's/#926wF /926wF /' subsetted-illinois/*13.yaml
sed -i 's/#1392R /1392R /' subsetted-illinois/*14.yaml
sed -i 's/#F566Euk /F566Euk /' subsetted-illinois/*15.yaml
sed -i 's/#R1200Euk /R1200Euk /' subsetted-illinois/*16.yaml
