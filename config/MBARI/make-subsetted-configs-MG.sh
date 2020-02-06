#!/bin/bash

for item in {01..14} ; do 

	cp config-MBARI-MG.yaml subsetted-MG/config-MBARI-MG-$item.yaml

done

sed -i 's/#926R /926R /' subsetted-MG/*01.yaml
sed -i 's/#806RB /806RB /' subsetted-MG/*02.yaml
sed -i 's/#515Y /515Y /' subsetted-MG/*03.yaml
sed -i 's/#V4F /V4F /' subsetted-MG/*04.yaml
sed -i 's/#V4R /V4R /' subsetted-MG/*05.yaml
sed -i 's/#V4RB /V4RB /' subsetted-MG/*06.yaml
sed -i 's/#338R /338R /' subsetted-MG/*07.yaml
sed -i 's/#341F /341F /' subsetted-MG/*08.yaml
sed -i 's/#785R /785R /' subsetted-MG/*09.yaml
sed -i 's/#27F /27F /' subsetted-MG/*10.yaml
sed -i 's/#1389F /1389F /' subsetted-MG/*11.yaml
sed -i 's/#1510R /1510R /' subsetted-MG/*12.yaml
sed -i 's/#926wF /926wF /' subsetted-MG/*13.yaml
sed -i 's/#1392R /1392R /' subsetted-MG/*14.yaml
