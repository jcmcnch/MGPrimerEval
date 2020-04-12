#!/bin/bash

for item in {01..14} ; do 

	cp config-MBARI-MT.yaml subsetted-MT/config-MBARI-MT-$item.yaml

done

sed -i 's/#926R /926R /' subsetted-MT/*01.yaml
sed -i 's/#806RB /806RB /' subsetted-MT/*02.yaml
sed -i 's/#515Y /515Y /' subsetted-MT/*03.yaml
sed -i 's/#V4F /V4F /' subsetted-MT/*04.yaml
sed -i 's/#V4R /V4R /' subsetted-MT/*05.yaml
sed -i 's/#V4RB /V4RB /' subsetted-MT/*06.yaml
sed -i 's/#338R /338R /' subsetted-MT/*07.yaml
sed -i 's/#341F /341F /' subsetted-MT/*08.yaml
sed -i 's/#785R /785R /' subsetted-MT/*09.yaml
sed -i 's/#27F /27F /' subsetted-MT/*10.yaml
sed -i 's/#1389F /1389F /' subsetted-MT/*11.yaml
sed -i 's/#1510R /1510R /' subsetted-MT/*12.yaml
sed -i 's/#926wF /926wF /' subsetted-MT/*13.yaml
sed -i 's/#1392R /1392R /' subsetted-MT/*14.yaml
