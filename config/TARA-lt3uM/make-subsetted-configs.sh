#!/bin/bash

mkdir -p subsetted

for item in {01..14} ; do 

	cp config-TARA-lt3uM.yaml subsetted/config-TARA-lt3uM-$item.yaml

done

sed -i 's/#926R /926R /' subsetted/*01.yaml
sed -i 's/#806RB /806RB /' subsetted/*02.yaml
sed -i 's/#515Y /515Y /' subsetted/*03.yaml
sed -i 's/#V4F /V4F /' subsetted/*04.yaml
sed -i 's/#V4R /V4R /' subsetted/*05.yaml
sed -i 's/#V4RB /V4RB /' subsetted/*06.yaml
sed -i 's/#338R /338R /' subsetted/*07.yaml
sed -i 's/#341F /341F /' subsetted/*08.yaml
sed -i 's/#785R /785R /' subsetted/*09.yaml
sed -i 's/#27F /27F /' subsetted/*10.yaml
sed -i 's/#1389F /1389F /' subsetted/*11.yaml
sed -i 's/#1510R /1510R /' subsetted/*12.yaml
sed -i 's/#926wF /926wF /' subsetted/*13.yaml
sed -i 's/#1392R /1392R /' subsetted/*14.yaml
