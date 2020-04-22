#!/bin/bash
for item in tsv-summaries/* ; do 

	(head -n 1 $item && tail -n +2 $item | sort -t$'\t' -r -k6 ) | sponge $item

done
