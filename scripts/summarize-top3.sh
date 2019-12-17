#!/bin/bash
tail -n+2 $1 | sort -k6 -gr -t$'\t' | head -n3 | cut -d$'\t' -f1,6 | sed -E 's/.*o:(.*)\t(.+)/\1\t\2/g' | awk '{printf $1 "(""%0.3f",$2")" ; printf "), "}' | sed 's/, /\n/3' > $2
