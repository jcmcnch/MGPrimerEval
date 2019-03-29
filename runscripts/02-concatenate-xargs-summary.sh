#!/bin/bash

find ./11-summary/ -type f -name "*.tsv" -print0 | xargs -0 cat > 190325_all_summary_less_27F.tsv
