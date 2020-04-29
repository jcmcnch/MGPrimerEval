#!/bin/bash

find ./11-summary/ -type f -name "*.tsv" -print0 | xargs -0 cat > 190317_all_summary_less_27F.tsv
