#!/bin/bash

find ./12-full-fastas -type f -name "*.fasta" -print0 | xargs -0 cat > tmp.mismatches.concatenated.fasta 
