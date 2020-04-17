#!/bin/bash

find ./14-subsampled-matched-fastas -type f -name "*.fasta" -print0 | xargs -0 cat > tmp.mismatches.subsampled.concatenated.fasta 
