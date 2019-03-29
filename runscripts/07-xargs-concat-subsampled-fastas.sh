#!/bin/bash

find ./14-subsampled-matched-fastas -type f -name "*.fasta" -print0 | xargs -0 cat > tmp.matches.subsampled.concatenated.fasta
