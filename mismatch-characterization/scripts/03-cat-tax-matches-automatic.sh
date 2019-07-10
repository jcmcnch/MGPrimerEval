#!/bin/bash

mkdir intermediate

group=$1
primer=$2
mismatches=$3

find ../15-matches-classified/individual/ -type f -name "*${group}*${primer}*${mismatches}*tax" -print0 | xargs -0 cat > intermediate/bioGEOTRACES.$group.$primer.$mismatches.hits.all.tax

