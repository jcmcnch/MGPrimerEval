#!/bin/bash

group=$1
primer=$2
mismatches=$3

find ../13-classified/individual/ -type f -name "*${group}*${primer}*${mismatches}*tax" -print0 | xargs -0 cat > intermediate/bioGEOTRACES.$group.$primer.$mismatches.nohits.all.tax

