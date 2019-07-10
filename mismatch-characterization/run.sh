#!/bin/bash

for item in `cat primers.txt`; do

	primer=$item

	for item in `cat groups.txt`; do

		group=$item

		./scripts/03-cat-tax-matches-automatic.sh $group $primer 0-mismatch

		for item in `cat mismatches.txt`; do

			mismatches=$item

			./scripts/00-cat-tax-mismatches-automatic.sh $group $primer $mismatches

		done

	done

done

./scripts/01-count-tax-mismatches.sh

./scripts/02-filter-mismatches-by-abundance.sh

./scripts/04-count-tax-matches.sh

./scripts/05-filter-matches-by-abundance.sh

./scripts/06-normalize-match-counts-by-total-seqs.sh

./scripts/07-compute-frac-mismatched.sh
