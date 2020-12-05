#!/bin/bash -i

#any old pandas-containing environment will do
conda activate opedia-env

#cleanup previous failed runs, if files exist
rm classify-workflow-intermediate/11-taxa-with-many-mismatches/*info
rm classify-workflow-intermediate/11-taxa-with-many-mismatches/*ids
rm -rf output-classify-workflow/taxa-mismatch-summaries

mkdir -p output-classify-workflow/taxa-mismatch-summaries

for primer in 926R 806RB 515Y V4F V4R V4RB 341F 785R 27F 1389F 1510R ; do

	for group in ARCH BACT-CYANO BACT-NON-CYANO EUK ; do

		for taxon in `cat classify-workflow-intermediate/11-taxa-with-many-mismatches/*.$group.$primer.0-mismatch.gt50mm-and-5pcmismatched-taxa.txt`; do

			cleantaxon=`echo $taxon | sed 's/://g' | sed 's/,/-/g'`
			input=`ls classify-workflow-intermediate/11-taxa-with-many-mismatches/*.$group.$primer.0-mismatch.gt50mm-and-5pcmismatched-taxa.headers.tsv`
			study=`basename $input | cut -f1 -d\.`
			output=classify-workflow-intermediate/11-taxa-with-many-mismatches/$study.$group.$primer.$cleantaxon.mismatched.ids
			info=classify-workflow-intermediate/11-taxa-with-many-mismatches/$study.$group.$primer.$cleantaxon.mismatched.info
		
			grep "$taxon" $input | cut -f1 > $output
			
			#add loop because grep -f takes forever and eats all of kraken's RAM!
			while read line ; do

				grep "$line" classify-workflow-intermediate/10-concatenated-info-files/*.$group.$primer.2-mismatch.nohits.all.info >> $info

			done < $output

			./scripts/make-mismatch-alignments-and-summarize.py --info $info --alignmentout output-classify-workflow/taxa-mismatch-summaries/$study.$group.$primer.$cleantaxon.ali.fasta --summaryout output-classify-workflow/taxa-mismatch-summaries/$study.$group.$primer.$cleantaxon.summary.tsv 

		done

	done

done
