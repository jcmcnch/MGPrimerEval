ruleorder: count_tax_matches_and_mismatches_cyano > count_tax_matches_and_mismatches
CUTOFF = config["cutoff"]
VSEARCHudbPath=config["VSEARCHudbPath"]
PhytoRefUdbPath=config["PhytoRefUdbPath"]

rule all:
	input:
		expand("{study}.EUK.pdf", study=config["study"]),
		expand("intermediate/classify-workflow/01-mismatches-classified/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.nohit.filtered.VSEARCHsintax-SILVA132.tax", sample=config["samples"], study=config["study"], group=["ARCH","BACT-NON-CYANO","EUK"], primer=config["primer"], mismatches=["0-mismatch", "1-mismatch", "2-mismatch"], direction=['fwd','rev']),
		expand("intermediate/classify-workflow/03-matches-classified/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.sub5k.hit.filtered.VSEARCHsintax-SILVA132.tax", sample=config["samples"], study=config["study"], group=["ARCH","BACT-NON-CYANO","EUK"], primer=config["primer"], mismatches=["0-mismatch", "1-mismatch", "2-mismatch"], direction=['fwd','rev']),
		#only one target necessary for phytoRef because mismatches and matches are classified in a single rule (no subsampling)
		expand("intermediate/classify-workflow/03-matches-classified/{sample}.SSU.{direction}.BACT-CYANO.{primer}.{mismatches}.sub5k.hit.filtered.VSEARCHsintax-PhytoRef.tax", sample=config["samples"], study=config["study"], primer=config["primer"], mismatches=["0-mismatch", "1-mismatch", "2-mismatch"], direction=['fwd','rev']),
		expand("intermediate/classify-workflow/07-normalized-counts/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.normalized.tsv", study=config["study"], group=config["groups"], primer=config["primer"], mismatches=["0-mismatch", "1-mismatch", "2-mismatch"]),
		expand("intermediate/classify-workflow/11-taxa-with-many-mismatches/{study}.{group}.{primer}.0-mismatch.gt50mm-and-5pcmismatched-taxa.headers.tsv", study=config["study"], group=config["groups"], primer=config["primer"]),
		expand("output/classify-workflow/overall-summaries/{study}.{group}.{primer}.{mismatches}.summary.tsv", sample=config["samples"], study=config["study"], group=config["groups"], primer=config["primer"], mismatches=["0-mismatch", "1-mismatch", "2-mismatch"]),
		expand("output/classify-workflow/overall-summaries/{study}.{group}.{primer}.{mismatches}.aln.summary.tsv", sample=config["samples"], study=config["study"], group=config["groups"], primer=config["primer"], mismatches=["0-mismatch", "1-mismatch", "2-mismatch"]),
		expand("output/classify-workflow/overall-summaries/{study}.{group}.{primer}.taxonFracMismatched.0-2mm.tsv", study=config["study"], group=config["groups"], primer=config["primer"]),
		expand("output/classify-workflow/plots/matchVSmismatch-barplots/{study}.{group}.{primer}.taxonFracMismatched.{mismatches}.pdf", study=config["study"], group=config["groups"], primer=config["primer"], mismatches=["0-mismatch", "1-mismatch", "2-mismatch"]),
		expand("output/classify-workflow/summary-mismatch-overlap-primer-pairs/{study}.{group}.{primer_pair}.avgCase.tsv", primer_pair=config["primer_pairs"], study=config["study"], group=config["groups"]),
		expand("output/classify-workflow/pasted-summaries/{study}.{group}.{primer_pair}.pasted.tsv", primer_pair=config["primer_pairs"], study=config["study"], group=config["groups"]),
		expand("output/classify-workflow/normalized-summaries/{study}.{group}.{primer_pair}.normalized.tsv", primer_pair=config["primer_pairs"], study=config["study"], group=config["groups"])

#this rule assumes your compute pipeline is done! remove and rerun if you add more samples to your pipeline
rule concatenate_compute_results:
	output:
		"{study}.compute-results.tsv"
	shell:
		"find ./output/compute-workflow/09-summary/ -type f -name \"*.tsv\" -print0 | xargs -0 cat > {output}"

rule plot_compute_results:
	input:
		"{study}.compute-results.tsv"
	params:
		"{study}"
	output:
		"{study}.EUK.pdf",
		"{study}.BACT-NON-CYANO.pdf",
		"{study}.BACT-CYANO.pdf",
		"{study}.ARCH.pdf"
	conda:
		"envs/ggplot2.yaml"
	shell:
		"scripts/plot-individual-primers-shell.R {params} {input} {output} || touch {output}"

rule classify_mismatches:
	input:
		"intermediate/compute-workflow/08-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.nohit.filtered.fastq",
	output:
		"intermediate/classify-workflow/01-mismatches-classified/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.nohit.filtered.VSEARCHsintax-SILVA132.tax"
	conda:
		"envs/vsearch.yaml"
	params:
		db=VSEARCHudbPath
	shell:
		#Double pipe is OR operator and will only be executed if vsearch returns an error. Necessary otherwise empty files will cause vsearch to fail.
		"""
		vsearch --sintax {input} --db {params.db} \
		--tabbedout {output} --threads 1 --sintax_cutoff 0 --top_hits_only --topn 1 --notrunclabels \
		|| touch {output}
		"""

#Take up to 5000 reads from the matched files
rule subsample_matched_fastqs:
	input:
		"intermediate/compute-workflow/08-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.hit.filtered.fastq",
	output:
		"intermediate/classify-workflow/02-matches-subsampled/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.sub5k.hit.filtered.fastq"
	conda:
		"envs/bbmap.yaml"
	shell:
		"reformat.sh samplereadstarget=5000 in={input} out={output} "

rule classify_matches_subsample:
	input:
		"intermediate/classify-workflow/02-matches-subsampled/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.sub5k.hit.filtered.fastq"
	output:
		"intermediate/classify-workflow/03-matches-classified/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.sub5k.hit.filtered.VSEARCHsintax-SILVA132.tax"
	conda:
		"envs/vsearch.yaml"
	params:
                db=VSEARCHudbPath
	shell:
		#Double pipe is OR operator and will only be executed if vsearch returns an error. Necessary otherwise empty files will cause vsearch to fail.
		"""
                vsearch --sintax {input} --db {params.db} \
		--tabbedout {output} --threads 1 --sintax_cutoff 0 --top_hits_only --topn 1 --notrunclabels \
		|| touch {output}
		"""
#Don't worry about subsampling, since there aren't that many chloroplast sequences
rule classify_cyano_fraction_phytoRef:
	input:
		nohits="intermediate/compute-workflow/08-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.BACT-CYANO.{primer}.{mismatches}.nohit.filtered.fastq",
		hits="intermediate/classify-workflow/02-matches-subsampled/{sample}.SSU.{direction}.BACT-CYANO.{primer}.{mismatches}.sub5k.hit.filtered.fastq"
	output:
		nohits="intermediate/classify-workflow/01-mismatches-classified/{sample}.SSU.{direction}.BACT-CYANO.{primer}.{mismatches}.nohit.filtered.VSEARCHsintax-PhytoRef.tax",
		hits="intermediate/classify-workflow/03-matches-classified/{sample}.SSU.{direction}.BACT-CYANO.{primer}.{mismatches}.sub5k.hit.filtered.VSEARCHsintax-PhytoRef.tax"
	conda:
		"envs/vsearch.yaml"
	params:
		db=PhytoRefUdbPath
	shell:
                #Double pipe is OR operator and will only be executed if vsearch returns an error. Necessary otherwise empty files will cause vsearch to fail.
                """
                vsearch --sintax {input.nohits} --db {params.db} \
                --tabbedout {output.nohits} --threads 1 --sintax_cutoff 0 --top_hits_only --topn 1 --notrunclabels \
                || touch {output.nohits}

		vsearch --sintax {input.hits} --db {params.db} \
                --tabbedout {output.hits} --threads 1 --sintax_cutoff 0 --top_hits_only --topn 1 --notrunclabels \
                || touch {output.hits}
                """


"""
Rules below comprise a workflow for generating summaries of
which taxa are most discriminated against by a particular primer set.
Implemented using common bash tools and tested on Ubuntu 16.04, not tested on other systems.
Results are not per-sample, but rather across an entire dataset.
IMPORTANT NOTE: it would be great to have a checkpoint at this stage to tell whether all the above files have been generated.
However, I haven't implemented this yet, so am just using the --until flag in a bash wrapper script to get this behaviour (see runscripts folder for example).
"""

#Concatenate taxonomy files, make sure everything is generated from above steps
rule cat_tax_for_all_samples_matches_and_mismatches:
	output:
		mismatches="intermediate/classify-workflow/04-tax-concatenated/{study}.{group}.{primer}.{mismatches}.nohits.all.tax",
		matches="intermediate/classify-workflow/04-tax-concatenated/{study}.{group}.{primer}.{mismatches}.hits.all.tax"
	shell:
		"find intermediate/classify-workflow/01-mismatches-classified/ -type f -name "
		"\"*{wildcards.group}*{wildcards.primer}*{wildcards.mismatches}*tax\" -print0 | "
		"xargs -0 cat > {output.mismatches} ; "
		"find intermediate/classify-workflow/03-matches-classified/ -type f -name "
		"\"*{wildcards.group}*{wildcards.primer}*{wildcards.mismatches}*tax\" -print0 | "
		"xargs -0 cat > {output.matches}"

#Counting order-level groupings (can adjust level with the "cut -d, -f1-4" parameter below)
rule count_tax_matches_and_mismatches:
	input:
		mismatches="intermediate/classify-workflow/04-tax-concatenated/{study}.{group}.{primer}.{mismatches}.nohits.all.tax",
		matches="intermediate/classify-workflow/04-tax-concatenated/{study}.{group}.{primer}.{mismatches}.hits.all.tax"
	output:
		matches="intermediate/classify-workflow/05-tax-counts/{study}.{group}.{primer}.{mismatches}.hits.all.order.counts.tsv",
		mismatches="intermediate/classify-workflow/05-tax-counts/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.tsv",
		taxTableMatches="intermediate/classify-workflow/05-tax-table/{study}.{group}.{primer}.{mismatches}.hits.all.order.counts.taxtable",
		taxTableMismatches="intermediate/classify-workflow/05-tax-table/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.taxtable"
	shell:
		"sed -re 's/\([0-9]{{1}}\.[0-9]{{2}}\)//g' {input.matches} | tee {output.taxTableMatches} |" #Remove confidence estimations from VSEARCH output, keep a copy for later steps but also pipe to subsequent commands
		"cut -f2 | sort | cut -d, -f1-4 | sort | uniq -c | " #Take only tax column, collapse to order level, then count unique occurrences
		"tail -f -n +2 | awk '{{print $1,\"\t\",$2}}' > {output.matches} ; " #Process output into tsv format to stdout
		"sed -re 's/\([0-9]{{1}}\.[0-9]{{2}}\)//g' {input.mismatches} | tee {output.taxTableMismatches} |" #Remove confidence estimations from VSEARCH output, keep a copy for later steps but also pipe to subsequent commands
		"cut -f2 | sort | cut -d, -f1-4 | sort | uniq -c | " #Take only tax column, collapse to order level, then count unique occurrences
		"tail -f -n +2 | awk '{{print $1,\"\t\",$2}}' > {output.mismatches}" #Process output into tsv format to stdout

#Counting order-level groupings for cyano at level 6
rule count_tax_matches_and_mismatches_cyano:
        input:
                mismatches="intermediate/classify-workflow/04-tax-concatenated/{study}.BACT-CYANO.{primer}.{mismatches}.nohits.all.tax",
                matches="intermediate/classify-workflow/04-tax-concatenated/{study}.BACT-CYANO.{primer}.{mismatches}.hits.all.tax"
        output:
                matches="intermediate/classify-workflow/05-tax-counts/{study}.BACT-CYANO.{primer}.{mismatches}.hits.all.order.counts.tsv",
                mismatches="intermediate/classify-workflow/05-tax-counts/{study}.BACT-CYANO.{primer}.{mismatches}.nohits.all.order.counts.tsv",
                taxTableMatches="intermediate/classify-workflow/05-tax-table/{study}.BACT-CYANO.{primer}.{mismatches}.hits.all.order.counts.taxtable",
                taxTableMismatches="intermediate/classify-workflow/05-tax-table/{study}.BACT-CYANO.{primer}.{mismatches}.nohits.all.order.counts.taxtable"
        shell:
                "sed -re 's/\([0-9]{{1}}\.[0-9]{{2}}\)//g' {input.matches} | tee {output.taxTableMatches} |" #Remove confidence estimations from VSEARCH output, keep a copy for later steps but also pipe to subsequent commands
                "cut -f2 | sort | cut -d, -f1-6 | sort | uniq -c | " #Take only tax column, collapse to 6th level (appropriate for chloroplasts), then count unique occurrences
                "tail -f -n +2 | awk '{{print $1,\"\t\",$2}}' > {output.matches} ; " #Process output into tsv format to stdout
                "sed -re 's/\([0-9]{{1}}\.[0-9]{{2}}\)//g' {input.mismatches} | tee {output.taxTableMismatches} |" #Remove confidence estimations from VSEARCH output, keep a copy for later steps but also pipe to subsequent commands
                "cut -f2 | sort | cut -d, -f1-6 | sort | uniq -c | " #Take only tax column, collapse to order level, then count unique occurrences
                "tail -f -n +2 | awk '{{print $1,\"\t\",$2}}' > {output.mismatches}" #Process output into tsv format to stdout

rule transform_tax_matches_to_proportions:
	input:
		matches="intermediate/classify-workflow/05-tax-counts/{study}.{group}.{primer}.{mismatches}.hits.all.order.counts.tsv",
		mismatches="intermediate/classify-workflow/05-tax-counts/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.tsv"
	output:
		matches="intermediate/classify-workflow/06-tax-counts-fractions/{study}.{group}.{primer}.{mismatches}.hits.all.order.counts.frac.tsv",
		mismatches="intermediate/classify-workflow/06-tax-counts-fractions/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.frac.tsv"
	shell:
		"scripts/mismatch-characterization/transform-to-fractional-abundance.py {input.matches} > {output.matches} ; "
		"scripts/mismatch-characterization/transform-to-fractional-abundance.py {input.mismatches} > {output.mismatches} "

#Since the matches were subsampled, they need to be normalized before calculating fractions to make them equivalent to the mismatches which were not subsampled
#This is a rough estimate, since it calculates the total fraction subsampled across the whole dataset
rule normalize_match_counts_by_total_seqs:
	input:
		path="intermediate/compute-workflow/08-checked/{primer}/{mismatches}/",
		countTableHits="intermediate/classify-workflow/05-tax-counts/{study}.{group}.{primer}.{mismatches}.hits.all.order.counts.tsv",
		countTableNoHits="intermediate/classify-workflow/05-tax-counts/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.tsv"
	output:
		hits="intermediate/classify-workflow/07-normalized-counts/{study}.{group}.{primer}.{mismatches}.hits.all.order.counts.normalized.tsv",
		nohits="intermediate/classify-workflow/07-normalized-counts/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.normalized.tsv"
	shell:
		"scripts/mismatch-characterization/normalize_match_counts_by_total_seqs.sh {input.path} "
		"{wildcards.group}.{wildcards.primer}.{wildcards.mismatches} " #Pattern for matching, not an input file
		"{wildcards.group}.{wildcards.primer}.{wildcards.mismatches}.sub5k.hit.filtered.fastq " #Pattern for matching, not an input file
		"{input.countTableHits} > {output.hits} ; "
		"scripts/mismatch-characterization/normalize_match_counts_by_total_seqs.sh {input.path} "
		"{wildcards.group}.{wildcards.primer}.{wildcards.mismatches} " #Pattern for matching, not an input file
		"{wildcards.group}.{wildcards.primer}.{wildcards.mismatches}.sub5k.hit.filtered.fastq " #Pattern for matching, not an input file
		"{input.countTableNoHits} > {output.nohits}"

#For a given mismatch threshold, counts the total number of sequences of both matches and classify_mismatches
#Used to calculate the quantitative importance of each mismatch in terms of the whole dataset
rule count_total_filtered_hits:
	input:
		"intermediate/compute-workflow/08-checked/"
	output:
		"intermediate/classify-workflow/08-total-filtered-seqs/{study}.{group}.{primer}.{mismatches}.totalFilteredSeqs.tsv"
	shell:
		"totalFilteredSeqs=`cat {input}/{wildcards.primer}/{wildcards.mismatches}/*{wildcards.group}*filtered.fastq | grep -c \"^@\"` || totalFilteredSeqs=0 ; "
		"printf \"{wildcards.primer}.{wildcards.group}.{wildcards.mismatches}\t$totalFilteredSeqs\n\" >> {output}"

#identify target taxonomies to quantify; choose only the abundant things
rule generate_target_files:
	input:
		expand("intermediate/classify-workflow/06-tax-counts-fractions/{study}.{group}.{primer}.{mismatches}.hits.all.order.counts.frac.tsv", study=config["study"], group=config["groups"], primer=config["primer"], mismatches=["0-mismatch", "1-mismatch", "2-mismatch"]),
		expand("intermediate/classify-workflow/06-tax-counts-fractions/{study}.{group}.{primer}.0-mismatch.nohits.all.order.counts.frac.tsv", study=config["study"], group=config["groups"], primer=config["primer"], mismatches=["0-mismatch", "1-mismatch", "2-mismatch"])
	output:
		"intermediate/classify-workflow/09-target-taxa/{study}.{group}.{primer}.targets"
	shell:
		"cat {input} | cut -f1 | sort | uniq > {output}"

#Make an overall taxonomic summary for each group identified in target files
#Summary shows the proportion of a taxon that is missed by the primer, as well as what fraction of the total dataset these taxon mismatches represent
rule compute_frac_mismatched:
	input:
		targets="intermediate/classify-workflow/09-target-taxa/{study}.{group}.{primer}.targets", #target taxonomies
		totalHits="intermediate/classify-workflow/08-total-filtered-seqs/{study}.{group}.{primer}.{mismatches}.totalFilteredSeqs.tsv", #A count of total filtered hits
		normalized="intermediate/classify-workflow/07-normalized-counts/{study}.{group}.{primer}.{mismatches}.hits.all.order.counts.normalized.tsv", #Subsampled hits normalized by total
		counts="intermediate/classify-workflow/05-tax-counts/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.tsv" #mismatched hits
	output:
		"output/classify-workflow/overall-summaries/{study}.{group}.{primer}.{mismatches}.summary.tsv" #A summary file that tells how quantitatively significant the mismatches are for the group in question and for the whole dataset
	shell:
		"scripts/mismatch-characterization/compute-frac-mismatched.sh " #bash script that takes the 4 input arguments above
		"{input.targets} {input.totalHits} {input.counts} {input.normalized} > {output}"


#Concatenate info files across whole study
rule concatenate_info_files_mismatches:
	output:
		"intermediate/classify-workflow/10-concatenated-info-files/{study}.{group}.{primer}.{mismatches}.nohits.all.info"
	shell:
		"find intermediate/compute-workflow/08-checked/{wildcards.primer}/{wildcards.mismatches} -type f -name "
		"\"*{wildcards.group}*{wildcards.primer}*{wildcards.mismatches}.info\" -print0 | "
		"xargs -0 cat > {output}"

#Make summaries of mismatches
#CAVEAT: Note that the most relevant output is not the 0-mismatch, since cutadapt does not recognize the primer region in the output if there are many mismatches vs. the cutoff. These more distant mismatches will be regarded as a negative result (-1 in the info file), and therefore some relevant mismatches will be omitted from the summary.
#So when you look at the output, the variants of the 0-mismatch will only be the things that match at 0-mismatch
#Script requires pandas
rule make_mismatch_alignments_and_summarize:
	output:
		summary="output/classify-workflow/overall-summaries/{study}.{group}.{primer}.{mismatches}.aln.summary.tsv",
		mismatchAlignment="output/classify-workflow/overall-summaries/{study}.{group}.{primer}.{mismatches}.aln.fasta"
	input:
		"intermediate/classify-workflow/10-concatenated-info-files/{study}.{group}.{primer}.{mismatches}.nohits.all.info"
	conda:
		"envs/biopython.yaml"
	shell:
		"scripts/make-mismatch-alignments-and-summarize.py --info {input} "
		"--summaryout {output.summary} --alignmentout {output.mismatchAlignment}"

rule summarize_mismatch_info:
	input:
		"output/classify-workflow/overall-summaries/{study}.{group}.{primer}.0-mismatch.summary.tsv",
		"output/classify-workflow/overall-summaries/{study}.{group}.{primer}.1-mismatch.summary.tsv",
		"output/classify-workflow/overall-summaries/{study}.{group}.{primer}.2-mismatch.summary.tsv"
	output:
		pastedSummaries=temp("output/classify-workflow/{study}.{group}.{primer}.0-2mm.pasted.tsv"),
		comparisonOutput="output/classify-workflow/overall-summaries/{study}.{group}.{primer}.taxonFracMismatched.0-2mm.tsv"
	shell:
		"paste {input} > {output.pastedSummaries} ; "
		"scripts/mismatch-characterization/summarize-taxa-mismatches.py {output.pastedSummaries} > {output.comparisonOutput}"

#greater than 50 total mismatches AND > 5% of total sequences mismatched for that taxon
rule get_taxa_with_many_mismatches:
	input:
		"output/classify-workflow/overall-summaries/{study}.{group}.{primer}.0-mismatch.summary.tsv"
	output:
		"intermediate/classify-workflow/11-taxa-with-many-mismatches/{study}.{group}.{primer}.0-mismatch.gt50mm-and-5pcmismatched-taxa.txt"	
	shell:
		"./scripts/printTaxaIfGt5pcMismatch100obs.py {input} > {output}"	

rule get_fastq_headers_for_mismatched_taxa:
	input:
		taxa="intermediate/classify-workflow/11-taxa-with-many-mismatches/{study}.{group}.{primer}.0-mismatch.gt50mm-and-5pcmismatched-taxa.txt",
		taxtable="intermediate/classify-workflow/05-tax-table/{study}.{group}.{primer}.0-mismatch.nohits.all.order.counts.taxtable"
	output:
		"intermediate/classify-workflow/11-taxa-with-many-mismatches/{study}.{group}.{primer}.0-mismatch.gt50mm-and-5pcmismatched-taxa.headers.tsv"
	shell:
		"grep -f {input.taxa} {input.taxtable} | cut -f1-2 > {output} || touch {output}"

rule generate_barplot_input:
	input:
                hits="intermediate/classify-workflow/07-normalized-counts/{study}.{group}.{primer}.{mismatches}.hits.all.order.counts.normalized.tsv",
                nohits="intermediate/classify-workflow/07-normalized-counts/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.normalized.tsv"		
	output:
                "output/classify-workflow/plots/{study}.{group}.{primer}.{mismatches}.barplot-input.tsv"
	shell:
                "scripts/mismatch-characterization/generate-barplot-input.py {input.hits} {input.nohits} {output} "

rule make_tax_matchVSmismatch_barplots:
	input:
		"output/classify-workflow/plots/{study}.{group}.{primer}.{mismatches}.barplot-input.tsv"
	params:
		"{study}.{group}.{primer}"
	conda:
		"envs/ggplot2.yaml"
	output:
		"output/classify-workflow/plots/matchVSmismatch-barplots/{study}.{group}.{primer}.taxonFracMismatched.{mismatches}.pdf"
	script:
		"scripts/mismatch-characterization/make-taxa-barplots-match-vs-mismatch.R"

#if a group has more than 10 reads in the summary, and represents at least 1% of the total mismatches, then print it out for further consideration
rule filter_summary_taxa:
        input:
                "output/classify-workflow/overall-summaries/{study}.{group}.{primer}.0-mismatch.summary.tsv"
	output:
		"output/classify-workflow/filtered-0-mismatches/{study}.{group}.{primer}.0-mismatch.gt1pc.gt10obs.tsv"
	shell:
		"./scripts/printIfGt1pc.py {input} | sort -r -k2 > {output}"

rule calc_primer_pair_mismatch_overlap:
	input:
		fwdprimer=lambda wildcards: "output/classify-workflow/filtered-0-mismatches/" + config["study"] + "." + config["groups"][wildcards.group] + "." + config["primer_pairs"][wildcards.primer_pair][0] + ".0-mismatch.gt1pc.gt10obs.tsv",
                revprimer=lambda wildcards: "output/classify-workflow/filtered-0-mismatches/" + config["study"] + "." + config["groups"][wildcards.group] + "." + config["primer_pairs"][wildcards.primer_pair][1] + ".0-mismatch.gt1pc.gt10obs.tsv"
	output:
		"output/classify-workflow/summary-mismatch-overlap-primer-pairs/{study}.{group}.{primer_pair}.avgCase.tsv"
	shell:
		"./scripts/primer-pair-subtract-overlap.py {input.fwdprimer} {input.revprimer} > {output}"

#assumes the compute pipeline is completed
rule paste_summaries:
	params:
		fwdprimer=lambda wildcards: config["primer_pairs"][wildcards.primer_pair][0],
                revprimer=lambda wildcards: config["primer_pairs"][wildcards.primer_pair][1]
	output:
                "output/classify-workflow/pasted-summaries/{study}.{group}.{primer_pair}.pasted.tsv"
	shell:
		"tmpfwdprimer=`mktemp /tmp/fwdprimer.summary.sorted.XXXXXXXXXXXXXXXX` ; "
		"tmprevprimer=`mktemp /tmp/revprimer.summary.sorted.XXXXXXXXXXXXXXXX` ; "
		"find ./output/compute-workflow/09-summary/ -type f -name \"*{wildcards.group}.{params.fwdprimer}.0-mismatch.summary.tsv\" -print0 | xargs -0 cat |  sort -t$'\\t' -k1,1 -k2,2 > $tmpfwdprimer ; "
		"find ./output/compute-workflow/09-summary/ -type f -name \"*{wildcards.group}.{params.revprimer}.0-mismatch.summary.tsv\" -print0 | xargs -0 cat |  sort -t$'\\t' -k1,1 -k2,2 > $tmprevprimer ; "
		"paste $tmpfwdprimer $tmprevprimer > {output} ; "
		"rm $tmpfwdprimer $tmprevprimer"
		

rule normalize_summaries:
	input:
                pasted="output/classify-workflow/pasted-summaries/{study}.{group}.{primer_pair}.pasted.tsv",
		normFactor="output/classify-workflow/summary-mismatch-overlap-primer-pairs/{study}.{group}.{primer_pair}.avgCase.tsv"
	output:
		"output/classify-workflow/normalized-summaries/{study}.{group}.{primer_pair}.normalized.tsv"
	shell:
		"./scripts/normalize-for-master-figure.py {input.pasted} `cat {input.normFactor}` {wildcards.primer_pair} > {output}"

