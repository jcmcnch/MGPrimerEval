CUTOFF = config["cutoff"]

rule all:
	input:
		#expand("classify-workflow-intermediate/01-mismatches-classified/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.nohit.filtered.VSEARCHsintax-SILVA132.tax", sample=config["samples"], study=config["study"], group=config["groups"], primer=config["primer"], mismatches=["0-mismatch", "1-mismatch", "2-mismatch"], direction=['fwd','rev']),
		#expand("classify-workflow-intermediate/03-matches-classified/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.sub5k.hit.filtered.VSEARCHsintax-SILVA132.tax", sample=config["samples"], study=config["study"], group=config["groups"], primer=config["primer"], mismatches=["0-mismatch", "1-mismatch", "2-mismatch"], direction=['fwd','rev']),
		#BELOW STEP NEEDS TO BE DEBUGGED - normalization does not seem to be working
		expand("output-classify-workflow/{study}.{group}.{primer}.{mismatches}.summary.tsv", sample=config["samples"], study=config["study"], group=config["groups"], primer=config["primer"], mismatches=["0-mismatch", "1-mismatch", "2-mismatch"], direction=['fwd','rev'])

rule classify_mismatches:
	input:
		"compute-workflow-intermediate/08-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.nohit.filtered.fastq",
	output:
		"classify-workflow-intermediate/01-mismatches-classified/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.nohit.filtered.VSEARCHsintax-SILVA132.tax"
	threads:
		20
	conda:
		"envs/vsearch.yaml"
	shell:
		#Double pipe is OR operator and will only be executed if vsearch returns an error. Necessary otherwise empty files will cause vsearch to fail.
		"""
		vsearch --sintax {input} \
		--db /home/db/VSEARCH/silva132_99_sintax.udb \
		--tabbedout {output} --threads 40 --sintax_cutoff 0 --top_hits_only --topn 1 --notrunclabels \
		|| touch {output}
		"""

#Take up to 5000 reads from the matched files
rule subsample_matched_fastas:
	input:
		"compute-workflow-intermediate/08-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.hit.filtered.fastq",
	output:
		"classify-workflow-intermediate/02-matches-subsampled/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.sub5k.hit.filtered.fastq"
	conda:
		"envs/bbmap.yaml"
	shell:
		"reformat.sh samplereadstarget=5000 in={input} out={output} "

rule classify_matches_subsample:
	input:
		"classify-workflow-intermediate/02-matches-subsampled/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.sub5k.hit.filtered.fastq"
	output:
		"classify-workflow-intermediate/03-matches-classified/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.sub5k.hit.filtered.VSEARCHsintax-SILVA132.tax"
	threads:
		20
	conda:
		"envs/vsearch.yaml"
	shell:
		#Double pipe is OR operator and will only be executed if vsearch returns an error. Necessary otherwise empty files will cause vsearch to fail.
		"""
		vsearch --sintax {input} \
		--db /home/db/VSEARCH/silva132_99_sintax.udb \
		--tabbedout {output} --threads 40 --sintax_cutoff 0 --top_hits_only --topn 1 --notrunclabels \
		|| touch {output}
		"""

"""
Rules below comprise a workflow for generating summaries of
which taxa are most discriminated against by a particular primer set.
Implemented using common bash tools and tested on Ubuntu 16.04, not tested on other systems.
Results are not per-sample, but rather across an entire dataset.
"""

#Concatenate taxonomy files
rule cat_tax_for_all_samples_matches_and_mismatches:
	output:
		mismatches="classify-workflow-intermediate/04-tax-concatenated/{study}.{group}.{primer}.{mismatches}.nohits.all.tax",
		matches="classify-workflow-intermediate/04-tax-concatenated/{study}.{group}.{primer}.{mismatches}.hits.all.tax"
	shell:
		"find classify-workflow-intermediate/01-mismatches-classified/ -type f -name "
		"\"*{wildcards.group}*{wildcards.primer}*{wildcards.mismatches}*tax\" -print0 | "
		"xargs -0 cat > {output.mismatches} ; "
		"find classify-workflow-intermediate/03-matches-classified/ -type f -name "
		"\"*{wildcards.group}*{wildcards.primer}*{wildcards.mismatches}*tax\" -print0 | "
		"xargs -0 cat > {output.matches}"

#Counting order-level groupings (can adjust level with the "cut -d, -f1-4" parameter below)
rule count_tax_matches_and_mismatches:
	input:
		mismatches="classify-workflow-intermediate/04-tax-concatenated/{study}.{group}.{primer}.{mismatches}.nohits.all.tax",
		matches="classify-workflow-intermediate/04-tax-concatenated/{study}.{group}.{primer}.{mismatches}.hits.all.tax"
	output:
		matches="classify-workflow-intermediate/05-tax-counts/{study}.{group}.{primer}.{mismatches}.hits.all.order.counts.tsv",
		mismatches="classify-workflow-intermediate/05-tax-counts/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.tsv",
		taxTableMatches="classify-workflow-intermediate/05-tax-table/{study}.{group}.{primer}.{mismatches}.hits.all.order.counts.taxtable",
		taxTableMismatches="classify-workflow-intermediate/05-tax-table/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.taxtable"
	shell:
		"sed -re 's/\([0-9]{{1}}\.[0-9]{{2}}\)//g' {input.matches} | tee {output.taxTableMatches} |" #Remove confidence estimations from VSEARCH output, keep a copy for later steps but also pipe to subsequent commands
		"cut -f2 | sort | cut -d, -f1-4 | sort | uniq -c | " #Take only tax column, collapse to order level, then count unique occurrences
		"tail -f -n +2 | awk '{{print $1,\"\t\",$2}}' > {output.matches} ; " #Process output into tsv format to stdout
		"sed -re 's/\([0-9]{{1}}\.[0-9]{{2}}\)//g' {input.mismatches} | tee {output.taxTableMismatches} |" #Remove confidence estimations from VSEARCH output, keep a copy for later steps but also pipe to subsequent commands
		"cut -f2 | sort | cut -d, -f1-4 | sort | uniq -c | " #Take only tax column, collapse to order level, then count unique occurrences
		"tail -f -n +2 | awk '{{print $1,\"\t\",$2}}' > {output.mismatches}" #Process output into tsv format to stdout

#Now take only those with greater than 1 % abundance (among mismatches) using basic python script (can change abundance cutoff if you desire)
rule filter_tax_matches_by_abundance:
	input:
		matches="classify-workflow-intermediate/05-tax-counts/{study}.{group}.{primer}.{mismatches}.hits.all.order.counts.tsv",
		mismatches="classify-workflow-intermediate/05-tax-counts/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.tsv"
	params:
		minAbund = CUTOFF  #Change fractional value in config file if desired, default 0.01
	output:
		matches="classify-workflow-intermediate/06-tax-counts-abund-filtered/{study}.{group}.{primer}.{mismatches}.hits.all.order.counts.min" + str(CUTOFF) + ".tsv",
		mismatches="classify-workflow-intermediate/06-tax-counts-abund-filtered/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.min" + str(CUTOFF) + ".tsv"
	shell:
		"scripts/mismatch-characterization/filter-by-fractional-abundance.py {input.matches} {params.minAbund} > {output.matches} ; "
		"scripts/mismatch-characterization/filter-by-fractional-abundance.py {input.mismatches} {params.minAbund} > {output.mismatches} "

#Since the matches were subsampled, they need to be normalized before calculating fractions to make them equivalent to the mismatches which were not subsampled
#This is a rough estimate, since it calculates the total fraction subsampled across the whole dataset
rule normalize_match_counts_by_total_seqs:
	input:
		path="compute-workflow-intermediate/08-checked/{primer}/{mismatches}/",
		countTable="classify-workflow-intermediate/05-tax-counts/{study}.{group}.{primer}.{mismatches}.hits.all.order.counts.tsv"
	output:
		"classify-workflow-intermediate/07-normalized-counts/{study}.{group}.{primer}.{mismatches}.hits.all.order.counts.normalized.tsv"
	shell:
		"scripts/mismatch-characterization/normalize_match_counts_by_total_seqs.sh {input.path} "
		"{wildcards.group}.{wildcards.primer}.{wildcards.mismatches} " #Pattern for matching, not an input file
		"{wildcards.group}.{wildcards.primer}.{wildcards.mismatches}.sub5k.hit.filtered.fastq " #Pattern for matching, not an input file
		"{input.countTable} > {output}"

#For a given mismatch threshold, counts the total number of sequences of both matches and classify_mismatches
#Used to calculate the quantitative importance of each mismatch in terms of the whole dataset
rule count_total_filtered_hits:
	input:
		"compute-workflow-intermediate/08-checked/"
	output:
		"classify-workflow-intermediate/08-total-filtered-seqs/{study}.{group}.{primer}.{mismatches}.totalFilteredSeqs.tsv"
	shell:
		"totalFilteredSeqs=`cat {input}/{wildcards.primer}/{wildcards.mismatches}/SRR*{wildcards.group}*filtered.fastq | grep -c \"^@\"` || totalFilteredSeqs=0 ; "
		"printf \"{wildcards.primer}.{wildcards.group}.{wildcards.mismatches}\t$totalFilteredSeqs\n\" >> {output}"

#identify target taxonomies to quantify; choose only the abundant things
rule generate_target_files:
	input:
		"classify-workflow-intermediate/06-tax-counts-abund-filtered/{study}.{group}.{primer}.{mismatches}.hits.all.order.counts.min" + str(CUTOFF) + ".tsv",
		"classify-workflow-intermediate/06-tax-counts-abund-filtered/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.min" + str(CUTOFF) + ".tsv"
	output:
		"classify-workflow-intermediate/09-target-taxa/{study}.{group}.{primer}.{mismatches}.targets"
	shell:
		"cat {input} | cut -f1 | sort | uniq > {output}"

#Make an overall taxonomic summary for each group identified in target files
#Summary shows the proportion of a taxon that is missed by the primer, as well as what fraction of the total dataset these taxon mismatches represent
rule compute_frac_mismatched:
	input:
		targets="classify-workflow-intermediate/09-target-taxa/{study}.{group}.{primer}.{mismatches}.targets", #target taxonomies
		totalHits="classify-workflow-intermediate/08-total-filtered-seqs/{study}.{group}.{primer}.{mismatches}.totalFilteredSeqs.tsv", #A count of total filtered hits
		normalized="classify-workflow-intermediate/07-normalized-counts/{study}.{group}.{primer}.{mismatches}.hits.all.order.counts.normalized.tsv", #Subsampled hits normalized by total
		counts="classify-workflow-intermediate/05-tax-counts/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.tsv" #mismatched hits
	output:
		"output-classify-workflow/{study}.{group}.{primer}.{mismatches}.summary.tsv" #A summary file that tells how quantitatively significant the mismatches are for the group in question and for the whole dataset
	shell:
		"scripts/mismatch-characterization/compute-frac-mismatched.sh " #bash script that takes the 4 input arguments above
		"{input.targets} {input.totalHits} {input.counts} {input.normalized} > {output}"

"""
#Concatenate info files across whole study
rule concatenate_info_files_mismatches:
	output:
		"classify-workflow-intermediate/10-concatenated-info-files/{study}.{group}.{primer}.{mismatches}.nohits.all.info"
	shell:
		"find compute-workflow-intermediate/08-checked/{wildcards.primer}/{wildcards.mismatches} -type f -name "
		"\"*{wildcards.group}*{wildcards.primer}*{wildcards.mismatches}.info\" -print0 | "
		"xargs -0 cat > {output}"

#Make summaries of mismatches
#CAVEAT: Note that the most relevant output is not the 0-mismatch, since cutadapt does not recognize the primer region in the output if there are many mismatches vs. the cutoff. These more distant mismatches will be regarded as a negative result (-1 in the info file), and therefore some relevant mismatches will be omitted from the summary.
rule make_mismatch_alignments_and_summarize:
	output:
		summary="output/{study}.{group}.{primer}.{mismatches}.aln.summary.tsv",
		mismatchAlignment="output/{study}.{group}.{primer}.{mismatches}.aln.fasta"
	input:
		"intermediate/{study}.{group}.{primer}.{mismatches}.nohits.all.info"
	conda:
		"envs/graftm.yaml"
	shell:
		"scripts/make-mismatch-alignments-and-summarize.py --info {input} "
		"--summaryout {output.summary} --alignmentout {output.mismatchAlignment}"
"""

"""
#Get ids for each taxon identified in target files
rule get_ids_for_target_files:
	input:
		targets="intermediate/{study}.{group}.{primer}.{mismatches}.targets", #target taxonomies
		taxtable="intermediate/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.taxtable" #A table with fastq headers and taxonomic assignments
	output:

"""
