CUTOFF = config["cutoff"]

rule all:
	input:
		#expand("15-matches-classified/individual/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.hit.filtered.VSEARCHsintax-SILVA132.tax", sample=config["samples"], group=config["groups"], primer=config["primer"], mismatches=config["mismatches"], direction=['fwd','rev'])
		#expand("13-classified/individual/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.nohit.filtered.VSEARCHsintax-SILVA132.tax", sample=config["samples"], group=config["groups"], primer=config["primer"], mismatches=config["mismatches"], direction=['fwd','rev'])
		#expand("intermediate/{study}.{group}.{primer}.0-mismatch.hits.all.order.counts.normalized.tsv", study=config["study"], sample=config["samples"], group=config["groups"], primer=config["primer"])
		#expand("intermediate/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.minAbund-0.01.tsv", study=config["study"], sample=config["samples"], group=config["groups"], primer=config["primer"], mismatches=["0-mismatch", "1-mismatch", "2-mismatch"])
		#"intermediate/totalFilteredHits.tsv"
		#expand("intermediate/{study}.{group}.{primer}.targets", study=config["study"], group=config["groups"], primer=config["primer"]),
		#expand("intermediate/{study}.{group}.{primer}.totalFilteredHits.tsv", study=config["study"], group=config["groups"], primer=config["primer"])
		#expand("intermediate/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.taxtable", study=config["study"], group=config["groups"], primer=config["primer"], mismatches=["0-mismatch", "1-mismatch", "2-mismatch"])
		#expand("intermediate/{study}.{group}.{primer}.{mismatches}.nohits.all.info", study=config["study"], group=config["groups"], primer=config["primer"], mismatches=["0-mismatch", "1-mismatch", "2-mismatch"])
		#expand("output/{study}.{group}.{primer}.{mismatches}.aln.summary.tsv", study=config["study"], group=config["groups"], primer=config["primer"], mismatches=["0-mismatch", "1-mismatch", "2-mismatch"])
		#expand("compareTags2MG/01-concatenated/{sample}.{group}.concat.515Y-926R.fasta", sample=config["samples"], group=config["groups"])
		expand("compareTags2MG/05-MG-blasted-against-ASVs/{sample}.PROK.nonzero.ASV.blastout.tsv", sample=config["samples"])

#Generate concatenated alignments that will be used for subsetting
rule concat_fwd_and_reverse_alignments:
	input:
		fwd="05-pyNAST-aligned/{sample}.fwd.{group}/{sample}.fwd.SSU.{group}_pynast_aligned.fasta",
		rev="05-pyNAST-aligned/{sample}.rev.{group}/{sample}.rev.SSU.{group}_pynast_aligned.fasta"
	output:
		temp("compareTags2MG/00-concatenated/{sample}.{group}.concat.fasta")
	shell:
		"cat {input.fwd} {input.rev} > {output}"

rule subset_to_primer_region:
	input:
		"compareTags2MG/00-concatenated/{sample}.{group}.concat.fasta"
	output:
		"compareTags2MG/01-concatenated/{sample}.{group}.concat.515Y-926R.fasta"
	params:
		start=lambda wildcards: config["primerROI"][wildcards.group]["515Y"][0],
		end=lambda wildcards: config["primerROI"][wildcards.group]["926R"][1]
	shell:
		"scripts/filter-pyNAST-for-ROI-v2.py --input {input} --start {params.start} --end {params.end} --padding 0 --output {output}  ; " #subset
		"sed -i 's/-//g' {output}" #remove gaps that may have been introduced by pyNAST

#concatenate the 3 different PROK categories
rule merge_PROK:
	input:
		"compareTags2MG/01-concatenated/{sample}.ARCH.concat.515Y-926R.fasta",
		"compareTags2MG/01-concatenated/{sample}.BACT-NON-CYANO.concat.515Y-926R.fasta",
		"compareTags2MG/01-concatenated/{sample}.BACT-CYANO.concat.515Y-926R.fasta"
	output:
		"compareTags2MG/01-concatenated/{sample}.PROK.concat.515Y-926R.fasta"
	shell:
		"cat {input} > {output}"

#get non-zero abundance ASVs from table, and output the ids
rule parse_16S_ASV_table:
	input:
		"GP13-sample-SRA.tsv",
		"191004_ASV_info/PROKs/191004_GP13_all-16S-seqs.with-tax.proportions.tsv"
	params:
		"{sample}"
	output:
		"compareTags2MG/02-ASV-ids/{sample}.PROK.nonzero.ASV.ids"
	script:
		"scripts/make-db-of-non-zero-abund-ASV.py"

rule get_sample_nonzero_16S_ASV_fastas:
	input:
		fasta="191004_ASV_info/PROKs/191004_GP13_all-16S-seqs.dna-sequences.fasta",
		ids="compareTags2MG/02-ASV-ids/{sample}.PROK.nonzero.ASV.ids"
	output:
		"compareTags2MG/03-ASV-fastas/{sample}.PROK.nonzero.ASV.fasta"
	shell:
		"seqtk subseq {input.fasta} {input.ids} > {output}"

rule make_blast_dbs_16S:
	input:
		"compareTags2MG/03-ASV-fastas/{sample}.PROK.nonzero.ASV.fasta"
	output:
		expand("compareTags2MG/04-ASV-blastdbs/{{sample}}.PROK.nonzero.ASV.db.{ext}", ext=["nhr", "nin", "nsq"])
	params:
		filestem="compareTags2MG/04-ASV-blastdbs/{sample}.PROK.nonzero.ASV.db"
	conda:
		"envs/qiime1.yaml"
	shell:
		"makeblastdb -in {input} -dbtype nucl -out {params.filestem} ; touch {output}"

rule blast_MG_vs_tags:
	input:
		database_files=lambda wildcards: expand("compareTags2MG/04-ASV-blastdbs/{sample}.PROK.nonzero.ASV.db.{ext}", ext=["nhr", "nin", "nsq"], sample=wildcards.sample),
		query="compareTags2MG/01-concatenated/{sample}.PROK.concat.515Y-926R.fasta"
	output:
		"compareTags2MG/05-MG-blasted-against-ASVs/{sample}.PROK.nonzero.ASV.blastout.tsv"
	params:
		dbname="compareTags2MG/04-ASV-blastdbs/{sample}.PROK.nonzero.ASV.db"
	conda:
		"envs/qiime1.yaml"
	shell:
		"blastn -qcov_hsp_perc 100 -perc_identity 100 -outfmt 6 -query {input.query} -db {params.dbname} > {output}"

"""
rule get_ASV_subsets:

rule BLAST_MG_SSU_rRNA_against_ASVs:

rule sift_BLAST_results:

#get proportion
rule count_positive_BLAST_hits:

rule plot_proportions:

rule plot_ASV_vs_BLAST_results:

#The metagenomic SSU rRNA that does NOT match the ASVs
rule classify_negative_BLAST_hits:


Need to compare sample to sample. To do:
1. Make lookup table between SRAid and sampleID
2. Get ASV subset for EACH sample. Take only non-zero abundance ASVs based on ASV table.
"""

"""
rule grab_full_fastas:
	input:
		fwdFQ="10-checked/{primer}/{mismatches}/{sample}.SSU.fwd.{group}.{primer}.{mismatches}.nohit.filtered.fastq",
		revFQ="10-checked/{primer}/{mismatches}/{sample}.SSU.rev.{group}.{primer}.{mismatches}.nohit.filtered.fastq",
		fwdFA="02-graftm_sifted/{sample}.fwd.SSU.hits.fa",
		revFA="02-graftm_sifted/{sample}.rev.SSU.hits.fa",
	output:
		fwd="12-full-fastas/{sample}.SSU.fwd.{group}.{primer}.{mismatches}.nohit.filtered.fasta",
		rev="12-full-fastas/{sample}.SSU.rev.{group}.{primer}.{mismatches}.nohit.filtered.fasta"
	conda:
		"envs/bbmap.yaml"
	shell:
		"filterbyname.sh names={input.fwdFQ} include=t in={input.fwdFA} out={output.fwd} ; "
		"filterbyname.sh names={input.revFQ} include=t in={input.revFA} out={output.rev} "

#Concatenating 0-mismatch misses because 1-mismatch and 2-mismatch are a subset of the 0-mismatch misses, and can be parsed out later
rule concatenate_mismatched_fastas:
#Note: concatenating to speed up classification step; using xargs for large numbers of files (otherwise bash will complain argument list too long)
	output:
		"tmp.mismatches.concatenated.fasta"
	shell:
		"find ./12-full-fastas -type f -name \"*0-mismatch*.fasta\" -print0 | xargs -0 cat > {output}"

rule classify_mismatches:
	input:
		"tmp.mismatches.concatenated.fasta"
	output:
		"13-classified/all/all.mismatches.VSEARCH.classified.tsv"
	threads:
		20
	conda:
		"envs/vsearch.yaml"
	shell:
		"vsearch --sintax {input} "
		"--db /home/db/VSEARCH/silva132_99_sintax.udb "
		"--tabbedout {output} --threads 40 --sintax_cutoff 0 --top_hits_only --topn 1 --notrunclabels"

rule deconcat_classifications:
	input:
		fasta="12-full-fastas/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.nohit.filtered.fasta",
		concatenated="13-classified/all/all.mismatches.VSEARCH.classified.tsv"
	output:
		"13-classified/individual/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.nohit.filtered.VSEARCHsintax-SILVA132.tax"
	shell:
		"tmpfile=`mktemp /tmp/fastq-ids.XXXXXXXXXXXXXXXX` ; seqmagick extract-ids {input.fasta} >> $tmpfile ; "
		"grep -f $tmpfile {input.concatenated} > {output} || touch {output} "


rule subsample_matched_fastas:
	input:
		fwdFQ="10-checked/{primer}/{mismatches}/{sample}.SSU.fwd.{group}.{primer}.{mismatches}.hit.filtered.fastq",
		revFQ="10-checked/{primer}/{mismatches}/{sample}.SSU.rev.{group}.{primer}.{mismatches}.hit.filtered.fastq",
		fwdFA="02-graftm_sifted/{sample}.fwd.SSU.hits.fa",
		revFA="02-graftm_sifted/{sample}.rev.SSU.hits.fa",
	output:
		fwd="14-subsampled-matched-fastas/{sample}.SSU.fwd.{group}.{primer}.{mismatches}.hit.filtered.fasta",
		rev="14-subsampled-matched-fastas/{sample}.SSU.rev.{group}.{primer}.{mismatches}.hit.filtered.fasta"
	conda:
		"envs/bbmap.yaml"
	shell:
		"filterbyname.sh names={input.fwdFQ} reads=5000 include=t in={input.fwdFA} out={output.fwd} ; "
		"filterbyname.sh names={input.revFQ} reads=5000 include=t in={input.revFA} out={output.rev} "

#Concatenating 2-mismatch hits because 1-mismatch and 0-mismatch are a subset of the 2-mismatch hits, and can be parsed out later
rule concatenate_matched_fastas:
	output:
		"tmp.matches.subsampled.concatenated.fasta"
	shell:
		"find ./14-subsampled-matched-fastas -type f -name \"*2-mismatch*.fasta\" -print0 | xargs -0 cat > {output}"


rule classify_matches_subsample:
	input:
		"tmp.matches.subsampled.concatenated.fasta"
	output:
		"15-matches-classified/all/all.matches.filtered.VSEARCH.classified.tsv"
	threads:
		20
	conda:
		"envs/vsearch.yaml"
	shell:
		"vsearch --sintax {input} "
                "--db /home/db/VSEARCH/silva132_99_sintax.udb "
                "--tabbedout {output} --threads 40 --sintax_cutoff 0 --top_hits_only --topn 1 --notrunclabels"

rule deconcat_matches_classifications:
	input:
		fasta="14-subsampled-matched-fastas/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.hit.filtered.fasta",
		concatenated="15-matches-classified/all/all.matches.filtered.VSEARCH.classified.tsv"
	output:
		"15-matches-classified/individual/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.hit.filtered.VSEARCHsintax-SILVA132.tax"
	shell:
		"tmpfile=`mktemp /tmp/fastq-ids.XXXXXXXXXXXXXXXX` ; seqmagick extract-ids {input.fasta} >> $tmpfile ; "
		"grep -f $tmpfile {input.concatenated} > {output} || touch {output} ; "


#Concatenate taxonomy files
rule cat_tax_for_all_samples_matches_and_mismatches:
	output:
		mismatches="intermediate/{study}.{group}.{primer}.{mismatches}.nohits.all.tax",
		matches="intermediate/{study}.{group}.{primer}.{mismatches}.hits.all.tax"
	shell:
		"find 13-classified/individual/ -type f -name "
		"\"*{wildcards.group}*{wildcards.primer}*{wildcards.mismatches}*tax\" -print0 | "
		"xargs -0 cat > {output.mismatches} ; "
		"find 15-matches-classified/individual/ -type f -name "
		"\"*{wildcards.group}*{wildcards.primer}*{wildcards.mismatches}*tax\" -print0 | "
		"xargs -0 cat > {output.matches}"

#Counting order-level groupings (can adjust level with the "cut -d, -f1-4" parameter below)
rule count_tax_matches_and_mismatches:
	input:
		matches="intermediate/{study}.{group}.{primer}.{mismatches}.hits.all.tax",
		mismatches="intermediate/{study}.{group}.{primer}.{mismatches}.nohits.all.tax"
	output:
		matches="intermediate/{study}.{group}.{primer}.{mismatches}.hits.all.order.counts.tsv",
		mismatches="intermediate/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.tsv",
		taxTableMatches="intermediate/{study}.{group}.{primer}.{mismatches}.hits.all.order.counts.taxtable",
		taxTableMismatches="intermediate/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.taxtable"
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
		matches="intermediate/{study}.{group}.{primer}.{mismatches}.hits.all.order.counts.tsv",
		mismatches="intermediate/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.tsv"
	params:
		minAbund = CUTOFF  #Change fractional value in config file if desired, default 0.01
	output:
		matches="intermediate/{study}.{group}.{primer}.{mismatches}.hits.all.order.counts.min" + str(CUTOFF) + ".tsv",
		mismatches="intermediate/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.min" + str(CUTOFF) + ".tsv"
	shell:
		"scripts/mismatch-characterization/filter-by-fractional-abundance.py {input.matches} {params.minAbund} > {output.matches} ; "
		"scripts/mismatch-characterization/filter-by-fractional-abundance.py {input.mismatches} {params.minAbund} > {output.mismatches} "

#Since the matches were subsampled, they need to be normalized before calculating fractions to make them equivalent to the mismatches which were not subsampled
#This is a rough estimate, since it calculates the total fraction subsampled across the whole dataset
rule normalize_match_counts_by_total_seqs:
	input:
		path="10-checked/{primer}/{mismatches}/",
		countTable="intermediate/{study}.{group}.{primer}.{mismatches}.hits.all.order.counts.tsv"
	output:
		"intermediate/{study}.{group}.{primer}.{mismatches}.hits.all.order.counts.normalized.tsv"
	shell:
		"scripts/mismatch-characterization/normalize_match_counts_by_total_seqs.sh {input.path} "
		"{wildcards.group}.{wildcards.primer}.{wildcards.mismatches} " #Pattern for matching, not an input file
		"{wildcards.group}.{wildcards.primer}.{wildcards.mismatches}.hit.filtered.fasta " #Pattern for matching, not an input file
		"{input.countTable} > {output}"

#For a given mismatch threshold, counts the total number of sequences of both matches and classify_mismatches
#Used to calculate the quantitative importance of each mismatch in terms of the whole dataset
rule count_total_filtered_hits:
	input:
		"10-checked"
	output:
		"intermediate/{study}.{group}.{primer}.{mismatches}.totalFilteredSeqs.tsv"
	shell:
		"totalFilteredSeqs=`cat {input}/{wildcards.primer}/{wildcards.mismatches}/SRR*{wildcards.group}*filtered.fastq | grep -c \"^@\"` || totalFilteredSeqs=0 ; "
		"printf \"{wildcards.primer}.{wildcards.group}.{wildcards.mismatches}\t$totalFilteredSeqs\n\" >> {output}"

#identify target taxonomies to quantify; choose only the abundant things
rule generate_target_files:
	input:
		"intermediate/{study}.{group}.{primer}.{mismatches}.hits.all.order.counts.min" + str(CUTOFF) + ".tsv",
		"intermediate/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.min" + str(CUTOFF) + ".tsv"
	output:
		"intermediate/{study}.{group}.{primer}.{mismatches}.targets"
	shell:
		"cat {input} | cut -f1 | sort | uniq > {output}"

#Make an overall taxonomic summary for each group identified in target files
#Summary shows the proportion of a taxon that is missed by the primer, as well as what fraction of the total dataset these taxon mismatches represent
rule compute_frac_mismatched:
	input:
		targets="intermediate/{study}.{group}.{primer}.{mismatches}.targets", #target taxonomies
		totalHits="intermediate/{study}.{group}.{primer}.{mismatches}.totalFilteredSeqs.tsv", #A count of total filtered hits
		normalized="intermediate/{study}.{group}.{primer}.{mismatches}.hits.all.order.counts.normalized.tsv", #Subsampled hits normalized by total
		counts="intermediate/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.tsv" #mismatched hits
	output:
		"output/{study}.{group}.{primer}.{mismatches}.summary.tsv" #A summary file that tells how quantitatively significant the mismatches are for the group in question and for the whole dataset
	shell:
		"scripts/mismatch-characterization/compute-frac-mismatched.sh " #bash script that takes the 4 input arguments above
		"{input.targets} {input.totalHits} {input.counts} {input.normalized} > {output}"

#Concatenate info files across whole study
rule concatenate_info_files_mismatches:
	output:
		"intermediate/{study}.{group}.{primer}.{mismatches}.nohits.all.info"
	shell:
		"find 10-checked/{wildcards.primer}/{wildcards.mismatches} -type f -name "
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


#Get ids for each taxon identified in target files
rule get_ids_for_target_files:
	input:
		targets="intermediate/{study}.{group}.{primer}.{mismatches}.targets", #target taxonomies
		taxtable="intermediate/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.taxtable" #A table with fastq headers and taxonomic assignments
	output:

"""
