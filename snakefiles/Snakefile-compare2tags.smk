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
		#"compareTags2MG/14-MG-not-in-ASVdb-classified-cat-parsed-fractions/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1.class.cat.counts.frac.min" + str(CUTOFF) + ".tsv"
		"compareTags2MG/15-MG-not-in-ASVdb-classified-kronas/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1andR2.krona.html"

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
		start=lambda wildcards: config["primerROI"][wildcards.group]["515Y"][1],
		end=lambda wildcards: config["primerROI"][wildcards.group]["926R"][0]
	shell:
		"scripts/filter-pyNAST-for-ROI-v2.py --input {input} --start {params.start} --end {params.end} --padding 0 --output {output}"

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
		"blastn -qcov_hsp_perc 100 -perc_identity 97 -outfmt 6 -query {input.query} -db {params.dbname} > {output}"

rule compare_MG_SSU_rRNA_with_ASVs:
	input:
		"GP13-sample-SRA.tsv",
		"191004_ASV_info/PROKs/191004_GP13_all-16S-seqs.with-tax.proportions.tsv",
		"compareTags2MG/05-MG-blasted-against-ASVs/{sample}.PROK.nonzero.ASV.blastout.tsv"
	params:
		"{sample}"
	conda:
		"envs/networkx.yaml"
	output:
		"compareTags2MG/06-MG-vs-ASV-tsv/{sample}.PROK.nonzero.ASV.comparison.tsv"
	script:
		"scripts/compare-MG-eASV-abund-from-blast-output.py"

rule plot_ASV_vs_BLAST_results:
	input:
		"compareTags2MG/06-MG-vs-ASV-tsv/{sample}.PROK.nonzero.ASV.comparison.tsv"
	params:
		"{sample}"
	output:
		"compareTags2MG/07-MG-vs-ASV-plots/{sample}.PROK.nonzero.ASV.comparison.svg"
	script:
		"scripts/seaborn-plot-correlations.py"

rule sift_unmatched:
	input:
		"compareTags2MG/05-MG-blasted-against-ASVs/{sample}.PROK.nonzero.ASV.blastout.tsv",
		"compareTags2MG/01-concatenated/{sample}.PROK.concat.515Y-926R.fasta"
	output:
		"compareTags2MG/08-MG-not-in-ASVdb/{sample}.PROK.515Y-926R.not-matching.ASVs.fasta"
	script:
		"scripts/get-non-matching.py"

#Align against E. coli for now, since just to get visual of locations
rule align_unmatched:
	input:
		seqs="compareTags2MG/08-MG-not-in-ASVdb/{sample}.PROK.515Y-926R.not-matching.ASVs.fasta",
		ref="SSU_refs/Ecoli_16s.fna"
	output:
		"compareTags2MG/09-MG-not-in-ASVdb-aligned/{sample}.PROK.515Y-926R.not-matching.ASVs.aligned.fasta"
	conda:
		"envs/pynast.yaml"
	shell:
		"pynast -p 10 -l 1 -i {input.seqs} -t {input.ref} -a {output}"

rule subset_non_matching:
	input:
		"compareTags2MG/09-MG-not-in-ASVdb-aligned/{sample}.PROK.515Y-926R.not-matching.ASVs.aligned.fasta"
	output:
		R1="compareTags2MG/10-MG-not-in-ASVdb-subsetted/{sample}.PROK.515Y-926R.not-matching.ASVs.aligned.R1.fasta",
		R2="compareTags2MG/10-MG-not-in-ASVdb-subsetted/{sample}.PROK.515Y-926R.not-matching.ASVs.aligned.R2.fasta",
		R1andR2="compareTags2MG/10-MG-not-in-ASVdb-subsetted/{sample}.PROK.515Y-926R.not-matching.ASVs.aligned.R1andR2.fasta"
	params:
		R1start=533,
		R1end=677,
		R2start=754,
		R2end=906,
		R1andR2start=533,
		R1andR2end=906
	shell:
		"scripts/filter-pyNAST-for-ROI-v2.py --input {input} --start {params.R1start} --end {params.R1end} --padding 0 --output {output.R1} ; "
		"scripts/filter-pyNAST-for-ROI-v2.py --input {input} --start {params.R2start} --end {params.R2end} --padding 0 --output {output.R2} ; "
		"scripts/filter-pyNAST-for-ROI-v2.py --input {input} --start {params.R1andR2start} --end {params.R1andR2end} --padding 0 --output {output.R1andR2} "

rule classify_unmatched:
	input:
		R1="compareTags2MG/10-MG-not-in-ASVdb-subsetted/{sample}.PROK.515Y-926R.not-matching.ASVs.aligned.R1.fasta",
		R2="compareTags2MG/10-MG-not-in-ASVdb-subsetted/{sample}.PROK.515Y-926R.not-matching.ASVs.aligned.R2.fasta",
		R1andR2="compareTags2MG/10-MG-not-in-ASVdb-subsetted/{sample}.PROK.515Y-926R.not-matching.ASVs.aligned.R1andR2.fasta"
	output:
		R1="compareTags2MG/11-MG-not-in-ASVdb-classified/{sample}.PROK.515Y-926R.not-matching.ASVs.aligned.R1.class.tsv",
		R2="compareTags2MG/11-MG-not-in-ASVdb-classified/{sample}.PROK.515Y-926R.not-matching.ASVs.aligned.R2.class.tsv",
		R1andR2="compareTags2MG/11-MG-not-in-ASVdb-classified/{sample}.PROK.515Y-926R.not-matching.ASVs.aligned.R1andR2.class.tsv"
	params:
		db="/home/db/VSEARCH/silva132_99_sintax.udb",
		options="--sintax_cutoff 0 --top_hits_only --topn 1 --notrunclabels",
		threads=1
	conda:
		"envs/vsearch.yaml"
	shell:
		"vsearch --sintax {input.R1} --db {params.db} --tabbedout {output.R1} --threads {params.threads} {params.options} ; "
		"vsearch --sintax {input.R2} --db {params.db} --tabbedout {output.R2} --threads {params.threads} {params.options} ; "
		"vsearch --sintax {input.R1andR2} --db {params.db} --tabbedout {output.R1andR2} --threads {params.threads} {params.options}"

#Concatenate taxonomy files
rule cat_tax_for_all_samples_matches_and_mismatches:
	output:
		R1="compareTags2MG/12-MG-not-in-ASVdb-classified-cat/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1.class.cat.tsv",
		R2="compareTags2MG/12-MG-not-in-ASVdb-classified-cat/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R2.class.cat.tsv",
		R1andR2="compareTags2MG/12-MG-not-in-ASVdb-classified-cat/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1andR2.class.cat.tsv"
	shell:
		"find compareTags2MG/11-MG-not-in-ASVdb-classified/ -type f -name "
		"\"*.R1.class.tsv\" -print0 | "
		"xargs -0 cat > {output.R1} ; "
		"find compareTags2MG/11-MG-not-in-ASVdb-classified/ -type f -name "
		"\"*.R2.class.tsv\" -print0 | "
		"xargs -0 cat > {output.R2} ; "
		"find compareTags2MG/11-MG-not-in-ASVdb-classified/ -type f -name "
		"\"*.R1andR2.class.tsv\" -print0 | "
		"xargs -0 cat > {output.R1andR2} ; "


#Counting order-level groupings (can adjust level with the "cut -d, -f1-4" parameter below)
rule count_tax_matches_and_mismatches:
	input:
		R1="compareTags2MG/12-MG-not-in-ASVdb-classified-cat/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1.class.cat.tsv",
		R2="compareTags2MG/12-MG-not-in-ASVdb-classified-cat/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R2.class.cat.tsv",
		R1andR2="compareTags2MG/12-MG-not-in-ASVdb-classified-cat/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1andR2.class.cat.tsv"
	output:
		R1counts="compareTags2MG/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1.class.cat.counts.tsv",
		R2counts="compareTags2MG/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R2.class.cat.counts.tsv",
		R1andR2counts="compareTags2MG/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1andR2.class.cat.counts.tsv",
		R1taxtable="compareTags2MG/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1.class.cat.taxtable.tsv",
		R2taxtable="compareTags2MG/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R2.class.cat.taxtable.tsv",
		R1andR2taxtable="compareTags2MG/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1andR2.class.cat.taxtable.tsv"
	shell:
		"sed -re 's/\([0-9]{{1}}\.[0-9]{{2}}\)//g' {input.R1} | tee {output.R1taxtable} |" #Remove confidence estimations from VSEARCH output, keep a copy for later steps but also pipe to subsequent commands
		"cut -f2 | sort | cut -d, -f1-4 | sort | uniq -c | " #Take only tax column, collapse to order level, then count unique occurrences
		"tail -f -n +2 | awk '{{print $1,\"\t\",$2}}' > {output.R1counts} ; " #Process output into tsv format to stdout
		"sed -re 's/\([0-9]{{1}}\.[0-9]{{2}}\)//g' {input.R2} | tee {output.R2taxtable} |" #Remove confidence estimations from VSEARCH output, keep a copy for later steps but also pipe to subsequent commands
		"cut -f2 | sort | cut -d, -f1-4 | sort | uniq -c | " #Take only tax column, collapse to order level, then count unique occurrences
		"tail -f -n +2 | awk '{{print $1,\"\t\",$2}}' > {output.R2counts} ; " #Process output into tsv format to stdout
		"sed -re 's/\([0-9]{{1}}\.[0-9]{{2}}\)//g' {input.R1andR2} | tee {output.R1andR2taxtable} |" #Remove confidence estimations from VSEARCH output, keep a copy for later steps but also pipe to subsequent commands
		"cut -f2 | sort | cut -d, -f1-4 | sort | uniq -c | " #Take only tax column, collapse to order level, then count unique occurrences
		"tail -f -n +2 | awk '{{print $1,\"\t\",$2}}' > {output.R1andR2counts} ; " #Process output into tsv format to stdout

#Now take only those with greater than 1 % abundance (among mismatches) using basic python script (can change abundance cutoff if you desire)
rule filter_tax_matches_by_abundance:
	input:
		R1="compareTags2MG/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1.class.cat.counts.tsv",
		R2="compareTags2MG/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R2.class.cat.counts.tsv",
		R1andR2="compareTags2MG/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1andR2.class.cat.counts.tsv"
	params:
		minAbund = 0.001  #Change fractional value in config file if desired, default 0.01
	output:
		R1="compareTags2MG/14-MG-not-in-ASVdb-classified-cat-parsed-fractions/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1.class.cat.counts.frac.min" + str(CUTOFF) + ".tsv",
		R2="compareTags2MG/14-MG-not-in-ASVdb-classified-cat-parsed-fractions/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R2.class.cat.counts.frac.min" + str(CUTOFF) + ".tsv",
		R1andR2="compareTags2MG/14-MG-not-in-ASVdb-classified-cat-parsed-fractions/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1andR2.class.cat.counts.frac.min" + str(CUTOFF) + ".tsv"
	shell:
		"scripts/mismatch-characterization/filter-by-fractional-abundance.py {input.R1} {params.minAbund} > {output.R1} ; "
		"scripts/mismatch-characterization/filter-by-fractional-abundance.py {input.R2} {params.minAbund} > {output.R2} ; "
		"scripts/mismatch-characterization/filter-by-fractional-abundance.py {input.R1andR2} {params.minAbund} > {output.R1andR2} "

rule make_kronas:
	input:
		R1="compareTags2MG/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1.class.cat.counts.tsv",
		R2="compareTags2MG/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R2.class.cat.counts.tsv",
		R1andR2="compareTags2MG/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1andR2.class.cat.counts.tsv"
	output:
		R1="compareTags2MG/15-MG-not-in-ASVdb-classified-kronas/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1.krona.input",
		R2="compareTags2MG/15-MG-not-in-ASVdb-classified-kronas/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R2.krona.input",
		R1andR2="compareTags2MG/15-MG-not-in-ASVdb-classified-kronas/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1andR2.krona.input",
		R1html="compareTags2MG/15-MG-not-in-ASVdb-classified-kronas/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1.krona.html",
		R2html="compareTags2MG/15-MG-not-in-ASVdb-classified-kronas/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R2.krona.html",
		R1andR2html="compareTags2MG/15-MG-not-in-ASVdb-classified-kronas/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1andR2.krona.html"
	conda:
		"envs/krona.yaml"
	shell:
		"sed 's/,/\t/g' {input.R1} | sed 's/\w://g' > {output.R1} ; ktImportText -c -n \"R1 reads\" -o {output.R1html} {output.R1} ; "
		"sed 's/,/\t/g' {input.R2} | sed 's/\w://g' > {output.R2} ; ktImportText -c -n \"R2 reads\" -o {output.R2html} {output.R2} ; "
		"sed 's/,/\t/g' {input.R1andR2} | sed 's/\w://g' > {output.R1andR2} ; ktImportText -c -n \"R1 and R2 reads\" -o {output.R1andR2html} {output.R1andR2} ; "
"""
rule get_pc_matching:

rule sift_unmatched:

#get proportion
rule count_positive_BLAST_hits:

rule plot_proportions:


#The metagenomic SSU rRNA that does NOT match the ASVs
rule classify_negative_BLAST_hits:


Need to compare sample to sample. To do:
1. Make lookup table between SRAid and sampleID
2. Get ASV subset for EACH sample. Take only non-zero abundance ASVs based on ASV table.
"""
