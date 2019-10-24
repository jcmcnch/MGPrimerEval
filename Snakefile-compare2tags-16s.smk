from datetime import datetime
datestamp=datetime.today().strftime('%Y-%m-%d')
cutoff=config["cutoff"]
pcid=config["pcid"]
strCutoff="minAbund-" + str(config["cutoff"])
strPcid="blastnPcID-" + str(config["pcid"])
denoiser=config["denoiser"]
outdir="compareTags2MG/" + '_'.join([datestamp, denoiser, strPcid, strCutoff, "vs_MG"])
iLenDeblurTrunc=config["iLenDeblurTrunc"] #e.g. for the 515Y/926R amplicons, the merged reads are 373bp but I typically truncate with deblur to 363bp, so this value would be equal to 10
ASVtable=config["ASVtable"]
ASVseqs=config["ASVseqs"]

rule all:
	input:
		expand("{outdir}/07-MG-vs-ASV-plots/{sample}.PROK.nonzero.ASV.comparison.svg", sample=config["samples"], outdir="compareTags2MG/" + '_'.join([datestamp, denoiser, strPcid, strCutoff, "vs_MG"])),
		expand("{outdir}/09-MG-not-in-ASVdb-aligned/{sample}.PROK.515Y-926R.not-matching.ASVs.aligned.fasta", sample=config["samples"], outdir="compareTags2MG/" + '_'.join([datestamp, denoiser, strPcid, strCutoff, "vs_MG"])),
		expand("{outdir}/15-MG-not-in-ASVdb-classified-kronas/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1andR2.krona.html", outdir="compareTags2MG/" + '_'.join([datestamp, denoiser, strPcid, strCutoff, "vs_MG"]))

#Generate concatenated alignments that will be used for subsetting
rule concat_fwd_and_reverse_alignments:
	input:
		fwd="05-pyNAST-aligned/{sample}.fwd.{group}/{sample}.fwd.SSU.{group}_pynast_aligned.fasta",
		rev="05-pyNAST-aligned/{sample}.rev.{group}/{sample}.rev.SSU.{group}_pynast_aligned.fasta"
	output:
		temp("{outdir}/00-concatenated/{sample}.{group}.concat.fasta")
	shell:
		"cat {input.fwd} {input.rev} > {output}"

rule subset_to_primer_region:
	input:
		"{outdir}/00-concatenated/{sample}.{group}.concat.fasta"
	output:
		"{outdir}/01-concatenated/{sample}.{group}.concat.515Y-926R.fasta"
	params:
		start=lambda wildcards: config["primerROI"][wildcards.group]["515Y"][1],
		end=lambda wildcards: config["primerROI"][wildcards.group]["926R"][0] - iLenDeblurTrunc #Necessary if amplicon region truncated during default deblur pipeline
	shell:
		"scripts/filter-pyNAST-for-ROI-v2.py --input {input} --start {params.start} --end {params.end} --padding 0 --output {output}"

#concatenate the 3 different PROK categories
rule merge_PROK:
	input:
		"{outdir}/01-concatenated/{sample}.ARCH.concat.515Y-926R.fasta",
		"{outdir}/01-concatenated/{sample}.BACT-NON-CYANO.concat.515Y-926R.fasta",
		"{outdir}/01-concatenated/{sample}.BACT-CYANO.concat.515Y-926R.fasta"
	output:
		"{outdir}/01-concatenated/{sample}.PROK.concat.515Y-926R.fasta"
	shell:
		"cat {input} > {output}"

#get non-zero abundance ASVs from table, and output the ids
rule parse_16S_ASV_table:
	input:
		"GP13-sample-SRA.tsv",
		expand({ASVtable}, ASVtable=config["ASVtable"])
	params:
		"{sample}"
	output:
		"{outdir}/02-ASV-ids/{sample}.PROK.nonzero.ASV.ids"
	script:
		"scripts/make-db-of-non-zero-abund-ASV.py"

rule get_sample_nonzero_16S_ASV_fastas:
	input:
		fasta=expand("{ASVseqs}", ASVseqs=config["ASVseqs"]),
		ids="{outdir}/02-ASV-ids/{sample}.PROK.nonzero.ASV.ids"
	output:
		"{outdir}/03-ASV-fastas/{sample}.PROK.nonzero.ASV.fasta"
	shell:
		"seqtk subseq {input.fasta} {input.ids} > {output}"

rule make_blast_dbs_16S:
	input:
		"{outdir}/03-ASV-fastas/{sample}.PROK.nonzero.ASV.fasta"
	output:
		expand("{{outdir}}/04-ASV-blastdbs/{{sample}}.PROK.nonzero.ASV.db.{ext}", ext=["nhr", "nin", "nsq"])
	params:
		filestem="{outdir}/04-ASV-blastdbs/{sample}.PROK.nonzero.ASV.db"
	conda:
		"envs/qiime1.yaml"
	shell:
		"makeblastdb -in {input} -dbtype nucl -out {params.filestem} ; touch {output}"

rule blast_MG_vs_tags:
	input:
		database_files=lambda wildcards: expand("{{outdir}}/04-ASV-blastdbs/{sample}.PROK.nonzero.ASV.db.{ext}", ext=["nhr", "nin", "nsq"], sample=wildcards.sample),
		query="{outdir}/01-concatenated/{sample}.PROK.concat.515Y-926R.fasta"
	output:
		"{outdir}/05-MG-blasted-against-ASVs/{sample}.PROK.nonzero.ASV.blastout.tsv"
	params:
		dbname="{outdir}/04-ASV-blastdbs/{sample}.PROK.nonzero.ASV.db"
	conda:
		"envs/qiime1.yaml"
	shell:
		"blastn -qcov_hsp_perc 100 -perc_identity {pcid} -outfmt 6 -query {input.query} -db {params.dbname} > {output}"

rule compare_MG_SSU_rRNA_with_ASVs:
	input:
		"GP13-sample-SRA.tsv",
		expand({ASVtable}, ASVtable=config["ASVtable"]),
		"{outdir}/05-MG-blasted-against-ASVs/{sample}.PROK.nonzero.ASV.blastout.tsv"
	params:
		"{sample}",
		"{outdir}/"
	conda:
		"envs/networkx.yaml"
	output:
		"{outdir}/06-MG-vs-ASV-tsv/{sample}.PROK.nonzero.ASV.comparison.tsv"
	script:
		"scripts/compare-MG-eASV-abund-from-blast-output.py"

rule plot_ASV_vs_BLAST_results:
	input:
		"{outdir}/06-MG-vs-ASV-tsv/{sample}.PROK.nonzero.ASV.comparison.tsv"
	params:
		"{sample}"
	output:
		"{outdir}/07-MG-vs-ASV-plots/{sample}.PROK.nonzero.ASV.comparison.svg"
	script:
		"scripts/seaborn-plot-correlations.py"

rule sift_unmatched:
	input:
		"{outdir}/05-MG-blasted-against-ASVs/{sample}.PROK.nonzero.ASV.blastout.tsv",
		"{outdir}/01-concatenated/{sample}.PROK.concat.515Y-926R.fasta"
	output:
		"{outdir}/08-MG-not-in-ASVdb/{sample}.PROK.515Y-926R.not-matching.ASVs.fasta"
	script:
		"scripts/get-non-matching.py"

#Align against E. coli for now, since just to get visual of locations
rule align_unmatched:
	input:
		seqs="{outdir}/08-MG-not-in-ASVdb/{sample}.PROK.515Y-926R.not-matching.ASVs.fasta",
		ref="SSU_refs/Ecoli_16s.fna"
	output:
		"{outdir}/09-MG-not-in-ASVdb-aligned/{sample}.PROK.515Y-926R.not-matching.ASVs.aligned.fasta"
	conda:
		"envs/pynast.yaml"
	shell:
		"pynast -p 10 -l 1 -i {input.seqs} -t {input.ref} -a {output}"

rule subset_non_matching:
	input:
		"{outdir}/09-MG-not-in-ASVdb-aligned/{sample}.PROK.515Y-926R.not-matching.ASVs.aligned.fasta"
	output:
		R1="{outdir}/10-MG-not-in-ASVdb-subsetted/{sample}.PROK.515Y-926R.not-matching.ASVs.aligned.R1.fasta",
		R2="{outdir}/10-MG-not-in-ASVdb-subsetted/{sample}.PROK.515Y-926R.not-matching.ASVs.aligned.R2.fasta",
		R1andR2="{outdir}/10-MG-not-in-ASVdb-subsetted/{sample}.PROK.515Y-926R.not-matching.ASVs.aligned.R1andR2.fasta"
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
		R1="{outdir}/10-MG-not-in-ASVdb-subsetted/{sample}.PROK.515Y-926R.not-matching.ASVs.aligned.R1.fasta",
		R2="{outdir}/10-MG-not-in-ASVdb-subsetted/{sample}.PROK.515Y-926R.not-matching.ASVs.aligned.R2.fasta",
		R1andR2="{outdir}/10-MG-not-in-ASVdb-subsetted/{sample}.PROK.515Y-926R.not-matching.ASVs.aligned.R1andR2.fasta"
	output:
		R1="{outdir}/11-MG-not-in-ASVdb-classified/{sample}.PROK.515Y-926R.not-matching.ASVs.aligned.R1.class.tsv",
		R2="{outdir}/11-MG-not-in-ASVdb-classified/{sample}.PROK.515Y-926R.not-matching.ASVs.aligned.R2.class.tsv",
		R1andR2="{outdir}/11-MG-not-in-ASVdb-classified/{sample}.PROK.515Y-926R.not-matching.ASVs.aligned.R1andR2.class.tsv"
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
	input:
		R1=expand("{{outdir}}/11-MG-not-in-ASVdb-classified/{sample}.PROK.515Y-926R.not-matching.ASVs.aligned.R1.class.tsv", sample=config["samples"]),
		R2=expand("{{outdir}}/11-MG-not-in-ASVdb-classified/{sample}.PROK.515Y-926R.not-matching.ASVs.aligned.R2.class.tsv", sample=config["samples"]),
		R1andR2=expand("{{outdir}}/11-MG-not-in-ASVdb-classified/{sample}.PROK.515Y-926R.not-matching.ASVs.aligned.R1andR2.class.tsv", sample=config["samples"])
	output:
		R1="{outdir}/12-MG-not-in-ASVdb-classified-cat/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1.class.cat.tsv",
		R2="{outdir}/12-MG-not-in-ASVdb-classified-cat/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R2.class.cat.tsv",
		R1andR2="{outdir}/12-MG-not-in-ASVdb-classified-cat/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1andR2.class.cat.tsv"
	shell:
		"find {outdir}/11-MG-not-in-ASVdb-classified/ -type f -name "
		"\"*.R1.class.tsv\" -print0 | "
		"xargs -0 cat > {output.R1} ; "
		"find {outdir}/11-MG-not-in-ASVdb-classified/ -type f -name "
		"\"*.R2.class.tsv\" -print0 | "
		"xargs -0 cat > {output.R2} ; "
		"find {outdir}/11-MG-not-in-ASVdb-classified/ -type f -name "
		"\"*.R1andR2.class.tsv\" -print0 | "
		"xargs -0 cat > {output.R1andR2} ; "


#Counting order-level groupings (can adjust level with the "cut -d, -f1-4" parameter below)
rule count_tax_matches_and_mismatches:
	input:
		R1="{outdir}/12-MG-not-in-ASVdb-classified-cat/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1.class.cat.tsv",
		R2="{outdir}/12-MG-not-in-ASVdb-classified-cat/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R2.class.cat.tsv",
		R1andR2="{outdir}/12-MG-not-in-ASVdb-classified-cat/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1andR2.class.cat.tsv"
	output:
		R1counts="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1.class.cat.counts.tsv",
		R2counts="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R2.class.cat.counts.tsv",
		R1andR2counts="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1andR2.class.cat.counts.tsv",
		R1taxtable="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1.class.cat.taxtable.tsv",
		R2taxtable="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R2.class.cat.taxtable.tsv",
		R1andR2taxtable="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1andR2.class.cat.taxtable.tsv"
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
		R1="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1.class.cat.counts.tsv",
		R2="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R2.class.cat.counts.tsv",
		R1andR2="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1andR2.class.cat.counts.tsv"
	params:
		minAbund = 0.001  #Change fractional value in config file if desired, default 0.01
	output:
		R1="{outdir}/14-MG-not-in-ASVdb-classified-cat-parsed-fractions/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1.class.cat.counts.frac.min" + str(cutoff) + ".tsv",
		R2="{outdir}/14-MG-not-in-ASVdb-classified-cat-parsed-fractions/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R2.class.cat.counts.frac.min" + str(cutoff) + ".tsv",
		R1andR2="{outdir}/14-MG-not-in-ASVdb-classified-cat-parsed-fractions/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1andR2.class.cat.counts.frac.min" + str(cutoff) + ".tsv"
	shell:
		"scripts/mismatch-characterization/filter-by-fractional-abundance.py {input.R1} {params.minAbund} > {output.R1} ; "
		"scripts/mismatch-characterization/filter-by-fractional-abundance.py {input.R2} {params.minAbund} > {output.R2} ; "
		"scripts/mismatch-characterization/filter-by-fractional-abundance.py {input.R1andR2} {params.minAbund} > {output.R1andR2} "

rule make_kronas:
	input:
		R1="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1.class.cat.counts.tsv",
		R2="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R2.class.cat.counts.tsv",
		R1andR2="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1andR2.class.cat.counts.tsv"
	output:
		R1="{outdir}/15-MG-not-in-ASVdb-classified-kronas/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1.krona.input",
		R2="{outdir}/15-MG-not-in-ASVdb-classified-kronas/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R2.krona.input",
		R1andR2="{outdir}/15-MG-not-in-ASVdb-classified-kronas/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1andR2.krona.input",
		R1html="{outdir}/15-MG-not-in-ASVdb-classified-kronas/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1.krona.html",
		R2html="{outdir}/15-MG-not-in-ASVdb-classified-kronas/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R2.krona.html",
		R1andR2html="{outdir}/15-MG-not-in-ASVdb-classified-kronas/GP13-PROK.515Y-926R.not-matching.ASVs.aligned.R1andR2.krona.html"
	conda:
		"envs/krona.yaml"
	shell:
		"sed 's/,/\t/g' {input.R1} | sed 's/\w://g' > {output.R1} ; ktImportText -c -n \"R1 reads\" -o {output.R1html} {output.R1} ; "
		"sed 's/,/\t/g' {input.R2} | sed 's/\w://g' > {output.R2} ; ktImportText -c -n \"R2 reads\" -o {output.R2html} {output.R2} ; "
		"sed 's/,/\t/g' {input.R1andR2} | sed 's/\w://g' > {output.R1andR2} ; ktImportText -c -n \"R1 and R2 reads\" -o {output.R1andR2html} {output.R1andR2} "
