from datetime import datetime
datestamp=datetime.today().strftime('%Y-%m-%d')
cutoff=config["cutoff"]
pcid=config["pcid"]
strCutoff="minAbund-" + str(config["cutoff"])
strPcid="blastnPcID-" + str(config["pcid"])
denoiser=config["denoiser"]
outdir="compare-workflow-intermediate/" + '_'.join([datestamp, denoiser, strPcid, strCutoff, "vs_MG-18S"])
ASVtable=config["ASVtable"]
ASVseqs=config["ASVseqs"]
iLenR1Trunc=config["iLenR1Trunc"]
iLenR2Trunc=config["iLenR2Trunc"]

rule all:
	input:
		expand("{outdir}/01-subsetted/{sample}.EUK.cleaned.515Y-926R.revcomped.sliced.fasta", sample=config["samples"], outdir="compare-workflow-intermediate/" + '_'.join([datestamp, denoiser, strPcid, strCutoff, "vs_MG"])),
		expand("{outdir}/07-MG-vs-ASV-plots/{sample}.EUK.nonzero.ASV.comparison.svg", sample=config["samples"], outdir="compare-workflow-intermediate/" + '_'.join([datestamp, denoiser, strPcid, strCutoff, "vs_MG-18S"])),
		expand("{outdir}/09-MG-not-in-ASVdb-aligned/{sample}.EUK.515Y-926R.not-matching.ASVs.aligned.fasta", sample=config["samples"], outdir="compare-workflow-intermediate/" + '_'.join([datestamp, denoiser, strPcid, strCutoff, "vs_MG-18S"])),
		expand("{outdir}/15-MG-not-in-ASVdb-classified-kronas/GP13-EUK.515Y-926R.not-matching.ASVs.aligned.R1andR2.krona.html", outdir="compare-workflow-intermediate/" + '_'.join([datestamp, denoiser, strPcid, strCutoff, "vs_MG-18S"]))

#Generate concatenated alignments that will be used for subsetting
rule concat_fwd_and_reverse_alignments:
	input:
		fwd="compute-workflow-intermediate/05-pyNAST-aligned/{sample}.fwd.SSU.EUK_pynast_aligned.fasta",
		rev="compute-workflow-intermediate/05-pyNAST-aligned/{sample}.rev.SSU.EUK_pynast_aligned.fasta"
	output:
		"{outdir}/00-concatenated/{sample}.EUK.concat.fasta"
	shell:
		"cat {input.fwd} {input.rev} > {output}"

rule subset_to_primer_region:
	input:
		"{outdir}/00-concatenated/{sample}.EUK.concat.fasta"
	output:
		R1="{outdir}/01-subsetted/{sample}.EUK.concat.515Y-926R.R1.fasta",
		R2="{outdir}/01-subsetted/{sample}.EUK.concat.515Y-926R.R2.fasta",
		R1andR2="{outdir}/01-subsetted/{sample}.EUK.concat.515Y-926R.fasta"
	params:
		startR1=lambda wildcards: config["primerROI"]["EUK"]["515Y"][1],
		endR1=lambda wildcards: config["primerROI"]["EUK"]["515Y"][1] + iLenR1Trunc,
		startR2=lambda wildcards: config["primerROI"]["EUK"]["926R"][0] - iLenR2Trunc,
		endR2=lambda wildcards: config["primerROI"]["EUK"]["926R"][0]
	shell:
		"""
		scripts/filter-pyNAST-for-ROI-v2.py --input {input} --start {params.startR1} --end {params.endR1} --padding 0 --output {output.R1} ;
		scripts/filter-pyNAST-for-ROI-v2.py --input {input} --start {params.startR2} --end {params.endR2} --padding 0 --output {output.R2} ;
		cat {output.R1} {output.R2} > {output.R1andR2}
		"""

rule get_names:
	input:
		names="{outdir}/01-subsetted/{sample}.EUK.concat.515Y-926R.fasta"
	output:
		"{outdir}/01-subsetted/{sample}.EUK.concat.515Y-926R.ids"
	shell:
		"grep \">\" {input} | sed 's/>//' | awk '{{print $1\" \"$2\" \"$3}}' | sort | uniq > {output}"

rule get_fastq_by_group:
	input:
		names="{outdir}/01-subsetted/{sample}.EUK.concat.515Y-926R.ids",
		fwd="compute-workflow-intermediate/03-low-complexity-filtered/{sample}.fwd.SSU.keep.fastq",
		rev="compute-workflow-intermediate/03-low-complexity-filtered/{sample}.rev.SSU.keep.fastq"
	output:
		"{outdir}/01-subsetted/{sample}.EUK.concat.515Y-926R.fastq"
	conda:
		"envs/bbmap.yaml"
	shell:
		"""
		filterbyname.sh include=t substring=f app=t \
		names={input.names} in={input.fwd} out={output}
		filterbyname.sh include=t substring=f app=t \
		names={input.names} in={input.rev} out={output}
		"""

rule revcomp_fastq_according_to_pyNAST:
	input:
		fasta="{outdir}/00-concatenated/{sample}.EUK.concat.fasta",
		fastq="{outdir}/01-subsetted/{sample}.EUK.concat.515Y-926R.fastq"
	output:
		fastq_revcomp="{outdir}/01-subsetted/{sample}.EUK.concat.515Y-926R.revcomped.fastq"
	conda:
		"envs/biopython.yaml"
	shell:
		"./scripts/revcompfastq_according_to_pyNAST.py --inpynast {input.fasta} --infastq {input.fastq} --outfastq {output.fastq_revcomp}"

#Since the above fastq will have additional bases included that are not part of the amplicon region
rule get_sliced_fastq:
	input:
		fasta="{outdir}/01-subsetted/{sample}.EUK.concat.515Y-926R.fasta",
		fastq="{outdir}/01-subsetted/{sample}.EUK.concat.515Y-926R.revcomped.fastq"
	output:
		"{outdir}/01-subsetted/{sample}.EUK.concat.515Y-926R.revcomped.sliced.fastq"
	shell:
		"scripts/get-fastq-slice.py --inputfasta {input.fasta} "
		"--inputfastq {input.fastq} --output {output}"

#Additional QC to remove a few remaining homopolymer runs and other things komplexity did not catch
rule remove_low_complexity_bbmap:
	input:
		"{outdir}/01-subsetted/{sample}.EUK.concat.515Y-926R.revcomped.sliced.fastq"
	output:
		masked="{outdir}/01-subsetted/{sample}.EUK.masked.515Y-926R.revcomped.sliced.fastq",
		cleaned="{outdir}/01-subsetted/{sample}.EUK.cleaned.515Y-926R.revcomped.sliced.fasta"
	conda:
		"envs/bbmap.yaml"
	shell:
		"""
		bbmask.sh overwrite=t -Xmx4g in={input} out={output.masked} entropy=0.7 \
		minkr=4 maxkr=8 mr=t minlen=20 minke=4 maxke=8 fastawrap=0
		reformat.sh maxns=0 qtrim=t trimq=30 minlength=90 in={output.masked} out={output.cleaned}
		"""

#get non-zero abundance ASVs from table, and output the ids
rule parse_16S_ASV_table:
	input:
		"GP13-sample-SRA.tsv",
		expand("{ASVtable}", ASVtable=config["ASVtable"])
	params:
		"{sample}"
	output:
		"{outdir}/02-ASV-ids/{sample}.EUK.nonzero.ASV.ids"
	script:
		"scripts/make-db-of-non-zero-abund-ASV.py"

rule get_sample_nonzero_16S_ASV_fastas:
	input:
		fasta=expand("{ASVseqs}", ASVseqs=config["ASVseqs"]),
		ids="{outdir}/02-ASV-ids/{sample}.EUK.nonzero.ASV.ids"
	output:
		"{outdir}/03-ASV-fastas/{sample}.EUK.nonzero.ASV.fasta"
	shell:
		"seqtk subseq {input.fasta} {input.ids} > {output}"

rule make_blast_dbs_18S:
	input:
		"{outdir}/03-ASV-fastas/{sample}.EUK.nonzero.ASV.fasta"
	output:
		expand("{{outdir}}/04-ASV-blastdbs/{{sample}}.EUK.nonzero.ASV.db.{ext}", ext=["nhr", "nin", "nsq"])
	params:
		filestem="{outdir}/04-ASV-blastdbs/{sample}.EUK.nonzero.ASV.db"
	conda:
		"envs/qiime1.yaml"
	shell:
		"makeblastdb -in {input} -dbtype nucl -out {params.filestem} ; touch {output}"

rule blast_MG_vs_tags:
	input:
		database_files=lambda wildcards: expand("{{outdir}}/04-ASV-blastdbs/{sample}.EUK.nonzero.ASV.db.{ext}", ext=["nhr", "nin", "nsq"], sample=wildcards.sample),
		query="{outdir}/01-subsetted/{sample}.EUK.cleaned.515Y-926R.revcomped.sliced.fasta"
	output:
		"{outdir}/05-MG-blasted-against-ASVs/{sample}.EUK.nonzero.ASV.blastout.tsv"
	params:
		dbname="{outdir}/04-ASV-blastdbs/{sample}.EUK.nonzero.ASV.db"
	conda:
		"envs/qiime1.yaml"
	shell:
		"blastn -qcov_hsp_perc 100 -perc_identity {pcid} -outfmt 6 -query {input.query} -db {params.dbname} > {output}"

rule compare_MG_SSU_rRNA_with_ASVs:
	input:
		"GP13-sample-SRA.tsv",
		expand("{ASVtable}", ASVtable=config["ASVtable"]),
		"{outdir}/05-MG-blasted-against-ASVs/{sample}.EUK.nonzero.ASV.blastout.tsv"
	params:
		"{sample}",
		"{outdir}/"
	conda:
		"envs/networkx.yaml"
	output:
		"{outdir}/06-MG-vs-ASV-tsv/{sample}.EUK.nonzero.ASV.comparison.tsv"
	script:
		"scripts/compare-MG-eASV-abund-from-blast-output-18S.py"

rule plot_ASV_vs_BLAST_results:
	input:
		"{outdir}/06-MG-vs-ASV-tsv/{sample}.EUK.nonzero.ASV.comparison.tsv"
	params:
		"{sample}"
	output:
		"{outdir}/07-MG-vs-ASV-plots/{sample}.EUK.nonzero.ASV.comparison.svg"
	script:
		"scripts/seaborn-plot-correlations.py"

rule sift_unmatched:
	input:
		"{outdir}/05-MG-blasted-against-ASVs/{sample}.EUK.nonzero.ASV.blastout.tsv",
		"{outdir}/01-subsetted/{sample}.EUK.concat.515Y-926R.fasta"
	output:
		"{outdir}/08-MG-not-in-ASVdb/{sample}.EUK.515Y-926R.not-matching.ASVs.fasta",
		"{outdir}/08-MG-not-in-ASVdb/{sample}.EUK.515Y-926R.matching.ASVs.fasta"
	script:
		"scripts/get-non-matching.py"

#Align against Saccharomyces_cerevisiae_S288C_18s to get visual of locations
rule align_unmatched:
	input:
		missed="{outdir}/08-MG-not-in-ASVdb/{sample}.EUK.515Y-926R.not-matching.ASVs.fasta",
		matched="{outdir}/08-MG-not-in-ASVdb/{sample}.EUK.515Y-926R.matching.ASVs.fasta",
		ref="SSU_refs/Saccharomyces_cerevisiae_S288C_18s-1_NR_132213.1.fa"
	output:
		matching="{outdir}/09-MG-not-in-ASVdb-aligned/{sample}.EUK.515Y-926R.matching.ASVs.aligned.fasta",
		missing="{outdir}/09-MG-not-in-ASVdb-aligned/{sample}.EUK.515Y-926R.not-matching.ASVs.aligned.fasta"
	conda:
		"envs/pynast.yaml"
	shell:
		"pynast -p 10 -l 1 -i {input.missed} -t {input.ref} -a {output.missing} ; "
		"pynast -p 10 -l 1 -i {input.matched} -t {input.ref} -a {output.matching}"

rule subset_non_matching:
	input:
		matching="{outdir}/09-MG-not-in-ASVdb-aligned/{sample}.EUK.515Y-926R.matching.ASVs.aligned.fasta",
		missing="{outdir}/09-MG-not-in-ASVdb-aligned/{sample}.EUK.515Y-926R.not-matching.ASVs.aligned.fasta"
	output:
		#R1missing="{outdir}/10-MG-not-in-ASVdb-subsetted/{sample}.EUK.515Y-926R.not-matching.ASVs.aligned.R1.fasta",
		#R2missing="{outdir}/10-MG-not-in-ASVdb-subsetted/{sample}.EUK.515Y-926R.not-matching.ASVs.aligned.R2.fasta",
		R1andR2missing="{outdir}/10-MG-not-in-ASVdb-subsetted/{sample}.EUK.515Y-926R.not-matching.ASVs.aligned.R1andR2.fasta",
		#R1matching="{outdir}/10-MG-not-in-ASVdb-subsetted/{sample}.EUK.515Y-926R.matching.ASVs.aligned.R1.fasta",
		#R2matching="{outdir}/10-MG-not-in-ASVdb-subsetted/{sample}.EUK.515Y-926R.matching.ASVs.aligned.R2.fasta",
		R1andR2matching="{outdir}/10-MG-not-in-ASVdb-subsetted/{sample}.EUK.515Y-926R.matching.ASVs.aligned.R1andR2.fasta"
	params:
		R1start=581,
		R1end=581 + iLenR1Trunc,
		R2start=1130 - iLenR2Trunc,
		R2end=1130,
		R1andR2start=581,
		R1andR2end=1130
	shell:
		#"scripts/filter-pyNAST-for-ROI-v2.py --input {input.matching} --start {params.R1start} --end {params.R1end} --padding 0 --output {output.R1matching} ; "
		#"scripts/filter-pyNAST-for-ROI-v2.py --input {input.matching} --start {params.R2start} --end {params.R2end} --padding 0 --output {output.R2matching} ; "
		"scripts/filter-pyNAST-for-ROI-v2.py --input {input.matching} --start {params.R1andR2start} --end {params.R1andR2end} --padding 0 --output {output.R1andR2matching} ; "
		#"scripts/filter-pyNAST-for-ROI-v2.py --input {input.missing} --start {params.R1start} --end {params.R1end} --padding 0 --output {output.R1missing} ; "
		#"scripts/filter-pyNAST-for-ROI-v2.py --input {input.missing} --start {params.R2start} --end {params.R2end} --padding 0 --output {output.R2missing} ; "
		"scripts/filter-pyNAST-for-ROI-v2.py --input {input.missing} --start {params.R1andR2start} --end {params.R1andR2end} --padding 0 --output {output.R1andR2missing}"

rule classify_unmatched:
	input:
		#R1missing="{outdir}/10-MG-not-in-ASVdb-subsetted/{sample}.EUK.515Y-926R.not-matching.ASVs.aligned.R1.fasta",
		#R2missing="{outdir}/10-MG-not-in-ASVdb-subsetted/{sample}.EUK.515Y-926R.not-matching.ASVs.aligned.R2.fasta",
		R1andR2missing="{outdir}/10-MG-not-in-ASVdb-subsetted/{sample}.EUK.515Y-926R.not-matching.ASVs.aligned.R1andR2.fasta",
		#R1matching="{outdir}/10-MG-not-in-ASVdb-subsetted/{sample}.EUK.515Y-926R.matching.ASVs.aligned.R1.fasta",
		#R2matching="{outdir}/10-MG-not-in-ASVdb-subsetted/{sample}.EUK.515Y-926R.matching.ASVs.aligned.R2.fasta",
		R1andR2matching="{outdir}/10-MG-not-in-ASVdb-subsetted/{sample}.EUK.515Y-926R.matching.ASVs.aligned.R1andR2.fasta"
	output:
		#R1missing="{outdir}/11-MG-not-in-ASVdb-classified/{sample}.EUK.515Y-926R.not-matching.ASVs.aligned.R1.class.tsv",
		#R2missing="{outdir}/11-MG-not-in-ASVdb-classified/{sample}.EUK.515Y-926R.not-matching.ASVs.aligned.R2.class.tsv",
		R1andR2missing="{outdir}/11-MG-not-in-ASVdb-classified/{sample}.EUK.515Y-926R.not-matching.ASVs.aligned.R1andR2.class.tsv",
		#R1matching="{outdir}/11-MG-not-in-ASVdb-classified/{sample}.EUK.515Y-926R.matching.ASVs.aligned.R1.class.tsv",
		#R2matching="{outdir}/11-MG-not-in-ASVdb-classified/{sample}.EUK.515Y-926R.matching.ASVs.aligned.R2.class.tsv",
		R1andR2matching="{outdir}/11-MG-not-in-ASVdb-classified/{sample}.EUK.515Y-926R.matching.ASVs.aligned.R1andR2.class.tsv"
	params:
		db="/home/db/VSEARCH/silva132_99_sintax.udb",
		options="--sintax_cutoff 0 --top_hits_only --topn 1 --notrunclabels",
		threads=1
	conda:
		"envs/vsearch.yaml"
	shell:
		#"vsearch --sintax {input.R1missing} --db {params.db} --tabbedout {output.R1missing} --threads {params.threads} {params.options} ; "
		#"vsearch --sintax {input.R2missing} --db {params.db} --tabbedout {output.R2missing} --threads {params.threads} {params.options} ; "
		"vsearch --sintax {input.R1andR2missing} --db {params.db} --tabbedout {output.R1andR2missing} --threads {params.threads} {params.options} ; "
		#"vsearch --sintax {input.R1matching} --db {params.db} --tabbedout {output.R1missing} --threads {params.threads} {params.options} ; "
		#"vsearch --sintax {input.R2matching} --db {params.db} --tabbedout {output.R2missing} --threads {params.threads} {params.options} ; "
		"vsearch --sintax {input.R1andR2matching} --db {params.db} --tabbedout {output.R1andR2matching} --threads {params.threads} {params.options}"

#Concatenate taxonomy files
rule cat_tax_for_all_samples_matches_and_mismatches:
	input:
		#R1missing=expand("{{outdir}}/11-MG-not-in-ASVdb-classified/{sample}.EUK.515Y-926R.not-matching.ASVs.aligned.R1.class.tsv", sample=config["samples"]),
		#R2missing=expand("{{outdir}}/11-MG-not-in-ASVdb-classified/{sample}.EUK.515Y-926R.not-matching.ASVs.aligned.R2.class.tsv", sample=config["samples"]),
		R1andR2missing=expand("{{outdir}}/11-MG-not-in-ASVdb-classified/{sample}.EUK.515Y-926R.not-matching.ASVs.aligned.R1andR2.class.tsv", sample=config["samples"]),
		#R1matching=expand("{{outdir}}/11-MG-not-in-ASVdb-classified/{sample}.EUK.515Y-926R.matching.ASVs.aligned.R1.class.tsv", sample=config["samples"]),
		#R2matching=expand("{{outdir}}/11-MG-not-in-ASVdb-classified/{sample}.EUK.515Y-926R.matching.ASVs.aligned.R2.class.tsv", sample=config["samples"]),
		R1andR2matching=expand("{{outdir}}/11-MG-not-in-ASVdb-classified/{sample}.EUK.515Y-926R.matching.ASVs.aligned.R1andR2.class.tsv", sample=config["samples"])
	output:
		#R1missing="{outdir}/12-MG-not-in-ASVdb-classified-cat/GP13-EUK.515Y-926R.not-matching.ASVs.aligned.R1.class.cat.tsv",
		#R2missing="{outdir}/12-MG-not-in-ASVdb-classified-cat/GP13-EUK.515Y-926R.not-matching.ASVs.aligned.R2.class.cat.tsv",
		R1andR2missing="{outdir}/12-MG-not-in-ASVdb-classified-cat/GP13-EUK.515Y-926R.not-matching.ASVs.aligned.R1andR2.class.cat.tsv",
		#R1matching="{outdir}/12-MG-not-in-ASVdb-classified-cat/GP13-EUK.515Y-926R.matching.ASVs.aligned.R1.class.cat.tsv",
		#R2matching="{outdir}/12-MG-not-in-ASVdb-classified-cat/GP13-EUK.515Y-926R.matching.ASVs.aligned.R2.class.cat.tsv",
		R1andR2matching="{outdir}/12-MG-not-in-ASVdb-classified-cat/GP13-EUK.515Y-926R.matching.ASVs.aligned.R1andR2.class.cat.tsv"
	shell:
		#"find {outdir}/11-MG-not-in-ASVdb-classified/ -type f -name "
		#"\"*not-matching*.R1.class.tsv\" -print0 | "
		#"xargs -0 cat > {output.R1missing} ; "
		#"find {outdir}/11-MG-not-in-ASVdb-classified/ -type f -name "
		#"\"*not-matching*.R2.class.tsv\" -print0 | "
		#"xargs -0 cat > {output.R2missing} ; "
		"find {outdir}/11-MG-not-in-ASVdb-classified/ -type f -name "
		"\"*not-matching*.R1andR2.class.tsv\" -print0 | "
		"xargs -0 cat > {output.R1andR2missing} ; "
		#"find {outdir}/11-MG-not-in-ASVdb-classified/ -type f -name "
		#"\"*.matching*.R1.class.tsv\" -print0 | "
		#"xargs -0 cat > {output.R1matching} ; "
		#"find {outdir}/11-MG-not-in-ASVdb-classified/ -type f -name "
		#"\"*.matching*.R2.class.tsv\" -print0 | "
		#"xargs -0 cat > {output.R2matching} ; "
		"find {outdir}/11-MG-not-in-ASVdb-classified/ -type f -name "
		"\"*.matching*.R1andR2.class.tsv\" -print0 | "
		"xargs -0 cat > {output.R1andR2matching} "


#Counting order-level groupings (can adjust level with the "cut -d, -f1-4" parameter below)
rule count_tax_matches_and_mismatches:
	input:
		#R1missing="{outdir}/12-MG-not-in-ASVdb-classified-cat/GP13-EUK.515Y-926R.not-matching.ASVs.aligned.R1.class.cat.tsv",
		#R2missing="{outdir}/12-MG-not-in-ASVdb-classified-cat/GP13-EUK.515Y-926R.not-matching.ASVs.aligned.R2.class.cat.tsv",
		R1andR2missing="{outdir}/12-MG-not-in-ASVdb-classified-cat/GP13-EUK.515Y-926R.not-matching.ASVs.aligned.R1andR2.class.cat.tsv",
		#R1matching="{outdir}/12-MG-not-in-ASVdb-classified-cat/GP13-EUK.515Y-926R.matching.ASVs.aligned.R1.class.cat.tsv",
		#R2matching="{outdir}/12-MG-not-in-ASVdb-classified-cat/GP13-EUK.515Y-926R.matching.ASVs.aligned.R2.class.cat.tsv",
		R1andR2matching="{outdir}/12-MG-not-in-ASVdb-classified-cat/GP13-EUK.515Y-926R.matching.ASVs.aligned.R1andR2.class.cat.tsv"
	output:
		#R1missingcounts="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-EUK.515Y-926R.not-matching.ASVs.aligned.R1.class.cat.counts.tsv",
		#R2missingcounts="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-EUK.515Y-926R.not-matching.ASVs.aligned.R2.class.cat.counts.tsv",
		R1andR2missingcounts="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-EUK.515Y-926R.not-matching.ASVs.aligned.R1andR2.class.cat.counts.tsv",
		#R1missingtaxtable="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-EUK.515Y-926R.not-matching.ASVs.aligned.R1.class.cat.taxtable.tsv",
		#R2missingtaxtable="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-EUK.515Y-926R.not-matching.ASVs.aligned.R2.class.cat.taxtable.tsv",
		R1andR2missingtaxtable="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-EUK.515Y-926R.not-matching.ASVs.aligned.R1andR2.class.cat.taxtable.tsv",
		#R1matchingcounts="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-EUK.515Y-926R.matching.ASVs.aligned.R1.class.cat.counts.tsv",
		#R2matchingcounts="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-EUK.515Y-926R.matching.ASVs.aligned.R2.class.cat.counts.tsv",
		R1andR2matchingcounts="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-EUK.515Y-926R.matching.ASVs.aligned.R1andR2.class.cat.counts.tsv",
		#R1matchingtaxtable="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-EUK.515Y-926R.matching.ASVs.aligned.R1.class.cat.taxtable.tsv",
		#R2matchingtaxtable="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-EUK.515Y-926R.matching.ASVs.aligned.R2.class.cat.taxtable.tsv",
		R1andR2matchingtaxtable="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-EUK.515Y-926R.matching.ASVs.aligned.R1andR2.class.cat.taxtable.tsv"
	shell:
		#"sed -re 's/\([0-9]{{1}}\.[0-9]{{2}}\)//g' {input.R1missing} | tee {output.R1missingtaxtable} |" #Remove confidence estimations from VSEARCH output, keep a copy for later steps but also pipe to subsequent commands
		#"cut -f2 | sort | cut -d, -f1-4 | sort | uniq -c | " #Take only tax column, collapse to order level, then count unique occurrences
		#"tail -f -n +2 | awk '{{print $1,\"\t\",$2}}' > {output.R1missingcounts} ; " #Process output into tsv format to stdout
		#"sed -re 's/\([0-9]{{1}}\.[0-9]{{2}}\)//g' {input.R2missing} | tee {output.R2missingtaxtable} |" #Remove confidence estimations from VSEARCH output, keep a copy for later steps but also pipe to subsequent commands
		#"cut -f2 | sort | cut -d, -f1-4 | sort | uniq -c | " #Take only tax column, collapse to order level, then count unique occurrences
		#"tail -f -n +2 | awk '{{print $1,\"\t\",$2}}' > {output.R2missingcounts} ; " #Process output into tsv format to stdout
		"sed -re 's/\([0-9]{{1}}\.[0-9]{{2}}\)//g' {input.R1andR2missing} | tee {output.R1andR2missingtaxtable} |" #Remove confidence estimations from VSEARCH output, keep a copy for later steps but also pipe to subsequent commands
		"cut -f2 | sort | cut -d, -f1-4 | sort | uniq -c | " #Take only tax column, collapse to order level, then count unique occurrences
		"tail -f -n +2 | awk '{{print $1,\"\t\",$2}}' > {output.R1andR2missingcounts} ; " #Process output into tsv format to stdout
		#"sed -re 's/\([0-9]{{1}}\.[0-9]{{2}}\)//g' {input.R1matching} | tee {output.R1matchingtaxtable} |" #Remove confidence estimations from VSEARCH output, keep a copy for later steps but also pipe to subsequent commands
		#"cut -f2 | sort | cut -d, -f1-4 | sort | uniq -c | " #Take only tax column, collapse to order level, then count unique occurrences
		#"tail -f -n +2 | awk '{{print $1,\"\t\",$2}}' > {output.R1matchingcounts} ; " #Process output into tsv format to stdout
		#"sed -re 's/\([0-9]{{1}}\.[0-9]{{2}}\)//g' {input.R2matching} | tee {output.R2matchingtaxtable} |" #Remove confidence estimations from VSEARCH output, keep a copy for later steps but also pipe to subsequent commands
		#"cut -f2 | sort | cut -d, -f1-4 | sort | uniq -c | " #Take only tax column, collapse to order level, then count unique occurrences
		#"tail -f -n +2 | awk '{{print $1,\"\t\",$2}}' > {output.R2matchingcounts} ; " #Process output into tsv format to stdout
		"sed -re 's/\([0-9]{{1}}\.[0-9]{{2}}\)//g' {input.R1andR2matching} | tee {output.R1andR2matchingtaxtable} |" #Remove confidence estimations from VSEARCH output, keep a copy for later steps but also pipe to subsequent commands
		"cut -f2 | sort | cut -d, -f1-4 | sort | uniq -c | " #Take only tax column, collapse to order level, then count unique occurrences
		"tail -f -n +2 | awk '{{print $1,\"\t\",$2}}' > {output.R1andR2matchingcounts} " #Process output into tsv format to stdout

#Now take only those with greater than 1 % abundance (among mismatches) using basic python script (can change abundance cutoff if you desire)
rule filter_tax_matches_by_abundance:
	input:
		#R1missing="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-EUK.515Y-926R.not-matching.ASVs.aligned.R1.class.cat.counts.tsv",
		#R2missing="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-EUK.515Y-926R.not-matching.ASVs.aligned.R2.class.cat.counts.tsv",
		R1andR2missing="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-EUK.515Y-926R.not-matching.ASVs.aligned.R1andR2.class.cat.counts.tsv",
		#R1matching="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-EUK.515Y-926R.matching.ASVs.aligned.R1.class.cat.counts.tsv",
		#R2matching="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-EUK.515Y-926R.matching.ASVs.aligned.R2.class.cat.counts.tsv",
		R1andR2matching="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-EUK.515Y-926R.matching.ASVs.aligned.R1andR2.class.cat.counts.tsv"
	params:
		minAbund = 0.001  #Change fractional value in config file if desired, default 0.01
	output:
		#R1missing="{outdir}/14-MG-not-in-ASVdb-classified-cat-parsed-fractions/GP13-EUK.515Y-926R.not-matching.ASVs.aligned.R1.class.cat.counts.frac.min" + str(cutoff) + ".tsv",
		#R2missing="{outdir}/14-MG-not-in-ASVdb-classified-cat-parsed-fractions/GP13-EUK.515Y-926R.not-matching.ASVs.aligned.R2.class.cat.counts.frac.min" + str(cutoff) + ".tsv",
		R1andR2missing="{outdir}/14-MG-not-in-ASVdb-classified-cat-parsed-fractions/GP13-EUK.515Y-926R.not-matching.ASVs.aligned.R1andR2.class.cat.counts.frac.min" + str(cutoff) + ".tsv",
		#R1matching="{outdir}/14-MG-not-in-ASVdb-classified-cat-parsed-fractions/GP13-EUK.515Y-926R.matching.ASVs.aligned.R1.class.cat.counts.frac.min" + str(cutoff) + ".tsv",
		#R2matching="{outdir}/14-MG-not-in-ASVdb-classified-cat-parsed-fractions/GP13-EUK.515Y-926R.matching.ASVs.aligned.R2.class.cat.counts.frac.min" + str(cutoff) + ".tsv",
		R1andR2matching="{outdir}/14-MG-not-in-ASVdb-classified-cat-parsed-fractions/GP13-EUK.515Y-926R.matching.ASVs.aligned.R1andR2.class.cat.counts.frac.min" + str(cutoff) + ".tsv"
	shell:
		#"scripts/mismatch-characterization/filter-by-fractional-abundance.py {input.R1missing} {params.minAbund} > {output.R1missing} ; "
		#"scripts/mismatch-characterization/filter-by-fractional-abundance.py {input.R2missing} {params.minAbund} > {output.R2missing} ; "
		"scripts/mismatch-characterization/filter-by-fractional-abundance.py {input.R1andR2missing} {params.minAbund} > {output.R1andR2missing} ; "
		#"scripts/mismatch-characterization/filter-by-fractional-abundance.py {input.R1matching} {params.minAbund} > {output.R1matching} ; "
		#"scripts/mismatch-characterization/filter-by-fractional-abundance.py {input.R2matching} {params.minAbund} > {output.R2matching} ; "
		"scripts/mismatch-characterization/filter-by-fractional-abundance.py {input.R1andR2matching} {params.minAbund} > {output.R1andR2matching} "

rule make_kronas:
	input:
		#R1missing="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-EUK.515Y-926R.not-matching.ASVs.aligned.R1.class.cat.counts.tsv",
		#R2missing="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-EUK.515Y-926R.not-matching.ASVs.aligned.R2.class.cat.counts.tsv",
		R1andR2missing="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-EUK.515Y-926R.not-matching.ASVs.aligned.R1andR2.class.cat.counts.tsv",
		#R1matching="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-EUK.515Y-926R.matching.ASVs.aligned.R1.class.cat.counts.tsv",
		#R2matching="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-EUK.515Y-926R.matching.ASVs.aligned.R2.class.cat.counts.tsv",
		R1andR2matching="{outdir}/13-MG-not-in-ASVdb-classified-cat-parsed/GP13-EUK.515Y-926R.matching.ASVs.aligned.R1andR2.class.cat.counts.tsv"
	output:
		#R1missing="{outdir}/15-MG-not-in-ASVdb-classified-kronas/GP13-EUK.515Y-926R.not-matching.ASVs.aligned.R1.krona.input",
		#R2missing="{outdir}/15-MG-not-in-ASVdb-classified-kronas/GP13-EUK.515Y-926R.not-matching.ASVs.aligned.R2.krona.input",
		R1andR2missing="{outdir}/15-MG-not-in-ASVdb-classified-kronas/GP13-EUK.515Y-926R.not-matching.ASVs.aligned.R1andR2.krona.input",
		#R1htmlmissing="{outdir}/15-MG-not-in-ASVdb-classified-kronas/GP13-EUK.515Y-926R.not-matching.ASVs.aligned.R1.krona.html",
		#R2htmlmissing="{outdir}/15-MG-not-in-ASVdb-classified-kronas/GP13-EUK.515Y-926R.not-matching.ASVs.aligned.R2.krona.html",
		R1andR2htmlmissing="{outdir}/15-MG-not-in-ASVdb-classified-kronas/GP13-EUK.515Y-926R.not-matching.ASVs.aligned.R1andR2.krona.html",
		#R1matching="{outdir}/15-MG-not-in-ASVdb-classified-kronas/GP13-EUK.515Y-926R.matching.ASVs.aligned.R1.krona.input",
		#R2matching="{outdir}/15-MG-not-in-ASVdb-classified-kronas/GP13-EUK.515Y-926R.matching.ASVs.aligned.R2.krona.input",
		R1andR2matching="{outdir}/15-MG-not-in-ASVdb-classified-kronas/GP13-EUK.515Y-926R.matching.ASVs.aligned.R1andR2.krona.input",
		#R1htmlmatching="{outdir}/15-MG-not-in-ASVdb-classified-kronas/GP13-EUK.515Y-926R.matching.ASVs.aligned.R1.krona.html",
		#R2htmlmatching="{outdir}/15-MG-not-in-ASVdb-classified-kronas/GP13-EUK.515Y-926R.matching.ASVs.aligned.R2.krona.html",
		R1andR2htmlmatching="{outdir}/15-MG-not-in-ASVdb-classified-kronas/GP13-EUK.515Y-926R.matching.ASVs.aligned.R1andR2.krona.html"
	conda:
		"envs/krona.yaml"
	shell:
		#"sed 's/,/\t/g' {input.R1missing} | sed 's/\w://g' > {output.R1missing} ; ktImportText -c -n \"R1 reads\" -o {output.R1htmlmissing} {output.R1missing} ; "
		#"sed 's/,/\t/g' {input.R2missing} | sed 's/\w://g' > {output.R2missing} ; ktImportText -c -n \"R2 reads\" -o {output.R2htmlmissing} {output.R2missing} ; "
		"sed 's/,/\t/g' {input.R1andR2missing} | sed 's/\w://g' > {output.R1andR2missing} ; ktImportText -c -n \"R1 and R2 reads\" -o {output.R1andR2htmlmissing} {output.R1andR2missing} ; "
		#"sed 's/,/\t/g' {input.R1matching} | sed 's/\w://g' > {output.R1matching} ; ktImportText -c -n \"R1 reads\" -o {output.R1htmlmatching} {output.R1matching} ; "
		#"sed 's/,/\t/g' {input.R2matching} | sed 's/\w://g' > {output.R2matching} ; ktImportText -c -n \"R2 reads\" -o {output.R2htmlmatching} {output.R2matching} ; "
		"sed 's/,/\t/g' {input.R1andR2matching} | sed 's/\w://g' > {output.R1andR2matching} ; ktImportText -c -n \"R1 and R2 reads\" -o {output.R1andR2htmlmatching} {output.R1andR2matching} "
