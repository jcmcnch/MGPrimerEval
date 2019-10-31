rule all:
	input:
	    expand("compute-workflow-intermediate/06-subsetted/{sample}.SSU.{direction}.{group}_pyNAST_{primer}.fasta", sample=config["samples"], group=config["groups"], primer=config["primer"], direction=['fwd','rev'])
		#expand("compute-workflow-intermediate/02-phyloFlash_sifted/{sample}.fwd.SSU.hits.fastq", sample=config["samples"])
		#expand("10-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.info", sample=config["samples"], group=config["groups"], primer=config["primer"], mismatches=config["mismatches"], direction=['fwd','rev']),
		#expand("11-summary/{sample}.{direction}.{group}.{primer}.{mismatches}.summary.tsv", sample=config["samples"], group=config["groups"], primer=config["primer"], mismatches=config["mismatches"], direction=['fwd','rev'])

rule fastp_clean:
	input:
		R1="00-fastq/{sample}_repaired_1.fastq.gz",
		R2="00-fastq/{sample}_repaired_2.fastq.gz"
	output:
		R1clean="compute-workflow-intermediate/01-fastp-cleaned/{sample}_1_clean.fastq.gz",
		R2clean="compute-workflow-intermediate/01-fastp-cleaned/{sample}_2_clean.fastq.gz",
		log="compute-workflow-intermediate/logs/01-fastp_cleaning/{sample}.log.html"
	threads:
		8
	shell:
		"fastp -3 -t 3 -i {input.R1} -I {input.R2} -o {output.R1clean} -O {output.R2clean} --thread {threads} --html {output.log}"

rule phyloFlash_sift:
	input:
		R1clean="compute-workflow-intermediate/01-fastp-cleaned/{sample}_1_clean.fastq.gz",
		R2clean="compute-workflow-intermediate/01-fastp-cleaned/{sample}_2_clean.fastq.gz"
	output:
		R1hits="compute-workflow-intermediate/02-phyloFlash_sifted/{sample}.fwd.SSU.hits.fastq",
		R2hits="compute-workflow-intermediate/02-phyloFlash_sifted/{sample}.rev.SSU.hits.fastq"
	params:
		workdir="compute-workflow-intermediate/tmp/",
		phyloFlash_other="compute-workflow-intermediate/02-phyloFlash_sifted/phyloFlash-other/"
	threads:
		8
	conda:
		"envs/phyloflash-env.yaml"
	shell:
		"""
		mkdir -p {params.phyloFlash_other} ; mkdir -p {params.workdir} ; cd {params.workdir} ;
		phyloFlash.pl -lib {wildcards.sample} -read1 ../../{input.R1clean} -read2 ../../{input.R2clean} -dbhome /home/db/phyloFlash/132 -readlength 144 -id 50 -CPUs {threads} -log -skip_spades -nozip ;
		mv {wildcards.sample}.`basename {input.R1clean}`.SSU.1.fq ../../{output.R1hits} ;
		mv {wildcards.sample}.`basename {input.R1clean}`.SSU.2.fq ../../{output.R2hits} ;
		mv {wildcards.sample}* ../../{params.phyloFlash_other}
		"""

rule remove_repeats_komplexity:
	input:
		R1hits="compute-workflow-intermediate/02-phyloFlash_sifted/{sample}.fwd.SSU.hits.fastq",
		R2hits="compute-workflow-intermediate/02-phyloFlash_sifted/{sample}.rev.SSU.hits.fastq"
	output:
		R1="compute-workflow-intermediate/03-low-complexity-filtered/{sample}.fwd.SSU.keep.fastq",
		R2="compute-workflow-intermediate/03-low-complexity-filtered/{sample}.rev.SSU.keep.fastq"
	conda:
		"envs/komplexity.yaml"
	shell:
		"kz --filter < {input.R1hits} > {output.R1} ; "
		"kz --filter < {input.R2hits} > {output.R2} "

rule sort_EUK:
	input:
		"compute-workflow-intermediate/03-low-complexity-filtered/{sample}.{direction}.SSU.keep.fastq"
	output:
		PROK="compute-workflow-intermediate/04-sorted/{sample}.{direction}.SSU.PROK.fastq",
		EUK="compute-workflow-intermediate/04-sorted/{sample}.{direction}.SSU.EUK.fastq"
	threads:
		8
	conda:
		"envs/bbmap.yaml"
	params:
		workdir="compute-workflow-intermediate/tmp/"
	shell:
		"cd {params.workdir} ; "
		"bbsplit.sh threads={threads} -Xmx100g overwrite=t usequality=f qtrim=f minratio=0.30 minid=0.30 pairedonly=f "
		"path=/home/db/bbsplit-db/EUK-PROK-bbsplit-db/ "
		"in=../../{input} basename={input}_%.fastq ; "
		"mv {input}*EUK*.fastq ../../{output.EUK}; mv {input}*PROK*.fastq ../../{output.PROK}"


rule sort_PROK:
	input:
		"compute-workflow-intermediate/04-sorted/{sample}.{direction}.SSU.PROK.fastq"
	output:
		ARCH="compute-workflow-intermediate/04-sorted/{sample}.{direction}.SSU.ARCH.fastq",
		BACT="compute-workflow-intermediate/04-sorted/{sample}.{direction}.SSU.BACT.fastq"
	threads:
		8
	conda:
		"envs/bbmap.yaml"
	params:
		workdir="compute-workflow-intermediate/tmp/"
	shell:
		"cd {params.workdir} ; "
		"bbsplit.sh threads={threads} -Xmx100g overwrite=t usequality=f qtrim=f minratio=0.30 minid=0.30 pairedonly=f "
		"path=/home/db/bbsplit-db/BACT-ARCH-bbsplit-db/ "
		"in=../../{input} basename={input}_%.fastq ; "
		"mv {input}*ARCH*.fastq ../../{output.ARCH}; mv {input}*BACT*.fastq ../../{output.BACT}"


rule sort_CYANO:
	input:
		"compute-workflow-intermediate/04-sorted/{sample}.{direction}.SSU.BACT.fastq"
	output:
		CYANO="compute-workflow-intermediate/04-sorted/{sample}.{direction}.SSU.BACT-CYANO.fastq",
		NONCYANO="compute-workflow-intermediate/04-sorted/{sample}.{direction}.SSU.BACT-NON-CYANO.fastq"
	threads:
		8
	conda:
		"envs/bbmap.yaml"
	params:
		workdir="compute-workflow-intermediate/tmp/"
	shell:
		"cd {params.workdir} ; "
		"bbsplit.sh threads={threads} -Xmx100g overwrite=t usequality=f qtrim=f minratio=0.30 minid=0.30 pairedonly=f "
		"path=/home/db/bbsplit-db/BACT-CYANO-bbsplit-db/ "
		"in=../../{input} basename={input}_%.fastq ; "
		"mv {input}*NON-CYANO*.fastq ../../{output.NONCYANO}; mv {input}*CYANO*.fastq ../../{output.CYANO}"


rule align_ARCH:
	input:
		seqs="compute-workflow-intermediate/04-sorted/{sample}.{direction}.SSU.ARCH.fastq",
		ref="SSU_refs/Sulfolobus_acidocaldarius_N8_16s.fasta"
	output:
		aligned="compute-workflow-intermediate/05-pyNAST-aligned/{sample}.{direction}.SSU.ARCH_pynast_aligned.fasta",
		log="compute-workflow-intermediate/05-pyNAST-aligned/{sample}.{direction}.SSU.ARCH.pyNAST.log"
	conda:
		"envs/pynast.yaml"
	shell:
		"tmpfile=`mktemp /tmp/fastq-ids.XXXXXXXXXXXXXXXX` ; seqtk seq -A {input.seqs} > $tmpfile ; "
		"pynast -p 10 -l 1 -i $tmpfile -t {input.ref} -a {output.aligned} -g {output.log} ; "
		"rm $tmpfile"


rule align_BACT:
	input:
		seqs="compute-workflow-intermediate/04-sorted/{sample}.{direction}.SSU.BACT-NON-CYANO.fastq",
		ref="SSU_refs/Ecoli_16s.fna"
	output:
		aligned="compute-workflow-intermediate/05-pyNAST-aligned/{sample}.{direction}.SSU.BACT-NON-CYANO_pynast_aligned.fasta",
		log="compute-workflow-intermediate/05-pyNAST-aligned/{sample}.{direction}.SSU.BACT-NON-CYANO.pyNAST.log"
	conda:
		"envs/pynast.yaml"
	shell:
		"tmpfile=`mktemp /tmp/fastq-ids.XXXXXXXXXXXXXXXX` ; seqtk seq -A {input.seqs} > $tmpfile ; "
		"pynast -p 10 -l 1 -i $tmpfile -t {input.ref} -a {output.aligned} -g {output.log} ; "
		"rm $tmpfile"


rule align_CYANO:
	input:
		seqs="compute-workflow-intermediate/04-sorted/{sample}.{direction}.SSU.BACT-CYANO.fastq",
		ref="SSU_refs/longest-CYANO-with-27F.fasta"
	output:
		aligned="compute-workflow-intermediate/05-pyNAST-aligned/{sample}.{direction}.SSU.BACT-CYANO_pynast_aligned.fasta",
		log="compute-workflow-intermediate/05-pyNAST-aligned/{sample}.{direction}.SSU.BACT-CYANO.pyNAST.log"
	conda:
		"envs/pynast.yaml"
	shell:
		"tmpfile=`mktemp /tmp/fastq-ids.XXXXXXXXXXXXXXXX` ; seqtk seq -A {input.seqs} > $tmpfile ; "
		"pynast -p 10 -l 1 -i $tmpfile -t {input.ref} -a {output.aligned} -g {output.log} ; "
		"rm $tmpfile"


rule align_EUK:
	input:
		seqs="compute-workflow-intermediate/04-sorted/{sample}.{direction}.SSU.EUK.fastq",
		ref="SSU_refs/Saccharomyces_cerevisiae_S288C_18s-1_NR_132213.1.fa"
	output:
		aligned="compute-workflow-intermediate/05-pyNAST-aligned/{sample}.{direction}.SSU.EUK_pynast_aligned.fasta",
		log="compute-workflow-intermediate/05-pyNAST-aligned/{sample}.{direction}.SSU.EUK.pyNAST.log"
	conda:
		"envs/pynast.yaml"
	shell:
		"tmpfile=`mktemp /tmp/fastq-ids.XXXXXXXXXXXXXXXX` ; seqtk seq -A {input.seqs} > $tmpfile ; "
		"pynast -p 10 -l 1 -i $tmpfile -t {input.ref} -a {output.aligned} -g {output.log} ; "
		"rm $tmpfile"

rule subset_to_primer_region:
	input:
		"compute-workflow-intermediate/05-pyNAST-aligned/{sample}.{direction}.SSU.{group}_pynast_aligned.fasta"
	output:
		"compute-workflow-intermediate/06-subsetted/{sample}.SSU.{direction}.{group}_pyNAST_{primer}.fasta"
	conda:
		"envs/biopython.yaml"
	params:
		start=lambda wildcards: config["primerROI"][wildcards.group][wildcards.primer][0],
		end=lambda wildcards: config["primerROI"][wildcards.group][wildcards.primer][1]
	shell:
		"scripts/filter-pyNAST-for-ROI.py --input {input} --output {output} --start {params.start} --end {params.end}"

"""
Since graftM actually told us which strand the SSU rRNA is on, then we need to add in a rule that infers this from the alignment.
"""

#If you have many samples and/or primers, it may make sense to implement the --batch flag here
rule get_fastq_for_subset:
	#Get the whole fastq files associated with the matching reads
	input:
		fasta="compute-workflow-intermediate/06-subsetted/{sample}.SSU.{direction}.{group}_pyNAST_{primer}.fasta",
		fastq="compute-workflow-intermediate/04-sorted/{sample}.{direction}.SSU.{group}.fastq"
	output:
		fastq=temp("compute-workflow-intermediate/06-subsetted/{sample}.SSU.{direction}.{group}_pyNAST_{primer}.full.fastq"),
		fastq_revcomp="compute-workflow-intermediate/07-subsetted-fastq/"
	conda:
		"envs/bbmap.yaml"
	shell:
		"filterbyname.sh names={input.fasta} include=t in={input.fastq} out={output.fastq} ; "
		"revcompfastq_according_to_pyNAST.py --inpynast {input.fasta} --infastq {output.fastq} --outfastq {output.fastq_revcomp}"


rule grab_matching_cutadapt_full:
	input:
		"09-complemented-fastqs/{sample}.SSU.{direction}.{group}_pyNAST_{primer}.full.revcomped.fastq"
	output:
		mismatch="10-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.nohit.fastq",
		match="10-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.hit.fastq",
		info="10-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.info"
	log:
		"logs/09-cutadapt/{sample}.{direction}.{group}.{primer}.{mismatches}.cutadapt.log"
	params:
		pattern=lambda wildcards : config["primer"][wildcards.primer],
		errorRate=lambda wildcards : config["mismatches"][wildcards.mismatches] / len(config["primer"][wildcards.primer]),
		lengthPrimer=lambda wildcards : len(config["primer"][wildcards.primer])
	shell:
		"cutadapt -f fastq --info-file={output.info} --no-indels --no-trim --overlap={params.lengthPrimer} -b {params.pattern} --error-rate={params.errorRate} --untrimmed-output={output.mismatch} --output={output.match} {input}"


rule quality_filter_primer_region:
    #Keep only sequences that have >30 phred score across the whole primer + 5 leading/trailing bases (implicit in script)
	input:
		mismatch="10-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.nohit.fastq",
		match="10-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.hit.fastq",
		info="10-checked/{primer}/6-mismatch/{sample}.SSU.{direction}.{group}.{primer}.6-mismatch.info" #Assume anything with -1 value is false positive, ignoring those with > ~30% mismatches to primer (e.g. for a 20bp primer)
	output:
		mismatch="10-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.nohit.filtered.fastq",
		match="10-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.hit.filtered.fastq"
	conda:
		"envs/biopython.yaml"
	shell:
		"scripts/filter-fastq-by-primer-ROI-coverage.py --info {input.info} --fastq {input.match} > {output.match} ; "
		"scripts/filter-fastq-by-primer-ROI-coverage.py --info {input.info} --fastq {input.mismatch} > {output.mismatch} "


rule compute_percentages:
	input:
		mismatch="10-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.nohit.filtered.fastq",
		match="10-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.hit.filtered.fastq"
	output:
		"11-summary/{sample}.{direction}.{group}.{primer}.{mismatches}.summary.tsv"
	shell:
		"numMatch=`grep -cE \"^@\" {input.match} || numMatch=0` ; numMismatch=`grep -cE \"^@\" {input.mismatch}` || numMismatch=0 ; "
		"sumTotal=`expr $numMatch + $numMismatch || sumTotal=0` ; "
		"if `[ $sumTotal -eq 0 ]` ; then fracMatch=\"NA\"; elif `[ $numMatch -ne 0 ]` && `[ $numMismatch -eq 0 ]`; then fracMatch=1; elif `[ $numMatch -eq 0 ]` && `[ $numMismatch -ne 0 ]`; then fracMatch=0; elif `[ $numMismatch -ne 0 ]` && `[ $numMatch -ne 0 ]`; then fracMatch=`bc <<< \"scale=4; $numMatch/$sumTotal\"`; fi ; "
		"printf \"{wildcards.sample}\t{wildcards.direction}\t{wildcards.group}\t{wildcards.primer}\t{wildcards.mismatches}\t$sumTotal\t$numMatch\t$numMismatch\t$fracMatch\n\" > {output}"
