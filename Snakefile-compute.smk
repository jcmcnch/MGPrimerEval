#readLength=config["readLength"] #no longer needed for phyloFlash 3.4
readlimit=config["readlimit"]
suffixR1=config["suffixR1"]
suffixR2=config["suffixR2"]
phyloFlashDB=config["phyloFlashDB"]
bbsplitDBpath=config["bbsplitDBpath"]
#Prefix to add uclust executable to path
uclustpath=config["uclustpath"]
shell.prefix('PATH=' + uclustpath + ':$PATH ')

rule all:
	input:
	    expand("intermediate/compute-workflow/06-subsetted/{sample}.SSU.{direction}.{group}_pyNAST_{primer}.fasta", sample=config["samples"], group=config["groups"], primer=config["primer"], direction=['fwd','rev']),
		expand("intermediate/compute-workflow/07-subsetted-fastq/{sample}.SSU.{direction}.{group}_pyNAST_{primer}.full.revcomped.fastq", sample=config["samples"], group=config["groups"], primer=config["primer"], direction=['fwd','rev']),
		expand("intermediate/compute-workflow/08-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.info", sample=config["samples"], group=config["groups"], primer=config["primer"], mismatches=config["mismatches"], direction=['fwd','rev']),
		expand("output/compute-workflow/09-summary/{sample}.{direction}.{group}.{primer}.{mismatches}.summary.tsv", sample=config["samples"], group=config["groups"], primer=config["primer"], mismatches=config["mismatches"], direction=['fwd','rev'])

rule fastp_clean:
	input:
		R1="intermediate/compute-workflow/00-fastq/{sample}" + suffixR1,
		R2="intermediate/compute-workflow/00-fastq/{sample}" + suffixR2
	output:
		R1clean=temp("intermediate/compute-workflow/01-fastp-cleaned/{sample}_1_clean.fastq.gz"),
		R2clean=temp("intermediate/compute-workflow/01-fastp-cleaned/{sample}_2_clean.fastq.gz"),
		log="intermediate/compute-workflow/logs/01-fastp_cleaning/{sample}.log.html"
	conda:
		"envs/fastp.yaml"
	threads:
		8
	shell:
		"fastp -3 -t 3 -i {input.R1} -I {input.R2} -o {output.R1clean} -O {output.R2clean} --thread {threads} --html {output.log}"

rule phyloFlash_sift:
	input:
		R1clean="intermediate/compute-workflow/01-fastp-cleaned/{sample}_1_clean.fastq.gz",
		R2clean="intermediate/compute-workflow/01-fastp-cleaned/{sample}_2_clean.fastq.gz"
	output:
		R1hits="intermediate/compute-workflow/02-phyloFlash_sifted/{sample}.fwd.SSU.hits.fastq",
		R2hits="intermediate/compute-workflow/02-phyloFlash_sifted/{sample}.rev.SSU.hits.fastq"
	params:
		workdir="intermediate/compute-workflow/tmp/",
		phyloFlash_other="intermediate/compute-workflow/02-phyloFlash_sifted/phyloFlash-other/",
		readlimit=readlimit,
		db=phyloFlashDB
	threads:
		8
	conda:
		"envs/phyloflash-env.yaml"
	shell:
                """
                mkdir -p {params.phyloFlash_other} ; mkdir -p {params.workdir} ; cd {params.workdir} ;
                phyloFlash.pl -readlimit {params.readlimit} -lib {wildcards.sample} -read1 ../../../{input.R1clean} -read2 ../../../{input.R2clean} -dbhome {params.db} -id 50 -CPUs {threads} -log -skip_spades -nozip ;
                mv {wildcards.sample}.`basename {input.R1clean}`.SSU.1.fq ../../../{output.R1hits} ;
                mv {wildcards.sample}.`basename {input.R1clean}`.SSU.2.fq ../../../{output.R2hits} ;
                mv {wildcards.sample}* ../../../{params.phyloFlash_other}
                """

rule remove_repeats_komplexity:
	input:
		R1hits="intermediate/compute-workflow/02-phyloFlash_sifted/{sample}.fwd.SSU.hits.fastq",
		R2hits="intermediate/compute-workflow/02-phyloFlash_sifted/{sample}.rev.SSU.hits.fastq"
	output:
		R1="intermediate/compute-workflow/03-low-complexity-filtered/{sample}.fwd.SSU.keep.fastq",
		R2="intermediate/compute-workflow/03-low-complexity-filtered/{sample}.rev.SSU.keep.fastq"
	conda:
		"envs/komplexity.yaml"
	shell:
		"kz --filter < {input.R1hits} > {output.R1} ; "
		"kz --filter < {input.R2hits} > {output.R2} "

rule sort_EUK:
	input:
		"intermediate/compute-workflow/03-low-complexity-filtered/{sample}.{direction}.SSU.keep.fastq"
	output:
		PROK="intermediate/compute-workflow/04-sorted/{sample}.{direction}.SSU.PROK.fastq",
		EUK="intermediate/compute-workflow/04-sorted/{sample}.{direction}.SSU.EUK.fastq"
	threads:
		8
	conda:
		"envs/bbmap.yaml"
	params:
		workdir="intermediate/compute-workflow/tmp/",
		bbsplitdb=bbsplitDBpath + "/EUK-PROK-bbsplit-db/"
	shell:
		"cd {params.workdir} ; "
		"bbsplit.sh ow=f threads={threads} -Xmx100g usequality=f qtrim=f minratio=0.30 minid=0.30 pairedonly=f "
		"path={params.bbsplitdb} "
		"in=../../../{input} basename={input}_%.fastq ; "
		"mv {input}*EUK*.fastq ../../../{output.EUK}; mv {input}*PROK*.fastq ../../../{output.PROK}"


rule sort_PROK:
	input:
		"intermediate/compute-workflow/04-sorted/{sample}.{direction}.SSU.PROK.fastq"
	output:
		ARCH="intermediate/compute-workflow/04-sorted/{sample}.{direction}.SSU.ARCH.fastq",
		BACT="intermediate/compute-workflow/04-sorted/{sample}.{direction}.SSU.BACT.fastq"
	threads:
		8
	conda:
		"envs/bbmap.yaml"
	params:
		workdir="intermediate/compute-workflow/tmp/",
                bbsplitdb=bbsplitDBpath + "/BACT-ARCH-bbsplit-db/"
	shell:
		"cd {params.workdir} ; "
		"bbsplit.sh ow=f threads={threads} -Xmx100g usequality=f qtrim=f minratio=0.30 minid=0.30 pairedonly=f "
		"path={params.bbsplitdb} "
		"in=../../../{input} basename={input}_%.fastq ; "
		"mv {input}*ARCH*.fastq ../../../{output.ARCH}; mv {input}*BACT*.fastq ../../../{output.BACT}"


rule sort_CYANO:
	input:
		"intermediate/compute-workflow/04-sorted/{sample}.{direction}.SSU.BACT.fastq"
	output:
		CYANO="intermediate/compute-workflow/04-sorted/{sample}.{direction}.SSU.BACT-CYANO.fastq",
		NONCYANO="intermediate/compute-workflow/04-sorted/{sample}.{direction}.SSU.BACT-NON-CYANO.fastq"
	threads:
		8
	conda:
		"envs/bbmap.yaml"
	params:
		workdir="intermediate/compute-workflow/tmp/",
                bbsplitdb=bbsplitDBpath + "/BACT-CYANO-bbsplit-db/"
	shell:
		"cd {params.workdir} ; "
		"bbsplit.sh ow=f threads={threads} -Xmx100g usequality=f qtrim=f minratio=0.30 minid=0.30 pairedonly=f "
		"path={params.bbsplitdb} "
		"in=../../../{input} basename={input}_%.fastq ; "
		"mv {input}*NON-CYANO*.fastq ../../../{output.NONCYANO}; mv {input}*CYANO*.fastq ../../../{output.CYANO}"


rule align_ARCH:
	input:
		seqs="intermediate/compute-workflow/04-sorted/{sample}.{direction}.SSU.ARCH.fastq",
		ref="SSU_refs/Sulfolobus_acidocaldarius_N8_16s.fasta"
	output:
		aligned="intermediate/compute-workflow/05-pyNAST-aligned/{sample}.{direction}.SSU.ARCH_pynast_aligned.fasta",
		log="intermediate/compute-workflow/05-pyNAST-aligned/{sample}.{direction}.SSU.ARCH.pyNAST.log"
	conda:
		"envs/pynast.yaml"
	shell:
		"tmpfile=`mktemp /tmp/fastq-ids.XXXXXXXXXXXXXXXX` ; seqtk seq -A {input.seqs} > $tmpfile ; "
		"pynast -p 10 -l 1 -i $tmpfile -t {input.ref} -a {output.aligned} -g {output.log} ; "
		"rm $tmpfile"


rule align_BACT:
	input:
		seqs="intermediate/compute-workflow/04-sorted/{sample}.{direction}.SSU.BACT-NON-CYANO.fastq",
		ref="SSU_refs/Ecoli_16s.fna"
	output:
		aligned="intermediate/compute-workflow/05-pyNAST-aligned/{sample}.{direction}.SSU.BACT-NON-CYANO_pynast_aligned.fasta",
		log="intermediate/compute-workflow/05-pyNAST-aligned/{sample}.{direction}.SSU.BACT-NON-CYANO.pyNAST.log"
	conda:
		"envs/pynast.yaml"
	shell:
		"tmpfile=`mktemp /tmp/fastq-ids.XXXXXXXXXXXXXXXX` ; seqtk seq -A {input.seqs} > $tmpfile ; "
		"pynast -p 10 -l 1 -i $tmpfile -t {input.ref} -a {output.aligned} -g {output.log} ; "
		"rm $tmpfile"


rule align_CYANO:
	input:
		seqs="intermediate/compute-workflow/04-sorted/{sample}.{direction}.SSU.BACT-CYANO.fastq",
		ref="SSU_refs/longest-CYANO-with-27F.fasta"
	output:
		aligned="intermediate/compute-workflow/05-pyNAST-aligned/{sample}.{direction}.SSU.BACT-CYANO_pynast_aligned.fasta",
		log="intermediate/compute-workflow/05-pyNAST-aligned/{sample}.{direction}.SSU.BACT-CYANO.pyNAST.log"
	conda:
		"envs/pynast.yaml"
	shell:
		"tmpfile=`mktemp /tmp/fastq-ids.XXXXXXXXXXXXXXXX` ; seqtk seq -A {input.seqs} > $tmpfile ; "
		"pynast -p 10 -l 1 -i $tmpfile -t {input.ref} -a {output.aligned} -g {output.log} ; "
		"rm $tmpfile"


rule align_EUK:
	input:
		seqs="intermediate/compute-workflow/04-sorted/{sample}.{direction}.SSU.EUK.fastq",
		ref="SSU_refs/Saccharomyces_cerevisiae_S288C_18s-1_NR_132213.1.fa"
	output:
		aligned="intermediate/compute-workflow/05-pyNAST-aligned/{sample}.{direction}.SSU.EUK_pynast_aligned.fasta",
		log="intermediate/compute-workflow/05-pyNAST-aligned/{sample}.{direction}.SSU.EUK.pyNAST.log"
	conda:
		"envs/pynast.yaml"
	shell:
		"tmpfile=`mktemp /tmp/fastq-ids.XXXXXXXXXXXXXXXX` ; seqtk seq -A {input.seqs} > $tmpfile ; "
		"pynast -p 10 -l 1 -i $tmpfile -t {input.ref} -a {output.aligned} -g {output.log} ; "
		"rm $tmpfile"

rule subset_to_primer_region:
	input:
		"intermediate/compute-workflow/05-pyNAST-aligned/{sample}.{direction}.SSU.{group}_pynast_aligned.fasta"
	output:
		"intermediate/compute-workflow/06-subsetted/{sample}.SSU.{direction}.{group}_pyNAST_{primer}.fasta"
	conda:
		"envs/biopython.yaml"
	params:
		start=lambda wildcards: config["primerROI"][wildcards.group][wildcards.primer][0],
		end=lambda wildcards: config["primerROI"][wildcards.group][wildcards.primer][1]
	shell:
		"scripts/filter-pyNAST-for-ROI.py --input {input} --output {output} --start {params.start} --end {params.end}"

#If you have many samples and/or primers, it probably makes sense to implement the --batch flag here
rule get_fastq_for_subset:
	#Get the whole fastq files associated with the matching reads
	input:
		fasta="intermediate/compute-workflow/06-subsetted/{sample}.SSU.{direction}.{group}_pyNAST_{primer}.fasta",
		fastq="intermediate/compute-workflow/04-sorted/{sample}.{direction}.SSU.{group}.fastq"
	output:
		fastq="intermediate/compute-workflow/07-subsetted-fastq/{sample}.SSU.{direction}.{group}_pyNAST_{primer}.full.fastq",
	conda:
		"envs/bbmap.yaml"
	shell:
		#substring=t must be added because pyNAST adds additional information to sequence header
		"filterbyname.sh substring=t names={input.fasta} include=t in={input.fastq} out={output.fastq} "

"""
Since graftM actually told us which strand the SSU rRNA is on, then we need to add in a rule that infers this from the alignment.
"""

rule revcomp_fastqs:
	input:
		fasta="intermediate/compute-workflow/06-subsetted/{sample}.SSU.{direction}.{group}_pyNAST_{primer}.fasta",
		fastq="intermediate/compute-workflow/07-subsetted-fastq/{sample}.SSU.{direction}.{group}_pyNAST_{primer}.full.fastq"
	output:
		fastq_revcomp="intermediate/compute-workflow/07-subsetted-fastq/{sample}.SSU.{direction}.{group}_pyNAST_{primer}.full.revcomped.fastq"
	conda:
		"envs/biopython.yaml"
	shell:
		"./scripts/revcompfastq_according_to_pyNAST.py --inpynast {input.fasta} --infastq {input.fastq} --outfastq {output.fastq_revcomp}"


rule grab_matching_cutadapt_full:
	input:
		"intermediate/compute-workflow/07-subsetted-fastq/{sample}.SSU.{direction}.{group}_pyNAST_{primer}.full.revcomped.fastq"
	output:
		mismatch="intermediate/compute-workflow/08-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.nohit.fastq",
		match="intermediate/compute-workflow/08-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.hit.fastq",
		info="intermediate/compute-workflow/08-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.info"
	log:
		"logs/08-cutadapt/{sample}.{direction}.{group}.{primer}.{mismatches}.cutadapt.log"
	conda:
		"envs/cutadapt-env.yaml"
	params:
		pattern=lambda wildcards : config["primer"][wildcards.primer],
		errorRate=lambda wildcards : config["mismatches"][wildcards.mismatches] / len(config["primer"][wildcards.primer]),
		lengthPrimer=lambda wildcards : len(config["primer"][wildcards.primer])
	shell:
		"cutadapt --info-file={output.info} --no-indels --no-trim --overlap={params.lengthPrimer} -b {params.pattern} --error-rate={params.errorRate} --untrimmed-output={output.mismatch} --output={output.match} {input}"


rule quality_filter_primer_region:
    #Keep only sequences that have >30 phred score across the whole primer + 5 leading/trailing bases (implicit in script)
	input:
		mismatch="intermediate/compute-workflow/08-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.nohit.fastq",
		match="intermediate/compute-workflow/08-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.hit.fastq",
		info="intermediate/compute-workflow/08-checked/{primer}/6-mismatch/{sample}.SSU.{direction}.{group}.{primer}.6-mismatch.info" #Assume anything with -1 value is false positive, ignoring those with > ~30% mismatches to primer (e.g. for a 20bp primer)
	output:
		mismatch="intermediate/compute-workflow/08-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.nohit.filtered.fastq",
		match="intermediate/compute-workflow/08-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.hit.filtered.fastq"
	conda:
		"envs/biopython.yaml"
	shell:
		"scripts/filter-fastq-by-primer-ROI-coverage.py --info {input.info} --fastq {input.match} > {output.match} ; "
		"scripts/filter-fastq-by-primer-ROI-coverage.py --info {input.info} --fastq {input.mismatch} > {output.mismatch} "


rule compute_percentages:
	input:
		mismatch="intermediate/compute-workflow/08-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.nohit.filtered.fastq",
		match="intermediate/compute-workflow/08-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.hit.filtered.fastq"
	output:
		"output/compute-workflow/09-summary/{sample}.{direction}.{group}.{primer}.{mismatches}.summary.tsv"
	shell:
		"numMatch=`grep -cE \"^@\" {input.match} || numMatch=0` ; numMismatch=`grep -cE \"^@\" {input.mismatch}` || numMismatch=0 ; "
		"sumTotal=`expr $numMatch + $numMismatch || sumTotal=0` ; "
		"if `[ $sumTotal -eq 0 ]` ; then fracMatch=\"NA\"; elif `[ $numMatch -ne 0 ]` && `[ $numMismatch -eq 0 ]`; then fracMatch=1; elif `[ $numMatch -eq 0 ]` && `[ $numMismatch -ne 0 ]`; then fracMatch=0; elif `[ $numMismatch -ne 0 ]` && `[ $numMatch -ne 0 ]`; then fracMatch=`bc <<< \"scale=4; $numMatch/$sumTotal\"`; fi ; "
		"printf \"{wildcards.sample}\t{wildcards.direction}\t{wildcards.group}\t{wildcards.primer}\t{wildcards.mismatches}\t$sumTotal\t$numMatch\t$numMismatch\t$fracMatch\n\" > {output}"
