#configfile: "config-27F.yaml"

rule all:
	input:
		#expand("10-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.info", sample=config["samples"], group=config["groups"], primer=config["primer"], mismatches=config["mismatches"], direction=['fwd','rev'])
		expand("11-summary/{sample}.{direction}.{group}.{primer}.{mismatches}.summary.tsv", sample=config["samples"], group=config["groups"], primer=config["primer"], mismatches=config["mismatches"], direction=['fwd','rev'])

rule fastp_clean:
	input:
		R1="00-fastq/{sample}_repaired_1.fastq.gz",
		R2="00-fastq/{sample}_repaired_2.fastq.gz"
	output:
		R1clean="01-fastp-cleaned/{sample}_1_clean.fastq.gz",
		R2clean="01-fastp-cleaned/{sample}_2_clean.fastq.gz",
		log="logs/01-fastp_cleaning/{sample}.log.html"
	threads: 8
	shell:
		"fastp -3 -t 3 -i {input.R1} -I {input.R2} -o {output.R1clean} -O {output.R2clean} --thread {threads} --html {output.log}"


rule graftm_sift:
	input:
		R1clean="01-fastp-cleaned/{sample}_1_clean.fastq.gz",
		R2clean="01-fastp-cleaned/{sample}_2_clean.fastq.gz",
		package="graftm/7.40.2013_08_greengenes_97_otus.with_euks.gpkg"
	output:
		graftm=directory("02-graftm_sifted/{sample}/"),
		log="logs/02-graftM_sifting/{sample}.graftM_log.txt",
		R1hits="02-graftm_sifted/{sample}.fwd.SSU.hits.fa",
		R2hits="02-graftm_sifted/{sample}.rev.SSU.hits.fa",
		R1hmmout="02-graftm_sifted/{sample}.fwd.SSU.hmmout.txt",
		R2hmmout="02-graftm_sifted/{sample}.rev.SSU.hmmout.txt"
	threads:
		1
	conda:
		"envs/graftm.yaml"
	shell:
		"rm -f {output.log} ; export PATH=\"$PWD/executables:$PATH\" ; "
		"graftM graft --force --forward {input.R1clean} --reverse {input.R2clean} --search_and_align_only --threads {threads} "
		"--graftm_package {input.package} --input_sequence_type nucleotide --search_method hmmsearch "
		"--verbosity 5 --log {output.log} --output_directory {output.graftm} ; "
		"ln -s `ls $PWD/02-graftm_sifted/{wildcards.sample}/{wildcards.sample}*/forward/{wildcards.sample}*_hits.fa` {output.R1hits} ; "
		"ln -s `ls $PWD/02-graftm_sifted/{wildcards.sample}/{wildcards.sample}*/reverse/{wildcards.sample}*_hits.fa` {output.R2hits} ; "
		"ln -s `ls $PWD/02-graftm_sifted/{wildcards.sample}/{wildcards.sample}*/forward/{wildcards.sample}*hmmout.txt` {output.R1hmmout} ; "
		"ln -s `ls $PWD/02-graftm_sifted/{wildcards.sample}/{wildcards.sample}*/reverse/{wildcards.sample}*hmmout.txt` {output.R2hmmout}"

rule remove_repeats_komplexity:
	input:
		R1sifted="02-graftm_sifted/{sample}.fwd.SSU.hits.fa",
		R2sifted="02-graftm_sifted/{sample}.rev.SSU.hits.fa"
	output:
		R1="03-low-complexity-filtered/{sample}.fwd.SSU.keep.fa",
		R2="03-low-complexity-filtered/{sample}.rev.SSU.keep.fa"
	conda:
		"envs/komplexity.yaml"
	shell:
		"kz --filter --fasta < {input.R1sifted} > {output.R1} ; "
		"kz --filter --fasta < {input.R2sifted} > {output.R2} "

rule sort_EUK:
	input:
		"03-low-complexity-filtered/{sample}.{direction}.SSU.keep.fa"
	output:
		PROK="04-sorted/{sample}.{direction}.SSU.PROK.fa",
		EUK="04-sorted/{sample}.{direction}.SSU.EUK.fa",
	log:
		"logs/04-sorting/{sample}.{direction}.bbsplit.PROK-EUK.log"
	threads: 8
	conda:
		"envs/bbmap.yaml"
	shell:
		"bbsplit.sh threads={threads} -Xmx100g overwrite=t usequality=f qtrim=f minratio=0.30 minid=0.30 pairedonly=f "
		"path=/home/db/bbsplit-db/EUK-PROK-bbsplit-db/ "
		"in={input} basename={input}_%.fasta ; "
		"mv {input}*EUK*.fasta {output.EUK}; mv {input}*PROK*.fasta {output.PROK}"


rule sort_PROK:
	input:
		"04-sorted/{sample}.{direction}.SSU.PROK.fa"
	output:
		ARCH="04-sorted/{sample}.{direction}.SSU.ARCH.fa",
		BACT="04-sorted/{sample}.{direction}.SSU.BACT.fa",
	log:
		"logs/04-sorting/{sample}.{direction}.bbsplit.ARCH-BACT.log"
	threads: 8
	conda:
		"envs/bbmap.yaml"
	shell:
		"bbsplit.sh threads={threads} -Xmx100g overwrite=t usequality=f qtrim=f minratio=0.30 minid=0.30 pairedonly=f "
		"path=/home/db/bbsplit-db/BACT-ARCH-bbsplit-db/ "
		"in={input} basename={input}_%.fasta ; "
		"mv {input}*ARCH*.fasta {output.ARCH}; mv {input}*BACT*.fasta {output.BACT}"


rule sort_CYANO:
	input:
		"04-sorted/{sample}.{direction}.SSU.BACT.fa"
	output:
		CYANO="04-sorted/{sample}.{direction}.SSU.BACT-CYANO.fa",
		NONCYANO="04-sorted/{sample}.{direction}.SSU.BACT-NON-CYANO.fa",
	log:
		"logs/04-sorting/{sample}.{direction}.bbsplit.BACT-CYANO.log"
	threads:
		8
	conda:
		"envs/bbmap.yaml"
	shell:
		"bbsplit.sh threads={threads} -Xmx100g overwrite=t usequality=f qtrim=f minratio=0.30 minid=0.30 pairedonly=f "
		"path=/home/db/bbsplit-db/BACT-CYANO-bbsplit-db/ "
		"in={input} basename={input}_%.fasta ; "
		"mv {input}*NON-CYANO*.fasta {output.NONCYANO}; mv {input}*CYANO*.fasta {output.CYANO}"


rule align_ARCH:
	input:
		seqs="04-sorted/{sample}.{direction}.SSU.ARCH.fa",
		ref="SSU_refs/Sulfolobus_acidocaldarius_N8_16s.fasta"
	output:
		alignment="05-pyNAST-aligned/{sample}.{direction}.ARCH/{sample}.{direction}.SSU.ARCH_pynast_aligned.fasta",
		log="logs/05-pyNAST-aligning/{sample}.{direction}.SSU.ARCH.pyNAST.log"
	conda:
		"envs/pynast.yaml"
	shell:
		"pynast -p 10 -l 1 -i {input.seqs} -t {input.ref} ;"
		"echo \"pyNAST info: `cat 04-sorted/{wildcards.sample}.SSU.ARCH_pynast_fail.fasta | wc -l` sequences failed to align.\" ;"
		"mv 04-sorted/{wildcards.sample}.{wildcards.direction}.SSU.ARCH_pynast_log.txt {output.log} ; mv 04-sorted/{wildcards.sample}.{wildcards.direction}.SSU.ARCH_pynast_fail.fasta {output.log} ; "
		"mv 04-sorted/{wildcards.sample}.{wildcards.direction}.SSU.ARCH_pynast_aligned.fasta {output.alignment} "


rule align_BACT:
	input:
		seqs="04-sorted/{sample}.{direction}.SSU.BACT-NON-CYANO.fa",
		ref="SSU_refs/Ecoli_16s.fna"
	output:
		alignment="05-pyNAST-aligned/{sample}.{direction}.BACT-NON-CYANO/{sample}.{direction}.SSU.BACT-NON-CYANO_pynast_aligned.fasta",
		log="logs/05-pyNAST-aligning/{sample}.{direction}.SSU.BACT-NON-CYANO.pyNAST.log"
	conda:
		"envs/pynast.yaml"
	shell:
		"pynast -p 10 -l 1 -i {input.seqs} -t {input.ref} ;"
		"echo \"pyNAST info: `cat 04-sorted/{wildcards.sample}.SSU.BACT-NON-CYANO_pynast_fail.fasta | wc -l` sequences failed to align.\" ;"
		"mv 04-sorted/{wildcards.sample}.{wildcards.direction}.SSU.BACT-NON-CYANO_pynast_log.txt {output.log} ; mv 04-sorted/{wildcards.sample}.{wildcards.direction}.SSU.BACT-NON-CYANO_pynast_fail.fasta {output.log} ; "
		"mv 04-sorted/{wildcards.sample}.{wildcards.direction}.SSU.BACT-NON-CYANO_pynast_aligned.fasta {output.alignment} "


rule align_CYANO:
	input:
		seqs="04-sorted/{sample}.{direction}.SSU.BACT-CYANO.fa",
		ref="SSU_refs/longest-CYANO-with-27F.fasta"
	output:
		alignment="05-pyNAST-aligned/{sample}.{direction}.BACT-CYANO/{sample}.{direction}.SSU.BACT-CYANO_pynast_aligned.fasta",
		log="logs/05-pyNAST-aligning/{sample}.{direction}.SSU.BACT-CYANO.pyNAST.log"
	conda:
		"envs/pynast.yaml"
	shell:
		"pynast -p 10 -l 1 -i {input.seqs} -t {input.ref} ;"
		"echo \"pyNAST info: `cat 04-sorted/{wildcards.sample}.SSU.BACT-CYANO_pynast_fail.fasta | wc -l` sequences failed to align.\" ;"
		"mv 04-sorted/{wildcards.sample}.{wildcards.direction}.SSU.BACT-CYANO_pynast_log.txt {output.log} ; mv 04-sorted/{wildcards.sample}.{wildcards.direction}.SSU.BACT-CYANO_pynast_fail.fasta {output.log} ; "
		"mv 04-sorted/{wildcards.sample}.{wildcards.direction}.SSU.BACT-CYANO_pynast_aligned.fasta {output.alignment} "

rule align_EUK:
	input:
		seqs="04-sorted/{sample}.{direction}.SSU.EUK.fa",
		ref="SSU_refs/Saccharomyces_cerevisiae_S288C_18s-1_NR_132213.1.fa"
	output:
		alignment="05-pyNAST-aligned/{sample}.{direction}.EUK/{sample}.{direction}.SSU.EUK_pynast_aligned.fasta",
		log="logs/05-pyNAST-aligning/{sample}.{direction}.SSU.EUK.pyNAST.log"
	conda:
		"envs/pynast.yaml"
	shell:
		"pynast -p 10 -l 1 -i {input.seqs} -t {input.ref} ;"
		"echo \"pyNAST info: `cat 04-sorted/{wildcards.sample}.SSU.EUK.pynast_fail.fasta | wc -l` sequences failed to align.\" ;"
		"mv 04-sorted/{wildcards.sample}.{wildcards.direction}.SSU.EUK_pynast_log.txt {output.log} ; mv 04-sorted/{wildcards.sample}.{wildcards.direction}.SSU.EUK_pynast_fail.fasta {output.log} ; "
		"mv 04-sorted/{wildcards.sample}.{wildcards.direction}.SSU.EUK_pynast_aligned.fasta {output.alignment} "


rule get_strandedness:
	input:
		fwd="02-graftm_sifted/{sample}.fwd.SSU.hmmout.txt",
		rev="02-graftm_sifted/{sample}.rev.SSU.hmmout.txt"
	output:
		fwd="strand_info/{sample}.strand.fwd.tsv",
		rev="strand_info/{sample}.strand.rev.tsv"
	shell:
		"sed '/^#/d' {input.fwd} | awk '{{print $1,\"\t\",$12}}' > {output.fwd} ; "
		"sed '/^#/d' {input.rev} | awk '{{print $1,\"\t\",$12}}' > {output.rev}"


rule get_ids:
	input:
		"05-pyNAST-aligned/{sample}.{direction}.{group}/{sample}.{direction}.SSU.{group}_pynast_aligned.fasta"
	output:
		"aligned-seq-ids/{sample}.SSU.{group}.pyNAST.{direction}.ids",
	shell:
		"seqmagick extract-ids {input} > {output} "


#obtain fastq records from headers for clean, filtered SSU reads
rule get_fastq:
	input:
		fwdIDs="aligned-seq-ids/{sample}.SSU.{group}.pyNAST.fwd.ids",
		revIDs="aligned-seq-ids/{sample}.SSU.{group}.pyNAST.rev.ids",
		fastqR1="01-fastp-cleaned/{sample}_1_clean.fastq.gz",
		fastqR2="01-fastp-cleaned/{sample}_2_clean.fastq.gz"
	output:
		fwd="06-fastq/{sample}.SSU.{group}.keep.fwd.fastq",
		rev="06-fastq/{sample}.SSU.{group}.keep.rev.fastq"
	conda:
		"envs/bbmap.yaml"
	shell:
		"filterbyname.sh names={input.fwdIDs} overwrite=t include=t in={input.fastqR1} out={output.fwd} ; "
		"filterbyname.sh names={input.revIDs} overwrite=t include=t in={input.fastqR2} out={output.rev} "


rule subset_to_primer_region:
	input:
		"05-pyNAST-aligned/{sample}.{direction}.{group}/{sample}.{direction}.SSU.{group}_pynast_aligned.fasta"
	output:
		"07-subsetted/{sample}.SSU.{direction}.{group}_pyNAST_{primer}.fasta"
	conda:
		"envs/biopython.yaml"
	params:
		start=lambda wildcards: config["primerROI"][wildcards.group][wildcards.primer][0],
		end=lambda wildcards: config["primerROI"][wildcards.group][wildcards.primer][1]
	shell:
		"scripts/filter-pyNAST-for-ROI.py --input {input} --output {output} --start {params.start} --end {params.end}"


rule get_fastq_for_subset:
	#Get the whole fastq files associated with the matching reads
	input:
		fasta="07-subsetted/{sample}.SSU.{direction}.{group}_pyNAST_{primer}.fasta",
		fastq="06-fastq/{sample}.SSU.{group}.keep.{direction}.fastq"
	output:
		"08-fastq-primer-region/{sample}.SSU.{direction}.{group}_pyNAST_{primer}.full.fastq"
	conda:
		"envs/bbmap.yaml"
	shell:
		"filterbyname.sh names={input.fasta} include=t in={input.fastq} out={output}"


rule reverse_complement_fastqs:
	#Cutadapt doesn't search the reverse complement, so have to do this semi-manually using the hmm output
	input:
		fastq="08-fastq-primer-region/{sample}.SSU.{direction}.{group}_pyNAST_{primer}.full.fastq",
		strandedness="strand_info/{sample}.strand.{direction}.tsv"
	output:
		"09-complemented-fastqs/{sample}.SSU.{direction}.{group}_pyNAST_{primer}.full.revcomped.fastq"
	shell:
		"scripts/reverse-complement-fastq-according-to-HMM-output.py --fastqinput {input.fastq} --strand {input.strandedness} > {output}"


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


rule quality_filter_primer_region_27F:
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
