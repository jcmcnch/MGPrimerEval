configfile: "config.yaml"

rule all:
	input:
		expand("13-classified/individual/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.nohit.RDP-SILVA132.tax", sample=config["samples"], group=config["groups"], primer=config["primer"], mismatches=config["mismatches"], direction=['fwd','rev'])

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
		graftm="02-graftm_sifted/{sample}/",
		log="logs/02-graftM_sifting/{sample}.graftM_log.txt",
		R1hits="02-graftm_sifted/{sample}.fwd.SSU.hits.fa",
		R2hits="02-graftm_sifted/{sample}.rev.SSU.hits.fa"
	threads:
		8
	conda:
		"envs/graftm.yaml"
	shell:
		"rm -f {output.log} ; export PATH=\"$PWD/executables:$PATH\" ; "
		"graftM graft --force --forward {input.R1clean} --reverse {input.R2clean} --search_and_align_only --threads {threads} "
		"--graftm_package {input.package} --input_sequence_type nucleotide --search_method hmmsearch "
		"--verbosity 5 --log {output.log} --output_directory {output.graftm} ; "
		"ln -s $PWD/02-graftm_sifted/{sample}/{sample}*/forward/{sample}*_hits.fa {R1hits} ;"
		"ln -s $PWD/02-graftm_sifted/{sample}/{sample}*/reverse/{sample}*_hits.fa {R2hits}"

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
		"echo \"pyNAST info: `cat 04-sorted/{wildcards.sample}.ARCH.pynast_fail.fasta | wc -l` sequences failed to align.\" ;"
		"mv 04-sorted/{wildcards.sample}.{wildcards.direction}.ARCH_pynast_log.txt {output.log} ; mv 04-sorted/{wildcards.sample}.{wildcards.direction}.ARCH_pynast_fail.fasta {output.log} ; "
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
		"echo \"pyNAST info: `cat 04-sorted/{wildcards.sample}.BACT-NON-CYANO_pynast_fail.fasta | wc -l` sequences failed to align.\" ;"
		"mv 04-sorted/{wildcards.sample}.{wildcards.direction}.BACT-NON-CYANO_pynast_log.txt {output.log} ; mv 04-sorted/{wildcards.sample}.{wildcards.direction}.BACT-NON-CYANO_pynast_fail.fasta {output.log} ; "
		"mv 04-sorted/{wildcards.sample}.{wildcards.direction}.SSU.BACT-NON-CYANO_pynast_aligned.fasta {output.alignment} "


rule align_CYANO:
	input:
		seqs="04-sorted/{sample}.{direction}.SSU.BACT-CYANO.fa",
		ref="SSU_refs/SILVA_132_longest-CYANO-non-Chloroplast.fasta"
	output:
		alignment="05-pyNAST-aligned/{sample}.{direction}.BACT-CYANO/{sample}.{direction}.SSU.BACT-CYANO_pynast_aligned.fasta",
		log="logs/05-pyNAST-aligning/{sample}.{direction}.SSU.BACT-CYANO.pyNAST.log"
	conda:
		"envs/pynast.yaml"
	shell:
		"pynast -p 10 -l 1 -i {input.seqs} -t {input.ref} ;"
		"echo \"pyNAST info: `cat 04-sorted/{wildcards.sample}.BACT-CYANO_pynast_fail.fasta | wc -l` sequences failed to align.\" ;"
		"mv 04-sorted/{wildcards.sample}.{wildcards.direction}.BACT-CYANO_pynast_log.txt {output.log} ; mv 04-sorted/{wildcards.sample}.{wildcards.direction}.BACT-CYANO_pynast_fail.fasta {output.log} ; "
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
		"echo \"pyNAST info: `cat 04-sorted/{wildcards.sample}.EUK.pynast_fail.fasta | wc -l` sequences failed to align.\" ;"
		"mv 04-sorted/{wildcards.sample}.{wildcards.direction}.EUK_pynast_log.txt {output.log} ; mv 04-sorted/{wildcards.sample}.{wildcards.direction}.EUK_pynast_fail.fasta {output.log} ; "
		"mv 04-sorted/{wildcards.sample}.{wildcards.direction}.SSU.EUK_pynast_aligned.fasta {output.alignment} "


rule get_strandedness:
	input:
		fwd="02-graftm_sifted/{sample}/{sample}_repaired_1/forward/{sample}_repaired_1_forward.hmmout.txt",
		rev="02-graftm_sifted/{sample}/{sample}_repaired_1/reverse/{sample}_repaired_1_reverse.hmmout.txt"
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
		fastqR1="00-fastq/{sample}_repaired_1.fastq.gz",
		fastqR2="00-fastq/{sample}_repaired_2.fastq.gz"
	output:
		fwd="06-fastq/{sample}.SSU.{group}.keep.fwd.fastq",
		rev="06-fastq/{sample}.SSU.{group}.keep.rev.fastq"
	conda:
		"envs/bbmap.yaml"
	shell:
		"filterbyname.sh names={input.fwdIDs} app=t include=t in={input.fastqR1} out={output.fwd} ; "
		"filterbyname.sh names={input.revIDs} app=t include=t in={input.fastqR2} out={output.rev} "


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


rule get_fastq_region:
	#Get the region around the primer with 5 leading and trailing bases
	input:
		fasta="07-subsetted/{sample}.SSU.{direction}.{group}_pyNAST_{primer}.fasta",
		fastq="06-fastq/{sample}.SSU.{group}.keep.{direction}.fastq",
		strand="strand_info/{sample}.strand.{direction}.tsv"
	output:
		"08-fastq-primer-region/{sample}.SSU.{direction}.{group}_pyNAST_{primer}.slice.fastq"
	conda:
		"envs/biopython.yaml"
	params:
		lenROI=lambda wildcards: len(config["primer"][wildcards.primer]),
		start=lambda wildcards: config["primerROI"][wildcards.group][wildcards.primer][0],
		padding=5
	shell:
		"scripts/grab-primer-region-in-fastq.py --strand {input.strand} --lenROI {params.lenROI} --padding {params.padding} --start {params.start} --fasta {input.fasta} --fastq {input.fastq} > {output}"


rule quality_filter_primer_region:
    #Keep only sequences that have >30 phred score across the whole primer + 5 leading/trailing bases
	input:
		"08-fastq-primer-region/{sample}.SSU.{direction}.{group}_pyNAST_{primer}.slice.fastq"
	output:
		"09-qual-filtered-primer-region/{sample}.SSU.{direction}.{group}_pyNAST_{primer}.slice.gtQ30.fastq"
	params:
		minlen=lambda wildcards: len(config["primer"][wildcards.primer]) + 10 ,
		minqual=30
	conda:
		"envs/fastp.yaml"
	shell:
		"fastp -A -q {params.minqual} -l {params.minlen} -u 0 -i {input} -o {output}"


rule grab_matching_cutadapt:
	input:
		"09-qual-filtered-primer-region/{sample}.SSU.{direction}.{group}_pyNAST_{primer}.slice.gtQ30.fastq"
	output:
		mismatch="10-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.nohit.fastq",
		match="10-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.hit.fastq"
	log:
		"logs/05-checking/{sample}.{direction}.{group}.{primer}.{mismatches}.cutadapt.log"
	params:
		pattern=lambda wildcards : config["primer"][wildcards.primer],
		errorRate=lambda wildcards : config["mismatches"][wildcards.mismatches] / len(config["primer"][wildcards.primer])
	shell:
		"cutadapt --no-indels --no-trim -b {params.pattern} --error-rate {params.errorRate} --untrimmed-output={output.mismatch} --output={output.match} {input}"


rule compute_percentages:
	input:
		mismatch="10-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.nohit.fastq",
		match="10-checked/{primer}/{mismatches}/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.hit.fastq"
	output:
		"11-summary/{sample}.{direction}.{group}.{primer}.{mismatches}.summary.tsv"
	shell:
		"numMatch=`grep -cE \"^@\" {input.match} || numMatch=0` ; numMismatch=`grep -cE \"^@\" {input.mismatch}` || numMismatch=0 ; "
		"sumTotal=`expr $numMatch + $numMismatch || sumTotal=0` ; "
		"if `[ $sumTotal -eq 0 ]` ; then fracMatch=\"NA\"; elif `[ $numMatch -ne 0 ]` && `[ $numMismatch -eq 0 ]`; then fracMatch=1; elif `[ $numMatch -eq 0 ]` && `[ $numMismatch -ne 0 ]`; then fracMatch=0; elif `[ $numMismatch -ne 0 ]` && `[ $numMatch -ne 0 ]`; then fracMatch=`bc <<< \"scale=4; $numMatch/$sumTotal\"`; fi ; "
		"printf \"{wildcards.sample}\t{wildcards.direction}\t{wildcards.group}\t{wildcards.primer}\t{wildcards.mismatches}\t$sumTotal\t$numMatch\t$numMismatch\t$fracMatch\n\" > {output}"


rule grab_full_fastas:
	input:
		fwdFQ="10-checked/{primer}/{mismatches}/{sample}.SSU.fwd.{group}.{primer}.{mismatches}.nohit.fastq",
		revFQ="10-checked/{primer}/{mismatches}/{sample}.SSU.rev.{group}.{primer}.{mismatches}.nohit.fastq",
		fwdFA="02-graftm_sifted/{sample}/{sample}_repaired_1/forward/{sample}_repaired_1_forward_hits.fa",
		revFA="02-graftm_sifted/{sample}/{sample}_repaired_1/reverse/{sample}_repaired_1_reverse_hits.fa"
	output:
		fwd="12-full-fastas/{sample}.SSU.fwd.{group}.{primer}.{mismatches}.nohit.fasta",
		rev="12-full-fastas/{sample}.SSU.rev.{group}.{primer}.{mismatches}.nohit.fasta"
	conda:
		"envs/bbmap.yaml"
	shell:
		"filterbyname.sh names={input.fwdFQ} include=t in={input.fwdFA} out={output.fwd} ; "
		"filterbyname.sh names={input.revFQ} include=t in={input.revFA} out={output.rev} "

rule concatenate_fastas:
#Note: concatenating because otherwise RDP classifier gets retrained for each set of sequences (very slow!)
	input:
		expand("12-full-fastas/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.nohit.fasta", sample=config["samples"], group=config["groups"], primer=config["primer"], mismatches=config["mismatches"], direction=['fwd','rev'])
	output:
		"tmp.mismatches.concatenated.fasta"
	shell:
		"cat {input} >> {output}"


rule classify_mismatches:
	input:
		"tmp.mismatches.concatenated.fasta"
	output:
		directory("13-classified/all/all.mismatches.RDP.classified")
	threads:
		10
	conda:
		"envs/qiime1.yaml"
	shell:
		"export RDP_JAR_PATH=\"/home/anaconda/miniconda2/envs/qiime1/bin/rdp_classifier-2.2/rdp_classifier-2.2.jar\" ; "
		"assign_taxonomy.py --confidence 0 -m rdp -i {input} -o {output} --rdp_max_memory=500000 "
  		"-t /home/db/SILVA_132/qiime_db/SILVA_132_QIIME_release/taxonomy/taxonomy_all/90/taxonomy_7_levels.txt "
		"-r /home/db/SILVA_132/qiime_db/SILVA_132_QIIME_release/rep_set/rep_set_all/90/silva132_90.fna ;"
		"rm tmp.mismatches.concatenated.fasta"


rule deconcat_classifications:
	input:
		fasta="12-full-fastas/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.nohit.fasta",
		concatenated="13-classified/all/all.mismatches.RDP.classified"
	output:
		"13-classified/individual/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.nohit.RDP-SILVA132.tax"
	shell:
		"tmpfile=`mktemp /tmp/fastq-ids.XXXXXXXXXXXXXXXX` ; seqmagick extract-ids {input.fasta} > $tmpfile ; "
		"grep -f $tmpfile {input.concatenated}/tmp.mismatches.concatenated_tax_assignments.txt > {output} || touch {output} ; "
		"rm $tmpfile"
