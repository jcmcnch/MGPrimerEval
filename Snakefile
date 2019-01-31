configfile: "config.yaml"

rule all:
	input:
		expand("logs/05-pyNAST-aligning/{sample}.{direction}.SSU.BACT-NON-CYANO.pyNAST.log", sample=config["samples"], direction=['fwd','rev']),
		expand("logs/05-pyNAST-aligning/{sample}.{direction}.SSU.BACT-CYANO.pyNAST.log", sample=config["samples"], direction=['fwd','rev']),
		expand("logs/05-pyNAST-aligning/{sample}.{direction}.SSU.ARCH.pyNAST.log", sample=config["samples"], direction=['fwd','rev']),
		expand("logs/05-pyNAST-aligning/{sample}.{direction}.SSU.EUK.pyNAST.log", sample=config["samples"], direction=['fwd','rev'])

rule fastp_clean:
	input:
		R1="00-fastq/{sample}_1.fastq.gz",
		R2="00-fastq/{sample}_2.fastq.gz"
	output:
		R1="01-fastp-cleaned/{sample}_1_clean.fastq.gz",
		R2="01-fastp-cleaned/{sample}_2_clean.fastq.gz",
		log="logs/01-fastp_cleaning/{sample}.log.html"
	threads: 8
	shell:
		"fastp -3 -t 3 -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} --thread {threads} --html {output.log}"


rule graftm_sift:
	input:
		R1="01-fastp-cleaned/{sample}_1_clean.fastq.gz",
		R2="01-fastp-cleaned/{sample}_2_clean.fastq.gz",
		package="graftm/7.40.2013_08_greengenes_97_otus.with_euks.gpkg"
	output:
		results=directory("02-graftm_sifted/{sample}"),
		R1hits="02-graftm_hits/{sample}_1_clean_forward_hits.fa",
		R2hits="02-graftm_hits/{sample}_1_clean_reverse_hits.fa",
		log="logs/02-graftM_sifting/{sample}.graftM_log.txt"
	threads: 8
	conda:
		"envs/graftm.yaml"
	shell:
		"rm -f {output.log} ; export PATH=\"$PWD/executables:$PATH\" ; "
		"graftM graft --force --forward {input.R1} --reverse {input.R2} --search_and_align_only --threads {threads} "
		"--graftm_package {input.package} --input_sequence_type nucleotide --search_method hmmsearch "
		"--verbosity 5 --log {output.log} --output_directory {output.results} ; cp {output.results}/*/*/*forward_hits.fa {output.R1hits} ; cp {output.results}/*/*/*reverse_hits.fa {output.R2hits}"

rule remove_repeats_komplexity:
	input:
		R1="02-graftm_hits/{sample}_1_clean_forward_hits.fa",
		R2="02-graftm_hits/{sample}_1_clean_reverse_hits.fa"
	output:
		R1="03-low-complexity-filtered/{sample}.fwd.SSU.keep.fa",
		R2="03-low-complexity-filtered/{sample}.rev.SSU.keep.fa"
	conda:
		"envs/komplexity.yaml"
	shell:
		"kz --filter --fasta < {input.R1} > {output.R1} ; "
		"kz --filter --fasta < {input.R2} > {output.R2} "


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
	threads: 8
	conda:
		"envs/bbmap.yaml"
	shell:
		"bbsplit.sh threads={threads} -Xmx100g overwrite=t usequality=f qtrim=f minratio=0.30 minid=0.30 pairedonly=f "
		"path=/home/db/bbsplit-db/BACT-CYANO-bbsplit-db/ "
		"in={input} basename={input}_%.fasta ; "
		"mv {input}*NON-CYANO*.fasta {output.NONCYANO}; mv {input}*CYANO*.fasta {output.CYANO}"


rule align_EUK:
	input:
		seqs="04-sorted/{sample}.{direction}.SSU.EUK.fa",
		ref="SSU_refs/Saccharomyces_cerevisiae_S288C_18s-1_NR_132213.1.fa"
	output:
		dir=directory("05-pyNAST-aligned/{sample}.{direction}.EUK/"),
		log="logs/05-pyNAST-aligning/{sample}.{direction}.SSU.EUK.pyNAST.log"
	conda:
		"envs/pynast.yaml"
	shell:
		"pynast -p 10 -l 1 -i {input.seqs} -t {input.ref} ; "
		"echo \"pyNAST info: `cat 04-sorted/{wildcards.sample}*pynast_fail.fasta | wc -l` sequences failed to align.\" ; "
		"mv 04-sorted/{wildcards.sample}.{wildcards.direction}*EUK_pynast_log.txt {output.log} ; mv 04-sorted/{wildcards.sample}.{wildcards.direction}*EUK*pynast* {output.dir}"


rule align_ARCH:
	input:
		seqs="04-sorted/{sample}.{direction}.SSU.ARCH.fa",
		ref="SSU_refs/Sulfolobus_acidocaldarius_N8_16s.fasta"
	output:
		dir=directory("05-pyNAST-aligned/{sample}.{direction}.ARCH/"),
		log="logs/05-pyNAST-aligning/{sample}.{direction}.SSU.ARCH.pyNAST.log"
	conda:
		"envs/pynast.yaml"
	shell:
		"pynast -p 10 -l 1 -i {input.seqs} -t {input.ref} ;"
		"echo \"pyNAST info: `cat 04-sorted/{wildcards.sample}*pynast_fail.fasta | wc -l` sequences failed to align.\" ;"
		"mv 04-sorted/{wildcards.sample}.{wildcards.direction}*ARCH_pynast_log.txt {output.log} ; mv 04-sorted/{wildcards.sample}.{wildcards.direction}*ARCH*pynast* {output.dir} "


rule align_BACT:
	input:
		seqs="04-sorted/{sample}.{direction}.SSU.BACT-NON-CYANO.fa",
		ref="SSU_refs/Ecoli_16s.fna"
	output:
		dir=directory("05-pyNAST-aligned/{sample}.{direction}.NON-CYANO/"),
		log="logs/05-pyNAST-aligning/{sample}.{direction}.SSU.BACT-NON-CYANO.pyNAST.log"
	conda:
		"envs/pynast.yaml"
	shell:
		"pynast -p 10 -l 1 -i {input.seqs} -t {input.ref} ;"
		"echo \"pyNAST info: `cat 04-sorted/{wildcards.sample}*pynast_fail.fasta | wc -l` sequences failed to align.\" ;"
		"mv 04-sorted/{wildcards.sample}.{wildcards.direction}*BACT-NON-CYANO_pynast_log.txt {output.log} ; mv 04-sorted/{wildcards.sample}.{wildcards.direction}*BACT-NON-CYANO_*pynast* {output.dir} "


rule align_CYANO:
	input:
		seqs="04-sorted/{sample}.{direction}.SSU.BACT-CYANO.fa",
		ref="SSU_refs/SILVA_132_longest-CYANO-non-Chloroplast.fasta"
	output:
		dir=directory("05-pyNAST-aligned/{sample}.{direction}.CYANO/"),
		log="logs/05-pyNAST-aligning/{sample}.{direction}.SSU.BACT-CYANO.pyNAST.log"
	conda:
		"envs/pynast.yaml"
	shell:
		"pynast -p 10 -l 1 -i {input.seqs} -t {input.ref} ;"
		"echo \"pyNAST info: `cat 04-sorted/{wildcards.sample}*pynast_fail.fasta | wc -l` sequences failed to align.\" ;"
		"mv 04-sorted/{wildcards.sample}.{wildcards.direction}*BACT-CYANO_pynast_log.txt {output.log} ; mv 04-sorted/{wildcards.sample}.{wildcards.direction}*BACT-CYANO_*pynast* {output.dir} "
