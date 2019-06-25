rule all:
	input:
		expand("15-matches-classified/individual/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.hit.filtered.VSEARCHsintax-SILVA132.tax", sample=config["samples"], group=config["groups"], primer=config["primer"], mismatches=config["mismatches"], direction=['fwd','rev']),
		expand("13-classified/individual/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.nohit.filtered.VSEARCHsintax-SILVA132.tax", sample=config["samples"], group=config["groups"], primer=config["primer"], mismatches=config["mismatches"], direction=['fwd','rev'])

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


rule concatenate_mismatched_fastas:
#Note: concatenating to speed up classification step; using xargs for large numbers of files (otherwise bash will complain argument list too long)
	output:
		"tmp.mismatches.concatenated.fasta"
	shell:
		"find ./12-full-fastas -type f -name \"*.fasta\" -print0 | xargs -0 cat > {output}"

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
		"grep -f $tmpfile {input.concatenated} > {output} || touch {output} ; "


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


rule concatenate_matched_fastas:
	output:
		"tmp.matches.subsampled.concatenated.fasta"
	shell:
		"find ./14-subsampled-matched-fastas -type f -name \"*.fasta\" -print0 | xargs -0 cat > {output}"


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
