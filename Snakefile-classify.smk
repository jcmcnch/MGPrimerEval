rule all:
	input:
		#expand("15-matches-classified/individual/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.hit.filtered.VSEARCHsintax-SILVA132.tax", sample=config["samples"], group=config["groups"], primer=config["primer"], mismatches=config["mismatches"], direction=['fwd','rev'])
		#expand("13-classified/individual/{sample}.SSU.{direction}.{group}.{primer}.{mismatches}.nohit.filtered.VSEARCHsintax-SILVA132.tax", sample=config["samples"], group=config["groups"], primer=config["primer"], mismatches=config["mismatches"], direction=['fwd','rev'])
		expand("intermediate/{study}.{group}.{primer}.0-mismatch.hits.all.order.counts.normalized.tsv", study=config["study"], sample=config["samples"], group=config["groups"], primer=config["primer"])
		#expand("intermediate/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.minAbund-0.01.tsv", study=config["study"], sample=config["samples"], group=config["groups"], primer=config["primer"], mismatches=["0-mismatch", "1-mismatch", "2-mismatch"])

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
		"find ./12-full-fastas -type f -name \"*2-mismatch*.fasta\" -print0 | xargs -0 cat > {output}"

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

"""
Rules below comprise a workflow for generating summaries of
which taxa are most discriminated against by a particular primer set.
Implemented using common bash tools and tested on Ubuntu 16.04, YMMV.
"""

rule cat_tax_for_all_samples_matches:
	output:
		"intermediate/{study}.{group}.{primer}.0-mismatch.hits.all.tax"
	shell:
		"find 15-matches-classified/individual/ -type f -name "
		"\"*{wildcards.group}*{wildcards.primer}*0-mismatch*tax\" -print0 | "
		"xargs -0 cat > {output}"

rule cat_tax_for_all_samples_mismatches:
	output:
		"intermediate/{study}.{group}.{primer}.{mismatches}.nohits.all.tax"
	shell:
		"find 13-classified/individual/ -type f -name "
		"\"*{wildcards.group}*{wildcards.primer}*{wildcards.mismatches}*tax\" -print0 | "
		"xargs -0 cat > {output}"

#Counting order-level groupings (can adjust level with the "cut -d, -f1-4" parameter below)
rule count_tax_matches:
	input:
		"intermediate/{study}.{group}.{primer}.0-mismatch.hits.all.tax"
	output:
		"intermediate/{study}.{group}.{primer}.0-mismatch.hits.all.order.counts.tsv"
	shell:
		"sed -re 's/\([0-9]{{1}}\.[0-9]{{2}}\)//g' {input} |" #Remove confidence estimations from VSEARCH output
		"cut -f2 | sort | cut -d, -f1-4 | sort | uniq -c | " #Take only tax column, collapse to order level, then count unique occurrences
		"tail -f -n +2 | awk '{{print $1,\"\t\",$2}}' > {output}" #Process output into tsv format to stdout

rule count_tax_mismatches:
	input:
		"intermediate/{study}.{group}.{primer}.{mismatches}.nohits.all.tax"
	output:
		"intermediate/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.tsv"
	shell:
		"sed -re 's/\([0-9]{{1}}\.[0-9]{{2}}\)//g' {input} |" #Remove confidence estimations from VSEARCH output
		"cut -f2 | sort | cut -d, -f1-4 | sort | uniq -c | " #Take only tax column, collapse to order level, then count unique occurrences
		"tail -f -n +2 | awk '{{print $1,\"\t\",$2}}' > {output}" #Process output into tsv format to stdout

#Now take only those with greater than 1 % abundance using basic python script (can change abundance cutoff if you desire)
rule filter_tax_matches_by_abundance:
	input:
		"intermediate/{study}.{group}.{primer}.0-mismatch.hits.all.order.counts.tsv"
	output:
		"intermediate/{study}.{group}.{primer}.0-mismatch.hits.all.order.counts.min0.01.tsv" #Should probably rename if you change fractional value
	shell:
		"scripts/mismatch-characterization/filter-by-fractional-abundance.py {input} 0.01 > {output} " #Can change fractional value here

rule filter_tax_mismatches_by_abundance:
	input:
		"intermediate/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.tsv"
	output:
		"intermediate/{study}.{group}.{primer}.{mismatches}.nohits.all.order.counts.min0.01.tsv" #Should probably rename if you change fractional value
	shell:
		"scripts/mismatch-characterization/filter-by-fractional-abundance.py {input} 0.01 > {output} " #Can change fractional value here

#Since the matches were subsampled, they need to be normalized before calculating fractions to make them equivalent to the mismatches which were not subsampled (took all of them)
rule normalize_match_counts_by_total_seqs:
	input:
		"intermediate/{study}.{group}.{primer}.0-mismatch.hits.all.order.counts.tsv"
	output:
		"intermediate/{study}.{group}.{primer}.0-mismatch.hits.all.order.counts.normalized.tsv"
	shell:
		"totalSeqs=`find 07-subsetted/ -type f -name \"*{wildcards.group}*{wildcards.primer}*\" -print0 | "
		"xargs -0 cat | grep -c \">\"` ; " #Count number of sequence records in file corresponding to group and primer
		"subsampledSeqs=`cat {input} | wc -l` ; " #Count number of subsampled seqs
		"fracSubsampled=`bc <<< \"scale=4; $subsampledSeqs/$totalSeqs\" ; " #Calculate fraction subsampled
		"while read line ; do ; "
		"num=`echo $line | awk '{{print $1}}'` ; "
		"tax=`echo $line | awk '{{print $2}}'` ; "
		"normalizedNum=`bc <<< \"scale=4; $num/$fracSubsampled\"` ; "
		"printf \"$tax\t$normalizedNum\n\" ; done < {input} > {output}"

#need to generate summary file for following script
rule count_total_filtered_hits:


#fixing bash script so can run with snakemake
rule compute_frac_mismatched:
	input:
		"intermediate/{study}.{group}.{primer}.0-mismatch.hits.all.order.counts.min0.01.tsv"
	output:
		targets="intermediate/{study}.{group}.{primer}.targets"
	shell:
		"scripts/mismatch-characterization/compute-frac-mismatched.sh {wildcards.primer}.{wildcards.group} {input} {output.targets}"
