from pathlib import Path

samples = [p.stem for p in Path('bams').glob('*.bam')]

rule all:
	input:
		expand("reports/{sample}/summary.txt", sample=samples)
rule find_vocs:
	input:
		"bams/{sample}.bam"
	output:
		"reports/{sample}/summary.txt"
	shell:
		"python mmmvi.py "
		"--bam {input} "
		"--reference data/nCoV-2018.fasta "
		"--mutations data/mutations.tsv "
		"--outdir reports/{wildcards.sample}/ "
