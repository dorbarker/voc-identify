from pathlib import Path

samples = [p.stem for p in Path('bams').glob('*.bam')]

rule all:
	input:
		expand("reports/{sample}/summary.txt", sample=samples)

rule bam_index:
		input:
			"bams/{sample}.bam"
		output:
			"bams/{sample}.bam.bai"
		shell:
			"samtools index -b {input}"

rule find_vocs:
	input:
		"bams/{sample}.bam",
		"bams/{sample}.bam.bai"
	output:
		"reports/{sample}/summary.txt"
	shell:
		"mmmvi "
		"--bam {input[0]} "
		"--reference data/reference.fasta "
		"--mutations data/mutations.tsv "
		"--outdir reports/{wildcards.sample}/ "
