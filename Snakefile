rule index_reference:
	input:
		'reference/reference.fasta'
	output:
		'reference/reference.mmi'

	shell:
		'minimap2 -d {output} {input}'

rule convert_bam_to_fasta:
	input:
		'samples/{name}.bam'
	output:
		'inputs/{name}.fasta'

	shell:
		'samtools fasta {input} > {output}'

rule link_fastq:
	input:
		'samples/{name}.fastq'
	output:
		'inputs/{name}.fastq'
	group:
		'fastq_symlink'
	shell:
		'ln -sr {input} {output}'

rule link_fastqgz:
	input:
		'samples/{name}.fastq.gz'
	output:
		'inputs/{name}.fastq.gz'
	group:
		'fastqgz_symlink'
	shell:
		'ln -sr {input} {output}'

rule align:
	input:
		rules.index_reference.output,
		'inputs/{name}.{ext}'
	output:
		'outputs/{name}.sam'

	threads:
		13

	shell:
		'minimap2 -t 12 -a {input[0]} {input[1]} > {output}'

