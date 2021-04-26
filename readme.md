# Identifying Variants of Concern in Metagenomic Samples

`mmmvi` is an in-progress tool to identify potential SARS-CoV-2 Variants of Concern
(VoC) in metagenomic samples. A manuscript is in preparation, but in the interim,
please cite this repository.

A VoC is identified by the presence of nucleotide differences relative to a
reference strain. In a heterogeneous, metagenomic sample, sequencing reads may
originate from multiple SARS-CoV-2 genomes. To improve confidence that a VoC is
present, we can look for cases where multiple nucleotide polymorphisms are
present on the same read.

# Install
```sh
conda install -c conda-forge -c bioconda -c dorbarker mmmvi
```

# Running mmmvi

## Basic Usage

```sh
mmmvi --bam your_sample.bam --mutations mutations.tsv --reference nCoV-2019.fasta --outdir reports/your_sample
```

## Command Line Arguments

```sh
usage: mmmvi.py [-h] --bam BAM --reference FASTA --mutations TABULAR --outdir
                DIR [--voc-column COLUMN] [--mutation-column COLUMN] [-d CHAR]
                [--only-vocs [VOC [VOC ...]]] [-v]

optional arguments:
  -h, --help            show this help message and exit
  --bam BAM             Path to a BAM file aligned against the reference
  --reference FASTA     Path to FASTA-formatted complete reference genome
  --mutations TABULAR   Path to tabular file describing Variants of Concern
  --outdir DIR          Output directory; will be created if it does not
                        already exist
  --voc-column COLUMN   Header for the column containing Variant of Concern
                        names [PangoLineage]
  --mutation-column COLUMN
                        Header for the column containing mutation descriptions
                        [NucName]
  -d CHAR, --delimiter CHAR
                        Delimiter character for tabular input and output [TAB]
  --only-vocs [VOC [VOC ...]]
                        Look only for these variants, plus the reference
  -v, --version         show program's version number and exit
```

## The Mutations File

The mutations file is a delimited file (tabs, by default) which describe which
mutations belong to which variants of concern. The default headers are
"PangoLineage" for VOC and "NucName" for the mutation description. However,
these may be overridden by the user with `--voc-column` and `--mutation-column`,
respectively. Additional columns are ignored.

All positions are relative to the reference genome and 1-based _i.e.,_ the first
position in the genome is position 1.

```sh
# mutations.tsv

voc        mutation
B.1.1.7    G24914C          # 1 base substitution from G to C at 24914
B.1.1.7    GAT28280CTA      # 3 base substitution from GAT to CTA at 282280-282282
B.1.1.7    [21765-21770]del # 21765-21770, inclusive, are deleted
P.1        28262AACA        # 4 base insertion between 28262 and 28263
```

# Automated Workflow
A [snakemake](https://snakemake.readthedocs.io/en/stable/) workflow,
`voc-identify.smk` is provided for convenience and easier parallel processing.

It assumes your project directory is structured thusly:

```
.
├── bams
│   ├── sample1.bam 
│   └── sample2.bam
│   |   ...
│   └── sampleN.bam
├── data
    ├── mutations.tsv
    └── reference.fasta
```

The workflow will place output reports for each sample in a corresponding
directory under `./reports/`:

```
reports/
└── sample1
    ├── cooccurence_matrices
    │   ├── absolute
    │   │   └── B.1.1.7.txt
    │   └── relative
    │       └── B.1.1.7.txt
    ├── read_report.txt
    ├── read_species.txt
    └── summary.txt
```
