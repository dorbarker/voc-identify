# Identifying Variants of Concern in Metagenomic Samples

`mmmvi` is a tool to identify potential SARS-CoV-2 Variants of Concern (VoC) in
metagenomic samples.

A VoC is identified by the presence of nucleotide differences relative to a
reference strain. In a heterogeneous, metagenomic sample, sequencing reads may
originate from multiple SARS-CoV-2 genomes. To improve confidence that a VoC is
present, we can look for cases where multiple nucleotide polymorphisms are
present on the same read.

# Install

## conda
```sh
conda install -c conda-forge -c bioconda -c dorbarker mmmvi
```

## pip
```sh
pip install git+https://github.com/dorbarker/voc-identify.git
```

# Running mmmvi

## Basic Usage

```sh
mmmvi --bam your_sample.bam        \
      --mutations mutations.tsv    \
      --reference reference.fasta  \
      --outdir reports/your_sample 
```

## Restrict your search to specific variants

```sh
mmmvi --bam your_sample.bam        \
      --mutations mutations.tsv    \
      --reference reference.fasta  \
      --outdir reports/your_sample \
      --only-vocs B.1.1.7 P.1
```

## Loading variant definitions from Public Health England YAML files

```sh
# Get the definitions for the first time
git clone https://github.com/phe-genomics/variant_definitions/

# Update existing definitions
git -C path/to/variant_definitions/ pull

mmmvi --mutations path/to/variant_definitions/variant_yaml/ \
      --bam your_sample.bam                                 \
      --reference reference.fasta                           \
      --outdir reports/your_sample                          \
      --only-vocs unloved-crouton denture-daughter
```

## Automation with Snakemake

See below for the mandatory directory structure required by this workflow.

A config file can be used to restrict which variants mmmvi will search for.
Omit `--configfile` and its argument if you wish to search for all variants.

Use the `--jobs` to control the number of parallel jobs.

```sh
snakemake --jobs 2 -s voc-identity.smk -d path/to/your/project --configfile example-config.yaml 
```

## Command Line Arguments

```sh
usage: mmmvi [-h] --bam BAM --reference FASTA --mutations PATH --outdir DIR
             [--voc-column COLUMN] [--mutation-column COLUMN] [-d CHAR]
             [--only-vocs [VOC [VOC ...]]] [-v] [--cite]

optional arguments:
  -h, --help            show this help message and exit
  --bam BAM             Path to a BAM file aligned against the reference
  --reference FASTA     Path to FASTA-formatted complete reference genome
  --mutations PATH      Either a path to tabular file describing all variant
                        mutations or a path to a directory containing PHE-
                        formatted YAML files which each provide the signature
                        mutations for a given variant
  --outdir DIR          Output directory; will be created if it does not
                        already exist
  --voc-column COLUMN   Header for the column of a tabular mutations file
                        containing variant names, and ignored if a YAML
                        directory is provided [PangoLineage]
  --mutation-column COLUMN
                        Header for the column of a tabular mutations file
                        containing mutation descriptions, and ignored if a
                        YAML directory is provided [NucName]
  -d CHAR, --delimiter CHAR
                        Delimiter character for tabular input and output [TAB]
  --only-vocs [VOC [VOC ...]]
                        Look only for these variants, plus the reference
  -v, --version         show program's version number and exit
  --cite                Print the citation in BibTeX format and exit

```

## The Mutations File

The mutations file is a delimited file (tabs, by default) which describe which
mutations belong to which variants of concern. The default headers are
"PangoLineage" for VOC and "NucName" for the mutation description. However,
these may be overridden by the user with `--voc-column` and `--mutation-column`,
respectively. Additional columns are ignored.

All positions are relative to the reference genome and 1-based _i.e.,_ the first
position in the genome is position 1.

As of version `0.10.0`, a directory containing YAML files using the format 
described by Public Health England may be used in place of the mutations file.
The specification for this format can be found
[here](https://github.com/phe-genomics/variant_definitions/). To use these variant
definitions, provide a path to the directory containing the YAML files.

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

It requires that your project directory is structured thusly:

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

# Citing MMMVI

If you have used `mmmvi` in your work, please cite us. A [preprint has been uploaded to bioRxiv](https://www.biorxiv.org/content/10.1101/2021.06.14.448421v1) in
preparation for peer review. If accepted for publication in a peer-reviewed journal, the citation will be updated to reflect that. 

```bibtex
@article {Barker2021.06.14.448421,
	author = {Barker, Dillon O.R. and Buchanan, Cody J. and Landgraff, Chrystal and Taboada, Eduardo N},
	title = {MMMVI: Detecting SARS-CoV-2 Variants of Concern in Metagenomic Samples},
	elocation-id = {2021.06.14.448421},
	year = {2021},
	doi = {10.1101/2021.06.14.448421},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2021/06/15/2021.06.14.448421},
	eprint = {https://www.biorxiv.org/content/early/2021/06/15/2021.06.14.448421.full.pdf},
	journal = {bioRxiv}
}
```

