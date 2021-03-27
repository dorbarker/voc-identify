# Identifying Variants of Concern in Metagenomic Samples

This is an in-progress tool to identify potential SARS-CoV-2 Variants of Concern
(VoC) in metagenomic samples.

A VoC is identified by the presence of nucleotide differences relative to a
reference strain. In a heterogeneous, metagenomic sample, sequencing reads may
originate from multiple SARS-CoV-2 genomes. To improve confidence that a VoC is
present, we can look for cases where multiple nucleotide polymorphisms are
present on the same read.

# Install
```sh
conda install -c conda-forge -c bioconda -c dorbarker mmmvi
```

# Running it

A [snakemake](https://snakemake.readthedocs.io/en/stable/) workflow,
`voc-identify.smk` is provided for convenience and easier parallel processing.

It assumes your project directory is structured thusly:

```
.
├── bams
│   ├── sample1.bam 
│   └── sample2.bam
|   |   ...
│   └── sampleN.bam
├── data
    └── mutations.tsv
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
