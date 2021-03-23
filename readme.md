# Identifying Variants of Concern in Metagenomic Samples

This is an in-progress tool to identify potential SARS-CoV-2 Variants of Concern
(VoC) in metagenomic samples.

A VoC is identified by the presence of nucleotide differences relative to a
reference strain. In a heterogeneous, metagenomic sample, sequencing reads may
originate from multiple SARS-CoV-2 genomes. To improve confidence that a VoC is
present, we can look for cases where multiple nucleotide polymorphisms are
present on the same read.

# Install
A proper `conda` installer will be available shortly. 

For the moment, clone the repository and ensure that `mmmvi.py` is available in
your project directory. This is, of course, terrible and will be fixed at the
first opportunity.

A snakemake workflow is provided and is described below. 

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
.
├── reports
│   ├── sample1 
│   │   ├── cooccurence_matrices
│   │   │   ├── B.1.1.7.csv
│   │   │   ├── B.1.351.csv
│   │   │   ├── B.1.525.csv
│   │   │   └── P.1.csv
│   │   ├── read_report.csv
│   │   └── summary.csv
```
