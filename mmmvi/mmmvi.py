import argparse
import logging

from pathlib import Path

from mmmvi import __version__
from mmmvi.lib import load_data, reporting, search

logging.basicConfig(
    format="%(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S", level=logging.INFO
)


def arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--bam",
        required=True,
        type=Path,
        metavar="BAM",
        help="Path to a BAM file aligned against the reference",
    )

    parser.add_argument(
        "--reference",
        required=True,
        type=Path,
        metavar="FASTA",
        help="Path to FASTA-formatted complete reference genome",
    )

    parser.add_argument(
        "--mutations",
        required=True,
        type=Path,
        metavar="TABULAR",
        help="Path to tabular file describing Variants of Concern",
    )

    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        metavar="DIR",
        help="Output directory; will be created if it does not already exist",
    )

    parser.add_argument(
        "--voc-column",
        default="PangoLineage",
        metavar="COLUMN",
        help="Header for the column containing Variant of Concern names [PangoLineage]",
    )

    parser.add_argument(
        "--mutation-column",
        default="NucName",
        metavar="COLUMN",
        help="Header for the column containing mutation descriptions [NucName]",
    )

    parser.add_argument(
        "-d",
        "--delimiter",
        default="\t",
        metavar="CHAR",
        help="Delimiter character for tabular input and output [TAB]",
    )

    parser.add_argument(
        "--only-vocs",
        nargs="*",
        metavar="VOC",
        default=[],
        help="Look only for these variants, plus the reference",
    )

    parser.add_argument(
        "-v", "--version", action="version", version=f"{parser.prog} {__version__}"
    )

    return parser.parse_args()


def main():

    args = arguments()
    logging.info("Begin")

    vocs = load_data.load_mutations(
        args.mutations,
        args.reference,
        args.voc_column,
        args.mutation_column,
        args.delimiter,
        args.only_vocs,
    )

    reads = load_data.load_reads(args.bam, args.reference)

    mutation_results = search.find_mutations(reads, vocs)

    reporting.write_reports(mutation_results, vocs, reads, args.outdir, args.delimiter)

    logging.info("Complete")


if __name__ == "__main__":
    main()
