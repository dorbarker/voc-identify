import argparse
import logging
import sys

from pathlib import Path

from mmmvi import __version__, CitationAction
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
        metavar="PATH",
        help="""Either a path to tabular file describing all variant 
                mutations or a path to a directory containing PHE-formatted
                YAML files which each provide the signature mutations for a
                given variant""",
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
        help="""Header for the column of a tabular mutations file containing 
                variant names, and ignored if a YAML directory is provided 
                [PangoLineage]""",
    )

    parser.add_argument(
        "--mutation-column",
        default="NucName",
        metavar="COLUMN",
        help="""Header for the column of a tabular mutations file containing 
                mutation descriptions, and ignored if a YAML directory is 
                provided [NucName]""",
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

    parser.add_argument(
        "--cite",
        nargs=0,
        action=CitationAction,
        help="Print the citation in BibTeX format and exit",
    )

    args = parser.parse_args()

    return args


def main():

    args = arguments()
    logging.info("Begin")

    try:
        vocs = load_data.load_mutations(
            args.mutations,
            args.reference,
            args.voc_column,
            args.mutation_column,
            args.delimiter,
            args.only_vocs,
        )
    except FileNotFoundError as fne:
        logging.error(str(fne))
        sys.exit(1)

    reads = load_data.load_reads(args.bam, args.reference)

    mutation_results = search.find_mutations(reads, vocs)

    reporting.write_reports(mutation_results, vocs, reads, args.outdir, args.delimiter)

    logging.info("Complete")


if __name__ == "__main__":
    main()
