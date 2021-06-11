import argparse
import itertools
import logging
import re
import string

from pathlib import Path
from typing import List, Tuple, Optional

import pysam
import pandas as pd

from mmmvi import __version__
from mmmvi import reporting
from mmmvi.types import Reads, VoCs, VoCResults, MutationResults, Mutations, Position

complements = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N", None: None}


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

    vocs = load_mutations(
        args.mutations,
        args.reference,
        args.voc_column,
        args.mutation_column,
        args.delimiter,
        args.only_vocs,
    )

    reads = load_reads(args.bam, args.reference)

    mutation_results = find_mutations(reads, vocs)

    reports = reporting.format_reports(reads, mutation_results, vocs)

    reporting.write_reports(reports, args.outdir, args.delimiter)

    logging.info("Complete")


def load_reference(reference: Path) -> str:
    # reference needs to be complete and in a single contig anyway
    #
    # written to avoid having to pull in all of biopython

    lines = []
    with reference.open("r") as f:
        for line in f:
            if not line.startswith(">"):
                lines.append(line.strip())

    seq = "".join(lines)

    return seq


def parse_mutation(s: str):

    if s.endswith("del"):

        position_range, wt, mutation = parse_deletion(s)

    elif s[0] in string.ascii_uppercase:

        position_range, wt, mutation = parse_substitution(s)

    else:
        position_range, wt, mutation = parse_insertion(s)

    return position_range, wt, mutation


def parse_deletion(s: str):
    # [123-125]del means that reference positions 123, 124, and 125 are deleted in the read

    _, *start_stop, _ = re.split(r"[\[\-\]]", s)

    try:
        start, stop = start_stop
    except ValueError:
        start = stop = start_stop[0]

    if stop < start:
        raise ValueError(f"stop is less than start in {s}")

    start = int(start)
    stop = int(stop)

    position_range = tuple(range(start - 1, stop))
    mutation = tuple(None for _ in position_range)
    wt = (None,)

    return position_range, wt, mutation


def parse_substitution(s: str):
    # A123T means A in the reference has been substituted by T in read
    # CAT123GTA means C, A, T at positions 123, 124, 125 have been substituted by G, T,  A

    wt, mutation = (tuple(x) for x in re.findall(r"[ATCG]+", s))

    if len(wt) != len(mutation):
        wt_str = "".join(wt)
        mut_str = "".join(mutation)
        raise ValueError(
            f"Mismatch between length of wild type '{wt_str}' and mutant '{mut_str}' in '{s}'"
        )

    start = int(re.search(r"\d+", s).group()) - 1
    position_range = tuple(range(start, start + len(wt)))

    return position_range, wt, mutation


def parse_insertion(s: str):
    # 123CAT indicates CAT has been inserted betwixt reference positions 123 and 124

    position = int("".join(itertools.takewhile(lambda x: x in string.digits, s)))

    # Exactly 1 None will get the whole insertion by exploiting pd.Series indexing
    position_range = (position - 1, None, position)

    mutation = tuple("".join(itertools.dropwhile(lambda x: x in string.digits, s)))

    wt = tuple(None for _ in mutation)

    return position_range, wt, mutation


def load_mutations(
    mutations_path: Path,
    reference_path: Path,
    voc_col: str,
    mut_col: str,
    delimiter: str,
    selected_vocs: List[str],
) -> VoCs:
    data = pd.read_csv(mutations_path, sep=delimiter)

    reference_seq = load_reference(reference_path)

    vocs = {"reference": {}}

    for idx, row in data.iterrows():

        voc = row[voc_col]

        if selected_vocs and voc not in selected_vocs:
            continue

        position_range, wt, mutant = parse_mutation(row[mut_col])

        if voc not in vocs:
            vocs[voc] = {}

        try:
            vocs[voc][position_range].add(mutant)
        except KeyError:
            vocs[voc][position_range] = {mutant}

        if wt == (None,):
            wt = tuple(reference_seq[position] for position in position_range)

        vocs["reference"][position_range] = [wt]

    return vocs


def load_reads(bam_path: Path, ref_path: Path) -> Reads:

    logging.info(f"Loading reads from {bam_path}")

    with pysam.AlignmentFile(
        bam_path, reference_filename=str(ref_path), mode="rb"
    ) as aln:
        return list(aln)


def find_mutations(reads: Reads, vocs: VoCs) -> VoCResults:
    results = {}

    for variant, mutations in vocs.items():
        logging.info(f"Searching for {variant} signature mutations")

        results[variant] = find_variant_mutations(reads, mutations)

    return results


def find_variant_mutations(reads: Reads, mutations: Mutations) -> MutationResults:
    results = {}

    for read in reads:

        orientation_tag = "rev" if read.is_reverse else "fwd"

        read_name = f"{read.query_name}:{orientation_tag}"

        seq = read.query_sequence

        pairs = read.get_aligned_pairs()

        results[read_name] = find_mutation_positions(seq, pairs, mutations)

    return results


def pad_seq_with_ambiguous(
    seq: str, query_positions: List[Optional[int]]
) -> List[Optional[str]]:
    new_seq = [None for _ in query_positions]

    for seq_element, query_position in zip(range(len(new_seq)), query_positions):
        try:
            nt = seq[query_position]
        except TypeError:  # None in the query positions
            continue
        new_seq[seq_element] = nt

    return new_seq


def is_insertion(position_range: Tuple[Optional[int], ...]) -> bool:

    try:
        result = position_range[1] is None and isinstance(position_range[0], int)
    except IndexError:
        result = False

    return result


def find_mutation_positions(seq: str, pairs, mutations) -> List[Position]:

    mutated_regions = []

    query_positions, subject_positions = zip(*pairs)

    aln = pd.Series(
        pad_seq_with_ambiguous(seq, query_positions), index=subject_positions
    )

    for mutation_positions, mutation_seqs in mutations.items():

        for mutation_seq in mutation_seqs:
            # has all of the mutations in the current group
            relevant = all(p in aln for p in mutation_positions)

            if not relevant:
                continue

            if is_insertion(mutation_positions):
                try:

                    start, _, stop = mutation_positions
                    has_mutation = aln.loc[start:stop][[None]]

                # This read spans the insertion locus, but doesn't actually have the insertion
                except KeyError:

                    # if the current 'VOC' is wild type
                    if all(x is None for x in mutation_seq):

                        mutated_regions.append((mutation_positions, mutation_seq))

                    continue
            else:
                has_mutation = aln.loc[list(mutation_positions)]

            try:

                is_mutated = has_mutation.equals(
                    pd.Series(mutation_seq, index=has_mutation.index)
                )

                if is_mutated:
                    mutated_regions.append((mutation_positions, mutation_seq))

            except ValueError:

                # spurious insertions, especially in Nanopore data
                if len(has_mutation) != len(mutation_seq):
                    continue

    return mutated_regions


if __name__ == "__main__":
    main()
