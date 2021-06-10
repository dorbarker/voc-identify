import argparse
import itertools

import numpy as np
import pysam
from collections import Counter
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import re
import string
import logging
import statistics

from . import __version__

complements = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N", None: None}

Position = Tuple[int, ...]
Mutation = Tuple[Optional[str], ...]
Mutations = Dict[Position, Mutation]
VoCs = Dict[str, Mutations]
Reads = List[pysam.AlignedSegment]
MutationResults = Dict[str, List[Position]]
VoCResults = Dict[str, MutationResults]


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

    reports = format_reports(reads, mutation_results, vocs)

    write_reports(reports, args.outdir, args.delimiter)


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


def one_index_range(position_mutation):

    position_range, mutation = position_mutation

    if None in position_range:
        result = [position_range[0] + 1]
    else:
        result = [pos + 1 for pos in position_range]

    return result


def one_index_results(voc_results: VoCResults) -> VoCResults:

    oir = (
        pd.DataFrame(voc_results)
        .applymap(lambda cell: [one_index_range(group) for group in cell])
        .to_dict()
    )

    return oir


def format_read_report(oir_results: VoCResults) -> pd.DataFrame:
    read_report = pd.DataFrame(oir_results)

    has_any_results = read_report.applymap(len).apply(sum, axis="columns") > 0

    return read_report[has_any_results]


def format_summary(voc_results: VoCResults, vocs: VoCs, reads) -> pd.DataFrame:

    mutation_df = pd.DataFrame(voc_results)

    count_of_reads_with_n_snps = mutation_df.applymap(len).agg(Counter)

    summary = (
        pd.DataFrame(count_of_reads_with_n_snps.to_dict(),)
        .transpose()
        .fillna(0)
        .applymap(int)
    )

    max_coverage = theoretical_maximum(reads, vocs)
    signature_counts = mutation_coverage(voc_results, vocs)

    return summary.join(max_coverage).join(signature_counts)


def shannon_entropy():
    pass


def mutation_coverage(voc_results, vocs):

    nonexclusive = {}
    exclusive = {}
    report = {}

    for voc, mutation_results in voc_results.items():

        signatures = set(vocs[voc])

        present_mutations = set(
            itertools.chain.from_iterable(mutation_results.values())
        )

        nonexclusive[voc] = {
            "maximum": len(signatures),
            "present": present_mutations,
            "signatures": signatures,
        }

    for voc, data in nonexclusive.items():

        other_sigs = set()
        for v in nonexclusive:
            if v not in {voc, "reference"}:
                other_sigs.update(nonexclusive[v]["signatures"])

        exclusive_sigs = data["signatures"].difference(other_sigs)

        # janky - fix
        exclusive_present = set(p for (p, m) in data["present"]).intersection(
            exclusive_sigs
        )

        exclusive[voc] = {
            "maximum": len(exclusive_sigs),
            "present": exclusive_present,
        }

    for voc in nonexclusive.keys():

        ne_num = len(nonexclusive[voc]["present"])
        ne_denom = nonexclusive[voc]["maximum"]

        e_num = len(exclusive[voc]["present"])
        e_denom = exclusive[voc]["maximum"]

        report[voc] = {
            "complete_signature_mutations": f"{ne_num}/{ne_denom}",
            "exclusive_signature_mutations": f"{e_num}/{e_denom}",
        }

    return pd.DataFrame.from_dict(report, orient="index")


def theoretical_maximum(reads: Reads, vocs: VoCs):

    median_read_length = statistics.median(
        [read.query_alignment_length for read in reads]
    )

    voc_max = {}
    for voc in vocs:

        position_ranges = sorted(vocs[voc].keys())

        max_covered = 0
        for position_range in position_ranges:

            start = position_range[0]
            end = start + median_read_length

            n_covered = sum([p[0] >= start and p[-1] <= end for p in position_ranges])
            max_covered = max(n_covered, max_covered)

        voc_max[voc] = {"median_maximum_coverage": max_covered}

    return pd.DataFrame.from_dict(voc_max, orient="index")


def format_mutation_string(position_range, mutation, wt):

    start = position_range[0] + 1  # adjust to 1-based counting for reporting

    # insertion
    if None in position_range:

        insertion_nt = "".join(mutation)
        s = f"{start}{insertion_nt}"

    # deletion
    elif None in mutation:

        stop = position_range[-1] + 1

        # point deletion
        if start == stop:
            s = f"[{start}]del"

        # multi-base deletion
        else:
            s = f"[{start}-{stop}]del"

    # substitution
    else:

        wildtype = wt[position_range][0]

        wildtype_nt = "".join(wildtype)
        variant_nt = "".join(mutation)

        s = f"{wildtype_nt}{start}{variant_nt}"

    return s


def initialize_matrix(voc, wt):

    lookup = {}
    mutation_strings = []

    for position_range, mutations in voc.items():

        lookup[position_range] = {}

        for mutation in mutations:

            mutation_string = format_mutation_string(position_range, mutation, wt)

            lookup[position_range][mutation] = mutation_string

            mutation_strings.append(mutation_string)

    mx = pd.DataFrame(data=0, index=mutation_strings, columns=mutation_strings)

    return lookup, mx


def format_cooccurrence_matrix(mutation_result, voc, wt):
    # For one VoC at a time

    lookup, mx = initialize_matrix(voc, wt)

    for read_mutations in mutation_result.values():

        for position, mutation in read_mutations:

            name = lookup[position][mutation]

            mx.loc[name, name] += 1

        for (row_pos, row_mut), (col_pos, col_mut) in itertools.permutations(
            read_mutations, r=2
        ):

            row_name = lookup[row_pos][row_mut]
            col_name = lookup[col_pos][col_mut]

            mx.loc[row_name, col_name] += 1

    return mx


def format_relative_cooccurrence_matrix(
    cooccurrence_matrix: pd.DataFrame,
) -> pd.DataFrame:

    rows = []

    for denominator_name in cooccurrence_matrix.columns:

        denominator = cooccurrence_matrix.loc[denominator_name, denominator_name]

        for numerator_name in cooccurrence_matrix.index:

            numerator = cooccurrence_matrix.loc[numerator_name, denominator_name]

            try:
                quotient = int(numerator) / int(denominator)
            except ZeroDivisionError:
                quotient = 0.0

            rows.append(
                {
                    "denominator": denominator_name,
                    "numerator": numerator_name,
                    "ratio": quotient,
                }
            )

    return pd.DataFrame(rows)


def format_relative_cooccurrence_matrices(absolute_cooccurrence_matrices):
    return {
        v: format_relative_cooccurrence_matrix(mx)
        for v, mx in absolute_cooccurrence_matrices.items()
    }


def format_cooccurrence_matrices(voc_results: VoCResults, vocs: VoCs):
    *variants, wt = sorted(vocs.keys(), key=lambda x: x == "reference")

    return {
        v: format_cooccurrence_matrix(voc_results[v], vocs[v], vocs[wt])
        for v in variants
    }


def format_read_species(voc_results, vocs, reads):

    species = {}
    total_reads = len(voc_results["reference"].keys())

    # get non-redundant set of positions across VOCs

    for key, species_data in nonredundant_read_species(voc_results):

        positions_mutations = species_data["positions_mutations"]

        species_positions, species_mutations = format_positions_mutations(
            positions_mutations
        )

        species[key] = {
            "positions": species_positions,
            "nucleotides": species_mutations,
            "count": species_data["count"],
        }

        bitarrays = make_voc_bitarray(positions_mutations, vocs)
        species[key].update(bitarrays)

    read_species = pd.DataFrame.from_dict(species, orient="index")

    read_species["proportion_total"] = read_species["count"] / total_reads

    overlapping_counts = read_species_overlap(read_species["positions"], reads)

    read_species["reads_overlapping"] = [
        overlapping_counts[positions] for positions in read_species["positions"]
    ]

    read_species["proportion_overlapping"] = (
        read_species["count"] / read_species["reads_overlapping"]
    )

    return read_species.sort_values(by="count", ascending=False)


def nonredundant_read_species(voc_results):

    nonredundant_reads = set()
    for read_results in voc_results.values():
        for read in read_results.keys():
            nonredundant_reads.add(read)

    nonredundant = {}
    for read in nonredundant_reads:
        positions_mutations = []
        for voc in voc_results:
            positions_mutations.extend(voc_results[voc][read])
        positions_mutations.sort()

        if not positions_mutations:
            continue

        key = str(positions_mutations)

        try:
            nonredundant[key]["count"] += 1

        except KeyError:
            nonredundant[key] = {
                "positions_mutations": positions_mutations,
                "count": 1,
            }

        yield from nonredundant.items()


def format_positions_mutations(positions_mutations):

    species_positions = []
    species_mutations = []

    for p, m in positions_mutations:

        # insertion
        if None in p:
            species_positions.append((p[0],))
            species_mutations.append(tuple("del" if x is None else x for x in m))

        # deletion
        elif None in m:
            species_positions.append(p)
            species_mutations.append(tuple("del" for _ in m))

        # substitution
        else:
            species_positions.append(p)
            species_mutations.append(m)

    species_positions = tuple(
        tuple(p + 1 for p in group) for group in species_positions
    )

    species_mutations = tuple(species_mutations)

    return species_positions, species_mutations


def make_voc_bitarray(positions_mutations, vocs: VoCs) -> Dict[str, Tuple[int, ...]]:
    # Returns a dictionary with a value of a bit array for each VOC key
    #
    # The bit array show for each position in the read species if it belongs (1) or not (0)
    # to each VOC or the reference.
    #
    # Consider a read species which covers a potential of 5 mutations which are diagnostic
    # for VOCs. In this case, the first position happens to be wild type (reference),
    # but the following 4 positions are diagnostic of the "A" variant. However, the "B" variant
    # shares the mutations at the 4th and 5th positions with "A". The results would be:
    #
    # {"reference": (1, 0, 0, 0, 0), "A": (0, 1, 1, 1, 1,), "B": (0, 0, 0, 1, 1)}

    bitarrays = {}
    for (position, nts), voc in itertools.product(positions_mutations, vocs):
        match = int(position in vocs[voc] and nts in vocs[voc][position])
        try:
            bitarrays[voc].append(match)
        except KeyError:
            bitarrays[voc] = [match]

    return {k: tuple(v) for k, v in bitarrays.items()}


def read_species_overlap(
    positions: pd.Series, reads: Reads
) -> Dict[Tuple[int, ...], int]:
    # Calculates the number of reads which overlap a read species.
    #
    # To be considered overlapping, the read must contain
    # all of the positions in the species.
    overlapping_counts = {species: 0 for species in positions}

    for read in reads:

        read_start, *_, read_end = sorted(read.get_reference_positions())

        for species_positions in overlapping_counts:

            sorted_positions = sorted(itertools.chain.from_iterable(species_positions))

            try:
                start, *_, stop = sorted_positions
                is_overlapping = start >= read_start and stop <= read_end
            # a species with a single point mutation
            except ValueError:
                start, *_ = sorted_positions
                is_overlapping = read_end >= start >= read_start

            overlapping_counts[species_positions] += is_overlapping

    return overlapping_counts


def format_reports(reads: Reads, voc_results: VoCResults, vocs: VoCs):

    logging.info("Formatting reports")

    oir_results = one_index_results(voc_results)

    reports = {
        "read_report": format_read_report(oir_results),
        "summary": format_summary(voc_results, vocs, reads),
        "absolute_cooccurrence_matrices": format_cooccurrence_matrices(
            voc_results, vocs
        ),
        "read_species": format_read_species(voc_results, vocs, reads),
    }

    reports["relative_cooccurrence_matrices"] = format_relative_cooccurrence_matrices(
        reports["absolute_cooccurrence_matrices"]
    )
    return reports


def write_cooccurrence_matrix(
    variant: str, directory: Path, data: pd.DataFrame, delimiter: str
) -> None:
    variant_out_name = variant.replace("/", "_")
    p = directory.joinpath(f"{variant_out_name}.txt")
    data.to_csv(p, sep=delimiter)


def write_reports(reports, outdir: Path, delimiter: str):

    logging.info("Writing reports")

    matrices_path = outdir.joinpath("cooccurrence_matrices")

    absolute_matrices = matrices_path.joinpath("absolute")
    absolute_matrices.mkdir(parents=True, exist_ok=True)

    relative_matrices = matrices_path.joinpath("relative")
    relative_matrices.mkdir(parents=True, exist_ok=True)

    reports["read_report"].to_csv(outdir / "read_report.txt", sep=delimiter)

    reports["summary"].to_csv(outdir / "summary.txt", sep=delimiter)

    reports["read_species"].to_csv(
        outdir / "read_species.txt", sep=delimiter, index=False
    )

    for variant, data in reports["absolute_cooccurrence_matrices"].items():
        write_cooccurrence_matrix(variant, absolute_matrices, data, delimiter)

    for variant, data in reports["relative_cooccurrence_matrices"].items():
        write_cooccurrence_matrix(variant, relative_matrices, data, delimiter)


if __name__ == "__main__":
    main()
