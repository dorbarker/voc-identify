import argparse
import itertools
import pysam
from collections import Counter
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union
import re
import string

from . import __version__

complements = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N", None: None}

Position = Tuple[int, ...]
Mutation = Tuple[Optional[str], ...]
Mutations = Dict[Position, Mutation]
VoCs = Dict[str, Mutations]
Reads = List[pysam.AlignedSegment]
MutationResults = Dict[str, List[Position]]
VoCResults = Dict[str, MutationResults]


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


def is_illumina(reads: Reads) -> bool:
    # Heuristically determine if the reads are paired or not.
    #
    # If duplicated read names outnumber singleton read names by
    # a factor of at least 10:1, then it's Illumina

    counts = Counter(Counter(read.query_name for read in reads).values())

    return (counts[2] / counts[1]) >= 10


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
            vocs[voc][position_range] = set([mutant])

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
        results[variant] = find_variant_mutations(reads, mutations)

    return results


def find_variant_mutations(reads: Reads, mutations: Mutations) -> MutationResults:
    if is_illumina(reads):
        result = find_variant_mutations_illumina(reads, mutations)

    else:
        result = find_variant_mutations_nanopore(reads, mutations)

    return result


def find_variant_mutations_nanopore(
    reads: Reads, mutations: Mutations
) -> MutationResults:
    results = {}

    for read in reads:
        read_name = read.query_name

        seq = read.query_sequence

        pairs = read.get_aligned_pairs()

        results[read_name] = find_mutation_positions(seq, pairs, mutations)

    return results


def find_variant_mutations_illumina(
    reads: Reads, mutations: Mutations
) -> MutationResults:
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


def format_summary(voc_results: VoCResults) -> pd.DataFrame:

    mutation_df = pd.DataFrame(voc_results)

    count_of_reads_with_n_snps = mutation_df.applymap(len).agg(Counter)

    summary = pd.DataFrame(count_of_reads_with_n_snps.to_dict()).transpose()

    return summary


def format_mutation_string(position_range: Position, mutations, wt):

    # filter the Nones in insertions
    no_missing_range = [x for x in position_range if x]

    start = min(no_missing_range)
    stop = max(no_missing_range)

    # deletions
    if mutations[position_range][0] is None:

        if start == stop:
            s = f"[{start}]del"
        else:
            s = f"[{start}-{stop}]del"

    # insertions
    elif is_insertion(position_range):

        insertion_nt = "".join(mutations[position_range])
        s = f"{start + 1}{insertion_nt}"

    # substitutions
    else:

        wildtype_nt = "".join(wt[position_range])
        variant_nt = "".join(mutations[position_range])

        s = f"{wildtype_nt}{start + 1}{variant_nt}"

    return s


def format_cooccurrence_matrix(mutation_result, mutations, wt) -> pd.DataFrame:
    # For one VoC at a time

    lookup = {
        position_range: format_mutation_string(position_range, mutations, wt)
        for position_range in mutations.keys()
    }

    mx = pd.DataFrame(data=0, index=lookup.values(), columns=lookup.values())

    # for each mutation position
    for positions in mutation_result.values():
        # self-vs-self
        for position in positions:
            name = lookup[position]

            mx.loc[name, name] += 1

        for row, col in itertools.permutations(positions, r=2):
            row_name = lookup[row]
            col_name = lookup[col]

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


def format_read_species(
    reads: Reads, mutation_results: VoCResults, read_report: pd.DataFrame, vocs: VoCs,
) -> pd.DataFrame:
    species = {}

    total_reads = len(mutation_results["reference"].keys())

    for _, variant_positions in read_report.iterrows():

        position_nts = set()

        for variant, positions in zip(read_report.columns, variant_positions):

            for position_range in positions:
                # convert between 1-based read_report and 0-based vocs
                voc_pos = tuple(position - 1 for position in position_range)

                try:
                    voc_nt = vocs[variant][voc_pos]
                except KeyError:
                    voc_pos = (voc_pos[0], None, voc_pos[0] + 1)
                    voc_nt = (vocs[variant][voc_pos],)

                    position_range = (position_range[0],)

                position_nt = (tuple(position_range), voc_nt)

                position_nts.add(position_nt)

        if not position_nts:  # reads matching no VoCs
            continue

        pairs = sorted(position_nts, key=lambda x: x[0])
        key = str(tuple(f"{p}{nt}" for p, nt in pairs))

        try:
            species[key]["count"] += 1

        except KeyError:

            locations, nucleotides = zip(*pairs)

            # a temporary fix to a weird regression introduced by multibase subs/dels
            locations = tuple(itertools.chain.from_iterable(locations))
            nucleotides = tuple(itertools.chain.from_iterable(nucleotides))

            species[key] = {
                "positions": locations,
                "nucleotides": tuple("del" if nt is None else nt for nt in nucleotides),
                "count": 1,  # the number of times this combination of mutations has been observed
            }

            # add the bitarrays
            species[key].update(make_voc_bitarray(locations, nucleotides, vocs))

    read_species = pd.DataFrame.from_dict(species, orient="index")

    read_species["proportion_total"] = read_species["count"] / total_reads

    overlapping_counts = read_species_overlap(read_species, reads)

    read_species["reads_overlapping"] = [
        overlapping_counts[positions] for positions in read_species["positions"]
    ]

    read_species["proportion_overlapping"] = (
        read_species["count"] / read_species["reads_overlapping"]
    )

    return read_species


def make_voc_bitarray(
    locations: Tuple[int, ...],
    nucleotides: Tuple[Union[str, Tuple[str, ...]], ...],
    vocs: VoCs,
) -> Dict[str, Tuple[int, ...]]:
    # locations is 1-indexed
    voc_bitarrays = {}

    locs = tuple(p - 1 for p in locations)  # back to 0-index

    are_insertions = [isinstance(x, tuple) for x in nucleotides]

    for variant in vocs:

        bitarray = []

        voc_positions, voc_nts = [], []

        for positions, mutations in vocs[variant].items():

            if len(positions) == len(mutations):

                for p, m in zip(positions, mutations):
                    voc_positions.append(p)
                    voc_nts.append(m)

            # insertion handling; otherwise we get off-by N errors for the reference
            # because a location is never None, this doesn't mess with .index() below
            else:
                for m in mutations:
                    voc_positions.append(None)
                    voc_nts.append(m)

        for loc, nt, insertion in zip(locs, nucleotides, are_insertions):

            if insertion:

                try:
                    match = int(vocs[variant][(loc, None, loc + 1)] == nt)
                except KeyError:
                    match = 0

            else:

                try:

                    idx = voc_positions.index(loc)
                    match = int(voc_nts[idx] == nt)

                except ValueError:
                    match = 0

            bitarray.append(match)

        voc_bitarrays[variant] = tuple(bitarray)

    return voc_bitarrays


def read_species_overlap(
    read_species: pd.DataFrame, reads: Reads
) -> Dict[Tuple[int, ...], int]:
    overlapping_counts = {species: 0 for species in read_species["positions"]}

    for read in reads:

        ref_positions = set(read.get_reference_positions())

        for species_positions in overlapping_counts:

            is_overlapping = all(p in ref_positions for p in species_positions)

            overlapping_counts[species_positions] += is_overlapping

    return overlapping_counts


def format_reports(reads: Reads, voc_results: VoCResults, vocs: VoCs):
    oir_results = one_index_results(voc_results)

    reports = {
        "read_report": format_read_report(oir_results),
        "summary": format_summary(voc_results),
        "absolute_cooccurrence_matrices": format_cooccurrence_matrices(
            voc_results, vocs
        ),
    }
    reports["read_species"] = format_read_species(
        reads, voc_results, reports["read_report"], vocs
    )

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
