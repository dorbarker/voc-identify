import argparse
import itertools
import pysam
from collections import Counter
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple
import re

from . import __version__

complements = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N", None: None}

Mutations = Dict[Tuple[int], Tuple[str]]
VoCs = Dict[str, Mutations]
Reads = List[pysam.AlignedSegment]
MutationResults = Dict[str, List[int]]


def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument("--bam", required=True, type=Path)

    parser.add_argument("--reference", required=True, type=Path)

    parser.add_argument("--mutations", required=True, type=Path)

    parser.add_argument("--outdir", required=True, type=Path)

    parser.add_argument(
        "-d",
        "--delimiter",
        default="\t",
        help="Delimiter character for tabular input and output [TAB]",
    )

    parser.add_argument(
        "-v", "--version", action="version", version=f"{parser.prog} {__version__}"
    )

    return parser.parse_args()


def main():

    args = arguments()

    vocs = load_mutations(args.mutations, args.reference, args.delimiter)

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
                lines.append(line)

    seq = "".join(lines)

    return seq


def parse_mutation(s):

    if s.endswith("del"):
        _, start, stop, _ = re.split("[\[\-\]]", s)
        position_range = tuple(range(int(start) - 1, int(stop)))
        mutation = tuple(None for _ in position_range)
        wt = (None,)

    else:
        wt, mutation = (tuple(x) for x in re.findall("[ATCG]+", s))

        start = int(re.search("\d+", s).group()) - 1
        position_range = tuple(range(start, start + len(wt)))

    return position_range, wt, mutation


def load_mutations(mutations_path: Path, reference_path: Path, delimiter: str) -> VoCs:

    data = pd.read_csv(mutations_path, sep=delimiter)

    reference_seq = load_reference(reference_path)

    vocs = {"reference": {}}

    for idx, row in data.iterrows():

        voc = row["PangoLineage"]

        position_range, wt, mutant = parse_mutation(row["NucName"])

        if voc not in vocs:
            vocs[voc] = {}

        vocs[voc][position_range] = mutant

        if wt == (None,):
            wt = tuple(reference_seq[position] for position in position_range)

        vocs["reference"][position_range] = wt

    return vocs


def load_reads(bam_path: Path, ref_path: Path) -> Reads:

    with pysam.AlignmentFile(
        bam_path, reference_filename=str(ref_path), mode="rb"
    ) as aln:
        return list(aln)


def find_mutations(reads: Reads, vocs: VoCs) -> Dict[str, MutationResults]:

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

        results[read_name] = find_mutation_positions(seq, pairs, False, mutations)

    return results


def find_variant_mutations_illumina(
    reads: Reads, mutations: Mutations
) -> MutationResults:

    results = {}

    for read in reads:

        revcomp = read.is_reverse
        orientation_tag = "rev" if revcomp else "fwd"
        read_name = f"{read.query_name}:{orientation_tag}"

        seq = read.get_forward_sequence()

        pairs = read.get_aligned_pairs()

        results[read_name] = find_mutation_positions(seq, pairs, revcomp, mutations)

    return results


def pad_seq_with_ambiguous(seq, query_positions):

    new_seq = [None for _ in query_positions]

    for seq_element, query_position in zip(range(len(new_seq)), query_positions):
        try:
            nt = seq[query_position]
        except TypeError:  # None in the query positions
            continue
        new_seq[seq_element] = nt

    return new_seq


def find_mutation_positions(seq, pairs, revcomp, mutations):

    mutated_regions = []

    original_orientation = True

    query_positions, subject_positions = zip(*pairs)

    aln = pd.DataFrame(
        {
            "seq": pad_seq_with_ambiguous(seq, query_positions),
            "subject_positions": subject_positions,
        }
    )

    for mutation_positions, mutation_seq in mutations.items():

        has_mutation = aln["subject_positions"].isin(mutation_positions)

        # has all of the mutations in the current group
        relevant = has_mutation.sum() == len(mutation_positions)

        if not relevant:
            continue

        # punt the reverse complementation of the sequence down here,
        # it's relatively costly do on tens of thousands of reads
        # if you don't need to

        # so we don't keep flipping orientations if multiple mutations hit the read
        if revcomp and original_orientation:
            aln["seq"] = [complements[nt] for nt in reversed(aln["seq"])]
            original_orientation = False

        current_seq = aln.loc[has_mutation, "seq"]

        is_mutated = current_seq.equals(
            pd.Series(mutation_seq, index=current_seq.index)
        )

        if is_mutated:
            mutated_regions.append(mutation_positions)

    return mutated_regions


def one_index_results(
    mutation_results: Dict[str, MutationResults]
) -> Dict[str, MutationResults]:

    oir = (
        pd.DataFrame(mutation_results)
        .applymap(lambda cell: [[pos + 1 for pos in group] for group in cell])
        .to_dict()
    )
    return oir


def format_read_report(oir_results: Dict[str, MutationResults]) -> pd.DataFrame:

    read_report = pd.DataFrame(oir_results)

    has_any_results = read_report.applymap(len).apply(sum, axis="columns") > 0

    return read_report[has_any_results]


def format_summary(mutation_results):

    mutation_df = pd.DataFrame(mutation_results)

    count_of_reads_with_N_snps = mutation_df.applymap(len).agg(Counter)

    return pd.DataFrame(count_of_reads_with_N_snps.to_dict()).transpose()


def format_mutation_string(position_range, mutations, wt):

    start = min(position_range)
    stop = max(position_range)

    if mutations[position_range][0] == None:

        if start == stop:
            s = f"[{start}]del"
        else:
            s = f"[{start}-{stop}]del"

    else:

        wildtype_nt = "".join(wt[position_range])
        variant_nt = "".join(mutations[position_range])

        s = f"{wildtype_nt}{start + 1}{variant_nt}"

    return s


def format_cooccurence_matrix(mutation_result, mutations, wt) -> pd.DataFrame:
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


def format_relative_coocurence_matrix(coocurence_matrix: pd.DataFrame) -> pd.DataFrame:

    rows = []

    for denominator_name in coocurence_matrix.columns:

        denominator = coocurence_matrix.loc[denominator_name, denominator_name]

        for numerator_name in coocurence_matrix.index:

            numerator = coocurence_matrix.loc[numerator_name, denominator_name]

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


def format_relative_coocurence_matrices(absolute_coocurrence_matrices):

    return {
        v: format_relative_coocurence_matrix(mx)
        for v, mx in absolute_coocurrence_matrices.items()
    }


def format_cooccurence_matrices(mutation_results, vocs):

    *variants, wt = sorted(vocs.keys(), key=lambda x: x == "reference")

    return {
        v: format_cooccurence_matrix(mutation_results[v], vocs[v], vocs[wt])
        for v in variants
    }


def format_read_species(
    reads: Reads,
    mutation_results: MutationResults,
    read_report: pd.DataFrame,
    vocs: VoCs,
) -> pd.DataFrame:

    species = {}

    total_reads = len(mutation_results["reference"].keys())

    for _, variant_positions in read_report.iterrows():

        position_nts = set()

        matching_variant = {v: 0 for v in vocs}

        for variant, positions in zip(read_report.columns, variant_positions):

            for position_range in positions:

                matching_variant[variant] += 1

                # convert between 1-based read_report and 0-based vocs
                voc_pos = tuple(position - 1 for position in position_range)

                position_nt = (tuple(position_range), vocs[variant][voc_pos])

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
                "nucleotides": nucleotides,
                "count": 1,
            }
            species[key].update(matching_variant)

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


def read_species_overlap(
    read_species: pd.DataFrame, reads: Reads
) -> Dict[Tuple[int, ...], int]:

    overlapping_counts = {species: 0 for species in read_species["positions"]}

    for read in reads:

        ref_positions = set(read.get_reference_positions())
        # print("-" * 10)
        # print(ref_positions)
        for species_positions in overlapping_counts:
            # print(species_positions)
            is_overlapping = all(p in ref_positions for p in species_positions)

            overlapping_counts[species_positions] += is_overlapping

    return overlapping_counts


def format_reports(reads: Reads, mutation_results, vocs):

    oir_results = one_index_results(mutation_results)

    reports = {
        "read_report": format_read_report(oir_results),
        "summary": format_summary(mutation_results),
        "absolute_cooccurence_matrices": format_cooccurence_matrices(
            mutation_results, vocs
        ),
    }
    reports["read_species"] = format_read_species(
        reads, mutation_results, reports["read_report"], vocs
    )

    reports["relative_cooccurence_matrices"] = format_relative_coocurence_matrices(
        reports["absolute_cooccurence_matrices"]
    )
    return reports


def write_cooccurence_matrix(
    variant: str, directory: Path, data: pd.DataFrame, delimiter: str
) -> None:

    variant_out_name = variant.replace("/", "_")
    p = directory.joinpath(f"{variant_out_name}.txt")
    data.to_csv(p, sep=delimiter)


def write_reports(reports, outdir: Path, delimiter: str):

    matrices_path = outdir.joinpath("cooccurence_matrices")

    absolute_matrices = matrices_path.joinpath("absolute")
    absolute_matrices.mkdir(parents=True, exist_ok=True)

    relative_matrices = matrices_path.joinpath("relative")
    relative_matrices.mkdir(parents=True, exist_ok=True)

    reports["read_report"].to_csv(outdir / "read_report.txt", sep=delimiter)

    reports["summary"].to_csv(outdir / "summary.txt", sep=delimiter)

    reports["read_species"].to_csv(
        outdir / "read_species.txt", sep=delimiter, index=False
    )

    for variant, data in reports["absolute_cooccurence_matrices"].items():

        write_cooccurence_matrix(variant, absolute_matrices, data, delimiter)

    for variant, data in reports["relative_cooccurence_matrices"].items():

        write_cooccurence_matrix(variant, relative_matrices, data, delimiter)


if __name__ == "__main__":
    main()
