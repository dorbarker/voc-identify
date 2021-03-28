import argparse
import itertools
import pysam
from collections import Counter
import pandas as pd
from pathlib import Path
from typing import Dict, List

from . import __version__

complements = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

Mutations = Dict[int, str]
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

    vocs = load_mutations(args.mutations, args.delimiter)

    reads = load_reads(args.bam, args.reference)

    mutation_results = find_mutations(reads, vocs)

    reports = format_reports(mutation_results, vocs)

    write_reports(reports, args.outdir, args.delimiter)


def is_illumina(reads: Reads) -> bool:
    # Heuristically determine if the reads are paired or not.
    #
    # If duplicated read names outnumber singleton read names by
    # a factor of at least 10:1, then it's Illumina

    counts = Counter(Counter(read.query_name for read in reads).values())

    return (counts[2] / counts[1]) >= 10


def load_mutations(mutations_path: Path, delimiter: str) -> VoCs:

    data = pd.read_csv(mutations_path, sep=delimiter)

    vocs = {"reference": {}}

    for idx, row in data.iterrows():

        # Currently only single-base substitutions are supported
        if row["Type"] == "Del" or len(row["Alt"]) > 1:
            continue

        voc = row["PangoLineage"]
        position = int(row["Position"]) - 1
        mutant = row["Alt"]

        if voc not in vocs:
            vocs[voc] = {}

        vocs[voc][position] = mutant

        vocs["reference"][position] = row["Ref"]

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

        results[read_name] = [
            s for q, s in pairs if is_mutant(q, s, seq, False, mutations)
        ]

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

        results[read_name] = [
            s for q, s in pairs if is_mutant(q, s, seq, revcomp, mutations)
        ]

    return results


def is_mutant(
    read_position: int,
    reference_position: int,
    read_sequence: str,
    revcomp: bool,
    mutations: Mutations,
) -> bool:

    if reference_position in mutations:

        if revcomp:
            read_sequence = "".join([complements[nt] for nt in reversed(read_sequence)])

        try:

            read_nt = read_sequence[read_position]
            mut_nt = mutations[reference_position]

            return read_nt == mut_nt

        except TypeError:
            return False


def one_index_results(
    mutation_results: Dict[str, MutationResults]
) -> Dict[str, MutationResults]:

    oir = {
        voc_name: {
            read: [i + 1 for i in positions] for read, positions in voc_results.items()
        }
        for voc_name, voc_results in mutation_results.items()
    }

    return oir


def format_read_report(oir_results: Dict[str, MutationResults]) -> pd.DataFrame:

    read_report = pd.DataFrame(oir_results)

    has_any_results = read_report.applymap(len).apply(sum, axis="columns") > 0

    return read_report[has_any_results]


def format_summary(mutation_results):

    mutation_df = pd.DataFrame(mutation_results)

    count_of_reads_with_N_snps = mutation_df.applymap(len).agg(Counter)

    return pd.DataFrame(count_of_reads_with_N_snps.to_dict()).transpose()


def format_cooccurence_matrix(mutation_result, mutations, wt) -> pd.DataFrame:
    # For one VoC at a time

    lookup = {pos: f"{wt[pos]}{pos+1}{mutations[pos]}" for pos in mutations.keys()}

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
    mutation_results: MutationResults, read_report: pd.DataFrame, vocs: VoCs
) -> pd.DataFrame:

    species = {}

    total_reads = len(mutation_results["reference"].keys())

    for _, variant_positions in read_report.iterrows():

        position_nts = set()

        matching_variant = {v: 0 for v in vocs}

        for variant, positions in zip(read_report.columns, variant_positions):

            for position in positions:

                matching_variant[variant] += 1

                # convert between 1-based read_report and 0-based vocs
                voc_pos = position - 1

                position_nt = (position, vocs[variant][voc_pos])
                position_nts.add(position_nt)

        if not position_nts:  # reads matching no VoCs
            continue

        pairs = sorted(position_nts, key=lambda x: x[0])
        key = str(tuple(f"{p}{nt}" for p, nt in pairs))

        try:
            species[key]["count"] += 1

        except KeyError:

            locations, nucleotides = zip(*pairs)
            species[key] = {
                "positions": locations,
                "nucleotides": nucleotides,
                "count": 1,
            }
            species[key].update(matching_variant)

    read_species = pd.DataFrame.from_dict(species, orient="index")

    read_species["proportion"] = read_species["count"] / total_reads

    return read_species


def format_reports(mutation_results, vocs):

    oir_results = one_index_results(mutation_results)

    reports = {
        "read_report": format_read_report(oir_results),
        "summary": format_summary(mutation_results),
        "absolute_cooccurence_matrices": format_cooccurence_matrices(
            mutation_results, vocs
        ),
    }
    reports["read_species"] = format_read_species(
        mutation_results, reports["read_report"], vocs
    )

    reports["relative_cooccurence_matrices"] = format_relative_coocurence_matrices(
        reports["absolute_cooccurence_matrices"]
    )
    return reports


def write_cooccurence_matrix(
    variant: str, directory: Path, data: pd.DataFrame, delimiter: str
) -> None:

    p = directory.joinpath(f"{variant}.txt")
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
