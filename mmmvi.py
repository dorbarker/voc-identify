import argparse
import itertools
import pysam
from collections import Counter
import pandas as pd
from pathlib import Path


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
        help="Delimiter character for tabular inpu and output [TAB]",
    )

    return parser.parse_args()


def main():

    args = arguments()

    vocs = load_mutations(args.mutations, args.delimiter)

    reads = load_reads(args.bam, args.reference)

    mutation_results = find_mutations(reads, vocs)

    reports = format_reports(mutation_results, vocs)

    write_reports(reports, args.outdir, args.delimiter)


def load_mutations(mutations_path: Path, delimiter: str):

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


def load_reads(bam_path: Path, ref_path: Path):

    with pysam.AlignmentFile(bam_path, reference_filename=ref_path, mode="rb") as aln:
        return list(aln)


def find_mutations(reads, vocs):

    results = {}

    for variant, mutations in vocs.items():

        results[variant] = find_variant_mutations(reads, mutations)

    return results


def find_variant_mutations(reads, mutations):

    results = {}

    for read in reads:

        read_name = read.query_name

        seq = read.get_forward_sequence()

        pairs = read.get_aligned_pairs()

        results[read_name] = [s for q, s in pairs if is_mutant(q, s, seq, mutations)]

    return results


def is_mutant(read_position, reference_position, read_sequence, mutations):

    if reference_position in mutations:

        try:

            read_nt = read_sequence[read_position]
            mut_nt = mutations[reference_position]

            return read_nt == mut_nt

        except TypeError:
            return False


def one_index_results(mutation_results):

    oir = {
        voc_name: {
            read: [i + 1 for i in positions] for read, positions in voc_results.items()
        }
        for voc_name, voc_results in mutation_results.items()
    }

    return oir


def format_read_report(oir_results):

    return pd.DataFrame(oir_results)


def format_summary(mutation_results):

    mutation_df = pd.DataFrame(mutation_results)

    count_of_reads_with_N_snps = mutation_df.applymap(len).agg(Counter)

    return pd.DataFrame(count_of_reads_with_N_snps.to_dict()).transpose()


def format_cooccurence_matrix(mutation_result, mutations, wt):
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


def format_cooccurence_matrices(mutation_results, vocs):

    *variants, wt = sorted(vocs.keys(), key=lambda x: x == "reference")

    return {
        v: format_cooccurence_matrix(mutation_results[v], vocs[v], vocs[wt])
        for v in variants
    }


def format_reports(mutation_results, vocs):

    oir_results = one_index_results(mutation_results)

    reports = {
        "read_report": format_read_report(oir_results),
        "summary": format_summary(mutation_results),
        "cooccurence_matrices": format_cooccurence_matrices(mutation_results, vocs),
    }

    return reports


def write_reports(reports, outdir: Path, delimiter):

    outdir.joinpath("cooccurence_matrices").mkdir(parents=True, exist_ok=True)

    reports["read_report"].to_csv(outdir / "read_report.txt", sep=delimiter)

    reports["summary"].to_csv(outdir / "summary.txt", sep=delimiter)

    for variant, data in reports["cooccurence_matrices"].items():

        p = outdir.joinpath("cooccurence_matrices", f"{variant}.txt")

        data.to_csv(p, sep=delimiter)


if __name__ == "__main__":
    main()
