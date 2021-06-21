import itertools
import logging
import statistics
from collections import Counter
from pathlib import Path
from typing import Dict, Tuple

import pandas as pd

from mmmvi.lib.types import VoCs, VoCResults, Reads


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


def expand_dataframe_by_reads(df: pd.DataFrame, reads: Reads):
    def _expand_rows():
        for seq, row in df.iterrows():
            for read_name in reads[seq]["reads"]:
                yield row._set_name(read_name, inplace=False)

    return pd.DataFrame(_expand_rows())


def format_read_report(oir_results: VoCResults, reads: Reads) -> pd.DataFrame:

    read_report = pd.DataFrame(oir_results)

    has_any_results = read_report.applymap(len).apply(sum, axis="columns") > 0

    reads_with_results = read_report[has_any_results]

    expanded = expand_dataframe_by_reads(reads_with_results, reads)

    expanded["first_pos"] = expanded.apply(
        lambda row: sorted(itertools.chain.from_iterable(row))[0], axis=1
    )

    expanded_sorted = expanded.sort_values(by="first_pos", ascending=True).drop(
        "first_pos", axis=1
    )

    return expanded_sorted


def format_summary(voc_results: VoCResults, vocs: VoCs, reads: Reads) -> pd.DataFrame:

    mutation_df = expand_dataframe_by_reads(pd.DataFrame(voc_results), reads)

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


def mutation_coverage(voc_results, vocs) -> pd.DataFrame:
    # For each variant, determines how many its signature mutations
    # are represented in the sample.
    #
    # The exclusive signature mutations are those those mutations
    # which are not shared between with any other variant.

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


def theoretical_maximum(reads: Reads, vocs: VoCs) -> pd.DataFrame:
    # Finds the maximum number of mutations for a given variant
    # that will fit on the median length read.
    #
    # The result depends both on read length and on the particular
    # genomic positions of the mutations for each variant

    all_lengths = []
    for seq, read_data in reads.items():
        seq_length = len(seq)
        for _ in read_data["reads"]:
            all_lengths.append(seq_length)

    median_read_length = statistics.median(all_lengths)

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


def format_cooccurrence_matrix(mutation_result, voc, wt, reads: Reads):
    # For one VoC at a time

    lookup, mx = initialize_matrix(voc, wt)

    for seq, read_mutations in mutation_result.items():

        for position, mutation in read_mutations:

            name = lookup[position][mutation]

            mx.loc[name, name] += len(reads[seq]["reads"])

        for (row_pos, row_mut), (col_pos, col_mut) in itertools.permutations(
            read_mutations, r=2
        ):

            row_name = lookup[row_pos][row_mut]
            col_name = lookup[col_pos][col_mut]

            mx.loc[row_name, col_name] += len(reads[seq]["reads"])

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


def format_cooccurrence_matrices(voc_results: VoCResults, vocs: VoCs, reads: Reads):
    *variants, wt = sorted(vocs.keys(), key=lambda x: x == "reference")

    return {
        v: format_cooccurrence_matrix(voc_results[v], vocs[v], vocs[wt], reads)
        for v in variants
    }


def format_read_species(voc_results, vocs, reads):

    species = {}
    total_reads = len(voc_results["reference"].keys())

    # get non-redundant set of positions across VOCs

    for key, species_data in nonredundant_read_species(voc_results, reads):

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


def nonredundant_read_species(voc_results, reads: Reads):
    # old -> {Voc: {name: mutations}}
    # new -> {Voc: {seq: mutations}}

    nonredundant_reads = set()  # full of sequences
    for read_results in voc_results.values():
        for seq in read_results.keys():
            nonredundant_reads.add(seq)

    nonredundant = {}
    for seq in nonredundant_reads:
        positions_mutations = set()
        for voc in voc_results:
            positions_mutations.update(voc_results[voc][seq])
        positions_mutations = sorted(positions_mutations)

        if not positions_mutations:
            continue

        key = str(positions_mutations)
        read_count = len(reads[seq]["reads"])

        try:
            nonredundant[key]["count"] += read_count

        except KeyError:
            nonredundant[key] = {
                "positions_mutations": positions_mutations,
                "count": read_count,
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

    for seq, read_data in reads.items():

        read_start, *_, read_end = sorted(read_data["reference_positions"])

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
        "read_report": format_read_report(oir_results, reads),
        "summary": format_summary(voc_results, vocs, reads),
        "absolute_cooccurrence_matrices": format_cooccurrence_matrices(
            voc_results, vocs, reads
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
