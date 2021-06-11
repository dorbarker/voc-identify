import logging

import pandas as pd

from typing import List, Optional, Tuple
from mmmvi.lib.types import (
    Reads,
    VoCs,
    VoCResults,
    Mutations,
    MutationResults,
    Position,
)


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
