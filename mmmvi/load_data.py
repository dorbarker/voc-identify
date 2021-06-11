import itertools
import logging
import re
import string
from pathlib import Path
from typing import List

import pandas as pd
import pysam

from .types import VoCs, Reads


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
    # Parses the mutation string from the mutations file
    #
    # The mutation string can be in one of 4 formats:
    # A123T         # point substitution
    # CAT123GTA     # multiple base substitution
    # [123-125]del  # deletion
    # 123CAT        # insertion

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
    # Loads the mutations file
    #
    # The mutations file is a tabular delimited file, which must
    # contain at least the following two columns, with other
    # columns ignored:
    #   1) a column containing the names of each variant (voc_col)
    #   2) a column containing the mutation strings (mut_col)

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
    # Loads reads from a BAM file on disk and returns them in a list
    #
    # The default pysam.AlignmentFile is an iterator which is consumed by use,
    # and storing the reads comes at a modest memory cost

    logging.info(f"Loading reads from {bam_path}")

    with pysam.AlignmentFile(
        bam_path, reference_filename=str(ref_path), mode="rb"
    ) as aln:
        return list(aln)
