import itertools
import logging
import re
import string
import yaml
from pathlib import Path
from typing import List, Tuple

import pandas as pd
import pysam

from mmmvi.lib.types import VoCs, Reads, Mutations


def load_mutations(
    mutations_path: Path,
    reference_path: Path,
    voc_col: str,
    mut_col: str,
    delimiter: str,
    selected_vocs: List[str],
) -> VoCs:
    # Decides whether to load variant definitions from a tabular file or from a directory
    # containing Public Health England-formatted YAML files.
    #
    # If the path provided by the --mutations command line argument is a file,
    # load_tabular_mutations is attempted. If the path instead refers to a directory,
    # load_mutations_phe is attempted.

    if mutations_path.is_file():
        vocs = load_tabular_mutations(
            mutations_path, reference_path, voc_col, mut_col, delimiter, selected_vocs
        )

    elif mutations_path.is_dir():
        vocs = load_mutations_phe(mutations_path, reference_path, selected_vocs)

    return vocs


def load_reference(reference: Path) -> str:
    # Loads the FASTA-formatted reference genome
    #
    # The reference genome *must* be a single complete sequence
    # in FASTA format

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


def load_tabular_mutations(
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

        mutation_string = row[mut_col].strip()

        try:
            position_range, wt, mutant = parse_mutation(mutation_string)

        # catch *all* exceptions from parsing,
        # because any problems here should stop the program
        except Exception:
            msg = f"Invalid mutation string: '{mutation_string}'"
            raise InvalidMutation(msg)

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


def load_mutations_phe(
    mutations_dir: Path, reference_path: Path, selected_vocs: List[str]
) -> VoCs:
    # Manages loading variant definitions from a directory full of YAML files
    # using the schema described by https://github.com/phe-genomics/variant_definitions/

    vocs = {"reference": {}}

    # the spec explicitly states the extension will be .yml, and so we can rely on it
    variant_files = mutations_dir.glob("*.yml")

    for variant in variant_files:

        voc = variant.stem  # per the spec, the file name matches its 'unique-id' value
        reference, mutations = load_variant_from_phe_yaml(variant, reference_path)

        if selected_vocs and voc not in selected_vocs:
            continue

        vocs["reference"].update(reference)
        vocs[voc] = mutations

    return vocs


def load_variant_from_phe_yaml(
    yaml_variant: Path, reference_path: Path
) -> Tuple[Mutations, Mutations]:
    # Loads VOC signature mutations from a YAML file using
    # Public Health England's format for SARS-CoV-2 variants:
    # https://github.com/phe-genomics/variant_definitions/

    reference = {}
    voc = {}

    data = yaml.safe_load(yaml_variant.read_text())
    reference_seq = load_reference(reference_path)

    for mutation in data["variants"]:

        start = mutation["one-based-reference-position"] - 1

        if mutation["type"] == "SNP":

            wt = (mutation["reference-base"],)
            mutant = (mutation["variant-base"],)
            position_range = (start,)

        elif mutation["type"] == "MNP":

            wt = tuple(mutation["reference-base"])
            mutant = tuple(mutation["variant-base"])
            position_range = tuple(range(start, start + len(wt)))

        elif mutation["type"] == "insertion":

            mutant = tuple(mutation["variant-base"][1:])
            wt = tuple(None for _ in mutant)
            position_range = (start, None, start + 1)

        elif mutation["type"] == "deletion":

            position_range = tuple(
                range(start, start + len(mutation["reference-base"]) - 1)
            )
            wt = (None,)
            mutant = tuple(None for _ in position_range)

        else:
            msg = "Mutation type '{}' is not implemented".format(mutation["type"])
            raise NotImplementedError(msg)

        if wt == (None,):
            wt = tuple(reference_seq[position] for position in position_range)

        try:
            voc[position_range].add(mutant)
        except KeyError:
            voc[position_range] = {mutant}

        reference[position_range] = [wt]

    return reference, voc


def load_reads(bam_path: Path, ref_path: Path) -> Reads:
    # Loads reads from a BAM file on disk and returns the unique reads.
    #
    # The the sequence is used as the key. The dictionary keeps track of the
    # set of read names which share that sequence, as well as a representative
    # pysam.AlignedSegment object

    logging.info(f"Loading reads from {bam_path}")

    reads = {}

    with pysam.AlignmentFile(
        bam_path, reference_filename=str(ref_path), mode="rb"
    ) as bam:

        for read in bam:

            seq = read.query_sequence
            orientation_tag = "rev" if read.is_reverse else "fwd"
            read_name = f"{read.query_name}:{orientation_tag}"

            try:
                reads[seq]["reads"].add(read_name)

            except KeyError:

                reads[seq] = {"reads": {read_name}, "read_obj": read}
    return reads


class InvalidMutation(Exception):
    pass
