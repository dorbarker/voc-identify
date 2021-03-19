import argparse
import pysam
from collections import Counter


def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument("--bam", required=True)

    parser.add_argument("--reference", required=True)

    parser.add_argument("--mutations", required=True)

    return parser.parse_args()


def main():

    args = arguments()

    vocs = load_mutations(args.mutations)

    reads = get_reads(args.bam, args.reference)

    mutant_reads = find_mutations(reads, vocs)

    print(format_report(mutant_reads))


def load_mutations(mutations_path: str):

    groups = {}

    with open(mutations_path, "r") as f:
        for idx, line in enumerate(f):

            l = line.strip()

            if l:
                groups[idx] = {}
                snps = l.split()

                for snp in snps:

                    # Convert 1-based reference to
                    # 0-based for BAM and pysam
                    position = int(snp[:-1]) - 1
                    nucleotide = snp[-1]

                groups[idx][position] = nucleotide

    return groups


def load_reads(bam_path: str, ref_path: str):

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

    count_of_reads_with_N_snps = mutation_results.applymap(len).agg(Counter)

    return pd.DataFrame(count_of_reads_with_N_snps.to_dict()).transpose()


def format_report(read_results):

    only_mutant_reads = filter(lambda x: read_results[x], read_results)

    sorted_all_reads = sorted(
        read_results, key=lambda x: len(read_results[x]), reverse=True
    )


def write_report():
    pass


if __name__ == "__main__":
    main()