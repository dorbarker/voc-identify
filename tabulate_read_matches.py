import argparse
import pysam


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


def get_reads(bam_path: str, ref_path: str):

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

        results[read.query_name] = set()

        seq = read.get_forward_sequence()

        start = read.reference_start

        for position in mutations:

            read_offset = position - start

            try:
                has_mutation = seq[read_offset] == mutations[position]

            except IndexError:
                has_mutation = False

            if has_mutation:
                results[read.query_name].add(position)

    return results


def format_report(read_results):

    only_mutant_reads = filter(lambda x: read_results[x], read_results)

    sorted_all_reads = sorted(
        read_results, key=lambda x: len(read_results[x]), reverse=True
    )


def write_report():
    pass


if __name__ == "__main__":
    main()
