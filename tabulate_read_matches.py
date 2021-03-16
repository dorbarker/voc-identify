# mutations.txt format
# 1234T 1422G

# notes
# samfile.fetch() will return all reads overlapping the requested region,
# including partial coverage
#
# Within pysam, coordinates are 0-based, half-open intervals
# BAM is also 0-based, so maybe should be converted to that before this script
# for least confusion


import pysam
import argparse
import itertools

def arguments():

    parser = argparse.ArgumentParser()

    parser = parser.add_argument('input')

    return parser.parse_args()


def main():
    pass


def load_mutations(mutations_path: str):

    groups = {}

    with open(mutations_path, 'r') as f:
        for idx, line in enumerate(f):

            l = line.strip()

            if l:
                groups[idx] = {}
                snps = l.split()

                for snp in snps:

                    position = int(snp[:-1])
                    nucleotide = snp[-1]

                groups[idx][position] = nucleotide

    return groups


def match_reads_to_region():
    pass


def is_mutant():
    pass


def has_all_mutations():
    pass


def format_report():
    pass

def write_report():
    pass


if __name__ == "__main__":
    main()
