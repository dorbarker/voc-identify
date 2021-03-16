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

def arguments():

    parser = argparse.ArgumentParser()

    parser = parser.add_argument('input')

    return parser.parse_args()


def main():
    pass


def load_mutations():
    pass


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
