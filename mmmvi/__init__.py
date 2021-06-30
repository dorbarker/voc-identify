import argparse

__version__ = "0.10.1"


class CitationAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        print(__citation__)
        parser.exit()


__citation__ = """
@article{Barker2021.06.14.448421,
	author = {Barker, Dillon O.R. and Buchanan, Cody J. and Landgraff, Chrystal and Taboada, Eduardo N},
	title = {MMMVI: Detecting SARS-CoV-2 Variants of Concern in Metagenomic Samples},
	elocation-id = {2021.06.14.448421},
	year = {2021},
	doi = {10.1101/2021.06.14.448421},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2021/06/15/2021.06.14.448421},
	eprint = {https://www.biorxiv.org/content/early/2021/06/15/2021.06.14.448421.full.pdf},
	journal = {bioRxiv}
}
"""
