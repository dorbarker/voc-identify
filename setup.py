from setuptools import setup, find_packages
from mmmvi import __version__

setup(
    name="mmmvi",
    version=__version__,
    packages=find_packages(),
    install_requires=["pandas", "pysam"],
    author="Dillon O.R. Barker",
    author_email="dillon.barker@canada.ca",
    entry_points={"console_scripts": ["mmmvi=mmmvi:main"]},
)
