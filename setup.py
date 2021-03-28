from setuptools import setup, find_packages
from src import __version__

setup(
    name="mmmvi",
    version=__version__,
    packages=find_packages(),
    author="Dillon O.R. Barker",
    author_email="dillon.barker@canada.ca",
    entry_points={"console_scripts": ["mmmvi=src.mmmvi:main",]},
)
