from setuptools import setup, find_packages

setup(
    name="mmmvi",
    version="0.1.1",
    packages=find_packages(),
    author="Dillon O.R. Barker",
    author_email="dillon.barker@canada.ca",
    entry_points={"console_scripts": ["mmmvi=src.mmmvi:main",]},
)
