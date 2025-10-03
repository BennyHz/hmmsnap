# setup.py
from setuptools import setup, find_packages
import os

setup(
    name="hmmsnap",
    version="0.1.0",
    description="Hidden Markov Model with signal smoothing for Nuclear Architecture-associated domains peak calling",
    long_description=open("README.md", encoding="utf-8").read() if os.path.exists("README.md") else "",
    long_description_content_type="text/markdown",
    author="BennyHz",
    author_email="haozhe_zhu@163.com",
    url="https://github.com/bennyhz/hmmsnap",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "hmmsnap=hmmsnap.cli:main",
        ],
    },
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.19",
        "pandas>=1.3",
        "pyBigWig>=0.3.18",
        "pysam>=0.18",
        "hmmlearn>=0.2.8",
        "pyranges>=0.0.120",
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
    ],
    keywords="chip-seq nad lad broad-domain hmm smoothing",
)