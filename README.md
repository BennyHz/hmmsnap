# HMMsnap

**Hidden Markov Model with signal smoothing for Nuclear Architecture-associated domains peak calling**

HMMsnap is a command-line tool for calling broad chromatin domains (e.g.,NADs, LADs, Histone modifications) from NAD-seq, DamID-seq, ChIP-seq or other genomic signal data using a multi-state Gaussian HMM with signal smoothing.

## Key Features
- Input from BAM, bigWig, or precomputed ratio files
- Automatic genome normalization (hg38/hg19/mm10/mm39)
- Signal smoothing (median or convolution) + optional rebinning
- 2-5 state HMM with adaptive initialization
- Gap merging and minimum length filtering

## Installation
```bash
pip install git+https://github.com/yourname/hmmsnap.git
