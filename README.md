# Neurodevelopmental Disorders (NDD) Gene Analysis Pipeline

![Python](https://img.shields.io/badge/python-3.8%2B-blue)
![Conda](https://img.shields.io/badge/conda-supported-brightgreen)
![License: code](https://img.shields.io/badge/code-MIT-green)
![License: data](https://img.shields.io/badge/data-CC%20BY%204.0-lightgrey)

A reproducible pipeline for identifying candidate genes in neurodevelopmental disorders through multi-database enrichment analysis.

## Table of Contents
- [Quick Start](#quick-start)
- [Detailed Setup](#detailed-setup)
- [Environment Management](#environment-management)
- [Running the Pipeline](#running-the-pipeline)
- [Output Interpretation](#output-interpretation)
- [Troubleshooting](#troubleshooting)


## Quick Start

1. Clone repository:
   ```bash
   git clone https://gitlab.com/biofold/ndd-analysis.git
   cd ndd-analysis
   ```

2. Generate Conda environment 
   ```bash
   conda config --set channel_priority flexible
   conda env create -f environment_[arch].yml
   conda activate ndd_analysis
   ```

3. Run the enrichment pipeline
   ```bash
   python run_ndd_analysis.py -c config.yml
   ```

## License

- **Code** (scripts, server, clients): MIT License, see [`LICENSE`](LICENSE).
- **Data and documentation** (gene classification, MOE scores, derived tables, mappings,
  documentation text): Creative Commons Attribution 4.0 International (CC BY 4.0), see
  [`LICENSE-DATA.md`](LICENSE-DATA.md). Third-party values remain under their sources' terms.

Please cite: Rivi C, Turina P, Capriotti E. *Identifying top candidate genes associated with
neurodevelopmental disorders*. Int. J. Mol. Sci. (manuscript ijms-4575333, in revision).
