# Neurodevelopmental Disorders (NDD) Gene Analysis Pipeline

![Python](https://img.shields.io/badge/python-3.8%2B-blue)
![Conda](https://img.shields.io/badge/conda-supported-brightgreen)
![License](https://img.shields.io/badge/license-MIT-green)

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

3. Run the whole pipeline
   ```bash
   python run_ndd_analysis.py
   ```
