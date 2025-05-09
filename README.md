# VIRGO2 Mapping and Taxonomy Workflow

This workflow is used to map reads to VIRGO2 and annotate the reads with taxonomic and functional information. It processes paired-end sequencing data and generates comprehensive taxonomic and functional annotations using the VIRGO2 database.

VIRGO2 is a resource developed by the Ravel lab at University of Maryland School of Medicine. The manuscript is available here: [doi: 10.1101/2025.03.04.641479](https://doi.org/10.1101/2025.03.04.641479)

## Overview

This workflow assumes you have already quality filtered, adapter trimmed, and host filtered your reads.

The workflow performs the following steps:
1. Maps sequencing reads to the VIRGO2 database
2. Generates taxonomic annotations
3. Calculates relative abundances
4. Produces summary reports and annotated output files

## Directory Structure

```
virgo2_mapping_and_taxonomy/
├── config/             # Configuration files
│   └── config.yaml     # Main configuration file
└── workflow/           # Snakemake workflow files
    ├── rules/         # Modular rule files
    ├── envs/          # Conda environment files
    ├── scripts/       # Analysis scripts
    ├── resources/     # Reference data and test files
    │   └── test_data/ # Test dataset
    └── Snakefile      # Main workflow definition
```

Note: The VIRGO2 database is not included in this repository and should be downloaded separately as described in the Database Setup section.

## Prerequisites

- conda
- Snakemake (version 8.20.0 or later)   
To install snakemake and conda, follow the instructions [on the kwondry website](https://kwondry.github.io/documentation/materials/getting-started/installation/#mambaminimamba)

### Installation

To install the workflow, run:

```bash
curl -L https://github.com/kwondry/virgo2_mapping_and_taxonomy/archive/refs/heads/main.zip -o main.zip
unzip main.zip && rm main.zip
```

### Database Setup

The VIRGO2 database is currently available via Dropbox. After publication, the files will be available from the Ravel lab on Zenodo.

To set up the database:

1. Get the Dropbox link from Michael France
2. Download and extract the database files to your preferred location
3. Update the `virgo2` section in `config/config.yaml` with the **absolute** path to your database location

Note: The database location can be anywhere on your system - it does not need to be within this workflow directory.

## Configuration

The workflow must be configured through `config/config.yaml`. 
- Input/output paths
- VIRGO2 script locations
- Resource allocations (threads, memory, runtime)


## Running the Workflow

### Test Data

A test dataset is provided in `resources/test_data/` containing three sample pairs. To run the workflow with test data:

1. Ensure you're in the workflow directory:
   ```bash
   cd virgo2_mapping_and_taxonomy
   ```

2. Run the workflow with test data:
   ```bash
   snakemake --use-conda --configfile config/config.yaml
   ```

The test configuration is already set up in `config/config.yaml`

### Running with your own data

To run the workflow with your own data:

1. Prepare a samplesheet in CSV format with the following columns:
   - sample: Sample identifier
   - fastq_1: Path to first read file
   - fastq_2: Path to second read file

2. Update the configuration in `config/config.yaml`:
   - Set the path to your samplesheet
   - Adjust resource requirements as needed

3. Run the workflow:
   ```bash
   snakemake --use-conda --configfile config/config.yaml
   ```

### Running on a Cluster

To submit the workflow to the cluster, use the provided submission script:

```bash
sbatch ./submit_jobs.sh
```

This is currently configured for the kwon lab on the O2 cluster, you may need to edit this script for your system.

## Output Files

The workflow generates several output files:
- `*.summary.NR.txt`: Summary of mapping results at the gene level
- `*_virgo2_NR_anno.csv`: results with the gene lengths and annotations added
- `*_virgo2_metagenomic_taxa.csv`: taxonomic relative abundances calculated from the gene counts corrected for gene length
