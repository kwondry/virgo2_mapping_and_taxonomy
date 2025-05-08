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
├── resources/          # Reference data and scripts
│   ├── VIRGO2_20250507/  # VIRGO2 database and scripts
│   ├── test_data/      # Test dataset
│   └── setup_virgo2_db.sh  # Database setup script
└── workflow/           # Snakemake workflow files
```

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

In the current version of this pipeline, the exact database structure is very picky, and you need to download the database from Dropbox. After publication, the files will be available from the Ravel lab on zenodo.

As of now, you need to:

1. get the dropbox link from Michael France  
2. download the folder and move the contents into the `resources/VIRGO2_20250507` folder
3. `gunzip` the `.txt.gz` files 
4. replace the python scripts with the ones from this repository
5. then do the following to setup the bowtie2 indicies:

Before running the workflow, you need to build the Bowtie2 index for the VIRGO2 database. A setup script is provided to help with this process:

1. Navigate to the resources directory:
   ```bash
   cd virgo2_mapping_and_taxonomy/resources
   ```

2. Make the setup script executable:
   ```bash
   chmod +x setup_virgo2_db.sh
   ```

3. Run the setup script:
   ```bash
   ./setup_virgo2_db.sh
   ```
Note: This step only needs to be performed once after downloading the workflow. The index files will be reused in subsequent runs.

If you want to store the VIRGO2 database in a different location (e.g., a shared reference directory), you can move the entire VIRGO2_20250507 folder and update `virgo2` section of the `config.yaml` to have the **absolute** paths to your database.


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

To run the workflow on a SLURM cluster:
```bash
snakemake --executor slurm --use-conda --configfile config/config.yaml --jobs 100
```

## Output Files

The workflow generates several output files:
- `*.summary.NR.txt`: Summary of mapping results
- `*_virgo2_NR_anno.csv`: Annotated results with taxonomic information
- `*_virgo2_metagenomic_taxa.csv`: Final output with taxonomic classifications
