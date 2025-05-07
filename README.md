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

- Snakemake (version 8.20.0 or later)
- Python 3.10 or later
- VIRGO2 database (included in resources/VIRGO2_20250507)

## Configuration

The workflow must be configured through `config/config.yaml`. 
- Input/output paths
- VIRGO2 script locations
- Resource allocations (threads, memory, runtime)


## Database Setup

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

The script will:
- Check if Bowtie2 is installed
- Verify the VIRGO2 database file exists
- Build the Bowtie2 index using up to 12 threads
- Create the index files in the VIRGO2_20250507 directory

Note: This step only needs to be performed once after downloading the workflow. The index files will be reused in subsequent runs.

If you want to store the VIRGO2 database in a different location (e.g., a shared reference directory), you can move the entire VIRGO2_20250507 folder and update `virgo2` section of the `config.yaml` to have the **absolute** paths to your database.

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

The test configuration is already set up in `config/config.yaml` with appropriate resource allocations for testing.

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
