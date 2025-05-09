#!/usr/bin/env python3

import os
import sys
import logging
import subprocess
import glob
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def gunzip_txt_files(directory):
    """
    Find and gunzip all .txt.gz files in the given directory.
    
    Args:
        directory (str): Directory to search for .txt.gz files
    """
    gz_files = glob.glob(os.path.join(directory, "*.txt.gz"))
    if gz_files:
        logger.info(f"Found {len(gz_files)} .txt.gz files to decompress")
        for gz_file in gz_files:
            logger.info(f"Decompressing {gz_file}")
            try:
                subprocess.run(["gunzip", "-f", gz_file], check=True)
                logger.info(f"Successfully decompressed {gz_file}")
            except subprocess.CalledProcessError as e:
                logger.error(f"Error decompressing {gz_file}: {e}")
                raise
    else:
        logger.info("No .txt.gz files found to decompress")

def run_bowtie2_build(database_dir):
    """
    Run bowtie2-build on the VIRGO2.fa.gz file if the index doesn't already exist.
    
    Args:
        database_dir (str): Directory containing the VIRGO2.fa.gz file
    """
    input_file = os.path.join(database_dir, "VIRGO2.fa.gz")
    output_prefix = os.path.join(database_dir, "VIRGO2")
    
    # Check for existing index files
    required_index_files = [
        f"{output_prefix}.1.bt2",
        f"{output_prefix}.2.bt2",
        f"{output_prefix}.3.bt2",
        f"{output_prefix}.4.bt2",
        f"{output_prefix}.rev.1.bt2",
        f"{output_prefix}.rev.2.bt2"
    ]
    
    # Check if all index files exist
    if all(os.path.exists(f) for f in required_index_files):
        logger.info("Bowtie2 index files already exist, skipping build")
        return
    
    logger.info("Building bowtie2 index...")
    try:
        subprocess.run(
            ["bowtie2-build", input_file, output_prefix],
            check=True,
            capture_output=True,
            text=True
        )
        logger.info("Successfully built bowtie2 index")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error building bowtie2 index: {e.stderr}")
        raise

def validate_database_location(config):
    """
    Validate that the database directory exists and contains necessary files.
    
    Args:
        config (dict): Configuration dictionary containing Virgo2 paths
        
    Returns:
        str: Path to the database directory
        
    Raises:
        FileNotFoundError: If required files or directories are missing
    """
    database_dir = config["virgo2"]["database_dir"]
    
    # Check if database directory exists
    if not os.path.isdir(database_dir):
        raise FileNotFoundError(
            f"Virgo2 database directory not found at: {database_dir}\n"
            "Please ensure the path in config.yaml is correct and the directory exists."
        )
    
    # Check for required database files
    required_files = [
        "VIRGO2.fa.gz",  # Main database file
    ]
    
    missing_files = []
    for file in required_files:
        file_path = os.path.join(database_dir, file)
        if not os.path.isfile(file_path):
            missing_files.append(file)
    
    if missing_files:
        raise FileNotFoundError(
            f"Missing required database files in {database_dir}:\n"
            f"- {'\n- '.join(missing_files)}\n"
            "Please ensure all required database files are present."
        )
    
    return database_dir

def setup_virgo2_db(config):
    """
    Validate the Virgo2 database location and required files,
    then run bowtie2-build and decompress any .txt.gz files.
    
    Args:
        config (dict): Configuration dictionary containing Virgo2 paths
    """
    try:
        # Validate database location and get path
        database_dir = validate_database_location(config)
        logger.info(f"Validated Virgo2 database in: {database_dir}")
        
        # Decompress any .txt.gz files
        gunzip_txt_files(database_dir)
        
        # Run bowtie2-build
        run_bowtie2_build(database_dir)
        
    except FileNotFoundError as e:
        logger.error(str(e))
        raise
    except subprocess.CalledProcessError as e:
        logger.error(f"Error during setup: {e}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error during Virgo2 setup: {e}")
        raise

if __name__ == "__main__":
    # Read config from snakemake
    if 'snakemake' in globals():
        config = snakemake.config
    else:
        # For testing outside of snakemake
        import yaml
        with open("config/config.yaml", 'r') as f:
            config = yaml.safe_load(f)
    
    setup_virgo2_db(config) 