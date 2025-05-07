#!/bin/bash

# Exit on error
set -e

# Check if bowtie2 is installed
if ! command -v bowtie2 &> /dev/null; then
    echo "Error: bowtie2 is not installed. Please activate an environment with bowtie2 installed."
    exit 1
fi

# Check if bowtie2-build is installed
if ! command -v bowtie2-build &> /dev/null; then
    echo "Error: bowtie2-build is not installed. Please activate an environment with bowtie2-build installed."
    exit 1
fi

# Set paths
VIRGO2_DIR="VIRGO2_20250507"
VIRGO2_FA="${VIRGO2_DIR}/VIRGO2.fa.gz"
OUTPUT_PREFIX="${VIRGO2_DIR}/VIRGO2"

# Check if input file exists
if [ ! -f "${VIRGO2_FA}" ]; then
    echo "Error: ${VIRGO2_FA} not found!"
    exit 1
fi

# Get number of available CPU cores
NUM_CORES=$(nproc)
if [ $NUM_CORES -gt 12 ]; then
    NUM_CORES=12
fi

echo "Building Bowtie2 index for VIRGO2 database..."
echo "Using ${NUM_CORES} threads"

# Build the index
bowtie2-build --threads ${NUM_CORES} -f "${VIRGO2_FA}" "${OUTPUT_PREFIX}"

echo "Done! Bowtie2 index has been created at ${OUTPUT_PREFIX}.*" 