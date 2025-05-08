#!/usr/bin/env python3
import pandas as pd
import numpy as np
import subprocess
import os

def covCorr(samFile, coverageFile):
    coverage = pd.read_csv(coverageFile, sep='\t', header=0, usecols=[0, 5], names=['Gene','PercentCovered'])
    coverage = dict(zip(coverage['Gene'], coverage['PercentCovered']))

    #open the .sam file as a pandas dataframe
    sam = pd.read_csv(samFile, sep='\t', header=None, usecols=[0,1], names=['Read','Gene'],quotechar='@')
    sam = sam.groupby('Read')['Gene'].apply(list).to_dict()

    #a dictionary that keeps track of best gene for each read
    best_gene_per_read = {}

    for Read, Genes in sam.items():
        if len(Genes) == 1 and Genes[0] == '*':
            continue
        best_gene_per_read[Read] = max(Genes, key=lambda x : coverage[x])

    samOut = pd.DataFrame(data={'Read':list(best_gene_per_read.keys()), 'Gene':list(best_gene_per_read.values())})

    return samOut

def noCovCorr(samFile):
    samOut = pd.read_csv(samFile, delimiter='\t', header=None, usecols=[0,1], names=['Read','Gene'], low_memory=False)
    samOut = samOut[samOut['Gene']!='*']
    return samOut

def geneCounts(samOut):
    samOut = samOut.drop(['Read'], axis=1)
    samOut['Count'] = 1
    samOut = samOut.groupby("Gene").sum()
    samOut = samOut.sort_values(by='Count', ascending=False)
    return samOut

def garbageCollect(file2del):
    try:
        os.remove(file2del)
    except OSError as e:
        print(f"Expected temporary file not available to delete: {file2del}")

def main():
    # Get database location from config
    database_dir = snakemake.config["virgo2"]["database_dir"]
    VIRGO2DBLoc = os.path.join(database_dir, 'VIRGO2')

    # Check if Bowtie2 index files exist
    required_index_files = [f"{VIRGO2DBLoc}.{i}.bt2" for i in range(1, 4)]
    missing_files = [f for f in required_index_files if not os.path.exists(f)]

    if missing_files:
        raise RuntimeError(f"Bowtie2 index files not found! Missing: {', '.join(missing_files)}")

    # Get parameters from snakemake
    reads = snakemake.input.reads
    output_prefix = snakemake.params.output_prefix
    threads = snakemake.threads
    use_coverage = snakemake.params.use_coverage

    if use_coverage:
        # Mapping with coverage correction
        try:
            subprocess.run(["bowtie2", "-N", "0", "-L", "20", "-D", "20", "-R", "3", "-k", "10",
                          "--sam-no-qname-trunc", "--end-to-end", "--seed", "343", "--no-unal",
                          "-p", str(threads), "-x", VIRGO2DBLoc, "-U", reads,
                          "-S", f"{output_prefix}.sam"], check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Bowtie2 mapping failed: {e}")

        # Sort SAM file
        try:
            subprocess.run(['samtools', 'sort', f'{output_prefix}.sam',
                          '-o', f'{output_prefix}.sorted.sam'], check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"SAM sorting failed: {e}")

        garbageCollect(f'{output_prefix}.sam')

        # Generate coverage values
        try:
            subprocess.run(['samtools', 'coverage', f"{output_prefix}.sorted.sam",
                          '-o', f"{output_prefix}.cov"], check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Coverage calculation failed: {e}")

        # Remove header from SAM file
        try:
            subprocess.run(f'samtools view -S {output_prefix}.sorted.sam | cut -f 1,3 > {output_prefix}.noHead.sam',
                         shell=True, check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"SAM header removal failed: {e}")

        garbageCollect(f'{output_prefix}.sorted.sam')

        # Process the mapping results
        samOut = covCorr(f"{output_prefix}.noHead.sam", f"{output_prefix}.cov")
        geneCounts = geneCounts(samOut)

        # Clean up temporary files
        garbageCollect(f'{output_prefix}.noHead.sam')
        garbageCollect(f'{output_prefix}.cov')

    else:
        # Mapping without coverage correction
        try:
            subprocess.run(["bowtie2", "-N", "0", "-L", "20", "-D", "20", "-R", "3",
                          "--sam-no-qname-trunc", "--end-to-end", "--seed", "343",
                          "--no-unal", "-p", str(threads), "-x", VIRGO2DBLoc,
                          "-U", reads, "-S", f"{output_prefix}.sam"], check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Bowtie2 mapping failed: {e}")

        # Sort SAM file
        try:
            subprocess.run(['samtools', 'sort', f'{output_prefix}.sam',
                          '-o', f'{output_prefix}.sorted.sam'], check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"SAM sorting failed: {e}")

        garbageCollect(f'{output_prefix}.sam')

        # Remove header from SAM file
        try:
            subprocess.run(['samtools', 'view', f"{output_prefix}.sorted.sam",
                          '-o', f"{output_prefix}.noHead.sam"], check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"SAM header removal failed: {e}")

        garbageCollect(f'{output_prefix}.sorted.sam')

        # Process the mapping results
        samOut = noCovCorr(f"{output_prefix}.noHead.sam")
        geneCounts = geneCounts(samOut)
        garbageCollect(f'{output_prefix}.noHead.sam')

    # Write output
    geneCounts.to_csv(snakemake.output[0], sep="\t")

if __name__ == "__main__":
    main() 