#!/usr/bin/env python3
import pandas as pd
import numpy as np

def main():
    # Read input data
    data = pd.read_csv(snakemake.input[0], sep=",")
    
    # Divide each gene's read counts by its length to correct
    glcorr_data = data[data.columns[4:]].div(data['Length'], axis=0)
    
    # Combine with metadata columns
    data = pd.concat([data[data.columns[0:4]], glcorr_data], axis=1)
    
    # Filter out MultiGenera and NaN categories
    data = data[data['Cat'] != 'MultiGenera']
    data = data[data['Cat'].notna()]
    
    # Drop metadata columns
    data = data.drop(['Cat', 'Length', 'Gene', 'GeneProduct'], axis=1)
    
    # Group by taxa and sum
    data_taxa = data.groupby('Taxa').sum()
    
    # Calculate relative abundances
    data_taxa_rel = data_taxa.div(data_taxa.sum(axis=0), axis=1)
    
    # Sort by total abundance
    data_taxa_rel['taxa_abundance'] = data_taxa_rel.sum(axis=1)
    data_taxa_rel = data_taxa_rel.sort_values(by=["taxa_abundance"], ascending=False)
    data_taxa_rel = data_taxa_rel.drop(['taxa_abundance'], axis=1)
    
    # Write output
    data_taxa_rel.to_csv(snakemake.output[0], sep=",")

if __name__ == "__main__":
    main() 