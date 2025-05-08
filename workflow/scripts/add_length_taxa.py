#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os

def main():
    # Read input summary file
    data = pd.read_csv(snakemake.input.summary, sep="\t")
    data['Gene'] = data['Gene'].astype('str')

    # Read gene length data
    len_data = pd.read_csv(os.path.join(snakemake.input.database_dir, "0.VIRGO2.geneLength.txt"), 
                          sep="\t", header=None)
    len_data.columns = ['Gene', 'Length']
    len_data['Gene'] = len_data['Gene'].astype('str')

    # Read taxa annotation
    taxa_data = pd.read_csv(os.path.join(snakemake.input.database_dir, "1.VIRGO2.taxon.txt"), 
                           sep="\t")
    taxa_data = taxa_data.drop(['Cluster'], axis=1)
    taxa_data['Gene'] = taxa_data['Gene'].astype('str')

    # Read gene product annotation
    gp_data = pd.read_csv(os.path.join(snakemake.input.database_dir, "6.VIRGO2.geneProduct.txt"), 
                         sep="\t")
    gp_data['Gene'] = gp_data['Gene'].astype('str')

    # Merge all data
    data = pd.merge(left=len_data, right=data, left_on="Gene", right_on='Gene', how='right')
    data = pd.merge(left=taxa_data, right=data, left_on="Gene", right_on='Gene', how='right')
    data = pd.merge(left=gp_data, right=data, on="Gene", how='right')

    # Write output
    data.to_csv(snakemake.output[0], sep=",", index=None)

if __name__ == "__main__":
    main() 