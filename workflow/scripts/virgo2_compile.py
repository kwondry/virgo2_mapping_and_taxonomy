#!/usr/bin/env python3
import pandas as pd
import os

def main():
    # Collect all input files and combine them
    output2concat = []
    for input_file in snakemake.input:
        try:
            # Get sample name from the input file path
            sample_name = os.path.basename(input_file).split('.')[0]
            mappingDF = pd.read_csv(input_file, sep='\t', index_col=0)
            mappingDF.columns = [sample_name]
            output2concat.append(mappingDF)
        except Exception as e:
            raise RuntimeError(f"Could not process file {input_file}: {e}")
    
    if not output2concat:
        raise RuntimeError('No VIRGO2 mapping output files to compile')
    
    # Combine all dataframes
    compiledOutput = pd.concat(output2concat, axis=1)
    compiledOutput = compiledOutput.fillna(0)
    
    # Write output
    compiledOutput.to_csv(snakemake.output[0], sep="\t")

if __name__ == "__main__":
    main() 