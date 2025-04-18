#!/usr/bin/env python
"""
This script prepares the mRNA expression data for unsupervised learning.
"""

import os
import pandas as pd
import numpy as np

# File paths
DATA_DIR = '/home/forde/PrecisionOncology/data/tcga_data/coadread_tcga_pan_can_atlas_2018'
OUTPUT_DIR = '/home/forde/PrecisionOncology/data/adapted_data'
MRNA_FILE = os.path.join(DATA_DIR, 'data_mrna_seq_v2_rsem.txt')
CLINICAL_SAMPLE_FILE = os.path.join(DATA_DIR, 'data_clinical_sample.txt')
OUTPUT_FILE = os.path.join(OUTPUT_DIR, 'mrna_matrix.csv')

# Make sure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

def main():
    print("Loading clinical sample data...")
    # Load clinical sample data with correct separator and skipping comment rows
    clinical_sample_df = pd.read_csv(CLINICAL_SAMPLE_FILE, sep='\t', comment='#')
    
    # Create mapping from sample ID to patient ID
    sample_to_patient = dict(zip(clinical_sample_df['SAMPLE_ID'], clinical_sample_df['PATIENT_ID']))
    
    print("Loading mRNA data...")
    # Load mRNA data
    # For large files, we can use chunksize to process in batches
    chunk_size = 1000
    chunks = pd.read_csv(MRNA_FILE, sep='\t', chunksize=chunk_size)
    
    first_chunk = next(chunks)
    # Initialize the DataFrame with the first chunk
    mrna_df = first_chunk.copy()
    
    # Process the rest of the chunks
    for chunk in chunks:
        mrna_df = pd.concat([mrna_df, chunk])
    
    # Set index to gene symbols
    mrna_df = mrna_df.set_index(['Hugo_Symbol', 'Entrez_Gene_Id'])
    
    # Rename columns to match patient IDs
    new_columns = {}
    for col in mrna_df.columns:
        if col in sample_to_patient:
            new_columns[col] = sample_to_patient[col]
    
    mrna_df = mrna_df.rename(columns=new_columns)
    
    # Drop the Entrez_Gene_Id level from the index to keep just Hugo_Symbol
    mrna_df.index = mrna_df.index.droplevel(1)
    
    # Save to CSV
    print(f"Saving mRNA matrix to {OUTPUT_FILE}...")
    mrna_df.to_csv(OUTPUT_FILE)
    
    # Print some statistics
    print(f"Number of genes: {mrna_df.shape[0]}")
    print(f"Number of patients: {mrna_df.shape[1]}")
    
    # Calculate basic statistics of mRNA expression values
    mrna_mean = mrna_df.mean().mean()
    mrna_std = mrna_df.std().mean()
    mrna_min = mrna_df.min().min()
    mrna_max = mrna_df.max().max()
    
    print(f"\nBasic statistics of mRNA expression values:")
    print(f"Mean: {mrna_mean:.2f}")
    print(f"Standard deviation: {mrna_std:.2f}")
    print(f"Minimum: {mrna_min:.2f}")
    print(f"Maximum: {mrna_max:.2f}")
    
    print("\nDone!")

if __name__ == "__main__":
    main()