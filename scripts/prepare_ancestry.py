#!/usr/bin/env python
"""
This script prepares the genetic ancestry data for unsupervised learning.
"""

import os
import pandas as pd
import numpy as np

# File paths
DATA_DIR = '/home/forde/PrecisionOncology/data/tcga_data/coadread_tcga_pan_can_atlas_2018'
OUTPUT_DIR = '/home/forde/PrecisionOncology/data/adapted_data'
ANCESTRY_FILE = os.path.join(DATA_DIR, 'data_genetic_ancestry.txt')
CLINICAL_SAMPLE_FILE = os.path.join(DATA_DIR, 'data_clinical_sample.txt')
OUTPUT_FILE = os.path.join(OUTPUT_DIR, 'ancestry_matrix.csv')

# Make sure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

def main():
    print("Loading clinical sample data...")
    # Load clinical sample data with correct separator and skipping comment rows
    clinical_sample_df = pd.read_csv(CLINICAL_SAMPLE_FILE, sep='\t', comment='#')
    
    # Create mapping from sample ID to patient ID
    sample_to_patient = dict(zip(clinical_sample_df['SAMPLE_ID'], clinical_sample_df['PATIENT_ID']))
    
    print("Loading genetic ancestry data...")
    # Load ancestry data
    ancestry_df = pd.read_csv(ANCESTRY_FILE, sep='\t')
    
    # Set index to ancestry types
    ancestry_df = ancestry_df.set_index(['ENTITY_STABLE_ID', 'NAME'])
    
    # Rename columns to match patient IDs
    new_columns = {}
    for col in ancestry_df.columns:
        if col in sample_to_patient:
            new_columns[col] = sample_to_patient[col]
    
    ancestry_df = ancestry_df.rename(columns=new_columns)
    
    # Drop the ENTITY_STABLE_ID level from the index to keep just NAME
    ancestry_df.index = ancestry_df.index.droplevel(0)
    
    # Save to CSV
    print(f"Saving ancestry matrix to {OUTPUT_FILE}...")
    ancestry_df.to_csv(OUTPUT_FILE)
    
    # Print some statistics
    print(f"Number of ancestry types: {ancestry_df.shape[0]}")
    print(f"Number of patients: {ancestry_df.shape[1]}")
    
    # Calculate statistics for each ancestry type
    print("\nStatistics for each ancestry type:")
    for ancestry_type in ancestry_df.index:
        data = ancestry_df.loc[ancestry_type]
        print(f"{ancestry_type}:")
        print(f"  Mean: {data.mean():.4f}")
        print(f"  Std: {data.std():.4f}")
        print(f"  Min: {data.min():.4f}")
        print(f"  Max: {data.max():.4f}")
    
    print("\nDone!")

if __name__ == "__main__":
    main()