#!/usr/bin/env python
"""
This script prepares the copy number alteration (CNA) data for unsupervised learning.
"""

import os
import pandas as pd
import numpy as np
import pathlib

# Get the project root directory (2 levels up from this script)
PROJECT_ROOT = pathlib.Path(__file__).parent.parent.parent.absolute()

# File paths
DATA_DIR = os.path.join(PROJECT_ROOT, 'data/tcga_data/coadread_tcga_pan_can_atlas_2018')
OUTPUT_DIR = os.path.join(PROJECT_ROOT, 'data/adapted_data')
CNA_FILE = os.path.join(DATA_DIR, 'data_cna.txt')
CLINICAL_SAMPLE_FILE = os.path.join(DATA_DIR, 'data_clinical_sample.txt')
OUTPUT_FILE = os.path.join(OUTPUT_DIR, 'cna_matrix.csv')

# Make sure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

def main():
    print("Loading clinical sample data...")
    # Load clinical sample data with correct separator and skipping comment rows
    clinical_sample_df = pd.read_csv(CLINICAL_SAMPLE_FILE, sep='\t', comment='#')
    
    # Create mapping from sample ID to patient ID
    sample_to_patient = dict(zip(clinical_sample_df['SAMPLE_ID'], clinical_sample_df['PATIENT_ID']))
    
    print("Loading CNA data...")
    # Load CNA data
    cna_df = pd.read_csv(CNA_FILE, sep='\t')
    
    # Set index to gene symbols
    cna_df = cna_df.set_index(['Hugo_Symbol', 'Entrez_Gene_Id'])
    
    # Rename columns to match patient IDs
    new_columns = {}
    for col in cna_df.columns:
        if col in sample_to_patient:
            new_columns[col] = sample_to_patient[col]
    
    cna_df = cna_df.rename(columns=new_columns)
    
    # Drop the Entrez_Gene_Id level from the index to keep just Hugo_Symbol
    cna_df.index = cna_df.index.droplevel(1)
    
    # Save to CSV
    print(f"Saving CNA matrix to {OUTPUT_FILE}...")
    cna_df.to_csv(OUTPUT_FILE)
    
    # Print some statistics
    print(f"Number of genes: {cna_df.shape[0]}")
    print(f"Number of patients: {cna_df.shape[1]}")
    
    # Calculate distribution of CNA values
    cna_values = cna_df.values.flatten()
    unique_values, counts = np.unique(cna_values, return_counts=True)
    
    print("\nDistribution of CNA values:")
    for val, count in zip(unique_values, counts):
        print(f"Value {val}: {count} occurrences ({count/len(cna_values)*100:.2f}%)")
    
    print("\nDone!")

if __name__ == "__main__":
    main()