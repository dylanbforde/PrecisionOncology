#!/usr/bin/env python
"""
This script merges the prepared data files and prepares the final dataset for unsupervised learning.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import pathlib

# Get the project root directory (2 levels up from this script)
PROJECT_ROOT = pathlib.Path(__file__).parent.parent.parent.absolute()

# File paths
DATA_DIR = os.path.join(PROJECT_ROOT, 'data/adapted_data')
OUTPUT_DIR = os.path.join(PROJECT_ROOT, 'data/adapted_data')
MUTATIONS_FILE = os.path.join(DATA_DIR, 'mutations_matrix.csv')
CNA_FILE = os.path.join(DATA_DIR, 'cna_matrix.csv')
MRNA_FILE = os.path.join(DATA_DIR, 'mrna_matrix.csv')
ANCESTRY_FILE = os.path.join(DATA_DIR, 'ancestry_matrix.csv')
OUTPUT_FILE = os.path.join(OUTPUT_DIR, 'merged_data.csv')

# Make sure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

def main():
    print("Loading data...")
    
    # Load prepared data
    try:
        # For mutations, load and filter by frequency
        print("Loading and filtering mutations data...")
        mutations_df = pd.read_csv(MUTATIONS_FILE, index_col=0)
        
        # Calculate mutation frequency for each gene
        mutation_freq = mutations_df.mean(axis=1)
        
        # Filter genes based on frequency (between 5% and 95%)
        filtered_genes = mutation_freq[(mutation_freq >= 0.05) & (mutation_freq <= 0.95)].index
        print(f"Filtered mutations from {len(mutation_freq)} to {len(filtered_genes)} genes")
        
        # Keep only the filtered genes
        mutations_df = mutations_df.loc[filtered_genes]
        
        # Transpose to have patients as rows and genes as columns
        mutations_df = mutations_df.T
        
        # Add prefix to column names
        mutations_df.columns = [f"mutation_{col}" for col in mutations_df.columns]
        
        # Load CNA data
        print("Loading and processing CNA data...")
        cna_df = pd.read_csv(CNA_FILE, index_col=0)
        
        # Replace NaN with 0
        cna_df = cna_df.fillna(0)
        
        # Filter out genes with low variance
        cna_var = cna_df.var(axis=1)
        cna_df = cna_df.loc[cna_var.sort_values(ascending=False).index[:1000]]
        print(f"Filtered CNA to top 1000 most variable genes")
        
        # Transpose
        cna_df = cna_df.T
        
        # Add prefix to column names
        cna_df.columns = [f"cna_{col}" for col in cna_df.columns]
        
        # Load mRNA data
        print("Loading and processing mRNA data...")
        mrna_df = pd.read_csv(MRNA_FILE, index_col=0)
        
        # Filter out genes with low variance
        mrna_var = mrna_df.var(axis=1)
        mrna_df = mrna_df.loc[mrna_var.sort_values(ascending=False).index[:1000]]
        print(f"Filtered mRNA to top 1000 most variable genes")
        
        # Transpose
        mrna_df = mrna_df.T
        
        # Add prefix to column names
        mrna_df.columns = [f"mrna_{col}" for col in mrna_df.columns]
        
        # Load ancestry data
        print("Loading ancestry data...")
        ancestry_df = pd.read_csv(ANCESTRY_FILE, index_col=0)
        
        # Transpose
        ancestry_df = ancestry_df.T
        
        # Add prefix to column names
        ancestry_df.columns = [f"ancestry_{col}" for col in ancestry_df.columns]
        
        # Merge all datasets
        print("Merging datasets...")
        # Find common patients across all datasets
        common_patients = list(set(mutations_df.index) & set(cna_df.index) & set(mrna_df.index) & set(ancestry_df.index))
        print(f"Number of common patients across all datasets: {len(common_patients)}")
        
        # Filter each dataset to include only common patients
        mutations_df = mutations_df.loc[common_patients]
        cna_df = cna_df.loc[common_patients]
        mrna_df = mrna_df.loc[common_patients]
        ancestry_df = ancestry_df.loc[common_patients]
        
        # Merge the dataframes
        merged_df = pd.concat([mutations_df, cna_df, mrna_df, ancestry_df], axis=1)
        
        # Save to CSV
        print(f"Saving merged data to {OUTPUT_FILE}...")
        merged_df.to_csv(OUTPUT_FILE)
        
        # Print some statistics
        print(f"\nMerged data statistics:")
        print(f"Number of patients: {merged_df.shape[0]}")
        print(f"Number of features: {merged_df.shape[1]}")
        print(f"Number of mutation features: {mutations_df.shape[1]}")
        print(f"Number of CNA features: {cna_df.shape[1]}")
        print(f"Number of mRNA features: {mrna_df.shape[1]}")
        print(f"Number of ancestry features: {ancestry_df.shape[1]}")
        
        print("\nDone!")
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()