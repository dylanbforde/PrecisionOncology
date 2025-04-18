#!/usr/bin/env python
"""
This script prepares the mutation data for unsupervised learning.
It creates a binary matrix where:
- Columns are the patient IDs
- Rows are gene symbols
- Values are 1 if the patient has a mutation in that gene, 0 otherwise
"""

import os
import pandas as pd
import numpy as np
from collections import defaultdict

# File paths
DATA_DIR = '/home/forde/PrecisionOncology/data/tcga_data/coadread_tcga_pan_can_atlas_2018'
OUTPUT_DIR = '/home/forde/PrecisionOncology/data/adapted_data'
MUTATIONS_FILE = os.path.join(DATA_DIR, 'data_mutations.txt')
CLINICAL_SAMPLE_FILE = os.path.join(DATA_DIR, 'data_clinical_sample.txt')
OUTPUT_FILE = os.path.join(OUTPUT_DIR, 'mutations_matrix.csv')

# Make sure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

def main():
    print("Loading clinical sample data...")
    # Load clinical sample data with correct separator and skipping comment rows
    clinical_sample_df = pd.read_csv(CLINICAL_SAMPLE_FILE, sep='\t', comment='#')
    
    # Create mapping from sample ID to patient ID
    sample_to_patient = dict(zip(clinical_sample_df['SAMPLE_ID'], clinical_sample_df['PATIENT_ID']))
    
    print("Loading mutations data...")
    # Load mutations data
    mutations_df = pd.read_csv(MUTATIONS_FILE, sep='\t')
    
    # Map sample IDs to patient IDs
    mutations_df['PATIENT_ID'] = mutations_df['Tumor_Sample_Barcode'].map(sample_to_patient)
    
    # Create a dictionary to store mutation data
    # Initialize with all patients having no mutations
    patient_ids = clinical_sample_df['PATIENT_ID'].unique()
    gene_mutation_dict = defaultdict(lambda: {patient_id: 0 for patient_id in patient_ids})
    
    print("Processing mutations...")
    # Populate the dictionary with mutation data
    for _, row in mutations_df.iterrows():
        if pd.notna(row['PATIENT_ID']):
            gene_symbol = row['Hugo_Symbol']
            patient_id = row['PATIENT_ID']
            gene_mutation_dict[gene_symbol][patient_id] = 1
    
    # Convert the dictionary to a DataFrame
    print("Creating mutations matrix...")
    mutations_matrix = pd.DataFrame.from_dict(gene_mutation_dict, orient='index')
    
    # Save to CSV
    print(f"Saving mutations matrix to {OUTPUT_FILE}...")
    mutations_matrix.to_csv(OUTPUT_FILE)
    
    # Print some statistics
    print(f"Number of genes: {mutations_matrix.shape[0]}")
    print(f"Number of patients: {mutations_matrix.shape[1]}")
    print(f"Number of mutations: {int(mutations_matrix.sum().sum())}")
    
    # Calculate mutation frequency per gene
    gene_mutation_counts = mutations_matrix.sum(axis=1)
    gene_mutation_freq = gene_mutation_counts / mutations_matrix.shape[1]
    
    # Find top 20 most frequently mutated genes
    top_genes = gene_mutation_freq.sort_values(ascending=False).head(20)
    print("\nTop 20 most frequently mutated genes:")
    for gene, freq in top_genes.items():
        print(f"{gene}: {freq:.4f} ({int(gene_mutation_counts[gene])} patients)")
    
    print("\nDone!")

if __name__ == "__main__":
    main()