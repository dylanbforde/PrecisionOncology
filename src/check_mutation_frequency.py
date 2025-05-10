#!/usr/bin/env python
"""
Script to check mutation frequencies in the dataset.
"""

import pandas as pd
import os
import pathlib

# Get the project root directory
PROJECT_ROOT = pathlib.Path(__file__).parent.parent.absolute()

# File paths
DATA_DIR = os.path.join(PROJECT_ROOT, 'data/adapted_data')
MERGED_DATA_FILE = os.path.join(DATA_DIR, 'merged_data.csv')

# Load the merged data
df = pd.read_csv(MERGED_DATA_FILE, index_col=0)

# Get mutation columns
mutation_cols = [col for col in df.columns if col.startswith('mutation_')]

# Print overall frequencies for key genes
key_genes = ['mutation_APC', 'mutation_TP53', 'mutation_KRAS', 'mutation_TTN', 'mutation_PIK3CA', 
             'mutation_FBXW7', 'mutation_SMAD4', 'mutation_BRAF']

print("Overall mutation frequencies:")
for gene in key_genes:
    if gene in df.columns:
        count = (df[gene] == 1).sum()
        total = len(df)
        print(f"{gene.replace('mutation_', '')} mutation frequency: {count} out of {total} samples ({count/total*100:.2f}%)")
    else:
        print(f"{gene.replace('mutation_', '')} not found in the dataset")

print("\nTop 10 most frequently mutated genes:")
mutation_freqs = {}
for gene in mutation_cols:
    count = (df[gene] == 1).sum()
    mutation_freqs[gene] = count / len(df) * 100

for gene, freq in sorted(mutation_freqs.items(), key=lambda x: x[1], reverse=True)[:10]:
    print(f"{gene.replace('mutation_', '')}: {freq:.2f}%")

# Print mutation frequencies by cluster
print("\nLoading cluster assignments...")
cluster_df = pd.read_csv(os.path.join(PROJECT_ROOT, 'results/cluster_assignments.csv'), index_col=0)

print("\nMutation frequencies by cluster:")
for cluster_id in sorted(cluster_df['K-means Cluster'].unique()):
    print(f"\nCluster {cluster_id}:")
    cluster_samples = cluster_df[cluster_df['K-means Cluster'] == cluster_id].index
    cluster_mutations = df.loc[cluster_samples, mutation_cols]
    
    for gene in key_genes:
        if gene in df.columns:
            count = (cluster_mutations[gene] == 1).sum()
            total = len(cluster_samples)
            print(f"{gene.replace('mutation_', '')} mutation frequency: {count} out of {total} samples ({count/total*100:.2f}%)")