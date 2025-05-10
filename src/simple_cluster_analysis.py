#!/usr/bin/env python
"""
Simplified script to analyze clusters and display key differences.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import silhouette_score
from scipy import stats
import pathlib

# Get the project root directory
PROJECT_ROOT = pathlib.Path(__file__).parent.parent.absolute()

# File paths
DATA_DIR = os.path.join(PROJECT_ROOT, 'data/adapted_data')
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results')
MERGED_DATA_FILE = os.path.join(DATA_DIR, 'merged_data.csv')
CLUSTER_ASSIGNMENTS_FILE = os.path.join(RESULTS_DIR, 'cluster_assignments.csv')
TSNE_RESULTS_FILE = os.path.join(RESULTS_DIR, 'tsne_results.csv')

def main():
    # Load data
    print("Loading data...")
    merged_df = pd.read_csv(MERGED_DATA_FILE, index_col=0)
    merged_df = merged_df.fillna(0)
    
    cluster_df = pd.read_csv(CLUSTER_ASSIGNMENTS_FILE, index_col=0)
    tsne_df = pd.read_csv(TSNE_RESULTS_FILE, index_col=0)
    
    print(f"Merged data shape: {merged_df.shape}")
    print(f"Cluster assignments shape: {cluster_df.shape}")
    print(f"t-SNE results shape: {tsne_df.shape}")
    
    # Check column names
    print("\nCluster assignment columns:")
    print(cluster_df.columns.tolist())
    
    # Split data by type
    mutation_cols = [col for col in merged_df.columns if col.startswith('mutation_')]
    cna_cols = [col for col in merged_df.columns if col.startswith('cna_')]
    mrna_cols = [col for col in merged_df.columns if col.startswith('mrna_')]
    ancestry_cols = [col for col in merged_df.columns if col.startswith('ancestry_')]
    
    mutations_df = merged_df[mutation_cols]
    cna_df = merged_df[cna_cols]
    mrna_df = merged_df[mrna_cols]
    ancestry_df = merged_df[ancestry_cols]
    
    print(f"\nNumber of mutation features: {len(mutation_cols)}")
    print(f"Number of CNA features: {len(cna_cols)}")
    print(f"Number of mRNA features: {len(mrna_cols)}")
    print(f"Number of ancestry features: {len(ancestry_cols)}")
    
    # Create a combined DataFrame with cluster assignments
    combined = pd.concat([tsne_df, cluster_df], axis=1)
    
    # Analyze the K-means clusters
    kmeans_clusters = cluster_df['K-means Cluster'].value_counts().sort_index()
    print("\nK-means cluster distribution:")
    print(kmeans_clusters)
    
    # Plot t-SNE with cluster coloring
    plt.figure(figsize=(12, 10))
    
    scatter = plt.scatter(
        tsne_df['t-SNE1'],
        tsne_df['t-SNE2'],
        c=cluster_df['K-means Cluster'],
        cmap='viridis',
        alpha=0.7,
        s=50
    )
    
    plt.colorbar(scatter, label='Cluster')
    plt.title('t-SNE Plot Colored by K-means Clusters', fontsize=14)
    plt.xlabel('t-SNE1', fontsize=12)
    plt.ylabel('t-SNE2', fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, 'simple_tsne_clusters.png'), dpi=300)
    plt.close()
    
    # Analyze mutation differences between clusters
    mutation_summary = {}
    for cluster_id in sorted(cluster_df['K-means Cluster'].unique()):
        # Get samples in this cluster
        cluster_samples = cluster_df[cluster_df['K-means Cluster'] == cluster_id].index
        
        # Calculate mutation frequency in this cluster
        cluster_mutation_freq = mutations_df.loc[cluster_samples].mean()
        
        # Store the top mutations
        top_mutations = cluster_mutation_freq.sort_values(ascending=False).head(5)
        mutation_summary[cluster_id] = top_mutations
    
    # Plot top mutations by cluster
    plt.figure(figsize=(15, 10))
    
    for i, (cluster_id, mutations) in enumerate(mutation_summary.items()):
        plt.subplot(2, 2, i+1)
        mutations.plot(kind='bar', color=plt.cm.viridis(i/4))
        plt.title(f'Top Mutations in Cluster {cluster_id}')
        plt.ylabel('Frequency')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
    
    plt.savefig(os.path.join(RESULTS_DIR, 'simple_top_mutations.png'), dpi=300)
    plt.close()
    
    # Analyze ancestry differences
    ancestry_by_cluster = {}
    
    for cluster_id in sorted(cluster_df['K-means Cluster'].unique()):
        # Get samples in this cluster
        cluster_samples = cluster_df[cluster_df['K-means Cluster'] == cluster_id].index
        
        # Calculate ancestry composition in this cluster
        cluster_ancestry = ancestry_df.loc[cluster_samples].mean()
        ancestry_by_cluster[cluster_id] = cluster_ancestry
    
    # Convert to DataFrame
    ancestry_summary = pd.DataFrame(ancestry_by_cluster).T
    
    # Plot ancestry composition by cluster
    plt.figure(figsize=(12, 8))
    
    ancestry_summary.plot(kind='bar', stacked=True, cmap='viridis')
    plt.title('Ancestry Composition by Cluster', fontsize=14)
    plt.xlabel('Cluster', fontsize=12)
    plt.ylabel('Average Value', fontsize=12)
    plt.legend(title='Ancestry Component')
    plt.tight_layout()
    
    plt.savefig(os.path.join(RESULTS_DIR, 'simple_ancestry_by_cluster.png'), dpi=300)
    plt.close()
    
    # Create a composite analysis of top features for each cluster
    plt.figure(figsize=(20, 15))
    
    # Plot 1: t-SNE with clusters
    plt.subplot(2, 2, 1)
    scatter = plt.scatter(
        tsne_df['t-SNE1'],
        tsne_df['t-SNE2'],
        c=cluster_df['K-means Cluster'],
        cmap='viridis',
        alpha=0.7,
        s=50
    )
    plt.colorbar(scatter, label='Cluster')
    plt.title('t-SNE Plot Colored by K-means Clusters', fontsize=14)
    plt.xlabel('t-SNE1', fontsize=12)
    plt.ylabel('t-SNE2', fontsize=12)
    
    # Plot 2: Mutation differences
    plt.subplot(2, 2, 2)
    mutation_avg = pd.DataFrame({f'Cluster {k}': v for k, v in mutation_summary.items()})
    mutation_avg.plot(kind='bar', ax=plt.gca())
    plt.title('Top Mutations by Cluster', fontsize=14)
    plt.ylabel('Frequency', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    
    # Plot 3: CNA differences
    plt.subplot(2, 2, 3)
    cna_summary = {}
    
    for cluster_id in sorted(cluster_df['K-means Cluster'].unique()):
        # Get samples in this cluster
        cluster_samples = cluster_df[cluster_df['K-means Cluster'] == cluster_id].index
        
        # Calculate mean CNA values
        cluster_cna = cna_df.loc[cluster_samples].mean()
        
        # Get top CNAs by absolute value (both amplifications and deletions)
        top_cnas = cluster_cna.reindex(cluster_cna.abs().sort_values(ascending=False).index).head(5)
        cna_summary[cluster_id] = top_cnas
    
    cna_avg = pd.DataFrame({f'Cluster {k}': v for k, v in cna_summary.items()})
    cna_avg.plot(kind='bar', ax=plt.gca())
    plt.title('Top CNAs by Cluster', fontsize=14)
    plt.ylabel('Average CNA Value', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    
    # Plot 4: Ancestry composition
    plt.subplot(2, 2, 4)
    ancestry_summary.plot(kind='bar', stacked=True, cmap='viridis', ax=plt.gca())
    plt.title('Ancestry Composition by Cluster', fontsize=14)
    plt.xlabel('Cluster', fontsize=12)
    plt.ylabel('Average Value', fontsize=12)
    plt.legend(title='Ancestry Component')
    
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, 'simple_cluster_composite.png'), dpi=300)
    plt.close()
    
    # Highlight key genes in a t-SNE plot
    for mutation in ['mutation_TP53', 'mutation_KRAS', 'mutation_APC']:
        if mutation in mutations_df.columns:
            plt.figure(figsize=(12, 10))
            
            # Color by cluster
            scatter = plt.scatter(
                tsne_df['t-SNE1'],
                tsne_df['t-SNE2'],
                c=cluster_df['K-means Cluster'],
                cmap='viridis',
                alpha=0.7,
                s=50
            )
            
            # Highlight mutations
            mutated_samples = mutations_df.index[mutations_df[mutation] == 1]
            plt.scatter(
                tsne_df.loc[mutated_samples, 't-SNE1'],
                tsne_df.loc[mutated_samples, 't-SNE2'],
                s=80,
                facecolors='none',
                edgecolors='black',
                linewidths=1.5,
                alpha=0.8,
                label=f"{mutation.replace('mutation_', '')} Mutated"
            )
            
            plt.colorbar(scatter, label='Cluster')
            plt.title(f't-SNE Plot with {mutation.replace("mutation_", "")} Mutations Highlighted', fontsize=14)
            plt.xlabel('t-SNE1', fontsize=12)
            plt.ylabel('t-SNE2', fontsize=12)
            plt.legend()
            plt.tight_layout()
            
            plt.savefig(os.path.join(RESULTS_DIR, f'simple_tsne_{mutation}.png'), dpi=300)
            plt.close()
    
    print(f"All results have been saved to the '{RESULTS_DIR}' directory.")

if __name__ == "__main__":
    main()