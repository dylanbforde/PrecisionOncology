#!/usr/bin/env python
"""
This script performs unsupervised learning on the merged data.
It applies dimensionality reduction and clustering algorithms and generates visualizations.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans, DBSCAN, AgglomerativeClustering
from sklearn.metrics import silhouette_score
import seaborn as sns
import pathlib

# Get the project root directory
PROJECT_ROOT = pathlib.Path(__file__).parent.parent.absolute()

# File paths
DATA_DIR = os.path.join(PROJECT_ROOT, 'data/adapted_data')
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results')
INPUT_FILE = os.path.join(DATA_DIR, 'merged_data.csv')

# Make sure results directory exists
os.makedirs(RESULTS_DIR, exist_ok=True)

def load_and_preprocess_data():
    """Load and preprocess the merged data."""
    print("Loading data...")

    # Load the merged data
    merged_df = pd.read_csv(INPUT_FILE, index_col=0)

    print(f"Data loaded with shape: {merged_df.shape}")

    # Check for missing values
    missing_values = merged_df.isna().sum().sum()
    print(f"Number of missing values: {missing_values}")

    # Fill missing values
    print("Filling missing values...")
    merged_df = merged_df.fillna(0)  # Replace NaN with 0

    # Split the data into different types
    mutation_cols = [col for col in merged_df.columns if col.startswith('mutation_')]
    cna_cols = [col for col in merged_df.columns if col.startswith('cna_')]
    mrna_cols = [col for col in merged_df.columns if col.startswith('mrna_')]
    ancestry_cols = [col for col in merged_df.columns if col.startswith('ancestry_')]

    print(f"Number of mutation features: {len(mutation_cols)}")
    print(f"Number of CNA features: {len(cna_cols)}")
    print(f"Number of mRNA features: {len(mrna_cols)}")
    print(f"Number of ancestry features: {len(ancestry_cols)}")

    # Create separate DataFrames for each data type
    mutations_df = merged_df[mutation_cols]
    cna_df = merged_df[cna_cols]
    mrna_df = merged_df[mrna_cols]
    ancestry_df = merged_df[ancestry_cols]

    # Standardize the data (except for mutations, which are already binary)
    print("Standardizing data...")
    scaler = StandardScaler()

    cna_scaled = pd.DataFrame(
        scaler.fit_transform(cna_df),
        index=cna_df.index,
        columns=cna_df.columns
    )

    mrna_scaled = pd.DataFrame(
        scaler.fit_transform(mrna_df),
        index=mrna_df.index,
        columns=mrna_df.columns
    )

    ancestry_scaled = pd.DataFrame(
        scaler.fit_transform(ancestry_df),
        index=ancestry_df.index,
        columns=ancestry_df.columns
    )

    # Combine the scaled data
    data_scaled = pd.concat([mutations_df, cna_scaled, mrna_scaled, ancestry_scaled], axis=1)

    return data_scaled, mutations_df, cna_scaled, mrna_scaled, ancestry_scaled

def perform_dimensionality_reduction(data_scaled):
    """Perform PCA and t-SNE for dimensionality reduction."""
    print("Performing PCA...")
    
    # PCA
    pca = PCA(n_components=50)
    pca_result = pca.fit_transform(data_scaled)
    
    # Plot PCA variance explained
    plt.figure(figsize=(10, 6))
    plt.plot(np.cumsum(pca.explained_variance_ratio_))
    plt.xlabel('Number of Components')
    plt.ylabel('Cumulative Explained Variance')
    plt.title('PCA Explained Variance')
    plt.grid(True)
    plt.savefig(os.path.join(RESULTS_DIR, 'pca_variance_explained.png'), dpi=300)
    plt.close()
    
    # Create a DataFrame with PCA results
    pca_df = pd.DataFrame(
        pca_result,
        index=data_scaled.index,
        columns=[f'PC{i+1}' for i in range(pca_result.shape[1])]
    )
    
    print("Performing t-SNE...")
    
    # t-SNE on PCA results
    tsne = TSNE(n_components=2, random_state=42, perplexity=30, n_iter=1000)
    tsne_result = tsne.fit_transform(pca_result[:, :20])  # Use first 20 PCs for t-SNE
    
    # Create a DataFrame with t-SNE results
    tsne_df = pd.DataFrame(
        tsne_result,
        index=data_scaled.index,
        columns=['t-SNE1', 't-SNE2']
    )
    
    return pca_df, tsne_df

def perform_clustering(pca_df, tsne_df):
    """Perform various clustering algorithms and evaluate them."""
    print("Performing clustering...")
    
    # Use PCA results for clustering (first 20 components)
    X_cluster = pca_df.iloc[:, :20]
    
    # K-means clustering for different numbers of clusters
    silhouette_scores = []
    for n_clusters in range(2, 11):
        kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        cluster_labels = kmeans.fit_predict(X_cluster)
        silhouette_avg = silhouette_score(X_cluster, cluster_labels)
        silhouette_scores.append(silhouette_avg)
        print(f"K-means with {n_clusters} clusters: Silhouette Score = {silhouette_avg:.3f}")
    
    # Plot silhouette scores
    plt.figure(figsize=(10, 6))
    plt.plot(range(2, 11), silhouette_scores, 'o-')
    plt.xlabel('Number of Clusters')
    plt.ylabel('Silhouette Score')
    plt.title('K-means Clustering: Silhouette Score vs. Number of Clusters')
    plt.grid(True)
    plt.savefig(os.path.join(RESULTS_DIR, 'kmeans_silhouette_scores.png'), dpi=300)
    plt.close()
    
    # Find optimal number of clusters based on silhouette score
    optimal_n_clusters = silhouette_scores.index(max(silhouette_scores)) + 2
    print(f"Optimal number of clusters based on silhouette score: {optimal_n_clusters}")
    
    # Perform K-means with optimal number of clusters
    kmeans = KMeans(n_clusters=optimal_n_clusters, random_state=42, n_init=10)
    kmeans_labels = kmeans.fit_predict(X_cluster)
    
    # Add cluster labels to t-SNE results
    tsne_df['K-means Cluster'] = kmeans_labels
    
    # Also try DBSCAN
    dbscan = DBSCAN(eps=3, min_samples=5)
    dbscan_labels = dbscan.fit_predict(X_cluster)
    
    # Add DBSCAN labels to t-SNE results
    tsne_df['DBSCAN Cluster'] = dbscan_labels
    
    # Try hierarchical clustering
    hc = AgglomerativeClustering(n_clusters=optimal_n_clusters)
    hc_labels = hc.fit_predict(X_cluster)
    
    # Add hierarchical clustering labels to t-SNE results
    tsne_df['Hierarchical Cluster'] = hc_labels
    
    return tsne_df, optimal_n_clusters

def visualize_clusters(tsne_df, data_scaled, optimal_n_clusters):
    """Create visualizations of the clusters."""
    print("Creating visualizations...")

    # Set a color palette
    palette = sns.color_palette("viridis", optimal_n_clusters)

    # K-means clustering visualization
    plt.figure(figsize=(10, 8))
    sns.scatterplot(
        x='t-SNE1', y='t-SNE2',
        hue='K-means Cluster',
        palette=palette,
        data=tsne_df,
        legend='full',
        alpha=0.8
    )
    plt.title('t-SNE Visualization with K-means Clustering')
    plt.savefig(os.path.join(RESULTS_DIR, 'tsne_kmeans_clusters.png'), dpi=300)
    plt.close()

    # DBSCAN clustering visualization
    plt.figure(figsize=(10, 8))
    sns.scatterplot(
        x='t-SNE1', y='t-SNE2',
        hue='DBSCAN Cluster',
        palette='Set1',
        data=tsne_df,
        legend='full',
        alpha=0.8
    )
    plt.title('t-SNE Visualization with DBSCAN Clustering')
    plt.savefig(os.path.join(RESULTS_DIR, 'tsne_dbscan_clusters.png'), dpi=300)
    plt.close()

    # Hierarchical clustering visualization
    plt.figure(figsize=(10, 8))
    sns.scatterplot(
        x='t-SNE1', y='t-SNE2',
        hue='Hierarchical Cluster',
        palette=palette,
        data=tsne_df,
        legend='full',
        alpha=0.8
    )
    plt.title('t-SNE Visualization with Hierarchical Clustering')
    plt.savefig(os.path.join(RESULTS_DIR, 'tsne_hierarchical_clusters.png'), dpi=300)
    plt.close()

    # Analyze clusters - average values for each feature type by cluster
    print("Analyzing clusters...")
    kmeans_cluster_df = tsne_df.copy()
    kmeans_cluster_df = pd.concat([kmeans_cluster_df, data_scaled], axis=1)

    # Get feature types
    feature_types = ['mutation', 'cna', 'mrna', 'ancestry']

    for feature_type in feature_types:
        # Get columns of this feature type
        cols = [col for col in data_scaled.columns if col.startswith(f'{feature_type}_')]

        if cols:
            # For large feature sets, select only the most important features
            if len(cols) > 50:
                print(f"Selecting top features for {feature_type}...")
                # Calculate variance of each feature across all samples
                feature_variance = data_scaled[cols].var()
                # Sort features by variance in descending order
                sorted_features = feature_variance.sort_values(ascending=False)
                # Select top 30 features
                top_features = sorted_features.head(30).index.tolist()
                cols = top_features

            # Calculate average value for each cluster
            cluster_means = kmeans_cluster_df.groupby('K-means Cluster')[cols].mean()

            # Heatmap of cluster means
            plt.figure(figsize=(16, 10))
            ax = sns.heatmap(cluster_means, cmap='viridis', linewidths=.5)

            # Improve x-axis labels readability
            if len(cols) > 10:
                plt.xticks(rotation=90, fontsize=8)
                plt.tight_layout()

            plt.title(f'Average {feature_type.upper()} Values by Cluster')
            plt.savefig(os.path.join(RESULTS_DIR, f'cluster_means_{feature_type}.png'), dpi=300)
            plt.close()

            # Save the cluster means for this feature type
            cluster_means.to_csv(os.path.join(RESULTS_DIR, f'cluster_means_{feature_type}.csv'))

    # Save cluster assignments
    cluster_assignments = tsne_df[['K-means Cluster', 'DBSCAN Cluster', 'Hierarchical Cluster']]
    cluster_assignments.to_csv(os.path.join(RESULTS_DIR, 'cluster_assignments.csv'))

    # Save t-SNE results with cluster assignments
    tsne_df.to_csv(os.path.join(RESULTS_DIR, 'tsne_results.csv'))

def main():
    # Load and preprocess the data
    data_scaled, mutations_df, cna_scaled, mrna_scaled, ancestry_scaled = load_and_preprocess_data()
    
    # Perform dimensionality reduction
    pca_df, tsne_df = perform_dimensionality_reduction(data_scaled)
    
    # Perform clustering
    tsne_df, optimal_n_clusters = perform_clustering(pca_df, tsne_df)
    
    # Create visualizations
    visualize_clusters(tsne_df, data_scaled, optimal_n_clusters)
    
    print(f"All results have been saved to the '{RESULTS_DIR}' directory.")

if __name__ == "__main__":
    main()