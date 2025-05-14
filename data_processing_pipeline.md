# Precision Oncology Data Processing Pipeline

This document outlines the complete data processing pipeline used in this project, from raw data preparation to clustering analysis and visualization. It provides a step-by-step explanation of how data was processed, transformed, and visualized to generate the analysis outputs.

## 1. Data Preparation

### 1.1 Raw Data Sources

The pipeline processes four distinct types of cancer data:
- **Mutation data**: Somatic mutations for each patient
- **Copy Number Alteration (CNA) data**: Gene-level copy number changes
- **mRNA expression data**: Gene expression levels
- **Ancestry data**: Genetic ancestry components

### 1.2 Data Processing Workflow

The data preparation is orchestrated by the `run_data_preparation.sh` script, which executes:

```bash
# Step 1: Run individual data preparation scripts
python src/data_preparation/prepare_all_data.py

# Step 2: Merge all prepared data
python src/data_preparation/merge_data.py
```

#### 1.2.1 Individual Data Type Processing

Each data type is processed by a dedicated script:

- **Mutation Processing** (`prepare_mutations.py`):
  - Creates a binary matrix (1 = mutation present, 0 = absent)
  - Maps sample IDs to patient IDs
  - Outputs `mutations_matrix.csv`

- **CNA Processing** (`prepare_cna.py`):
  - Processes copy number alteration data
  - Outputs `cna_matrix.csv`

- **mRNA Processing** (`prepare_mrna.py`):
  - Processes mRNA expression data
  - Outputs `mrna_matrix.csv`

- **Ancestry Processing** (`prepare_ancestry.py`):
  - Processes ancestry/genetic background data
  - Outputs `ancestry_matrix.csv`

#### 1.2.2 Data Merging

The `merge_data.py` script combines all prepared datasets:

1. **Mutation Filtering**:
   - Calculates frequency for each mutation
   - Retains mutations with 5%-95% frequency
   - Discards very rare or very common mutations

2. **Feature Selection for CNA and mRNA**:
   - Calculates variance for each gene
   - Selects top 1000 most variable genes
   - Reduces dimensionality while preserving informative features

3. **Data Standardization**:
   - Transposes matrices to have patients as rows and features as columns
   - Adds data type prefixes to column names (e.g., "mutation_", "cna_")

4. **Data Integration**:
   - Identifies patients common to all datasets
   - Creates a unified matrix with all data types
   - Outputs `merged_data.csv`

## 2. Clustering Analysis

The clustering analysis is performed by the `clustering_analysis.py` script.

### 2.1 Data Preprocessing

1. **Loading and Cleaning**:
   - Loads the merged dataset
   - Fills missing values with zeros
   - Separates data by type (mutations, CNA, mRNA, ancestry)

2. **Standardization**:
   - Applies `StandardScaler` to CNA, mRNA, and ancestry data
   - Keeps mutation data as binary (0/1)
   - Standardizes each data type separately to preserve relative importance

### 2.2 Dimensionality Reduction

1. **Principal Component Analysis (PCA)**:
   - Reduces data to 50 principal components
   - Calculates and visualizes cumulative explained variance
   - Generates `pca_variance_explained.png`

2. **t-SNE (t-Distributed Stochastic Neighbor Embedding)**:
   - Takes first 20 PCA components as input
   - Reduces to 2 dimensions for visualization
   - Parameters: perplexity=30, n_iter=1000
   - Preserves local structure for better visualization

### 2.3 Clustering Algorithms

1. **K-means Clustering**:
   - Tests cluster counts from 2 to 10
   - Calculates silhouette scores for each
   - Selects optimal number of clusters based on highest silhouette score
   - Generates `kmeans_silhouette_scores.png`

2. **DBSCAN (Density-Based Spatial Clustering)**:
   - Parameters: eps=3, min_samples=5
   - Identifies clusters based on density

3. **Hierarchical Clustering**:
   - Uses agglomerative clustering
   - Sets cluster count to match K-means optimal count

### 2.4 Cluster Analysis and Visualization

1. **t-SNE Visualizations**:
   - Colored by K-means clusters: `tsne_kmeans_clusters.png`
   - Colored by DBSCAN clusters: `tsne_dbscan_clusters.png`
   - Colored by hierarchical clusters: `tsne_hierarchical_clusters.png`

2. **Cluster Characterization**:
   - Calculates mean values for each feature type by cluster
   - For large feature sets, selects top 30 most variable features
   - Generates heatmaps for each data type:
     - `cluster_means_mutation.png`
     - `cluster_means_cna.png`
     - `cluster_means_mrna.png`
     - `cluster_means_ancestry.png`

3. **Result Storage**:
   - Saves cluster assignments: `cluster_assignments.csv`
   - Saves t-SNE coordinates with cluster labels: `tsne_results.csv`
   - Saves cluster means for each data type as CSVs

## 3. Simplified Cluster Analysis

The `simple_cluster_analysis.py` script provides focused analysis and more accessible visualizations.

### 3.1 Key Analysis Components

1. **Cluster Distribution**:
   - Counts samples per cluster
   - Visualizes t-SNE plot with cluster coloring: `simple_tsne_clusters.png`

2. **Mutation Analysis**:
   - Calculates mutation frequency by cluster
   - Identifies top 5 mutations per cluster
   - Creates bar charts: `simple_top_mutations.png`

3. **Ancestry Analysis**:
   - Calculates average ancestry composition by cluster
   - Creates stacked bar chart: `simple_ancestry_by_cluster.png`

4. **Composite Visualization**:
   - Combines t-SNE, mutations, CNA, and ancestry in one figure
   - Provides comprehensive overview of cluster characteristics
   - Generates `simple_cluster_composite.png`

5. **Key Gene Highlighting**:
   - Creates t-SNE plots highlighting key cancer genes:
     - TP53: `simple_tsne_mutation_TP53.png`
     - KRAS: `simple_tsne_mutation_KRAS.png`
     - APC: `simple_tsne_mutation_APC.png`

## 4. Summary of Analysis Outputs

### 4.1 Clustering Results

- **Optimal number of clusters**: Determined by silhouette score
- **Cluster assignments**: Available in `cluster_assignments.csv`
- **Cluster profiles**: Characterized by mutations, CNA, mRNA, and ancestry

### 4.2 Key Visualizations

- **Dimensionality Reduction**:
  - PCA explained variance: `pca_variance_explained.png`
  - t-SNE cluster plots: `tsne_*_clusters.png`

- **Cluster Characteristics**:
  - Heatmaps of feature means: `cluster_means_*.png`
  - Top mutations by cluster: `simple_top_mutations.png`
  - Ancestry composition: `simple_ancestry_by_cluster.png`

- **Feature Highlights**:
  - Key mutation visualizations: `simple_tsne_mutation_*.png`
  - Composite analysis: `simple_cluster_composite.png`

This pipeline integrates multiple data types to identify patient clusters with distinct molecular profiles, potentially informing precision oncology approaches.