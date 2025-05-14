#!/usr/bin/env python
"""
This script performs survival analysis on the identified clusters.
It generates Kaplan-Meier survival curves for each cluster.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import seaborn as sns
import pathlib
from matplotlib.lines import Line2D

# Get the project root directory
PROJECT_ROOT = pathlib.Path(__file__).parent.parent.absolute()

# File paths
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results')
CLINICAL_FILE = os.path.join(PROJECT_ROOT, 'data/tcga_data/coadread_tcga_pan_can_atlas_2018/data_clinical_patient.txt')
CLUSTER_FILE = os.path.join(RESULTS_DIR, 'cluster_assignments.csv')

# Make sure results directory exists
os.makedirs(RESULTS_DIR, exist_ok=True)

def load_data():
    """Load clinical data and cluster assignments."""
    print("Loading clinical data...")
    
    # Load clinical data (skipping metadata rows)
    clinical_df = pd.read_csv(CLINICAL_FILE, sep='\t', skiprows=4)
    
    # Standardize the patient ID format to match cluster assignments
    clinical_df['PATIENT_ID'] = clinical_df['PATIENT_ID'].str.replace('-', '-', regex=False)
    
    # Load cluster assignments
    cluster_df = pd.read_csv(CLUSTER_FILE, index_col=0)
    cluster_df.index.name = 'PATIENT_ID'
    
    # Reset index to have PATIENT_ID as a column
    cluster_df = cluster_df.reset_index()
    
    # Inspect the loaded data
    print(f"Clinical data shape: {clinical_df.shape}")
    print(f"Cluster data shape: {cluster_df.shape}")
    
    return clinical_df, cluster_df

def merge_data(clinical_df, cluster_df):
    """Merge clinical data with cluster assignments."""
    print("Merging clinical data with cluster assignments...")
    
    # Merge based on patient ID
    merged_df = pd.merge(cluster_df, clinical_df, on='PATIENT_ID', how='inner')
    
    print(f"Merged data shape: {merged_df.shape}")
    print(f"Number of unique patients: {merged_df['PATIENT_ID'].nunique()}")
    
    return merged_df

def prepare_survival_data(merged_df):
    """Prepare data for survival analysis."""
    print("Preparing survival data...")
    
    # Extract relevant columns
    # First check what columns are available
    print("\nAvailable columns:", merged_df.columns.tolist())
    
    # Extract columns based on what's available
    try:
        survival_df = merged_df[['PATIENT_ID', 'K-means Cluster', 'Hierarchical Cluster', 'OS_STATUS', 'OS_MONTHS']].copy()
    except KeyError:
        # Fallback to just the essential columns
        survival_df = merged_df[['PATIENT_ID', 'K-means Cluster', 'OS_STATUS', 'OS_MONTHS']].copy()
    
    # Convert OS_STATUS to binary (1 for event, 0 for censored)
    survival_df['event'] = survival_df['OS_STATUS'].apply(lambda x: 1 if x == '1:DECEASED' else 0)
    
    # Rename for clarity
    survival_df.rename(columns={'OS_MONTHS': 'time'}, inplace=True)
    
    # Handle missing values
    survival_df = survival_df.dropna(subset=['time', 'event'])
    
    # Convert time to numeric
    survival_df['time'] = pd.to_numeric(survival_df['time'])
    
    # Print summary
    event_count = survival_df['event'].sum()
    total_count = len(survival_df)
    print(f"Total patients for analysis: {total_count}")
    print(f"Events (deaths): {event_count} ({event_count/total_count*100:.1f}%)")
    print(f"Censored: {total_count - event_count} ({(total_count - event_count)/total_count*100:.1f}%)")
    
    return survival_df

def generate_survival_curves(survival_df):
    """Generate Kaplan-Meier survival curves for each cluster."""
    print("Generating survival curves...")
    
    # Set up the plot
    plt.figure(figsize=(12, 8))
    
    # Initialize the Kaplan-Meier fitter
    kmf = KaplanMeierFitter()
    
    # Define a color palette for clusters
    cluster_ids = sorted(survival_df['K-means Cluster'].unique())
    palette = sns.color_palette("viridis", len(cluster_ids))
    
    # Create survival curves for each cluster
    for i, cluster_id in enumerate(cluster_ids):
        # Filter data for this cluster
        cluster_data = survival_df[survival_df['K-means Cluster'] == cluster_id]
        
        # Fit the KM model
        label = f'Cluster {cluster_id} (n={len(cluster_data)})'
        kmf.fit(cluster_data['time'], cluster_data['event'], label=label)
        
        # Plot the curve
        kmf.plot(ax=plt.gca(), ci_show=True, color=palette[i])
    
    # Add labels and title
    plt.xlabel('Time (Months)', fontsize=12)
    plt.ylabel('Survival Probability', fontsize=12)
    plt.title('Kaplan-Meier Survival Curves by Cluster', fontsize=14)
    plt.grid(alpha=0.3)
    
    # Add legend
    plt.legend(loc='best', frameon=True, framealpha=0.5)
    
    # Save the figure
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, 'survival_curves_by_cluster.png'), dpi=300)
    plt.close()
    
    # Create a version with 95% confidence intervals for better visualization
    plt.figure(figsize=(12, 8))
    
    for i, cluster_id in enumerate(cluster_ids):
        cluster_data = survival_df[survival_df['K-means Cluster'] == cluster_id]
        label = f'Cluster {cluster_id} (n={len(cluster_data)})'
        kmf.fit(cluster_data['time'], cluster_data['event'], label=label)
        kmf.plot(ax=plt.gca(), ci_show=False, color=palette[i])
    
    plt.xlabel('Time (Months)', fontsize=12)
    plt.ylabel('Survival Probability', fontsize=12)
    plt.title('Kaplan-Meier Survival Curves by Cluster (without CIs)', fontsize=14)
    plt.grid(alpha=0.3)
    plt.legend(loc='best', frameon=True, framealpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, 'survival_curves_clean.png'), dpi=300)
    plt.close()
    
    return cluster_ids, palette

def perform_logrank_test(survival_df, cluster_ids):
    """Perform log-rank test to compare survival between clusters."""
    print("Performing log-rank test...")
    
    # Create a dataframe to store results
    results = []
    
    # Perform pairwise log-rank tests
    for i, cluster1 in enumerate(cluster_ids):
        for cluster2 in cluster_ids[i+1:]:
            # Get data for each cluster
            group1 = survival_df[survival_df['K-means Cluster'] == cluster1]
            group2 = survival_df[survival_df['K-means Cluster'] == cluster2]
            
            # Perform log-rank test
            lr_test = logrank_test(
                group1['time'], group2['time'],
                group1['event'], group2['event']
            )
            
            # Store results
            results.append({
                'Cluster 1': cluster1,
                'Cluster 2': cluster2,
                'p-value': lr_test.p_value,
                'Test Statistic': lr_test.test_statistic
            })
    
    # Convert to DataFrame
    results_df = pd.DataFrame(results)
    
    # Print results
    print("\nLog-rank test results (pairwise):")
    print(results_df)
    
    # Save to CSV
    results_df.to_csv(os.path.join(RESULTS_DIR, 'logrank_results.csv'), index=False)
    
    # Create a heatmap of p-values
    if len(cluster_ids) > 2:
        p_value_matrix = np.zeros((len(cluster_ids), len(cluster_ids)))
        
        # Fill the matrix with p-values
        for result in results:
            i = cluster_ids.index(result['Cluster 1'])
            j = cluster_ids.index(result['Cluster 2'])
            p_value_matrix[i, j] = result['p-value']
            p_value_matrix[j, i] = result['p-value']  # Mirror
        
        # Set diagonal to 1 (no difference)
        np.fill_diagonal(p_value_matrix, 1)
        
        # Create heatmap
        plt.figure(figsize=(10, 8))
        sns.heatmap(
            p_value_matrix,
            annot=True,
            fmt='.3f',
            cmap='viridis_r',
            xticklabels=[f'Cluster {c}' for c in cluster_ids],
            yticklabels=[f'Cluster {c}' for c in cluster_ids]
        )
        plt.title('Log-rank Test p-values Between Clusters', fontsize=14)
        plt.tight_layout()
        plt.savefig(os.path.join(RESULTS_DIR, 'logrank_pvalues_heatmap.png'), dpi=300)
        plt.close()
    
    return results_df

def create_composite_plot(survival_df, cluster_ids, palette):
    """Create a composite plot with survival curves and clinical characteristics."""
    print("Creating composite survival visualization...")
    
    # Set up the figure
    fig, axes = plt.subplots(1, 2, figsize=(20, 10))
    
    # 1. Survival curves
    kmf = KaplanMeierFitter()
    
    for i, cluster_id in enumerate(cluster_ids):
        cluster_data = survival_df[survival_df['K-means Cluster'] == cluster_id]
        label = f'Cluster {cluster_id} (n={len(cluster_data)})'
        kmf.fit(cluster_data['time'], cluster_data['event'], label=label)
        kmf.plot(ax=axes[0], ci_show=False, color=palette[i])
    
    axes[0].set_xlabel('Time (Months)', fontsize=12)
    axes[0].set_ylabel('Survival Probability', fontsize=12)
    axes[0].set_title('Kaplan-Meier Survival Curves by Cluster', fontsize=14)
    axes[0].grid(alpha=0.3)
    axes[0].legend(loc='best', frameon=True, framealpha=0.5)
    
    # 2. Median survival time by cluster
    median_survival = []
    
    for cluster_id in cluster_ids:
        cluster_data = survival_df[survival_df['K-means Cluster'] == cluster_id]
        kmf.fit(cluster_data['time'], cluster_data['event'])
        
        try:
            median = kmf.median_survival_time_
        except:
            # If median not reached, use the last timepoint
            median = cluster_data['time'].max()
            
        median_survival.append({
            'Cluster': f'Cluster {cluster_id}',
            'Median Survival (Months)': median
        })
    
    median_df = pd.DataFrame(median_survival)
    
    # Bar plot of median survival
    bars = axes[1].bar(
        median_df['Cluster'],
        median_df['Median Survival (Months)'],
        color=[palette[cluster_ids.index(int(c.split()[-1]))] for c in median_df['Cluster']]
    )
    
    # Add text labels
    for bar in bars:
        height = bar.get_height()
        axes[1].text(
            bar.get_x() + bar.get_width()/2.,
            height + 1,
            f'{height:.1f}',
            ha='center',
            va='bottom',
            fontsize=10
        )
    
    axes[1].set_xlabel('Cluster', fontsize=12)
    axes[1].set_ylabel('Median Survival (Months)', fontsize=12)
    axes[1].set_title('Median Survival Time by Cluster', fontsize=14)
    axes[1].grid(axis='y', alpha=0.3)
    
    # Add p-value annotation if significant differences exist
    try:
        logrank_results = pd.read_csv(os.path.join(RESULTS_DIR, 'logrank_results.csv'))
        significant_results = logrank_results[logrank_results['p-value'] < 0.05]
        
        if not significant_results.empty:
            min_pvalue = logrank_results['p-value'].min()
            axes[0].annotate(
                f'Log-rank p-value: {min_pvalue:.4f}',
                xy=(0.5, 0.05),
                xycoords='axes fraction',
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8)
            )
    except:
        pass
    
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, 'survival_composite.png'), dpi=300)
    plt.close()

def main():
    # Load data
    clinical_df, cluster_df = load_data()
    
    # Merge clinical data with cluster assignments
    merged_df = merge_data(clinical_df, cluster_df)
    
    # Prepare data for survival analysis
    survival_df = prepare_survival_data(merged_df)
    
    # Generate survival curves
    cluster_ids, palette = generate_survival_curves(survival_df)
    
    # Perform log-rank test
    logrank_results = perform_logrank_test(survival_df, cluster_ids)
    
    # Create composite visualization
    create_composite_plot(survival_df, cluster_ids, palette)
    
    print(f"All results have been saved to the '{RESULTS_DIR}' directory.")

if __name__ == "__main__":
    main()