#!/usr/bin/env python3
"""
Cluster significant genes from CELINA results to find co-regulated gene modules.
For each cell type and sample, identifies groups of genes that change together.
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path
import glob
import argparse
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist
import warnings
warnings.filterwarnings('ignore')

def find_optimal_clusters_elbow(data, max_k=10):
    """
    Find optimal number of clusters using elbow method.
    
    Args:
        data: Gene expression data (genes x cells)
        max_k: Maximum number of clusters to test
    
    Returns:
        optimal_k: Optimal number of clusters
        wcss_values: Within-cluster sum of squares for each k
    """
    wcss_values = []
    k_range = range(1, min(max_k + 1, len(data) + 1))
    
    for k in k_range:
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        kmeans.fit(data)
        wcss_values.append(kmeans.inertia_)  # inertia_ is WCSS
    
    # Find elbow point using second derivative
    if len(wcss_values) > 2:
        # Calculate second derivative to find elbow
        second_derivative = np.diff(wcss_values, 2)
        optimal_k = np.argmax(second_derivative) + 2  # +2 because of double diff
    else:
        optimal_k = 2  # Default fallback
    
    return optimal_k, wcss_values

def find_optimal_clusters_silhouette(data, max_k=10):
    """
    Find optimal number of clusters using silhouette method.
    
    Args:
        data: Gene expression data (genes x cells)
        max_k: Maximum number of clusters to test
    
    Returns:
        optimal_k: Optimal number of clusters
        silhouette_scores: Silhouette scores for each k
    """
    silhouette_scores = []
    k_range = range(2, min(max_k + 1, len(data) + 1))
    
    for k in k_range:
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        cluster_labels = kmeans.fit_predict(data)
        silhouette_avg = silhouette_score(data, cluster_labels)
        silhouette_scores.append(silhouette_avg)
    
    if silhouette_scores:
        optimal_k = k_range[np.argmax(silhouette_scores)]
    else:
        optimal_k = 2
    
    return optimal_k, silhouette_scores

def plot_elbow_curve(wcss_values, optimal_k, output_path):
    """Plot elbow curve for cluster selection."""
    plt.figure(figsize=(10, 6))
    plt.plot(range(1, len(wcss_values) + 1), wcss_values, 'bo-')
    plt.axvline(x=optimal_k, color='red', linestyle='--', 
                label=f'Optimal k = {optimal_k}')
    plt.xlabel('Number of Clusters (k)')
    plt.ylabel('Within-Cluster Sum of Squares (WCSS)')
    plt.title('Elbow Method for Optimal Cluster Selection')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def interactive_cluster_selection(wcss_values, silhouette_scores, max_k, sample_name, cell_type):
    """Interactive cluster selection with elbow plot display."""
    print(f"\n{'='*60}")
    print(f"INTERACTIVE CLUSTER SELECTION")
    print(f"Sample: {sample_name}")
    print(f"Cell Type: {cell_type}")
    print(f"{'='*60}")
    
    # Create the elbow plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Elbow plot
    k_range = range(1, len(wcss_values) + 1)
    ax1.plot(k_range, wcss_values, 'bo-', linewidth=2, markersize=8)
    ax1.set_xlabel('Number of Clusters (k)')
    ax1.set_ylabel('Within-Cluster Sum of Squares (WCSS)')
    ax1.set_title('Elbow Method')
    ax1.grid(True, alpha=0.3)
    
    # Silhouette plot
    if silhouette_scores:
        sil_k_range = range(2, len(silhouette_scores) + 2)
        ax2.plot(sil_k_range, silhouette_scores, 'go-', linewidth=2, markersize=8)
        ax2.set_xlabel('Number of Clusters (k)')
        ax2.set_ylabel('Silhouette Score')
        ax2.set_title('Silhouette Method')
        ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    # Get user input
    while True:
        try:
            print(f"\nRecommended number of clusters:")
            print(f"  - Elbow method suggests: {np.argmax(np.diff(wcss_values, 2)) + 2 if len(wcss_values) > 2 else 2}")
            if silhouette_scores:
                print(f"  - Silhouette method suggests: {range(2, len(silhouette_scores) + 2)[np.argmax(silhouette_scores)]}")
            
            user_k = int(input(f"\nEnter number of clusters (2-{max_k}): "))
            
            if 2 <= user_k <= max_k:
                print(f"Selected {user_k} clusters for {sample_name} - {cell_type}")
                return user_k
            else:
                print(f"Please enter a number between 2 and {max_k}")
                
        except ValueError:
            print("Please enter a valid number")
        except KeyboardInterrupt:
            print("\nExiting...")
            return None

def plot_silhouette_curve(silhouette_scores, optimal_k, output_path):
    """Plot silhouette curve for cluster selection."""
    plt.figure(figsize=(10, 6))
    k_range = range(2, len(silhouette_scores) + 2)
    plt.plot(k_range, silhouette_scores, 'go-')
    plt.axvline(x=optimal_k, color='red', linestyle='--', 
                label=f'Optimal k = {optimal_k}')
    plt.xlabel('Number of Clusters (k)')
    plt.ylabel('Silhouette Score')
    plt.title('Silhouette Method for Optimal Cluster Selection')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def perform_clustering_analysis(expression_data, method='kmeans', n_clusters=None):
    """
    Perform clustering analysis on gene expression data.
    
    Args:
        expression_data: Gene expression matrix (genes x cells)
        method: Clustering method ('kmeans' or 'hierarchical')
        n_clusters: Number of clusters (if None, will be determined automatically)
    
    Returns:
        cluster_labels: Cluster assignments for each gene
        cluster_centers: Cluster centers (for k-means)
    """
    # Standardize the data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(expression_data)
    
    if method == 'kmeans':
        if n_clusters is None:
            n_clusters, _ = find_optimal_clusters_elbow(scaled_data)
        
        kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        cluster_labels = kmeans.fit_predict(scaled_data)
        cluster_centers = kmeans.cluster_centers_
        
    elif method == 'hierarchical':
        if n_clusters is None:
            n_clusters, _ = find_optimal_clusters_silhouette(scaled_data)
        
        hierarchical = AgglomerativeClustering(n_clusters=n_clusters)
        cluster_labels = hierarchical.fit_predict(scaled_data)
        cluster_centers = None
    
    return cluster_labels, cluster_centers

def plot_cluster_heatmap(expression_data, cluster_labels, gene_names, output_path):
    """Plot heatmap of clustered genes."""
    # Create DataFrame with cluster assignments
    df = pd.DataFrame(expression_data, index=gene_names)
    df['cluster'] = cluster_labels
    
    # Sort by cluster
    df_sorted = df.sort_values('cluster')
    
    # Remove cluster column for plotting
    plot_data = df_sorted.drop('cluster', axis=1)
    
    # Create the plot
    plt.figure(figsize=(15, 10))
    
    # Create color map for clusters
    n_clusters = len(set(cluster_labels))
    cluster_colors = plt.cm.Set3(np.linspace(0, 1, n_clusters))
    
    # Plot heatmap
    sns.heatmap(plot_data, cmap='viridis', cbar_kws={'label': 'Expression Level'})
    
    # Add cluster boundaries
    cluster_boundaries = []
    current_cluster = df_sorted['cluster'].iloc[0]
    for i, cluster in enumerate(df_sorted['cluster']):
        if cluster != current_cluster:
            cluster_boundaries.append(i)
            current_cluster = cluster
    
    for boundary in cluster_boundaries:
        plt.axhline(y=boundary, color='red', linewidth=2)
    
    plt.title('Gene Expression Heatmap by Clusters')
    plt.xlabel('Cells')
    plt.ylabel('Genes (sorted by cluster)')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def plot_pca_clusters(expression_data, cluster_labels, gene_names, output_path):
    """Plot PCA visualization of clustered genes."""
    # Standardize data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(expression_data)
    
    # Perform PCA
    pca = PCA(n_components=2)
    pca_data = pca.fit_transform(scaled_data)
    
    # Create the plot
    plt.figure(figsize=(12, 8))
    
    # Plot each cluster with different colors
    n_clusters = len(set(cluster_labels))
    colors = plt.cm.Set3(np.linspace(0, 1, n_clusters))
    
    for i in range(n_clusters):
        mask = cluster_labels == i
        plt.scatter(pca_data[mask, 0], pca_data[mask, 1], 
                   c=[colors[i]], label=f'Cluster {i}', alpha=0.7, s=50)
    
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)')
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)')
    plt.title('PCA Visualization of Gene Clusters')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def analyze_cluster_characteristics(expression_data, cluster_labels, gene_names):
    """Analyze characteristics of each cluster."""
    cluster_stats = []
    
    for cluster_id in sorted(set(cluster_labels)):
        mask = cluster_labels == cluster_id
        cluster_genes = [gene_names[i] for i in range(len(gene_names)) if mask[i]]
        cluster_data = expression_data[mask]
        
        stats = {
            'cluster_id': cluster_id,
            'n_genes': len(cluster_genes),
            'genes': cluster_genes,
            'mean_expression': np.mean(cluster_data),
            'std_expression': np.std(cluster_data),
            'max_expression': np.max(cluster_data),
            'min_expression': np.min(cluster_data)
        }
        cluster_stats.append(stats)
    
    return cluster_stats

def load_celina_results(sample_dir):
    """Load CELINA results for a sample."""
    # Use the FDR-corrected p-values file with global Bonferroni correction
    pvalues_file = os.path.join(sample_dir, 'all_p_values_FDR_global_bonferroni.csv')
    if os.path.exists(pvalues_file):
        return pd.read_csv(pvalues_file)
    return None

def load_celltype_matrix(sample_dir, cell_type):
    """Load cell-type-specific gene expression matrix."""
    matrix_file = os.path.join(sample_dir, 'celltype_matrices', 
                              f'gene_cell_matrix_{cell_type}.csv')
    if os.path.exists(matrix_file):
        return pd.read_csv(matrix_file, index_col=0)
    return None

def load_gene_filter_file(filter_file):
    """
    Load gene filter file (CSV with gene and cell_type columns).
    Returns a dictionary mapping cell_type -> set of genes.
    """
    if filter_file is None or not os.path.exists(filter_file):
        return None
    
    try:
        df = pd.read_csv(filter_file)
        if 'gene' not in df.columns or 'cell_type' not in df.columns:
            print(f"Warning: Filter file must have 'gene' and 'cell_type' columns. Ignoring filter.")
            return None
        
        filter_dict = {}
        for _, row in df.iterrows():
            cell_type = row['cell_type']
            gene = row['gene']
            if cell_type not in filter_dict:
                filter_dict[cell_type] = set()
            filter_dict[cell_type].add(gene)
        
        print(f"Loaded gene filter: {sum(len(genes) for genes in filter_dict.values())} gene-cell_type pairs")
        return filter_dict
    except Exception as e:
        print(f"Error loading filter file: {e}. Ignoring filter.")
        return None

def get_significant_genes(celina_results, cell_type=None, pvalue_threshold=0.01, gene_filter=None):
    """Get significant genes from CELINA results for a specific cell type."""
    if celina_results is None:
        return []
    
    # The FDR-corrected file has genes as rows and cell types as columns
    if 'gene' in celina_results.columns:
        if cell_type and cell_type in celina_results.columns:
            # Filter for significant genes in specific cell type
            # Handle NA values by converting to numeric and dropping NaN
            numeric_col = pd.to_numeric(celina_results[cell_type], errors='coerce')
            sig_genes = celina_results[numeric_col < pvalue_threshold]['gene'].tolist()
            
            # Apply gene filter if provided
            if gene_filter is not None and cell_type in gene_filter:
                sig_genes = [gene for gene in sig_genes if gene in gene_filter[cell_type]]
            
            return sig_genes
        else:
            # If no cell type specified, return all genes
            sig_genes = celina_results['gene'].tolist()
            
            # Apply gene filter if provided (for all cell types)
            if gene_filter is not None:
                all_filtered_genes = set()
                for genes in gene_filter.values():
                    all_filtered_genes.update(genes)
                sig_genes = [gene for gene in sig_genes if gene in all_filtered_genes]
            
            return sig_genes
    
    return []

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Cluster significant genes from CELINA results to find co-regulated gene modules.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with default settings
  python cluster_significant_genes.py --input-dir ./example_data --output-dir ./clustering_results
  
  # Interactive mode to choose number of clusters
  python cluster_significant_genes.py --input-dir ./example_data --output-dir ./clustering_results --interactive
  
  # Custom p-value threshold and output directory
  python cluster_significant_genes.py --input-dir ./example_data --output-dir my_results --pvalue-threshold 0.05
  
  # Use hierarchical clustering with more clusters
  python cluster_significant_genes.py --input-dir ./example_data --output-dir ./clustering_results --method hierarchical --max-clusters 15
        """
    )
    
    parser.add_argument('--input-dir', '-i', required=True,
                       help='Directory containing sample folders with CELINA results')
    parser.add_argument('--output-dir', '-o', required=True,
                       help='Output directory for clustering results')
    parser.add_argument('--pvalue-threshold', type=float, default=0.01,
                       help='P-value threshold for significant genes (default: 0.01)')
    parser.add_argument('--max-clusters', type=int, default=10,
                       help='Maximum number of clusters to test (default: 10)')
    parser.add_argument('--method', type=str, choices=['kmeans', 'hierarchical'], 
                       default='kmeans',
                       help='Clustering method (default: kmeans)')
    parser.add_argument('--min-genes', type=int, default=5,
                       help='Minimum number of genes needed for clustering (default: 5)')
    parser.add_argument('--samples', nargs='*', default=None,
                       help='Specific samples to analyze (default: all samples)')
    parser.add_argument('--verbose', action='store_true',
                       help='Enable verbose output')
    parser.add_argument('--interactive', action='store_true',
                       help='Interactive mode: show elbow plot and let user choose number of clusters')
    parser.add_argument('--gene-filter-file', type=str, default=None,
                       help='CSV file with gene and cell_type columns to filter genes')
    
    return parser.parse_args()

def main():
    """Main function to perform clustering analysis."""
    # Parse command line arguments
    args = parse_arguments()
    
    # Load gene filter file if provided
    gene_filter = load_gene_filter_file(args.gene_filter_file)
    
    # Set up output directory
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    
    if args.verbose:
        print(f"Clustering analysis starting with the following settings:")
        print(f"  Input directory: {args.input_dir}")
        print(f"  Output directory: {output_dir}")
        print(f"  P-value threshold: {args.pvalue_threshold}")
        print(f"  Max clusters: {args.max_clusters}")
        print(f"  Clustering method: {args.method}")
        print(f"  Min genes per analysis: {args.min_genes}")
        print(f"  Interactive mode: {'Yes' if args.interactive else 'No'}")
        if args.gene_filter_file:
            print(f"  Gene filter file: {args.gene_filter_file}")
        if args.samples:
            print(f"  Specific samples: {', '.join(args.samples)}")
        else:
            print(f"  Samples: all available")
        print()
    
    # Get all sample directories
    all_sample_dirs = []
    for item in os.listdir(args.input_dir):
        item_path = os.path.join(args.input_dir, item)
        if os.path.isdir(item_path):
            all_sample_dirs.append(item_path)
    
    # Filter samples if specific ones were requested
    if args.samples:
        sample_dirs = []
        for sample_dir in all_sample_dirs:
            sample_name = os.path.basename(sample_dir)
            if sample_name in args.samples:
                sample_dirs.append(sample_dir)
        
        if not sample_dirs:
            print(f"Error: No matching samples found. Available samples:")
            for sample_dir in all_sample_dirs:
                print(f"  - {os.path.basename(sample_dir)}")
            return
    else:
        sample_dirs = all_sample_dirs
    
    results_summary = []
    
    for sample_dir in sample_dirs:
        sample_name = os.path.basename(sample_dir)
        print(f"Processing sample: {sample_name}")
        
        # Load CELINA results
        celina_results = load_celina_results(sample_dir)
        if celina_results is None:
            print(f"  No CELINA results found for {sample_name}")
            continue
        
        # Get available cell types for this sample
        celltype_matrices_dir = os.path.join(sample_dir, 'celltype_matrices')
        if not os.path.exists(celltype_matrices_dir):
            print(f"  No cell-type matrices found for {sample_name}")
            continue
        
        celltype_files = glob.glob(os.path.join(celltype_matrices_dir, 'gene_cell_matrix_*.csv'))
        cell_types = [os.path.basename(f).replace('gene_cell_matrix_', '').replace('.csv', '') 
                     for f in celltype_files]
        
        for cell_type in cell_types:
            print(f"    Processing cell type: {cell_type}")
            
            # Get significant genes for this specific cell type (filtered by gene_filter if provided)
            significant_genes = get_significant_genes(celina_results, cell_type, args.pvalue_threshold, gene_filter)
            if not significant_genes:
                if gene_filter is not None and cell_type in gene_filter:
                    print(f"      No genes found for {cell_type} in {sample_name} after filtering")
                else:
                    print(f"      No significant genes found for {cell_type} in {sample_name}")
                continue
            
            print(f"      Found {len(significant_genes)} significant genes for {cell_type}")
            
            # Load cell-type-specific matrix
            matrix = load_celltype_matrix(sample_dir, cell_type)
            if matrix is None:
                continue
            
            # Filter for significant genes that are available in the matrix
            available_genes = [gene for gene in significant_genes if gene in matrix.index]
            if len(available_genes) < args.min_genes:
                print(f"      Not enough significant genes ({len(available_genes)}) for clustering (minimum: {args.min_genes})")
                continue
            
            # Extract expression data for significant genes
            expression_data = matrix.loc[available_genes].values
            
            # Perform clustering analysis
            try:
                # Determine optimal number of clusters
                optimal_k_elbow, wcss_values = find_optimal_clusters_elbow(expression_data, args.max_clusters)
                optimal_k_silhouette, silhouette_scores = find_optimal_clusters_silhouette(expression_data, args.max_clusters)
                
                # Choose number of clusters
                if args.interactive:
                    # Interactive mode: show plots and get user input
                    user_k = interactive_cluster_selection(wcss_values, silhouette_scores, 
                                                         min(args.max_clusters, len(available_genes) // 2), 
                                                         sample_name, cell_type)
                    if user_k is None:  # User cancelled
                        print(f"      Skipping {cell_type} (user cancelled)")
                        continue
                    optimal_k = user_k
                else:
                    # Automatic mode: use the average of both methods
                    optimal_k = int(np.round((optimal_k_elbow + optimal_k_silhouette) / 2))
                    optimal_k = max(2, min(optimal_k, len(available_genes) // 2))  # Reasonable bounds
                
                # Perform clustering
                cluster_labels, cluster_centers = perform_clustering_analysis(
                    expression_data, method=args.method, n_clusters=optimal_k)
                
                # Create output directory for this sample and cell type
                sample_output_dir = os.path.join(output_dir, sample_name, cell_type)
                os.makedirs(sample_output_dir, exist_ok=True)
                
                # Plot elbow curve
                elbow_plot_path = os.path.join(sample_output_dir, 'elbow_curve.png')
                plot_elbow_curve(wcss_values, optimal_k_elbow, elbow_plot_path)
                
                # Plot silhouette curve
                silhouette_plot_path = os.path.join(sample_output_dir, 'silhouette_curve.png')
                plot_silhouette_curve(silhouette_scores, optimal_k_silhouette, silhouette_plot_path)
                
                # Plot cluster heatmap
                heatmap_path = os.path.join(sample_output_dir, 'cluster_heatmap.png')
                plot_cluster_heatmap(expression_data, cluster_labels, available_genes, heatmap_path)
                
                # Plot PCA visualization
                pca_path = os.path.join(sample_output_dir, 'pca_clusters.png')
                plot_pca_clusters(expression_data, cluster_labels, available_genes, pca_path)
                
                # Analyze cluster characteristics
                cluster_stats = analyze_cluster_characteristics(expression_data, cluster_labels, available_genes)
                
                # Save cluster results
                cluster_results = []
                for i, gene in enumerate(available_genes):
                    cluster_results.append({
                        'gene': gene,
                        'cluster': cluster_labels[i]
                    })
                
                cluster_df = pd.DataFrame(cluster_results)
                cluster_df.to_csv(os.path.join(sample_output_dir, 'cluster_assignments.csv'), index=False)
                
                # Save cluster statistics
                stats_df = pd.DataFrame(cluster_stats)
                stats_df.to_csv(os.path.join(sample_output_dir, 'cluster_statistics.csv'), index=False)
                
                # Create cluster-specific gene lists
                cluster_genes = {}
                for cluster_id in range(optimal_k):
                    cluster_genes[f'cluster_{cluster_id}_genes'] = []
                
                # Read cluster assignments to get genes per cluster
                cluster_assignments_file = os.path.join(sample_output_dir, 'cluster_assignments.csv')
                if os.path.exists(cluster_assignments_file):
                    cluster_df = pd.read_csv(cluster_assignments_file)
                    for _, row in cluster_df.iterrows():
                        gene = row['gene']
                        cluster_id = row['cluster']
                        cluster_genes[f'cluster_{cluster_id}_genes'].append(gene)
                
                # Convert lists to comma-separated strings
                cluster_genes_str = {}
                for cluster_id in range(optimal_k):
                    genes_list = cluster_genes[f'cluster_{cluster_id}_genes']
                    cluster_genes_str[f'cluster_{cluster_id}_genes'] = ', '.join(genes_list)
                
                # Add to summary
                summary_entry = {
                    'sample': sample_name,
                    'cell_type': cell_type,
                    'n_significant_genes': len(available_genes),
                    'n_clusters': optimal_k,
                    'optimal_k_elbow': optimal_k_elbow,
                    'optimal_k_silhouette': optimal_k_silhouette,
                    'genes': available_genes
                }
                # Add cluster-specific gene columns
                summary_entry.update(cluster_genes_str)
                results_summary.append(summary_entry)
                
                if args.verbose:
                    print(f"      Clustered {len(available_genes)} genes into {optimal_k} clusters using {args.method}")
                else:
                    print(f"      Clustered {len(available_genes)} genes into {optimal_k} clusters")
                
            except Exception as e:
                print(f"      Error clustering {cell_type}: {str(e)}")
                continue
    
    # Save overall summary
    if results_summary:
        summary_df = pd.DataFrame(results_summary)
        summary_df.to_csv(os.path.join(output_dir, 'clustering_summary.csv'), index=False)
    
    print(f"\nClustering analysis complete! Results saved to {output_dir}")
    print(f"Processed {len(results_summary)} sample-cell type combinations")

if __name__ == "__main__":
    main()
