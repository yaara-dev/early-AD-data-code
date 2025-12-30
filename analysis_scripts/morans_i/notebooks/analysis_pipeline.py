# -*- coding: utf-8 -*-
"""
Moran's I Statistical Analysis Pipeline

This notebook performs comprehensive statistical analysis on Moran's I spatial autocorrelation results:
1. Loads concatenated regional results (from concatenate_regional_csvs.py)
2. Applies quantile normalization to both counts and Moran's I data
3. Performs permutation tests to compare WT vs 5X samples
4. Applies FDR correction (Benjamini-Hochberg method)
5. Generates visualizations (heatmaps, volcano plots, histograms)
6. Identifies significant genes and analyzes direction of effects

This is part of the Moran's I spatial autocorrelation analysis publication pipeline.

Use #%% cells for interactive execution in VS Code or PyCharm
"""

#%% Imports and Setup
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats
from statsmodels.stats.multitest import multipletests
from scipy.stats import norm
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

# Define the base path for your regional data
# ============================================================================
# CONFIGURATION: Choose which data to use
# ============================================================================

# Option 1: Use EXAMPLE DATA (for testing)
USE_EXAMPLE_DATA = True  # Set to True to use example data, False for real data

if USE_EXAMPLE_DATA:
    # Example data path (relative to notebook location)
    script_dir = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
    base_path = os.path.join(script_dir, "..", "example_output", "concatenated_results")
    base_path = os.path.abspath(base_path)
    regions = ['CA1']  # Example data only has CA1
else:
    # REAL DATA path - Update this to your actual data location
    # Default assumes data is in the parent morans_i directory
    script_dir = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
    # Go up two levels from notebooks/ to get to morans_i_publication/, then up one more to morans_i/
    parent_dir = os.path.dirname(os.path.dirname(script_dir))
    base_path = os.path.join(parent_dir, "results", "fovs", "grid_based", "binary_6", "concatenated_results")
    base_path = os.path.abspath(base_path)
    
    # All available regions in real data
    regions = ['CA1', 'CA3', 'SLM', 'SM', 'DG', 'upper_CA1', 'inner_DG', 'under_DG']

# Alternative: Specify absolute path directly (uncomment and modify if needed)
# base_path = "/Users/yaarakarasik/morans_i/results/fovs/grid_based/binary_6/concatenated_results"

print("Starting analysis with the following configuration:")
print(f"Base path: {base_path}")
print(f"Regions to analyze: {regions}")

# Flag to skip statistical analysis (only generate heatmaps)
SKIP_STATISTICAL_ANALYSIS = False  # Set to True to skip statistical analysis

# Number of permutations for statistical tests
NUM_PERMUTATIONS = 1000  # Increase for more precision (e.g., 10000), decrease for faster testing (e.g., 100)

#%% Data Loading Functions
# Function to load data from your CSV files
def load_regional_data():
    """
    Load Moran's I and counts data from the separate CSV files for each region
    """
    all_regions_morans = {}
    all_regions_counts = {}
    
    print(f"Base path: {base_path}")
    print(f"Base path exists: {os.path.exists(base_path)}")
    print(f"Regions to process: {regions}")
    
    for region in regions:
        morans_file = os.path.join(base_path, f"{region}_all_samples_morans_concatenated.csv")
        counts_file = os.path.join(base_path, f"{region}_all_samples_counts_concatenated.csv")
        
        if os.path.exists(morans_file) and os.path.exists(counts_file):
            print(f"Loading data for {region}")
            
            # Load Moran's I data
            morans_df = pd.read_csv(morans_file, index_col=0)
            all_regions_morans[region] = morans_df
            
            # Load counts data
            counts_df = pd.read_csv(counts_file, index_col=0)
            all_regions_counts[region] = counts_df
            
            print(f"  {region} Moran's I: {morans_df.shape}")
            print(f"  {region} Counts: {counts_df.shape}")
        else:
            print(f"Warning: Files not found for {region}")
    
    return all_regions_morans, all_regions_counts

# Load the data
print("\nLoading regional data...")
all_regions_morans, all_regions_counts = load_regional_data()

print(f"\nLoaded data for {len(all_regions_morans)} regions")

#%% Sample Group Identification
# Function to identify WT vs 5X samples based on column names
def identify_sample_groups(df):
    """
    Identify WT and 5X sample columns based on column naming patterns
    """
    wt_cols = [col for col in df.columns if any(x in col.upper() for x in ['WT', 'FEM2_WT', 'FEM3_WTE', 'FEM4_WT'])]
    five_x_cols = [col for col in df.columns if any(x in col.upper() for x in ['5X', 'FEM2_5X', 'FEM3_5X', 'FEM4_5X'])]
    
    return wt_cols, five_x_cols

# Test the function on one region
if regions and regions[0] in all_regions_morans:
    test_df = all_regions_morans[regions[0]]
    wt_cols, five_x_cols = identify_sample_groups(test_df)
    print(f"\nSample identification test for {regions[0]}:")
    print(f"  WT samples: {len(wt_cols)}")
    print(f"  5X samples: {len(five_x_cols)}")
    if wt_cols:
        print(f"  Example WT column: {wt_cols[0]}")
    if five_x_cols:
        print(f"  Example 5X column: {five_x_cols[0]}")
        
    # Show a few sample column names
    print(f"\nFirst 10 column names in {regions[0]}:")
    print(list(test_df.columns[:10]))

#%% Data Validation and Cleaning - Before Normalization (SKIPPED - using pre-normalized data)
# Skip data cleaning since we're loading pre-normalized matrices
if not SKIP_STATISTICAL_ANALYSIS:
    print("\n=== DATA VALIDATION AND CLEANING (Before Normalization) ===")
else:
    print("\nSkipping data validation/cleaning - using pre-normalized matrices")

def validate_and_clean_data(all_regions_morans, all_regions_counts):
    """
    Validate and clean data before normalization:
    - Moran's I: Keep NaN as missing data (exclude from analysis)
    - Counts: Convert NaN to 0 (no expression detected)
    """
    print("Validating and cleaning data...")
    
    cleaned_morans = {}
    cleaned_counts = {}
    
    for region in all_regions_morans.keys():
        if region in all_regions_counts:
            print(f"\nProcessing {region}...")
            
            # Get the data
            morans_df = all_regions_morans[region].copy()
            counts_df = all_regions_counts[region].copy()
            
            # Check initial state
            morans_nan_before = morans_df.isnull().sum().sum()
            counts_nan_before = counts_df.isnull().sum().sum()
            
            print(f"  Moran's I - NaN values before cleaning: {morans_nan_before:,}")
            print(f"  Counts - NaN values before cleaning: {counts_nan_before:,}")
            
            # For Moran's I: Keep NaN values (they represent missing data)
            # Just ensure they're properly formatted as NaN (not empty strings)
            morans_df = morans_df.replace('', np.nan)
            morans_df = morans_df.replace(' ', np.nan)
            
            # For Counts: Convert NaN to 0 (no expression detected)
            counts_df = counts_df.fillna(0)
            
            # Check final state
            morans_nan_after = morans_df.isnull().sum().sum()
            counts_nan_after = counts_df.isnull().sum().sum()
            
            print(f"  Moran's I - NaN values after cleaning: {morans_nan_after:,}")
            print(f"  Counts - NaN values after cleaning: {counts_nan_after:,}")
            
            # Store cleaned data
            cleaned_morans[region] = morans_df
            cleaned_counts[region] = counts_df
            
            print(f"  ✓ {region} cleaned")
    
    return cleaned_morans, cleaned_counts

# Clean and validate the data
print("Cleaning and validating data before normalization...")
cleaned_morans, cleaned_counts = validate_and_clean_data(all_regions_morans, all_regions_counts)

# Update the main data dictionaries with cleaned data
print("\nUpdating main data dictionaries with cleaned data...")
all_regions_morans = cleaned_morans
all_regions_counts = cleaned_counts

print("✓ Data cleaning and validation completed!")

# Show summary of cleaning results
print("\n=== CLEANING RESULTS SUMMARY ===")
total_morans_nan = sum(df.isnull().sum().sum() for df in all_regions_morans.values())
total_counts_nan = sum(df.isnull().sum().sum() for df in all_regions_counts.values())

print(f"Total NaN values after cleaning:")
print(f"  Moran's I: {total_morans_nan:,} (missing data - will be excluded from analysis)")
print(f"  Counts: {total_counts_nan:,} (should be 0 - all converted to 0)")

if total_counts_nan == 0:
    print("✓ Counts data properly cleaned - no NaN values remain")
else:
    print(f"⚠ Warning: Counts still has {total_counts_nan} NaN values")

    print("\nData is now ready for normalization!")

#%% Quantile Normalization (SKIPPED - loading pre-normalized data)
# Quantile normalization function (keep the same as original)
def quantileNormalize(df_input):
    df = df_input.copy()
    #compute rank
    dic = {}
    for col in df:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis = 1).tolist()
    #sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df

#%% Quantile Normalization
# Load existing normalized matrices OR compute normalization if they don't exist
print("\n=== QUANTILE NORMALIZATION ===")
normalized_matrices_dir = f"{base_path}/normalized_matrices"
all_regions_normalized_moransi = {}
all_regions_normalized_counts = {}

# Try to load existing normalized matrices first
print("Checking for existing normalized matrices...")
regions_to_normalize = []

for region in regions:
    morans_normalized_file = f"{normalized_matrices_dir}/{region}_morans_quantile_normalized.csv"
    counts_normalized_file = f"{normalized_matrices_dir}/{region}_counts_quantile_normalized.csv"
    
    if os.path.exists(morans_normalized_file) and os.path.exists(counts_normalized_file):
        print(f"Loading normalized matrices for {region}...")
        all_regions_normalized_moransi[region] = pd.read_csv(morans_normalized_file, index_col=0)
        all_regions_normalized_counts[region] = pd.read_csv(counts_normalized_file, index_col=0)
        print(f"  ✓ {region} loaded from existing files")
    else:
        print(f"  ⚠ Normalized matrices not found for {region}, will compute normalization")
        regions_to_normalize.append(region)

# If any regions need normalization, compute it now
if regions_to_normalize:
    print(f"\nComputing quantile normalization for {len(regions_to_normalize)} region(s)...")
    os.makedirs(normalized_matrices_dir, exist_ok=True)
    
    for region in regions_to_normalize:
        if region in all_regions_morans and region in all_regions_counts:
            print(f"  Normalizing {region}...")
            
            # Get the data
            morans_df = all_regions_morans[region].copy()
            counts_df = all_regions_counts[region].copy()
            
            # Convert N/A strings to NaN for Moran's I
            morans_df = morans_df.replace('N/A', np.nan)
            morans_df = morans_df.replace('', np.nan)
            
            # Convert to numeric (will convert NaN strings to NaN)
            morans_df = morans_df.apply(pd.to_numeric, errors='coerce')
            counts_df = counts_df.apply(pd.to_numeric, errors='coerce')
            
            # Fill NaN counts with 0
            counts_df = counts_df.fillna(0)
            
            # Perform quantile normalization
            morans_normalized = quantileNormalize(morans_df)
            counts_normalized = quantileNormalize(counts_df)
            
            # Store normalized data
            all_regions_normalized_moransi[region] = morans_normalized
            all_regions_normalized_counts[region] = counts_normalized
            
            # Save normalized matrices for future use
            morans_normalized_file = f"{normalized_matrices_dir}/{region}_morans_quantile_normalized.csv"
            counts_normalized_file = f"{normalized_matrices_dir}/{region}_counts_quantile_normalized.csv"
            morans_normalized.to_csv(morans_normalized_file)
            counts_normalized.to_csv(counts_normalized_file)
            print(f"    ✓ {region} normalized and saved")

print(f"\n✓ Normalized matrices ready for {len(all_regions_normalized_moransi)} region(s)")

#%% Skip Boxplot Visualization (not needed for heatmap generation)
# Create boxplots to visualize the effect of quantile normalization
def create_normalization_boxplots(region='CA1'):
    """
    Create boxplots showing data distribution before and after quantile normalization
    for both Moran's I and counts data
    """
    if region not in all_regions_morans or region not in all_regions_normalized_moransi:
        print(f"Region {region} not available for visualization")
        return
    
    # Create figure with 2x2 subplots
    fig, axes = plt.subplots(2, 2, figsize=(20, 16))
    
    # Moran's I - Before normalization
    morans_before = all_regions_morans[region]
    wt_cols, five_x_cols = identify_sample_groups(morans_before)
    
    # Organize columns: 5X first, then WT
    organized_cols = five_x_cols + wt_cols
    morans_before_organized = morans_before[organized_cols]
    
    # Create boxplot for Moran's I before normalization
    sns.boxplot(data=morans_before_organized, ax=axes[0,0])
    axes[0,0].set_title(f"Moran's I Distribution - {region} (Before Quantile Normalization)")
    axes[0,0].set_ylabel("Moran's I")
    axes[0,0].set_xlabel("Samples")
    axes[0,0].tick_params(axis='x', rotation=45)
    
    # Add vertical line to separate 5X and WT samples
    if five_x_cols and wt_cols:
        axes[0,0].axvline(x=len(five_x_cols)-0.5, color='red', linestyle='--', alpha=0.7)
        axes[0,0].text(len(five_x_cols)/2, axes[0,0].get_ylim()[1]*0.9, '5X Samples', 
                       ha='center', va='top', fontsize=12, color='red')
        axes[0,0].text(len(five_x_cols) + len(wt_cols)/2, axes[0,0].get_ylim()[1]*0.9, 'WT Samples', 
                       ha='center', va='top', fontsize=12, color='blue')
    
    # Moran's I - After normalization
    morans_after = all_regions_normalized_moransi[region]
    morans_after_organized = morans_after[organized_cols]
    
    sns.boxplot(data=morans_after_organized, ax=axes[0,1])
    axes[0,1].set_title(f"Moran's I Distribution - {region} (After Quantile Normalization)")
    axes[0,1].set_ylabel("Normalized Moran's I")
    axes[0,1].set_xlabel("Samples")
    axes[0,1].tick_params(axis='x', rotation=45)
    
    # Add vertical line to separate 5X and WT samples
    if five_x_cols and wt_cols:
        axes[0,1].axvline(x=len(five_x_cols)-0.5, color='red', linestyle='--', alpha=0.7)
        axes[0,1].text(len(five_x_cols)/2, axes[0,1].get_ylim()[1]*0.9, '5X Samples', 
                       ha='center', va='top', fontsize=12, color='red')
        axes[0,1].text(len(five_x_cols) + len(wt_cols)/2, axes[0,1].get_ylim()[1]*0.9, 'WT Samples', 
                       ha='center', va='top', fontsize=12, color='blue')
    
    # Counts - Before normalization
    counts_before = all_regions_counts[region]
    counts_before_organized = counts_before[organized_cols]
    
    sns.boxplot(data=counts_before_organized, ax=axes[1,0])
    axes[1,0].set_title(f"Counts Distribution - {region} (Before Quantile Normalization)")
    axes[1,0].set_ylabel("Counts")
    axes[1,0].set_xlabel("Samples")
    axes[1,0].tick_params(axis='x', rotation=45)
    
    # Add vertical line to separate 5X and WT samples
    if five_x_cols and wt_cols:
        axes[1,0].axvline(x=len(five_x_cols)-0.5, color='red', linestyle='--', alpha=0.7)
        axes[1,0].text(len(five_x_cols)/2, axes[1,0].get_ylim()[1]*0.9, '5X Samples', 
                       ha='center', va='top', fontsize=12, color='red')
        axes[1,0].text(len(five_x_cols) + len(wt_cols)/2, axes[1,0].get_ylim()[1]*0.9, 'WT Samples', 
                       ha='center', va='top', fontsize=12, color='blue')
    
    # Counts - After normalization
    counts_after = all_regions_normalized_counts[region]
    counts_after_organized = counts_after[organized_cols]
    
    sns.boxplot(data=counts_after_organized, ax=axes[1,1])
    axes[1,1].set_title(f"Counts Distribution - {region} (After Quantile Normalization)")
    axes[1,1].set_ylabel("Normalized Counts")
    axes[1,1].set_xlabel("Samples")
    axes[1,1].tick_params(axis='x', rotation=45)
    
    # Add vertical line to separate 5X and WT samples
    if five_x_cols and wt_cols:
        axes[1,1].axvline(x=len(five_x_cols)-0.5, color='red', linestyle='--', alpha=0.7)
        axes[1,1].text(len(five_x_cols)/2, axes[1,1].get_ylim()[1]*0.9, '5X Samples', 
                       ha='center', va='top', fontsize=12, color='red')
        axes[1,1].text(len(five_x_cols) + len(wt_cols)/2, axes[1,1].get_ylim()[1]*0.9, 'WT Samples', 
                       ha='center', va='top', fontsize=12, color='blue')
    
    plt.tight_layout()
    
    # Save the plot
    output_dir = "normalization_plots"
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(f'{output_dir}/normalization_comparison_{region}.png', dpi=300, bbox_inches='tight')
    
    plt.show()
    plt.close()

#%% Configuration Flags for Visualization
# Option to create normalization boxplots to visualize the effect of quantile normalization
CREATE_NORMALIZATION_PLOTS = True  # Set to False to skip normalization plots

# Option to generate heatmaps
CREATE_HEATMAPS = True  # Set to False to skip heatmap generation

if CREATE_NORMALIZATION_PLOTS:
    print("\nCreating normalization plots...")
    # Create plots for first region as example (you can modify to plot all regions)
    if regions and regions[0] in all_regions_morans:
        create_normalization_boxplots(region=regions[0])
        print(f"  ✓ Normalization plots created for {regions[0]}")
else:
    print("\nSkipping normalization plots (set CREATE_NORMALIZATION_PLOTS = True to enable)")

#%% Z-Scoring for Heatmap Visualization
# Function to apply z-scoring for heatmap visualization
def apply_zscore_to_rows(df):
    """
    Apply z-scoring to each row (gene) for heatmap visualization
    """
    df_zscore = df.copy()
    for idx in df.index:
        row_data = df.loc[idx]
        # Remove NaN values for z-scoring
        valid_data = row_data.dropna()
        if len(valid_data) > 1:  # Need at least 2 values for z-scoring
            mean_val = valid_data.mean()
            std_val = valid_data.std()
            if std_val > 0:
                df_zscore.loc[idx] = (row_data - mean_val) / std_val
            else:
                df_zscore.loc[idx] = 0
        else:
            df_zscore.loc[idx] = 0
    return df_zscore

# Test z-scoring on a small dataset
if regions and regions[0] in all_regions_normalized_moransi:
    test_data = all_regions_normalized_moransi[regions[0]].iloc[:5, :5]  # First 5 genes, first 5 samples
    test_zscore = apply_zscore_to_rows(test_data)
    print(f"Z-scoring test for {regions[0]}:")
    print(f"  Original data shape: {test_data.shape}")
    print(f"  Z-scored data shape: {test_zscore.shape}")
    print(f"  Z-scored data mean: {test_zscore.mean().mean():.3f}")
    print(f"  Z-scored data std: {test_zscore.std().mean():.3f}")

#%% VERIFICATION:- Clustering function added (only visualization change)
def cluster_genes_hierarchical(data, method='ward', metric='euclidean'):
    """
    Perform hierarchical clustering on rows (genes) and return reordered indices.
    
    Parameters:
    - data: DataFrame with genes as rows, samples as columns (already z-scored/log-normalized)
    - method: linkage method ('ward' by default)
    - metric: distance metric ('euclidean' by default)
    
    Returns:
    - List of gene indices in clustered order
    """
    # Remove any NaN values for clustering
    data_clean = data.fillna(0)
    
    # Calculate pairwise distances between rows (genes)
    distances = pdist(data_clean.values, metric=metric)
    
    # Perform hierarchical clustering
    linkage_matrix = linkage(distances, method=method)
    
    # Get the order of leaves (genes) from the dendrogram
    clustered_order = leaves_list(linkage_matrix)
    
    # Map back to gene names
    clustered_gene_names = data_clean.index[clustered_order].tolist()
    
    return clustered_gene_names

#%% Heatmap Generation
# Function to create heatmaps with z-scoring applied to each row
def heatmap_by_region(region, value, save_png=True, normalize=False, show=False):
    """
    Create heatmap for a specific region and value type
    value: 'Count' or 'Morans_i'
    normalize: if True, use normalized data; if False, use raw data
    """
    # Select the right data based on the value and normalization
    if value == 'Count':
        if normalize:
            data = all_regions_normalized_counts[region]
            normalize_str = 'quantile normalized'
        else:
            data = all_regions_counts[region]
            normalize_str = 'raw'
    elif value == 'Morans_i':
        if normalize:
            data = all_regions_normalized_moransi[region]
            normalize_str = 'quantile normalized'
        else:
            data = all_regions_morans[region]
            normalize_str = 'raw'
    else:
        raise ValueError("value must be 'Count' or 'Morans_i'")
    
    # Organize columns: WT samples first, then 5X samples (if available)
    wt_cols, five_x_cols = identify_sample_groups(data)
    if len(wt_cols) > 0 and len(five_x_cols) > 0:
        organized_cols = wt_cols + five_x_cols
    else:
        # If no WT/5X groups found, use all columns in original order
        print(f"  Note: No WT/5X groups found for {region}, using all columns")
        organized_cols = list(data.columns)
    data_organized = data[organized_cols]
    
    # VERIFICATION: CHANGED - Replaced mean-difference sorting with hierarchical clustering
    # Apply z-scoring first for clustering
    data_zscore_for_clustering = apply_zscore_to_rows(data_organized)
    
    # Cluster genes using hierarchical clustering
    if len(data_organized) > 1:  # Need at least 2 genes to cluster
        clustered_gene_order = cluster_genes_hierarchical(data_zscore_for_clustering, method='ward', metric='euclidean')
        data_organized_sorted = data_organized.loc[clustered_gene_order]
        diff_info = "\nGenes ordered by hierarchical clustering (Ward linkage, Euclidean distance)"
        sorted_indices = clustered_gene_order  # For return value consistency
    else:
        data_organized_sorted = data_organized
        diff_info = ""
        sorted_indices = None
    
    # Apply z-scoring to each row (gene) for visualization
    data_zscore = apply_zscore_to_rows(data_organized_sorted)
    
    # Create the heatmap
    plt.figure(figsize=(30, 30))
    sns.heatmap(data_zscore, annot=False, cmap='coolwarm', linewidths=.5)
    
    # Set title based on whether WT/5X groups exist
    if five_x_cols and wt_cols:
        title_suffix = 'WT Samples | 5X Samples'
    else:
        title_suffix = 'All Samples'
    plt.title(f'{value} Heatmap for {region} ({normalize_str}) - Z-scored by gene\n{title_suffix}{diff_info}')
    # VERIFICATION: CHANGED - Updated y-axis label
    plt.ylabel('Genes (ordered by hierarchical clustering)')
    plt.xlabel('Samples')
    
    # Add a vertical line to separate WT and 5X samples (if groups exist)
    if five_x_cols and wt_cols:
        plt.axvline(x=len(wt_cols), color='black', linewidth=3, alpha=0.8)
        # Add text labels for the two groups (WT left, 5X right)
        plt.text(len(wt_cols)/2, -len(data_zscore)*0.02, 'WT Samples', 
                ha='center', va='top', fontsize=16, color='blue', weight='bold')
        plt.text(len(wt_cols) + len(five_x_cols)/2, -len(data_zscore)*0.02, '5X Samples', 
                ha='center', va='top', fontsize=16, color='red', weight='bold')
    
    plt.tight_layout()

    # Save as PNG if needed
    if save_png:
        # VERIFICATION: CHANGED - Output directory name changed to avoid overwriting originals
        output_dir = f"{base_path}/heatmaps_zscore_clustered"
        os.makedirs(output_dir, exist_ok=True)
        plt.savefig(f'{output_dir}/heatmap_{value}_{region}_{normalize_str}_zscored.png', dpi=300, bbox_inches='tight')
        
        # Also save as TIFF in tiff_format subfolder
        tiff_dir = f"{output_dir}/tiff_format"
        os.makedirs(tiff_dir, exist_ok=True)
        plt.savefig(f'{tiff_dir}/heatmap_{value}_{region}_{normalize_str}_zscored.tiff', dpi=300, bbox_inches='tight')

    # Show plot if required
    if show:
        plt.show()
    plt.close()
    
    # Return the sorted data and ordering for consistency across heatmaps
    if sorted_indices is not None:
        return data_organized_sorted, sorted_indices
    else:
        return data_organized_sorted, None

#%% Heatmap Generation
# Generate heatmaps for all regions
if CREATE_HEATMAPS:
    print("\nGenerating heatmaps...")
    for region in regions:
        if region in all_regions_morans:
            print(f"Creating heatmaps for {region}")
            
            # First create Moran's I heatmap to get the gene ordering
            morans_data, gene_ordering = heatmap_by_region(region, 'Morans_i', normalize=True)
            
            # Create counts heatmap with the same gene ordering
            if gene_ordering is not None:
                # Reorder counts data to match Moran's I ordering
                counts_data = all_regions_normalized_counts[region]
                wt_cols, five_x_cols = identify_sample_groups(counts_data)
                organized_cols = wt_cols + five_x_cols
                counts_organized = counts_data[organized_cols]
                counts_organized_sorted = counts_organized.loc[gene_ordering]
                
                # Apply z-scoring and create counts heatmap with same ordering
                counts_zscore = apply_zscore_to_rows(counts_organized_sorted)
                
                # Create the counts heatmap with same gene order
                plt.figure(figsize=(30, 30))
                sns.heatmap(counts_zscore, annot=False, cmap='coolwarm', linewidths=.5)
                
                # Set title based on whether WT/5X groups exist
                if five_x_cols and wt_cols:
                    title_suffix = 'WT Samples | 5X Samples'
                else:
                    title_suffix = 'All Samples'
                plt.title(f'Counts Heatmap for {region} (quantile normalized) - Z-scored by gene\n{title_suffix}\nSame gene ordering as Moran\'s I')
                # VERIFICATION: CHANGED - Updated y-axis label
                plt.ylabel('Genes (ordered by hierarchical clustering)')
                plt.xlabel('Samples')
                
                # Add vertical line to separate WT and 5X samples (if groups exist)
                if five_x_cols and wt_cols:
                    plt.axvline(x=len(wt_cols), color='black', linewidth=3, alpha=0.8)
                    plt.text(len(wt_cols)/2, -len(counts_zscore)*0.02, 'WT Samples', 
                            ha='center', va='top', fontsize=16, color='blue', weight='bold')
                    plt.text(len(wt_cols) + len(five_x_cols)/2, -len(counts_zscore)*0.02, '5X Samples', 
                            ha='center', va='top', fontsize=16, color='red', weight='bold')
                
                plt.tight_layout()
                
                # Save counts heatmap
                # VERIFICATION: CHANGED - Output directory name changed to avoid overwriting originals
                output_dir = f"{base_path}/heatmaps_zscore_clustered"
                plt.savefig(f'{output_dir}/heatmap_Counts_{region}_quantile normalized_zscored_ordered.png', dpi=300, bbox_inches='tight')
                
                # Also save as TIFF in tiff_format subfolder
                tiff_dir = f"{output_dir}/tiff_format"
                os.makedirs(tiff_dir, exist_ok=True)
                plt.savefig(f'{tiff_dir}/heatmap_Counts_{region}_quantile normalized_zscored_ordered.tiff', dpi=300, bbox_inches='tight')
                plt.close()
                
                print(f"  Created ordered heatmaps for {region}")
            else:
                # Fallback to original method if no gene ordering available
                heatmap_by_region(region, 'Count', normalize=True)
            
            # Create raw data heatmaps (not normalized)
            heatmap_by_region(region, 'Morans_i', normalize=False)
            heatmap_by_region(region, 'Count', normalize=False)

    print("Heatmap generation completed!")
else:
    print("\nSkipping heatmap generation (set CREATE_HEATMAPS = True to enable)")

#%% Statistical Analysis - Permutation Tests
if SKIP_STATISTICAL_ANALYSIS:
    print("\n" + "="*60)
    print("SKIPPING STATISTICAL ANALYSIS")
    print("="*60)
    print("Statistical results already exist in:")
    print(f"  {base_path}/statistical_results/")
    print("Only heatmap visualization with clustering is being regenerated.")
    print("="*60 + "\n")
else:
    # Function to run permutation tests (keep the same as original)
    def permutation_test(df, wt_cols, five_x_cols, num_permutations=NUM_PERMUTATIONS):
        wt_mean = df[wt_cols].mean(axis=1)
        five_x_mean = df[five_x_cols].mean(axis=1)
        observed_mean_diff = wt_mean - five_x_mean
        n_genes = df.shape[0]

        # Store results for each gene
        p_values_direct = np.zeros(n_genes)
        p_values_fit = np.zeros(n_genes)

        for i in range(n_genes):
            gene_data = pd.concat([df[wt_cols].iloc[i], df[five_x_cols].iloc[i]])
            
            # Skip genes with all NaN values
            if gene_data.isna().all():
                p_values_direct[i] = np.nan
                p_values_fit[i] = np.nan
                continue
            
            # Remove NaN values for permutation test
            gene_data_clean = gene_data.dropna()
            
            # Skip if not enough valid values for comparison
            if len(gene_data_clean) < 2:
                p_values_direct[i] = np.nan
                p_values_fit[i] = np.nan
                continue

            # Generate permutations
            perm_mean_diffs = np.zeros(num_permutations)
            for perm in range(num_permutations):
                # Shuffle labels and divide into two groups (WT and 5x)
                shuffled = np.random.permutation(gene_data_clean.values)
                wt_perm = shuffled[:len(wt_cols)]
                five_x_perm = shuffled[len(wt_cols):]

                # Compute the difference in means for this permutation
                perm_mean_diffs[perm] = wt_perm.mean() - five_x_perm.mean()

            # Compute the direct p-value: proportion of permuted diffs that are more extreme than observed
            observed_diff = observed_mean_diff.iloc[i]
            if pd.isna(observed_diff):
                p_values_direct[i] = np.nan
                p_values_fit[i] = np.nan
                continue
                
            more_extreme = np.sum(np.abs(perm_mean_diffs) >= np.abs(observed_diff))
            p_values_direct[i] = more_extreme  / num_permutations

            # Fit a normal distribution to the permuted mean differences
            # Remove any NaN or inf values before fitting
            perm_mean_diffs_clean = perm_mean_diffs[np.isfinite(perm_mean_diffs)]
            
            if len(perm_mean_diffs_clean) < 2:
                p_values_fit[i] = np.nan
                continue
                
            mu, std = norm.fit(perm_mean_diffs_clean)
            
            if std == 0 or not np.isfinite(mu) or not np.isfinite(std):
                p_values_fit[i] = np.nan
                continue

            # Compute z-score for observed mean difference
            z_score = (observed_diff - mu) / std

            # Calculate p-value based on normal distribution fit (two-tailed test)
            p_values_fit[i] = 2 * (1 - norm.cdf(np.abs(z_score)))

        # Add p-values to the dataframe
        df['wt_mean'] = wt_mean
        df['five_x_mean'] = five_x_mean
        df['mean_diff'] = observed_mean_diff
        df['permutation_p_value_direct'] = p_values_direct
        df['permutation_p_value_fit'] = p_values_fit

        return df

    # Function to run permutation tests for all regions
    def run_permutation_test(region_data, num_permutations=NUM_PERMUTATIONS):
        wt_cols, five_x_cols = identify_sample_groups(region_data)
        
        if len(wt_cols) == 0 or len(five_x_cols) == 0:
            print(f"  ⚠ Warning: No WT or 5X columns found for statistical analysis")
            print(f"     Found columns: {list(region_data.columns[:5])}...")
            print(f"     Statistical tests require WT and 5X sample groups")
            print(f"     Skipping permutation tests for this region")
            # Return dataframe with NaN values for statistical columns
            region_data = region_data.copy()
            region_data['wt_mean'] = np.nan
            region_data['five_x_mean'] = np.nan
            region_data['mean_diff'] = np.nan
            region_data['permutation_p_value_direct'] = np.nan
            region_data['permutation_p_value_fit'] = np.nan
            return region_data
        else:
            print(f"  Running permutation test with {len(wt_cols)} WT and {len(five_x_cols)} 5X samples")
            permutation_test_df = permutation_test(region_data, wt_cols, five_x_cols, num_permutations = num_permutations)
            return permutation_test_df

    # Run permutation tests for all regions
    print("\nRunning permutation tests...")
    for region in regions:
        if region in all_regions_normalized_moransi:
            print(f"Running permutation test for Moran's I in {region}")
            all_regions_normalized_moransi[region] = run_permutation_test(all_regions_normalized_moransi[region])
            
            print(f"Running permutation test for counts in {region}")
            all_regions_normalized_counts[region] = run_permutation_test(all_regions_normalized_counts[region])
    
    print("Permutation tests completed!")

#%% FDR Correction
if SKIP_STATISTICAL_ANALYSIS:
    pass  # Skip FDR correction
else:
    # Apply FDR correction with rate 0.05
    print("\nApplying FDR correction (rate = 0.05)...")

def apply_fdr_correction(df, p_value_col='permutation_p_value_fit', alpha=0.05):
    """
    Apply FDR correction using Benjamini-Hochberg method on the fit p-values
    
    Args:
        df: DataFrame containing p-values
        p_value_col: Column name containing fit p-values (default: 'permutation_p_value_fit')
        alpha: FDR rate (default 0.05)
    
    Returns:
        DataFrame with FDR-corrected p-values added as new column
    """
    if p_value_col not in df.columns:
        print(f"  Warning: Column '{p_value_col}' not found, skipping FDR correction")
        return df
    
    # Get valid fit p-values (drop NaN)
    valid_mask = df[p_value_col].notna()
    valid_pvals = df.loc[valid_mask, p_value_col]
    
    if len(valid_pvals) == 0:
        print(f"  Warning: No valid fit p-values found, skipping FDR correction")
        df['fdr_corrected_fit'] = np.nan
        df['significant_fit'] = False
        return df
    
    print(f"  Processing {len(valid_pvals)} valid fit p-values...")
    
    # Apply FDR correction to fit p-values
    try:
        # multipletests returns (rejected, pvals_corrected, alphacSidak, alphacBonf)
        _, pvals_corrected, _, _ = multipletests(valid_pvals, method='fdr_bh', alpha=alpha)
        
        # Add corrected fit p-values to dataframe as new column
        df.loc[valid_mask, 'fdr_corrected_fit'] = pvals_corrected
        
        # Add significance flag for fit p-values (corrected p < alpha)
        df['significant_fit'] = df['fdr_corrected_fit'] < alpha
        
        # Count significant results
        significant_count = df['significant_fit'].sum()
        print(f"  ✓ FDR correction completed: {significant_count} significant results at α={alpha}")
        
    except Exception as e:
        print(f"  Error in FDR correction: {e}")
        df['fdr_corrected_fit'] = np.nan
        df['significant_fit'] = False
    
    return df

# Apply FDR correction to Moran's I data
print("\nMoran's I data:")
for region in regions:
    if region in all_regions_normalized_moransi:
        print(f"  {region}:")
        all_regions_normalized_moransi[region] = apply_fdr_correction(
            all_regions_normalized_moransi[region], 
            alpha=0.05
        )

# Apply FDR correction to Counts data
print("\nCounts data:")
for region in regions:
    if region in all_regions_normalized_counts:
        print(f"  {region}:")
        all_regions_normalized_counts[region] = apply_fdr_correction(
            all_regions_normalized_counts[region], 
            alpha=0.05
        )

        print("\n✓ FDR correction completed for all regions!")

#%% Save Statistical Results to CSV
if SKIP_STATISTICAL_ANALYSIS:
    pass  # Skip saving statistical results
else:
    # Save all statistical results to CSV files for further analysis
    print("\n=== SAVING STATISTICAL RESULTS TO CSV ===")

def save_statistical_results_to_csv(all_regions_normalized_moransi, all_regions_normalized_counts, output_dir=f"{base_path}/statistical_results"):
    """
    Save all statistical results to CSV files for further analysis
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}")
    
    # Save Moran's I results
    print("\nSaving Moran's I statistical results...")
    for region in all_regions_normalized_moransi.keys():
        df = all_regions_normalized_moransi[region]
        
        # Create filename
        filename = f"{region}_morans_statistical_results.csv"
        filepath = os.path.join(output_dir, filename)
        
        # Save to CSV
        df.to_csv(filepath)
        print(f"  ✓ {region}: {filename} ({df.shape[0]} genes, {df.shape[1]} columns)")
    
    # Save Counts results
    print("\nSaving Counts statistical results...")
    for region in all_regions_normalized_counts.keys():
        df = all_regions_normalized_counts[region]
        
        # Create filename
        filename = f"{region}_counts_statistical_results.csv"
        filepath = os.path.join(output_dir, filename)
        
        # Save to CSV
        df.to_csv(filepath)
        print(f"  ✓ {region}: {filename} ({df.shape[0]} genes, {df.shape[1]} columns)")
    
    # Create summary file
    print("\nCreating summary file...")
    summary_data = []
    
    for region in all_regions_normalized_moransi.keys():
        if region in all_regions_normalized_counts:
            morans_df = all_regions_normalized_moransi[region]
            counts_df = all_regions_normalized_counts[region]
            
            # Count significant genes (if FDR correction was applied)
            morans_sig = 0
            counts_sig = 0
            
            if 'significant_fit' in morans_df.columns:
                morans_sig = morans_df['significant_fit'].sum()
            
            if 'significant_fit' in counts_df.columns:
                counts_sig = counts_df['significant_fit'].sum()
            
            # Get total genes
            total_genes = len(morans_df)
            
            summary_data.append({
                'Region': region,
                'Total_Genes': total_genes,
                'Morans_Significant_Genes': morans_sig,
                'Counts_Significant_Genes': counts_sig,
                'Morans_Significance_Rate': f"{(morans_sig/total_genes)*100:.2f}%" if total_genes > 0 else "N/A",
                'Counts_Significance_Rate': f"{(counts_sig/total_genes)*100:.2f}%" if total_genes > 0 else "N/A"
            })
    
    # Save summary
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_filepath = os.path.join(output_dir, "statistical_analysis_summary.csv")
        summary_df.to_csv(summary_filepath, index=False)
        print(f"  ✓ Summary saved: statistical_analysis_summary.csv")
    
    # Create comprehensive summary with gene names
    print("\nCreating comprehensive summary with gene names...")
    comprehensive_summary = []
    
    for region in all_regions_normalized_moransi.keys():
        if region in all_regions_normalized_counts:
            morans_df = all_regions_normalized_moransi[region]
            counts_df = all_regions_normalized_counts[region]
            
            # Get significant gene names
            morans_sig_genes = []
            counts_sig_genes = []
            
            if 'significant_fit' in morans_df.columns:
                morans_sig_genes = morans_df[morans_df['significant_fit']].index.tolist()
            
            if 'significant_fit' in counts_df.columns:
                counts_sig_genes = counts_df[counts_df['significant_fit']].index.tolist()
            
            # Analyze direction for significant genes
            morans_direction_summary = []
            counts_direction_summary = []
            
            # Analyze Moran's I direction
            for gene in morans_sig_genes:
                if gene in morans_df.index:
                    mean_diff = morans_df.loc[gene, 'mean_diff']
                    if mean_diff > 0:
                        direction = "5X > WT"
                    elif mean_diff < 0:
                        direction = "WT > 5X"
                    else:
                        direction = "No diff"
                    morans_direction_summary.append(f"{gene}({direction})")
            
            # Analyze Counts direction
            for gene in counts_sig_genes:
                if gene in counts_df.index:
                    mean_diff = counts_df.loc[gene, 'mean_diff']
                    if mean_diff > 0:
                        direction = "5X > WT"
                    elif mean_diff < 0:
                        direction = "WT > 5X"
                    else:
                        direction = "No diff"
                    counts_direction_summary.append(f"{gene}({direction})")
            
            # Create comprehensive entry with direction information
            comprehensive_summary.append({
                'Region': region,
                'Total_Genes': len(morans_df),
                'Morans_Significant_Count': len(morans_sig_genes),
                'Counts_Significant_Count': len(counts_sig_genes),
                'Morans_Significant_Genes': '; '.join(morans_sig_genes) if morans_sig_genes else 'None',
                'Counts_Significant_Genes': '; '.join(counts_sig_genes) if counts_sig_genes else 'None',
                'Morans_Direction_Summary': '; '.join(morans_direction_summary) if morans_direction_summary else 'None',
                'Counts_Direction_Summary': '; '.join(counts_direction_summary) if counts_direction_summary else 'None'
            })
    
    # Save comprehensive summary
    if comprehensive_summary:
        comp_df = pd.DataFrame(comprehensive_summary)
        comp_filepath = os.path.join(output_dir, "comprehensive_summary_with_gene_names.csv")
        comp_df.to_csv(comp_filepath, index=False)
        print(f"  ✓ Comprehensive summary saved: comprehensive_summary_with_gene_names.csv")
    
    print(f"\n✓ All statistical results saved to: {output_dir}/")
    print("Files created:")
    print("  - {region}_morans_statistical_results.csv (for each region)")
    print("  - {region}_counts_statistical_results.csv (for each region)")
    print("  - statistical_analysis_summary.csv (overall summary)")
    print("  - comprehensive_summary_with_gene_names.csv (summary with gene names + direction)")
    print("  - direction_analysis_results.csv (detailed direction analysis)")

    # Save the results
    save_statistical_results_to_csv(all_regions_normalized_moransi, all_regions_normalized_counts)

#%% Significant Genes Visualization
if SKIP_STATISTICAL_ANALYSIS:
    pass  # Skip visualization
else:
    # Create visualizations to show significant genes and p-value thresholds
    print("\n=== CREATING SIGNIFICANT GENES VISUALIZATIONS ===")

def create_significant_genes_plots(all_regions_normalized_moransi, all_regions_normalized_counts, max_regions=4):
    """
    Create plots showing significant genes with p-value thresholds
    """
    if not all_regions_normalized_moransi or not all_regions_normalized_counts:
        print("No data available for visualization")
        return
    
    # Limit to first few regions to avoid too many plots
    regions_to_plot = list(all_regions_normalized_moransi.keys())[:max_regions]
    
    fig, axes = plt.subplots(2, 2, figsize=(20, 16))
    axes = axes.flatten()
    
    for i, region in enumerate(regions_to_plot):
        if i >= len(axes):
            break
            
        if region not in all_regions_normalized_moransi or region not in all_regions_normalized_counts:
            continue
        
        morans_df = all_regions_normalized_moransi[region]
        counts_df = all_regions_normalized_counts[region]
        
        # Check if FDR correction was applied
        if 'fdr_corrected_fit' not in morans_df.columns or 'fdr_corrected_fit' not in counts_df.columns:
            print(f"  Warning: FDR correction not found for {region}, skipping visualization")
            continue
        
        # Get p-values
        morans_pvals = morans_df['fdr_corrected_fit']
        counts_pvals = counts_df['fdr_corrected_fit']
        
        # Remove NaN values for plotting
        morans_valid = morans_pvals.dropna()
        counts_valid = counts_pvals.dropna()
        
        if len(morans_valid) == 0 or len(counts_valid) == 0:
            print(f"  Warning: No valid p-values found for {region}")
            continue
        
        # Create subplot
        ax = axes[i]
        
        # Plot Moran's I p-values
        if len(morans_valid) > 0:
            ax.scatter(range(len(morans_valid)), -np.log10(morans_valid), 
                      alpha=0.6, s=20, label='Moran\'s I', color='blue')
        
        # Plot Counts p-values
        if len(counts_valid) > 0:
            ax.scatter(range(len(counts_valid)), -np.log10(counts_valid), 
                      alpha=0.6, s=20, label='Counts', color='red')
        
        # Add significance threshold line (α = 0.05)
        threshold = -np.log10(0.05)
        ax.axhline(y=threshold, color='red', linestyle='--', alpha=0.8, 
                  label=f'Significance threshold (α=0.05)')
        
        # Add FDR threshold line (α = 0.05)
        ax.axhline(y=threshold, color='orange', linestyle=':', alpha=0.8, 
                  label=f'FDR threshold (α=0.05)')
        
        # Customize plot
        ax.set_title(f'{region} - Significant Genes (FDR Corrected P-values)')
        ax.set_xlabel('Genes (ranked by p-value)')
        ax.set_ylabel('-log10(FDR corrected p-value)')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Add text showing counts
        morans_sig = morans_df.get('significant_fit', pd.Series([False])).sum()
        counts_sig = counts_df.get('significant_fit', pd.Series([False])).sum()
        
        ax.text(0.02, 0.98, f'Moran\'s I significant: {morans_sig}\nCounts significant: {counts_sig}', 
                transform=ax.transAxes, verticalalignment='top',
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    # Hide unused subplots
    for i in range(len(regions_to_plot), len(axes)):
        axes[i].set_visible(False)
    
    plt.tight_layout()
    
    # Save plot as PNG
    output_dir = "statistical_plots"
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(f'{output_dir}/significant_genes_plots_{region}.png', dpi=300, bbox_inches='tight')
    print(f"  ✓ Plot saved: significant_genes_plots_{region}.png")
    
    plt.show()
    plt.close()

def create_volcano_plots(all_regions_normalized_moransi, all_regions_normalized_counts, max_regions=4):
    """
    Create volcano plots showing effect size vs significance
    """
    if not all_regions_normalized_moransi or not all_regions_normalized_counts:
        print("No data available for volcano plots")
        return
    
    # Limit to first few regions
    regions_to_plot = list(all_regions_normalized_moransi.keys())[:max_regions]
    
    fig, axes = plt.subplots(2, 2, figsize=(20, 16))
    axes = axes.flatten()
    
    for i, region in enumerate(regions_to_plot):
        if i >= len(axes):
            break
            
        if region not in all_regions_normalized_moransi or region not in all_regions_normalized_counts:
            continue
        
        morans_df = all_regions_normalized_moransi[region]
        counts_df = all_regions_normalized_counts[region]
        
        # Check if we have the necessary columns
        if 'fdr_corrected_fit' not in morans_df.columns:
            continue
        
        # Get effect size (mean difference between 5X and WT)
        wt_cols, five_x_cols = identify_sample_groups(morans_df)
        
        if len(wt_cols) > 0 and len(five_x_cols) > 0:
            # Calculate effect size for Moran's I
            five_x_mean = morans_df[five_x_cols].mean(axis=1)
            wt_mean = morans_df[wt_cols].mean(axis=1)
            effect_size = five_x_mean - wt_mean
            
            # Get p-values
            pvals = morans_df['fdr_corrected_fit']
            
            # Remove NaN values
            valid_mask = pvals.notna() & effect_size.notna()
            valid_pvals = pvals[valid_mask]
            valid_effect = effect_size[valid_mask]
            
            if len(valid_pvals) > 0:
                ax = axes[i]
                
                # Create volcano plot
                ax.scatter(valid_effect, -np.log10(valid_pvals), alpha=0.6, s=20)
                
                # Add significance threshold
                threshold = -np.log10(0.05)
                ax.axhline(y=threshold, color='red', linestyle='--', alpha=0.8)
                
                # Add effect size threshold lines
                ax.axvline(x=0, color='gray', linestyle='-', alpha=0.5)
                
                # Customize plot
                ax.set_title(f'{region} - Volcano Plot (Moran\'s I)')
                ax.set_xlabel('Effect Size (5X - WT)')
                ax.set_ylabel('-log10(FDR corrected p-value)')
                ax.grid(True, alpha=0.3)
                
                # Count significant genes
                sig_count = (valid_pvals < 0.05).sum()
                ax.text(0.02, 0.98, f'Significant genes: {sig_count}', 
                       transform=ax.transAxes, verticalalignment='top',
                       bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    # Hide unused subplots
    for i in range(len(regions_to_plot), len(axes)):
        axes[i].set_visible(False)
    
    plt.tight_layout()
    
    # Save plot as PNG
    output_dir = f"{base_path}/statistical_plots"
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(f'{output_dir}/volcano_plots_{region}.png', dpi=300, bbox_inches='tight')
    print(f"  ✓ Plot saved: volcano_plots_{region}.png")
    
    plt.show()
    plt.close()

    # Create visualizations
    print("Creating significant genes plots...")
    create_significant_genes_plots(all_regions_normalized_moransi, all_regions_normalized_counts)

    print("Creating volcano plots...")
    create_volcano_plots(all_regions_normalized_moransi, all_regions_normalized_counts)

#%% Moran's I Distribution Histograms per Region
if SKIP_STATISTICAL_ANALYSIS:
    pass  # Skip histograms
else:
    # Create histograms to visualize Moran's I distribution for each region
    def create_morans_histograms():
        """
        Create histograms showing Moran's I distribution for each region
        """
        print("\nCreating Moran's I distribution histograms...")
        
        # Calculate how many regions we have to determine subplot layout
        n_regions = len(regions)
        n_cols = 3  # 3 columns per row
        n_rows = (n_regions + n_cols - 1) // n_cols  # Ceiling division
        
        # Create figure for before normalization
        fig_before, axes_before = plt.subplots(n_rows, n_cols, figsize=(20, 6*n_rows))
        if n_rows == 1:
            axes_before = axes_before.reshape(1, -1)
        
        # Create figure for after normalization
        fig_after, axes_after = plt.subplots(n_rows, n_cols, figsize=(20, 6*n_rows))
        if n_rows == 1:
            axes_after = axes_after.reshape(1, -1)
        
        # Flatten axes for easier iteration
        axes_before_flat = axes_before.flatten()
        axes_after_flat = axes_after.flatten()
        
        for i, region in enumerate(regions):
            if region in all_regions_morans and region in all_regions_normalized_moransi:
                # Get data for this region
                morans_before = all_regions_morans[region]
                morans_after = all_regions_normalized_moransi[region]
                
                # Flatten the data for histogram
                morans_before_flat = morans_before.values.flatten()
                morans_after_flat = morans_after.values.flatten()
                
                # Convert to numeric and handle non-numeric values
                try:
                    morans_before_flat = pd.to_numeric(morans_before_flat, errors='coerce')
                    morans_after_flat = pd.to_numeric(morans_after_flat, errors='coerce')
                except Exception as e:
                    print(f"  Warning: Could not convert data to numeric for {region}: {e}")
                    continue
                
                # Remove NaN values
                morans_before_flat = morans_before_flat[~np.isnan(morans_before_flat)]
                morans_after_flat = morans_after_flat[~np.isnan(morans_after_flat)]
                
                # Check if we have valid data after cleaning
                if len(morans_before_flat) == 0 or len(morans_after_flat) == 0:
                    print(f"  Warning: No valid numeric data found for {region} after cleaning")
                    continue
                
                # Create histogram for before normalization
                axes_before_flat[i].hist(morans_before_flat, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
                axes_before_flat[i].set_title(f'{region} - Before Normalization')
                axes_before_flat[i].set_xlabel('Moran\'s I')
                axes_before_flat[i].set_ylabel('Frequency')
                axes_before_flat[i].grid(True, alpha=0.3)
                
                # Add statistics
                mean_val = np.mean(morans_before_flat)
                std_val = np.std(morans_before_flat)
                axes_before_flat[i].axvline(mean_val, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_val:.3f}')
                axes_before_flat[i].legend()
                axes_before_flat[i].text(0.02, 0.98, f'Mean: {mean_val:.3f}\nStd: {std_val:.3f}\nN: {len(morans_before_flat)}', 
                                        transform=axes_before_flat[i].transAxes, verticalalignment='top',
                                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                
                # Create histogram for after normalization
                axes_after_flat[i].hist(morans_after_flat, bins=50, alpha=0.7, color='lightgreen', edgecolor='black')
                axes_after_flat[i].set_title(f'{region} - After Quantile Normalization')
                axes_after_flat[i].set_xlabel('Normalized Moran\'s I')
                axes_after_flat[i].set_ylabel('Frequency')
                axes_after_flat[i].grid(True, alpha=0.3)
                
                # Add statistics
                mean_val_norm = np.mean(morans_after_flat)
                std_val_norm = np.std(morans_after_flat)
                axes_after_flat[i].axvline(mean_val_norm, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_val_norm:.3f}')
                axes_after_flat[i].legend()
                axes_after_flat[i].text(0.02, 0.98, f'Mean: {mean_val_norm:.3f}\nStd: {std_val_norm:.3f}\nN: {len(morans_after_flat)}', 
                                        transform=axes_after_flat[i].transAxes, verticalalignment='top',
                                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                
                print(f"  Created histograms for {region}")
            else:
                # Hide unused subplots
                axes_before_flat[i].set_visible(False)
                axes_after_flat[i].set_visible(False)
        
        # Hide any remaining unused subplots
        for i in range(n_regions, len(axes_before_flat)):
            axes_before_flat[i].set_visible(False)
            axes_after_flat[i].set_visible(False)
        
        # Adjust layout and save
        fig_before.suptitle('Moran\'s I Distribution Before Quantile Normalization', fontsize=16, y=0.98)
        fig_before.tight_layout()
        
        fig_after.suptitle('Moran\'s I Distribution After Quantile Normalization', fontsize=16, y=0.98)
        fig_after.tight_layout()
        
        # Display plots (not saving to disk)
        print("Displaying Moran's I distribution histograms...")
        plt.show()
        
        # Close figures to free memory
        plt.close(fig_before)
        plt.close(fig_after)

    # Create the histograms
    create_morans_histograms()

#%% Enhanced Statistical Analysis and Direction Analysis
if SKIP_STATISTICAL_ANALYSIS:
    pass  # Skip enhanced analysis
else:
    # Analyze direction of effect and unique gene patterns for significant genes
    print("\n=== ENHANCED STATISTICAL ANALYSIS ===")

def analyze_direction_and_uniqueness(all_regions_normalized_moransi, all_regions_normalized_counts):
    """
    Analyze the direction of effect (WT vs 5X) and identify unique patterns
    """
    print("Analyzing direction of effect and unique gene patterns...")
    
    direction_analysis = []
    
    for region in all_regions_normalized_moransi.keys():
        if region in all_regions_normalized_counts:
            print(f"\nAnalyzing {region}...")
            
            morans_df = all_regions_normalized_moransi[region]
            counts_df = all_regions_normalized_counts[region]
            
            # Get significant genes (using the FDR-corrected results we already calculated)
            morans_sig_genes = []
            counts_sig_genes = []
            
            if 'significant_fit' in morans_df.columns:
                morans_sig_genes = morans_df[morans_df['significant_fit']].index.tolist()
            
            if 'significant_fit' in counts_df.columns:
                counts_sig_genes = counts_df[counts_df['significant_fit']].index.tolist()
            
            # Analyze direction for Moran's I significant genes
            morans_direction_analysis = []
            for gene in morans_sig_genes:
                if gene in morans_df.index:
                    wt_mean = morans_df.loc[gene, 'wt_mean']
                    five_x_mean = morans_df.loc[gene, 'five_x_mean']
                    mean_diff = morans_df.loc[gene, 'mean_diff']
                    
                    # Determine direction
                    if mean_diff > 0:
                        direction = "WT > 5X (higher spatial autocorrelation in WT)"
                    elif mean_diff < 0:
                        direction = "5X > WT (higher spatial autocorrelation in 5X)"
                    else:
                        direction = "No difference"
                    
                    morans_direction_analysis.append({
                        'Gene': gene,
                        'Region': region,
                        'Analysis_Type': 'Moran\'s I',
                        'WT_Mean': wt_mean,
                        '5X_Mean': five_x_mean,
                        'Mean_Difference': mean_diff,
                        'Direction': direction,
                        'FDR_Corrected_P': morans_df.loc[gene, 'fdr_corrected_fit']
                    })
            
            # Analyze direction for Counts significant genes
            counts_direction_analysis = []
            for gene in counts_sig_genes:
                if gene in counts_df.index:
                    wt_mean = counts_df.loc[gene, 'wt_mean']
                    five_x_mean = counts_df.loc[gene, 'five_x_mean']
                    mean_diff = counts_df.loc[gene, 'mean_diff']
                    
                    # Determine direction
                    if mean_diff > 0:
                        direction = "WT > 5X (higher expression in WT)"
                    elif mean_diff < 0:
                        direction = "5X > WT (higher expression in 5X)"
                    else:
                        direction = "No difference"
                    
                    counts_direction_analysis.append({
                        'Gene': gene,
                        'Region': region,
                        'Analysis_Type': 'Gene Counts',
                        'WT_Mean': wt_mean,
                        '5X_Mean': five_x_mean,
                        'Mean_Difference': mean_diff,
                        'Direction': direction,
                        'FDR_Corrected_P': counts_df.loc[gene, 'fdr_corrected_fit']
                    })
            
            # Identify unique patterns
            morans_only = set(morans_sig_genes) - set(counts_sig_genes)
            counts_only = set(counts_sig_genes) - set(morans_sig_genes)
            both_significant = set(morans_sig_genes) & set(counts_sig_genes)
            
            print(f"  Moran's I only significant: {len(morans_only)} genes")
            print(f"  Counts only significant: {len(counts_only)} genes")
            print(f"  Both significant: {len(both_significant)} genes")
            
            # Add to overall analysis
            direction_analysis.extend(morans_direction_analysis)
            direction_analysis.extend(counts_direction_analysis)
    
    return direction_analysis

# Run the enhanced analysis
print("Running enhanced statistical analysis...")
direction_analysis_results = analyze_direction_and_uniqueness(all_regions_normalized_moransi, all_regions_normalized_counts)

print(f"✓ Enhanced analysis completed! Analyzed {len(direction_analysis_results)} significant gene patterns")

#%% Results Summary and Export
# Save results
if len(direction_analysis_results) > 0:
    output_dir = f"{base_path}/analysis_results"
    os.makedirs(output_dir, exist_ok=True)
    
    # Convert to DataFrame and save
    direction_df = pd.DataFrame(direction_analysis_results)
    direction_df.to_csv(f"{output_dir}/direction_analysis_results.csv", index=False)
    print(f"\nResults saved to {output_dir}/direction_analysis_results.csv")
    
    # Summary statistics
    print(f"\nTotal significant gene patterns analyzed: {len(direction_df)}")
    
    if 'Analysis_Type' in direction_df.columns:
        print("\nBreakdown by analysis type:")
        print(direction_df['Analysis_Type'].value_counts())
    
    if 'Region' in direction_df.columns:
        print("\nBreakdown by region:")
        print(direction_df['Region'].value_counts())
        
    if 'Direction' in direction_df.columns:
        print("\nBreakdown by direction:")
        print(direction_df['Direction'].value_counts())
    else:
        print("No significant genes found for direction analysis")

#%% Add Alternative FDR Corrected P-values Column
if SKIP_STATISTICAL_ANALYSIS:
    pass  # Skip alternative FDR
else:
    # Add a new column with FDR correction using a different threshold
    print("\n=== ADDING ALTERNATIVE FDR CORRECTED P-VALUES COLUMN ===")

    # Set the alternative FDR threshold (you can change this value)
    alternative_fdr_threshold = 0.1  
    print(f"Adding FDR correction with threshold: {alternative_fdr_threshold}")

    def add_alternative_fdr_column(df, p_value_col='permutation_p_value_fit', alpha=0.1):
        """
        Add a new column with FDR correction using alternative threshold
        """
        if p_value_col not in df.columns:
            print(f"  Warning: Column '{p_value_col}' not found, skipping FDR correction")
            return df
        
        # Get valid fit p-values (drop NaN)
        valid_mask = df[p_value_col].notna()
        valid_pvals = df.loc[valid_mask, p_value_col]
        
        if len(valid_pvals) == 0:
            print(f"  Warning: No valid fit p-values found, skipping FDR correction")
            df[f'fdr_corrected_fit_{alpha}'] = np.nan
            df[f'significant_fit_{alpha}'] = False
            return df
        
        # Apply FDR correction to fit p-values
        try:
            # multipletests returns (rejected, pvals_corrected, alphacSidak, alphacBonf)
            _, pvals_corrected, _, _ = multipletests(valid_pvals, method='fdr_bh', alpha=alpha)
            
            # Add corrected fit p-values to dataframe as new column
            df.loc[valid_mask, f'fdr_corrected_fit_{alpha}'] = pvals_corrected
            
            # Add significance flag for fit p-values (corrected p < alpha)
            df[f'significant_fit_{alpha}'] = df[f'fdr_corrected_fit_{alpha}'] < alpha
            
            # Count significant results
            significant_count = df[f'significant_fit_{alpha}'].sum()
            print(f"  ✓ FDR correction completed: {significant_count} significant results at α={alpha}")
            
        except Exception as e:
            print(f"  Error in FDR correction: {e}")
            df[f'fdr_corrected_fit_{alpha}'] = np.nan
            df[f'significant_fit_{alpha}'] = False
        
        return df

    # Add alternative FDR columns to Moran's I data
    print("\nAdding alternative FDR columns to Moran's I data:")
    for region in regions:
        if region in all_regions_normalized_moransi:
            print(f"  {region}:")
            all_regions_normalized_moransi[region] = add_alternative_fdr_column(
                all_regions_normalized_moransi[region], 
                alpha=alternative_fdr_threshold
            )

    # Add alternative FDR columns to Counts data
    print("\nAdding alternative FDR columns to Counts data:")
    for region in regions:
        if region in all_regions_normalized_counts:
            print(f"  {region}:")
            all_regions_normalized_counts[region] = add_alternative_fdr_column(
                all_regions_normalized_counts[region], 
                alpha=alternative_fdr_threshold
            )

    print(f"\n✓ Alternative FDR columns added to all regions!")
    print(f"New columns added:")
    print(f"  - fdr_corrected_fit_{alternative_fdr_threshold}")
    print(f"  - significant_fit_{alternative_fdr_threshold}")

    #%% Save Updated Results with New Columns
    # Save the updated results with the new FDR columns
    print(f"\n=== SAVING UPDATED RESULTS WITH ALTERNATIVE FDR COLUMNS ===")
    
    # Create output directory if it doesn't exist
    output_dir = os.path.join(base_path, "statistical_results")
    os.makedirs(output_dir, exist_ok=True)

    # Save Moran's I results with new columns
    print("\nSaving Moran's I results with alternative FDR columns...")
    for region in all_regions_normalized_moransi.keys():
        df = all_regions_normalized_moransi[region]
        
        # Create filename
        filename = f"{region}_morans_statistical_results_with_alt_fdr.csv"
        filepath = os.path.join(output_dir, filename)
        
        # Save to CSV
        df.to_csv(filepath)
        print(f"  ✓ {region}: {filename} ({df.shape[0]} genes, {df.shape[1]} columns)")

    # Save Counts results with new columns
    print("\nSaving Counts results with alternative FDR columns...")
    for region in all_regions_normalized_counts.keys():
        df = all_regions_normalized_counts[region]
        
        # Create filename
        filename = f"{region}_counts_statistical_results_with_alt_fdr.csv"
        filepath = os.path.join(output_dir, filename)
        
        # Save to CSV
        df.to_csv(filepath)
        print(f"  ✓ {region}: {filename} ({df.shape[0]} genes, {df.shape[1]} columns)")

    print(f"\n✓ Updated results saved with alternative FDR columns!")
    print(f"Now you can compare:")
    print(f"  - fdr_corrected_fit (original α=0.05)")
    print(f"  - fdr_corrected_fit_{alternative_fdr_threshold} (alternative α={alternative_fdr_threshold})")
    print(f"  - significant_fit (original α=0.05)")
    print(f"  - significant_fit_{alternative_fdr_threshold} (alternative α={alternative_fdr_threshold})")

# print("\n" + "="*60)
# if SKIP_STATISTICAL_ANALYSIS:
#     print("SCRIPT COMPLETE - Heatmaps with clustering generated")
#     print("Statistical analysis was skipped (set SKIP_STATISTICAL_ANALYSIS = False to run)")
# else:
#     print("SCRIPT COMPLETE - Full analysis pipeline completed")
#     print("Generated: normalization, heatmaps, statistical tests, and FDR correction")
# print("="*60)

# %%
