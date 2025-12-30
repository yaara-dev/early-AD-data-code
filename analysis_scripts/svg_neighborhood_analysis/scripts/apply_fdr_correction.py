#!/usr/bin/env python3
"""
FDR Correction Script for CELINA Results

This script applies FDR correction to p-values in CELINA results across multiple samples.
It provides options for different correction strategies to handle multiple comparisons.
"""

import os
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
import argparse
from pathlib import Path

def apply_fdr_correction(data, method='fdr_bh', alpha=0.05, strategy='per_celltype'):
    """
    Apply FDR correction to p-values with different strategies.
    
    Parameters:
    -----------
    data : pandas.DataFrame
        DataFrame with genes as index and cell types as columns
    method : str
        FDR correction method ('fdr_bh', 'fdr_by', 'bonferroni')
    alpha : float
        Significance level
    strategy : str
        Correction strategy:
        - 'per_celltype': Correct within each cell type separately
        - 'global': Correct across all p-values from all cell types
        - 'per_sample': Correct within each sample (all cell types combined)
    
    Returns:
    --------
    pandas.DataFrame
        DataFrame with FDR-corrected p-values
    """
    
    # Remove NaN values for correction
    data_clean = data.copy()
    
    if strategy == 'per_celltype':
        # Apply FDR correction within each cell type
        corrected_data = data_clean.copy()
        for col in data_clean.columns:
            p_values = data_clean[col].dropna()
            if len(p_values) > 0:
                # Apply FDR correction
                rejected, p_corrected, _, _ = multipletests(
                    p_values.values, 
                    alpha=alpha, 
                    method=method
                )
                
                # Create a series with corrected p-values
                corrected_series = pd.Series(index=p_values.index, data=p_corrected)
                corrected_data[col] = corrected_series
        
    elif strategy == 'global':
        # Apply FDR correction across all p-values
        all_p_values = data_clean.values.flatten()
        valid_indices = ~pd.isna(all_p_values)
        
        if np.sum(valid_indices) > 0:
            valid_p_values = all_p_values[valid_indices]
            rejected, p_corrected, _, _ = multipletests(
                valid_p_values, 
                alpha=alpha, 
                method=method
            )
            
            # Reconstruct the corrected matrix
            corrected_data = data_clean.copy()
            corrected_values = np.full_like(all_p_values, np.nan)
            corrected_values[valid_indices] = p_corrected
            corrected_data.iloc[:, :] = corrected_values.reshape(data_clean.shape)
        else:
            corrected_data = data_clean.copy()
            
    elif strategy == 'per_sample':
        # Apply FDR correction across all cell types in the sample
        all_p_values = data_clean.values.flatten()
        valid_indices = ~pd.isna(all_p_values)
        
        if np.sum(valid_indices) > 0:
            valid_p_values = all_p_values[valid_indices]
            rejected, p_corrected, _, _ = multipletests(
                valid_p_values, 
                alpha=alpha, 
                method=method
            )
            
            # Reconstruct the corrected matrix
            corrected_data = data_clean.copy()
            corrected_values = np.full_like(all_p_values, np.nan)
            corrected_values[valid_indices] = p_corrected
            corrected_data.iloc[:, :] = corrected_values.reshape(data_clean.shape)
        else:
            corrected_data = data_clean.copy()
    
    return corrected_data

def process_sample_folder(sample_folder, method='fdr_bh', alpha=0.05, strategy='per_celltype', output_suffix='_FDR_corrected'):
    """
    Process a single sample folder and apply FDR correction.
    
    Parameters:
    -----------
    sample_folder : str
        Path to the sample folder
    method : str
        FDR correction method
    alpha : float
        Significance level
    strategy : str
        Correction strategy
    output_suffix : str
        Suffix for output filename
    
    Returns:
    --------
    dict
        Summary statistics for the sample
    """
    
    # Find the all_p_values file
    p_values_file = os.path.join(sample_folder, 'all_p_values_CombinedPvals.csv')
    
    if not os.path.exists(p_values_file):
        print(f"⚠️  No all_p_values_CombinedPvals.csv found in {sample_folder}")
        return None
    
    print(f"📊 Processing {os.path.basename(sample_folder)}...")
    
    # Read the p-values data
    try:
        data = pd.read_csv(p_values_file, index_col=0)
        print(f"   Loaded {data.shape[0]} genes × {data.shape[1]} cell types")
    except Exception as e:
        print(f"❌ Error reading {p_values_file}: {e}")
        return None
    
    # Apply FDR correction
    try:
        corrected_data = apply_fdr_correction(data, method=method, alpha=alpha, strategy=strategy)
        
        # Save corrected data
        output_file = os.path.join(sample_folder, f'all_p_values{output_suffix}.csv')
        corrected_data.to_csv(output_file)
        print(f"   ✅ Saved corrected p-values to {os.path.basename(output_file)}")
        
        # Calculate summary statistics
        summary = {
            'sample': os.path.basename(sample_folder),
            'total_genes': data.shape[0],
            'total_cell_types': data.shape[1],
            'total_tests': data.count().sum(),
            'significant_before': (data < alpha).sum().sum(),
            'significant_after': (corrected_data < alpha).sum().sum(),
            'correction_strategy': strategy,
            'fdr_method': method,
            'alpha': alpha
        }
        
        # Per cell type statistics
        for col in data.columns:
            before_sig = (data[col] < alpha).sum()
            after_sig = (corrected_data[col] < alpha).sum()
            summary[f'{col}_before'] = before_sig
            summary[f'{col}_after'] = after_sig
        
        return summary
        
    except Exception as e:
        print(f"❌ Error applying FDR correction: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(description='Apply FDR correction to CELINA p-values')
    parser.add_argument('--input-dir', '-i', 
                       required=True,
                       help='Directory containing sample folders (or single sample folder)')
    parser.add_argument('--method', '-m', 
                       choices=['fdr_bh', 'fdr_by', 'bonferroni'],
                       default='fdr_bh',
                       help='FDR correction method (default: fdr_bh)')
    parser.add_argument('--alpha', '-a', 
                       type=float, default=0.05,
                       help='Significance level (default: 0.05)')
    parser.add_argument('--strategy', '-s',
                       choices=['per_celltype', 'global', 'per_sample'],
                       default='global',
                       help='Correction strategy: per_celltype, global, or per_sample (default: global)')
    parser.add_argument('--output-suffix', '-o',
                       default='_FDR_global_bonferroni',
                       help='Suffix for output files (default: _FDR_global_bonferroni)')
    parser.add_argument('--summary-file', 
                       default='fdr_correction_summary.csv',
                       help='Output file for summary statistics (default: fdr_correction_summary.csv)')
    
    args = parser.parse_args()
    
    print("🧬 FDR Correction for CELINA Results")
    print("=" * 50)
    print(f"Input directory: {args.input_dir}")
    print(f"Method: {args.method}")
    print(f"Alpha: {args.alpha}")
    print(f"Strategy: {args.strategy}")
    print(f"Output suffix: {args.output_suffix}")
    print()
    
    # Check if input is a single file or directory
    if os.path.isfile(args.input_dir):
        # Single file mode
        sample_folder = os.path.dirname(args.input_dir)
        if sample_folder == '':
            sample_folder = '.'
        summaries = []
        summary = process_sample_folder(
            sample_folder, 
            method=args.method, 
            alpha=args.alpha, 
            strategy=args.strategy,
            output_suffix=args.output_suffix
        )
        if summary:
            summaries.append(summary)
    else:
        # Directory mode - find all sample folders
        sample_folders = []
        for item in os.listdir(args.input_dir):
            item_path = os.path.join(args.input_dir, item)
            if os.path.isdir(item_path):
                sample_folders.append(item_path)
        
        if not sample_folders:
            print("❌ No sample folders found!")
            return
        
        print(f"Found {len(sample_folders)} sample folders")
        print()
        
        # Process each sample
        summaries = []
        for folder in sample_folders:
            summary = process_sample_folder(
                folder, 
                method=args.method, 
                alpha=args.alpha, 
                strategy=args.strategy,
                output_suffix=args.output_suffix
            )
            if summary:
                summaries.append(summary)
            print()
    
    # Save summary
    if summaries:
        summary_df = pd.DataFrame(summaries)
        summary_file = os.path.join(args.input_dir, args.summary_file)
        summary_df.to_csv(summary_file, index=False)
        print(f"📋 Summary saved to {summary_file}")
        
        # Print overall statistics
        print("\n📊 Overall Statistics:")
        print(f"   Total samples processed: {len(summaries)}")
        print(f"   Total significant before correction: {summary_df['significant_before'].sum()}")
        print(f"   Total significant after correction: {summary_df['significant_after'].sum()}")
        if summary_df['significant_before'].sum() > 0:
            print(f"   Average reduction: {((summary_df['significant_before'] - summary_df['significant_after']) / summary_df['significant_before'] * 100).mean():.1f}%")
    
    print("\n✅ FDR correction completed!")

if __name__ == "__main__":
    main()
