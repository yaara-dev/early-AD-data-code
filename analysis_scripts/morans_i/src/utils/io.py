#!/usr/bin/env python3
"""
I/O Utilities

This module contains utilities for data loading, saving, and progress tracking
for Moran's I analysis workflows.
"""

import pandas as pd
import numpy as np
import time
from typing import Dict, List, Optional, Any
from pathlib import Path


def load_spatial_data(filename: str, 
                     coordinate_conversion: bool = True,
                     gene_mapping: bool = True,
                     filter_unassigned: bool = True) -> pd.DataFrame:
    """
    Load and preprocess spatial sequencing data.
    
    Parameters:
        filename: Path to CSV file with spatial data
        coordinate_conversion: Whether to convert coordinates to microns
        gene_mapping: Whether to apply gene name standardization
        filter_unassigned: Whether to filter unassigned cells
        
    Returns:
        Preprocessed DataFrame
        
    Example:
        >>> # data = load_spatial_data('sample.csv')
        >>> # len(data.columns) >= 5  # X, Y, Z, Gene, Cell
        True
    """
    from .spatial import (convert_coordinates_to_microns, apply_gene_name_mapping, 
                         filter_unassigned_cells)
    
    print(f"Loading spatial data from: {filename}")
    data = pd.read_csv(filename)
    print(f"Original data shape: {data.shape}")
    
    if coordinate_conversion:
        print("Converting coordinates to microns...")
        data = convert_coordinates_to_microns(data)
        
        # Create clean DataFrame with standardized columns
        processed_data = pd.DataFrame({
            'X': data['X'],
            'Y': data['Y'], 
            'Z': data['Z'],
            'Gene': data['gene'],
            'Cell': data['cell']
        })
        data = processed_data
    else:
        # Even without coordinate conversion, standardize column names
        print("Standardizing column names...")
        processed_data = pd.DataFrame({
            'X': data['x'],
            'Y': data['y'], 
            'Z': data['z'],
            'Gene': data['gene'],
            'Cell': data['cell']
        })
        data = processed_data
    
    if gene_mapping:
        print("Applying gene name mapping...")
        # Determine the correct gene column name
        gene_col = 'Gene' if 'Gene' in data.columns else 'gene'
        data = apply_gene_name_mapping(data, gene_col=gene_col)
    
    if filter_unassigned:
        original_count = len(data)
        data = filter_unassigned_cells(data)
        print(f"Filtered unassigned cells: {original_count:,} → {len(data):,} points")
    
    print(f"Final data shape: {data.shape}")
    return data


def save_morans_results(results: Dict[str, float], 
                       output_filename: str,
                       sample_id: str,
                       gene_counts: Dict[str, int] = None) -> None:
    """
    Save Moran's I results to CSV file.
    
    Parameters:
        results: Dictionary mapping gene names to Moran's I values
        output_filename: Output CSV filename
        sample_id: Sample identifier for column naming
        gene_counts: Optional dictionary mapping gene names to counts for zero count filtering
        
    Example:
        >>> results = {'ACTB': 0.5, 'GAPDH': 0.3}
        >>> save_morans_results(results, 'test.csv', 'SAMPLE')
    """
    from .spatial import get_global_gene_list
    
    print(f"Saving results to: {output_filename}")
    
    # Get global gene list and create results array
    global_genes = get_global_gene_list()
    global_results = []
    
    # Map results to global gene list, replacing zero Moran's I with N/A for zero count genes
    found_genes = 0
    for gene in global_genes:
        if gene in results:
            value = results[gene]
            # If gene counts provided and gene has zero count, show N/A in Moran's I
            if gene_counts and gene in gene_counts and gene_counts[gene] == 0:
                global_results.append('N/A')
            else:
                global_results.append(value)
            found_genes += 1
        else:
            # Gene not found in results, set to N/A
            global_results.append('N/A')
    
    print(f"Mapped {found_genes}/{len(global_genes)} genes from global list")
    
    # Create DataFrame and save
    results_df = pd.DataFrame({
        'GeneNames': global_genes,
        sample_id: global_results
    })
    
    results_df.to_csv(output_filename, index=False)
    print(f"✅ Results saved to: {output_filename}")
    
    # Log zero count filtering if applicable
    if gene_counts:
        zero_count_genes = [gene for gene, count in gene_counts.items() if count == 0]
        if zero_count_genes:
            print(f"   ⚠️  Genes with zero counts show N/A in Moran's I: {', '.join(zero_count_genes)}")


def print_progress(current: int, total: int, start_time: float, 
                  label: str = "Processing", frequency: int = 250) -> None:
    """
    Print progress updates with timing information.
    
    Parameters:
        current: Current item number (0-indexed)
        total: Total number of items
        start_time: Start time from time.time()
        label: Description label
        frequency: Print frequency (every N items)
        
    Example:
        >>> import time
        >>> start = time.time()
        >>> print_progress(500, 1000, start, "Processing cells")
    """
    if current % frequency == 0 or current == total - 1:
        elapsed = time.time() - start_time
        percent = 100.0 * (current + 1) / total
        print(f"  📊 {label}: {current + 1:,}/{total:,} ({percent:.1f}%) - {elapsed:.1f}s elapsed", 
              flush=True)


def validate_gene_coverage(data: pd.DataFrame, 
                         gene_col: str = 'Gene') -> Dict[str, Any]:
    """
    Validate gene coverage against global gene list.
    
    Parameters:
        data: DataFrame with gene column
        gene_col: Name of gene column
        
    Returns:
        Dictionary with validation statistics
        
    Example:
        >>> data = pd.DataFrame({'Gene': ['ACTB', 'GAPDH', 'Unknown']})
        >>> stats = validate_gene_coverage(data)
        >>> 'found_genes' in stats
        True
    """
    from .spatial import get_global_gene_list
    
    global_genes = set(get_global_gene_list())
    data_genes = set(data[gene_col].unique())
    
    found_genes = global_genes.intersection(data_genes)
    missing_genes = global_genes - data_genes
    extra_genes = data_genes - global_genes
    
    stats = {
        'total_global_genes': len(global_genes),
        'found_genes': len(found_genes),
        'missing_genes': len(missing_genes),
        'extra_genes': len(extra_genes),
        'coverage_percentage': 100.0 * len(found_genes) / len(global_genes),
        'found_gene_list': sorted(list(found_genes)),
        'missing_gene_list': sorted(list(missing_genes)),
        'extra_gene_list': sorted(list(extra_genes))
    }
    
    return stats


def print_gene_validation_report(stats: Dict[str, Any]) -> None:
    """
    Print a formatted gene validation report.
    
    Parameters:
        stats: Dictionary from validate_gene_coverage()
        
    Example:
        >>> stats = {'found_genes': 95, 'total_global_genes': 101}
        >>> print_gene_validation_report(stats)
    """
    print(f"\n" + "="*50)
    print(f"GENE VALIDATION REPORT")
    print(f"="*50)
    print(f"Global genes: {stats['total_global_genes']}")
    print(f"Found genes: {stats['found_genes']}")
    print(f"Missing genes: {stats['missing_genes']}")
    print(f"Extra genes: {stats['extra_genes']}")
    print(f"Coverage: {stats['coverage_percentage']:.1f}%")
    
    if stats['missing_genes'] > 0:
        print(f"\nMissing genes: {', '.join(stats['missing_gene_list'][:10])}")
        if stats['missing_genes'] > 10:
            print(f"... and {stats['missing_genes'] - 10} more")
    
    if stats['extra_genes'] > 0:
        print(f"\nExtra genes: {', '.join(stats['extra_gene_list'][:10])}")
        if stats['extra_genes'] > 10:
            print(f"... and {stats['extra_genes'] - 10} more")


def ensure_output_directory(filepath: str) -> None:
    """
    Ensure the output directory exists for a given filepath.
    
    Parameters:
        filepath: Output file path
        
    Example:
        >>> ensure_output_directory('results/analysis/output.csv')
    """
    output_dir = Path(filepath).parent
    output_dir.mkdir(parents=True, exist_ok=True)

