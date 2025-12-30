#!/usr/bin/env python3

"""
Grid-Based Moran's I Pipeline

This pipeline calculates Moran's I spatial autocorrelation for gene expression
data using a regular spatial grid as the spatial units.
"""

import sys
import time
import argparse
import numpy as np
import pandas as pd
import os
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from algorithms.morans_i import (
    process_morans_i_analysis
)
from utils.io import (
    load_spatial_data
)
from utils.spatial import (
    create_3d_grid,
    assign_to_3d_grid,
    create_3d_grid_centroids,
    get_global_gene_list
)


def run_grid_based_morans_analysis(input_file: str,
                                 output_file: str,
                                 sample_id: str,
                                 grid_size: float = 35.0,
                                 weight_scheme: str = "inverse_square") -> None:
    """
    Run complete grid-based Moran's I analysis.
    
    Parameters:
        input_file: Path to input CSV file
        output_file: Path to output CSV file
        sample_id: Sample identifier for output columns
        grid_size: Grid cell size in microns
        weight_scheme: 'inverse_square' (default), 'inverse', 'binary_4', or 'binary_6'
    """
    start_time = time.time()
    
    print(f"🏗️ 3D GRID-BASED MORAN'S I ANALYSIS")
    print(f"="*60)
    print(f"Input: {input_file}")
    print(f"Output: {output_file}")
    print(f"Sample ID: {sample_id}")
    print(f"Grid size: {grid_size}μm × {grid_size}μm (3D centroids)")
    print(f"Weight scheme: {weight_scheme}")
    print(f"="*60)
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file) if os.path.dirname(output_file) else '.', exist_ok=True)
    
    # Step 1: Load and preprocess data
    print(f"\n📥 STEP 1: Loading and preprocessing data")
    data = load_spatial_data(
        input_file,
        coordinate_conversion=False,
        gene_mapping=True, 
        filter_unassigned=False  # Grid-based: use all points
    )
    
    # Step 2: Data loaded successfully
    print(f"\n🔍 STEP 2: Data loaded successfully")
    print(f"Data shape: {data.shape}")
    print(f"Columns: {list(data.columns)}")
    
    # Step 3: Create 3D spatial grid
    print(f"\n🏗️ STEP 3: Creating 3D spatial grid")
    x_min, x_max = data['X'].min(), data['X'].max()
    y_min, y_max = data['Y'].min(), data['Y'].max()
    z_min, z_max = data['Z'].min(), data['Z'].max()
    
    print(f"Data extent: X=[{x_min:.1f}, {x_max:.1f}], Y=[{y_min:.1f}, {y_max:.1f}], Z=[{z_min:.1f}, {z_max:.1f}]")
    
    x_edges, y_edges, z_edges, x_centers, y_centers, z_centers = create_3d_grid(
        x_min, x_max, y_min, y_max, z_min, z_max, grid_size
    )
    
    print(f"3D Grid dimensions: {len(x_centers)} × {len(y_centers)} × {len(z_centers)} = {len(x_centers) * len(y_centers) * len(z_centers):,} total cells")
    
    # Step 4: Assign points to 3D grid
    print(f"\n📍 STEP 4: Assigning points to 3D grid cells")
    max_x_cells = len(x_centers)  # Number of cells in X direction
    max_y_cells = len(y_centers)  # Number of cells in Y direction
    
    grid_ids = assign_to_3d_grid(
        data['X'].values, data['Y'].values, data['Z'].values,
        x_min, y_min, z_min, grid_size, max_x_cells, max_y_cells
    )
    data['GridCell'] = grid_ids
    
    unique_grid_cells = np.unique(grid_ids)
    print(f"Points assigned to {len(unique_grid_cells):,} 3D grid cells")
    
    # Step 5: Create 3D grid centroids
    print(f"\n🎯 STEP 5: Creating 3D grid centroids")
    
    # Calculate centroids only for active grid cells (those with data)
    centroid_matrix = np.zeros((len(unique_grid_cells), 3))
    
    for idx, grid_id in enumerate(unique_grid_cells):
        # Get all points in this grid cell
        cell_data = data[data['GridCell'] == grid_id]
        
        # Calculate actual 3D centroid (mean of all points in this grid cell)
        mean_x = cell_data['X'].mean()
        mean_y = cell_data['Y'].mean() 
        mean_z = cell_data['Z'].mean()
        
        centroid_matrix[idx, :] = [mean_x, mean_y, mean_z]
    
    print(f"✅ Created 3D centroids for {len(unique_grid_cells):,} active grid cells")
    
    # Step 6: Build expression matrix
    print(f"\n🧮 STEP 6: Building expression matrix")
    global_genes = get_global_gene_list()
    
    # Optimized: use pandas pivot_table for much faster expression matrix creation
    expression_df = data.groupby(['GridCell', 'Gene']).size().reset_index(name='count')
    expression_matrix = expression_df.pivot_table(
        index='GridCell', 
        columns='Gene', 
        values='count', 
        fill_value=0
    ).reindex(columns=global_genes, fill_value=0).values
    
    print(f"Expression matrix shape: {expression_matrix.shape}")
    
    # Step 7: Calculate Moran's I
    print(f"\n🔢 STEP 7: Calculating Moran's I for {len(global_genes)} genes")
    morans_start = time.time()
    
    results = process_morans_i_analysis(
        centroid_matrix,
        expression_matrix, 
        global_genes,
        weight_scheme=weight_scheme,
        grid_ids=unique_grid_cells if weight_scheme in ["binary_4", "binary_6"] else None,
        max_x_cells=max_x_cells if weight_scheme in ["binary_4", "binary_6"] else None,
        max_y_cells=max_y_cells if weight_scheme in ["binary_4", "binary_6"] else None,
    )
    
    morans_time = time.time() - morans_start
    print(f"✅ Moran's I calculation completed in {morans_time:.1f}s")
    
    # Step 8: Save results
    print(f"\n💾 STEP 8: Saving results")
    
    # Calculate gene counts using pandas operations
    gene_counts = dict(zip(global_genes, expression_matrix.sum(axis=0).astype(int)))
    
    # Create DataFrame with both Moran's I and gene counts
    output_data = []
    for gene in global_genes:
        morans_i = results.get(gene, 'N/A')
        gene_count = gene_counts.get(gene, 0)
        output_data.append([gene, morans_i, gene_count])
    
    results_df = pd.DataFrame(output_data, columns=[
        'GeneNames', 
        f'MoransI_{sample_id}', 
        f'GeneCount_{sample_id}'
    ])
    results_df.to_csv(output_file, index=False)
    print(f"✅ Results saved to: {output_file}")
    
    # Summary statistics
    print(f"\n" + "="*60)
    print(f"ANALYSIS SUMMARY")
    print(f"="*60)
    
    # Show all genes with results (positive, negative, and zero Moran's I)
    genes_with_results = {k: v for k, v in results.items() if k in gene_counts and gene_counts[k] > 0}
    if genes_with_results:
        sorted_results = sorted(genes_with_results.items(), key=lambda x: abs(x[1]), reverse=True)
        mean_morans = np.mean(list(genes_with_results.values()))
        
        print(f"Grid size: {grid_size}μm × {grid_size}μm (3D analysis)")
        print(f"Total grid cells: {len(x_centers) * len(y_centers) * len(z_centers):,}")
        print(f"Active grid cells: {len(unique_grid_cells):,}")
        print(f"Total genes analyzed: {len(global_genes)}")
        print(f"Genes with expression data: {len(genes_with_results)}")
        print(f"Mean Moran's I: {mean_morans:.6f}")
        print(f"Total processing time: {time.time() - start_time:.1f}s")
        
        print(f"\nTop 5 genes by absolute Moran's I (showing clustering and dispersion):")
        for i, (gene, score) in enumerate(sorted_results[:5]):
            pattern = "clustering" if score > 0 else "dispersion" if score < 0 else "random"
            print(f"  {i+1}. {gene:<15} {score:.6f} ({pattern})")
    else:
        print("⚠️ No genes had expression data")
    
    print(f"\n🎉 Analysis completed successfully!")
    print(f"Results saved to: {output_file}")


def main():
    """Main entry point for grid-based Moran's I analysis."""
    parser = argparse.ArgumentParser(
        description="Grid-based Moran's I spatial autocorrelation analysis"
    )
    parser.add_argument(
        '--input', 
        required=True,
        help='Input CSV file with spatial data'
    )
    parser.add_argument(
        '--output',
        required=True, 
        help='Output CSV file for results'
    )
    parser.add_argument(
        '--sample-id',
        required=True,
        help='Sample identifier for output column naming'
    )
    parser.add_argument(
        '--grid-size',
        type=float,
        default=35.0,
        help='Grid cell size in microns (default: 35.0)'
    )
    parser.add_argument(
        '--weight-scheme',
        choices=['inverse_square', 'inverse', 'binary_4', 'binary_6'],
        default='inverse_square',
        help="Spatial weight scheme: 'inverse_square' (1/d^2, default), 'inverse' (1/d), 'binary_4' (XY adjacency), or 'binary_6' (3D adjacency)"
    )
    
    args = parser.parse_args()
    
    try:
        run_grid_based_morans_analysis(
            args.input,
            args.output,
            args.sample_id,
            args.grid_size,
            weight_scheme=args.weight_scheme,
        )
    except Exception as e:
        print(f"❌ Error: {e}")
        raise


if __name__ == "__main__":
    main()

