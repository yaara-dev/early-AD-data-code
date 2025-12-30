#!/usr/bin/env python3
"""
Create separate gene-cell matrices for each cell type within each sample.
Each matrix will have the same format as gene_cell_matrix.csv but only include cells of a specific type.
"""

import pandas as pd
import numpy as np
import os
import argparse
from pathlib import Path

def load_cell_type_mapping(sample_dir):
    """
    Load cell type mapping from RNA_with_cells.csv file.
    Returns a dictionary mapping cell_id to cell_type.
    """
    rna_file = os.path.join(sample_dir, 'RNA_with_cells.csv')
    if not os.path.exists(rna_file):
        print(f"Warning: RNA_with_cells.csv not found in {sample_dir}")
        return {}
    
    df = pd.read_csv(rna_file)
    # Create mapping from cell_id to cell_type
    cell_type_mapping = df.set_index('cell_id')['cell_type'].to_dict()
    return cell_type_mapping

def load_expression_matrix(sample_dir):
    """
    Load gene expression matrix from gene_cell_matrix.csv file.
    Returns a pandas DataFrame with genes as rows and cells as columns.
    """
    matrix_file = os.path.join(sample_dir, 'gene_cell_matrix.csv')
    if not os.path.exists(matrix_file):
        print(f"Warning: gene_cell_matrix.csv not found in {sample_dir}")
        return None
    
    # Load the expression matrix
    df = pd.read_csv(matrix_file, index_col=0)  # First column is gene names
    return df

def create_celltype_matrices_for_sample(sample_dir, sample_name, output_dir):
    """
    Create separate gene-cell matrices for each cell type in a sample.
    """
    print(f"Processing sample: {sample_name}")
    
    # Load cell type mapping
    cell_type_mapping = load_cell_type_mapping(sample_dir)
    if not cell_type_mapping:
        print(f"  No cell type mapping found for {sample_name}")
        return
    
    # Load expression matrix
    expression_matrix = load_expression_matrix(sample_dir)
    if expression_matrix is None:
        print(f"  No expression matrix found for {sample_name}")
        return
    
    print(f"  Found {len(expression_matrix)} genes and {len(expression_matrix.columns)} cells")
    
    # Get all unique cell types
    cell_types = set(cell_type_mapping.values())
    print(f"  Cell types: {sorted(cell_types)}")
    
    # Create output directory for this sample
    sample_output_dir = os.path.join(output_dir, sample_name)
    os.makedirs(sample_output_dir, exist_ok=True)
    
    # Create a matrix for each cell type
    for cell_type in cell_types:
        # Get cells of this type
        cells_of_type = [cell_id for cell_id, ct in cell_type_mapping.items() 
                        if ct == cell_type and cell_id in expression_matrix.columns]
        
        if cells_of_type:
            # Extract expression matrix for this cell type
            celltype_matrix = expression_matrix[cells_of_type].copy()
            
            # Clean cell type name for filename
            clean_cell_type = cell_type.replace(' ', '_').replace('/', '_')
            output_file = os.path.join(sample_output_dir, f"gene_cell_matrix_{clean_cell_type}.csv")
            
            # Save the matrix
            celltype_matrix.to_csv(output_file)
            
            print(f"    Created {clean_cell_type}: {len(celltype_matrix)} genes × {len(celltype_matrix.columns)} cells")
        else:
            print(f"    No cells found for cell type: {cell_type}")

def main():
    """
    Main function to process all samples and create cell-type-specific matrices.
    """
    parser = argparse.ArgumentParser(description='Create cell-type-specific gene expression matrices')
    parser.add_argument('--input-dir', '-i', 
                       required=True,
                       help='Directory containing sample folders (or single sample folder)')
    parser.add_argument('--output-dir', '-o',
                       required=True,
                       help='Output directory for celltype matrices')
    
    args = parser.parse_args()
    
    print("🧬 Creating Cell-Type-Specific Matrices")
    print("=" * 50)
    print(f"Input directory: {args.input_dir}")
    print(f"Output directory: {args.output_dir}")
    print()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Check if input is a single sample folder or directory of samples
    if os.path.isdir(args.input_dir):
        # Check if it's a sample folder (has gene_cell_matrix.csv) or a directory of samples
        if os.path.exists(os.path.join(args.input_dir, 'gene_cell_matrix.csv')):
            # Single sample folder
            sample_name = os.path.basename(args.input_dir)
            if sample_name == '':
                sample_name = os.path.basename(os.path.dirname(args.input_dir))
            create_celltype_matrices_for_sample(args.input_dir, sample_name, args.output_dir)
        else:
            # Directory of sample folders
            sample_dirs = [d for d in os.listdir(args.input_dir) 
                          if os.path.isdir(os.path.join(args.input_dir, d))]
            
            print(f"Found {len(sample_dirs)} sample directories")
            print()
            
            # Process each sample
            for sample_dir in sample_dirs:
                sample_path = os.path.join(args.input_dir, sample_dir)
                sample_name = sample_dir
                create_celltype_matrices_for_sample(sample_path, sample_name, args.output_dir)
                print()
    else:
        print(f"❌ Input directory not found: {args.input_dir}")
        return
    
    print("✅ All cell-type-specific matrices created!")
    print(f"Files saved to: {args.output_dir}")
    print("\nDirectory structure:")
    print("celltype_matrices/")
    print("├── [sample_name]/")
    print("│   ├── gene_cell_matrix_[cell_type_1].csv")
    print("│   ├── gene_cell_matrix_[cell_type_2].csv")
    print("│   └── ...")
    print("└── ...")

if __name__ == "__main__":
    main()
