#!/usr/bin/env python3
"""
Add cell type column to regions_genes CSV files by merging with cell type assignment CSVs.
"""

import pandas as pd
from pathlib import Path

# Sample name mapping: regions_genes sample name -> cell_type file name
SAMPLE_MAPPING = {
    'fem3_5x_E7_A_left': '5x_E7_A_left',
    'fem2_5x_F5_B_left': '5x_F5_B_left',
    'fem2_5x_F5_B_right': '5x_F5_B_right_cut',
    'fem2_WT_F3_B_left': 'WT_F3_B_left_updated',
    'fem3_WTE1_B_L': 'WTE1_B_L_UPDATED',
    'fem3_WTE1_B_R': 'WTE1_B_R',
    'fem4_5x_F8_A_R': '5x_F8_A_R',
    'fem4_WT_F11': 'WT_F11',
}

# Paths
CELL_TYPING_DIR = Path('/Users/yaarakarasik/data_for_publication/cell_typing_all_samples')
REGIONS_GENES_DIR = Path('/Users/yaarakarasik/assign_regions/data/outputs/output_correct')
OUTPUT_DIR = Path('/Users/yaarakarasik/data_for_publication/explore-ExSeq-brain-AD')


def add_cell_types(sample_name, cell_type_name):
    """Add cell_type column to a regions_genes CSV file."""
    print(f"\nProcessing {sample_name}...")
    
    # Construct file paths
    regions_genes_file = REGIONS_GENES_DIR / sample_name / f"{sample_name}_regions_genes.csv"
    cell_type_file = CELL_TYPING_DIR / f"cell_type_{cell_type_name}.csv"
    output_file = OUTPUT_DIR / f"{sample_name}_regions_genes_with_cell_types.csv"
    
    # Check if files exist
    if not regions_genes_file.exists():
        print(f"  ERROR: Regions genes file not found: {regions_genes_file}")
        return False
    
    if not cell_type_file.exists():
        print(f"  ERROR: Cell type file not found: {cell_type_file}")
        return False
    
    # Load data
    print(f"  Loading regions_genes from: {regions_genes_file}")
    regions_df = pd.read_csv(regions_genes_file)
    print(f"    Loaded {len(regions_df):,} rows")
    
    print(f"  Loading cell types from: {cell_type_file}")
    cell_type_df = pd.read_csv(cell_type_file)
    print(f"    Loaded {len(cell_type_df):,} cell type assignments")
    
    # Create a mapping dictionary from cell_index to cell_type
    # Handle potential string formatting differences (quotes, whitespace)
    cell_type_map = {}
    for _, row in cell_type_df.iterrows():
        cell_index = str(row['cell_index']).strip().strip('"')
        cell_type = str(row['cell_type']).strip().strip('"')
        cell_type_map[cell_index] = cell_type
    
    print(f"    Created mapping for {len(cell_type_map):,} cells")
    
    # Map cell_type to regions_df
    # Clean cell column values (remove quotes if present)
    regions_df['cell_clean'] = regions_df['cell'].astype(str).str.strip().str.strip('"')
    
    # Add cell_type column using map, defaulting to 'unassigned' for unmatched cells
    regions_df['cell_type'] = regions_df['cell_clean'].map(cell_type_map).fillna('unassigned')
    
    # Drop the temporary clean column
    regions_df = regions_df.drop(columns=['cell_clean'])
    
    # Count statistics
    assigned_count = (regions_df['cell_type'] != 'unassigned').sum()
    unassigned_count = (regions_df['cell_type'] == 'unassigned').sum()
    print(f"  Statistics:")
    print(f"    Assigned cells: {assigned_count:,} rows")
    print(f"    Unassigned cells: {unassigned_count:,} rows")
    
    # Save output
    print(f"  Saving to: {output_file}")
    regions_df.to_csv(output_file, index=False)
    print(f"  ✓ Successfully saved {len(regions_df):,} rows")
    
    return True


def main():
    """Process all sample pairs."""
    print("=" * 60)
    print("Adding Cell Types to Regions Genes CSVs")
    print("=" * 60)
    
    success_count = 0
    total_count = len(SAMPLE_MAPPING)
    
    for sample_name, cell_type_name in SAMPLE_MAPPING.items():
        if add_cell_types(sample_name, cell_type_name):
            success_count += 1
    
    print("\n" + "=" * 60)
    print(f"Processing complete: {success_count}/{total_count} samples processed successfully")
    print("=" * 60)


if __name__ == '__main__':
    main()


