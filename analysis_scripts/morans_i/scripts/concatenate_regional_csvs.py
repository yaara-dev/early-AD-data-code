#!/usr/bin/env python3
"""
Script to concatenate CSV files from different samples within each region.
For each region, combines Moran's I values and counts from different samples side by side,
adding sample name prefixes to column names.

This script is part of the Moran's I spatial autocorrelation analysis pipeline.
It should be run after the FOV grid-based analysis to prepare data for statistical analysis.
"""

import os
import pandas as pd
import glob
import argparse
from pathlib import Path

def get_sample_name_from_filename(filename):
    """Extract sample name from filename (everything before the region name)."""
    # Remove file extension and get basename
    basename = os.path.basename(filename).replace('.csv', '')
    
    # Find the region name in the filename
    regions = ['CA1', 'CA3', 'DG', 'SM', 'SLM', 'upper_CA1', 'inner_DG', 'under_DG']
    
    for region in regions:
        if region in basename:
            # Extract everything before the region name
            sample_name = basename.split(f'_{region}')[0]
            return sample_name
    
    # Fallback: return the filename without extension
    return basename

def concatenate_region_csvs(base_path, output_dir):
    """
    Concatenate CSV files from different samples within each region.
    
    Args:
        base_path: Path to the directory containing region folders
        output_dir: Directory to save the concatenated results
    """
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Get all region directories
    region_dirs = [d for d in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, d))]
    
    print(f"Found {len(region_dirs)} regions: {region_dirs}")
    
    for region in region_dirs:
        region_path = os.path.join(base_path, region)
        print(f"\nProcessing region: {region}")
        
        # Find all Moran's I files for this region (FOV-combined outputs)
        moran_files = glob.glob(os.path.join(region_path, f"*_{region}_morans_combined.csv"))
        
        # Find all counts files for this region (FOV-combined outputs)
        count_files = glob.glob(os.path.join(region_path, f"*_{region}_counts_combined.csv"))
        
        # Also check for direct FOV output files (format: SAMPLE_REGION_morans.csv)
        moran_files_direct = glob.glob(os.path.join(region_path, f"*_{region}_morans.csv"))
        count_files_direct = glob.glob(os.path.join(region_path, f"*_{region}_counts.csv"))
        
        # Combine both file types
        moran_files.extend(moran_files_direct)
        count_files.extend(count_files_direct)
        
        print(f"  Found {len(moran_files)} Moran's I files")
        print(f"  Found {len(count_files)} counts files")
        
        if not moran_files and not count_files:
            print(f"  No CSV files found for region {region}")
            continue
        
        # Process Moran's I files
        if moran_files:
            print(f"  Processing Moran's I files...")
            moran_dfs = []
            
            for file_path in sorted(moran_files):
                sample_name = get_sample_name_from_filename(file_path)
                print(f"    Reading {sample_name}...")
                
                try:
                    df = pd.read_csv(file_path)
                    
                    # Handle different file formats
                    # If file has GeneNames column, use it; otherwise use first column
                    if 'GeneNames' in df.columns:
                        gene_col = 'GeneNames'
                    else:
                        gene_col = df.columns[0]
                        df = df.rename(columns={gene_col: 'GeneNames'})
                    
                    # Add sample prefix to all columns except GeneNames
                    new_columns = ['GeneNames']
                    for col in df.columns:
                        if col != 'GeneNames':
                            new_columns.append(f"{sample_name}_{col}")
                    
                    df.columns = new_columns
                    moran_dfs.append(df)
                    
                except Exception as e:
                    print(f"    Error reading {file_path}: {e}")
                    continue
            
            if moran_dfs:
                # Merge all Moran's I dataframes
                print(f"    Merging {len(moran_dfs)} Moran's I dataframes...")
                merged_morans = moran_dfs[0]
                
                for df in moran_dfs[1:]:
                    merged_morans = merged_morans.merge(df, on='GeneNames', how='outer')
                
                # Save merged Moran's I results
                output_file = os.path.join(output_dir, f"{region}_all_samples_morans_concatenated.csv")
                merged_morans.to_csv(output_file, index=False)
                print(f"    Saved merged Moran's I to: {output_file}")
                print(f"    Shape: {merged_morans.shape}")
        
        # Process counts files
        if count_files:
            print(f"  Processing counts files...")
            count_dfs = []
            
            for file_path in sorted(count_files):
                sample_name = get_sample_name_from_filename(file_path)
                print(f"    Reading {sample_name}...")
                
                try:
                    df = pd.read_csv(file_path)
                    
                    # Handle different file formats
                    # If file has GeneNames column, use it; otherwise use first column
                    if 'GeneNames' in df.columns:
                        gene_col = 'GeneNames'
                    else:
                        gene_col = df.columns[0]
                        df = df.rename(columns={gene_col: 'GeneNames'})
                    
                    # Add sample prefix to all columns except GeneNames
                    new_columns = ['GeneNames']
                    for col in df.columns:
                        if col != 'GeneNames':
                            new_columns.append(f"{sample_name}_{col}")
                    
                    df.columns = new_columns
                    count_dfs.append(df)
                    
                except Exception as e:
                    print(f"    Error reading {file_path}: {e}")
                    continue
            
            if count_dfs:
                # Merge all counts dataframes
                print(f"    Merging {len(count_dfs)} counts dataframes...")
                merged_counts = count_dfs[0]
                
                for df in count_dfs[1:]:
                    merged_counts = merged_counts.merge(df, on='GeneNames', how='outer')
                
                # Save merged counts results
                output_file = os.path.join(output_dir, f"{region}_all_samples_counts_concatenated.csv")
                merged_counts.to_csv(output_file, index=False)
                print(f"    Saved merged counts to: {output_file}")
                print(f"    Shape: {merged_counts.shape}")

def main():
    """Main function to run the concatenation script."""
    
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(
        description="Concatenate CSV files from different samples within each region",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scripts/concatenate_regional_csvs.py example_output/fov_grid_binary6_results
  python scripts/concatenate_regional_csvs.py example_output/fov_grid_binary6_results --output example_output/concatenated_results
        """
    )
    
    parser.add_argument(
        'base_path',
        help='Path to the directory containing region folders'
    )
    
    parser.add_argument(
        '--output', '-o',
        default=None,
        help='Output directory for concatenated results (default: base_path/concatenated_results)'
    )
    
    args = parser.parse_args()
    
    # Set output directory
    if args.output:
        output_dir = args.output
    else:
        output_dir = os.path.join(args.base_path, "concatenated_results")
    
    print("Starting CSV concatenation process...")
    print(f"Base path: {args.base_path}")
    print(f"Output directory: {output_dir}")
    
    if not os.path.exists(args.base_path):
        print(f"Error: Base path {args.base_path} does not exist!")
        return
    
    try:
        concatenate_region_csvs(args.base_path, output_dir)
        print("\n✅ CSV concatenation completed successfully!")
        print(f"Results saved to: {output_dir}")
        
    except Exception as e:
        print(f"\n❌ Error during concatenation: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()

