#!/usr/bin/env python3
"""
FOV-Specific Grid-Based Moran's I Analysis

For each region in a region assignment CSV, group by FOV, filter rows for each region-FOV combination,
standardize columns to x,y,z,gene,cell, run the grid-based pipeline, and combine per-FOV results
into per-region wide CSV files with FOV size information and flags.
"""

import os
import sys
import argparse
import tempfile
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))


def load_region_assignments(region_file: str) -> pd.DataFrame:
    """Load and validate region assignment data with FOV information."""
    if not os.path.exists(region_file):
        raise FileNotFoundError(f"Region file not found: {region_file}")
    
    df = pd.read_csv(region_file)
    required = ['region_name', 'gene', 'global_x', 'global_y', 'Z', 'cell', 'fov']
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    
    print(f"📥 Loaded {len(df):,} gene-region-FOV assignments")
    print(f"🧠 Found {df['region_name'].nunique()} regions")
    print(f"🔍 Found {df['fov'].nunique()} FOVs")
    print(f"🧬 Total unique genes: {df['gene'].nunique()}")
    
    # Apply gene name mapping to ensure consistency with Moran's I analysis
    print("🔄 Applying gene name mapping for consistency...")
    from utils.spatial import apply_gene_name_mapping
    df = apply_gene_name_mapping(df, gene_col='gene')
    
    print(f"✅ Gene names standardized. Total unique genes after mapping: {df['gene'].nunique()}")
    
    return df


def calculate_fov_coverage(fov_data: pd.DataFrame) -> tuple[float, int, str]:
    """
    Calculate FOV coverage area and size flag.
    
    Returns:
        coverage_area_um2: Area in square microns
        point_count: Number of data points
        size_flag: Size classification (FULL, PARTIAL, SMALL)
    """
    if len(fov_data) == 0:
        return 0.0, 0, "EMPTY"
    
    # Calculate coverage area from coordinate ranges
    x_range = fov_data['global_x'].max() - fov_data['global_x'].min()
    y_range = fov_data['global_y'].max() - fov_data['global_y'].min()
    coverage_area_um2 = x_range * y_range
    point_count = len(fov_data)
    
    # Expected full FOV size: 350x350 = 122,500 μm²
    expected_full_area = 122500.0
    coverage_ratio = coverage_area_um2 / expected_full_area
    
    # Determine size flag
    if coverage_ratio >= 0.8:
        size_flag = "FULL"
    elif coverage_ratio >= 0.5:
        size_flag = "PARTIAL"
    else:
        size_flag = "SMALL"
    
    return coverage_area_um2, point_count, size_flag


def calculate_gene_counts(fov_df: pd.DataFrame) -> pd.Series:
    """Calculate gene counts for a single FOV."""
    return fov_df['gene'].value_counts()


def prepare_fov_grid_input(df: pd.DataFrame, region_name: str, fov_id: str) -> pd.DataFrame:
    """Prepare region-FOV specific data for grid-based analysis."""
    # Filter to this region and FOV
    sub = df[(df['region_name'] == region_name) & (df['fov'] == fov_id)].copy()
    
    if sub.empty:
        raise ValueError(f"No data for region {region_name}, FOV {fov_id}")
    
    # The grid_based pipeline expects lowercase column names when coordinate_conversion=False
    # Rename to match what load_spatial_data expects
    out = sub[['gene', 'global_x', 'global_y', 'Z', 'cell']].rename(columns={
        'gene': 'gene',
        'global_x': 'x',         # Rename to lowercase for pipeline compatibility
        'global_y': 'y',         # Rename to lowercase for pipeline compatibility
        'Z': 'z',                # Rename to lowercase for pipeline compatibility
        'cell': 'cell',
    })
    
    return out


def run_grid_for_fov(fov_df: pd.DataFrame,
                     region_name: str,
                     fov_id: str,
                     sample_id: str,
                     output_dir: str,
                     grid_size: float,
                     weight_scheme: str,
                     coverage_area: float,
                     point_count: int,
                     size_flag: str) -> tuple[dict, dict]:
    """Run grid-based Moran's I analysis for a single region-FOV combination."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Create temporary input file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as tmp:
        fov_df.to_csv(tmp.name, index=False)
        tmp_input = tmp.name
    
    try:
        # Create temporary output file (we won't keep this)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as tmp_out:
            tmp_output = tmp_out.name
        
        fov_sample_id = f"{sample_id}_{region_name}_{fov_id}"
        
        # Run grid-based analysis
        cmd = [
            sys.executable,
            "src/pipelines/grid_based_morans.py",
            "--input", tmp_input,
            "--output", tmp_output,
            "--sample-id", fov_sample_id,
            "--grid-size", str(grid_size),
            "--weight-scheme", weight_scheme,
        ]
        
        print(f"   🚀 Running grid-based for {region_name} FOV{fov_id}...")
        print(f"      Grid size: {grid_size} μm")
        print(f"      Weight scheme: {weight_scheme}")
        res = subprocess.run(cmd, capture_output=True, text=True, cwd=".")
        
        if res.returncode != 0:
            print(f"   ❌ Error running grid-based for {region_name} FOV{fov_id}")
            print(f"   STDOUT: {res.stdout}")
            print(f"   STDERR: {res.stderr}")
            return None, {}
        
        print(f"   ✅ {region_name} FOV{fov_id} grid-based completed")
        
        # Load results and calculate gene counts
        results_df = pd.read_csv(tmp_output)
        morans_col = [c for c in results_df.columns if c.startswith('MoransI_')][0]
        
        # Create results dictionary
        results = {
            'morans_i': results_df.set_index('GeneNames')[morans_col],
            'gene_counts': calculate_gene_counts(fov_df)
        }
        
        # Return results and metadata
        metadata = {
            'coverage_area_um2': coverage_area,
            'point_count': point_count,
            'size_flag': size_flag,
            'fov_id': fov_id
        }
        
        return results, metadata
        
    finally:
        try:
            os.unlink(tmp_input)
            os.unlink(tmp_output)
        except Exception:
            pass


def save_individual_fov_csvs(per_fov_results: list[tuple], region_name: str, sample_id: str, output_dir: str) -> tuple[str, str] | None:
    """
    Save individual FOV CSV files and create sample-region combined files.
    
    Args:
        per_fov_results: List of (results, metadata) tuples from FOV analysis
        region_name: Name of the region being processed
        sample_id: Sample identifier
        output_dir: Output directory for results
        
    Returns:
        Tuple of (morans_path, counts_path) for the sample-region combined files
    """
    if not per_fov_results:
        print(f"❌ No FOV results to save for region {region_name}")
        return None
    
    print(f"   🔄 Saving individual FOV CSV files for {len(per_fov_results)} FOVs...")
    
    # Create sample-specific FOV folder structure
    # The sample_id might include the region name, so we need to extract just the sample part
    # For example: "FEM2_5X_F5_B_LEFT_CA1" -> sample: "FEM2_5X_F5_B_LEFT", region: "CA1"
    if sample_id.endswith(f"_{region_name}"):
        # Extract the sample name without the region suffix
        actual_sample_id = sample_id.replace(f"_{region_name}", "")
    else:
        actual_sample_id = sample_id
    
    # Create folders under the FOV output directory
    # Sample-specific FOV folders: output_dir/sample_name/region_name_fovs/
    sample_fov_dir = os.path.join(output_dir, actual_sample_id, f"{region_name}_fovs")
    Path(sample_fov_dir).mkdir(parents=True, exist_ok=True)
    
    # Region folders for combined files: output_dir/region_name/
    region_dir = os.path.join(output_dir, region_name)
    Path(region_dir).mkdir(parents=True, exist_ok=True)
    
    morans_map = {}
    counts_map = {}
    gene_index = None
    
    for results, metadata in per_fov_results:
        if not results:
            continue
            
        fov_id = metadata.get('fov_id', 'unknown_fov')
        
        # Store for combined file
        morans_map[fov_id] = results['morans_i']
        counts_map[fov_id] = results['gene_counts']
        
        # Set gene index
        if gene_index is None:
            gene_index = morans_map[fov_id].index
        
        # Create individual FOV CSV with metadata
        fov_morans_df = pd.DataFrame({
            'GeneNames': results['morans_i'].index,
            'morans_i': results['morans_i'].values,
            'coverage_um2': metadata['coverage_area_um2'],
            'point_count': metadata['point_count'],
            'size_flag': metadata['size_flag']
        })
        
        fov_counts_df = pd.DataFrame({
            'GeneNames': results['gene_counts'].index,
            'gene_count': results['gene_counts'].values,
            'coverage_um2': metadata['coverage_area_um2'],
            'point_count': metadata['point_count'],
            'size_flag': metadata['size_flag']
        })
        
        # Save individual FOV files
        fov_morans_path = os.path.join(sample_fov_dir, f"{fov_id}_morans.csv")
        fov_counts_path = os.path.join(sample_fov_dir, f"{fov_id}_counts.csv")
        
        fov_morans_df.to_csv(fov_morans_path, index=False)
        fov_counts_df.to_csv(fov_counts_path, index=False)
        
        print(f"      💾 Saved FOV {fov_id}: {os.path.basename(fov_morans_path)}, {os.path.basename(fov_counts_path)}")
    
    if not morans_map:
        print(f"❌ No valid FOV results to save for region {region_name}")
        return None
    
    # Create sample-region combined files (data only, no metadata)
    morans_df = pd.DataFrame(morans_map, index=gene_index).reset_index()
    morans_df = morans_df.rename(columns={'index': 'GeneNames'})
    
    counts_df = pd.DataFrame(counts_map, index=gene_index).reset_index()
    counts_df = counts_df.rename(columns={'index': 'GeneNames'})
    
    # Save sample-region combined files
    morans_path = os.path.join(region_dir, f"{actual_sample_id}_{region_name}_morans.csv")
    counts_path = os.path.join(region_dir, f"{actual_sample_id}_{region_name}_counts.csv")
    
    morans_df.to_csv(morans_path, index=False)
    counts_df.to_csv(counts_path, index=False)
    
    print(f"✅ Saved individual FOV files and sample-region combined files for {region_name}:")
    print(f"   Individual FOVs: {sample_fov_dir}")
    print(f"   Sample-region combined: {os.path.basename(morans_path)}")
    print(f"   FOVs included: {', '.join(sorted(morans_map.keys()))}")
    
    return morans_path, counts_path


def main():
    """Main entry point for FOV-specific grid-based Moran's I analysis."""
    ap = argparse.ArgumentParser(description="FOV-Specific Grid-Based Moran's I Analysis")
    ap.add_argument('--region-file', required=True, help='Region assignments CSV with FOV column')
    ap.add_argument('--output-dir', required=True, help='Output directory for results')
    ap.add_argument('--sample-id', required=True, help='Sample ID')
    ap.add_argument('--grid-size', type=float, default=35.0, help='Grid cell size in microns')
    ap.add_argument('--min-points-per-fov', type=int, default=50, help='Minimum points per FOV')
    ap.add_argument('--weight-scheme', choices=['inverse_square', 'inverse', 'binary_4', 'binary_6'], 
                    default='binary_6', help='Weight scheme: binary_6 (default), binary_4, inverse_square, or inverse')
    
    args = ap.parse_args()
    
    # Create output directory structure based on weight scheme
    fov_output_dir = os.path.join(args.output_dir, 'fovs', 'grid_based', args.weight_scheme)
    Path(fov_output_dir).mkdir(parents=True, exist_ok=True)
    
    print(f"🔬 FOV-SPECIFIC GRID-BASED MORAN'S I ANALYSIS")
    print(f"="*60)
    print(f"Region file: {args.region_file}")
    print(f"Output directory: {fov_output_dir}")
    print(f"Sample ID: {args.sample_id}")
    print(f"Grid size: {args.grid_size} μm")
    print(f"Weight scheme: {args.weight_scheme}")
    print(f"Min points per FOV: {args.min_points_per_fov}")
    print(f"="*60)
    
    # Load data
    df = load_region_assignments(args.region_file)
    
    # Get unique regions and FOVs
    regions = sorted(df['region_name'].unique())
    
    print(f"\n🎯 Processing {len(regions)} regions...")
    
    for region_name in regions:
        print(f"\n🧠 REGION: {region_name}")
        print("-" * 40)
        
        # Get FOVs for this region
        region_data = df[df['region_name'] == region_name]
        fovs = sorted(region_data['fov'].unique())
        
        print(f"   Found {len(fovs)} FOVs: {', '.join(map(str, fovs))}")
        
        # Process each FOV
        per_fov_results = []
        
        for fov_id in fovs:
            fov_data = region_data[region_data['fov'] == fov_id]
            
            # Check minimum points requirement
            if len(fov_data) < args.min_points_per_fov:
                print(f"   ⚠️  FOV{fov_id}: {len(fov_data)} points (below threshold, skipping)")
                continue
            
            # Calculate coverage info
            coverage_area, point_count, size_flag = calculate_fov_coverage(fov_data)
            print(f"   📊 FOV{fov_id}: {coverage_area:.0f} μm² coverage, {point_count} points ({size_flag})")
            
            # Prepare and run analysis
            try:
                fov_input = prepare_fov_grid_input(df, region_name, fov_id)
                results, metadata = run_grid_for_fov(
                    fov_input, region_name, fov_id, args.sample_id, fov_output_dir, 
                    args.grid_size, args.weight_scheme, coverage_area, point_count, size_flag
                )
                
                if results:
                    per_fov_results.append((results, metadata))
                    
            except Exception as e:
                print(f"   ❌ Error processing FOV{fov_id}: {e}")
                continue
        
        # Save individual FOV CSVs and create sample-region combined files
        if per_fov_results:
            save_individual_fov_csvs(per_fov_results, region_name, args.sample_id, fov_output_dir)
        else:
            print(f"   ❌ No valid FOVs for region {region_name}")
    
    print(f"\n🎉 FOV-specific grid-based analysis completed!")
    print(f"Results saved to: {fov_output_dir}")


if __name__ == '__main__':
    main()

