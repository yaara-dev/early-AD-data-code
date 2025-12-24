#!/usr/bin/env python3
"""
Generate manifest.json file listing all CSV files in data/csvs/ directory.
This is run at build time to create a manifest for the web interface.
"""

import json
import os
from pathlib import Path

def generate_manifest():
    """Generate manifest.json from CSV files in data/all_samples/"""
    csv_dir = Path("data/all_samples")
    manifest_path = Path("data/csvs/manifest.json")
    
    if not csv_dir.exists():
        print(f"Warning: {csv_dir} does not exist")
        return
    
    # Find normalized CSV files (_points.csv files)
    # Prefer normalized files over original files
    points_files = sorted(csv_dir.glob("*_points.csv"))
    
    # Also check for original files (for backward compatibility)
    original_files = sorted(csv_dir.glob("*.csv"))
    # Filter out backup files and normalized files
    original_files = [f for f in original_files 
                     if not f.name.endswith('.backup') 
                     and not f.name.endswith('_points.csv')
                     and not f.name.endswith('_regions.csv')
                     and not f.name.endswith('_points_regions.csv')
                     and not f.name.endswith('_cells.csv')
                     and not f.name.endswith('_points_cells.csv')]
    
    if not points_files and not original_files:
        print(f"Warning: No CSV files found in {csv_dir}")
        # Create empty manifest
        manifest = []
    else:
        # Generate manifest entries
        manifest = []
        
        # Add normalized files first (preferred)
        for csv_file in points_files:
            # Get relative path from index.html (which is at repo root)
            relative_path = str(csv_file)
            
            # Generate display name from filename
            # Remove _points suffix
            name = csv_file.stem.replace("_points", "")
            # Remove other common suffixes (for backward compatibility with old names)
            if "_with_regions_celltypes" in name:
                name = name.replace("_with_regions_celltypes", "")
            elif "_regions_genes_with_cell_types" in name:
                name = name.replace("_regions_genes_with_cell_types", "")
            elif "_regions_genes" in name:
                name = name.replace("_regions_genes", "")
            
            manifest.append({
                "path": relative_path,
                "name": name
            })
        
        # Add original files for backward compatibility
        for csv_file in original_files:
            relative_path = str(csv_file)
            name = csv_file.stem
            # Remove common suffixes
            if "_regions_genes_with_cell_types" in name:
                name = name.replace("_regions_genes_with_cell_types", "")
            elif "_regions_genes" in name:
                name = name.replace("_regions_genes", "")
            
            manifest.append({
                "path": relative_path,
                "name": name + " (original)"
            })
    
    # Write manifest
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)
    
    print(f"Generated manifest with {len(manifest)} CSV files")
    print(f"Manifest saved to: {manifest_path}")
    
    return manifest

if __name__ == "__main__":
    generate_manifest()

