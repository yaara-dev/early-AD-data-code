#!/usr/bin/env python3
"""
Generate manifest.json file listing all CSV files in data/csvs/ directory.
This is run at build time to create a manifest for the web interface.
"""

import json
import os
from pathlib import Path

def generate_manifest():
    """Generate manifest.json from CSV files in data/csvs/"""
    csv_dir = Path("data/csvs")
    manifest_path = csv_dir / "manifest.json"
    
    if not csv_dir.exists():
        print(f"Warning: {csv_dir} does not exist")
        return
    
    # Find all CSV files
    csv_files = sorted(csv_dir.glob("*.csv"))
    
    if not csv_files:
        print(f"Warning: No CSV files found in {csv_dir}")
        # Create empty manifest
        manifest = []
    else:
        # Generate manifest entries
        manifest = []
        for csv_file in csv_files:
            # Get relative path from index.html (which is at repo root)
            relative_path = str(csv_file)
            
            # Generate display name from filename
            # Remove common suffixes and format nicely
            name = csv_file.stem
            # Remove common suffixes like "_regions_genes_with_cell_types"
            if "_regions_genes_with_cell_types" in name:
                name = name.replace("_regions_genes_with_cell_types", "")
            elif "_regions_genes" in name:
                name = name.replace("_regions_genes", "")
            
            manifest.append({
                "path": relative_path,
                "name": name
            })
    
    # Write manifest
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)
    
    print(f"Generated manifest with {len(manifest)} CSV files")
    print(f"Manifest saved to: {manifest_path}")
    
    return manifest

if __name__ == "__main__":
    generate_manifest()

