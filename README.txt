ExSeq Brain AD - Spatial Genomics Visualization

================================================================================
DATA EXPLORER
================================================================================

The full interactive data explorer is available via GitHub Pages:

https://yaara-dev.github.io/explore-ExSeq-brain-AD/

Click the link above to explore the spatial genomics data with interactive 2D 
and 3D visualizations, filtering options, and comprehensive dashboards.

================================================================================
FEATURES
================================================================================

- Multi-sample support: Browse between different samples
- Interactive 2D scatter plot: Visualize spatial data with x/y coordinates
- Interactive 3D scatter plot: Explore data in three dimensions with rotation 
  and zoom
- Comprehensive dashboard: View statistics, heatmaps, distributions, and more
- Filtering: Filter by region, gene, cell type, and Z-slice
- Statistics: View total records, unique genes, regions, and visible points
- Tooltips: Hover over points to see detailed information

================================================================================
SAMPLES AVAILABLE
================================================================================

- WT_1
- WT_2.1
- WT_2.2
- WT_3
- 5xFAD_1.1
- 5xFAD_1.2
- 5xFAD_2
- 5xFAD_3

All samples use normalized CSV format with 5 files per sample located in 
data/all_samples/.

DATA FORMAT
-----------

Normalized CSV structure with 5 files per sample:
- *_points.csv - Point data (point_id, gene, Z, global_x, global_y, fov)
- *_regions.csv - Region information (region_id, region_name, region_area, 
  region_proportion)
- *_points_regions.csv - Point-to-region mapping (point_id, region_id)
- *_cells.csv - Cell information (cell_id, cell_type)
- *_points_cells.csv - Point-to-cell mapping (point_id, cell_id)

================================================================================
FOR DEVELOPERS
================================================================================

The following information is only relevant if you want to modify or extend the 
visualization code.

PROJECT STRUCTURE
-----------------

- index.html - Main visualization file (contains 2D view, 3D view, and 
  dashboard)
- data/csvs/manifest.json - Sample manifest file (auto-generated)
- data/all_samples/ - Normalized CSV files (5 files per sample)
- data_explorer_scripts/ - Data generation and normalization scripts
  - generate_manifest.py - Generates manifest.json from CSV files
  - add_cell_types.py - Adds cell type information to CSV files
  - server.py - Local development server with cache control

LOCAL DEVELOPMENT
-----------------

1. Generate the manifest file (if data changes):
   python3 data_explorer_scripts/generate_manifest.py

2. Run a local server:
   python3 data_explorer_scripts/server.py
   
   Or use Python's built-in server:
   python3 -m http.server 8000

3. Open in browser:
   - Go to http://localhost:8000
   - The visualization will automatically load the manifest and first sample

DEPLOYMENT
----------

The project uses GitHub Actions for automatic deployment. The workflow 
(.github/workflows/deploy.yml) automatically:
- Generates the manifest file
- Deploys to GitHub Pages on changes to index.html, data/all_samples/, or 
  data_explorer_scripts/generate_manifest.py

REQUIREMENTS
-----------

- A modern web browser (Chrome, Firefox, Safari, Edge)
- Python 3 (for running scripts and local development)
- Normalized CSV files in data/all_samples/ folder

