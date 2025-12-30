# ExSeq Brain AD - Spatial Genomics Data Explorer & Analysis Scripts

This repository contains both an interactive web-based data explorer and computational analysis pipelines for spatial transcriptomics data from the ExSeq Brain AD study.

---

## 🚀 Data Explorer

**The full interactive data explorer is available via GitHub Pages:**

👉 **[https://yaara-dev.github.io/explore-ExSeq-brain-AD/](https://yaara-dev.github.io/explore-ExSeq-brain-AD/)**

Click the link above to explore the spatial genomics data with interactive 2D and 3D visualizations, filtering options, and comprehensive dashboards.

**📋 For detailed information about data structure, column descriptions, and file relationships, see [data/README.md](data/README.md)**

### Features

- **Multi-sample support:** Browse between different samples
- **Interactive 2D scatter plot:** Visualize spatial data with x/y coordinates
- **Interactive 3D scatter plot:** Explore data in three dimensions with rotation and zoom
- **Comprehensive dashboard:** View statistics, heatmaps, distributions, and more
- **Filtering:** Filter by region, gene, cell type, and Z-slice
- **Statistics:** View total records, unique genes, regions, and visible points
- **Tooltips:** Hover over points to see detailed information

### Samples Available

- WT_1
- WT_2.1
- WT_2.2
- WT_3
- 5xFAD_1.1
- 5xFAD_1.2
- 5xFAD_2
- 5xFAD_3

All samples use normalized CSV format with 5 files per sample located in `data/all_samples/`.

---

## 📊 Analysis Scripts

This repository includes four computational analysis pipelines for spatial transcriptomics data:

1. **Cell Typing** - Automated cell type annotation using Seurat and marker gene catalogs
2. **Moran's I Spatial Autocorrelation** - Spatial pattern analysis using grid-based Moran's I statistics
3. **RNA Velocity** - Analysis of RNA velocity and spatial cell-state transitions
4. **SVG Neighborhood Analysis** - Post-CELINA analysis of spatially variable genes and co-regulated modules

Each analysis pipeline is self-contained with its own documentation, requirements, and example data.

**📖 For detailed information about each analysis, see [analysis_scripts/README.md](analysis_scripts/README.md)**

Each analysis directory contains its own README with:
- Installation instructions
- Input data format requirements
- Step-by-step usage examples
- Output format descriptions

---

## For Developers

The following information is only relevant if you want to modify or extend the visualization code.

### Project Structure

**Data Explorer Components:**
- `index.html` - Main visualization file (contains 2D view, 3D view, and dashboard)
- `data/csvs/manifest.json` - Sample manifest file (auto-generated)
- `data/all_samples/` - Normalized CSV files (5 files per sample)
- `data_explorer_scripts/` - Data generation and normalization scripts
  - `generate_manifest.py` - Generates manifest.json from CSV files
  - `add_cell_types.py` - Adds cell type information to CSV files
  - `server.py` - Local development server with cache control

**Analysis Scripts:**
- `analysis_scripts/` - Computational analysis pipelines
  - `cell_typing/` - Cell type annotation pipeline (R/Seurat)
  - `morans_i/` - Spatial autocorrelation analysis (Python)
  - `rna_velocity/` - RNA velocity analysis (Python)
  - `svg_neighborhood_analysis/` - SVG neighborhood analysis (Python)

### Data Format

The visualization uses a normalized CSV structure with 5 files per sample:
- `*_points.csv` - Point data (point_id, gene, Z, X, Y, fov)
- `*_regions.csv` - Region information (region_id, region_name, region_area, region_proportion)
- `*_points_regions.csv` - Point-to-region mapping (point_id, region_id)
- `*_cells.csv` - Cell information (cell_id, cell_type)
- `*_points_cells.csv` - Point-to-cell mapping (point_id, cell_id)

**📋 For detailed information about data structure, column descriptions, and file relationships, see [data/README.md](data/README.md)**

### Local Development

1. **Generate the manifest file** (if data changes):
   ```bash
   python3 data_explorer_scripts/generate_manifest.py
   ```

2. **Run a local server:**
   ```bash
   python3 data_explorer_scripts/server.py
   ```
   Or use Python's built-in server:
   ```bash
   python3 -m http.server 8000
   ```

3. **Open in browser:**
   - Go to `http://localhost:8000`
   - The visualization will automatically load the manifest and first sample

### Deployment

The project uses GitHub Actions for automatic deployment. The workflow (`.github/workflows/deploy.yml`) automatically:
- Generates the manifest file
- Deploys to GitHub Pages on changes to `index.html`, `data/all_samples/`, or `data_explorer_scripts/generate_manifest.py`

### Requirements

**For Data Explorer:**
- A modern web browser (Chrome, Firefox, Safari, Edge)
- Python 3 (for running scripts and local development)
- Normalized CSV files in `data/all_samples/` folder

**For Analysis Scripts:**
- See individual analysis README files in `analysis_scripts/` for specific requirements
- Most Python analyses require Python 3.8+ with scientific packages (pandas, numpy, scipy, etc.)
- Cell typing analysis requires R with Seurat package
