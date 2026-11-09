# ExSeq Brain AD - Spatial Genomics Data Explorer & Analysis Scripts

This repository contains an interactive web-based data explorer, Xenium HTML viewers, and computational analysis pipelines for spatial transcriptomics data from the ExSeq Brain AD study.

---

## 🚀 Data Explorer

**The full interactive data explorer is available via GitHub Pages:**

👉 **[https://yaara-dev.github.io/early-AD-data-code/](https://yaara-dev.github.io/early-AD-data-code/)**

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

**WT:** `WT_1`, `WT_2.1`, `WT_2.2`, `WT_3`, `WT_4`, `WT_5`, `WT_6.1`, `WT_6.2`

**5xFAD:** `5xFAD_1.1`, `5xFAD_1.2`, `5xFAD_2`, `5xFAD_3`, `5xFAD_4`, `5xFAD_5`, `5xFAD_6.1`, `5xFAD_6.2`

All samples use a **separated/normalized** CSV format with 5 files per sample in `data/viewer_normalized/` (for the web viewer).

---

## 📊 Analysis Scripts

This repository includes four computational analysis pipelines for spatial transcriptomics data:

1. **Cell Typing** - Supervised STELLAR encoder label transfer for ExSeq cell-type annotation
2. **Expression, Spatial Distribution & Moran's I** - MATLAB scripts for differential expression, localization/spatial distribution, and Moran’s I
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

## 🔬 Xenium HTML viewers

Self-contained **Xenium Ranger** analysis-summary reports for the Xenium samples. Open any `.html` file in a web browser (no server required) to inspect run quality, transcript and cell metrics, and spatial summaries.

Files live in [`data/xenium_html/`](data/xenium_html/) in this repository. The same reports are also included in both Zenodo deposits ([code](https://doi.org/10.5281/zenodo.17568649) and [analysis data](https://doi.org/10.5281/zenodo.21903277)).

**WT:** `WT_5`, `WT_6_section1`, `WT_6_section2`

**5xFAD:** `5xFAD_5_section1`, `5xFAD_5_section2`, `5xFAD_6`

| File |
|------|
| `WT_5_analysis_summary.html` |
| `WT_6_section1_analysis_summary.html` |
| `WT_6_section2_analysis_summary.html` |
| `5xFAD_5_section1_analysis_summary.html` |
| `5xFAD_5_section2_analysis_summary.html` |
| `5xFAD_6_analysis_summary.html` |

These reports are for visual exploration of the Xenium runs. They are not used by the analysis scripts.

---

## Supplementary figures

- **[Fig. S2](figures/Fig.S2.pdf)** — Qualitative comparison of ExSeq spatial gene patterns in sample `WT_2` with Allen Mouse Brain Atlas in situ hybridization data. See [`figures/README.md`](figures/README.md) for a more detailed description.

---

## For Developers

The following information is only relevant if you want to modify or extend the visualization code.

### Project Structure

**Data Explorer Components:**
- `index.html` - Main visualization file (contains 2D view, 3D view, and dashboard)
- `data/csvs/manifest.json` - Sample manifest file (auto-generated)
- `data/viewer_normalized/` - Separated normalized CSVs for the web viewer (5 files per sample)
- `data/exseq/`, `data/xenium/`, `data/xenium_full_sections/` - Analysis tables (download from the [data deposit](https://doi.org/10.5281/zenodo.21903277); see [data/README.md](data/README.md))
- `data/xenium_html/` - Xenium Ranger HTML reports (in this repository, and also in both Zenodo deposits)
- `figures/` - Supplementary figures (e.g. Fig. S2 Allen Atlas comparison)
- `data_explorer_scripts/` - Data generation and normalization scripts
  - `generate_manifest.py` - Generates manifest.json from CSV files
  - `add_cell_types.py` - Adds cell type information to CSV files
  - `server.py` - Local development server with cache control

**Analysis Scripts:**
- `analysis_scripts/` - Computational analysis pipelines
  - `cell_typing/` - Supervised STELLAR encoder label transfer (Python)
  - `expression_spatial_distribution_morans_i/` - Expression, spatial distribution/localization, and Moran’s I (MATLAB)
  - `rna_velocity/` - RNA velocity analysis (Python)
  - `svg_neighborhood_analysis/` - SVG neighborhood analysis (Python)

### Zenodo archives

- **Code** (this repository release, including Xenium Ranger HTML reports): [https://doi.org/10.5281/zenodo.17568649](https://doi.org/10.5281/zenodo.17568649)
- **Analysis data** (ExSeq / Xenium tables, and the same Xenium Ranger HTML reports): [https://doi.org/10.5281/zenodo.21903277](https://doi.org/10.5281/zenodo.21903277)

Large ExSeq/Xenium tables are **not** stored in git. To run the MATLAB analyses:

1. Clone this repository (or download the [code Zenodo archive](https://doi.org/10.5281/zenodo.17568649)).
2. Download the [data deposit](https://doi.org/10.5281/zenodo.21903277).
3. Unpack so these folders are populated:
   - `data/exseq/`
   - `data/xenium/`
   - `data/xenium_full_sections/`
4. Run the scripts under `analysis_scripts/expression_spatial_distribution_morans_i/` — they already point at those paths.

Xenium Ranger HTML reports are in `data/xenium_html/` in this repository (open in a browser). They are also included in both Zenodo deposits. They are not used by the MATLAB scripts.

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
- Deploys to GitHub Pages on changes to `index.html`, `data/viewer_normalized/`, or `data_explorer_scripts/generate_manifest.py`

### Requirements

**For Data Explorer:**
- A modern web browser (Chrome, Firefox, Safari, Edge)
- Python 3 (for running scripts and local development)
- Separated normalized CSV files in `data/viewer_normalized/` folder

**For Analysis Scripts:**
- See individual analysis README files in `analysis_scripts/` for specific requirements
- Most Python analyses require Python 3.8+ with scientific packages (pandas, numpy, scipy, etc.)
- Cell typing (`cell_typing/`) requires Python packages listed in `analysis_scripts/cell_typing/environment.yml`
- Expression / spatial distribution / Moran’s I (`expression_spatial_distribution_morans_i/`) requires MATLAB (see that folder’s README for toolboxes)
