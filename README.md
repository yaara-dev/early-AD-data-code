# ExSeq Brain AD - Spatial Genomics Visualization

## 🚀 Data Explorer

**The full interactive data explorer is available via GitHub Pages:**

👉 **[https://yaara-dev.github.io/explore-ExSeq-brain-AD/](https://yaara-dev.github.io/explore-ExSeq-brain-AD/)**

Click the link above to explore the spatial genomics data with interactive 2D and 3D visualizations, filtering options, and comprehensive dashboards.

---

## Features

- **Multi-sample support:** Browse between different samples
- **Interactive 2D scatter plot:** Visualize spatial data with x/y coordinates
- **Interactive 3D scatter plot:** Explore data in three dimensions with rotation and zoom
- **Comprehensive dashboard:** View statistics, heatmaps, distributions, and more
- **Filtering:** Filter by region, gene, cell type, and Z-slice
- **Statistics:** View total records, unique genes, regions, and visible points
- **Tooltips:** Hover over points to see detailed information


## Samples Available

- fem2_5x_F5_B_left
- fem2_5x_F5_B_right
- fem2_WT_F3_B_left_updated
- fem3_5x_E7_A_left
- WTE1_B_L_UPDATED
- WTE1_B_R
- fem4_5x_F8_A_R
- fem4_WT_F11

All samples use normalized CSV format with 5 files per sample located in `data/all_samples/`.

---

## For Developers

The following information is only relevant if you want to modify or extend the visualization code.

### Project Structure

- `index.html` - Main visualization file (contains 2D view, 3D view, and dashboard)
- `data/csvs/manifest.json` - Sample manifest file (auto-generated)
- `data/all_samples/` - Normalized CSV files (5 files per sample)
- `scripts/` - Data generation and normalization scripts
  - `generate_manifest.py` - Generates manifest.json from CSV files
  - `add_cell_types.py` - Adds cell type information to CSV files
  - `server.py` - Local development server with cache control

### Data Format

The visualization uses a normalized CSV structure with 5 files per sample:
- `*_points.csv` - Point data (point_id, gene, Z, global_x, global_y, fov)
- `*_regions.csv` - Region information (region_id, region_name, region_area, region_proportion)
- `*_points_regions.csv` - Point-to-region mapping (point_id, region_id)
- `*_cells.csv` - Cell information (cell_id, cell_type)
- `*_points_cells.csv` - Point-to-cell mapping (point_id, cell_id)

### Local Development

1. **Generate the manifest file** (if data changes):
   ```bash
   python3 scripts/generate_manifest.py
   ```

2. **Run a local server:**
   ```bash
   python3 scripts/server.py
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
- Deploys to GitHub Pages on changes to `index.html`, `data/all_samples/`, or `scripts/generate_manifest.py`

### Requirements

- A modern web browser (Chrome, Firefox, Safari, Edge)
- Python 3 (for running scripts and local development)
- Normalized CSV files in `data/all_samples/` folder
