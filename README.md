# ExSeq Brain AD - Spatial Transcriptomics Visualization

Interactive visualization of spatial transcriptomics data from ExSeq brain AD analysis. This project uses static HTML files with Plotly for fast, serverless visualization on GitHub Pages.

## Data Overview

- **Total data points**: 522,723
- **Genes**: 104
- **Brain regions**: 5 (CA3, DG, SM, inner_DG, under_DG)
- **Fields of view (FOV)**: 19

## Visualizations

Three interactive visualizations are available:

1. **index.html** - 2D spatial view with region-based coloring
2. **view_3d.html** - 3D spatial view including z-coordinate depth
3. **dashboard.html** - Combined dashboard with 2D/3D views, statistics table, and gene count analysis

### Features

- Interactive zoom, pan, and rotation (3D)
- **Gene filtering**: Click legend items to show/hide specific genes in 2D view
- Region-based coloring for easy identification in dashboard views
- Larger, more visible point sizes (size 5 for genes, size 3-4 for regions)
- Hover tooltips showing gene, region, cell ID, FOV, and coordinates
- Statistics table showing cell counts, unique genes, and point counts per region
- Cell count bar chart comparing regions (instead of redundant gene counts)

## Installation & Setup

### Prerequisites

- Python 3.10 or higher
- [uv](https://github.com/astral-sh/uv) package manager

### Installation

1. Install dependencies:
```bash
uv pip install -e .
```

2. Generate visualizations:
```bash
python generate_visualization.py
```

This will create:
- `index.html` - 2D visualization
- `view_3d.html` - 3D visualization
- `dashboard.html` - Combined dashboard
- `data/overview.json` - Downsampled data (52K points)
- `data/stats.json` - Region statistics

## Usage

### Local Development

Simply open the HTML files in your browser:
```bash
open index.html
```

### GitHub Pages Deployment

#### Automatic Deployment (Recommended)

The repository includes a GitHub Actions workflow that automatically:
- Generates visualizations when data is updated
- Deploys to GitHub Pages

To set up:
1. Push the repository to GitHub
2. Enable GitHub Pages in repository settings
3. Select "GitHub Actions" as the publishing source

#### Manual Deployment

1. Commit the generated HTML files to your repository
2. Enable GitHub Pages in repository settings
3. Select branch and `/` (root) as source
4. Access at: `https://username.github.io/repository-name/`

## Data Structure

The visualization processes data from `fem3_WTE1_B_R_regions_genes.csv`:

- **region**: Brain region (CA3, DG, SM, inner_DG, under_DG)
- **x_coordinate, y_coordinate, z_coordinate**: Spatial coordinates
- **gene**: Gene name (104 unique genes)
- **cell_id**: Cell identifier
- **fov**: Field of view identifier
- **region_area**: Area of the region
- **region_proportion**: Proportion of total area

## Performance Optimization

- **Downsampling**: Displays every 10th point for fast initial loading (52K points from 522K)
- **Progressive detail**: Each view samples ~6-8K points for smooth rendering
- **CDN delivery**: Plotly.js loaded from CDN to reduce file size

## Regeneration

To regenerate visualizations after data updates:

```bash
python generate_visualization.py
```

## Repository Structure

```
explore-ExSeq-brain-AD/
├── index.html                    # Main 2D visualization
├── view_3d.html                 # 3D visualization
├── dashboard.html                # Combined dashboard
├── fem3_WTE1_B_R_regions_genes.csv  # Source data (63MB)
├── generate_visualization.py    # Generation script
├── pyproject.toml               # uv package configuration
├── README.md                    # This file
├── LICENSE                      # License file
├── .gitignore                   # Git ignore rules
├── .nojekyll                    # GitHub Pages configuration
├── .github/
│   └── workflows/
│       └── deploy.yml          # GitHub Actions workflow
└── data/
    ├── overview.json            # Downsampled data (13MB)
    └── stats.json              # Region statistics
```

## Files to Commit

For GitHub Pages deployment, commit:
- All `*.html` files
- `data/` directory with JSON files
- Source CSV file (or exclude if too large)
- `README.md`, `LICENSE`, `pyproject.toml`
- `.nojekyll` file (important for GitHub Pages)
- `.github/workflows/` directory

See `.gitignore` for Python artifacts that should NOT be committed.

## License

See LICENSE file for details.

