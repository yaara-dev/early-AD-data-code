# Analysis Scripts

This directory contains computational analysis pipelines for spatial transcriptomics data from the ExSeq Brain AD study. Each analysis pipeline is self-contained with its own documentation, requirements, and example data.

## Overview

The analysis scripts provide four complementary approaches to analyzing spatial gene expression patterns:

1. **Cell Typing** - Automated cell type annotation using Seurat and marker gene catalogs
2. **Moran's I Spatial Autocorrelation** - Spatial pattern analysis using grid-based Moran's I statistics
3. **RNA Velocity** - Analysis of RNA velocity and spatial cell-state transitions
4. **SVG Neighborhood Analysis** - Post-CELINA analysis of spatially variable genes and co-clustered genes

## Available Analyses

### 1. Cell Typing (`cell_typing/`)

**Purpose**: Automated cell type annotation from single-cell expression data using marker gene catalogs.

**Key Features**:
- Seurat-based clustering and dimensionality reduction
- Marker catalog-driven cell type assignment
- Quality control visualizations (PCA, UMAP, support histograms)
- Handles multiple cell types simultaneously

**Language**: R  
**Main Script**: `code/celltyping.R`  
**Documentation**: See [cell_typing/README.md](cell_typing/README.md)

---

### 2. Moran's I Spatial Autocorrelation (`morans_i/`)

**Purpose**: Quantify spatial clustering and dispersion patterns of gene expression using spatial autocorrelation statistics.

**Key Features**:
- FOV-specific analysis with 3D grid-based approach
- Binary_6 weight scheme for 3D spatial adjacency
- Per-FOV and combined regional results
- Statistical comparison between experimental groups
- Quantile normalization and FDR correction

**Language**: Python  
**Main Scripts**: `src/pipelines/fov_grid_morans.py`, `notebooks/analysis_pipeline.py`  
**Documentation**: See [morans_i/README.md](morans_i/README.md)

---

### 3. RNA Velocity (`rna_velocity/`)

**Purpose**: Analyze RNA velocity and spatial cell-state transitions to understand directional changes in gene expression.

**Key Features**:
- Per-cell normalization of spliced/unspliced RNA counts
- KNN pooling in expression/PCA space
- Gene-wise γ parameter estimation
- Distance-stratified analysis of cell-state changes
- Region-stratified analysis within anatomical regions
- Permutation testing for phase analysis

**Language**: Python  
**Main Script**: `src/rna_velocity_all_tissues.py`  
**Documentation**: See [rna_velocity/README.md](rna_velocity/README.md)

---

### 4. SVG Neighborhood Analysis (`svg_neighborhood_analysis/`)

**Purpose**: Post-CELINA analysis to identify and cluster spatially variable genes (SVGs) within spatial neighborhoods.

**Key Features**:
- FDR correction on CELINA p-values (multiple strategies)
- Cell-type-specific expression matrix creation
- Gene clustering to identify co-regulated modules
- K-means and hierarchical clustering methods
- Comprehensive visualizations and statistics

**Language**: Python  
**Main Scripts**: `scripts/apply_fdr_correction.py`, `scripts/create_celltype_matrices.py`, `scripts/cluster_significant_genes.py`  
**Documentation**: See [svg_neighborhood_analysis/README.md](svg_neighborhood_analysis/README.md)

---

## Quick Start

Each analysis has its own detailed README with installation instructions, input requirements, and usage examples. Navigate to the individual analysis directories for:

- Installation and requirements
- Input data format specifications
- Step-by-step usage instructions
- Example workflows
- Output format descriptions

## Common Requirements

Most Python-based analyses require:
- Python 3.8 or higher
- Common scientific Python packages (pandas, numpy, scipy, matplotlib, seaborn)
- See individual `requirements.txt` files for specific dependencies

The R-based cell typing analysis requires:
- R (with Seurat package)

## Directory Structure

```
analysis_scripts/
├── README.md                    # This file
├── cell_typing/                 # Cell type annotation pipeline
│   ├── README.md
│   ├── code/
│   ├── input_example/
│   └── output/
├── morans_i/                    # Spatial autocorrelation analysis
│   ├── README.md
│   ├── src/
│   ├── notebooks/
│   ├── scripts/
│   ├── example_data/
│   └── example_output/
├── rna_velocity/                # RNA velocity analysis
│   ├── README.md
│   ├── src/
│   ├── data/
│   └── RESULTS/
└── svg_neighborhood_analysis/   # SVG neighborhood analysis
    ├── README.md
    ├── scripts/
    ├── example_data/
    └── example_output/
```

## Getting Help

For detailed information about each analysis:
1. Navigate to the specific analysis directory
2. Read the individual README.md file
3. Check example data and output directories for format references
4. Review the source code comments for implementation details

