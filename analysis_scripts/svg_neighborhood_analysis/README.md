# SVG Neighborhood Analysis

## Overview

This toolkit provides a computational pipeline for post-CELINA analysis of spatial gene expression patterns. The pipeline performs FDR correction on CELINA p-values, creates cell-type-specific expression matrices, and clusters significant genes to identify co-regulated gene modules within spatial neighborhoods.

## Method Description

The analysis pipeline processes CELINA (spatial gene expression analysis) results to identify and cluster spatially variable genes (SVGs) that show significant spatial patterns. The pipeline consists of three main steps:

1. **FDR Correction**: Applies multiple testing correction to CELINA p-values
2. **Cell-type Matrix Creation**: Generates cell-type-specific gene expression matrices
3. **Gene Clustering**: Clusters significant genes to identify co-regulated modules

### Key Features

- **FDR Correction**: Multiple correction strategies (per-cell-type, global, per-sample)
- **Cell-type Analysis**: Separates expression patterns by cell type
- **Gene Clustering**: Identifies co-regulated gene modules using k-means or hierarchical clustering
- **Comprehensive Output**: Generates cluster assignments, statistics, and visualizations

## Installation

### Requirements

- Python 3.8 or higher
- Required Python packages (see `requirements.txt`):
  - pandas >= 1.3.0
  - numpy >= 1.21.0
  - scipy >= 1.7.0
  - matplotlib >= 3.5.0
  - seaborn >= 0.11.0
  - statsmodels >= 0.14.0
  - scikit-learn >= 1.0.0

### Setup

1. Clone or download this repository:
```bash
git clone [repository_url]
cd svg_neighborhood_analysis
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

## Prerequisites

This pipeline requires **CELINA** (an R package) to be run first. CELINA analyzes spatial gene expression patterns and generates p-values for spatial variability. Please refer to the CELINA documentation for installation and usage instructions.

**Input Requirements**: Your CELINA output should include:
- `all_p_values_CombinedPvals.csv`: Raw p-values from CELINA (genes × cell types)
- `gene_cell_matrix.csv`: Gene expression matrix (genes × cells)
- `RNA_with_cells.csv`: RNA location data with cell type annotations (for cell-type matrix creation)

## Complete Analysis Workflow

The post-CELINA analysis pipeline consists of three main steps:

### Step 1: Apply FDR Correction

Apply FDR correction to CELINA p-values:

```bash
python scripts/apply_fdr_correction.py \
    --input-dir <sample_directory> \
    --method fdr_bh \
    --strategy global \
    --output-suffix _FDR_global_bonferroni
```

This generates `all_p_values_FDR_global_bonferroni.csv` with corrected p-values.

### Step 2: Create Cell-type Matrices

Create cell-type-specific gene expression matrices:

```bash
python scripts/create_celltype_matrices.py \
    --input-dir <sample_directory> \
    --output-dir <output_directory>/celltype_matrices
```

This creates separate matrices for each cell type in `celltype_matrices/<sample_name>/gene_cell_matrix_<cell_type>.csv`.

### Step 3: Cluster Significant Genes

Cluster significant genes to identify co-regulated modules:

```bash
python scripts/cluster_significant_genes.py \
    --input-dir <sample_directory> \
    --output-dir <output_directory>/clustering_results \
    --pvalue-threshold 0.01 \
    --method kmeans \
    --max-clusters 10
```

This generates cluster assignments, statistics, and visualizations for each cell type.

## Usage

### FDR Correction

```bash
python scripts/apply_fdr_correction.py \
    --input-dir <input_directory> \
    --method fdr_bh \
    --strategy global \
    --alpha 0.05 \
    --output-suffix _FDR_global_bonferroni
```

**Arguments**:
- `--input-dir` (required): Directory containing sample folders with CELINA results
- `--method`: FDR correction method (`fdr_bh`, `fdr_by`, `bonferroni`; default: `fdr_bh`)
- `--strategy`: Correction strategy (`per_celltype`, `global`, `per_sample`; default: `global`)
- `--alpha`: Significance level (default: 0.05)
- `--output-suffix`: Suffix for output files (default: `_FDR_global_bonferroni`)

### Create Cell-type Matrices

```bash
python scripts/create_celltype_matrices.py \
    --input-dir <input_directory> \
    --output-dir <output_directory>
```

**Arguments**:
- `--input-dir` (required): Directory containing sample folders with `gene_cell_matrix.csv` and `RNA_with_cells.csv`
- `--output-dir` (required): Output directory for cell-type matrices

### Cluster Significant Genes

```bash
python scripts/cluster_significant_genes.py \
    --input-dir <input_directory> \
    --output-dir <output_directory> \
    --pvalue-threshold 0.01 \
    --method kmeans \
    --max-clusters 10
```

**Arguments**:
- `--input-dir` (required): Directory containing sample folders with FDR-corrected p-values and cell-type matrices
- `--output-dir` (required): Output directory for clustering results
- `--pvalue-threshold`: P-value threshold for significant genes (default: 0.01)
- `--method`: Clustering method (`kmeans`, `hierarchical`; default: `kmeans`)
- `--max-clusters`: Maximum number of clusters to test (default: 10)
- `--min-genes`: Minimum number of genes required for clustering (default: 5)
- `--interactive`: Interactive mode to choose number of clusters (optional)

## Example Workflow

Using the provided example data:

```bash
# Step 1: FDR Correction
python scripts/apply_fdr_correction.py \
    --input-dir example_data \
    --method fdr_bh \
    --strategy global \
    --output-suffix _FDR_global_bonferroni

# Step 2: Create Cell-type Matrices (requires RNA_with_cells.csv)
# Note: This step requires RNA_with_cells.csv which is not included in example_data
# For full workflow, ensure your sample directory contains this file

# Step 3: Cluster Significant Genes
python scripts/cluster_significant_genes.py \
    --input-dir example_data \
    --output-dir example_output/clustering_results \
    --pvalue-threshold 0.01
```

## Input Data Format

### CELINA Output Files

1. **all_p_values_CombinedPvals.csv**: CSV file with genes as rows and cell types as columns, containing raw p-values from CELINA analysis.

2. **gene_cell_matrix.csv**: CSV file with genes as rows and cells as columns, containing gene expression counts.

3. **RNA_with_cells.csv** (for cell-type matrices): CSV file containing RNA location data with columns:
   - `cell_id`: Cell identifier
   - `cell_type`: Cell type annotation
   - Additional columns (x, y, z, gene, etc.)

### Example Input Format

**all_p_values_CombinedPvals.csv**:
```csv
gene,activated microglia,astrocytes,endothelial,excitatory neurons,inhibitory neurons
ACTB,8.99e-21,9.03e-08,4.08e-09,1.13e-21,0.0093
APP,6.60e-09,0.659,0.310,0.710,
```

**gene_cell_matrix.csv**:
```csv
gene,cell_001,cell_002,cell_003,cell_004
ACTB,95,7,0,0
APP,2,0,0,0
```

## Output Format

### FDR Correction Output

- `all_p_values_FDR_global_bonferroni.csv`: FDR-corrected p-values (same format as input)
- `fdr_correction_summary.csv`: Summary statistics for each sample

### Cell-type Matrices Output

- `celltype_matrices/<sample_name>/gene_cell_matrix_<cell_type>.csv`: Expression matrix for each cell type

### Clustering Output

For each sample-cell type combination:
- `cluster_assignments.csv`: Gene-to-cluster assignments
- `cluster_statistics.csv`: Statistics for each cluster
- `elbow_curve.png`: Elbow method plot for cluster selection
- `silhouette_curve.png`: Silhouette score plot
- `cluster_heatmap.png`: Heatmap of clustered genes
- `pca_clusters.png`: PCA visualization of clusters

Overall summary:
- `clustering_summary.csv`: Summary of all clustering analyses

## File Structure

```
svg_neighborhood_analysis/
├── README.md                    # This file
├── requirements.txt             # Python dependencies
├── LICENSE                      # BSD 3-Clause License
├── .gitignore                   # Git ignore rules
│
├── scripts/                     # Analysis scripts
│   ├── apply_fdr_correction.py  # FDR correction utility
│   ├── create_celltype_matrices.py  # Cell-type matrix creation
│   └── cluster_significant_genes.py  # Gene clustering analysis
│
├── example_data/                # Example input data
│   ├── all_p_values_CombinedPvals.csv  # Example CELINA p-values
│   └── gene_cell_matrix.csv     # Example gene expression matrix
│
└── example_output/              # Example output (user-generated)
    ├── celltype_matrices/       # Cell-type matrices
    └── clustering_results/       # Clustering results
```

## Method Details

### FDR Correction Strategies

- **per_celltype**: Corrects p-values within each cell type separately
- **global**: Corrects all p-values across all cell types together
- **per_sample**: Corrects p-values across all cell types within a sample

### Clustering Methods

- **k-means**: Partitions genes into k clusters based on expression similarity
- **hierarchical**: Creates a tree-based clustering structure

Both methods use optimal cluster number selection via elbow method and silhouette analysis.

## Citation

If you use this software in your research, please cite:

[Your Publication Citation]

## Contact

For questions, issues, or contributions, please contact:
- [Your Name/Email]
- [GitHub Issues: repository_url/issues]

## License

This software is licensed under the BSD 3-Clause License. See LICENSE file for details.

## Acknowledgments

[Any acknowledgments or funding information]

