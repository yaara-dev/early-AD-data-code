# FOV Grid-Based Moran's I Spatial Autocorrelation Analysis

## Overview

This toolkit provides a computational pipeline for calculating Moran's I spatial autocorrelation statistics on spatial transcriptomics data using a Field-of-View (FOV) grid-based approach with binary_6 weight scheme. The method analyzes spatial patterns of gene expression by partitioning data into spatial grid cells and computing spatial autocorrelation for each FOV separately.

## Method Description

The analysis employs a **3D grid-based approach** where spatial transcriptomics data points are assigned to regular 3D grid cells. For each FOV, Moran's I statistics are calculated using a **binary_6 weight scheme**, which defines spatial adjacency based on face-adjacent neighbors in 3D space (6 neighbors: ±x, ±y, ±z directions). This approach captures spatial clustering and dispersion patterns of gene expression within each FOV.

### Key Features

- **FOV-specific analysis**: Processes each Field-of-View independently
- **3D spatial gridding**: Uses 3D grid cells for spatial unit definition
- **Binary_6 adjacency**: Defines spatial weights based on 3D face-adjacent neighbors
- **Gene expression analysis**: Calculates Moran's I for multiple genes simultaneously
- **Comprehensive output**: Generates both per-FOV and combined regional results

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

### Setup

1. Clone or download this repository:
```bash
git clone [repository_url]
cd morans_i_publication
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

## Complete Analysis Workflow

The Moran's I analysis pipeline consists of three main steps:

### Step 1: Run FOV Grid-Based Analysis

First, run the FOV grid-based Moran's I analysis on your spatial transcriptomics data:

```bash
python src/pipelines/fov_grid_morans.py \
    --region-file <input_file.csv> \
    --output-dir <output_directory> \
    --sample-id <sample_identifier> \
    --weight-scheme binary_6 \
    --grid-size 35.0
```

This generates per-FOV and combined regional results for each sample.

### Step 2: Concatenate Regional Results

After running the analysis for multiple samples, concatenate the results across samples for each region:

```bash
python scripts/concatenate_regional_csvs.py <output_directory>/fovs/grid_based/binary_6 \
    --output <output_directory>/concatenated_results
```

This script combines results from different samples side-by-side, adding sample name prefixes to column names. The output files will be named `{region}_all_samples_morans_concatenated.csv` and `{region}_all_samples_counts_concatenated.csv`.

**Important**: Run this step before proceeding to statistical analysis. The concatenation script expects the directory structure created by Step 1, with region folders containing sample-specific CSV files.

### Step 3: Statistical Analysis

Run the analysis notebook to perform statistical tests, normalization, and visualization:

```bash
# Open the notebook in your preferred environment (VS Code, PyCharm, Jupyter)
# The notebook uses #%% cells for interactive execution
python notebooks/analysis_pipeline.py
```

Or execute interactively using the `#%%` cell markers in VS Code or PyCharm.

The analysis notebook performs:
- **Quantile normalization** of Moran's I and counts data
- **Permutation tests** comparing WT vs 5X samples (or your experimental groups)
- **FDR correction** using Benjamini-Hochberg method
- **Visualization** including heatmaps, volcano plots, and distribution histograms
- **Significant gene identification** with direction of effect analysis

**Configuration**: Before running, update the `base_path` variable in the notebook to point to your concatenated results directory (default: `example_output/concatenated_results`).

## Usage

### Basic Command

```bash
python src/pipelines/fov_grid_morans.py \
    --region-file <input_file.csv> \
    --output-dir <output_directory> \
    --sample-id <sample_identifier> \
    --weight-scheme binary_6 \
    --grid-size 35.0
```

### Required Arguments

- `--region-file`: Path to input CSV file containing spatial transcriptomics data with region and FOV assignments
- `--output-dir`: Directory where results will be saved
- `--sample-id`: Sample identifier used for output file naming

### Optional Arguments

- `--grid-size`: Grid cell size in microns (default: 35.0)
- `--weight-scheme`: Spatial weight scheme - use `binary_6` for 3D adjacency (default: binary_6)
- `--min-points-per-fov`: Minimum number of data points required per FOV (default: 50)

### Example

Using the provided example script:

```bash
bash scripts/run_fov_grid_binary6.sh
```

Or run directly:

```bash
python src/pipelines/fov_grid_morans.py \
    --region-file example_data/sample_data.csv \
    --output-dir example_output \
    --sample-id EXAMPLE_SAMPLE \
    --weight-scheme binary_6 \
    --grid-size 35.0 \
    --min-points-per-fov 10
```

## Input Data Format

The input CSV file must contain the following columns:

- `region_name`: Brain region or spatial region identifier (e.g., "CA1", "SM", "DG")
- `gene`: Gene name (e.g., "ACTB", "GAPDH", "SNAP25")
- `global_x`: X coordinate in microns
- `global_y`: Y coordinate in microns
- `Z`: Z coordinate in microns
- `cell`: Cell identifier
- `fov`: Field-of-View identifier (e.g., "FOV_01", "FOV_02")

### Example Input Format

```csv
region_name,gene,global_x,global_y,Z,cell,fov
CA1,ACTB,1250.5,2450.2,12.3,C01,FOV_01
CA1,GAPDH,1350.8,2550.1,11.8,C02,FOV_01
SM,SNAP25,1400.9,2600.2,12.4,C03,FOV_02
```

## Output Format

The analysis generates two types of output files:

### 1. Per-FOV Results

Individual CSV files for each FOV containing:
- `GeneNames`: Gene identifier
- `morans_i`: Moran's I statistic value
- `coverage_um2`: FOV coverage area in square microns
- `point_count`: Number of data points in the FOV
- `size_flag`: FOV size classification (FULL, PARTIAL, SMALL)

**Location**: `output_dir/fovs/grid_based/binary_6/<sample_id>/<region_name>_fovs/<fov_id>_morans.csv`

### 2. Combined Regional Results

Combined CSV files per region containing:
- `GeneNames`: Gene identifier
- Columns for each FOV: `morans_i` values for each FOV in the region

**Location**: `output_dir/fovs/grid_based/binary_6/<region_name>/<sample_id>_<region_name>_morans.csv`

### Output Interpretation

- **Positive Moran's I**: Indicates spatial clustering (similar expression values are spatially adjacent)
- **Negative Moran's I**: Indicates spatial dispersion (dissimilar expression values are spatially adjacent)
- **Moran's I ≈ 0**: Indicates random spatial distribution
- **N/A**: Indicates insufficient data for calculation (zero gene counts)

## Method Details

### Spatial Gridding

Data points are assigned to 3D grid cells based on their spatial coordinates. The grid cell size (default: 35μm) determines the spatial resolution of the analysis.

### Binary_6 Weight Scheme

The binary_6 weight scheme creates a binary adjacency matrix where:
- Weight = 1 for face-adjacent grid cells in 3D space (±x, ±y, ±z directions)
- Weight = 0 for all other pairs
- Diagonal = 0 (no self-influence)

This approach captures immediate spatial neighbors in 3D space, making it suitable for analyzing 3D spatial transcriptomics data.

### Moran's I Calculation

For each gene, Moran's I is calculated as:

```
I = (n/W) × ΣᵢΣⱼ wᵢⱼ(xᵢ - x̄)(xⱼ - x̄) / Σᵢ(xᵢ - x̄)²
```

Where:
- `n`: Number of spatial units (grid cells)
- `W`: Sum of all spatial weights
- `wᵢⱼ`: Spatial weight between units i and j
- `xᵢ`: Expression value at unit i
- `x̄`: Mean expression value

## File Structure

```
morans_i_publication/
├── README.md                    # This file
├── requirements.txt             # Python dependencies
├── LICENSE                      # BSD 3-Clause License
│
├── src/                        # Source code
│   ├── algorithms/
│   │   └── morans_i.py         # Core Moran's I algorithm
│   ├── pipelines/
│   │   ├── grid_based_morans.py      # Grid-based pipeline
│   │   └── fov_grid_morans.py        # FOV grid analysis (main entry point)
│   └── utils/
│       ├── spatial.py           # Spatial utilities
│       └── io.py                # I/O utilities
│
├── example_data/               # Example input data
│   └── sample_data.csv          # Sample dataset
│
├── notebooks/                  # Analysis notebooks
│   └── analysis_pipeline.py       # Statistical analysis pipeline
│
├── example_data/               # Example input data
│   └── sample_data.csv          # Sample dataset
│
├── example_output/             # Example results
│   ├── fov_grid_binary6_results/  # FOV analysis outputs
│   └── concatenated_results/      # Concatenated results (input for analysis notebook)
│
└── scripts/                    # Utility scripts
    ├── run_fov_grid_binary6.sh    # Example execution script
    └── concatenate_regional_csvs.py  # Concatenation script
```

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

