# Cell typing (Seurat)

This repository contains a single R script that performs clustering and cell type annotation using a marker catalog,
and writes figures + tables to an output folder.

## Quick start (input example)
From the repository root:

```bash
Rscript code/celltyping.R
```

Outputs are written to:
- `output/` (default)

## Files
- `code/celltyping.R` — main pipeline + an input example block (runs when executed via `Rscript`)
- `input_example/input_counts.csv` — example counts matrix (genes x cells)
- `input_example/input_markers.csv` — example marker catalog (Marker, cell_type)
- `output/` — default output folder (created automatically; included for convenience)

## Viewing outputs
After running, open the PDF files in `output/`:
- `umap__input_example.pdf`
- `PCA_QC__input_example.pdf`
- `support_hist__input_example.pdf`
- `support_density__input_example.pdf`

and the tables:
- `cell_type__input_example.csv`
- `ignored_summary.csv`

## Running on your own data
In R (RStudio Console):

```r
source("code/celltyping.R")

res <- run_celltyping_pipeline(
  counts_csv  = "PATH/TO/YOUR_COUNTS.csv",
  markers_csv = "PATH/TO/YOUR_MARKERS.csv",
  out_dir     = "output",
  sample_id   = "my_sample",
  run_jackstraw = FALSE
)
```

### Windows paths
Use forward slashes:

```r
counts_csv = "C:/Users/.../my_counts.csv"
```

(or escape backslashes: `C:\\Users\\...`).
