# ExSeq Brain AD - Spatial Genomics Visualization

Interactive serverless HTML visualization for exploring spatial genomics data from multiple samples.

## Quick Start

### Option 1: Automatic Loading (Recommended)

1. **Generate the manifest file** (first time only):
   ```bash
   python3 generate_manifest.py
   ```
   This creates `data/csvs/manifest.json` listing all available CSV files.

2. **Run the startup script:**
   - **Mac/Linux:** Double-click `start_server.sh` or run `./start_server.sh` in terminal
   - **Windows:** Double-click `start_server.bat`
   - **Or manually:** Run `python3 -m http.server 8000` in the project directory

3. **Open in browser:**
   - Go to `http://localhost:8000`
   - The visualization will automatically load the manifest and first sample
   - Use the dropdown to switch between samples

### Option 2: Manual File Selection

1. **Open `index.html` directly** in your browser (double-click the file)
2. **Use the file input** to browse and select a CSV file from the `data/csvs/` folder
3. The visualization will load automatically

**Note:** When opening directly (file://), the manifest won't load due to browser security, so use the file upload option.

## Features

- **Multi-sample support:** Browse between 8 different samples
- **Interactive 2D scatter plot:** Visualize spatial data with x/y coordinates
- **Filtering:** Filter by region, gene, cell type, and Z-slice
- **Statistics:** View total records, unique genes, regions, and visible points
- **Tooltips:** Hover over points to see detailed information
- **Column normalization:** Automatically handles different CSV formats

## Samples Available

- fem2_5x_F5_B_left
- fem2_5x_F5_B_right
- fem2_WT_F3_B_left
- fem3_5x_E7_A_left
- fem3_WTE1_B_L
- fem3_WTE1_B_R
- fem4_5x_F8_A_R
- fem4_WT_F11

All samples are located in the `data/csvs/` folder.

## Requirements

- A modern web browser (Chrome, Firefox, Safari, Edge)
- Python 3 (for automatic loading via server) - or use manual file selection
- CSV files in `data/csvs/` folder

## Troubleshooting

**Files don't load automatically:**
- Make sure you're using the web server (Option 1) instead of opening the file directly
- Or use the file input to manually select a CSV file

**Server won't start:**
- Make sure Python 3 is installed: `python3 --version`
- Or use the manual file selection method (Option 2)

## Data Format

The visualization supports CSV files with the following columns:
- `region` or `region_name`: Brain region
- `gene`: Gene name
- `x_coordinate` or `global_x`: X coordinate
- `y_coordinate` or `global_y`: Y coordinate
- `z_coordinate` or `Z`: Z coordinate
- `cell_id` or `cell`: Cell identifier
- `cell_type`: Cell type (optional)
- `fov`: Field of view

The visualization automatically normalizes column names to handle different formats.
