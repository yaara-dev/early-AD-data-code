#!/bin/bash
# Example script to run FOV grid-based Moran's I analysis with binary_6 weight scheme

# Set paths
REGION_FILE="example_data/sample_data.csv"
OUTPUT_DIR="example_output"
SAMPLE_ID="EXAMPLE_SAMPLE"
GRID_SIZE=35.0
WEIGHT_SCHEME="binary_6"

# Run the analysis
python src/pipelines/fov_grid_morans.py \
    --region-file "$REGION_FILE" \
    --output-dir "$OUTPUT_DIR" \
    --sample-id "$SAMPLE_ID" \
    --grid-size "$GRID_SIZE" \
    --weight-scheme "$WEIGHT_SCHEME" \
    --min-points-per-fov 10

echo "Analysis complete! Results saved to $OUTPUT_DIR"

