#!/bin/bash

# Script to run the analysis notebook with example data
# This script runs the notebook as a Python script

set -e

echo "Running Moran's I Analysis Pipeline with Example Data"
echo "======================================================"

# Change to the publication directory
cd "$(dirname "$0")/.."

# Check if required files exist
if [ ! -f "notebooks/analysis_pipeline.py" ]; then
    echo "Error: notebooks/analysis_pipeline.py not found!"
    exit 1
fi

if [ ! -d "example_output/concatenated_results" ]; then
    echo "Error: example_output/concatenated_results directory not found!"
    echo "Please run the concatenation script first or ensure example data exists."
    exit 1
fi

# Run the notebook
echo ""
echo "Starting analysis..."
echo ""

python notebooks/analysis_pipeline.py

echo ""
echo "Analysis complete!"
echo "Check the output directories for results."

