#!/bin/bash
# Test script to verify GitHub Actions deployment workflow locally

set -e  # Exit on error

echo "=========================================="
echo "Testing GitHub Actions Deployment Workflow"
echo "=========================================="
echo ""

# Step 1: Check Python version
echo "Step 1: Checking Python version..."
python3 --version
if [ $? -ne 0 ]; then
    echo "❌ ERROR: Python 3 not found"
    exit 1
fi
echo "✓ Python found"
echo ""

# Step 2: Generate CSV manifest
echo "Step 2: Generating CSV manifest..."
python3 scripts/generate_manifest.py
if [ $? -ne 0 ]; then
    echo "❌ ERROR: Manifest generation failed"
    exit 1
fi
echo "✓ Manifest generated successfully"
echo ""

# Step 3: Verify manifest.json exists and is valid
echo "Step 3: Verifying manifest.json..."
if [ ! -f "data/csvs/manifest.json" ]; then
    echo "❌ ERROR: manifest.json not found"
    exit 1
fi

# Check if it's valid JSON
python3 -m json.tool data/csvs/manifest.json > /dev/null
if [ $? -ne 0 ]; then
    echo "❌ ERROR: manifest.json is not valid JSON"
    exit 1
fi
echo "✓ manifest.json is valid"
echo ""

# Step 4: Verify required files exist
echo "Step 4: Verifying required files..."
REQUIRED_FILES=("index.html" "data/csvs/manifest.json")
for file in "${REQUIRED_FILES[@]}"; do
    if [ ! -f "$file" ]; then
        echo "❌ ERROR: Required file not found: $file"
        exit 1
    fi
    echo "✓ Found: $file"
done
echo ""

# Step 5: Count CSV files
echo "Step 5: Checking CSV files..."
CSV_COUNT=$(find data/all_samples -name "*_points.csv" | wc -l | tr -d ' ')
echo "Found $CSV_COUNT sample files (*_points.csv)"
if [ "$CSV_COUNT" -lt 1 ]; then
    echo "⚠️  WARNING: No CSV files found"
else
    echo "✓ CSV files found"
fi
echo ""

echo "=========================================="
echo "✅ All checks passed! Deployment should work."
echo "=========================================="
echo ""
echo "Note: This test doesn't actually deploy to GitHub Pages."
echo "The actual deployment happens on GitHub Actions when you push to main."

