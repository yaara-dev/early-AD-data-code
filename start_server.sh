#!/bin/bash
# Simple script to start a local web server for the visualization
# This allows the HTML to automatically load CSV files from the data/csvs folder

echo "=========================================="
echo "  ExSeq Brain AD Visualization Server"
echo "=========================================="
echo ""
echo "Starting local web server..."
echo ""

# Detect OS and open browser
if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS
    (sleep 2 && open "http://localhost:8000") &
elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
    # Linux
    (sleep 2 && xdg-open "http://localhost:8000" 2>/dev/null) &
elif [[ "$OSTYPE" == "msys" || "$OSTYPE" == "cygwin" ]]; then
    # Windows (Git Bash)
    (sleep 2 && start "http://localhost:8000") &
fi

echo "The visualization will be available at:"
echo "  http://localhost:8000"
echo ""
echo "Your browser should open automatically in a few seconds."
echo "Press Ctrl+C to stop the server"
echo ""
echo "=========================================="
echo ""

# Try Python 3 first (most common)
if command -v python3 &> /dev/null; then
    python3 -m http.server 8000
# Fall back to Python 2
elif command -v python &> /dev/null; then
    python -m SimpleHTTPServer 8000
# Try Node.js http-server if available
elif command -v npx &> /dev/null; then
    npx http-server -p 8000
else
    echo "Error: Could not find Python or Node.js to start a web server."
    echo "Please install Python 3 or Node.js, or use the file input in the HTML."
    echo ""
    echo "To install Python 3:"
    echo "  macOS: brew install python3"
    echo "  Linux: sudo apt-get install python3"
    exit 1
fi

