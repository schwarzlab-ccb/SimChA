#!/bin/bash
cd ..
set -euo pipefail

if [ $# -gt 0 ]; then out=$1; else out="./out"; fi

echo "Selected data folder $out"
echo "Creating consistent segmentation"
python3 scripts/consistent_segmentation.py -i $out/copynumbers.tsv -o $out/copynumbers_consistent.out

echo "Plotting CN Tracks"
# check if tree file exists
if [ -f $out/parent_graph.new ]; then
    echo "Plotting with evolutionary tree"
    python3 scripts/plot_copynumbers.py --input $out/copynumbers_consistent.out --tree $out/parent_graph.new --output-folder $out --type heatmap
else
    echo "Plotting without evolutionary tree"
    python3 scripts/plot_copynumbers.py --input $out/copynumbers_consistent.out --output-folder $out --type heatmap
fi
