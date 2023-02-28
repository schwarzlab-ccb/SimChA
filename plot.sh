#!/bin/bash
set -euo pipefail

dir="./out"

echo "Selected data folder $dir"
echo "Creating consistent segmentation"
python3 scripts/consistent_segmentation.py -i $dir/copynumbers.out -o $dir/copynumbers_consistent.out

echo "Plotting CN Tracks"
# check if tree file exists
echo "$dir/clone_tree.new"
if [ -f $dir/clone_tree.new ]; then
    echo "Found .new tree, plotting it"
    python3 scripts/plot_copynumbers.py --input $dir/copynumbers_consistent.out --tree $dir/clone_tree.new --output-folder $dir --type heatmap
elif  [ -f $dir/clone_tree.dot ]; then
    echo "Found .dot tree, plotting it"
    python3 scripts/plot_copynumbers.py --input $dir/copynumbers_consistent.out --tree $dir/clone_tree.dot --output-folder $dir --type heatmap
else
    echo "No tree found"
    python3 scripts/plot_copynumbers.py --input $dir/copynumbers_consistent.out --output-folder $dir --type heatmap
fi
