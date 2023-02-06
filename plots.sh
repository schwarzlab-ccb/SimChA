out="./out"$1

echo "Selected data folder $out"
echo "Creating consistent segmentation"
python3 scripts/consistent_segmentation.py -i $out/copynumbers.out -o $out/copynumbers_consistent.out

echo "Plotting CN Tracks"
# check if tree file exists
if [ -f $out/parent_graph.new ]; then
    echo "Found .new tree, plotting it"
    python3 scripts/plot_copynumbers.py --input $out/copynumbers_consistent.out --tree $out/parent_graph.new --output-folder $out --type heatmap
elif  [ -f $out/parent_graph.dot ]; then
    echo "Found .dot tree, plotting it"
    python3 scripts/plot_copynumbers.py --input $out/copynumbers_consistent.out --tree $out/parent_graph.dot --output-folder $out --type heatmap
else
    echo "No tree found"
    python3 scripts/plot_copynumbers.py --input $out/copynumbers_consistent.out --output-folder $out --type heatmap
fi
