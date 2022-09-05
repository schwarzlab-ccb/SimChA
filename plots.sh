out="./out"$1
echo "Writing to $out" 
echo "Plotting Parent Graph"
dot -Tpng $out/parent_graph.dot > $out/parent_tree.png
echo "Plotting CN Tracks"
python3 scripts/plot_cn_tracks.py --input_folder $out --output_folder $out --fraction
