echo "Plotting Parent Graph"
dot -Tpng out/parent_graph.dot > out/parent_tree.png
echo "Plotting Fish Plot"
python3 scripts/fish.py out/populations.csv out/parent_tree.csv out/fish.png -S 42
python3 scripts/fish.py out/populations.csv out/parent_tree.csv out/fish_abs.png -S 42 --absolute
echo "Plotting CN Tracks"
python3 scripts/plot_cn_tracks.py --input_folder out --output_folder out --fraction
echo "Plotting Population Dynamics"
python3 scripts/plot_population_dynamics.py --input_folder out --output_folder out
