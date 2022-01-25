dot -Tjpg out/parent_graph.dot > out/parent_tree.jpg
python scripts/fish_py.py out/populations.csv out/parent_tree.csv out/fish.png
python scripts/plot_cn_tracks.py --input_folder out --output_folder out --fraction
