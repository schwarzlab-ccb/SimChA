out="out"$1
echo "Writing to $out" 
echo "Plotting Parent Graph"
dot -Tpng $out/parent_graph.dot > $out/parent_tree.png
