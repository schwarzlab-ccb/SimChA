out="./out"$1
echo "Writing to $out" 
echo "Plotting CN Tracks"
python3 scripts/plot_cn_tracks.py --input_folder $out --output_folder $out
