#!/bin/bash
run_cmd="dotnet run --project SimChA -T in/refphase.csv -C my_params.json"

seg_cmd="python scripts/consistent_segmentation.py -i out/copynumbers.tsv -o out/segments.tsv"

filename="segment_distribution.tsv"
rm $filename
touch $filename

#perl -le 'print "sample_id", "\t", "mean_seg_l", "\t", "n_segs", "\t", "ploidy"' > $filename

for i in $(seq 1 1200):
do
  echo "Event $i"
  $run_cmd
  #sed -i "1s/chrom/chr/" out/copynumbers.tsv
  #$seg_cmd
  # Run MSAI finder
  python scripts/count_segments.py -i out/copynumbers.tsv -p out/clones.tsv -o $filename -e $i
done
