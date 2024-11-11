#!/bin/bash
simcha_path="SimChA"
params_file="armless.json"
data_path="data/hg19"
cnps="out/tetraploid_lin/cnps.tsv"

mkdir -p out/tetraploid_sqr
for i in {1..400}; do
  stress=$(echo "scale=2; ($i - 1) * 0.005" | bc)
  task_params_file="out/tetraploid_sqr/tmp_params.json"

  python update_params.py "$params_file" "$task_params_file" "$stress"
  dotnet run --no-build --project $simcha_path -P $cnps -O "out/tetraploid_sqr/$i" -C $task_params_file -D $data_path
done
