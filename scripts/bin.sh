#!/bin/bash

threads=24

dir=$1
input="${dir}/copynumbers.tsv"
common_input="$input --samples ${dir}/clones.tsv"
#cns bin $common_input --out "${dir}/bin_10MB.tsv" --bins 10000000 --remove gaps --filter 1000000 --threads $threads --verbose
cns bin $common_input --out "${dir}/bin_3MB.tsv" --bins 3000000 --remove gaps --filter 300000 --threads $threads --verbose
#cns bin $common_input --out "${dir}/bin_1MB.tsv" --bins 1000000 --remove gaps --filter 100000 --threads $threads --verbose
#cns bin $common_input --out "${dir}/bin_300KB.tsv" --bins 300000 --remove gaps --filter 30000 --threads $threads --verbose
