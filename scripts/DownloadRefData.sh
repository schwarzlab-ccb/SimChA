#!/usr/bin/bash

# Download the hg19 dataset
baseAPI="https://storage.googleapis.com/genomics-public-data/references/GRCh37/"
chrs=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")
# hg19 stores chromosomes individually.
for i in ${!chrs[@]}; do
  # Shorthand for the chromosome ID
  chr=${chrs[$i]}
  # Download and unzip
  wget $baseAPI$chr.fa.gz
  gzip -d $chr.fa.gz
  # Append the chromosomes to the final reference FASTA file
  cat $chr.fa >> genome.fa
  # Delete individual chromosome FASTA
  rm $chr.fa
done

# Fix chromosome headers from NCBI format to UCSC format
# Convert ">gi|...|gb|...| Homo sapiens chromosome N,..." to ">chrN"
sed -E 's/>gi\|[0-9]+\|gb\|[^|]+\| Homo sapiens chromosome ([0-9XY]+),.*/>\1/' genome.fa > genome.tmp
sed 's/^>\([0-9XY]\)/>chr\1/' genome.tmp > genome.fa
rm genome.tmp

mv genome.fa ../data/hg19/
 
# Download hg38 dataset
wget https://storage.googleapis.com/genomics-public-data/references/GRCh38_Verily/GRCh38_Verily_v1.genome.fa -O genome.fa
mv genome.fa ../data/hg38/