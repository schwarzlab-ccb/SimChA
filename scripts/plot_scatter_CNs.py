#!/usr/bin/env python3

import pandas as pd
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from os.path import join
from utils import load_dataset, hg19_chr_lengths, chromosome_colors, chromosome_names, hg38_chr_lengths

# Define a function to convert a row to a list
def row_to_list(row, hg19_chr_cum_starts, hg38_chr_cum_starts, step_size=1000000):
    # Initialize an empty list
    res = []
    
    start = row['start'] + hg19_chr_cum_starts[row['chrom']]
    if start % step_size != 0:
        start = start - (start % step_size) + step_size
    end = row['end'] + hg19_chr_cum_starts[row['chrom']]

    # Loop over each position between start and end, divisible by step
    for pos in range(start, end, step_size):
        cn = row['cn']
        res.append((pos, cn))
    
    # Return the list
    return res

def sample_to_CNs(sample, step_size, hg19_chr_cum_starts, hg38_chr_cum_starts):
    rows = sample.apply(lambda x: row_to_list(x, hg19_chr_cum_starts, hg38_chr_cum_starts, step_size), axis=1)
    filtered = rows[rows.apply(lambda x: len(x) > 0)]
    return pd.Series(dict(np.concatenate(filtered.values)))

def genome_length_init():
    # Initialize the cumulative start position variable
    cum_start = 0
    hg19_chr_cum_starts = {}
    # Loop over each chromosome in the chromosome lengths dictionary
    for chrom, length in hg19_chr_lengths.items():
        hg19_chr_cum_starts[chrom] = cum_start
        cum_start += length
    cum_start = 0
    hg38_chr_cum_starts = {}
    # Loop over each chromosome in the chromosome lengths dictionary
    for chrom, length in hg38_chr_lengths.items():
        hg38_chr_cum_starts[chrom] = cum_start
        cum_start += length    

    # total length of the genome
    return hg19_chr_cum_starts, hg38_chr_cum_starts

def plot_scatter_CNs(data, use_hg19=True, output="../out"):
    # Initialize the cumulative start position variable
    hg19_chr_cum_starts, hg38_chr_cum_starts = genome_length_init() 
    cns = pd.read_csv(join(data,"copynumbers.tsv"), index_col=0, sep="\t")
    cns['cn'] = cns['cn_a'] + cns['cn_b']
    # list of unique indices in cns
    samples = cns.index.unique()
    positions = {}
    step_size = 1000000
    for sample in samples:
        positions[sample] = sample_to_CNs(cns.loc[sample, :], step_size, hg19_chr_cum_starts, hg38_chr_cum_starts) 
    df_CNs = pd.DataFrame(positions)

    # print(df_CNs)
    # Create a figure and ax[0]is object
    fig, ax = plt.subplots(1, figsize=(32, 9))

    max_cn = 11
    # Set the y-ax[0]is limits
    ax.set_ylim(0, max_cn)

    # Loop over each chromosome in the hg19_chrom_lengths dictionary
    x_pos = 0
    for chrom, length in hg19_chr_lengths.items():
        # Get the color for the current chromosome
        color = chromosome_colors[chrom]
        
        # Add a rectangle to the plot with the appropriate color and width
        rect = Rectangle((x_pos, 0), length, max_cn, color=color, alpha=0.2)
        ax.add_patch(rect)
        
        # Update the x position for the next chromosome
        x_pos += length

    for sample in df_CNs.columns:
        sample_coords = df_CNs[sample]
        sample_count = len(sample_coords)
        jitter = np.random.normal(0, 0.1, sample_count)
        ax.scatter(df_CNs.index, sample_coords+jitter, s=.1, alpha=0.1, color="gray")

    ax.plot(df_CNs.index, df_CNs.mean(axis=1), color="red", linewidth=1.5, alpha=0.7)
    # plot the a band representing the standard deviation
    ax.fill_between(df_CNs.index, df_CNs.mean(axis=1) - df_CNs.std(axis=1), df_CNs.mean(axis=1) + df_CNs.std(axis=1), color="red", alpha=0.15)

    # Set the x-ax[0]is limits to match the total genome length
    ax.set_xlim(0, x_pos)
    # Set the x-ax[0]is ticks to be at the middle of each chromosome
    ax.set_xticks([hg19_chr_cum_starts[chrom] + (length / 2) for chrom, length in hg19_chr_lengths.items()])
    # set tick labels to be chromosome names
    ax.set_xticklabels(list(hg19_chr_cum_starts.keys()))

    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Copy number")
    ax.set_title("Copy number variation")

    # plot a line at y=2
    ax.plot(df_CNs.index, [2] * len(df_CNs.index), color="black", linewidth=1, alpha=0.7)

    fig.savefig(f"{output}/scatter_CNs.png", dpi=150, bbox_inches="tight")

def plot_scatter_CNs():
    datasets = {"PCAWG":"../out/PCAWG/", "PCAWG, Filtered, 95%": "../out/PCAWG_filtered_95_pc", "PCAWG, Filtered, 99%":"../out/PCAWG_filtered_99_pc"}#, "Neutral": "../out", "Fitness-Driven": "../mcmc_complex"}
    ls = ["-", "-.", "--"]
    # Initialize the cumulative start position variable
    hg19_chr_cum_starts, hg38_chr_cum_starts = genome_length_init()

    x_pos = 0
    for chrom, length in hg19_chr_lengths.items():
        # Update the x position for the next chromosome
        x_pos += length
     
    fig, ax = plt.subplots(1, figsize=(32, 9))
    for i, key in enumerate(datasets.keys()):
        data = datasets[key]
        cns = pd.read_csv(join(data,"copynumbers.tsv"), index_col=0, sep="\t")
        cns['cn'] = cns['cn_a'] + cns['cn_b']
        # list of unique indices in cns
        samples = cns.index.unique()
        positions = {}
        step_size = 1000000
        for sample in samples:
            positions[sample] = sample_to_CNs(cns.loc[sample, :], step_size, hg19_chr_cum_starts, hg38_chr_cum_starts) 
        df_CNs = pd.DataFrame(positions)

        # Set the y-ax[0]is limits
        ax.set_ylim(0.4, 3.5)

        ax.plot(df_CNs.index, df_CNs.mean(axis=1), label=key, linewidth=1.5, ls = ls[i])

        # Set the x-ax[0]is limits to match the total genome length
        ax.set_xlim(0, x_pos)

        if i == 0:
            # plot a line at y=2
            ax.plot(df_CNs.index, [2] * len(df_CNs.index), color="black", linewidth=1, alpha=0.7, ls=":")
    
    # Set the x-ax[0]is ticks to be at the middle of each chromosome
    ax.set_xticks([hg19_chr_cum_starts[chrom] + (length / 2) for chrom, length in hg19_chr_lengths.items()])
    # set tick labels to be chromosome names
    ax.set_xticklabels(list(hg19_chr_cum_starts.keys()))
    ax.legend()

    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Copy number")
    ax.set_title("Copy number variation")


    fig.savefig(f"./scatter_CNs.svg", dpi=300, bbox_inches="tight")
    




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot the average scatter CN plots for a given dataset')
    parser.add_argument('-I', "--input", type=str, default=".", help='The folder to input dataset to plot')
    parser.add_argument("--hg19",action="store_true", dest="use_hg19", default=False, help="Use hg19 instead of hg38")
    parser.add_argument('-O', "--output", type=str, default=".", help='The folder to output the plots to')
    args = parser.parse_args()

    plot_scatter_CNs()
    #plot_scatter_CNs(args.input, args.use_hg19, args.output)




