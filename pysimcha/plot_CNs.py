#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
import pandas as pd
import math
import sys
sys.path.append('../scripts')
from utils import chromosome_colors

# a plot with rectangles between start-end on a line given by the copy number
def plot_hap(ax, chr_data, chr, hap):    
    if (len(chr_data) == 0):
        return
    for index, row in chr_data.iterrows():
        start = row["start"] / 1000000
        end = row["end"] / 1000000
        height = .2
        width = end - start
        if hap == "a":
            y_pos = row["cn_a"] - .1  
            hatch = "///"
            alpha = 0.5
        elif hap == "b":
            y_pos = row["cn_b"] - .1  
            hatch = "\\\\\\"            
            alpha = 0.5
        else:
            y_pos = row["cn_both"] - .1        
            hatch = None            
            alpha = 1
        fc = chromosome_colors[chr]
        rect = mpatches.Rectangle((start, y_pos), width, height, fc=fc, alpha=alpha, hatch=hatch)   
        ax.add_patch(rect)           

def plot_chr(ax, data, chr, start = 0, end = 0, join_haps = False):
    ax.set_ylabel(f"CN {chr}")
    chr_data = data.loc[data["chrom"] == chr]
    haps = ["both"] if join_haps else ["a", "b"]
    max_end = chr_data.iloc[-1]["end"] if len(chr_data) > 0 else 1
    max_cna = int(math.ceil(chr_data["cn_both"].max() if join_haps else chr_data[["cn_a", "cn_b"]].max().max())) if len(chr_data) > 0 else 0
    end = max_end / 1000000 if end == 0 else end
    ax.set_xlim(start, end)
    ax.set_ylim(-0.5,  max_cna + .5)
    ax.set_yticks(np.arange(0, max_cna + 1, 1))

    # add a thin line for every y tick
    for y in range(max_cna + 1):
        ax.plot([start, end], [y, y], color="black", linewidth=0.1, alpha=0.5)
    
    for i in range(len(haps)):
        plot_hap(ax, chr_data, chr, haps[i])


def plot_CNs(data, sample, join_haps = False, dpi=150):
    chr_names = chromosome_colors.keys()
    # set 24 subplots vertically stacked
    fig, axs = plt.subplots(len(chr_names), 1)
    fig.set_size_inches(1920 / dpi, 1080 / dpi * 8)
    plt.tight_layout()
    axs[0].set_title(f"Haplotype-specific CNs of sample {sample}")
    data.loc[:,"cn_both"] = data[["cn_a", "cn_b"]].sum(axis=1)
    for i, chr_name in enumerate(chr_names):
        plot_chr(axs[i], data, chr_name, join_haps = join_haps)


def get_data(input_file, sample_name):
    df = pd.read_csv(input_file, sep="\t", index_col=0)
    sample_name = sample_name if sample_name != "" else df.index[0]
    sample_data = df.loc[[sample_name]]
    return sample_data, sample_name

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot haplotype specific Copy Numbers for a sample")
    parser.add_argument("-I", "--input", default="copynumbers.tsv", help="The file with copynumbers.")
    parser.add_argument("-S", "--sample", default="", help="Sample ID")
    parser.add_argument("-O", "--output", default="copy_numbers.png", help="Output file path")
    parser.add_argument("-J", "--joint", action="store_true", help="Plot both haplotypes jointly (default: False))")
    args = parser.parse_args()

    sample_data, sample_name = get_data(args.input, args.sample)
    plot_CNs(sample_data, sample_name, join_haps = args.joint, dpi=150)
    # save using tight layout
    plt.savefig(args.output, dpi=150)