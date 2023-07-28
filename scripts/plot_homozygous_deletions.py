import pandas as pd
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from os.path import join
from utils import load_dataset, hg19_chr_lengths, chromosome_colors, chromosome_names, hg38_chr_lengths


def genome_length_init():
    # Initialize the cumulative start position variable
    cum_start = 0
    hg19_chr_cum_starts = {}
    # Loop over each chromosome in the chromosome lengths dictionary
    for chrom, length in hg19_chr_lengths.items():
        hg19_chr_cum_starts[chrom] = cum_start
        cum_start += length

    # total length of the genome
    return hg19_chr_cum_starts

def plot_homozygous_deletions(output):
    fig, ax = plt.subplots(1, figsize=(16,9))
    inputs = ["../PCAWG", "../mcmc_simple", "../mcmc_complex", "../mcmc_simple_more_ess", "../mcmc_complex_more_ess"] #"../simple_no_fitness", "../complex_no_fitness"
    labels = ["PCAWG", "Simple, Ess = 0.1", "Complex, Ess = 0.1", "Simple, Ess = 0.5", "Complex, Ess = 0.5"]
    ls = ["-", "--", "--", "-.", "-."]
    for i, input in enumerate(inputs):
        if i == 1 or i == 2:
            df = pd.read_csv(join(input, "copynumbers.tsv"), index_col=0, sep="\t", names=["sample_id", "chrom", "start", "end", "cn_a", "cn_b", "n_snvs"])
        else:
            df = pd.read_csv(join(input, "copynumbers.tsv"), index_col=0, sep="\t")
        # get the unique ids of the samples
        sample_ids = df.index.unique()
        hist = []
        for sample in sample_ids:
            for index, row in df.loc[sample].iterrows():
                # Get the start and end positions
                start = row['start']# + hg19_chr_cum_starts[row['chrom']]
                end = row['end']# + hg19_chr_cum_starts[row['chrom']]
                if (row['chrom'] == "chrX" or row['chrom'] == "chrY"):
                    continue
                if (row['cn_a'] == 0 and row['cn_b'] == 0):
                    hist.append((end-start)/1000000)
        ax.hist(hist, 75 ,range=[0,75], density=True, histtype='step', label=labels[i], linestyle=ls[i])
    ax.set_xlabel("Length of homozygous deletions (MB)")
    ax.set_ylabel("Frequency")
    ax.set_yscale('log')
    ax.set_xbound(lower=0,upper=75)
    ax.legend()
    fig.savefig(join(output, "homozygous_deletions.png"), dpi=150, bbox_inches='tight')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot the average scatter CN plots for a given dataset')
    parser.add_argument('-i', "--input", type=str, help='The folder to input dataset to plot')
    parser.add_argument("--hg19",action="store_true", dest="use_hg19", default=False, help="Use hg19 instead of hg38")
    parser.add_argument('-o', "--output", type=str, help='The folder to output the plots to')
    args = parser.parse_args()

    plot_homozygous_deletions(args.output)