#!/usr/bin/env python3

import pandas as pd
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from os.path import join
from cns.process.binning import bin_by_break_type
from cns.utils.assemblies import hg19_chr_lengths, hg19_chr_starts, human_chr_colors, aut_names, chr_names

# Define a function to convert a row to a list
def row_to_list(row, bins, column):
    # Initialize an empty list
    res = []
    start_index = np.digitize(row['start']+hg19_chr_starts[row['chrom']], bins)
    start = bins[start_index]
    #start = row['start'] + hg19_chr_cum_starts[row['chrom']]
    # Round start to the nearest bin
    #if start % step_size != 0:
    #    start = start - (start % step_size) + step_size
    end_index = np.digitize(row['end']+hg19_chr_starts[row['chrom']]-1, bins)
    end = bins[end_index]
    # The rounding can cause the start to now be greater than the end
    # in which case, it is a single bin
    if start_index == end_index:
        return [(start, row[column])]

    # Loop over each position between start and end, divisible by step
    for index in range(start_index, end_index+1):
        cn = row[column]
        res.append((bins[index], cn))
    # Return the list
    return res


def chromosome_to_bins(chrom, step_size=1_000_000):
    bins = np.arange(step_size, hg19_chr_lengths[chrom], step_size)
    bins += hg19_chr_starts[chrom] - 1
    bins = np.append(bins, hg19_chr_lengths[chrom] + hg19_chr_starts[chrom])
    return bins


def sample_to_homozygous_spots(sample, bins):
    count = {b : 0 for _, b in enumerate(bins)}
    rows = sample.apply(lambda x: row_to_list(x, bins, 'homozygous_del'), axis=1)
    filtered = rows[rows.apply(lambda x: len(x) > 0)]
    for items in filtered.values:
        for pos, cn in items:
            count[pos] += cn
    return pd.Series(count)


def get_homozygous_deletion_locations(df, step_size=1_000_000):
    df['cn'] = df['cn_a'] + df['cn_b']
    # 1 = homozygously deleted region, 0 otherwise
    df['homozygous_del'] = (df['cn'] == 0).astype(int)
    # Ignore sex chromosomes
    df = df[(df['chrom'] != 'chrX') & (df['chrom'] != 'chrY')]
    #df.reset_index(inplace=True)
    #df.set_index(['sample_id', 'chrom'], inplace=True)
    samples = df.index.unique()
    positions = {}
    bins = get_chromosome_bins(step_size, includeSexChromosomes=False)
    # Loop over the samples
    for sample in samples:
        # For each chromosome, count the number of homozygous deletions
        positions[sample] = sample_to_homozygous_spots(df.loc[sample,:], bins)
        #positions[sample] = sample_to_homozygous_spots(df.loc[sample, :], bins, step_size)
    df_dels = pd.DataFrame(positions)
    return df_dels


def homozygous_length_distribution(df, bin_size=1_000_000):
    df["cn"] = df["cn_a"] + df["cn_b"]
    sample_ids = df["sample_id"].unique()
    hist = []
    for id in sample_ids:
        sample = df[df["sample_id"] == id]
        for _, row in sample.iterrows():
            # Skip segments that are not homozygously deleted or are sex chromosomes
            if row["cn"] != 0 or row['chrom'] in ["chrX", "chrY"]:
                continue
            if row["cn_a"] == 0 and row["cn_b"] == 0:
                length = (row["end"] - row["start"])/bin_size
                hist.append(length)
    return hist


def get_chromosome_bins(step_size=1_000_000, includeSexChromosomes=False):
    bins = np.array([])
    chroms = aut_names if not includeSexChromosomes else chr_names
    for chrom in chroms:
        chr_bins = chromosome_to_bins(chrom, step_size=step_size)
        bins = np.append(bins, chr_bins)
    return bins


def sample_to_SNPs(sample, bins, column):
    count = {b : 0 for _, b in enumerate(bins)}
    rows = sample.apply(lambda x: row_to_list(x, bins, column), axis=1)
    filtered = rows[rows.apply(lambda x: len(x) > 0)]
    for items in filtered.values:
        for pos, cn in items:
            count[pos] += cn
    return pd.Series(count)
    #return pd.Series(np.concatenate(filtered.values))


def samples_to_SNPs(cns, column='cn', step_size=1_000_000, includeSexChromosomes=False):
    # list of unique indices in cns
    samples = cns.index.unique()
    positions = {}
    bins = get_chromosome_bins(step_size, includeSexChromosomes)
    for sample in samples:
        positions[sample] = sample_to_SNPs(cns.loc[sample, :], bins, column)
    df_CNs = pd.DataFrame(positions)
    return df_CNs



def plot_scatter_CNs(data, use_hg19=True, output="../out"):
    cns = pd.read_csv(join(data, "copynumbers.tsv"), index_col=0, sep="\t")
    cns['cn'] = cns['cn_a'] + cns['cn_b']
    # list of unique indices in cns
    samples = cns.index.unique()
    positions = {}
    step_size = 10_000_000
    for sample in samples:
        positions[sample] = bin_by_break_type(cns.loc[sample, :], step_size) 
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
        color = human_chr_colors[chrom]
        
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
    ax.set_xticks([hg19_chr_starts[chrom] + (length / 2) for chrom, length in hg19_chr_lengths.items()])
    # set tick labels to be chromosome names
    ax.set_xticklabels(list(hg19_chr_starts.keys()))

    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Copy number")
    ax.set_title("Copy number variation")

    # plot a line at y=2
    ax.plot(df_CNs.index, [2] * len(df_CNs.index), color="black", linewidth=1, alpha=0.7)

    fig.savefig(f"{output}/scatter_CNs.png", dpi=150, bbox_inches="tight")


def plot_scatter_CNs():
    datasets = {"PCAWG":"../out/PCAWG_filtered_95_pc/"}#, "Stick-Break":"../out/simulated_stick_break/", "Limited Stick-Break": "../out/simulated_stick_break_limited/"}
    ls = ["-", "-.", "--", ":"]

    x_pos = 0
    for chrom, length in hg19_chr_lengths.items():
        # Update the x position for the next chromosome
        x_pos += length
     
    fig, ax = plt.subplots(1, figsize=(32, 9))
    for i, key in enumerate(datasets.keys()):
        data = datasets[key]
        cns = pd.read_csv(join(data,"copynumbers.tsv"), index_col = 0, sep="\t")
        cns['cn'] = cns['cn_a'] + cns['cn_b']
        step_size = 1_000_000
        df_CNs = samples_to_SNPs(cns, 'cn', step_size, includeSexChromosomes=True)
        # Set the y-ax[0]is limits
        ax.set_ylim(0.4, 3.5)

        ax.plot(df_CNs.index, df_CNs.mean(axis=1), label=key, linewidth=1.5, ls = ls[i])

        # Set the x-ax[0]is limits to match the total genome length
        ax.set_xlim(0, x_pos)

        if i == 0:
            # plot a line at y=2
            ax.plot(df_CNs.index, [2] * len(df_CNs.index), color="black", linewidth=1, alpha=0.7, ls=":")
    
    # Set the x-ax[0]is ticks to be at the middle of each chromosome
    ax.set_xticks([hg19_chr_starts[chrom] + (length / 2) for chrom, length in hg19_chr_lengths.items()])
    # set tick labels to be chromosome names
    ax.set_xticklabels(list(hg19_chr_starts.keys()))
    ax.legend()

    ax.set_xlabel("Genomic position (1MB bins)")
    ax.set_ylabel("Mean copy-number")
    ax.set_title("Mean Copy-Number Profile of PCAWG Dataset")


    fig.savefig(f"./scatter_CNs_pcawg.png", dpi=300, bbox_inches="tight")
    




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot the average scatter CN plots for a given dataset')
    parser.add_argument('-I', "--input", type=str, default=".", help='The folder to input dataset to plot')
    parser.add_argument("--hg19",action="store_true", dest="use_hg19", default=False, help="Use hg19 instead of hg38")
    parser.add_argument('-O', "--output", type=str, default=".", help='The folder to output the plots to')
    args = parser.parse_args()

    plot_scatter_CNs()
    #plot_scatter_CNs(args.input, args.use_hg19, args.output)




