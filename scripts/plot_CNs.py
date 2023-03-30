import argparse
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
import plot_karyotype as pk


def get_chr_data(data, chr_name, haplotype = "both"):
    chr_data = []
    for sample in data:
        for item in sample:
            if item["chr"] == chr_name and (item["haplotype"] == haplotype or haplotype == "both"):
                chr_data.append(item)
    return chr_data


def get_segmentation(chr_data):
    breakpoints = [seg["start"] for seg in chr_data] + [seg["end"] for seg in chr_data]
    breakpoints = list(set(breakpoints))
    breakpoints.sort()
    return breakpoints


# calculate the copy numbers
def get_copy_numbers(chr_data, breakpoints):
    copy_numbers = []
    for i in range(len(breakpoints) - 1):
        start = breakpoints[i]
        end = breakpoints[i + 1]
        copy_number = 0
        for seg in chr_data:
            if seg["start"] <= start and seg["end"] >= end:
                copy_number += 1
        copy_numbers.append(copy_number)
    return copy_numbers


# a plot with rectangles between start-end on a line given by the copy number
def plot_chr(ax, breakpoints, cnas, chr, hap, start = 0, end = 0):    
    ax.set_ylabel(f"CN {chr}")
    if (len(breakpoints) == 0):
        return
    end = breakpoints[-1] / 1000000 if end == 0 else end
    ax.set_xlim(start, end)
    ax.set_ylim(-0.5, max(cnas) + .5)
    ax.set_yticks(np.arange(0, max(cnas) + 1, 1))
    for i in range(len(breakpoints) - 1):
        start = breakpoints[i] / 1000000
        end = breakpoints[i + 1] / 1000000
        copy_number = cnas[i]
        if hap == "1":
            height = .25
            y_pos = copy_number
        elif hap == "2":
            height = .25
            y_pos = copy_number - .25
        else:
            height = .5
            y_pos = copy_number - .25
            
        ax.add_patch(mpatches.Rectangle((start, y_pos), end - start, height, fc=pk.chromosome_colors[chr], edgecolor="black"))
        # add a thin line for every y tick
        for y in range(max(cnas) + 1):
            ax.plot([start, end], [y, y], color="black", linewidth=0.1, alpha=0.5)


def plot_CNs(data, sample, dpi=200):
    chr_names = pk.chromosome_colors.keys()
    # set 24 subplots vertically stacked
    fig, axs = plt.subplots(len(chr_names), 1)
    fig.set_size_inches(1920 / dpi, 1080 / dpi * 8)
    axs[0].set_title(f"Haplotype-specific CNs of sample {sample}")
    for i, chr_name in enumerate(chr_names):
        # set figure to full hd, tight layout  
        for hap in ["1", "2"]:
            chr_data = get_chr_data(data, chr_name, hap)
            breakpoints = get_segmentation(chr_data)
            cnas = get_copy_numbers(chr_data, breakpoints)
            plot_chr(axs[i], breakpoints, cnas, chr_name, hap)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot karyotype for a sample")
    parser.add_argument("--input", default="./out/clones.tsv", help="Input file path")
    parser.add_argument("--sample", default="1", help="Sample ID")
    parser.add_argument("--output", default="./out/karyotype.png", help="Output file path")
    args = parser.parse_args()

    data = pk.get_data(args.input, args.sample)
    plot_CNs(data, args.sample, dpi=200)
    # save using tight layout
    plt.tight_layout()
    plt.savefig(args.output, dpi=200)