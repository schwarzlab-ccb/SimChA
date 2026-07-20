#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
from os.path import join, abspath, dirname
import matplotlib.patches as mpatches

human_chr_colors = {
    "chr1": "#FF0000",
    "chr2": "#00FF00",
    "chr3": "#0000FF",
    "chr4": "#FFFF00",
    "chr5": "#FF00FF",
    "chr6": "#00FFFF",
    "chr7": "#800000",
    "chr8": "#008000",
    "chr9": "#000080",
    "chr10": "#808000",
    "chr11": "#800080",
    "chr12": "#008080",
    "chr13": "#C0C0C0",
    "chr14": "#FFA500",
    "chr15": "#A52A2A",
    "chr16": "#800080",
    "chr17": "#008000",
    "chr18": "#0000FF",
    "chr19": "#FF1493",
    "chr20": "#1E90FF",
    "chr21": "#32CD32",
    "chr22": "#FFD700",
    "chrX": "#DA70D6",
    "chrY": "#A9A9A9"
}

CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def get_test_data():
    # Four datasets with length (float) and direction (bool)
    return [[
        {"chr": "chr3", "length": 3.2, "direction": True},
        {"chr": "chr1", "length": 4.6, "direction": False},
        {"chr": "chr2", "length": 2.8, "direction": True},
        {"chr": "chr3", "length": 5.1, "direction": False},
    ],[
        {"chr": "chr8", "length": 6.2, "direction": True},
        {"chr": "chr4", "length": 5.1, "direction": False},
    ]]


def parse_region_string(s):
    # Define the regular expression pattern to match the string
    pattern = r'H([1, 2])([<>])(chr[\d,Y,X]+)\[(\d+):(\d+)\)'
    
    # Use re.search() to find the pattern in the string and extract the groups
    match = re.search(pattern, s)
    if match:
        # Extract the groups from the match object
        haplotype = match.group(1)
        direction = True if match.group(2) == '>' else False
        chromosome = match.group(3)
        start = int(match.group(4))
        end = int(match.group(5))

        # Return the extracted values as a dictionary
        return {
            'chr': chromosome,
            'haplotype': haplotype,
            'direction': direction,
            'start': start,
            'end': end
        }
    else:
        return None
    

def parse_karyotype(kar):
    contigs = kar.split(';')
    return [[parse_region_string(seg) for seg in contig[1:-1].split('~')] for contig in contigs]


def draw_horizontal_arrow(ax, x, y, dx, dy, color, last_left):
    width = 0.1
    head_len = min(2.5, abs(dx / 2))
    arrow = mpatches.FancyArrow(x, y, dx, dy, fc=color, ec=color, width=width, head_width=0.5, head_length=head_len, length_includes_head=True, overhang=.5, alpha=0.5)
    if last_left and dx > 0:
        end_width = min(.5, abs(dx / 2))
        end_height = 0.5
        ending = mpatches.Rectangle((x - end_width / 2, y - end_height / 2), end_width, end_height, fc=color, ec=color, alpha=0.5, linewidth=0)        
        ax.add_patch(ending)    
    ax.add_patch(arrow)
    return dx < 0
    

def get_data(input_file, sample_name = ""):
    # check if input file exists
    try:
        open(input_file, 'r')
    except IOError:
        print(f"File {input_file} not found")
        raise SystemExit
    clones = pd.read_csv(input_file, sep='\t')

    clones["sample_id"] = clones["sample_id"].astype(str)    
    clones.set_index("sample_id", inplace=True)    

    sample_name = sample_name if sample_name != "" else clones.index[0]

    # clones row where ID is sample_name
    sample = clones.loc[sample_name]
    if (sample.empty):
        print(f"Sample {sample_name} not found")
        raise SystemExit
    
    return parse_karyotype(sample["karyotype"]), sample_name


def plot_karyotype(contigs, dpi=100):
    contig_count = len(contigs)

    fig, ax = plt.subplots()
    # set figure to full hd, tight layout
    fig.set_size_inches(1920 / dpi, 32 * len(contigs) / dpi)

    # Set the y-axis labels
    ax.set_yticks([i for i in range(contig_count)])
    ylabels = ["contig" + str(i) for i in range(0, contig_count)]
    ylabels.reverse()
    ax.set_yticklabels(ylabels)
    starts = [0.0] * contig_count

    # Loop through the data and draw horizontal arrows
    for i, sample in enumerate(contigs):
        last_left = False
        for j, item in enumerate(sample):
            chr_no = item["chr"]
            start = item["start"]
            end = item["end"]
            # Calculate the length in millions and round to 2 decimal places
            length =  round((end - start) / 1000000, 2)
            direction = item["direction"]
            x = starts[i] + (0 if direction else length)
            dx = length if direction else -length
            starts[i] += abs(dx)
            last_left = draw_horizontal_arrow(ax, x, contig_count - i - 1, dx, 0, human_chr_colors[chr_no], last_left)

    # Set the x-axis limits and labels
    ax.set_xlim(0, np.max(starts))
    ax.set_xlabel("MBase count")
    ax.set_ylim(-1, contig_count)
    # plot all chromosome colors as a legend
    for i, (chr_no, color) in enumerate(human_chr_colors.items()):
        ax.plot(0, 0, color=color, label=chr_no)
        ax.legend(loc="upper left", bbox_to_anchor=(1, 1))

    return fig, ax


def get_cn_data(input_file, sample_name=""):
    cn_df = pd.read_csv(input_file, sep='\t')
    cn_df["sample_id"] = cn_df["sample_id"].astype(str)
    sample_name = sample_name if sample_name != "" else cn_df["sample_id"].iloc[0]
    sample_df = cn_df[cn_df["sample_id"] == sample_name].copy()
    if sample_df.empty:
        raise ValueError(f"Sample {sample_name} not found in {input_file}")
    return sample_df, sample_name


def plot_cn_tracks(sample_df, sample_name="", dpi=100):
    offset = 0.3
    present_chroms = sample_df["chrom"].unique()
    chroms = [c for c in CHROMOSOMES if c in present_chroms]
    nrows = len(chroms)

    global_max_pos = sample_df["end"].max()

    # Compute per-chromosome CN range for height ratios
    height_ratios = []
    for chrom in chroms:
        chrom_df = sample_df[sample_df["chrom"] == chrom]
        if chrom_df.empty:
            height_ratios.append(2)
        else:
            max_cn = max(chrom_df["cn_a"].max(), chrom_df["cn_b"].max(), 2)
            height_ratios.append(max_cn + 1.0)

    total_height = sum(height_ratios) * 0.2
    fig, axes = plt.subplots(nrows, 1, figsize=(7, max(total_height, 6)), dpi=dpi,
                             sharex=True, gridspec_kw={"height_ratios": height_ratios})
    if nrows == 1:
        axes = [axes]

    for idx, chrom in enumerate(chroms):
        ax = axes[idx]
        color = human_chr_colors[chrom]
        chrom_df = sample_df[sample_df["chrom"] == chrom]

        if chrom_df.empty:
            ax.set_ylabel(chrom, fontsize=8, rotation=0, labelpad=40, va="center")
            ax.set_yticks([])
            continue

        for _, row in chrom_df.iterrows():
            start = row["start"]
            width = row["end"] - row["start"]
            cn_a = row["cn_a"]
            cn_b = row["cn_b"]

            # cn_a: rectangle on [cn_a, cn_a + offset]
            rect_a = mpatches.Rectangle(
                (start, cn_a), width, offset,
                fc=color, ec="black", alpha=0.7, linewidth=0.4)
            ax.add_patch(rect_a)

            # cn_b: rectangle on [cn_b - offset, cn_b]
            rect_b = mpatches.Rectangle(
                (start, cn_b - offset), width, offset,
                fc=color, ec="black", alpha=0.4, linewidth=0.4)
            ax.add_patch(rect_b)

        max_pos = chrom_df["end"].max()
        max_cn = max(chrom_df["cn_a"].max(), chrom_df["cn_b"].max(), 2)
        ax.set_xlim(0, global_max_pos)
        ax.set_ylim(-0.5, max_cn + 0.5)
        ax.set_yticks(range(0, int(max_cn) + 1))
        ax.axhline(0, color="black", linewidth=0.5, linestyle="--")
        ax.set_ylabel(chrom, fontsize=8, rotation=0, labelpad=40, va="center")
        ax.tick_params(labelsize=7)
        # Hide x ticks on all but the bottom subplot
        if idx < nrows - 1:
            ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)

    # Bottom subplot: show tick marks and labels at 20 Mb intervals
    import matplotlib.ticker as mticker
    axes[-1].tick_params(axis="x", which="both", bottom=True, labelbottom=True)
    axes[-1].xaxis.set_major_locator(mticker.MultipleLocator(10e6))
    axes[-1].xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x / 1e6)}"))
    axes[-1].set_xlabel("Mb", fontsize=9)

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor="gray", alpha=0.7, edgecolor="black", label="cn_a (above)"),
        mpatches.Patch(facecolor="gray", alpha=0.4, edgecolor="black", label="cn_b (below)"),
    ]
    fig.legend(handles=legend_elements, loc="lower right", fontsize=8, ncols=2)

    title = f"CN tracks – {sample_name}" if sample_name else "CN tracks"
    fig.suptitle(title, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.subplots_adjust(hspace=0.1)
    return fig, axes


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot karyotype or CN tracks for a sample")
    parser.add_argument("-I", "--input", default="karyotypes.tsv", help="The file with the karyotype data.")
    parser.add_argument("-S", "--sample", default="Sample_1", help="Sample ID")
    parser.add_argument("-O", "--output", default="", help="Output file path")
    parser.add_argument("--dpi", default=100, help="DPI for output", type=int)
    parser.add_argument("--cn", action="store_true", help="Plot CN tracks instead of karyotype")
    args = parser.parse_args()

    if args.cn:
        cn_file = join(abspath(dirname(args.input)), "copynumbers.tsv")
        sample_df, sample_name = get_cn_data(cn_file, args.sample)
        fig, axes = plot_cn_tracks(sample_df, sample_name, dpi=args.dpi)
    else:
        data, sample_name = get_data(args.input, args.sample)
        fig, ax = plot_karyotype(data, dpi=args.dpi)
        ax.set_title(f"Karyotype of sample {sample_name}")

    out_path = args.output if args.output != "" else join(abspath(dirname(args.input)), f"{sample_name}_karyotype.png")
    print(f"Saving plot to {out_path}")
    plt.savefig(out_path, dpi=args.dpi)
