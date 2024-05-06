#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import matplotlib.patches as mpatches
from utils import chr_colors

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
    

def get_data(input_file, sample_name):
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


def plot_karyotype(data, sample_name, dpi=150):
    sample_count = len(data)

    fig, ax = plt.subplots()
    # set figure to full hd, tight layout
    fig.set_size_inches(1920 / dpi, 1080 / dpi)

    # Set the y-axis labels
    ax.set_yticks([i for i in range(sample_count)])
    ylabels = ["contig" + str(i) for i in range(0, sample_count)]
    ylabels.reverse()
    ax.set_yticklabels(ylabels)
    starts = [0.0] * sample_count

    # Loop through the data and draw horizontal arrows
    for i, sample in enumerate(data):
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
            last_left = draw_horizontal_arrow(ax, x, sample_count - i - 1, dx, 0, chr_colors[chr_no], last_left)

    # Set the x-axis limits and labels
    ax.set_xlim(0, np.max(starts))
    ax.set_xlabel("MBase count")
    ax.set_ylim(-1, sample_count)
    ax.set_title(f"Karyotype of sample {sample_name}")

    # plot all chromosome colors as a legend
    for i, (chr_no, color) in enumerate(chr_colors.items()):
        ax.plot(0, 0, color=color, label=chr_no)
        ax.legend(loc="upper left", bbox_to_anchor=(1, 1))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot karyotype for a sample")
    parser.add_argument("-I", "--input", default="karyotypes.tsv", help="The file with the karyotype data.")
    parser.add_argument("-S", "--sample", default="sample_1", help="Sample ID")
    parser.add_argument("-O", "--output", default="karyotype.png", help="Output file path")
    parser.add_argument("--dpi", default=150, help="Output file path", type=int)
    args = parser.parse_args()

    data, sample_name = get_data(args.input, args.sample)
    plot_karyotype(data, sample_name, dpi=args.dpi)    
    plt.savefig(args.output, dpi=args.dpi)
