import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import matplotlib.patches as mpatches

chromosome_colors = {
    'chr1': 'red',
    'chr2': 'mediumblue',
    'chr3': 'forestgreen',
    'chr4': 'darkorange',
    'chr5': 'dodgerblue',
    'chr6': 'olivedrab',
    'chr7': 'purple',
    'chr8': 'gray',
    'chr9': 'gold',
    'chr10': 'salmon',
    'chr11': 'turquoise',
    'chr12': 'darkviolet',
    'chr13': 'green',
    'chr14': 'indianred',
    'chr15': 'steelblue',
    'chr16': 'sienna',
    'chr17': 'royalblue',
    'chr18': 'darkgoldenrod',
    'chr19': 'mediumvioletred',
    'chr20': 'teal',
    'chr21': 'peru',
    'chr22': 'navy',
    'chrX': 'chocolate',
    'chrY': 'darkslateblue'
}
   

def parse_region_string(s):
    # Define the regular expression pattern to match the string
    pattern = r'(chr[\d,Y,X]+)_H([1, 2])([<>])\[(\d+):(\d+)\)'
    
    # Use re.search() to find the pattern in the string and extract the groups
    match = re.search(pattern, s)
    if match:
        # Extract the groups from the match object
        chromosome = match.group(1)
        haplotype = match.group(2)
        direction = True if match.group(3) == '>' else False
        start = int(match.group(4))
        end = int(match.group(5))
        # Calculate the length in millions and round to 2 decimal places
        length = round((end - start) / 1000000, 2)

        # Return the extracted values as a dictionary
        return {
            'chr': chromosome,
            'haplotype': haplotype,
            'direction': direction,
            'length': length
        }
    else:
        return None
    

def parse_karyotype(kar):
    contigs = kar.split(';')
    return [[parse_region_string(seg) for seg in contig[1:-1].split('~')] for contig in contigs]


def draw_horizontal_arrow(ax, x, y, dx, dy, color, last_left):
    width = 0.1
    head_len = min(5, abs(dx / 2))
    arrow = mpatches.FancyArrow(x, y - width/2, dx, dy, fc=color, ec=color, width=width, head_width=0.5, head_length=head_len, length_includes_head=True, overhang=.5, alpha=0.5)
    if last_left and dx > 0:
        end_width = min(1, abs(dx / 2))
        end_height = 0.5
        ending = mpatches.Rectangle((x - end_width / 2, y - end_height / 2), end_width, end_height, fc=color, ec=color, alpha=0.3)        
        ax.add_patch(ending)    
    ax.add_patch(arrow)
    return dx < 0
    

def get_data(input_file, sample_name):
    # check if input file exists
    try:
        open(input_file, 'r')
    except IOError:
        print(f"File {input_file} not found")
        return
    clones = pd.read_csv(input_file, sep='\t')

    # set clones["ID"] to string
    clones["ID"] = clones["ID"].astype(str)    
    clones.set_index("ID", inplace=True)    

    # clones row where ID is sample_name
    sample = clones.loc[sample_name]
    if (sample.empty):
        print(f"Sample {sample_name} not found")
        return
    
    return parse_karyotype(sample["Karyotype"])


def plot_karyotype(data, sample_name, output_file):
    dpi = 150
    sample_count = len(data)

    fig, ax = plt.subplots()
    # set figure to full hd, tight layout
    fig.set_size_inches(1920 / dpi, 1080 / dpi)
    fig.tight_layout(pad=1.25, rect=[0.025, 0, .925, 1])    

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
            length = item["length"]
            direction = item["direction"]
            x = starts[i] + (0 if direction else length)
            dx = length if direction else -length
            starts[i] += abs(dx)
            last_left = draw_horizontal_arrow(ax, x, sample_count - i - 1, dx, 0, chromosome_colors[chr_no], last_left)

    # Set the x-axis limits and labels
    ax.set_xlim(0, np.max(starts))
    ax.set_xlabel("MBase count")
    ax.set_ylim(-1, sample_count)
    ax.set_title(f"Karyotype of sample {sample_name}")

    # plot all chromosome colors as a legend
    for i, (chr_no, color) in enumerate(chromosome_colors.items()):
        ax.plot(0, 0, color=color, label=chr_no)
        ax.legend(loc="upper left", bbox_to_anchor=(1, 1))

    # Show the plot
    plt.savefig(output_file, dpi=dpi)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot karyotype for a sample")
    parser.add_argument("--input", default="./out/clones.tsv", help="Input file path")
    parser.add_argument("--sample", default="1", help="Sample ID")
    parser.add_argument("--output", default="./out/karyotype.png", help="Output file path")
    args = parser.parse_args()

    data = get_data(args.input, args.sample)
    plot_karyotype(data, args.sample, args.output)