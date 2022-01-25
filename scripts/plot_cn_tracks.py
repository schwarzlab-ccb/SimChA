import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo

from dot_to_newick import convert_dot_to_newick as dot_to_newick
from plotting_functions import (_get_y_positions, plot_chrom_boundaries,
                                plot_single_copynumber, plot_tree)

TREE_WIDTH_SCALE = 1
TRACK_WIDTH_SCALE = 1
HEIGHT_SCALE = 1

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_folder", type=str, default=None, required=False)
    parser.add_argument("-o", "--output_folder", type=str, default=None, required=False)
    parser.add_argument("--fraction", action="store_true", dest="fraction")
    args = parser.parse_args()

    if args.input_folder is None:
        input_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), '../SimChA/out'))

    if args.output_folder is None:
        output_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), '../SimChA/out'))

    # load copynumbers
    copynumbers = pd.read_csv(os.path.join(input_folder, 'copynumbers.out'), sep='\t', index_col=0)
    samples = copynumbers.index.unique()

    # load tree
    if not 'parent_graph.new' in os.listdir(input_folder):
        if 'parent_graph.dot' in os.listdir(input_folder):
            dot_to_newick(os.path.join(input_folder, 'parent_graph.dot'),
                          os.path.join(input_folder, 'parent_graph.new'))
        else:
            raise BaseException(f"Did not find either 'parent_graph.new' or 'parent_graph.dot'"
                                 "in input folder {input_folder}")

    tree = Phylo.read(os.path.join(input_folder, 'parent_graph.new'), 'newick')
    tree.ladderize()  # Flip branches so deeper clades are displayed at top

    # Display the population size as fractions
    if args.fraction:
        total_pop = 0
        for clade in tree.find_clades():
            if clade.name is not None:
                total_pop += int(clade.name.split('-')[1])
        for clade in tree.find_clades():
            if clade.name is not None:
                clade.name = clade.name.split('-')[0] + '-' + \
                    str(np.round(float(clade.name.split('-')[1])/total_pop, 3))

    normal_sample = [x for x in tree.find_clades(
    ) if x.name is not None and x.name.split('-')[0] == '0']
    if len(normal_sample) == 1:
        normal_sample = normal_sample[0].name
    else:
        raise TypeError(f"multiple normal samples: '0' for normal name {normal_sample}")


    # Figure dimensions
    nsegs = len(copynumbers.loc[samples[0]])
    track_width = nsegs * 0.2 * TRACK_WIDTH_SCALE
    tree_width = 2.5 * TREE_WIDTH_SCALE  # in figsize
    nrows = len(samples)
    plotheight = 4 * 0.2 * nrows * HEIGHT_SCALE
    plotwidth = tree_width + track_width
    tree_width_ratio = tree_width / plotwidth
    fig = plt.figure(figsize=(plotwidth, plotheight), constrained_layout=True)

    # plot tree
    gs = fig.add_gridspec(nrows, 2, width_ratios=[tree_width_ratio, 1-tree_width_ratio])
    tree_ax = fig.add_subplot(gs[0:len(samples), 0])
    cn_axes = [fig.add_subplot(gs[i]) for i in range(1, (2*(nrows))+1, 2)]
    y_posns = _get_y_positions(tree, normal_sample, adjust=True)
    y_posns = {clade.name.split('-')[0]: y_pos for clade, y_pos in y_posns.items()
        if clade.name is not None and clade.name != 'root'}
    plot_tree(tree,
              ax=tree_ax,
              normal_name=normal_sample)

    # plot copy-number tracks
    ymax = copynumbers[['cn_a', 'cn_b']].max().max()
    for sample, group in copynumbers.groupby('sample_id'):
        if sample not in samples:
            continue
        index_to_plot = y_posns[str(sample)] - 1
        plot_single_copynumber(group,
                               ax=cn_axes[index_to_plot],
                               print_chrom=False,
                            #    print_chrom=index_to_plot==0,
                               ymax=ymax)
        cn_axes[index_to_plot].set_ylabel('')
    for ax in cn_axes[:-1]:
        ax.set_xticks([])
    plot_chrom_boundaries(cn_axes[0], print_chrom=True, outside=True)

    # save
    fig.savefig(os.path.join(input_folder, 'copy_numbers.pdf'))
