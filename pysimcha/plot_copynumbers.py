#!/usr/bin/env python3

import argparse
import os

from Bio import Phylo

from dot_to_newick import convert_dot_to_newick as dot_to_newick
from plotting_functions import plot_cn_bars, plot_cn_heatmap, load_data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True)
    parser.add_argument('-t', '--tree', type=str, default=None, required=False)
    parser.add_argument("-o", "--output-folder", dest='output_folder', type=str, default=None, required=False)
    parser.add_argument("--fraction", action="store_true", dest="fraction")

    parser.add_argument('--type', type=str, default='heatmap', required=False,
                        choices=["bars", "heatmap", "both"],
                        help="""Type of copy-number plot to save. 'bars' is recommended for <10 samples, 
                         heatmap for more samples, 'both' will plot both. (default: auto).""")

    args = parser.parse_args()

    if args.output_folder is not None:
        output_folder = args.output_folder
    else:
        output_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), '../out'))

    # load copynumbers
    copynumbers = load_data(args.input)
    assert copynumbers.unstack('sample_id').isna().sum().sum() == 0, 'no consistent segmentation'
    cmax = copynumbers[['cn_a', 'cn_b']].max().max()

    # load tree
    if args.tree is None:
        tree = None
    elif not os.path.exists(args.tree):
        raise FileNotFoundError(f"Did not find tree file {args.tree}")
    else:
        if args.tree.endswith('.dot'):
            dot_to_newick(args.tree, args.tree.replace('.dot', '.new'))
        elif not args.tree.endswith('.new'):
            raise BaseException("Tree must be either .dot or .new format.")

        tree = Phylo.read(args.tree.replace('dot', 'new'), 'newick')
        tree.ladderize()  # Flip branches so deeper clades are displayed at top

    if args.type == 'both':
        plot_type = ['heatmap', 'bars']
    else:
        plot_type = [args.type]

    if 'bars' in plot_type:
        fig = plot_cn_bars(copynumbers=copynumbers, tree=tree,
                           fraction=args.fraction, cmax=8, tree_width_ratio=1)
        fig.savefig(os.path.join(output_folder, 'copy_numbers_bars.png'))

    if 'heatmap' in plot_type:
        fig = plot_cn_heatmap(copynumbers=copynumbers,
                              tree=tree, figsize=(20, 15),
                              cmax=cmax,
                              hide_internal_nodes=True,
                              ignore_segment_lengths=False
                              )
        fig.savefig(os.path.join(output_folder, 'copy_numbers_heatmap.png'))
