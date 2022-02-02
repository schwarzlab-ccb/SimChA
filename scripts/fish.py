#!/usr/bin/env python

import argparse
import random
import sys
from datetime import timedelta
from timeit import default_timer as timer

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import cm
from PIL import Image


def create_ordering(tree, cur_clone):
    res = []
    if cur_clone in tree.keys():
        res += [cur_clone]
        for child in tree[cur_clone]:
            res += create_ordering(tree, child)
        res += [cur_clone]
    else:
        return 2 * [cur_clone]
    return res


# A dict-based tree where for each parent there is a list of children
def build_tree(parent_df, ids, root):
    tree = {cid: list(children) for cid in ids if len(
        children := parent_df.loc[parent_df["ParentId"] == cid, "ChildId"]) > 0}
    if -1 not in tree:
        tree[-1] = [root]
    return tree


if __name__ == '__main__':
    randint = random.randint(0, sys.maxsize)
    parser = argparse.ArgumentParser(
        description='Create a Fish (Muller) plot for the given evolutionary tree.')
    parser.add_argument("populations", type=str, help="A CSV file with the header \"Id,Step,Pop\".")
    parser.add_argument("parent_tree", type=str,
                        help="A CSV file with the header \"ParentId,ChildId\".")
    parser.add_argument("output", type=str,
                        help="Output image filepath. The format must support alpha channels.")
    parser.add_argument("-F", "--first", dest="first_step", type=int,
                        help="The step to start plotting from.")
    parser.add_argument("-L", "--last", dest="last_step", type=int,
                        help="The step to end the plotting at.")
    parser.add_argument("-S", "--seed", dest="seed", type=int,
                        help="Seed for colors", default=randint)
    parser.add_argument('-a', '--absolute', dest="absolute", action="store_true", default=False,
                        help='Plot the populations in absolute numbers rather than normalized')

    args = parser.parse_args()

    random.seed(args.seed)

    # parent tree
    parent_df = pd.read_csv(args.parent_tree)
    parents = parent_df["ParentId"].unique()
    children = parent_df["ChildId"].unique()
    root_list = np.setdiff1d(parents, children)
    if len(root_list) != 1:
        raise Exception(
            "Failed to determine root. There must be exactly one node with no parent in the parent tree.")
    root_id = root_list[0]
    ids = np.concatenate([[root_id], children])
    tree = build_tree(parent_df, ids, root_id)

    ordering = create_ordering(tree, 0)

    # populations dataframe
    populations_df = pd.read_csv(args.populations)

    min_step = populations_df["Step"].min()
    first_step = max(args.last_step, min_step) if args.first_step else min_step

    max_step = populations_df["Step"].max()
    last_step = min(args.last_step, max_step) if args.last_step else max_step

    populations_df = populations_df.loc[(populations_df["Step"] >= first_step)
                                        & (populations_df["Step"] <= last_step)]

    populations_df = populations_df.set_index(['Id', 'Step'])[['Pop']].unstack('Step').fillna(0)
    populations_df = populations_df.loc[ordering] / 2

    if not args.absolute:
        populations_df = populations_df / populations_df.sum(axis=0)

    # Create colors
    colors = np.array(cm.rainbow(np.linspace(0, 1, len(ids))).tolist())
    np.random.shuffle(colors)
    colors = pd.DataFrame(colors, index=ids)
    colors.loc[root_id] = [0.5, 0.5, 0.5, 1]
    colors = colors.loc[ordering].values

    # Plot
    plt.figure(figsize=(20, 10))
    plt.stackplot(np.arange(populations_df.shape[1]), populations_df, colors=colors)

    plt.xlim(first_step, last_step)
    if not args.absolute:
        plt.ylim(0, 1)

    plt.savefig(args.output)
