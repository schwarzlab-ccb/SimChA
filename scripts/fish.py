#!/usr/bin/env python

import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import cm
from scipy.ndimage import gaussian_filter


def stackplot(x, *args, ax=None, colors=None, labels=(), **kwargs):
    '''Taken largely from the matplotlib implemenation except that the keywords `edgecolor` and 
    `where` are provided to `fill_between`'''

    y = np.row_stack(args)

    labels = iter(labels)

    if ax is None:
        ax = plt.gca()

    if colors is not None:
        ax.set_prop_cycle(color=colors)

    # Assume data passed has not been 'stacked', so stack it here.
    # We'll need a float buffer for the upcoming calculations.
    stack = np.cumsum(y, axis=0, dtype=np.promote_types(y.dtype, np.float32))

    # Color between x = 0 and the first array.
    color = ax._get_lines.get_next_color()
    coll = ax.fill_between(x, 0, stack[0, :], where=(stack[0, :] != 0), interpolate=True,
                           facecolor=color, edgecolor=color, label=next(labels, None),
                           **kwargs)
    coll.sticky_edges.y[:] = [0]
    r = [coll]

    # Color between array i-1 and array i
    for i in range(len(y) - 1):
        color = ax._get_lines.get_next_color()
        r.append(ax.fill_between(x, stack[i, :], stack[i+1, :], 
                                 where=(stack[i, :] != stack[i+1, :]), facecolor=color,
                                 edgecolor=color, label=next(labels, None), interpolate=True,
                                 **kwargs))

    return r


def create_ordering(cur_tree, cur_clone):
    res = []
    if cur_clone in cur_tree.keys():
        res += [cur_clone]
        for child in cur_tree[cur_clone]:
            res += create_ordering(cur_tree, child)
        res += [cur_clone]
    else:
        return [cur_clone]
    return res


def build_tree(parent_tree):
    '''A dict-based tree where for each parent there is a list of children'''
    parent_df = pd.read_csv(parent_tree)
    parents = parent_df["ParentId"].unique()
    children = parent_df["ChildId"].unique()
    root_list = np.setdiff1d(parents, children)
    if len(root_list) != 1:
        raise Exception(
            "Failed to determine root. There must be exactly one node with no parent in the parent tree.")
    root_id = root_list[0]
    ids = np.concatenate([[root_id], children])

    tree = {cid: list(children) for cid in ids if len(
        children := parent_df.loc[parent_df["ParentId"] == cid, "ChildId"]) > 0}
    if -1 not in tree:
        tree[-1] = [root_id]

    return tree, ids, root_id


if __name__ == '__main__':
    randint = np.random.randint(0, 2**32 - 1)
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
    parser.add_argument("--smooth", type=int, default=None,
                        help="STD for Gaussian smoothing")
    parser.add_argument("--infer-empty", dest='infer_empty', action="store_true", default=False,
                        help="Whether to infer empty entries")
    parser.add_argument("-S", "--seed", dest="seed", type=int,
                        help="Seed for colors", default=randint)
    parser.add_argument('-a', '--absolute', dest="absolute", action="store_true", default=False,
                        help='Plot the populations in absolute numbers rather than normalized')

    args = parser.parse_args()

    np.random.seed(args.seed)

    # population dataframe
    populations_df = pd.read_csv(args.populations)

    min_step = populations_df["Step"].min()
    first_step = max(args.first_step, min_step) if args.first_step else min_step

    max_step = populations_df["Step"].max()
    last_step = min(args.last_step, max_step) if args.last_step else max_step

    populations_df = populations_df.loc[(populations_df["Step"] >= first_step) &
                                        (populations_df["Step"] <= last_step)]
    samples = populations_df['Id'].unique()

    # parent tree
    tree, ids, root_id = build_tree(args.parent_tree)

    ordering = create_ordering(tree, 0)
    ordering = [x for x in ordering if x in samples]

    # populations dataframe processing
    populations_df = populations_df.set_index(['Id', 'Step'])[['Pop']].unstack('Step')

    if args.infer_empty:
        populations_df[('Pop', 0)] = populations_df[('Pop', 0)].fillna(0)
        populations_df[('Pop', populations_df.shape[1] - 1)
                    ] = populations_df[('Pop', populations_df.shape[1] - 1)].fillna(0)
        populations_df = populations_df.interpolate(axis=1)
    else:
        populations_df = populations_df.fillna(0)

    populations_df = populations_df.loc[ordering]
    val, count = np.unique(populations_df.index, return_counts=True)
    doubles = val[count > 1]
    populations_df.loc[doubles] = populations_df.loc[doubles] / 2

    if not args.absolute:
        populations_df_sum = populations_df.sum(axis=0)
        populations_df_sum[populations_df_sum == 0] = 1
        populations_df = populations_df / populations_df_sum

    populations_df = populations_df.values

    if args.smooth:
        populations_df = gaussian_filter(populations_df, (0, args.smooth))

    # Create colors
    colors = np.array(cm.rainbow(np.linspace(0, 1, len(ids))).tolist())
    np.random.shuffle(colors)
    colors = pd.DataFrame(colors, index=ids)
    colors.loc[root_id] = [0.5, 0.5, 0.5, 1]
    colors = colors.loc[ordering].values

    # Plot
    plt.figure(figsize=(20, 10))
    stackplot(np.arange(first_step, last_step + 1), populations_df, colors=colors,
                  aa=True, lw=1)

    plt.xlim(first_step, last_step)

    if args.absolute:
        plt.yscale('log')
    else:
        plt.ylim(0, 1)

    plt.savefig(args.output)
