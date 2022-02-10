#!/usr/bin/env python
"""Function to create fishplots."""

import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import cm
from scipy.ndimage import gaussian_filter


def stackplot(x, *args, ax=None, colors=None, labels=(), **kwargs):
    """Draw a stacked area plot.

    Taken largely from the matplotlib implementation except that the keywords `edgecolor` and
    `where` are provided to `fill_between`
    """
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


def fish_plot(steps, pops_stack, colors=None, pop_max=None):
    """Plot the actual fish plot."""
    stackplot(steps, pops_stack, colors=colors)

    first_step = pops_stack.columns[0]
    last_step = pops_stack.columns[-1]
    plt.xlim(first_step, last_step)
    plt.ylim(0, 1)
    label_text = np.round(np.abs(np.arange(-0.5, 0.6, 0.1)), 1)
    plt.yticks(np.arange(0, 1.1, step=0.1),
               (label_text * pop_max).astype(int) if pop_max is not None else label_text)
    plt.xticks(np.arange(first_step, last_step + 1, step=(last_step - first_step) / 10).astype(int))


def _create_ordering(cur_tree, cur_clone):
    """Create index for population DataFrame.

    Recursively traverses the parent tree.
    Build a list such that children are listed between two instances of parent.
    """
    res = []
    if cur_clone in cur_tree.keys():
        res += [cur_clone]
        for child in cur_tree[cur_clone]:
            res += _create_ordering(cur_tree, child)
        res += [cur_clone]
    else:
        return [cur_clone]
    return res


def _build_tree(parent_tree):
    """Create a dict-based tree where for each parent there is a list of children."""
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


def _create_colors(ids, root_id, ordering, seed=0, cmap="rainbow"):
    np.random.seed(seed)
    try:
        cmap = cm.get_cmap(cmap)
    except ValueError:
        print("WARNING: colormap not recognized, setting to default")
        cmap = cm.get_cmap("rainbow")

    colors = np.array(cmap(np.linspace(0, 1, len(ids))))
    np.random.shuffle(colors)
    colors = pd.DataFrame(colors, index=ids)
    colors.loc[-1] = np.ones(4)
    colors.loc[root_id] = .5 * np.ones(4)
    colors = colors.loc[ordering].values

    return colors


def load_data(populations_file, parent_tree,
              first_step=None, last_step=None, interpolation=0, absolute=False, smooth=None, seed=0,
              cmap="rainbow"):
    """Load data required for plotting.

    Args:
        populations_file (str): Location of populations file
        parent_tree (str): Location of parents file
        first_step (int, optional): First step to plot. Defaults to None.
        last_step (int, optional): Last step to plot. Defaults to None.
        interpolation (int, optional): Order of interpolation. Defaults to 0.
        absolute (bool, optional): Does not normalize data if set True. Defaults to False.
        smooth (int, optional): Window for Gaussian smoothing. Defaults to None.
        seed (int, optional): Seed used for coloring. Defaults to 0.
        cmap (str, optional): Matplotlib colormap. Defaults to None.

    Returns:
        pd.DataFrame: DataFrame with population information
        pd.Series: Series used for x-axis
        list: List of colors
        int: Maximum population. None if not absolute
    """
    # population dataframe
    populations_df = pd.read_csv(populations_file)

    min_step = populations_df["Step"].min()
    first_step = max(first_step, min_step) if first_step else min_step

    max_step = populations_df["Step"].max()
    last_step = min(last_step, max_step) if last_step else max_step

    populations_df = populations_df.loc[(populations_df["Step"] >= first_step) &
                                        (populations_df["Step"] <= last_step)]

    # populations dataframe processing
    pops_table = populations_df.set_index(['Id', 'Step'])['Pop'].unstack('Step')

    if interpolation > 0:
        if pops_table.shape[1] > 50:
            print("WARNING: interpolation is not recommened for large data")
        pops_table[first_step].fillna(0, inplace=True)
        pops_table[last_step].fillna(0, inplace=True)

        if (~pops_table.isna()).sum(axis=1).min() - 1 < interpolation:
            raise ValueError(f"For interpolation order {interpolation}, the iput data has not "
                             f"enough datapoints (at least {interpolation + 1} per sample)")

        larger_pops_table = pd.DataFrame(index=pops_table.index,
                                         columns=np.arange(0, pops_table.shape[1] - 0.1, 0.1))
        larger_pops_table[pops_table.columns.astype(float)] = pops_table
        larger_pops_table = larger_pops_table.astype(float)

        linear_interpolation = larger_pops_table.interpolate(axis=1, method='linear')
        pops_table_interpolate = larger_pops_table.interpolate(axis=1, method='spline',
                                                               order=interpolation)
        pops_table_interpolate[linear_interpolation <= 0] = 0

        pops_table = pops_table_interpolate

        # pops_table = pops_table.interpolate(axis=1)
    else:
        pops_table = pops_table.fillna(0)

    steps = pops_table.columns

    if smooth:
        pops_table = pd.DataFrame(gaussian_filter(pops_table, (0, smooth)),
                                  index=pops_table.index, columns=pops_table.columns)

    if absolute:
        pops_sums = pops_table.sum(axis=0)
        pop_max = pops_sums.max()
        pops_rest = pop_max - pops_sums
        pops_table.loc[-1] = pops_rest
    else:
        pop_max = None

    # Build parental relationship
    tree, ids, root_id = _build_tree(parent_tree)
    samples = populations_df['Id'].unique()
    ordering = _create_ordering(tree, -1 if absolute else root_id)
    ordering = [x for x in ordering if x in samples or x == -1 and absolute]

    pops_stack = pops_table.loc[ordering]
    val, count = np.unique(pops_stack.index, return_counts=True)
    doubles = val[count > 1]
    pops_stack.loc[doubles] = pops_stack.loc[doubles] / 2

    pops_sum = pops_stack.sum(axis=0)
    pops_sum[pops_sum == 0] = 1
    pops_stack = pops_stack / (pops_sum)

    colors = _create_colors(ids, root_id, ordering, seed, cmap)
    return pops_stack, steps, colors, pop_max


def _parse_arguments():
    """Parse the input arguments."""
    randint = np.random.randint(0, 2**31 - 1)
    parser = argparse.ArgumentParser(description='Create a Fish (Muller) plot '
                                                 'for the given evolutionary tree.')
    parser.add_argument("populations", type=str,
                        help="A CSV file with the header \"Id,Step,Pop\".")
    parser.add_argument("parent_tree", type=str,
                        help="A CSV file with the header \"ParentId,ChildId\".")
    parser.add_argument("output", type=str,
                        help="Output image filepath. The format must support alpha channels.")
    parser.add_argument('-a', '--absolute', dest="absolute", action="store_true", default=False,
                        help='Plot the populations in absolute numbers rather than normalized.')
    parser.add_argument("-i", "--interpolation", dest='interpolation', type=int, default=0,
                        help="Order of interpolation for empty data (0 means no interpolation).")
    parser.add_argument("-S", "--smooth", type=float, default=None,
                        help="STDev for Gaussian convolutional filter. The higher the value "
                             "the smoother the resulting bands will be. Recommended is around 1.0.")
    parser.add_argument("-F", "--first", dest="first_step", type=int,
                        help="The step to start plotting from.")
    parser.add_argument("-L", "--last", dest="last_step", type=int,
                        help="The step to end the plotting at.")
    parser.add_argument("-R", "--seed", dest="seed", type=int,
                        help="Random seed for selection of colors.", default=randint)
    parser.add_argument("--cmap", type=str, default="rainbow",
                        help="Colormap to use. Has to be a matplotlib colormap Uses rainbow by default")
    parser.add_argument("-W", "--width", dest="width", type=int, default=1920,
                        help="Output image width")
    parser.add_argument("-H", "--height", dest="height", type=int, default=1080,
                        help="Output image height")

    args = parser.parse_args()

    return args


def setup_figure(width=1920, height=1080):
    """Create figure with specified height and width."""
    dpi = (width + height) // 20
    plt.figure(figsize=(width // dpi, height // dpi), dpi=dpi)


if __name__ == '__main__':

    args = _parse_arguments()

    pops_stack, steps, colors, pop_max = load_data(
        args.populations, args.parent_tree, args.first_step, args.last_step,
        args.interpolation, args.absolute, args.smooth, args.seed, args.cmap)

    # Plot
    setup_figure(args.width, args.height)
    fish_plot(steps, pops_stack, colors, pop_max)
    plt.savefig(args.output)
