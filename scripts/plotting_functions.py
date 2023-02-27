"""
Plotting functions largely taken from MEDICC2 (see https://bitbucket.org/schwarzlab/medicc2)
"""

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from utils import format_chromosomes, format_chromosomes_int

COL_ALLELE_A = mpl.colors.to_rgba('orange')
COL_ALLELE_B = mpl.colors.to_rgba('teal')
LINEWIDTH_COPY_NUMBERS = 4
SMALL_SEGMENTS_LIMIT = 1e7
COL_VLINES = '#1f77b4'
COL_MARKER_INTERNAL = COL_VLINES
COL_MARKER_TERMINAL = 'black'
COL_MARKER_NORMAL = 'green'
LINEWIDTH_COPY_NUMBERS = 2
TREE_MARKER_SIZE = 40
XLABEL_TICK_SIZE = 8
XLABEL_FONT_SIZE = 10
YLABEL_FONT_SIZE = 12
SMALL_SEGMENTS_LIMIT = 1e7
TREE_WIDTH_SCALE = 1
TRACK_WIDTH_SCALE = 1
HEIGHT_SCALE = 1

CHROM_SIZES = [247249719, 242951149, 199501827, 191273063, 180857866, 170899992,
    158821424, 146274826, 140273252, 135374737, 134452384, 132349534,
    114142980, 106368585, 100338915, 88827254, 78774742, 76117153,
    63811651, 62435964, 46944323, 49691432, 154913754]

# sum = 0
# CHROM_ABS_START = {list(CHROM_SIZE.keys())[i]: (sum := sum+([0]+list(CHROM_SIZE.values()))[i])
#                    for i in range(len(CHROM_SIZE))}


def load_data(filename):
    '''Load the data from a file'''
    data = pd.read_csv(filename, sep='\t')
    assert np.all(data.columns == ['sample_id', 'chr', 'start', 'end',
                  'cn_a', 'cn_b']), 'Data columns are not correct'
    data['chr'] = format_chromosomes_int(data['chr'])
    data = data.set_index(['sample_id', 'chr', 'start', 'end'])
    data = data.sort_index()

    assert data.index.is_unique, 'Index is not unique'
    assert len(data.columns) == 2, 'Data has more than two columns'

    return data

def add_chrom_offset(data):
    chrom_offset = pd.DataFrame([(f'chr{str(i+1).replace("23", "X")}', l)
                                    for i, l in enumerate(CHROM_SIZES)], columns=['chr', 'chrom_offset'])
    chrom_offset['chrom_offset'] = np.append(0, np.cumsum(chrom_offset['chrom_offset'])[:-1])
    chrom_offset['chr'] = format_chromosomes_int(chrom_offset['chr'])
    chrom_offset = chrom_offset.set_index('chr')[['chrom_offset']]
    data = data.join(chrom_offset, on='chr').copy()

    return data


def calc_total_pos(data, col='pos'):
    if 'chrom_offset' not in data.columns:
        data = add_chrom_offset(data)
    return data[col] + data['chrom_offset']


def plot_chrom_boundaries(ax, print_chrom=False, outside=False):
    chrom_offset = pd.DataFrame([(f'chr{str(i+1).replace("23", "X")}', l)
                                for i, l in enumerate(CHROM_SIZES)], columns=['chr', 'chrom_offset'])
    chrom_offset['chrom_offset'] = np.append(0, np.cumsum(chrom_offset['chrom_offset'])[:-1])
    chrom_offset = chrom_offset.set_index('chr')['chrom_offset']
    for chrom, chrom_offset in chrom_offset.items():
        ax.axvline(chrom_offset, color='grey')
        if print_chrom:
            if outside:
                ax.text(chrom_offset + 1e7, ax.get_ylim()[1],
                        chrom.replace('chr', ''), fontsize=8, ha='left', va='bottom')
            else:
                ax.text(chrom_offset + 1e7, ax.get_ylim()[0] + (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.9,
                        chrom.replace('chr', ''), fontsize=8, ha='left')

    ax.set_xlim(0, np.sum(CHROM_SIZES))


def plot_single_baf(data, ax=None, print_chrom=False):
    if ax is None:
        fig, ax = plt.subplots(figsize=(20, 5))

    sample = data.columns[-1]
    data['total_pos'] = calc_total_pos(data)

    ax.axhline(0.5, linestyle='--', color='black')
    ax.set_ylim(-0.01, 1.01)
    ax.plot(data['total_pos'], data[sample], 'o', color='grey')
    plot_chrom_boundaries(ax, print_chrom=print_chrom)
    ax.set_ylabel('BAF')


def plot_single_logr(data, ax=None, print_chrom=False):
    if ax is None:
        fig, ax = plt.subplots(figsize=(20, 5))

    sample = data.columns[-1]
    data['total_pos'] = calc_total_pos(data)

    ax.axhline(0.0, linestyle='--', color='black')
    ax.plot(data['total_pos'], data[sample], 'o', color='grey')
    plot_chrom_boundaries(ax, print_chrom=print_chrom)
    ax.set_ylabel('logR')


def plot_single_copynumber(data, ax=None, print_chrom=False, ymax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(20, 5))
    data = data.copy()
    data['total_start'] = calc_total_pos(data, col='start')
    data['total_end'] = calc_total_pos(data, col='end')
    segment_lengths = data['end'] - data['start']
    data['small_segment'] = segment_lengths < SMALL_SEGMENTS_LIMIT

    lines_a = []
    lines_b = []
    circles_a = []
    circles_b = []
    lines_alpha = []

    for _, dat in data.iterrows():
        if dat['small_segment']:
            circles_a.append(
                (dat['total_start'] + 0.5*(dat['total_end'] - dat['total_start']), dat['cn_a']))
            circles_b.append(
                (dat['total_start'] + 0.5*(dat['total_end'] - dat['total_start']), dat['cn_b']))

        lines_alpha.append(0.5 if dat['cn_a'] == dat['cn_b'] else 1)
        lines_a.append([(dat['total_start'], dat['cn_a']), (dat['total_end'], dat['cn_a'])])
        lines_b.append([(dat['total_start'], dat['cn_b']), (dat['total_end'], dat['cn_b'])])

    colors_a = np.array([COL_ALLELE_A] * len(lines_a))
    # clonal mutations
    colors_a[:, 3] = np.array(lines_alpha)

    colors_b = np.array([COL_ALLELE_B] * len(lines_b))
    colors_b[:, 3] = np.array(lines_alpha)
    # a and b are overlapping
    colors_a[data['cn_a'] == data['cn_b'], 3] = 0.5
    colors_b[data['cn_a'] == data['cn_b'], 3] = 0.5
    colors = np.row_stack([colors_a, colors_b])

    lc = mpl.collections.LineCollection(lines_a + lines_b, colors=colors,
                                        linewidth=LINEWIDTH_COPY_NUMBERS)
    ax.add_collection(lc)

    if len(circles_a) > 0:
        a_b_overlap = (data.loc[data['small_segment']]['cn_a'] ==
                       data.loc[data['small_segment']]['cn_b']).values

        ax.plot(np.array(circles_a)[:, 0], np.array(circles_a)[:, 1],
                'o', ms=4, color='black', alpha=1., zorder=5)
        ax.plot(np.array(circles_b)[:, 0], np.array(circles_b)[:, 1],
                'o', ms=4, color='black', alpha=1., zorder=5)
        ax.plot(np.array(circles_a)[~a_b_overlap, 0], np.array(circles_a)[~a_b_overlap, 1],
                'o', ms=3, color=COL_ALLELE_A, alpha=1., zorder=6)
        ax.plot(np.array(circles_a)[a_b_overlap, 0], np.array(circles_a)[a_b_overlap, 1],
                'o', ms=3, color=COL_ALLELE_A, alpha=0.5, zorder=6)
        ax.plot(np.array(circles_b)[~a_b_overlap, 0], np.array(circles_b)[~a_b_overlap, 1],
                'o', ms=3, color=COL_ALLELE_B, alpha=1., zorder=6)
        ax.plot(np.array(circles_b)[a_b_overlap, 0], np.array(circles_b)[a_b_overlap, 1],
                'o', ms=3, color=COL_ALLELE_B, alpha=0.5, zorder=6)

    if ymax is None:
        ymax = data[['cn_a', 'cn_b']].max().max()
    ax.set_ylim(-0.1, ymax + 0.1)
    plot_chrom_boundaries(ax, print_chrom=print_chrom)
    ax.set_yticks(range(ymax+1))

    seg_bound_first = data['total_start'].values[0]
    seg_bounds = data['total_end'].values
    ax.vlines(np.append(seg_bound_first, seg_bounds), ymin=0, ymax=ymax, ls='--', alpha=0.5,
              color='grey', linewidth=1)


def _get_x_positions(tree):
    """Create a mapping of each clade to its horizontal position.
    Dict of {clade: x-coord}
    """
    depths = tree.depths()
    # If there are no branch lengths, assume unit branch lengths
    if not max(depths.values()):
        depths = tree.depths(unit_branch_lengths=True)
    return depths


def _get_y_positions(tree, normal_name, adjust=False):
    """Create a mapping of each clade to its vertical position.
    Dict of {clade: y-coord}.
    Coordinates are negative, and integers for tips.
    """
    maxheight = tree.count_terminals()
    heights = {tip: maxheight - 1 - i for i,
               tip in enumerate(reversed([x for x in tree.get_terminals() if x.name != normal_name]))}
    heights.update({list(tree.find_clades(normal_name))[0]: -1})

    # Internal nodes: place at midpoint of children
    def calc_row(clade):
        for subclade in clade:
            if subclade not in heights:
                calc_row(subclade)
        # Closure over heights
        heights[clade] = (heights[clade.clades[0]] + heights[clade.clades[-1]]) / 2.0

    if tree.root.clades:
        calc_row(tree.root)

    if adjust:
        pos = pd.DataFrame([(clade, val) for clade, val in heights.items()],
                           columns=['clade', 'pos']).sort_values('pos')
        pos['newpos'] = 0
        count = 0
        for i in pos.index:
            if pos.loc[i, 'clade'].name is not None and pos.loc[i, 'clade'].name != 'root':
                count = count+1
            pos.loc[i, 'newpos'] = count

        pos.set_index('clade', inplace=True)
        heights = pos.to_dict()['newpos']

    return heights


def plot_tree(input_tree,
              label_func=None,
              title='',
              ax=None,
              output_name=None,
              normal_name=0,
              show_branch_lengths=True,
              branch_labels=None,
              label_colors=None,
              marker_size=None,
              line_width=None,
              hide_internal_nodes=False,
              **kwargs):
    """Plot the given tree using matplotlib (or pylab).
    The graphic is a rooted tree, drawn with roughly the same algorithm as
    draw_ascii.
    Additional keyword arguments passed into this function are used as pyplot
    options. The input format should be in the form of:
    pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict), or
    pyplot_option_name=(dict).
    Example using the pyplot options 'axhspan' and 'axvline'::
        from Bio import Phylo, AlignIO
        from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
        constructor = DistanceTreeConstructor()
        aln = AlignIO.read(open('TreeConstruction/msa.phy'), 'phylip')
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)
        tree = constructor.upgma(dm)
        Phylo.draw(tree, axhspan=((0.25, 7.75), {'facecolor':'0.5'}),
        ... axvline={'x':0, 'ymin':0, 'ymax':1})
    Visual aspects of the plot can also be modified using pyplot's own functions
    and objects (via pylab or matplotlib). In particular, the pyplot.rcParams
    object can be used to scale the font size (rcParams["font.size"]) and line
    width (rcParams["lines.linewidth"]).
    :Parameters:
        label_func : callable
            A function to extract a label from a node. By default this is str(),
            but you can use a different function to select another string
            associated with each node. If this function returns None for a node,
            no label will be shown for that node.
        do_show : bool
            Whether to show() the plot automatically.
        show_support : bool
            Whether to display confidence values, if present on the tree.
        ax : matplotlib/pylab axes
            If a valid matplotlib.axes.Axes instance, the phylogram is plotted
            in that Axes. By default (None), a new figure is created.
        branch_labels : dict or callable
            A mapping of each clade to the label that will be shown along the
            branch leading to it. By default this is the confidence value(s) of
            the clade, taken from the ``confidence`` attribute, and can be
            easily toggled off with this function's ``show_support`` option.
            But if you would like to alter the formatting of confidence values,
            or label the branches with something other than confidence, then use
            this option.
        label_colors : dict or callable
            A function or a dictionary specifying the color of the tip label.
            If the tip label can't be found in the dict or label_colors is
            None, the label will be shown in black.
    """

    import matplotlib.collections as mpcollections

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 7))

    label_func = label_func if label_func is not None else lambda x: x

    # options for displaying label colors.
    if label_colors is not None:
        if callable(label_colors):
            def get_label_color(label):
                return label_colors(label)
        else:
            # label_colors is presumed to be a dict
            def get_label_color(label):
                return label_colors.get(label, "black")
    else:
        clade_colors = {}
        for sample in [x.name for x in list(input_tree.find_clades('')) if x.name is not None]:
            # determine if sample is terminal
            is_terminal = True
            matches = list(input_tree.find_clades(sample))
            if len(matches) > 0:
                clade = matches[0]
                is_terminal = clade.is_terminal()
            # determine if sample is normal
            clade_colors[sample] = COL_MARKER_TERMINAL
            if not is_terminal:
                clade_colors[sample] = COL_MARKER_INTERNAL
            if sample == normal_name:
                clade_colors[sample] = COL_MARKER_NORMAL

        def get_label_color(label):
            return clade_colors.get(label, "black")

    if marker_size is None:
        marker_size = TREE_MARKER_SIZE

    def marker_func(x): return (marker_size, get_label_color(
        x.name)) if x.name is not None else None

    ax.axes.get_yaxis().set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True, prune=None))
    ax.xaxis.set_tick_params(labelsize=XLABEL_TICK_SIZE)
    ax.xaxis.label.set_size(XLABEL_FONT_SIZE)
    ax.set_title(title, x=0.01, y=1.0, ha='left', va='bottom',
                 fontweight='bold', fontsize=16, zorder=10)
    x_posns = _get_x_positions(input_tree)
    y_posns = _get_y_positions(input_tree, adjust=not hide_internal_nodes, normal_name=normal_name)
    # Arrays that store lines for the plot of clades
    horizontal_linecollections = []
    vertical_linecollections = []

    # Options for displaying branch labels / confidence
    def value_to_str(value):
        if value is None or value == 0:
            return None
        elif int(value) == value:
            return str(int(value))
        else:
            return str(value)

    if not branch_labels:
        if show_branch_lengths:
            def format_branch_label(x):
                return value_to_str(np.round(x.branch_length, 1)) if x.name != 'root' and x.name is not None else None
        else:
            def format_branch_label(clade):
                return None

    elif isinstance(branch_labels, dict):
        def format_branch_label(clade):
            return branch_labels.get(clade)
    else:
        if not callable(branch_labels):
            raise TypeError(
                "branch_labels must be either a dict or a callable (function)"
            )

        def format_branch_label(clade):
            return value_to_str(branch_labels(clade))

    def draw_clade_lines(
        use_linecollection=False,
        orientation="horizontal",
        y_here=0,
        x_start=0,
        x_here=0,
        y_bot=0,
        y_top=0,
        color="black",
        lw=".1",
    ):
        """Create a line with or without a line collection object.
        Graphical formatting of the lines representing clades in the plot can be
        customized by altering this function.
        """
        if not use_linecollection and orientation == "horizontal":
            ax.hlines(y_here, x_start, x_here, color=color, lw=lw)
        elif use_linecollection and orientation == "horizontal":
            horizontal_linecollections.append(
                mpcollections.LineCollection(
                    [[(x_start, y_here), (x_here, y_here)]], color=color, lw=lw
                )
            )
        elif not use_linecollection and orientation == "vertical":
            ax.vlines(x_here, y_bot, y_top, color=color)
        elif use_linecollection and orientation == "vertical":
            vertical_linecollections.append(
                mpcollections.LineCollection(
                    [[(x_here, y_bot), (x_here, y_top)]], color=color, lw=lw
                )
            )

    def draw_clade(clade, x_start, color, lw):
        """Recursively draw a tree, down from the given clade."""
        x_here = x_posns[clade]
        y_here = y_posns[clade]
        # phyloXML-only graphics annotations
        if hasattr(clade, "color") and clade.color is not None:
            color = clade.color.to_hex()
        if hasattr(clade, "width") and clade.width is not None:
            lw = clade.width * plt.rcParams["lines.linewidth"]
        # Draw a horizontal line from start to here
        draw_clade_lines(
            use_linecollection=True,
            orientation="horizontal",
            y_here=y_here,
            x_start=x_start,
            x_here=x_here,
            color=color,
            lw=lw,
        )
        # Add node marker
        if marker_func is not None:
            marker = marker_func(clade)
            if marker is not None and clade is not None and \
                    not(hide_internal_nodes and not clade.is_terminal()):
                marker_size, marker_col = marker_func(clade)
                ax.scatter(x_here, y_here, s=marker_size, c=marker_col, zorder=3)

        # Add node/taxon labels
        label = label_func(str(clade))
        ax_scale = ax.get_xlim()[1] - ax.get_xlim()[0]

        if label not in (None, clade.__class__.__name__) and \
                not (hide_internal_nodes and not clade.is_terminal()):
            ax.text(
                x_here + min(0.5*ax_scale, 1),
                y_here,
                " %s" % label,
                verticalalignment="center",
                color=get_label_color(label),
            )
        # Add label above the branch
        conf_label = format_branch_label(clade)
        if conf_label:
            ax.text(
                0.5 * (x_start + x_here),
                y_here - 0.15,
                conf_label,
                fontsize="small",
                horizontalalignment="center",
            )

        if clade.clades:
            # Draw a vertical line connecting all children
            all_children_y = [y_posns[x] for x in clade.clades if x.name is not None] + \
                ([y_posns[clade]] if clade.name is not None else [])
            y_top = np.max(all_children_y)
            y_bot = np.min(all_children_y)
            # Only apply widths to horizontal lines, like Archaeopteryx
            draw_clade_lines(
                use_linecollection=True,
                orientation="vertical",
                x_here=x_here,
                y_bot=y_bot,
                y_top=y_top,
                color=color,
                lw=lw,
            )
            # Draw descendents
            for child in clade:
                draw_clade(child, x_here, color, lw)

    if line_width is None:
        line_width = plt.rcParams["lines.linewidth"]
    draw_clade(input_tree.root, 0, "k", line_width)

    # If line collections were used to create clade lines, here they are added
    # to the pyplot plot.
    for i in horizontal_linecollections:
        ax.add_collection(i)
    for i in vertical_linecollections:
        ax.add_collection(i)

    ax.set_xlabel("branch length")
    ax.set_ylabel("taxa")

    # Add margins around the `tree` to prevent overlapping the ax
    xmax = max(x_posns.values())
    #ax.set_xlim(-0.05 * xmax, 1.25 * xmax)
    ax.set_xlim(-0.05 * xmax, 1.05 * xmax)
    # Also invert the y-axis (origin at the top)
    # Add a small vertical margin, but avoid including 0 and N+1 on the y axis
    ax.set_ylim(max(y_posns.values()) + 0.5, 0.5)

    # Parse and process key word arguments as pyplot options
    for key, value in kwargs.items():
        try:
            # Check that the pyplot option input is iterable, as required
            list(value)
        except TypeError:
            raise ValueError(
                'Keyword argument "%s=%s" is not in the format '
                "pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict),"
                " or pyplot_option_name=(dict) " % (key, value)
            ) from None
        if isinstance(value, dict):
            getattr(plt, str(key))(**dict(value))
        elif not (isinstance(value[0], tuple)):
            getattr(plt, str(key))(*value)
        elif isinstance(value[0], tuple):
            getattr(plt, str(key))(*value[0], **dict(value[1]))

    if output_name is not None:
        plt.savefig(output_name + ".png")

    return plt.gcf()


def plot_cn_bars(copynumbers, tree=None, fraction=False, y_posns=None, cmax=8, tree_width_ratio=1,
                 hide_internal_nodes=True):

    copynumbers = copynumbers.copy().reset_index().set_index('sample_id')
    # Display the population size as fractions
    if fraction:
        total_pop = 0
        for clade in tree.find_clades():
            if clade.name is not None:
                total_pop += int(clade.name.split('-')[1])
        for clade in tree.find_clades():
            if clade.name is not None:
                clade.name = clade.name.split('-')[0] + '-' + \
                    str(np.round(float(clade.name.split('-')[1])/total_pop, 3))

    samples = list(copynumbers.index.unique())

    # Figure dimensions
    nsegs = len(copynumbers.loc[samples[0]])
    track_width = nsegs * 0.2 * TRACK_WIDTH_SCALE
    tree_width = 0 if tree is None else 2.5 * TREE_WIDTH_SCALE  # in figsize
    nrows = len(samples)
    plotheight = 4 * 0.2 * nrows * HEIGHT_SCALE
    plotwidth = tree_width + track_width
    tree_width_ratio = tree_width / plotwidth
    fig = plt.figure(figsize=(plotwidth, plotheight), constrained_layout=True)

    if tree is None:
        gs = fig.add_gridspec(nrows, 1)
        cn_axes = [fig.add_subplot(gs[i]) for i in range(0, nrows)]
        y_posns = {sample: i for i, sample in enumerate(samples)}  # as they appear
    else:
        normal_sample = tree.root.name
        gs = fig.add_gridspec(nrows, 2, width_ratios=[tree_width_ratio, 1-tree_width_ratio])
        tree_ax = fig.add_subplot(gs[0:len(samples), 0])
        cn_axes = [fig.add_subplot(gs[i]) for i in range(1, (2*(nrows))+1, 2)]
        y_posns = _get_y_positions(tree, normal_sample, adjust=True)
        y_posns = {clade.name: y_pos for clade, y_pos in y_posns.items()
                   if clade.name is not None and clade.name != 'root'}
        plot_tree(tree,
                  ax=tree_ax,
                  label_func=lambda x: '',
                  normal_name=normal_sample,
                  hide_internal_nodes=hide_internal_nodes)

    for sample in samples:
        cur_data = copynumbers.loc[sample]
        index_to_plot = y_posns[sample] - 1
        plot_single_copynumber(cur_data,
                               ax=cn_axes[index_to_plot],
                               ymax=cmax)
        cn_axes[index_to_plot].set_ylabel(
            sample, fontsize=YLABEL_FONT_SIZE, rotation=0, va='center', ha='right')
    for ax in cn_axes[:-1]:
        ax.set_xticks([])

    return fig


def plot_cn_heatmap(copynumbers, tree=None, y_posns=None, cmax=8, total_copy_numbers=False,
                    alleles=['cn_a', 'cn_b'], tree_width_ratio=1, cbar_width_ratio=0.05, figsize=(20, 10),
                    tree_line_width=0.5, tree_marker_size=1, hide_internal_nodes=True, title='',
                    tree_label_colors=None, cmap='coolwarm', normal_name=0,
                    ignore_segment_lengths=False):

    copynumbers = copynumbers[alleles].copy()

    if not isinstance(alleles, list) and not isinstance(alleles, tuple):
        alleles = [alleles]
    nr_alleles = len(alleles)

    cmax = min(cmax, np.max(copynumbers[alleles].values.astype(int)))

    if tree is None:
        fig, axs = plt.subplots(figsize=figsize, ncols=1+nr_alleles, sharey=False,
                                gridspec_kw={'width_ratios': nr_alleles*[1] + [cbar_width_ratio]})

        cur_sample_labels = (copynumbers.index.get_level_values('sample_id').unique())
        if y_posns is None:
            y_posns = {s: i for i, s in enumerate(cur_sample_labels)}

        cn_axes = axs[:-1]
    else:
        fig, axs = plt.subplots(figsize=figsize, ncols=2+nr_alleles, sharey=False,
                                gridspec_kw={'width_ratios': [tree_width_ratio] + nr_alleles*[1] + [cbar_width_ratio]})
        tree_ax = axs[0]
        cn_axes = axs[1:-1]

        if not hide_internal_nodes:
            cur_sample_labels = np.array([x.name for x in list(
                tree.find_clades()) if x.name is not None])
        else:
            cur_sample_labels = np.array([x.name for x in tree.get_terminals()])

        y_posns = {k.name: v for k, v in _get_y_positions(
            tree, adjust=not hide_internal_nodes, normal_name=normal_name).items()}

        _ = plot_tree(tree, ax=tree_ax, normal_name=normal_name,
                      label_func=lambda x: '', hide_internal_nodes=hide_internal_nodes,
                      show_branch_lengths=False, line_width=tree_line_width, marker_size=tree_marker_size,
                      title=title, label_colors=tree_label_colors)
        tree_ax.set_axis_off()
        tree_ax.set_axis_off()
        fig.set_constrained_layout_pads(w_pad=0, h_pad=0, hspace=0.0, wspace=100)
        tree_ax.set_ylim(len(cur_sample_labels)-0.5, -0.5)

    cax = axs[-1]

    ind = [y_posns.get(x, -1) for x in cur_sample_labels]
    cur_sample_labels = cur_sample_labels[np.argsort(ind)]
    if cmax == (2 if total_copy_numbers else 1):
        cmax += 1
    color_norm = mcolors.TwoSlopeNorm(vmin=0, vcenter=(2 if total_copy_numbers else 1), vmax=cmax)

    chr_ends = copynumbers.loc[cur_sample_labels[0]].copy()
    chr_ends['end_pos'] = np.cumsum([1]*len(chr_ends))
    chr_ends = chr_ends.reset_index().groupby('chr').max()['end_pos']
    chr_ends = chr_ends.dropna()
    chr_ends = chr_ends.astype(int)

    if ignore_segment_lengths:
        x_pos = np.arange(
            len(copynumbers.loc[cur_sample_labels].astype(int).unstack('sample_id'))+1)
    else:
        x_pos = np.append([0], np.cumsum(copynumbers.loc[cur_sample_labels].astype(int).unstack(
            'sample_id').loc[:, (alleles[0])].loc[:, cur_sample_labels].eval('end-start').values))
    y_pos = np.arange(len(cur_sample_labels)+1)

    for ax, allele in zip(cn_axes, alleles):
        im = ax.pcolormesh(x_pos, y_pos,
                           copynumbers.loc[cur_sample_labels].astype(int).unstack(
                               'sample_id').loc[:, (allele)].loc[:, cur_sample_labels].values.T,
                           cmap=cmap,
                           norm=color_norm)

        for _, line in chr_ends.iteritems():
            ax.axvline(x_pos[line], color='black', linewidth=0.75)

        xtick_pos = np.append([0], x_pos[chr_ends.values][:-1])
        xtick_pos = (xtick_pos + np.roll(xtick_pos, -1))/2
        xtick_pos[-1] += x_pos[-1]/2
        ax.set_xticks(xtick_pos)
        ax.set_xticklabels([str(x).replace('chr', '').replace('23', 'X').replace('24', 'Y') for x in chr_ends.index], ha='center', rotation=90, va='bottom')
        ax.tick_params(width=0)
        ax.xaxis.set_tick_params(labelbottom=False, labeltop=True, bottom=False)
        ax.set_yticks([])

    # add sample names
    cn_axes[0].set_yticks(np.arange(len(cur_sample_labels)) + 0.5)
    cn_axes[0].set_yticklabels(cur_sample_labels, ha='right',
                               va='center', fontsize=YLABEL_FONT_SIZE)

    # add colorbar
    cax.pcolormesh([0, 1],
                   np.arange(0, cmax+2),
                   np.arange(0, cmax+1)[:, np.newaxis],
                   cmap=cmap,
                   norm=color_norm)

    cax.set_xticks([])
    cax.set_yticks(np.arange(0, cmax+1)+0.5)
    if np.max(copynumbers.values.astype(int)) > cmax:
        cax.set_yticklabels([str(x) + '+' if x == cmax else str(x)
                             for x in np.arange(0, cmax+1)], ha='left')
    else:
        cax.set_yticklabels(np.arange(0, cmax+1), ha='left')
    cax.yaxis.set_tick_params(left=False, labelleft=False, labelright=True)

    # for ax in axs[:-1]:
    #     for i in range(len(cur_sample_labels)+1):
    #         ax.axhline(i, color='black', linewidth=0.5)
    #         ax.axhline(i+0.5, ls='--', color='black', linewidth=0.5)
    
    for ax in cn_axes:
        ax.set_ylim(len(cur_sample_labels), 0)

    return fig
