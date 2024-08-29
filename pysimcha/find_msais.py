import argparse
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
import numpy as np

COL_ALLELE_A = mpl.colors.to_rgba('orange')
COL_ALLELE_B = mpl.colors.to_rgba('teal')
COL_NONMSAI = mpl.colors.to_rgba('gray')
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

end_clones_ids = [5, 6, 8, 9, 10, 14, 15, 16, 18, 19]
end_names = [f"refphase_{i}" for i in end_clones_ids]


def format_chromosomes_int(chroms):
    return chroms.astype(str).str.replace('chr', '').replace('X', '23').replace('Y', '24').astype(int)


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

def plot_single_copynumber(data, msai_locs, ax=None, print_chrom=False, ymax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(20, 5))
    data = data.copy()
    data['total_start'] = calc_total_pos(data, col='start')
    data['total_end'] = calc_total_pos(data, col='end')
    segment_lengths = data['end'] - data['start']
    data['small_segment'] = segment_lengths < SMALL_SEGMENTS_LIMIT
    msais = pd.DataFrame()
    msais['total_start'] = calc_total_pos(msai_locs, col='start')
    msais['total_end']   = calc_total_pos(msai_locs, col='end')

    lines_a = []
    lines_b = []
    circles_a = []
    col_circle_a = []
    circles_b = []
    col_circle_b = []
    lines_alpha = []
    col_lines_a = []
    col_lines_b = []

    for _, dat in data.iterrows():
        #print(dat['total_start'])
        #print(dat['total_start'] in msai_total_starts)
        if dat['total_start'] in msais['total_start'].values:
            if dat['small_segment']:
                col_circle_a.append(COL_ALLELE_A)
                col_circle_b.append(COL_ALLELE_B)
            col_lines_a.append(COL_ALLELE_A)
            col_lines_b.append(COL_ALLELE_B)
        else:
            if dat['small_segment']:
                col_circle_a.append(COL_NONMSAI)
                col_circle_b.append(COL_NONMSAI)
            col_lines_a.append(COL_NONMSAI)
            col_lines_b.append(COL_NONMSAI)
        if dat['small_segment']:
            
            circles_a.append(
                (dat['total_start'] + 0.5*(dat['total_end'] - dat['total_start']), dat['cn_a']))
            circles_b.append(
                (dat['total_start'] + 0.5*(dat['total_end'] - dat['total_start']), dat['cn_b']))
                
        lines_alpha.append(0.5 if dat['cn_a'] == dat['cn_b'] else 1)
        lines_a.append([(dat['total_start'], dat['cn_a']), (dat['total_end'], dat['cn_a'])])
        lines_b.append([(dat['total_start'], dat['cn_b']), (dat['total_end'], dat['cn_b'])])

    colors_a = np.array([COL_ALLELE_A] * len(lines_a))
    colors_a = np.array(col_lines_a)
    # clonal mutations
    colors_a[:, 3] = np.array(lines_alpha)

    colors_b = np.array([COL_ALLELE_B] * len(lines_b))
    colors_b = np.array(col_lines_b)
    colors_b[:, 3] = np.array(lines_alpha)
    # a and b are overlapping
    colors_a[data['cn_a'] == data['cn_b'], 3] = 0.5
    colors_b[data['cn_a'] == data['cn_b'], 3] = 0.5
    colors = np.row_stack([colors_a, colors_b])

    lc = mpl.collections.LineCollection(lines_a + lines_b, colors=colors,
                                        linewidth=LINEWIDTH_COPY_NUMBERS)
    ax.add_collection(lc)

    if len(circles_a) > 0:
        col_circle_a = np.array(col_circle_a)
        col_circle_b = np.array(col_circle_b)

        a_b_overlap = (data.loc[data['small_segment']]['cn_a'] ==
                       data.loc[data['small_segment']]['cn_b']).values

        ax.scatter(np.array(circles_a)[:, 0], np.array(circles_a)[:, 1],
                marker='o', color='black', alpha=1., zorder=5)
        ax.scatter(np.array(circles_b)[:, 0], np.array(circles_b)[:, 1],
                marker='o', color='black', alpha=1., zorder=5)
        ax.scatter(np.array(circles_a)[~a_b_overlap, 0], np.array(circles_a)[~a_b_overlap, 1],
                marker='o', color=np.array(col_circle_a)[~a_b_overlap], alpha=1., zorder=6)
        ax.scatter(np.array(circles_a)[a_b_overlap, 0], np.array(circles_a)[a_b_overlap, 1],
                marker='o', color=np.array(col_circle_a)[a_b_overlap], alpha=0.5, zorder=6)
        ax.scatter(np.array(circles_b)[~a_b_overlap, 0], np.array(circles_b)[~a_b_overlap, 1],
                marker='o', color=np.array(col_circle_b)[~a_b_overlap], alpha=1., zorder=6)
        ax.scatter(np.array(circles_b)[a_b_overlap, 0], np.array(circles_b)[a_b_overlap, 1],
                marker='o', color=np.array(col_circle_b)[a_b_overlap], alpha=0.5, zorder=6)

    if ymax is None:
        ymax = data[['cn_a', 'cn_b']].max().max()
    ax.set_ylim(-0.1, ymax + 0.1)
    plot_chrom_boundaries(ax, print_chrom=print_chrom)
    #ax.set_yticks(range(ymax+1))
    ax.set_ylabel(data["sample_id"][0][-2:])
    #seg_bound_first = data['total_start'].values[0]
    #seg_bounds = data['total_end'].values
    #ax.vlines(np.append(seg_bound_first, seg_bounds), ymin=0, ymax=ymax, ls='--', alpha=0.5,
    #          color='grey', linewidth=1)


# Plotting functionality to make it look like Refphase's paper
def plot_msais(data, msai_locs):
    n_samples = len(data) 
    fig, axes = plt.subplots(n_samples,1,figsize=(24, 60), constrained_layout=True)
    for i, df in enumerate(data):
        plot_single_copynumber(df, msai_locs, ax=axes[i])
    
    plt.show()
    fig.savefig("./out/cn.png")
    return

# Find all MSAIs in a given phylogenetic tree
def find_msais(data, output):
    # Select the refphase 0  segments
    seg_0 = data[data["sample_id"]=="refphase_0"]
    segments = []
    for name in end_names:
        segments.append(data[data["sample_id"] == name].reset_index())

    n_msais = 0
    where_msais = []

    for index, row in seg_0.iterrows():
        hierarchy = ["cn_a" for _ in range(len(end_names))]
        flipped = False
        # Compare the segments
        for j, clonal_seg in enumerate(segments):
            info = clonal_seg.iloc[index]
            if info["cn_a"] < info["cn_b"]:
                hierarchy[j] = "cn_b"
        if len(set(hierarchy)) == 2:
            n_msais += 1
            where_msais.append([row["chr"],row["start"],row["end"]])
    msai_locs = pd.DataFrame(np.array(where_msais), columns=["chr","start", "end"])

    sum_msai_lengths = (msai_locs["end"] - msai_locs["start"]).sum()

    total_seg_lengths = (segments[0]["end"] - segments[0]["start"]).sum()
    
    frac = sum_msai_lengths/total_seg_lengths

    print(f"Number of MSAI segments: {n_msais};   total number of segments: {len(segments[0])};   fractional length of MSAIs: {frac}")
    f = open(output, "a")
    f.write(f"{frac}\n")
    #plot_msais(segments, msai_locs)


        # Get all the relevant segments
        #all_segs = data[data["start"] == seg["start"]]
        #print(all_segs)
    
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find MSAIs, i.e. mirrored subclonal imbalances")
    parser.add_argument("-i","--input", default="../out/copynumbers.tsv", help="The file with the karyotype data.")
    parser.add_argument("--sample", default="sample_1", help="Sample ID")
    parser.add_argument("-o", "--output", default="./out/karyotype.png", help="Output file path")
    parser.add_argument("--joint", action="store_true", help="Plot both haplotypes jointly (default: False))")
    args = parser.parse_args()

    data = pd.read_csv(args.input, sep="\t")
    find_msais(data, args.output)
