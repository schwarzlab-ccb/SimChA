import re

import numpy as np
import pandas as pd
import json

from os.path import join

chromosome_names = [
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
    'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
    'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21',
    'chr22', 'chrX', 'chrY'
]

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


def format_chromosomes_int(chroms):
    return chroms.astype(str).str.replace('chr', '').replace('X', '23').replace('Y', '24').astype(int)


def format_chromosomes(ds):
    """ Expects pandas Series with chromosome names. 
    The goal is to take recognisalbe chromosome names, i.e. chr4 or chrom3 and turn them into chr3 format.
    If the chromosomes names are not recognized, return them unchanged."""
    ds = ds.astype('str')
    pattern = re.compile(r"(chr|chrom)?((\d+)|X|Y)", flags=re.IGNORECASE)
    matches = ds.apply(pattern.match)
    matchable = ~matches.isnull().any()
    if matchable:
        newchr = matches.apply(lambda x: "chr%s" % x[2].upper())
        numchr = matches.apply(lambda x: int(x[3]) if x[3] is not None else -1)
        chrlevels = np.sort(numchr.unique())
        chrlevels = np.setdiff1d(chrlevels, [-1])
        chrcats = ["chr%d" % i for i in chrlevels]
        if 'chrX' in list(newchr):
            chrcats += ['chrX', ]
        if 'chrY' in list(newchr):
            chrcats += ['chrY', ]
        newchr = pd.Categorical(newchr, categories=chrcats)
    else:
        print("Could not match the chromosome labels. Rename the chromosomes according chr1, "
              "chr2, ... to avoid potential errors.")
        newchr = pd.Categorical(ds, categories=ds.unique())

    return newchr


def get_min_maj_CNs(data, is_major):
    cols = data[['cn_a', 'cn_b']]
    return cols.max(axis=1) if is_major else cols.min(axis=1)

# calculate the average CN for each sample


def get_CNs_by_sample(data, is_major):
    CNs = data.groupby(['sample_id']).apply(get_min_maj_CNs, is_major)
    return CNs.groupby(['sample_id']).mean()


def get_hap_by_sample(data, hap):
    return data.groupby(['sample_id'])[hap].mean()


def load_dataset(dataset_path):
    loaded_data = {
        "CNs": pd.read_csv(join(dataset_path, "copynumbers.tsv"), index_col=0, sep="\t"),
        "clones": pd.read_csv(join(dataset_path, "clones.tsv"), index_col=0, sep="\t"),
        "karyotypes": pd.read_csv(join(dataset_path, "karyotypes.tsv"), index_col=0, sep="\t"),
        "samples": pd.read_csv(join(dataset_path, "samples.tsv"), index_col=0, sep="\t"),
    }
    config_file = join(dataset_path, "sim_params.json")
    # import json into a dict
    with open(config_file) as f:
        loaded_data["config"] = json.load(f)
    return loaded_data


def calc_CNs(dataset):
    cns = [get_CNs_by_sample(dataset, True), get_CNs_by_sample(dataset, False), get_hap_by_sample(dataset, "cn_a"), get_hap_by_sample(dataset, "cn_b")]
    df = pd.DataFrame(cns).T
    df.columns = ["major", "minor", "hap_a", "hap_b"]
    return df

def get_seg_lengths(data, include_cn_normal = True):
    sample_ids = pd.Series({c: data[c].unique() for c in data})["sample_id"]
    all_seg_length = []
    for _, id in enumerate(sample_ids):
        sample = data[data["sample_id"] == id]
        sample_seg_length = []
        for _, row in sample.iterrows():
            if row["cn_a"] == 1 and row["cn_b"] == 1 and not include_cn_normal:
                continue
            seg_length = row["end"] - row["start"]
            sample_seg_length.append(seg_length)
        all_seg_length.append(np.array(sample_seg_length).mean())
    return all_seg_length

def get_changepoints(data):
    sample_ids = pd.Series({c: data[c].unique() for c in data})["sample_id"]
    all_changepoints = []
    for _, id in enumerate(sample_ids):
        sample = data[data["sample_id"] == id]
        for chr in chromosome_names:
            chr_changepoints = []
            segs = sample[sample["chrom"] == chr]
            last_seg = 2
            for index, seg in segs.iterrows():
                this_seg = seg["cn_a"]+seg["cn_b"]
                chr_changepoints.append(abs(this_seg - last_seg))
                lastSeg = this_seg
            all_changepoints += chr_changepoints
    return all_changepoints

def get_BP_per_chromosome(data):
    sample_ids = pd.Series({c: data[c].unique() for c in data})["sample_id"]
    all_chr_bins = []
    for _, id in enumerate(sample_ids):
        sample = data[data["sample_id"] == id]
        for _, chr in enumerate(chromosome_names):
            segs = sample[sample["chrom"] == chr]
            all_chr_bins.append(len(segs)-1)
    return all_chr_bins

def get_BP_per_10MB(data):
    SIZE = 10000000
    sample_ids = pd.Series({c: data[c].unique() for c in data})["sample_id"]
    all_chr_bins = []
    for _, id in enumerate(sample_ids):
        sample = data[data["sample_id"] == id]
        for _, chr in enumerate(chromosome_names): 
            segs = sample[sample["chrom"] == chr]
            intervals = np.arange(0, hg19_chr_lengths[chr]+SIZE, SIZE)
            bins = [0 for _ in range(len(intervals)-1)]
            for _, seg in segs.iterrows():
                # To which bin does the start of the segment belong?
                [start_index,end_index] = np.searchsorted(intervals, [seg["start"], seg["end"]])
                bins[start_index - 1] += 1
                if (start_index != end_index):
                    bins[end_index - 1] += 1
            bins = [val - 1 if val >= 1 else val for val in bins]
            all_chr_bins += bins
    return all_chr_bins

def calc_hallmarks(dataset):
    hallmarks = [get_seg_lengths(dataset), get_changepoints(dataset), get_BP_per_chromosome(dataset), get_BP_per_10MB(dataset)]
    df = pd.DataFrame(hallmarks).T
    df.columns = ["seg_lengths", "changepoints", "bps_per_chromosome", "bps_per_10MB"]
    return df


# https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh37
hg19_chr_lengths = {
    "chr1": 249_250_621,
    "chr2": 243_199_373,
    "chr3": 198_022_430,
    "chr4": 191_154_276,
    "chr5": 180_915_260,
    "chr6": 171_115_067,
    "chr7": 159_138_663,
    "chr8": 146_364_022,
    "chr9": 141_213_431,
    "chr10": 135_534_747,
    "chr11": 135_006_516,
    "chr12": 133_851_895,
    "chr13": 115_169_878,
    "chr14": 107_349_540,
    "chr15": 102_531_392,
    "chr16": 90_354_753,
    "chr17": 81_195_210,
    "chr18": 78_077_248,
    "chr19": 59_128_983,
    "chr20": 63_025_520,
    "chr21": 48_129_895,
    "chr22": 51_304_566,
    "chrX": 155_270_560,
    "chrY": 59_373_566
}

# https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh38
hg38_chr_lengths = {
    "chr1": 248_956_422,
    "chr2": 242_193_529,
    "chr3": 198_295_559,
    "chr4": 190_214_555,
    "chr5": 181_538_259,
    "chr6": 170_805_979,
    "chr7": 159_345_973,
    "chr8": 145_138_636,
    "chr9": 138_394_717,
    "chr10": 133_797_422,
    "chr11": 135_086_622,
    "chr12": 133_275_309,
    "chr13": 114_364_328,
    "chr14": 107_043_718,
    "chr15": 101_991_189,
    "chr16": 90_338_345,
    "chr17": 83_257_441,
    "chr18": 80_373_285,
    "chr19": 58_617_616,
    "chr20": 64_444_167,
    "chr21": 46_709_983,
    "chr22": 50_818_468,
    "chrX": 156_040_895,
    "chrY": 57_227_415
}

hg19_chr_cum_starts = {
    'chr1': 0,
    'chr2': 249250621,
    'chr3': 492449994,
    'chr4': 690472424,
    'chr5': 881626700,
    'chr6': 1062541960,
    'chr7': 1233657027,
    'chr8': 1392795690,
    'chr9': 1539159712,
    'chr10': 1680373143,
    'chr11': 1815907890,
    'chr12': 1950914406,
    'chr13': 2084766301,
    'chr14': 2199936179,
    'chr15': 2307285719,
    'chr16': 2409817111,
    'chr17': 2500171864,
    'chr18': 2581367074,
    'chr19': 2659444322,
    'chr20': 2718573305,
    'chr21': 2781598825,
    'chr22': 2829728720,
    'chrX': 2881033286,
    'chrY': 3036303846
 }

hg38_chr_cum_starts = {
    'chr1': 0,
    'chr2': 248956422,
    'chr3': 491149951,
    'chr4': 689445510,
    'chr5': 879660065,
    'chr6': 1061198324,
    'chr7': 1232004303,
    'chr8': 1391350276,
    'chr9': 1536488912,
    'chr10': 1674883629,
    'chr11': 1808681051,
    'chr12': 1943767673,
    'chr13': 2077042982,
    'chr14': 2191407310,
    'chr15': 2298451028,
    'chr16': 2400442217,
    'chr17': 2490780562,
    'chr18': 2574038003,
    'chr19': 2654411288,
    'chr20': 2713028904,
    'chr21': 2777473071,
    'chr22': 2824183054,
    'chrX': 2875001522,
    'chrY': 3031042417
 }

hg19_genome_length = 3095677412
hg19_autosome_length = 2881033286
hg38_genome_length = 3088269832
hg38_autosome_length = 2875001522

# Define a function to convert a row to a list
def row_to_list(row, column, step_size):
    # Initialize an empty list
    res = []

    start = row['start'] + hg19_chr_cum_starts[row['chrom']]
    if start % step_size != 0:
        start = start - (start % step_size) + step_size
    end = row['end'] + hg19_chr_cum_starts[row['chrom']]

    # Loop over each position between start and end, divisible by step
    for pos in range(start, end, step_size):
        cn = row[column]
        res.append((pos, cn))

    # Return the list
    return res


def sample_to_SNPs(sample, column, step_size):
    rows = sample.apply(lambda x: row_to_list(x, column, step_size), axis=1)
    filtered = rows[rows.apply(lambda x: len(x) > 0)]
    return pd.Series(dict(np.concatenate(filtered.values)))


def samples_to_SNPs(cns, column='cn', step_size=1_000_000):
    # list of unique indices in cns
    samples = cns.index.unique()
    positions = {}
    for sample in samples:
        positions[sample] = sample_to_SNPs(cns.loc[sample, :], column, step_size)
    df_CNs = pd.DataFrame(positions)
    return df_CNs
