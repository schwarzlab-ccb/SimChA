import re

import numpy as np
import pandas as pd
import json
from collections import defaultdict


from os.path import join

chromosome_names = [
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
    'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
    'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21',
    'chr22', 'chrX', 'chrY'
]

autosome_names = chromosome_names[:-2]

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

# Get the segment lengths for a dataset (only non-diploid segments by default)
def get_seg_lengths(data, include_cn_normal = False, include_loh = False):
    sample_ids = pd.Series({c: data[c].unique() for c in data})["sample_id"]
    all_seg_lengths = []
    for _, id in enumerate(sample_ids):
        sample = data[data["sample_id"] == id]
        for _, row in sample.iterrows():
            if row["cn_a"] == 1 and row["cn_b"] == 1 and not include_cn_normal:
                continue
            elif row["cn_a"] != 1 and row["cn_a"]+row["cn_b"] == 2 and not include_loh:
                continue
            all_seg_lengths.append(row["end"] - row["start"])
    return all_seg_lengths

# Get the changepoint values for a dataset
def get_changepoints(data, include_cn_normal = False, include_loh = False):
    sample_ids = pd.Series({c: data[c].unique() for c in data})["sample_id"]
    all_changepoints = []
    data["cn"] = data["cn_a"] + data["cn_b"]
    for _, id in enumerate(sample_ids):
        sample = data[data["sample_id"] == id]
        for chr in chromosome_names:
            chr_changepoints = []
            segs = sample[sample["chrom"] == chr]
            last_seg = 2
            for _, seg in segs.iterrows():
                # By default we don't count diploid segments
                if seg["cn"] == 2:
                    # Update this segment
                    last_seg = 2
                    # Ordinary diploid segment
                    if seg["cn_a"] == 1 and seg["cn_b"] == 1 and not include_cn_normal:
                        continue
                    # LOH segments, only need to check if the the copy number of a is not 1
                    elif seg["cn_a"] != 1 and not include_loh:
                        continue
                this_seg = seg["cn"]
                chr_changepoints.append(abs(this_seg - last_seg))
                last_seg = this_seg
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

def get_BP_per_bin_size(data, bin_size=10_000_000):
    sample_ids = pd.Series({c: data[c].unique() for c in data})["sample_id"]
    all_chr_bins = []
    # Loop over the samples
    for i, id in enumerate(sample_ids):
        sample = data[data["sample_id"] == id]
        # Loop over the chromosomes
        for _, chr in enumerate(chromosome_names): 
            segs = sample[sample["chrom"] == chr]
            intervals = np.arange(0, hg19_chr_lengths[chr]+bin_size, bin_size)
            # Just want the counts of all the end points
            res = np.histogram(segs['end'].iloc[:-1], bins=intervals)[0]
            all_chr_bins.extend(res)
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

hg19_genome_xy_length = 3095677412
hg19_genome_xx_length = 3036303846
hg19_autosome_length = 2881033286
hg19_genome_length = hg19_genome_xy_length

hg38_genome_xy_length = 3088269832
hg38_genome_xx_length = 3031042417
hg38_autosome_length = 2875001522
hg38_genome_length = hg38_genome_xy_length

# Define a function to convert a row to a list
def row_to_list(row, bins, column):
    # Initialize an empty list
    res = []
    start_index = np.digitize(row['start']+hg19_chr_cum_starts[row['chrom']], bins)
    start = bins[start_index]
    #start = row['start'] + hg19_chr_cum_starts[row['chrom']]
    # Round start to the nearest bin
    #if start % step_size != 0:
    #    start = start - (start % step_size) + step_size
    end_index = np.digitize(row['end']+hg19_chr_cum_starts[row['chrom']]-1, bins)
    end = bins[end_index]
    # The rounding can cause the start to now be greater than the end
    # in which case, it is a single bin
    if start_index == end_index:
        return [(start, row[column])]

    # Loop over each position between start and end, divisible by step
    for index in range(start_index, end_index+1):
        cn = row[column]
        res.append((bins[index], cn))
    # Return the list
    return res

def chromosome_to_bins(chrom, step_size=1_000_000):
    bins = np.arange(step_size, hg19_chr_lengths[chrom], step_size)
    bins += hg19_chr_cum_starts[chrom] - 1
    bins = np.append(bins, hg19_chr_lengths[chrom] + hg19_chr_cum_starts[chrom])
    return bins

def get_chromosome_bins(step_size=1_000_000, includeSexChromosomes=False):
    bins = np.array([])
    chroms = autosome_names if not includeSexChromosomes else chromosome_names
    for chrom in chroms:
        chr_bins = chromosome_to_bins(chrom, step_size=step_size)
        bins = np.append(bins, chr_bins)
    return bins

def sample_to_SNPs(sample, bins, column):
    count = {b : 0 for _, b in enumerate(bins)}
    rows = sample.apply(lambda x: row_to_list(x, bins, column), axis=1)
    filtered = rows[rows.apply(lambda x: len(x) > 0)]
    for items in filtered.values:
        for pos, cn in items:
            count[pos] += cn
    return pd.Series(count)
    #return pd.Series(np.concatenate(filtered.values))

def samples_to_SNPs(cns, column='cn', step_size=1_000_000, includeSexChromosomes=False):
    # list of unique indices in cns
    samples = cns.index.unique()
    positions = {}
    bins = get_chromosome_bins(step_size, includeSexChromosomes)
    for sample in samples:
        positions[sample] = sample_to_SNPs(cns.loc[sample, :], bins, column)
    df_CNs = pd.DataFrame(positions)
    return df_CNs

def sample_to_homozygous_spots(sample, bins):
    count = {b : 0 for _, b in enumerate(bins)}
    rows = sample.apply(lambda x: row_to_list(x, bins, 'homozygous_del'), axis=1)
    filtered = rows[rows.apply(lambda x: len(x) > 0)]
    for items in filtered.values:
        for pos, cn in items:
            count[pos] += cn
    return pd.Series(count)

def get_homozygous_deletion_locations(df, step_size=1_000_000):
    df['cn'] = df['cn_a'] + df['cn_b']
    # 1 = homozygously deleted region, 0 otherwise
    df['homozygous_del'] = (df['cn'] == 0).astype(int)
    # Ignore sex chromosomes
    df = df[(df['chrom'] != 'chrX') & (df['chrom'] != 'chrY')]
    #df.reset_index(inplace=True)
    #df.set_index(['sample_id', 'chrom'], inplace=True)
    samples = df.index.unique()
    positions = {}
    bins = get_chromosome_bins(step_size, includeSexChromosomes=False)
    # Loop over the samples
    for sample in samples:
        # For each chromosome, count the number of homozygous deletions
        positions[sample] = sample_to_homozygous_spots(df.loc[sample,:], bins)
        #positions[sample] = sample_to_homozygous_spots(df.loc[sample, :], bins, step_size)
    df_dels = pd.DataFrame(positions)
    return df_dels

def homozygous_length_distribution(df, bin_size=1_000_000):
    df["cn"] = df["cn_a"] + df["cn_b"]
    sample_ids = df["sample_id"].unique()
    hist = []
    for id in sample_ids:
        sample = df[df["sample_id"] == id]
        for _, row in sample.iterrows():
            # Skip segments that are not homozygously deleted or are sex chromosomes
            if row["cn"] != 0 or row['chrom'] in ["chrX", "chrY"]:
                continue
            if row["cn_a"] == 0 and row["cn_b"] == 0:
                length = (row["end"] - row["start"])/bin_size
                hist.append(length)
    return hist

