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
        "CNs" : pd.read_csv(join(dataset_path, "copynumbers.tsv"), index_col=0, sep="\t"),
        "clones" : pd.read_csv(join(dataset_path, "clones.tsv"), index_col=0, sep="\t"),
        "karyotypes" : pd.read_csv(join(dataset_path, "karyotypes.tsv"), index_col=0, sep="\t"),
        "samples" : pd.read_csv(join(dataset_path, "samples.tsv"), index_col=0, sep="\t"),
    }
    config_file = join(dataset_path, "sim_params.json")
    # import json into a dict
    with open(config_file) as f:
        loaded_data["config"] = json.load(f)
    return loaded_data

def calc_CNs(dataset):
    cns = [get_CNs_by_sample(dataset, True), get_CNs_by_sample(dataset, False), get_hap_by_sample(dataset, "cn_a"), get_hap_by_sample(dataset, "cn_b")]
    df =  pd.DataFrame(cns).T
    df.columns = ["major", "minor", "hap_a", "hap_b"]
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