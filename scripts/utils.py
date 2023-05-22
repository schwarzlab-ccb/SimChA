import re

import numpy as np
import pandas as pd
import json

from os.path import join

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
    loaded_data["major"] = get_CNs_by_sample(loaded_data["CNs"], True)
    loaded_data["minor"] = get_CNs_by_sample(loaded_data["CNs"], False) 
    loaded_data["hap_a"] = get_hap_by_sample(loaded_data["CNs"], "cn_a")
    loaded_data["hap_b"] = get_hap_by_sample(loaded_data["CNs"], "cn_b")
    config_file = join(dataset_path, "sim_params.json")
    # import json into a dict
    with open(config_file) as f:
        loaded_data["config"] = json.load(f)
    return loaded_data