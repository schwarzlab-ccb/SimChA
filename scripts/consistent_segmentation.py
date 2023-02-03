# import all necessary packages for this file
import re
import argparse

import numpy as np
import pandas as pd
from numba import jit
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(x): return x

from plotting_functions import load_data


def consistend_segmentation(raw_data):
    '''Segment the data in a consistent way across all samples'''

    samples = raw_data.index.get_level_values('sample_id').unique()

    shared_segments = [[], [], []]  # chrom, start, end
    for chrom in raw_data.index.get_level_values('chr').unique():
        breakpoints = np.sort(np.unique(np.concatenate([raw_data.loc[(slice(None), chrom), :].eval('start'),
                                                        raw_data.loc[(slice(None), chrom), :].eval('end')])))

        shared_segments[0].extend(((len(breakpoints) - 1)) * [chrom])
        shared_segments[1].extend(breakpoints[:-1])
        shared_segments[2].extend(breakpoints[1:])

    shared_segments = [(x[0], x[1], x[2]) for x in zip(*shared_segments)]
    segmented_data = pd.DataFrame(columns=pd.MultiIndex.from_product([samples, ['cn_a', 'cn_b']], names=['sample_id', 'alleles']),
                                  index=pd.MultiIndex.from_tuples(shared_segments, names=['chr', 'start', 'end']))

    for chrom in tqdm(raw_data.index.get_level_values('chr').unique()):
        breakpoints = np.sort(np.unique(np.concatenate([raw_data.loc[(slice(None), chrom), :].eval('start'),
                                                        raw_data.loc[(slice(None), chrom), :].eval('end')])))

        all_sample_chrom_data = raw_data.loc[(
            slice(None), chrom), :].reset_index().values[:, 3:].astype(int)
        sample_index = np.array(
            raw_data.loc[(slice(None), chrom), :].index.get_level_values('sample_id'))
        res = _consistend_segmentation_step(
            breakpoints, samples, sample_index, all_sample_chrom_data)
        segmented_data.loc[chrom, breakpoints[:-1], breakpoints[1:]] = res

    segmented_data = segmented_data.stack('sample_id').reset_index()
    segmented_data['chr'] = format_chromosomes(segmented_data['chr'])
    segmented_data = (segmented_data
                      .set_index(['sample_id', 'chr', 'start', 'end'])
                      .sort_index()
                      .astype(int))

    return segmented_data


@jit(nopython=False, forceobj=True)
def _consistend_segmentation_step(breakpoints, samples, sample_index, all_sample_chrom_data):

    sample_cn_data = np.zeros((len(breakpoints) - 1, 2 * len(samples)))
    for n, sample in enumerate(samples):
        cur_segment_i = 0
        for i, (b1, b2) in enumerate(zip(breakpoints, breakpoints[1:])):
            sample_chrom_data = all_sample_chrom_data[sample_index == sample]
            if b2 > sample_chrom_data[cur_segment_i, 0]:
                cur_segment_i += 1
            sample_cn_data[i, (2*n):(2*n+2)] = sample_chrom_data[cur_segment_i, 1:]

    return sample_cn_data


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


def parse_args():
    parser = argparse.ArgumentParser(
        description='Segment the data in a consistent way across all samples')
    parser.add_argument('--i', '-input', dest='input', help='Input file')
    parser.add_argument('--o', '-output', dest='output', help='Output file')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    raw_data = load_data(args.input)
    segmented_data = consistend_segmentation(raw_data)
    segmented_data.to_csv(args.output, sep='\t')
