import argparse

import numpy as np
import pandas as pd
from numba import jit
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(x): return x

from plotting_functions import load_data
from utils import format_chromosomes


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
