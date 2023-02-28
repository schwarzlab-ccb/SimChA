import argparse

import numpy as np
import pandas as pd

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(x): return x

from plotting_functions import load_data
try:
    from segmentation_cython.consistent_segmentation_step import \
        _consistent_segmentation_step_cython
except ImportError:
    raise ImportError('Please run `python setup.py build_ext --inplace` in the '
                      'simcha/scripts/segmentation_cython directory with cython installed')


def consistent_segmentation(raw_data, alleles=['cn_a', 'cn_b'], chrom_col='chr'):
    '''Segment the data in a consistent way across all samples'''

    assert len(np.setdiff1d(raw_data.columns, ['sample_id', 'chr', 'start', 'end'])) == len(alleles)

    assert (raw_data.groupby(['sample_id', chrom_col], observed=True)['start'].min() == 0).all()
    assert (raw_data.groupby(['sample_id', chrom_col], observed=True)['end'].max().groupby(chrom_col, observed=True).std() == 0).all()

    samples = raw_data['sample_id'].unique()

    shared_segments = [[], [], []]  # chrom, start, end
    segmented_data = []
    for chrom in tqdm(raw_data[chrom_col].unique()):
        raw_data_chrom = raw_data.query('chr == @chrom')
        breakpoints = np.sort(np.unique(np.concatenate([raw_data_chrom['start'],
                                                        raw_data_chrom['end']])))

        shared_segments[0].extend(((len(breakpoints) - 1)) * [chrom])
        shared_segments[1].extend(breakpoints[:-1])
        shared_segments[2].extend(breakpoints[1:])

        all_sample_chrom_data = raw_data_chrom.values[:, 3:].astype(int)
        sample_index = raw_data_chrom['sample_id'].values
        
        label, count = np.unique(sample_index, return_counts=True)
        sample_index_count = count[np.argsort(label)]

        res = _consistent_segmentation_step_cython(
                    breakpoints[1:], sample_index_count,
                    all_sample_chrom_data, len(alleles))
        assert (res != -1).all()

        cur_segmented_data = pd.DataFrame(np.concatenate(res.swapaxes(0, 2).swapaxes(1, 2), axis=1),
                                          columns=pd.MultiIndex.from_product([alleles, samples], names=['alleles', 'sample_id']))
        cur_segmented_data[chrom_col] = chrom
        cur_segmented_data['start'] = breakpoints[:-1]
        cur_segmented_data['end'] = breakpoints[1:]
        cur_segmented_data = (cur_segmented_data
                              .set_index([chrom_col, 'start', 'end'])
                              .stack('sample_id')
                              .reorder_levels(['sample_id', chrom_col, 'start', 'end']))
        segmented_data.append(cur_segmented_data)

    segmented_data = pd.concat(segmented_data, axis=0)
    segmented_data = (segmented_data
                      .reorder_levels(['sample_id', chrom_col, 'start', 'end'])
                      .sort_index()) # this line is the most computationally expensive

    assert raw_data.eval(f'(end - start)').sum() == segmented_data.eval(f'(end - start)').sum()
    for allele in alleles:
        assert raw_data.eval(f'(end - start) * {allele}').sum() == segmented_data.eval(f'(end - start) * {allele}').sum()

    return segmented_data


def parse_args():
    parser = argparse.ArgumentParser(
        description='Segment the data in a consistent way across all samples')
    parser.add_argument('--i', '-input', dest='input', help='Input file')
    parser.add_argument('--o', '-output', dest='output', help='Output file')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    raw_data = load_data(args.input).reset_index()
    segmented_data = consistent_segmentation(raw_data)
    segmented_data.to_csv(args.output, sep='\t')

