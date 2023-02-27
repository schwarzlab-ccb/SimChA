cimport numpy as np
import numpy as np
import cython

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
def _consistent_segmentation_step_cython_single_allele(
    np.ndarray[np.long_t, ndim=1] ends,
    np.ndarray[np.long_t, ndim=1] sample_index_count, 
    np.ndarray[np.long_t, ndim=2] all_sample_chrom_data
    ):

    cdef int offset = 0
    cdef np.long_t cur_end = 0
    cdef Py_ssize_t n, i, cur_segment_i
    cdef Py_ssize_t n_max = sample_index_count.shape[0]
    cdef Py_ssize_t i_max = ends.shape[0]

    cdef np.ndarray[np.long_t, ndim=2] sample_cn_data = -1 * np.ones([i_max, n_max], dtype=np.int64)

    for n in range(n_max):
        cur_segment_i = 0
        for i in range(i_max):
            cur_end = ends[i]
            if cur_end > all_sample_chrom_data[offset + cur_segment_i, 0]:
                cur_segment_i += 1
            sample_cn_data[i, n] = all_sample_chrom_data[offset + cur_segment_i, 1]
        offset += sample_index_count[n]

    return sample_cn_data

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
def _consistent_segmentation_step_cython(
    np.ndarray[np.long_t, ndim=1] ends,
    np.ndarray[np.long_t, ndim=1] sample_index_count, 
    np.ndarray[np.long_t, ndim=2] all_sample_chrom_data,
    np.int32_t n_alleles
    ):

    cdef int offset = 0
    cdef np.long_t cur_end = 0
    cdef Py_ssize_t n, i, cur_segment_i, cur_allele
    cdef Py_ssize_t n_max = sample_index_count.shape[0]
    cdef Py_ssize_t i_max = ends.shape[0]

    cdef np.ndarray[np.long_t, ndim=3] sample_cn_data = -1 * np.ones([i_max, n_max, n_alleles], dtype=np.int64)

    for n in range(n_max):
        cur_segment_i = 0
        for i in range(i_max):
            cur_end = ends[i]
            if cur_end > all_sample_chrom_data[offset + cur_segment_i, 0]:
                cur_segment_i += 1
            for cur_allele in range(n_alleles):
                sample_cn_data[i, n, cur_allele] = all_sample_chrom_data[offset + cur_segment_i, 1 + cur_allele]
        offset += sample_index_count[n]

    return sample_cn_data
