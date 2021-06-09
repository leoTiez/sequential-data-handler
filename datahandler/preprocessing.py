#!/usr/bin/python3
"""Preprocessing module

Functions for more specific and data manipulation and analysis. They can be used in the preprocessing pipeline
before a sophisticated in depth analysis is carried out. Available functions are:
* peak_detect_smooth -  Find the relative peak values for the signal
* cancel_noise_cpd - Filter CPD signal for values where CPDs are possible
"""
import re
import numpy as np
from scipy.signal import argrelmax


def peak_detect_smooth(sig, peak_range=200, mode='wrap'):
    """
    Find the relative peak values for the signal. It uses a moving window for determining relative maxima.
    :param sig: Data values
    :type sig: numpy.array
    :param peak_range: Size of the the moving window. This means how many values to the left and the right are
    considered for determining the relative maximum
    :type peak_range: int
    :param mode: Behaviour at the edges. Possible are 'wrap' or 'clip.
    See scipy.signal.argrelmax documentation for more information
    :type mode: str
    :return: numpy.array with indices for the peak values
    """
    peaks, = argrelmax(sig, order=peak_range, mode=mode)
    return peaks


def cancel_noise_cpd(cpd_sig, chrom_start, dna_seq):
    """
    Set CPD signal to zero where there is no adjacent pyrimidines
    :param cpd_sig: Data array with the cpd signal
    :type cpd_sig: numpy.array
    :param chrom_start: Dictionary with chromosome names as indices and the indices where the chromosomes start in
    the cpd_sig array
    :type chrom_start: dict
    :param dna_seq: The DNA sequence as parsed fasta or fastq
    :type dna_seq: iterable
    :return: Filtered CPD signal
    """
    combinations = ['TT', 'CT', 'TC', 'CC', 'AA', 'GA', 'AG', 'GG']
    for start, seq in zip(chrom_start.values(), dna_seq):
        seq = seq.seq
        ind = set()
        for c in combinations:
            c_ind = set([m.start() for m in re.finditer(c, str(seq))])
            ind = ind.union(c_ind)

        sig_mask = np.zeros(len(seq)).astype('bool')
        sig_mask[np.asarray(list(ind))] = True
        sig_mask[np.asarray(list(ind)) + 1] = True
        cpd_sig[start:start + len(seq)][~sig_mask] = 0

    return cpd_sig


