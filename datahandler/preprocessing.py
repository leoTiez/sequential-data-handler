#!/usr/bin/python3
import re
import numpy as np
from scipy.signal import argrelmax


def peak_detect_smooth(sig, peak_range=200, mode='wrap'):
    peaks, = argrelmax(sig, order=peak_range, mode=mode)
    return peaks


def cancel_noise_cpd(cpd_sig, chrom_start, dna_seq):
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


