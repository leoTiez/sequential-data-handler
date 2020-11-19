#!/usr/bin/python3
import unittest
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from datahandler import reader, preprocessing, seqDataHandler


class TestPreprocessing(unittest.TestCase):
    def test_peak_detect_smooth(self):
        test_seq = np.asarray([0, 1, 0, 2, 3, 3, 3, 3, 4, 3, 2, 1, 0, 1, 1, 2, 0])
        peak_range = 3
        exp_peaks = [8, 15]

        peaks = preprocessing.peak_detect_smooth(test_seq, peak_range=peak_range)
        self.assertListEqual(peaks.tolist(), exp_peaks)

    def test_cancel_noise_cpd(self):
        genome = [SeqRecord(Seq('ACGACGTTA'))]
        cpd = np.arange(1, 10)
        exp_filtered = np.arange(1, 10)
        exp_filtered_idx = np.asarray([2, 3, 6, 7])
        mask = np.zeros(9)
        mask[exp_filtered_idx] = 1.
        exp_filtered[~mask.astype('bool')] = 0.
        chrom_start = {'chrI': 0}

        filtered_cpd = preprocessing.cancel_noise_cpd(cpd, chrom_start, genome)
        self.assertListEqual(exp_filtered.tolist(), filtered_cpd.tolist())


if __name__ == '__main__':
    unittest.main()
