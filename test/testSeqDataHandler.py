#!/usr/bin/python3
import unittest
import os
import numpy as np
import wget

from datahandler import reader

from datahandler import seqDataHandler as seq


class TestSeqDataHandler(unittest.TestCase):
    def test_center_norm_and_center_norm_all(self):
        data = np.arange(9)
        mean = 4.
        std = 2.58
        expect_data = (data - mean) / std
        norm_data = seq.center_norm(data)

        for norm, exp in zip(norm_data, expect_data):
            self.assertAlmostEqual(exp, norm, 2)

        data_2 = np.arange(2, 20)
        mean_2 = 10.5
        std_2 = 5.18
        expect_data_2 = (data_2 - mean_2) / std_2
        norm_data_all = seq.center_norm_all([data, data_2])

        for nd, ed in zip(norm_data_all, [expect_data, expect_data_2]):
            for norm, exp in zip(nd, ed):
                self.assertAlmostEqual(exp, norm, 2)

    def test_remap_norm_and_remap_norm_all(self):
        data = np.arange(9.)
        expect_data = data / 8.
        norm_data = seq.remap_norm(data)

        for norm, exp in zip(norm_data, expect_data):
            self.assertAlmostEqual(exp, norm, 2)

        data_2 = np.arange(2., 20.)
        expect_data_2 = data_2 - 2
        expect_data_2 /= 17.
        norm_data_all = seq.remap_norm_all([data, data_2])

        for nd, ed in zip(norm_data_all, [expect_data, expect_data_2]):
            for norm, exp in zip(nd, ed):
                self.assertAlmostEqual(exp, norm, 2)

    def test_smooth_and_smooth_all(self):
        test_data = np.arange(9.)
        smooth_size = 3
        exp_result = [0.33, 1., 2., 3., 4., 5., 6., 7., 5.]
        smooth_data = seq.smooth(test_data, smooth_size=smooth_size)

        for smooth, exp in zip(smooth_data, exp_result):
            self.assertAlmostEqual(exp, smooth, 2)

        data_2 = np.arange(2., 20.)
        smooth_size_2 = 4
        exp_result_2 = [1.25, 2.25, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5,
                        11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 13.5]
        smooth_data_all = seq.smooth_all([test_data, data_2], smooth_list=[smooth_size, smooth_size_2])
        for sd, ed in zip(smooth_data_all, [exp_result, exp_result_2]):
            for smooth, exp in zip(sd, ed):
                self.assertAlmostEqual(exp, smooth, 2)

    def test_annotate_and_annotate_all(self):
        name_bed = 'test_bed.bed'
        rel_path = 'data'
        anno_idx = 50
        test_start = 115941
        test_end = 119543

        if not os.path.isfile('data/bigWigExample.bw'):
            url = 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM1499626&format=file&file=GSM1499626%5Finput%' \
                  '5F2r%5Fdelta%5Fbromo%5Funstressed%5Frep2%5Fnterm%2Ebw'
            rel_path = 'data/'
            name_big = wget.download(url, out=rel_path + 'bigWigExample.bw')
        else:
            name_big = 'data/bigWigExample.bw'

        bed = reader.load_bam_bed_file(name_bed, rel_path=rel_path, is_abs_path=False)
        bigwig = reader.load_big_file(name_big, rel_path='', is_abs_path=False)
        all_values, chrom_dict = seq.get_values([bigwig])

        anno, _ = seq.annotate(all_values[0], bed, chrom_start=chrom_dict)
        self.assertListEqual(anno[anno_idx].tolist(), all_values[0][test_start:test_end].tolist())

        anno_l, _ = seq.annotate_all(all_values, bed, chrom_start=chrom_dict)
        self.assertListEqual(anno_l[0][anno_idx].tolist(), all_values[0][test_start:test_end].tolist())

    def test_rescaling_and_rescaling_all(self):
        data = np.arange(5)
        data_2 = np.arange(15)
        vec_size = 10
        exp_result = [0., 0.444, 0.888, 1.333, 1.777, 2.222, 2.666, 3.111, 3.555, 4.]
        exp_result_2 = [0., 1.555, 3.111, 4.666, 6.222, 7.777, 9.333, 10.888, 12.444, 14.]

        rescale = seq.rescale(data, vec_len=vec_size)
        for scale, exp in zip(rescale, exp_result):
            self.assertAlmostEqual(exp, scale, 2)

        rescale_all = seq.rescale_all([data, data_2], vec_len=vec_size)
        for rd, ed in zip(rescale_all, [exp_result, exp_result_2]):
            for scale, exp in zip(rd, ed):
                self.assertAlmostEqual(exp, scale, 2)

    def test_get_values(self):
        test_start = 230208
        test_end = 230210
        exp_result = [2.87, 2.87]
        test_chrom = 'chrIV'
        test_chrom_start = 230208 + 813178 + 316617
        if not os.path.isfile('data/bigWigExample.bw'):
            url = 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM1499626&format=file&file=GSM1499626%5Finput%' \
                  '5F2r%5Fdelta%5Fbromo%5Funstressed%5Frep2%5Fnterm%2Ebw'
            rel_path = 'data/'
            name = wget.download(url, out=rel_path + 'bigWigExample.bw')
        else:
            name = 'data/bigWigExample.bw'
        bigwig = reader.load_big_file(name, rel_path='', is_abs_path=False)

        all_values, chrom_dict = seq.get_values([bigwig])

        result = all_values[test_start:test_end]
        for exp, res in zip(exp_result, result):
            self.assertAlmostEqual(exp, res, 2)

        self.assertEqual(chrom_dict[test_chrom], test_chrom_start)


if __name__ == '__main__':
    unittest.main()
