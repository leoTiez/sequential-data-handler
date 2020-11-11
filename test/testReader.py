#!/usr/bin/python3
import unittest
import os
import src.reader as reader
import wget


class TestReader(unittest.TestCase):
    def test_set_path(self):
        name = 'test_name'
        rel_path = 'data/data'
        curr_dir = os.getcwd()

        path = reader.set_path(name, rel_path=rel_path, is_abs_path=True)
        self.assertEqual(path, name)

        path = reader.set_path(name, rel_path=rel_path, is_abs_path=False)
        self.assertEqual(path, curr_dir + '/' + rel_path + '/' + name)

    def test_load_bam_bed_file(self):
        name = 'test_bed.bed'
        rel_path = 'data'
        count = 4077
        bed = reader.load_bam_bed_file(name=name, rel_path=rel_path, is_abs_path=False)

        self.assertIsNot(bed, None)
        self.assertEqual(bed.count(), count)

    def test_load_big_file(self):
        test_chrom = 'chrX'
        test_start = 0
        test_end = 1
        exp_result = 3.55
        if not os.path.isfile('data/bigWigExample.bw'):
            url = 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM1499626&format=file&file=GSM1499626%5Finput%' \
                  '5F2r%5Fdelta%5Fbromo%5Funstressed%5Frep2%5Fnterm%2Ebw'
            rel_path = 'data/'
            name = wget.download(url, out=rel_path + 'bigWigExample.bw')
        else:
            name = 'data/bigWigExample.bw'

        bigwig = reader.load_big_file(name, rel_path='', is_abs_path=False)
        self.assertIsNot(bigwig, None)
        self.assertTrue(bigwig.isBigWig())
        self.assertAlmostEqual(bigwig.values(test_chrom, test_start, test_end)[0], exp_result, 2)

    def test_load_fast(self):
        if not os.path.isfile('data/NC_005816.fna'):
            url = 'https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.fna'
            rel_path = 'data'
            name = wget.download(url, out=rel_path)
        else:
            name = 'data/NC_005816.fna'
        fasta = reader.load_fast(name, rel_path='', is_abs_path=False, is_fastq=False)
        self.assertIsNot(fasta, None)


if __name__ == '__main__':
    unittest.main()
