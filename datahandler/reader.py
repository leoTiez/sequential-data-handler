#!/usr/bin/python3
import os
import numpy as np
import pyBigWig
from pybedtools import BedTool
from Bio import SeqIO


def set_path(name, rel_path='data', is_abs_path=False):
    if not is_abs_path:
        cur_dir = os.getcwd()
        return cur_dir + '/' + rel_path + '/' + name
    else:
        return name


def load_bam_bed_file(name, rel_path='data', is_abs_path=False):
    path = set_path(name, rel_path=rel_path, is_abs_path=is_abs_path)
    file = BedTool(path)
    return file


def load_big_file(name, rel_path='data', is_abs_path=False):
    path = set_path(name, rel_path=rel_path, is_abs_path=is_abs_path)
    file = pyBigWig.open(path)
    return file


def load_fast(name, rel_path='data', is_abs_path=False, is_fastq=True):
    path = set_path(name, rel_path=rel_path, is_abs_path=is_abs_path)
    if is_fastq:
        return list(SeqIO.parse(path, 'fastq'))
    else:
        return list(SeqIO.parse(path, 'fasta'))


def create_bed_random_fragments(bw, max_chunk=6000, name='random_fragments', path='/'):
    curr_dir = os.getcwd()
    bed = open('%s/%s/%s.bed' % (curr_dir, path, name), 'w+')
    for chrom, length in bw.chroms().items():
        count = 0
        while True:
            next_stop = count + int(np.random.random() * max_chunk)
            if next_stop >= length:
                next_stop = length - 1

            bed.write('%s\t%s\t%s\n' % (chrom, count, next_stop))
            if next_stop >= length - 1:
                break
            count = next_stop

    bed.close()

