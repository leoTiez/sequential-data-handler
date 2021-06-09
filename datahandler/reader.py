#!/usr/bin/python3
"""Reader module

Collection of functions for loading files and data values that are related to sequencing data and genomic analysis.
The provided functions are
* set_path - Create path to file
* load_bam_bed_file - Load bam or bed file
* load_gff - Load gff annotation file
* load_big_file - load bigwig file
* load_fast - Load fasta or fastq file
* create_bed_random_fragments - Create random fragments and save them in a bed file
"""
import os
import numpy as np
import pyBigWig
from pybedtools import BedTool
from Bio import SeqIO


def set_path(name, rel_path='data', is_abs_path=False):
    """
    Create path to file
    :param name: Name of the file or absolute path if is_abs_path is set to True
    :type name: str
    :param rel_path: Relative path without the name from current directory
    :type rel_path: str
    :param is_abs_path: If True, name is interpreted as absolute path.
    :type is_abs_path: bool
    :return: Path to file as string
    """
    if not is_abs_path:
        cur_dir = os.getcwd()
        return cur_dir + '/' + rel_path + '/' + name
    else:
        return name


def load_bam_bed_file(name, rel_path='data', is_abs_path=False):
    """
    Load bam or bed file
    :param name: Name of the file or absolute path if is_abs_path is set to True
    :type name: str
    :param rel_path: Relative path without the name from current directory
    :type rel_path: str
    :param is_abs_path: If True, name is interpreted as absolute path.
    :type is_abs_path: bool
    :return: BedTool file object
    """
    path = set_path(name, rel_path=rel_path, is_abs_path=is_abs_path)
    file = BedTool(path)
    return file


def load_gff(name, rel_path='data', is_abs_path=False):
    """
    Load gff annotation file
    :param name: Name of the file or absolute path if is_abs_path is set to True
    :type name: str
    :param rel_path: Relative path without the name from current directory
    :type rel_path: str
    :param is_abs_path: If True, name is interpreted as absolute path.
    :type is_abs_path: bool
    :return: Pipe to file
    """
    path = set_path(name, rel_path=rel_path, is_abs_path=is_abs_path)
    file = open(path)
    return file


def load_big_file(name, rel_path='data', is_abs_path=False):
    """
    Load bigwig file
    :param name: Name of the file or absolute path if is_abs_path is set to True
    :type name: str
    :param rel_path: Relative path without the name from current directory
    :type rel_path: str
    :param is_abs_path: If True, name is interpreted as absolute path.
    :type is_abs_path: bool
    :return: bigWigFile object
    """
    path = set_path(name, rel_path=rel_path, is_abs_path=is_abs_path)
    file = pyBigWig.open(path)
    return file


def load_fast(name, rel_path='data', is_abs_path=False, is_fastq=True):
    """
    Load fasta or fastq file
    :param name: Name of the file or absolute path if is_abs_path is set to True
    :type name: str
    :param rel_path: Relative path without the name from current directory
    :type rel_path: str
    :param is_abs_path: If True, name is interpreted as absolute path.
    :type is_abs_path: bool
    :param is_fastq: If true, the loaded file is interpreted as fastq
    :type is_fastq: bool
    :return: Parsed fasta/fastq file converted to a list to make it reusable
    """
    path = set_path(name, rel_path=rel_path, is_abs_path=is_abs_path)
    if is_fastq:
        return list(SeqIO.parse(path, 'fastq'))
    else:
        return list(SeqIO.parse(path, 'fasta'))


def create_bed_random_fragments(chrom_dict, max_chunk=6000, name='random_fragments', path='/'):
    """
    Create random fragments and save them in a bed file
    :param chrom_dict: Dictionary with chromosome name as key and chromosome size as value
    :type chrom_dict: dict
    :param max_chunk: Maximal size of a fragment
    :type max_chunk: int
    :param name: Name of the bed file that is to be created
    :type name: str
    :param path: Path where the bed file is to be saved
    :type path: str
    :return: None
    """
    curr_dir = os.getcwd()
    bed = open('%s/%s/%s.bed' % (curr_dir, path, name), 'w+')
    for chrom, length in chrom_dict.items():
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

