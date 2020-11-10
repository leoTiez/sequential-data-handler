#!/usr/bin/python3
import os

import roman
import numpy as np
import scipy.interpolate as interp
import pyBigWig
from pybedtools import BedTool


def center_norm(data):
    return (data - data.mean()) / data.std()


def remap_norm(data):
    data -= data.min()
    return data / data.max()


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


def load_bam_bed_file(name, rel_path='data'):
    cur_dir = os.getcwd()
    abs_path = cur_dir + '/' + rel_path + '/' + name
    file = BedTool(abs_path)
    return file


def load_big_file(name, rel_path='data'):
    cur_dir = os.getcwd()
    abs_path = cur_dir + '/' + rel_path + '/' + name

    file = pyBigWig.open(abs_path)
    return file


def get_values(bw_list):
    if len(bw_list) == 0:
        raise ValueError('List with bigwig objects must not be empty.')

    all_values = []
    for _ in bw_list:
        all_values.append([])

    counter = 0
    chrom_start = {}
    for chrom_num, _ in enumerate(bw_list[0].chroms().values()):
        roman_num = roman.toRoman(chrom_num+1)
        if chrom_num == 16:
            roman_num = 'M'
        length = int(bw_list[0].chroms('chr%s' % roman_num))
        chrom_start['chr%s' % roman_num] = counter
        counter += length
        for num, bw in enumerate(bw_list):
            values = np.asarray(bw.values('chr%s' % roman_num, 0, length))
            all_values[num].extend(np.nan_to_num(values, nan=0.0).tolist())

    return all_values, chrom_start


def normalise_over_annotation(bw_list, bed_ref, smoothing=None, normalise=None, trans_dict=False):
    all_values, chrom_start = get_values(bw_list=bw_list)
    bw_gen_mapping = [[] for _ in all_values]

    means = []
    stds = []
    if smoothing is None:
        smoothing = np.repeat(None, len(bw_list))
    for num, smooth in enumerate(smoothing):
        all_values[num] = np.asarray(all_values[num])
        if smooth is not None:
            all_values[num] = np.convolve(all_values[num], np.ones(smooth)/float(smooth), mode='same') \
                              + np.flip(np.convolve(np.flip(all_values[num]), np.ones(smooth)/float(smooth), mode='same')) / 2.
        means.append(all_values[num].mean())
        stds.append(all_values[num].std())
        if normalise is not None:
            if normalise.lower() == 'center':
                all_values[num] = center_norm(all_values[num])
            elif normalise.lower() == 'remap':
                all_values[num] = remap_norm(all_values[num])

    if trans_dict:
        t_dict = {}
    for int_num, interval in enumerate(bed_ref):
        # index 0: Chromosome, 1: start, 2: end, 3: name, 5: strand
        if trans_dict:
            t_v = []
        for num, (values, gen_mapping) in enumerate(zip(all_values, bw_gen_mapping)):
            start = chrom_start[interval[0]]
            frag_values = values[start + int(interval[1]): start + int(interval[2])]
            try:
                if interval[5] == '-':
                    frag_values = np.flip(frag_values)
            except IndexError:
                pass
            frag_values = np.nan_to_num(frag_values, copy=False, nan=0.)
            gen_mapping.append(frag_values)
            if trans_dict:
                t_v.append(frag_values)

        if trans_dict:
            t_dict[interval[3]] = t_v

    if not trans_dict:
        return bw_gen_mapping, means, stds, all_values, chrom_start
    else:
        return bw_gen_mapping, means, stds, all_values, chrom_start, t_dict


def data_scaling(transcript_data, vec_len=1000):
    vec_len = vec_len if vec_len is not None else len(max(transcript_data, key=len))
    data_array = np.zeros((len(transcript_data), vec_len))
    for num, td in enumerate(transcript_data):
        td = np.asarray(td)
        inter_td = interp.interp1d(np.arange(td.size), td)
        data_array[num] = inter_td(np.linspace(0, td.size - 1, vec_len))

    return data_array

