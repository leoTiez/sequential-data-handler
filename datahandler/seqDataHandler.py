#!/usr/bin/python3
"""Manipulation of sequential data

The seqDataHandler contains some few important functions to retrieve data from bigwig files and to
manipulate them. This includes the following functions

* center_norm - Normalise data by subtracting the mean and dividing the standard deviation
* remap_norm - Normalise data by remapping all data values between 0 (minimum) and 1 (maximum)
* smooth - Smooth data with a moving-average window
* annotate_gff_from_bw - Segment bigwig file according to gff annotation file
* annotate - Segment data according to bed annotation file
* rescale - Rescale data array
* center_norm_all - Wrapper function to center all data arrays in list
* remap_norm_all - Wrapper function to remap all data arrays in list
* smooth_all - Wrapper function to smooth all data arrays in list
* annotate_all - Wrapper function to segment all data arrays according to the passed bed file
* rescale_all - Wrapper function to rescale all data arrays in a list
* get_values  - Retrieve all data values from bigwig file and concatenate them together
"""
import warnings
import numpy as np
import scipy.interpolate as interp
from BCBio import GFF
from datahandler.reader import load_gff


def center_norm(data):
    """
    Normalise data by subtracting the mean and dividing the standard deviation
    :param data: Data array
    :type data: numpy.array
    :return: Normalised data array
    """
    return (data - data.mean()) / data.std()


def remap_norm(data):
    """
    Normalise data by remapping all data values between 0 (minimum) and 1 (maximum)
    :param data: Data array
    :type data: numpy array
    :return: Normalised data array
    """
    data -= data.min()
    return data / data.max()


def smooth(data, smooth_size=20):
    """
    Smooth data with a moving-average window
    :param data: Data array
    :type data: numpy.array
    :param smooth_size: Size of the moving average window
    :type smooth_size: int
    :return: Smoothed data
    """
    return np.convolve(data, np.ones(smooth_size) / float(smooth_size), mode='same')


def annotate_gff_from_bw(bw, gff_path, gff_source_type=[('ensembl_havana', 'gene')]):
    """
    Segment bigwig file according to gff annotation file
    :param bw: Bigwig file
    :type bw: bigWigFile
    :param gff_path: path to gff file
    :type gff_path: str
    :param gff_source_type: parameters for which the gff file is filtered
    :type gff_source_type: list(tuple(str))
    :return: Segmented data array
    """
    gff = load_gff(gff_path, rel_path='', is_abs_path=True)
    examiner = GFF.GFFExaminer()
    chrom_list = list(examiner.available_limits(gff)['gff_id'].keys())
    gff.close()

    gen_mapping = []
    gff = load_gff(gff_path, rel_path='', is_abs_path=True)
    for chrom in chrom_list:
        limit_info = dict(gff_id=chrom, gff_source_type=gff_source_type)
        for rec in GFF.parse(gff, limit_info=limit_info):
            for num, r in enumerate(rec.features):
                anno = bw.values('chr%s' % chrom[0], int(r.location.start), int(r.location.end))
                if int(r.location.strand) == -1:
                    anno = np.flip(anno)
                anno = np.nan_to_num(anno, copy=False, nan=0.)
                gen_mapping.append(anno)

    gff.close()
    return gen_mapping


def annotate(data, bed_ref, chrom_start):
    """
    Segment data according to bed annotation file
    :param data: Data array (it is assumed that all chromosomes are concatenated together)
    :type data: numpy.array
    :param bed_ref: reference to bed file
    :type bed_ref: str
    :param chrom_start: Dictionary with index positions where the chromosomes start
    :type chrom_start: dict
    :return: Segmented data array, dictionary which maps from gene name to data array
    """
    gen_mapping = []
    trans_dict = {}
    make_trans_dict = True
    if len(list(bed_ref[0])) < 4:
        warnings.warn('No annotation names found. Return empty dict for trans_dict', RuntimeWarning)
        make_trans_dict = False

    for int_num, interval in enumerate(bed_ref):
        # index 0: Chromosome, 1: start, 2: end, 3: name, 5: strand
        start = chrom_start[interval[0]]
        frag_values = data[start + int(interval[1]): start + int(interval[2])]
        try:
            if interval[5] == '-':
                frag_values = np.flip(frag_values)
        except IndexError:
            pass
        frag_values = np.nan_to_num(frag_values, copy=False, nan=0.)
        gen_mapping.append(frag_values)
        if make_trans_dict:
            trans_dict[interval[3]] = frag_values

    return gen_mapping, trans_dict


def rescale(data, vec_len=1000):
    """
    Rescale data array
    :param data: Data array
    :type data: numpy.array
    :param vec_len: Vector length to which the data is supposed to be rescaled
    :type vec_len: int
    :return: rescaled data array
    """
    data = np.asarray(data)
    inter_td = interp.interp1d(np.arange(data.size), data)
    return inter_td(np.linspace(0, data.size - 1, vec_len))


def center_norm_all(all_values):
    """
    Wrapper function to center all data arrays in list
    :param all_values: List with data arrays
    :type all_values: list(numpy.array)
    :return: List with normalised data arrays
    """
    return [center_norm(data) for data in all_values]


def remap_norm_all(all_values):
    """
    Wrapper function to remap all data arrays in list
    :param all_values: List with data arrays
    :type all_values: list(numpy.array)
    :return: List with normalised data arrays
    """
    return [remap_norm(data) for data in all_values]


def smooth_all(all_values, smooth_list):
    """
    Wrapper function to smooth all data arrays in list
    :param all_values: All data arrays in a list
    :type all_values: list(data.array)
    :param smooth_list: List or integer with sizes of the moving-average windows
    :type smooth_list: list(int) or int
    :return: List with smoothed data arrays
    """
    if type(smooth_list) == list:
        return [smooth(data, smooth_size=smooth_size) for data, smooth_size in zip(all_values, smooth_list)]
    elif type(smooth_list) == int:
        return [smooth(data, smooth_size=smooth_list) for data in all_values]
    else:
        raise ValueError('smooth_list must be either list or int')


def annotate_all(all_values, bed_ref, chrom_start):
    """
    Wrapper function to segment all data arrays according to the passed bed file
    :param all_values: List with data array (it is assumed that all chromosomes are concatenated together per data
    array)
    :type all_values: list(numpy.array)
    :param bed_ref: reference to bed file
    :type bed_ref: str
    :param chrom_start: Dictionary with index positions where the chromosomes start
    :type chrom_start: dict
    :return: List with segmented data arrays, list with dictionaries which map from gene name to data array
    """
    bw_gen_mapping = []
    t_dict_list = []
    for data in all_values:
        gm, t_d = annotate(data, bed_ref, chrom_start)
        bw_gen_mapping.append(gm)
        t_dict_list.append(t_d)

    return bw_gen_mapping, t_dict_list


def rescale_all(transcript_data, vec_len=1000):
    """
    Wrapper function to rescale all data arrays in a list
    :param transcript_data: List with data arrays
    :type transcript_data: list(numpy.array)
    :param vec_len: Length to which all arrays are rescaled
    :type vec_len: int
    :return: numpy.array with rescaled data arrays (two-dim numpy array)
    """
    vec_len = vec_len if vec_len is not None else len(max(transcript_data, key=len))
    data_array = np.zeros((len(transcript_data), vec_len))
    for num, td in enumerate(transcript_data):
        data_array[num] = rescale(td, vec_len=vec_len)

    return data_array


def get_values(bw_list):
    """
    Retrieve all data values from bigwig file and concatenate them together
    :param bw_list: List with bigwig files
    :type bw_list: list(bigWigFile)
    :return: List with one data array per bigwig file, dictionary with starting indices for chromosomes.
    """
    if len(bw_list) == 0:
        raise ValueError('List with bigwig objects must not be empty.')

    all_values = [[] for _ in bw_list]

    counter = 0
    chrom_start = {}
    for chrom, length in bw_list[0].chroms().items():
        chrom_start[chrom] = counter
        counter += length
        for num, bw in enumerate(bw_list):
            values = np.asarray(bw.values(chrom, 0, length))
            all_values[num].extend(np.nan_to_num(values, nan=0.0).tolist())

    all_values = [np.asarray(genome) for genome in all_values]
    return all_values, chrom_start

