# Data Handler for ChIP-seq Data
## Introduction
Libraries for analysing and processing ChIP-seq data are 
all over the place; while some libraries can handle some file types,
some others fail to load them but have advantageous additional
functionality. This small library aims to bundle these forces
to create an easy-to-use module for ChIP-seq signal manipulation and
pre-processing.

## Pre-Requirements and Installation
To install and use the `datahandler` you need to have Python3 (preferably
Python >= 3.6) and pip installed. If you want to contribute to the development
or add some functionality yourself, you can clone this repository and install
the necessary libraries with
```bash
python3 -m pip install .
```
We currently haven't released a pip package yet. 
If you've set this library as a dependency for your Python script you can install
the newest stable version via
```bash
python3 -m pip install -e git+https://git@github.com/leoTiez/seqDataHandler.git#egg=seqDataHandler
```
or if you rather install a particular version 
```bash
python3 -m pip install -e git+https://git@github.com/leoTiez/seqDataHandler.git@RELEASETAG#egg=seqDataHandler
```
where `RELEASETAG` is replaced by the version number. The releases can also be downloaded
separately from the GitHub release page. The experimental version that is
currently under development (which often coincide with the `main` or stable version) can
be installed with 
```bash
python3 -m pip install -e git+https://git@github.com/leoTiez/seqDataHandler.git@develop#egg=seqDataHandler
```

## Usage
### Load Files
To load files import the `reader` module via

```python
from datahandler import reader
```

We currently provide loading
- bigwig of bigbed files `load_big_file` (which returns a `bigWigFile` object, see https://github.com/deeptools/pyBigWig),
- bam or bed files `load_bam_bed_file` (which returns a `BedTool` object, see https://daler.github.io/pybedtools/topical-documentation-contents.html), and
- fasta or fastq `load_fast` (which returns a `SeqIO` object, see https://biopython.org/wiki/SeqIO).

They all follow the same layout

```python
data_object = reader.loader_function(name, rel_path='data', is_abs_path=False)
```

where `loader_function` is replaced by one of the function names given above;
`name` is the file name or the path to the file if `is_abs_path` set to `True`;
and `rel_path` is the relative path from your current directory. Naturally
`rel_path` is ignored when `is_abs_path` is set to `True`, as the absolute path is
then given already via `name`.

`load_fast` has one additional flag `is_fastq` to determine whether the file
that is to be loaded is a fastq file.


The `reader` module also provides a function to create a bed file with random 
chunks called 

```python
reader.create_bed_random_fragments(bw, max_chunk=6000, name='random_fragments', path='/')
```

which comes handy for creating negative control setups for determining trends.
`bw` denotes a loaded bigwig file; `max_chunk` defines the maximal length of a 
piece of DNA; `name` is the name of the newly created bed file (note that the
file suffix `.bed` is automatically attached and doesn't need to be included
in the `name` parameter); and the path to the directory where the bed file is to
be saved.

### Handle Sequential ChIP-seq Data
We provide some fundamental transformation functions that can be applied
to the bigwig data object which are implemented inthe `seqDataHandler`.
It can be imported via

```python
from datahandler import seqDataHandler
``` 

Due to convenience, most functions expect a `numpy` array as input parameter.
To retrieve the values--chromosomes are concatenated--from the bigwig data object,
run

```python
all_values, chrom_dict = seqDataHandler.get_values(bw_list)
``` 

where `bw_list` is a list with bigwig data objects. The function returns a list 
of `numpy.array`s and a dictionary with the indices where in the array the chromosomes
start. If you want to retrieve the values from a single bigwig file,
run `get_values([bw_file])[0]`.

Most of the other functions are available for single for either single numpy arrays
or for a list of numpy arrays. They can be discriminated via the `_all` suffix in 
the function name.

Firstly, we provide two normalisation functions.

```python
norm_data = seqDataHandler.center_norm[_all](data[_list])
```
and
```python
norm_data = seqDataHandler.remap_norm[_all](data[_list])
```

where the brackets represent the extension that is necessary to apply the function
to a list of numpy arrays. The `center` normalisation returns the data with 
zero mean and a standard deviation of one. The `remap` normalisation procedure
rescales the data values such that the data minimum is zero and the data maximum is one.
For most ChIP-seq data we recommend to use the remap normalisation; however
keep in mind that this normalisation method does not replace biologically plausible
normalisation methods (spike normalisation or normalisation with qPCR).

Secondly, the data can be smoothed through

```python
smoothed_data = seqDataHandler.smooth[_all](data[_list], smooth[_list])
```

where the brackets, again, represent the suffix for the list function. The `smooth`
parameter indicates the size of your sliding smoothing window. If you want to use
the list function (`_all`) then it is necessary to apply one smoothing factor per 
array in the list; if one array is not to be smoothed, pass `None`.

Thirdly, the data can be annotated if a suitable bed file and the dictionary
with the starting indices for the chromosomes (`chrom_start`) are passed. 

```python
transcript[_list], trans_dict[_list] = seqDataHandler.annotate[_all](data[_list], bed_ref, chrom_start)
```

As before, the suffix in the brackets indicate the extension for the function to 
apply for several data arrays simultaneously; the `transcript[_list]` is a list (of lists
when `[_all]`) with the transcripts in the same order as they're defined in the bed file;
`trans_dict[_list]` is a dictionary (a list of dictionaries) that map from the transcript
name to the sequence; `data[_list]` is the data array (list of data arrays); `bed_ref` denotes
the bed file; and `chrom_start` is the directory with the chromosome start indices.

Lastly, for sometimes it is beneficial to rescale data arrays (for example the transcripts
retrieved from the `annotate` method) to let them match sizse, which can be used, for
example, for plotting or for other machine learning approaches (e.g. neural networks).
Rescale the data via

```python
rescaled[_list] = seqDataHandler.rescale[_all](transcript[_list], vec_len=1000)
``` 

If the `[_all]` method is used, it returns a two-dimensional numpy array `n x m`,
with `n` representing the number of data arrays in the list (e.g. transcripts),
and `m` is equal to the `vec_lenght`; if the method is used without the naming
suffix, then a single array is passed with the length `m = vec_lenght`.

### Process the ChIP-seq data
The `preprocess` library is currently under development and contains in its recent
status only two functions. Import it via

```python
from datahandler import preprocessing
```

To detect peaks in the ChIP-seq signal, we provide
a straight forward method that searches for maxima within a certain range.
We recommend to set this range to the signal resolution. 

```python
peaks = preprocessing.peak_detect_smooth(sig, peak_range=200, mode='wrap')
```

The function returns the indices with the peaks; `sig` denotes the input signal;
`peak_range` is the area in which the value needs to be a maxima; and `mode`
defines how boundary affects are managed. Pass `'wrap'` if the the peak window continues
at the end (beginning) and involves the in the beginning (end). Choose this variant
if `sig` represents concatenated chromosomes. If transcripts are passed, use `'clip'`
(you can read more about the different modes here
https://numpy.org/devdocs/reference/generated/numpy.take.html#numpy.take).

CPD signals should be noise filtered, knowing that they can only occur at
two adjacent pyrimidines (both at the transcribing and non-transcribing strand).
Whenever there's a signal at a position where there are no adjacent pyrimidines,
the signal can be muted. Run the filtering method via

```python
filtered_cpd = preprocessing.cancel_noise_cpd(cpd_sig, chrom_start, dna_seq)
```

where `cpd_sig` denotes the ChIP-seq signal for CPDs; `chrom_start` is the
directory with the start indices for the chromosomes; and `dna_seq` represents
the DNA sequence loaded from a fasta or fastq file (see the reading functions above).
It's important to keep in mind that the method replaces the signal in place, despite
the fact that it has a return value.

